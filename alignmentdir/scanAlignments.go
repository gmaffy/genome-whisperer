package alignmentdir

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// ============================= 1. Scan all dirs to see the stage of each sample ==================================== //

type FileInfo struct {
	Path         string
	Size         int64
	Present      bool
	Valid        bool
	ValidateErr  error
	IndexPresent bool
	IndexSize    int64
}

type SampleBamState struct {
	Sample    string
	SortedBam FileInfo
	RgmdBam   FileInfo
	RgmdCram  FileInfo
	BqsrBam   FileInfo
	BqsrCram  FileInfo

	OtherFiles []FileInfo
}

func isUsable(f FileInfo) bool {
	return f.Present && f.Valid
}

func hasIndex(f FileInfo) bool {
	return f.IndexPresent && f.IndexSize > 0
}

// ------------------------------------ Shared completeness / plan decision ----------------------------------------- //

// SampleAssessment is what AlignReads (data-directory mode) would decide about a
// sample: whether it is finished, and if not, the steps it would run.
type SampleAssessment struct {
	Sample   string
	LongRead bool
	Complete bool
	// External is the externally produced alignment a long-read sample would
	// adopt, if any.
	External string
	// AmbiguousExternal lists the candidates when there is more than one;
	// AlignReads refuses to choose and skips the sample.
	AmbiguousExternal []string
	Plan              []ProcessingReason
	NeedsReads        bool
}

// assessSample decides whether a sample is complete and, if not, plans its
// remaining steps. RunAlignReadsDir and ScanAlignments both call it, so the
// scan report always matches what a run would do.
func assessSample(state SampleBamState, bqsr, gatkNeedsBam bool) SampleAssessment {
	a := SampleAssessment{Sample: state.Sample, LongRead: IsLongReadSample(state.Sample)}

	// Adoption is decided from names alone; the header reads and the full
	// validation belong in the worker, inside ReasonAdoptExternalAlignment.
	if a.LongRead {
		candidates := adoptionCandidates(state)
		if len(candidates) > 1 {
			// Two movies, or a movie and a stale partial. Picking one silently
			// drops half a sample's data, so this is a question for a person;
			// nothing in the directory is touched.
			a.AmbiguousExternal = candidates
			return a
		}
		if len(candidates) == 1 {
			a.External = candidates[0]
		}
	}

	rgmdCramOK := isUsable(state.RgmdCram) && hasIndex(state.RgmdCram)
	bqsrCramOK := isUsable(state.BqsrCram) && hasIndex(state.BqsrCram)
	noIntermediates := !state.SortedBam.Present && !state.RgmdBam.Present && !state.BqsrBam.Present

	if a.LongRead {
		// A long-read sample is finished at the rgmd.cram: BQSR is never run on
		// it, so requiring a bqsr.cram would mean no long-read sample was ever
		// complete.
		//
		// OtherFiles is not required to be empty either. An adopted alignment
		// stays in the directory permanently — it is an input, and cleanup will
		// not remove it — so demanding its absence would re-plan and
		// re-validate every adopted sample on every run, which on a
		// multi-hundred-gigabyte CRAM is not free.
		a.Complete = rgmdCramOK && noIntermediates
	} else {
		a.Complete = rgmdCramOK && bqsrCramOK && noIntermediates && len(state.OtherFiles) == 0
	}
	if a.Complete {
		return a
	}

	// Build the full ordered plan from the current state so processSample can
	// execute every step without re-scanning the filesystem.
	a.Plan, a.NeedsReads = buildPlan(state, planInputs{
		bqsr:         bqsr && !a.LongRead,
		gatkNeedsBam: gatkNeedsBam,
		longRead:     a.LongRead,
		external:     a.External != "",
	})
	return a
}

// ------------------------------------------------- Scan ----------------------------------------------------------- //

// SampleStatus is the one-word verdict on a sample in a scan report.
type SampleStatus string

const (
	StatusComplete  SampleStatus = "COMPLETE"
	StatusCleanup   SampleStatus = "CLEANUP"   // outputs are done; only cleanup remains
	StatusResume    SampleStatus = "RESUME"    // an alignment exists; later steps remain
	StatusAlign     SampleStatus = "ALIGN"     // nothing usable; must align from reads
	StatusAmbiguous SampleStatus = "AMBIGUOUS" // several external alignments; needs a person
	StatusError     SampleStatus = "ERROR"     // the bams directory could not be read
)

// SampleReport is one sample in a scan.
type SampleReport struct {
	SampleAssessment
	Project string
	BamDir  string
	State   SampleBamState
	Err     error
}

// Status summarises the sample in one word.
func (r SampleReport) Status() SampleStatus {
	switch {
	case r.Err != nil:
		return StatusError
	case len(r.AmbiguousExternal) > 0:
		return StatusAmbiguous
	case r.Complete:
		return StatusComplete
	case containsReason(r.Plan, ReasonAlignFromReads), containsReason(r.Plan, ReasonAlignLongReads):
		return StatusAlign
	case len(r.Plan) == 1 && r.Plan[0] == ReasonCleanup:
		return StatusCleanup
	default:
		return StatusResume
	}
}

// ScanReport is the result of ScanAlignments.
type ScanReport struct {
	DataDir      string
	Species      string
	RefVer       string
	Reference    string
	BQSR         bool
	GatkNeedsBam bool
	Samples      []SampleReport
}

// ScanAlignments inspects every sample's bams directory under
//
//	<dataDir>/<species>/<project>/<sample>/reference_genomes/<refVer>/bams/
//
// with the same file rules and validation AlignReads uses, and reports what
// AlignReads --bqsr=<bqsr> would do with each sample. It reads files but never
// changes them (apart from creating the reference's GATK dictionary if it is
// missing, as every other command does).
func ScanAlignments(dataDir, species, refVer, genomesDir, refFasta string, bqsr, verbose, quick bool) (ScanReport, error) {
	report := ScanReport{DataDir: dataDir, Species: species, RefVer: refVer, BQSR: bqsr}

	if dataDir == "" || species == "" || refVer == "" {
		return report, fmt.Errorf("--data-dir, --species and --ref-version are required")
	}
	dInfo, err := os.Stat(dataDir)
	if err != nil {
		return report, fmt.Errorf("data directory: %w", err)
	}
	if !dInfo.IsDir() {
		return report, fmt.Errorf("data directory %s is not a directory", dataDir)
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		return report, err
	}

	resolvedFasta, err := resolveRefFasta(refFasta, genomesDir, species, refVer)
	if err != nil {
		return report, fmt.Errorf("could not resolve reference fasta: %w", err)
	}
	if fi, err := os.Stat(resolvedFasta); err != nil || !fi.Mode().IsRegular() {
		return report, fmt.Errorf("reference fasta %s is not a regular file", resolvedFasta)
	}
	if err := utils.EnsureGatkDict(resolvedFasta); err != nil {
		return report, err
	}
	report.Reference = resolvedFasta
	report.GatkNeedsBam = gatkReadsNeedBam(resolvedFasta)

	pattern := filepath.Join(dataDirAbs, species, "*", "*", "reference_genomes", refVer, "bams")
	matches, _ := filepath.Glob(pattern)
	if len(matches) == 0 {
		return report, fmt.Errorf("no bams directories match %s", pattern)
	}
	sort.Strings(matches)

	for i, bamDir := range matches {
		// <dataDirAbs>/<species>/<project>/<sample>/reference_genomes/<refVer>/bams
		sampleDir := filepath.Dir(filepath.Dir(filepath.Dir(bamDir)))
		sample := filepath.Base(sampleDir)
		project := filepath.Base(filepath.Dir(sampleDir))
		fmt.Fprintf(os.Stderr, "[%d/%d] checking %s/%s\n", i+1, len(matches), project, sample)

		r := SampleReport{Project: project, BamDir: bamDir}
		r.Sample = sample
		state, scanErr := inspectSampleBamDir(sample, bamDir, resolvedFasta, verbose, quick)
		r.State = state
		if scanErr != nil {
			r.Err = scanErr
		} else {
			r.SampleAssessment = assessSample(state, bqsr, report.GatkNeedsBam)
		}
		report.Samples = append(report.Samples, r)
	}
	return report, nil
}

// -------------------------------------------------- Report -------------------------------------------------------- //

var reasonText = map[ProcessingReason]string{
	ReasonAlignFromReads:         "align reads from clean_reads",
	ReasonAlignLongReads:         "align long reads with pbmm2",
	ReasonAdoptExternalAlignment: "adopt the external alignment",
	ReasonMarkDupSortedBam:       "mark duplicates (sorted.bam to RGMD.bam)",
	ReasonIndexRgmdBam:           "index RGMD.bam",
	ReasonConvertRgmdBam:         "convert RGMD.bam to RGMD.cram",
	ReasonIndexRgmdCram:          "index RGMD.cram",
	ReasonMaterializeRgmdBam:     "decode RGMD.cram to RGMD.bam for GATK",
	ReasonRunBQSR:                "run BQSR",
	ReasonConvertBqsrBam:         "convert bqsr.bam to bqsr.cram",
	ReasonIndexBqsrCram:          "index bqsr.cram",
	ReasonCleanup:                "remove intermediates",
}

func describePlan(plan []ProcessingReason) string {
	parts := make([]string, len(plan))
	for i, p := range plan {
		if t, ok := reasonText[p]; ok {
			parts[i] = t
		} else {
			parts[i] = string(p)
		}
	}
	return strings.Join(parts, " → ")
}

func humanSize(n int64) string {
	const unit = 1024
	if n < unit {
		return fmt.Sprintf("%d B", n)
	}
	div, exp := int64(unit), 0
	for m := n / unit; m >= unit; m /= unit {
		div *= unit
		exp++
	}
	return fmt.Sprintf("%.1f %ciB", float64(n)/float64(div), "KMGTPE"[exp])
}

func describeFile(f FileInfo) string {
	if !f.Present {
		return "-"
	}
	var s string
	switch {
	case !f.Valid && f.ValidateErr != nil:
		s = "INVALID: " + strings.TrimSpace(f.ValidateErr.Error())
	case !f.Valid:
		s = "INVALID"
	case hasIndex(f):
		s = "ok, indexed"
	default:
		s = "ok, NO INDEX"
	}
	return fmt.Sprintf("%-40s %9s  %s", filepath.Base(f.Path), humanSize(f.Size), s)
}

var statusColor = map[SampleStatus]*color.Color{
	StatusComplete:  color.New(color.FgGreen, color.Bold),
	StatusCleanup:   color.New(color.FgCyan, color.Bold),
	StatusResume:    color.New(color.FgYellow, color.Bold),
	StatusAlign:     color.New(color.FgMagenta, color.Bold),
	StatusAmbiguous: color.New(color.FgRed, color.Bold),
	StatusError:     color.New(color.FgRed, color.Bold),
}

// WriteScanReport prints a per-sample breakdown followed by totals.
func WriteScanReport(w io.Writer, r ScanReport) {
	bqsr := "with BQSR"
	if !r.BQSR {
		bqsr = "without BQSR"
	}
	fmt.Fprintf(w, "\nAlignment scan: %s / %s (%s)\n", r.Species, r.RefVer, r.DataDir)
	fmt.Fprintf(w, "Reference: %s\n", r.Reference)
	fmt.Fprintf(w, "Plans show what AlignReads would do %s.\n", bqsr)
	if r.GatkNeedsBam {
		fmt.Fprintln(w, "Reference has contigs past the BAI limit: BQSR reads the RGMD.bam, not the RGMD.cram.")
	}
	fmt.Fprintln(w, strings.Repeat("─", 100))

	counts := map[SampleStatus]int{}
	var blocked []string
	for _, s := range r.Samples {
		st := s.Status()
		counts[st]++

		kind := "short read"
		if s.LongRead {
			kind = "long read"
		}
		fmt.Fprintf(w, "%s %s/%s (%s)\n", statusColor[st].Sprintf("%-9s", st), s.Project, s.Sample, kind)

		if s.Err != nil {
			fmt.Fprintf(w, "    error: %v\n\n", s.Err)
			continue
		}
		for _, f := range []struct {
			label string
			info  FileInfo
		}{
			{"sorted.bam", s.State.SortedBam},
			{"RGMD.bam", s.State.RgmdBam},
			{"RGMD.cram", s.State.RgmdCram},
			{"bqsr.bam", s.State.BqsrBam},
			{"bqsr.cram", s.State.BqsrCram},
		} {
			if f.info.Present {
				fmt.Fprintf(w, "    %-10s %s\n", f.label, describeFile(f.info))
			}
		}
		for _, o := range s.State.OtherFiles {
			note := "unrecognised"
			if o.Path == s.External {
				note = "external alignment (will be adopted, then kept)"
			}
			fmt.Fprintf(w, "    %-10s %-40s %9s  %s\n", "other", filepath.Base(o.Path), humanSize(o.Size), note)
		}

		switch st {
		case StatusAmbiguous:
			fmt.Fprintln(w, "    more than one external alignment; AlignReads will skip this sample until only one remains")
		case StatusComplete:
		default:
			fmt.Fprintf(w, "    next: %s\n", describePlan(s.Plan))
			if s.NeedsReads {
				fmt.Fprintln(w, "    needs the sample's reads in clean_reads/")
			}
			if !s.LongRead && st == StatusCleanup && len(s.State.OtherFiles) > 0 {
				// Cleanup removes only the pipeline's own intermediates, so
				// these files keep the sample from ever counting as complete.
				fmt.Fprintln(w, "    unrecognised files block completion; cleanup does not remove them — move or delete them by hand")
				blocked = append(blocked, s.Project+"/"+s.Sample)
			}
		}
		fmt.Fprintln(w)
	}

	fmt.Fprintln(w, strings.Repeat("─", 100))
	fmt.Fprintf(w, "Samples: %d\n", len(r.Samples))
	for _, st := range []struct {
		s    SampleStatus
		desc string
	}{
		{StatusComplete, "finished, nothing to do"},
		{StatusCleanup, "outputs done, only intermediates to remove"},
		{StatusResume, "alignment exists, later steps remain"},
		{StatusAlign, "nothing usable, align from reads"},
		{StatusAmbiguous, "several external alignments, needs a person"},
		{StatusError, "bams directory could not be read"},
	} {
		if counts[st.s] > 0 {
			fmt.Fprintf(w, "  %s %4d  %s\n", statusColor[st.s].Sprintf("%-9s", st.s), counts[st.s], st.desc)
		}
	}
	if len(blocked) > 0 {
		fmt.Fprintf(w, "\nBlocked by unrecognised files: %s\n", strings.Join(blocked, ", "))
	}
}
