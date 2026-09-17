package alignmentdir

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
)

// Long-read support for the data directory pipeline.
//
// A long-read sample reaches the same durable artefacts as a short-read one —
// <sample>.RGMD.cram plus its index — by one of two routes:
//
//	one FASTQ        --pbmm2-->            <sample>.sorted.bam --pbmarkdup--> ...
//	X.aligned.cram   --addreplacerg-->     <sample>.sorted.bam --pbmarkdup--> ...
//
// Both normalise onto <sample>.sorted.bam, which is what lets the existing
// markdup → cram → index chain serve all three sample shapes unchanged, and
// what makes an interrupted run resumable through the ordinary scan: that name
// classifies as slotSortedBam, so the next run picks up at markdup.

// isFastqRead reports whether a filename is a FASTQ.
//
// A duplicate of warehouse.IsFastq, which cannot be imported: warehouse/scan.go
// imports this package, so the dependency only runs one way. Keep the two in
// step — the original is warehouse/standard.go.
func isFastqRead(name string) bool {
	lower := strings.ToLower(name)
	return strings.HasSuffix(lower, ".fastq.gz") ||
		strings.HasSuffix(lower, ".fq.gz") ||
		strings.HasSuffix(lower, ".fastq") ||
		strings.HasSuffix(lower, ".fq")
}

// GetReadsLong returns the single long-read FASTQ in a clean_reads directory.
//
// Identified by being the only FASTQ present, not by role. Identifying it by
// elimination — the one file ClassifyRead calls ReadRoleUnknown — looks right
// and does not work: a PacBio movie name ends in a six-digit timestamp, so
// m64011_200521_073741.fastq.gz satisfies isFwd's "1.fastq.gz" test and
// classifies as a forward read. Every real HiFi sample would report no reads.
//
// Exactly one file is required. Two would be two movies, or two subsets of one,
// and choosing between them is not this function's decision to make.
func GetReadsLong(cleanReadsDir string) (string, error) {
	entries, err := os.ReadDir(cleanReadsDir)
	if err != nil {
		return "", err
	}

	var found []string
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if isFastqRead(entry.Name()) {
			found = append(found, entry.Name())
		}
	}

	switch len(found) {
	case 1:
		return filepath.Join(cleanReadsDir, found[0]), nil
	case 0:
		return "", fmt.Errorf("no FASTQ found in %s", cleanReadsDir)
	default:
		// Print the role each file classified as: where a movie name was
		// mistaken for a paired read, this is the line that says so.
		described := make([]string, 0, len(found))
		for _, name := range found {
			described = append(described, fmt.Sprintf("%s (classifies as %s)", name, ClassifyRead(name)))
		}
		return "", fmt.Errorf("expected one long-read FASTQ in %s, found %d: %s",
			cleanReadsDir, len(found), strings.Join(described, ", "))
	}
}

// adoptableAlignment reports whether a bams-directory filename is an externally
// produced alignment eligible for adoption: a BAM or CRAM this pipeline neither
// writes nor recognises.
//
// Defined in terms of classifyAlignmentFile rather than by listing suffixes, so
// that everything adoptable is by construction RoleOther to the warehouse
// scanner and the estate report's view of these files does not change.
//
// Names deliberately need not contain the sample: warehouse/inspect.go records
// that a sample directory named MO971_LR holds MENINA_LONG_READS.* files.
//
// Two kinds of name are excluded because they are this pipeline's own leavings
// rather than anybody's input: *.RG.bam, which the old single-sample pbmm2 path
// wrote, and anything carrying .partial, which marks a write that was
// interrupted. A partial is doubly wrong to adopt — it is ours, and it is
// truncated.
func adoptableAlignment(name string) bool {
	lower := strings.ToLower(name)
	if !strings.HasSuffix(lower, ".bam") && !strings.HasSuffix(lower, ".cram") {
		return false
	}
	if strings.HasSuffix(lower, ".rg.bam") || strings.Contains(lower, ".partial") {
		return false
	}
	return classifyAlignmentFile(name) == slotOther
}

// adoptionCandidates returns the external alignments eligible for adoption.
//
// Empty unless the sample has nothing usable of its own: the pipeline's own
// RGMD.cram, RGMD.bam or sorted.bam always wins over an external file, and it
// wins by slot rather than by directory order. That matters because
// inspectSampleBamDir assigns slots instead of accumulating them, so two files
// competing for one slot resolve to whichever sorts last — stable enough to
// look deliberate, arbitrary enough to be wrong. An external alignment and a
// .RGMD.cram are never in the same slot, so this precedence is decided by a
// switch on state and not by a filename.
func adoptionCandidates(state SampleBamState) []string {
	if isUsable(state.RgmdCram) || isUsable(state.RgmdBam) || isUsable(state.SortedBam) {
		return nil
	}

	// Names this pipeline writes for this sample are never adoptable, whatever
	// state they are in: they are leftovers from an interrupted run of our own,
	// and cleanup collects them.
	ours := make(map[string]struct{})
	for _, name := range pipelineIntermediates(state.Sample) {
		ours[name] = struct{}{}
	}

	var candidates []string
	for _, f := range state.OtherFiles {
		base := filepath.Base(f.Path)
		if _, mine := ours[base]; mine {
			continue
		}
		if adoptableAlignment(base) {
			candidates = append(candidates, f.Path)
		}
	}
	return candidates
}

// alignLongReads aligns one long-read FASTQ to sortedBam with pbmm2.
//
// Written under a temporary name and renamed, because the alternative is
// silent: utils.ValidateBam's deep form is "samtools view -h", which treats a
// missing BGZF end-of-file marker as a warning rather than an error. A pbmm2
// killed on a block boundary therefore leaves a file that the deep check
// passes, buildPlan accepts as a usable sorted.bam, and marks duplicates on —
// producing a complete-looking cram with reads missing. A file under a name no
// slot recognises cannot be mistaken for a finished one.
func alignLongReads(sePath, refIndex, sortedBam, sample, preset string, threads int, verbose bool) error {
	// The suffix order matters: the temporary name must end in .bam so pbmm2
	// and samtools sort both infer the format, and must not end in
	// "sorted.bam", which classifyAlignmentFile would put in the sorted slot.
	partial := strings.TrimSuffix(sortedBam, ".bam") + ".partial.bam"

	if err := clearAlignment(sortedBam); err != nil {
		return err
	}
	if err := clearAlignment(partial); err != nil {
		return err
	}

	lib := fmt.Sprintf("%s_1", sample)
	if err := alignment.Pbmm2Align(sePath, refIndex, partial, sample, lib, preset, threads, verbose); err != nil {
		_ = os.Remove(partial)
		return err
	}
	return os.Rename(partial, sortedBam)
}

// clearAlignment removes an alignment file and every index spelling that might
// be sitting beside it.
//
// The index list is indexCandidates rather than a hand-written one: BamIndex
// writes a CSI for a BAM, so removing only the .bai — as the paired-read path
// used to — leaves a stale .csi that findIndex will happily validate and reuse
// against the different data the re-alignment just wrote.
func clearAlignment(path string) error {
	targets := append([]string{path}, indexCandidates(path)...)
	return removeIfExists(targets...)
}

// verifyAdoptable checks an externally produced alignment hard enough to build
// a call set on.
//
// Deliberately stricter than checkAlignmentContigs, which passes a file whose
// @SQ names merely form a subset of the reference's. That leniency is right for
// its own purpose — diagnosing a renamed header — but here it would let a cram
// aligned to a different build of a same-named assembly through, and the
// resulting calls would be wrong with nothing on disk to show why.
func verifyAdoptable(path, refFasta string) error {
	out, err := exec.Command("samtools", "view", "-H", path).Output()
	if err != nil {
		return fmt.Errorf("reading the header of %s: %w", path, err)
	}
	header := string(out)

	if !headerIsCoordinateSorted(header) {
		return fmt.Errorf("%s is not coordinate sorted (@HD SO:coordinate); sort it before adopting it", path)
	}
	if !strings.Contains(header, "@RG\t") {
		return fmt.Errorf("%s has no @RG line; adoption cannot establish which sample its reads belong to", path)
	}

	refSQ, err := referenceSQ(refFasta)
	if err != nil {
		return fmt.Errorf("reading the sequence dictionary for %s: %w", refFasta, err)
	}
	fileSQ := parseSQ(header)
	if len(fileSQ) == 0 {
		return fmt.Errorf("%s has no @SQ lines; it is not aligned to anything", path)
	}

	byName := make(map[string]sqRecord, len(refSQ))
	for _, rec := range refSQ {
		byName[rec.name] = rec
	}

	for _, rec := range fileSQ {
		refRec, ok := byName[rec.name]
		if !ok {
			// sameSequences exists precisely because this is usually a rename
			// rather than the wrong assembly, and a rename is fixable without
			// re-aligning. Say which it is, and leave the file alone.
			if sameSequences(fileSQ, refSQ) {
				return fmt.Errorf("%s holds the same sequences as the reference under different names: "+
					"rename its header with scripts/reheader-crams.sh -r %s %s and re-run — the file has not been touched",
					path, refFasta, path)
			}
			return fmt.Errorf("%s was aligned to a different assembly (contig %q is not in the reference); "+
				"it has not been touched", path, rec.name)
		}
		if rec.length != refRec.length {
			return fmt.Errorf("%s contig %q is %s bases, the reference's is %s: a different build of the same assembly; "+
				"it has not been touched", path, rec.name, rec.length, refRec.length)
		}
		if rec.m5 != "" && refRec.m5 != "" && rec.m5 != refRec.m5 {
			return fmt.Errorf("%s contig %q has a different sequence checksum from the reference's: "+
				"a different build of the same assembly; it has not been touched", path, rec.name)
		}
	}
	return nil
}

// headerIsCoordinateSorted reports whether a SAM header declares a coordinate
// sort order.
func headerIsCoordinateSorted(header string) bool {
	for _, line := range strings.Split(header, "\n") {
		if !strings.HasPrefix(line, "@HD\t") {
			continue
		}
		for _, field := range strings.Split(line, "\t") {
			if field == "SO:coordinate" {
				return true
			}
		}
	}
	return false
}

// normaliseAdoptedAlignment writes an adopted alignment out as this sample's
// sorted.bam, with a read group naming this sample.
//
// One samtools pass does both jobs, and both are necessary:
//
//   - BAM, because pbmarkdup reads bam/xml/fa/fq and not CRAM.
//   - The read group, because the adopted file's SM is the movie or the
//     collaborator's own sample name — MENINA_LONG_READS where the sample
//     directory says MO971_LR — and per warehouse/inspect.go that mismatch is
//     the normal case, not the exception. pbmarkdup preserves @RG, so left
//     alone it would reach the gVCF's sample column and swap the sample's
//     identity in the joint call set.
//
// The adopted file is read and never written: it is an input, and where no
// long-read FASTQ was kept it is the only copy of the sample's reads.
func normaliseAdoptedAlignment(adopted, sortedBam, refFasta, sample string, threads int, verbose bool) error {
	partial := strings.TrimSuffix(sortedBam, ".bam") + ".partial.bam"

	if err := clearAlignment(sortedBam); err != nil {
		return err
	}
	if err := clearAlignment(partial); err != nil {
		return err
	}

	rg := fmt.Sprintf(`@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:PACBIO`, sample, sample, sample)
	cmdStr := fmt.Sprintf(`samtools addreplacerg -@ %d -m overwrite_all -r '%s' --reference %s -O bam -o %s %s`,
		samThreads(threads), rg, refFasta, partial, adopted)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	if err != nil {
		_ = os.Remove(partial)
		return fmt.Errorf("normalising %s into %s: %w", adopted, sortedBam, err)
	}
	return os.Rename(partial, sortedBam)
}

// checkMarkedDuplicates confirms a duplicate-marked long-read BAM is still an
// aligned, coordinate-sorted file.
//
// pbmarkdup accepts FASTA and FASTQ as well as BAM, and this pipeline hands it
// a BAM produced from a FASTQ, so its output is worth one header read. Without
// this the failure is silent: checkAlignmentContigs returns nil for a file with
// no @SQ lines at all — correct leniency for its own purpose — so an unaligned
// or re-ordered output would validate clean, convert to cram, index, and reach
// HaplotypeCaller looking finished.
func checkMarkedDuplicates(rgmdBam string) error {
	out, err := exec.Command("samtools", "view", "-H", rgmdBam).Output()
	if err != nil {
		return fmt.Errorf("reading the header of %s: %w", rgmdBam, err)
	}
	header := string(out)

	if len(parseSQ(header)) == 0 {
		return fmt.Errorf("%s has no @SQ lines after duplicate marking: the reads came back unaligned", rgmdBam)
	}
	if !headerIsCoordinateSorted(header) {
		return fmt.Errorf("%s is not coordinate sorted after duplicate marking", rgmdBam)
	}
	return nil
}

// samThreads bounds a samtools -@ value the same way alignment.samThreads does:
// the codecs stop scaling much past 8, and this is a per-sample budget several
// samples spend at once.
func samThreads(threads int) int {
	if threads < 1 {
		return 1
	}
	if threads > 8 {
		return 8
	}
	return threads
}

// pbmm2Presets are the alignment modes pbmm2 accepts, as reported by
// `pbmm2 align --help` (pbmm2 26.2.0).
var pbmm2Presets = []string{"SUBREAD", "CCS", "HIFI", "ISOSEQ", "UNROLLED"}

// validatePreset rejects an unknown --preset before a run starts.
//
// pbmm2 interpolates the value straight into its own choice check and fails,
// which without this happens after the short-read samples have already spent
// hours. Worth noting: the flag's help text long read "CSS", a transposition of
// CCS, so a value copied from the help itself was invalid.
func validatePreset(preset string) error {
	for _, p := range pbmm2Presets {
		if strings.EqualFold(preset, p) {
			return nil
		}
	}
	return fmt.Errorf("preset %q is not a pbmm2 alignment mode; choose one of %s",
		preset, strings.Join(pbmm2Presets, ", "))
}
