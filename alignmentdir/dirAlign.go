package alignmentdir

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"

	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
)

func isFwd(filename string) bool {
	return strings.HasSuffix(filename, "1.fastq.gz") ||
		strings.HasSuffix(filename, "1.fq.gz") ||
		strings.HasSuffix(filename, "1.fastq") ||
		strings.HasSuffix(filename, "1.fq") ||
		strings.HasPrefix(filename, "Forward") ||
		strings.Contains(filename, "_R1_") ||
		strings.Contains(filename, "_1_") ||
		strings.Contains(filename, "1.fq")
}

func isRev(filename string) bool {
	return strings.HasSuffix(filename, "2.fastq.gz") ||
		strings.HasSuffix(filename, "2.fq.gz") ||
		strings.HasSuffix(filename, "2.fastq") ||
		strings.HasSuffix(filename, "2.fq") ||
		strings.HasPrefix(filename, "Reverse") ||
		strings.Contains(filename, "_R2_") ||
		strings.Contains(filename, "_2_") ||
		strings.Contains(filename, "_2.fq")
}

func GetReadsPE(cleanReadsDir string) ([]string, []string, error) {
	var fwdReads, revReads []string

	entries, err := os.ReadDir(cleanReadsDir)
	if err != nil {
		return nil, nil, err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}

		name := entry.Name()
		fullPath := filepath.Join(cleanReadsDir, name)

		if isFwd(name) {
			fwdReads = append(fwdReads, fullPath)
		} else if isRev(name) {
			revReads = append(revReads, fullPath)
		}
	}

	sort.Strings(fwdReads)
	sort.Strings(revReads)

	return fwdReads, revReads, nil
}

type GenomeRef struct {
	RefVer    string
	FastaPath string
	DictPath  string
}

func GetValidGenomesFromDisk(genomesDir string) (map[string][]GenomeRef, error) {
	result := make(map[string][]GenomeRef)
	species, err := os.ReadDir(genomesDir)
	if err != nil {
		return nil, fmt.Errorf("reading genomes dir %s: %w", genomesDir, err)
	}

	for _, sp := range species {
		if !sp.IsDir() {
			continue
		}
		spKey := strings.ToUpper(sp.Name())
		spDir := filepath.Join(genomesDir, sp.Name())

		refs, err := os.ReadDir(spDir)
		if err != nil {
			continue
		}

		for _, ref := range refs {
			if !ref.IsDir() {
				continue
			}
			assemblyDir := filepath.Join(spDir, ref.Name(), "assembly")
			info, err := os.Stat(assemblyDir)
			if err != nil || !info.IsDir() {
				continue
			}

			entries, err := os.ReadDir(assemblyDir)
			if err != nil {
				continue
			}

			var fastaPath, gzFastaPath, dictPath string
			for _, f := range entries {
				if f.IsDir() {
					continue
				}
				name := f.Name()
				switch {
				case utils.IsBgzippedFasta(name):
					gzFastaPath = filepath.Join(assemblyDir, name)
				case utils.IsFasta(name):
					fastaPath = filepath.Join(assemblyDir, name)
				case strings.HasSuffix(name, ".dict"):
					dictPath = filepath.Join(assemblyDir, name)
				}
			}

			// An assembly directory holding both forms is used uncompressed:
			// GATK's recalibration tools refuse a block-compressed reference,
			// and ReadDir order alone would pick whichever name sorts last.
			if fastaPath == "" {
				fastaPath = gzFastaPath
			}

			if fastaPath != "" && dictPath != "" {
				result[spKey] = append(result[spKey], GenomeRef{
					RefVer:    strings.ToUpper(ref.Name()),
					FastaPath: fastaPath,
					DictPath:  dictPath,
				})
			}
		}
	}

	fmt.Printf("Discovered genomes: %d\n", len(result))

	return result, nil
}

func resolveRefFasta(refFasta, genomesDir, species, refVer string) (string, error) {
	if refFasta != "" {
		if _, err := os.Stat(refFasta); err != nil {
			return "", fmt.Errorf("provided refFasta %q not accessible: %w", refFasta, err)
		}
		return refFasta, nil
	}
	if genomesDir == "" {
		return "", fmt.Errorf("either refFasta or genomesDir must be provided")
	}
	genomes, err := GetValidGenomesFromDisk(genomesDir)
	if err != nil {
		return "", fmt.Errorf("auto-discovering genomes: %w", err)
	}
	spKey := strings.ToUpper(species)
	refs, ok := genomes[spKey]
	if !ok {
		return "", fmt.Errorf("no genomes found for species %q in %s", species, genomesDir)
	}
	verKey := strings.ToUpper(refVer)
	for _, r := range refs {
		if r.RefVer == verKey {
			color.Green("Auto-resolved reference: %s\n", r.FastaPath)
			return r.FastaPath, nil
		}
	}
	return "", fmt.Errorf("species %q found but no assembly matching refVer %q in %s", species, refVer, genomesDir)
}

// ProcessingReason describes an individual pipeline step for a sample.
type ProcessingReason string

const (
	// ReasonAlignFromReads — no usable intermediates; align fwd/rev reads from scratch.
	ReasonAlignFromReads ProcessingReason = "align_from_reads"
	// ReasonMarkDupSortedBam — sorted.bam present; run MarkDuplicates then convert to cram.
	ReasonMarkDupSortedBam ProcessingReason = "markdup_sorted_bam"
	// ReasonIndexRgmdBam — rgmd.bam present but has no index; index it before converting.
	ReasonIndexRgmdBam ProcessingReason = "index_rgmd_bam"
	// ReasonConvertRgmdBam — rgmd.bam present; convert to rgmd.cram.
	ReasonConvertRgmdBam ProcessingReason = "convert_rgmd_bam"
	// ReasonIndexRgmdCram — rgmd.cram present but missing index.
	ReasonIndexRgmdCram ProcessingReason = "index_rgmd_cram"
	// ReasonMaterializeRgmdBam — the rgmd.cram is the only surviving alignment
	// but GATK cannot read an indexed cram on this reference; decode it back to
	// a CSI-indexed rgmd.bam for the BQSR steps. ReasonCleanup removes it again.
	ReasonMaterializeRgmdBam ProcessingReason = "materialize_rgmd_bam"
	// ReasonRunBQSR — rgmd alignment complete; BQSR not yet done.
	ReasonRunBQSR ProcessingReason = "run_bqsr"
	// ReasonConvertBqsrBam — bqsr.bam present; convert to bqsr.cram.
	ReasonConvertBqsrBam ProcessingReason = "convert_bqsr_bam"
	// ReasonIndexBqsrCram — bqsr.cram present but missing index.
	ReasonIndexBqsrCram ProcessingReason = "index_bqsr_cram"
	// ReasonCleanup — all required crams present; remove intermediate files.
	ReasonCleanup ProcessingReason = "cleanup"
)

type sampleTask struct {
	sample        string
	bamDir        string
	cleanReadsDir string
	fwd           string
	rev           string
	state         SampleBamState
	plan          []ProcessingReason
}

type sampleResult struct {
	sample  string
	success bool
	err     error
}

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}

	if species == "" {
		fmt.Println("Please provide species name")
		return
	}
	if refVer == "" {
		fmt.Println("Please provide reference version name")
		return
	}

	resolvedFasta, resolveErr := resolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		color.Red("Could not resolve reference fasta: %v\n", resolveErr)
		return
	}

	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", resolvedFasta)
		return
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", resolvedFasta)
		return
	}

	dictFilePath := utils.DictPath(resolvedFasta)
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	if bqsr {
		if aligner == "pbmm2" {
			color.Red("BQSR is not supported for pbmm2. Use bwa-mem2 or disable BQSR.\n")
			return
		}
		if len(knownSites) == 0 && !bootstrap {
			color.Red("BQSR requested: either provide known-sites or enable bootstrap.\n")
			return
		}
		if len(knownSites) > 0 && bootstrap {
			color.Red("Choose either known-sites OR bootstrap, not both.\n")
			return
		}
		for _, ks := range knownSites {
			if _, err := os.Stat(ks); err != nil {
				color.Red("Known-sites file not found: %s\n", ks)
				return
			}
		}
	}

	color.Cyan("===================== Getting directories with the clean_reads directories =====================\n")
	cleanReadsDirs, err := getPathsWithCleanReadsDir(dataDir, species)
	if err != nil {
		color.Red("Error finding samples dirs: %v\n", err)
		return
	}
	if len(cleanReadsDirs) == 0 {
		color.Yellow("No samples dirs found under %s/%s\n", dataDir, species)
		return
	}
	color.Green("Sample dirs found: %d\n", len(cleanReadsDirs))
	var (
		tasks           []sampleTask
		alreadyComplete int
		scanFailures    []string
		missingReads    []string
		longReadSamples []string
	)

	// Decided once for the run: on a reference with contigs past the BAI limit
	// the GATK steps have to be given a bam, because htsjdk turns a .crai into a
	// BAI it cannot address. Said out loud because it changes both the plan and
	// the peak disk a sample needs.
	gatkNeedsBam := gatkReadsNeedBam(resolvedFasta)
	if gatkNeedsBam {
		color.Yellow("Reference has contigs longer than %d bases: the BQSR steps will read the rgmd.bam, not the rgmd.cram\n", baiMaxPosition)
		color.Yellow("(htsjdk builds a BAI from a .crai on open, and a BAI cannot address a position that far into a contig)\n\n")
	}

	color.Cyan("===================== Checking and validating alignment files currently in directory =====================\n")
	for _, cleanReadsDir := range cleanReadsDirs {
		sampleDir := filepath.Dir(cleanReadsDir)
		sample := filepath.Base(sampleDir)
		bamDir := filepath.Join(sampleDir, "reference_genomes", refVer, "bams")

		state, scanErr := inspectSampleBamDir(sample, bamDir, resolvedFasta, verbose, quick)
		if scanErr != nil {
			color.Red("[%s] Scanning the bams directory failed: %v\n", sample, scanErr)
			scanFailures = append(scanFailures, sample)
			continue
		}

		rgmdCramOK := state.RgmdCram.Present && state.RgmdCram.Valid && state.RgmdCram.IndexPresent && state.RgmdCram.IndexSize > 0
		bqsrCramOK := state.BqsrCram.Present && state.BqsrCram.Valid && state.BqsrCram.IndexPresent && state.BqsrCram.IndexSize > 0

		if rgmdCramOK && bqsrCramOK && !state.SortedBam.Present && !state.RgmdBam.Present && !state.BqsrBam.Present && len(state.OtherFiles) == 0 {
			alreadyComplete++
			color.Green("[%s]- ✅ PASS. Skipping ...\n", sample)
			continue
		}

		if strings.HasSuffix(strings.ToUpper(sample), "LR") {
			color.Blue("[%s] Long reads sample ..\n", sample)
			longReadSamples = append(longReadSamples, sample)
			continue
		}

		// Build the full ordered plan from the current state so processSample
		// can execute every step without re-scanning the filesystem.
		plan, needsReads := buildPlan(state, bqsr, gatkNeedsBam)

		var fwd, rev string
		if needsReads {
			fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
			if readsErr != nil || len(fwdReads) != 1 || len(revReads) != 1 {
				color.Red("[%s] Forward and reverse reads not found in %s\n", sample, cleanReadsDir)
				missingReads = append(missingReads, sample)
				continue
			}
			fwd = fwdReads[0]
			rev = revReads[0]
		}

		color.Yellow("[%s] Sample queued — %d steps: %v\n\n", sample, len(plan), plan)
		tasks = append(tasks, sampleTask{
			sample:        sample,
			bamDir:        bamDir,
			cleanReadsDir: cleanReadsDir,
			fwd:           fwd,
			rev:           rev,
			state:         state,
			plan:          plan,
		})
	}

	color.Green("Samples complete: %d\n", alreadyComplete)
	color.Green("Samples queued:   %d\n", len(tasks))
	fmt.Printf("\n-------------------------------------- Queued Samples --------------------------------------\n")
	for _, task := range tasks {
		color.Yellow("%s\n", task.sample)
		fmt.Printf("\nTo DO:\n-------------------------------------------------\n\n")
		for _, reason := range task.plan {
			fmt.Printf("%s\n", reason)
		}
		fmt.Printf("\n\n================================================================\n\n")
	}
	fmt.Printf("----------------------------------------------------------------------------------------------\n")

	if len(scanFailures) > 0 {
		color.Red("Scan failures:    %d\n", len(scanFailures))
	}
	if len(longReadSamples) > 0 {
		color.Yellow("Long read samples:      %d\n\n", len(longReadSamples))
	}

	if len(tasks) == 0 {
		return
	}

	totalCores := runtime.NumCPU()
	if threads <= 0 {
		threads = totalCores
	}
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}

	sem := make(chan struct{}, maxParallelJobs)
	results := make([]sampleResult, len(tasks))
	var wg sync.WaitGroup

	opts := newRuntimeOpts(threads, maxParallelJobs)

	color.Cyan("Processing %d samples using %d threads each (Max parallel jobs: %d)\n", len(tasks), threads, maxParallelJobs)
	color.Cyan("GATK: %d JVM slots, %s per sample step, %s per shard, %d interval shards, %d pair-HMM threads\n\n",
		cap(gatkSlots), opts.javaOpts, opts.shardJava, opts.shardCount, opts.pairHmm)
	for i, task := range tasks {
		wg.Add(1)
		go func(idx int, task sampleTask) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()

			fmt.Printf("\nProcessing sample: %s ....\n\n", task.sample)
			results[idx] = processSample(task, resolvedFasta, gatkLogLevel, aligner, verbose, quick, skipVer, bqsr, bootstrap, knownSites, opts)
			if results[idx].success {
				color.Green("[%s] done\n", task.sample)
				return
			}
			color.Red("[%s] failed: %v\n", task.sample, results[idx].err)
		}(i, task)
	}

	wg.Wait()

	var failed []string
	var successfulSampes []string
	successful := 0
	for _, result := range results {
		if result.success {
			successfulSampes = append(successfulSampes, result.sample)
			successful++
			continue
		}
		failed = append(failed, result.sample)
	}

	fmt.Println()
	color.Green("Processed successfully: %d\n", successful)
	color.Green("Successful samples: %s\n", strings.Join(successfulSampes, ", "))
	if len(failed) > 0 {
		color.Red("Failed:                %d\n", len(failed))
		color.Red("Failed samples:        %s\n", strings.Join(failed, ", "))
	}
	if len(missingReads) > 0 {
		color.Yellow("Skipped missing reads: %s\n", strings.Join(missingReads, ", "))
	}
	if len(longReadSamples) > 0 {
		color.Yellow("Skipped unsupported:   %s\n", strings.Join(longReadSamples, ", "))
	}
	if len(scanFailures) > 0 {
		color.Red("Scan failures:         %s\n", strings.Join(scanFailures, ", "))
	}
}

func getPathsWithCleanReadsDir(dataDir, species string) ([]string, error) {
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, err
	}

	var roots []string
	for _, match := range matches {
		info, err := os.Stat(match)
		if err != nil || !info.IsDir() {
			continue
		}
		roots = append(roots, match)
	}

	sort.Strings(roots)
	return roots, nil
}

// buildPlan inspects the current SampleBamState and returns the complete
// ordered list of steps still required to finish the sample, plus a bool
// indicating whether paired FASTQ reads must be located before queuing.
// Steps are determined once from the snapshot; processSample executes them
// in order without any further filesystem re-scans between steps.
func buildPlan(state SampleBamState, bqsr, gatkNeedsBam bool) ([]ProcessingReason, bool) {
	var plan []ProcessingReason
	needsReads := false

	rgmdCramValid := isUsable(state.RgmdCram)

	if rgmdCramValid {
		if !hasIndex(state.RgmdCram) {
			plan = append(plan, ReasonIndexRgmdCram)
		}
	} else {
		switch {
		case isUsable(state.RgmdBam):
			// rgmd.bam exists — index it if needed, then convert to cram.
			if !hasIndex(state.RgmdBam) {
				plan = append(plan, ReasonIndexRgmdBam)
			}
			plan = append(plan, ReasonConvertRgmdBam)
		case isUsable(state.SortedBam):
			plan = append(plan, ReasonMarkDupSortedBam, ReasonConvertRgmdBam)
		default:
			plan = append(plan, ReasonAlignFromReads, ReasonConvertRgmdBam)
			needsReads = true
		}
		// After producing the cram we always need to index it.
		plan = append(plan, ReasonIndexRgmdCram)
	}

	if bqsr {
		if !isUsable(state.BqsrCram) {
			if isUsable(state.BqsrBam) {
				plan = append(plan, ReasonConvertBqsrBam)
			} else {
				// Where GATK cannot read an indexed cram, the recalibration
				// steps read the rgmd.bam instead, so the plan has to guarantee
				// an indexed one is in place by the time BQSR starts.
				if gatkNeedsBam {
					switch {
					case containsReason(plan, ReasonAlignFromReads), containsReason(plan, ReasonMarkDupSortedBam):
						// A fresh rgmd.bam is about to be written. Converting it
						// to cram does not remove it — cleanup does, at the end
						// — so all it still needs is an index.
						if !containsReason(plan, ReasonIndexRgmdBam) {
							plan = append(plan, ReasonIndexRgmdBam)
						}
					case isUsable(state.RgmdBam):
						if !hasIndex(state.RgmdBam) && !containsReason(plan, ReasonIndexRgmdBam) {
							plan = append(plan, ReasonIndexRgmdBam)
						}
					default:
						// The cram is the only copy left; decode it back.
						plan = append(plan, ReasonMaterializeRgmdBam)
					}
				}
				plan = append(plan, ReasonRunBQSR, ReasonConvertBqsrBam)
			}
			plan = append(plan, ReasonIndexBqsrCram)
		} else if !hasIndex(state.BqsrCram) {
			plan = append(plan, ReasonIndexBqsrCram)
		}
	}

	plan = append(plan, ReasonCleanup)
	return plan, needsReads
}

// containsReason reports whether a step is already in a plan.
func containsReason(plan []ProcessingReason, want ProcessingReason) bool {
	for _, step := range plan {
		if step == want {
			return true
		}
	}
	return false
}

func processSample(task sampleTask, refFasta, gatkLogLevel, aligner string, verbose, quick, skipVer, bqsr, bootstrap bool, knownSites []string, opts runtimeOpts) sampleResult {
	if err := os.MkdirAll(task.bamDir, 0o755); err != nil {
		return sampleResult{sample: task.sample, err: err}
	}

	sn := task.sample

	// Same pure function on the same reference as the one buildPlan was given,
	// so the steps in the plan and the file this picks below cannot disagree.
	gatkNeedsBam := gatkReadsNeedBam(refFasta)

	// Live paths: updated in-place as each step produces a new file, so later
	// steps always have the correct path without any filesystem re-scan.
	rgmdBamPath := task.state.RgmdBam.Path
	rgmdCramPath := task.state.RgmdCram.Path
	bqsrBamPath := task.state.BqsrBam.Path
	bqsrCramPath := task.state.BqsrCram.Path

	color.Cyan("[%s] Plan (%d steps): %v\n", sn, len(task.plan), task.plan)

	for i, step := range task.plan {
		color.Cyan("[%s] Step %d/%d: %s\n", sn, i+1, len(task.plan), step)

		switch step {

		// ------------------------------------------------------------------ //
		// Align from raw reads → sorted.bam → rgmd.bam → rgmd.cram           //
		// ------------------------------------------------------------------ //
		case ReasonAlignFromReads:
			if task.fwd == "" || task.rev == "" {
				return sampleResult{sample: sn, err: fmt.Errorf("no usable alignment intermediates or paired reads found")}
			}
			if !skipVer {
				color.Cyan("[%s] Validating forward reads: %s\n", sn, task.fwd)
				if err := utils.ValidateFastqGz(task.fwd, verbose, quick); err != nil {
					color.Red("[%s] Forward read validation failed: %v\n", sn, err)
					return sampleResult{sample: sn, err: fmt.Errorf("forward read validation failed: %w", err)}
				}
				color.Green("[%s] Forward reads are valid\n", sn)

				color.Cyan("[%s] Validating reverse reads: %s\n", sn, task.rev)
				if err := utils.ValidateFastqGz(task.rev, verbose, quick); err != nil {
					color.Red("[%s] Reverse read validation failed: %v\n", sn, err)
					return sampleResult{sample: sn, err: fmt.Errorf("reverse read validation failed: %w", err)}
				}
				color.Green("[%s] Reverse reads are valid\n", sn)
			}

			sortedBam := filepath.Join(task.bamDir, sn+".sorted.bam")
			color.Cyan("[%s] Aligning PE reads → %s using %s ...\n", sn, sortedBam, aligner)
			if err := alignPairedReads(task.fwd, task.rev, refFasta, sortedBam, sn, aligner, opts.threads, verbose); err != nil {
				color.Red("[%s] Alignment failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] Alignment done: %s\n", sn, sortedBam)

			rgmdBamPath = alignment.RgmdBamPath(sortedBam)
			color.Cyan("[%s] MarkDuplicates: sorted.bam → rgmd.bam ...\n", sn)
			if err := alignment.MarkDuplicates(refFasta, sortedBam, verbose, aligner, gatkLogLevel, opts.javaOpts); err != nil {
				color.Red("[%s] MarkDuplicates failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] MarkDuplicates done: %s\n", sn, rgmdBamPath)

			// The next step (ReasonConvertRgmdBam) is always in the plan after
			// ReasonAlignFromReads; set rgmdCramPath so it has the right target.
			rgmdCramPath = strings.TrimSuffix(rgmdBamPath, filepath.Ext(rgmdBamPath)) + ".cram"

		// ------------------------------------------------------------------ //
		// sorted.bam → rgmd.bam → rgmd.cram (MarkDuplicates path)            //
		// ------------------------------------------------------------------ //
		case ReasonMarkDupSortedBam:
			rgmdBamPath = alignment.RgmdBamPath(task.state.SortedBam.Path)
			color.Cyan("[%s] MarkDuplicates: sorted.bam → rgmd.bam ...\n", sn)
			if err := alignment.MarkDuplicates(refFasta, task.state.SortedBam.Path, verbose, aligner, gatkLogLevel, opts.javaOpts); err != nil {
				color.Red("[%s] MarkDuplicates failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] MarkDuplicates done: %s\n", sn, rgmdBamPath)
			rgmdCramPath = strings.TrimSuffix(rgmdBamPath, filepath.Ext(rgmdBamPath)) + ".cram"

		// ------------------------------------------------------------------ //
		// Index rgmd.bam (needed before samtools view / BamToCram)           //
		// ------------------------------------------------------------------ //
		case ReasonIndexRgmdBam:
			color.Cyan("[%s] Indexing rgmd.bam: %s ...\n", sn, rgmdBamPath)
			if err := alignment.BamIndex(rgmdBamPath, opts.threads, verbose); err != nil {
				color.Red("[%s] rgmd.bam index failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] rgmd.bam indexed successfully\n", sn)

		// ------------------------------------------------------------------ //
		// rgmd.bam → rgmd.cram                                                //
		// ------------------------------------------------------------------ //
		case ReasonConvertRgmdBam:
			color.Cyan("[%s] Converting rgmd.bam → rgmd.cram ...\n", sn)
			if err := alignment.BamToCram(rgmdBamPath, refFasta, opts.threads, verbose); err != nil {
				color.Red("[%s] BamToCram (rgmd) failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			rgmdCramPath = strings.TrimSuffix(rgmdBamPath, filepath.Ext(rgmdBamPath)) + ".cram"
			color.Green("[%s] rgmd.cram created: %s\n", sn, rgmdCramPath)

		// ------------------------------------------------------------------ //
		// Index rgmd.cram                                                     //
		// ------------------------------------------------------------------ //
		case ReasonIndexRgmdCram:
			color.Cyan("[%s] Indexing rgmd.cram: %s ...\n", sn, rgmdCramPath)
			if err := alignment.BamIndex(rgmdCramPath, opts.threads, verbose); err != nil {
				color.Red("[%s] rgmd.cram index failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] rgmd.cram indexed successfully\n", sn)

		// ------------------------------------------------------------------ //
		// rgmd.cram → rgmd.bam (only where GATK cannot read the cram)        //
		// ------------------------------------------------------------------ //
		case ReasonMaterializeRgmdBam:
			color.Cyan("[%s] Decoding rgmd.cram → rgmd.bam for the GATK steps ...\n", sn)
			decoded, err := alignment.CramToBam(rgmdCramPath, refFasta, opts.threads, verbose)
			if err != nil {
				color.Red("[%s] CramToBam (rgmd) failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			rgmdBamPath = decoded
			color.Cyan("[%s] Indexing rgmd.bam: %s ...\n", sn, rgmdBamPath)
			if err := alignment.BamIndex(rgmdBamPath, opts.threads, verbose); err != nil {
				color.Red("[%s] rgmd.bam index failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] rgmd.bam ready: %s\n", sn, rgmdBamPath)

		// ------------------------------------------------------------------ //
		// BQSR: rgmd alignment → bqsr.bam                                    //
		// ------------------------------------------------------------------ //
		case ReasonRunBQSR:
			color.Cyan("[%s] BQSR route: %s\n", sn, describeRoute(refFasta))

			// The cram is the input everywhere it can be read; where it cannot,
			// buildPlan has already put the steps that produce the bam in front
			// of this one.
			bqsrInput := rgmdCramPath
			if gatkNeedsBam {
				if rgmdBamPath == "" {
					return sampleResult{sample: sn, err: fmt.Errorf("BQSR needs an rgmd.bam on this reference but none was produced")}
				}
				bqsrInput = rgmdBamPath
				color.Cyan("[%s] Reading %s — a .crai cannot be read as a BAI on this reference\n", sn, filepath.Base(bqsrInput))
			}

			sampleKnownSites := append([]string(nil), knownSites...)
			if len(sampleKnownSites) == 0 && bootstrap {
				color.Cyan("[%s] Bootstrapping BQSR known sites from %s ...\n", sn, bqsrInput)
				var err error
				sampleKnownSites, err = createBootstrapKnownSites(refFasta, bqsrInput, task.bamDir, gatkLogLevel, opts, verbose)
				if err != nil {
					color.Red("[%s] Bootstrapping BQSR known sites failed: %v\n", sn, err)
					return sampleResult{sample: sn, err: err}
				}
				color.Green("[%s] Bootstrap known sites created\n", sn)
			}
			color.Cyan("[%s] Running BQSR on %s ...\n", sn, bqsrInput)
			outBam, err := runBQSR(refFasta, bqsrInput, task.bamDir, sampleKnownSites, opts, verbose)
			if err != nil {
				color.Red("[%s] BQSR failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			bqsrBamPath = outBam
			color.Green("[%s] BQSR done: %s\n", sn, bqsrBamPath)

		// ------------------------------------------------------------------ //
		// bqsr.bam → bqsr.cram                                               //
		// ------------------------------------------------------------------ //
		case ReasonConvertBqsrBam:
			color.Cyan("[%s] Converting bqsr.bam → bqsr.cram ...\n", sn)
			if err := alignment.BamToCram(bqsrBamPath, refFasta, opts.threads, verbose); err != nil {
				color.Red("[%s] BamToCram (bqsr) failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			bqsrCramPath = strings.TrimSuffix(bqsrBamPath, filepath.Ext(bqsrBamPath)) + ".cram"
			color.Green("[%s] bqsr.cram created: %s\n", sn, bqsrCramPath)

		// ------------------------------------------------------------------ //
		// Index bqsr.cram                                                     //
		// ------------------------------------------------------------------ //
		case ReasonIndexBqsrCram:
			if bqsrCramPath == "" {
				bqsrCramPath = strings.TrimSuffix(bqsrBamPath, filepath.Ext(bqsrBamPath)) + ".cram"
			}
			color.Cyan("[%s] Indexing bqsr.cram: %s ...\n", sn, bqsrCramPath)
			if err := alignment.BamIndex(bqsrCramPath, opts.threads, verbose); err != nil {
				color.Red("[%s] bqsr.cram index failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] bqsr.cram indexed successfully\n", sn)

		// ------------------------------------------------------------------ //
		// Cleanup — remove all intermediates; keep final CRAMs + indexes      //
		// ------------------------------------------------------------------ //
		case ReasonCleanup:
			color.Cyan("[%s] Cleaning up intermediate files in %s ...\n", sn, task.bamDir)
			if err := cleanupSampleOutputs(task.bamDir, rgmdCramPath, bqsrCramPath, bqsr); err != nil {
				color.Red("[%s] Cleanup failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] Cleanup done\n", sn)

		default:
			return sampleResult{sample: sn, err: fmt.Errorf("unknown processing step: %q", step)}
		}
	}

	color.Green("[%s] All steps complete ✅\n", sn)
	return sampleResult{sample: sn, success: true}
}

func inspectSampleBamDir(sample, bamDir, refFasta string, verbose, quick bool) (SampleBamState, error) {
	state := SampleBamState{Sample: sample}

	info, err := os.Stat(bamDir)
	if err != nil {
		if os.IsNotExist(err) {
			return state, nil
		}
		return state, err
	}
	if !info.IsDir() {
		return state, fmt.Errorf("%s is not a directory", bamDir)
	}

	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return state, err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}

		name := entry.Name()
		lowerName := strings.ToLower(name)
		fullPath := filepath.Join(bamDir, name)

		switch {
		case strings.HasSuffix(lowerName, "sorted.bam"):
			state.SortedBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "rgmd.bam"):
			state.RgmdBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "rgmd.cram"):
			state.RgmdCram = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "bqsr.bam"):
			state.BqsrBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "bqsr.cram"):
			state.BqsrCram = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, ".bai"), strings.HasSuffix(lowerName, ".csi"), strings.HasSuffix(lowerName, ".crai"),
			strings.HasSuffix(lowerName, ".pdf"), strings.HasSuffix(lowerName, ".txt"), strings.HasSuffix(lowerName, ".list"):
			continue
		default:
			fileInfo, err := entry.Info()
			size := int64(0)
			if err == nil {
				size = fileInfo.Size()
			}
			state.OtherFiles = append(state.OtherFiles, FileInfo{
				Path:    fullPath,
				Size:    size,
				Present: true,
			})
		}
	}

	return state, nil
}

func getAlignmentFileInfo(path, refFasta string, verbose, quick bool) FileInfo {
	//color.Cyan("[%s] Checking %s\n", path, path)
	info := FileInfo{
		Path:    path,
		Present: true,
	}

	stat, err := os.Stat(path)
	if err != nil || !stat.Mode().IsRegular() {
		return info
	}

	info.Size = stat.Size()
	info.ValidateErr = utils.ValidateBam(path, refFasta, verbose, quick)
	if info.ValidateErr == nil {
		// A file can be perfectly well formed and still not belong to this
		// reference. quickcheck does not decode, so it passes such a file
		// straight through to GATK, which fails much later and far less
		// clearly. See checkAlignmentContigs.
		if err := checkAlignmentContigs(path, refFasta); err != nil {
			color.Red("%s: %v\n", filepath.Base(path), err)
			info.ValidateErr = err
		}
	}
	info.Valid = info.ValidateErr == nil

	lower := strings.ToLower(path)
	if strings.HasSuffix(lower, ".bam") || strings.HasSuffix(lower, ".cram") {
		// findIndex reads the index rather than only stat-ing it, so an index
		// left truncated by an interrupted run is reported as missing and gets
		// rebuilt instead of being handed to GATK.
		_, idxSize, found := findIndex(path, verbose)
		info.IndexPresent = found
		info.IndexSize = idxSize
	} else {
		info.IndexPresent = false
		info.IndexSize = 0
	}

	return info
}

func alignPairedReads(fwd, rev, refFasta, sortedBam, sample, aligner string, threads int, verbose bool) error {
	if err := removeIfExists(sortedBam); err != nil {
		return err
	}

	// Inline removal of index files for sortedBam
	var idxCandidates []string
	if strings.HasSuffix(strings.ToLower(sortedBam), ".bam") {
		idxCandidates = []string{strings.TrimSuffix(sortedBam, filepath.Ext(sortedBam)) + ".bai"}
	} else if strings.HasSuffix(strings.ToLower(sortedBam), ".cram") {
		stem := strings.TrimSuffix(sortedBam, filepath.Ext(sortedBam))
		idxCandidates = []string{stem + ".crai", sortedBam + ".crai"}
	}
	if err := removeIfExists(idxCandidates...); err != nil {
		return err
	}

	lib := fmt.Sprintf("%s_1", sample)
	var err error
	switch aligner {
	case "bwa-mem2":
		err = alignment.BwaMem2Align(fwd, rev, refFasta, sample, lib, threads, sortedBam, verbose)
	case "bwa-mem":
		err = alignment.BwaMemAlign(fwd, rev, refFasta, sample, lib, threads, sortedBam, verbose)
	case "bowtie2":
		err = alignment.Bowtie2Align(fwd, rev, refFasta, sortedBam, sample, lib, threads, verbose)
	default:
		return fmt.Errorf("aligner %s is not supported for paired-read directory alignment", aligner)
	}
	return err
}

func runBQSR(refFasta, input, bamDir string, knownSites []string, opts runtimeOpts, verbose bool) (string, error) {
	var output string
	switch {
	case strings.HasSuffix(strings.ToLower(input), ".bam"):
		output = strings.TrimSuffix(input, filepath.Ext(input)) + "_bqsr.bam"
	case strings.HasSuffix(strings.ToLower(input), ".cram"):
		output = strings.TrimSuffix(input, filepath.Ext(input)) + "_bqsr.bam"
	default:
		return "", fmt.Errorf("input file %s is not a BAM or CRAM file", input)
	}
	base := strings.TrimSuffix(input, filepath.Ext(input))
	sn := filepath.Base(base)
	recalTable := base + "recal_table.txt"
	recalTable2 := base + "recal_table2.txt"
	plots := base + "recal_table_plots.pdf"
	shardDir := filepath.Join(bamDir, "shards")

	if err := removeIfExists(output, recalTable, recalTable2, plots); err != nil {
		return "", err
	}
	if err := removeIfExists(indexCandidates(output)...); err != nil {
		return "", err
	}

	var args []string
	for _, site := range knownSites {
		args = append(args, "--known-sites "+shQuote(site))
	}
	knownSitesArgs := strings.Join(args, " ")

	// Each recalibration pass is one process over the whole reference, unless the
	// reference is large enough that useScatter sends it to the interval scatter
	// in scatter.go. Everything else in this function is the same either way.
	scatter := useScatter(refFasta)

	recalibrate := func(bam, table, label string) error {
		if scatter {
			return scatterBaseRecalibrator(refFasta, bam, table, knownSitesArgs, shardDir, label, opts, verbose)
		}
		cmd := fmt.Sprintf(`gatk --java-options "%s" BaseRecalibrator -R %s -I %s %s -O %s --maximum-cycle-value %d --tmp-dir %s`,
			opts.javaOpts, shQuote(refFasta), shQuote(bam), knownSitesArgs, shQuote(table), alignment.MaxCycleValue, shQuote(alignment.WorkTmpDir(table)))
		return runGatk(cmd, verbose)
	}

	if err := recalibrate(input, recalTable, "before"); err != nil {
		color.Red("[%s] BaseRecalibrator failed: %v\n", sn, err)
		return "", err
	}
	color.Green("[%s] BaseRecalibrator completed successfully\n", sn)

	// --create-output-bam-index false is not optional on a large genome: htsjdk
	// can only write a BAI, and a BAI cannot address positions past 2^29-1, so
	// ApplyBQSR would do the whole pass and then die writing the index. The BAM
	// is indexed as CSI by BamIndex below instead.
	cmd2 := fmt.Sprintf(`gatk --java-options "%s" ApplyBQSR -R %s -I %s -bqsr %s -O %s --create-output-bam-index false --tmp-dir %s`,
		opts.javaOpts, shQuote(refFasta), shQuote(input), shQuote(recalTable), shQuote(output), shQuote(alignment.WorkTmpDir(output)))
	if err := runGatk(cmd2, verbose); err != nil {
		color.Red("[%s] ApplyBQSR failed: %v\n", sn, err)
		return "", err
	}
	color.Green("[%s] ApplyBQSR completed successfully\n", sn)

	// ApplyBQSR was told not to write an index, so one is built here. The
	// scattered second pass additionally needs it, to query intervals.
	if err := alignment.BamIndex(output, opts.threads, verbose); err != nil {
		color.Red("[%s] BamIndex failed: %v\n", sn, err)
		return "", err
	}

	if err := recalibrate(output, recalTable2, "after"); err != nil {
		color.Red("[%s] Second BaseRecalibrator failed: %v\n", sn, err)
		return "", err
	}
	color.Green("[%s] Second BaseRecalibrator completed successfully\n", sn)

	// AnalyzeCovariates only draws the before/after QC plot, and it draws it
	// through Rscript, so it fails on a host where R or one of BQSR.R's packages
	// (gsalib, gplots) is missing. That is not a reason to throw away a
	// recalibration that has already finished, so it is reported and stepped
	// over. Both recalibration tables are kept either way.
	cmd4 := fmt.Sprintf(`gatk --java-options "%s" AnalyzeCovariates -before %s -after %s -plots %s --tmp-dir %s`,
		opts.javaOpts, shQuote(recalTable), shQuote(recalTable2), shQuote(plots), shQuote(alignment.WorkTmpDir(plots)))
	if err := runGatk(cmd4, verbose); err != nil {
		color.Yellow("[%s] AnalyzeCovariates failed: %v\n", sn, err)
		color.Yellow("[%s] Recalibration is complete; only the QC plot is missing. Check Rscript and its gsalib/gplots packages.\n", sn)
	} else {
		color.Green("[%s] AnalyzeCovariates completed successfully\n", sn)
	}

	return output, nil
}

func createBootstrapKnownSites(refFasta, input, bamDir, gatkLogLevel string, opts runtimeOpts, verbose bool) ([]string, error) {
	base := strings.TrimSuffix(input, filepath.Ext(input))
	// Bgzip everywhere it works, uncompressed where a .tbi cannot address the
	// contigs — the two filtered files at the end are read back by
	// BaseRecalibrator, which needs an index it can query. See variantSuffixes.
	vcfExt, idxExt := variantSuffixes(refFasta)
	rawVCF := base + ".raw" + vcfExt
	snpVCF := base + ".raw.SNP" + vcfExt
	indelVCF := base + ".raw.INDEL" + vcfExt
	snpColumns := base + ".raw.SNP.columns" + vcfExt
	indelColumns := base + ".raw.INDEL.columns" + vcfExt
	filteredSNP := base + ".raw.SNP.hard_filtered" + vcfExt
	filteredINDEL := base + ".raw.INDEL.hard_filtered" + vcfExt
	shardDir := filepath.Join(bamDir, "shards")

	// Only the calling step has two implementations, so it carries a func while
	// the filtering steps that follow it stay plain commands. Those read one VCF
	// and write another and are not worth scattering at any genome size.
	type bootstrapStep struct {
		name   string
		output string
		cmd    string
		run    func() error
	}

	// One HaplotypeCaller over the whole reference, unless the reference is large
	// enough that useScatter sends it to the interval scatter in scatter.go. The
	// pair-HMM thread count is left at GATK's own default here: this process has
	// the sample to itself, which is not true of a shard.
	wholeGenomeCall := fmt.Sprintf(`gatk --java-options "%s" HaplotypeCaller -R %s -I %s -O %s --tmp-dir %s`,
		opts.javaOpts, shQuote(refFasta), shQuote(input), shQuote(rawVCF), shQuote(alignment.WorkTmpDir(rawVCF)))

	steps := []bootstrapStep{
		{
			name:   "HaplotypeCaller",
			output: rawVCF,
			run: func() error {
				if useScatter(refFasta) {
					return scatterHaplotypeCaller(refFasta, input, rawVCF, shardDir, opts, verbose)
				}
				return runGatk(wholeGenomeCall, verbose)
			},
		},
		{
			name:   "SelectVariants (SNP)",
			output: snpVCF,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" SelectVariants -V %s --select-type-to-include SNP -O %s --verbosity %s`, opts.javaOpts, shQuote(rawVCF), shQuote(snpVCF), gatkLogLevel),
		},
		{
			name:   "SelectVariants (INDEL)",
			output: indelVCF,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" SelectVariants -V %s --select-type-to-include INDEL -O %s --verbosity %s`, opts.javaOpts, shQuote(rawVCF), shQuote(indelVCF), gatkLogLevel),
		},
		{
			name:   "VariantFiltration (SNP)",
			output: snpColumns,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O %s --verbosity %s`, opts.javaOpts, shQuote(snpVCF), shQuote(snpColumns), gatkLogLevel),
		},
		{
			name:   "SelectVariants (Filtered SNP)",
			output: filteredSNP,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" SelectVariants --exclude-filtered -V %s -O %s --verbosity %s`, opts.javaOpts, shQuote(snpColumns), shQuote(filteredSNP), gatkLogLevel),
		},
		{
			name:   "VariantFiltration (INDEL)",
			output: indelColumns,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O %s --verbosity %s`, opts.javaOpts, shQuote(indelVCF), shQuote(indelColumns), gatkLogLevel),
		},
		{
			name:   "SelectVariants (Filtered INDEL)",
			output: filteredINDEL,
			cmd:    fmt.Sprintf(`gatk --java-options "%s" SelectVariants --exclude-filtered -V %s -O %s --verbosity %s`, opts.javaOpts, shQuote(indelColumns), shQuote(filteredINDEL), gatkLogLevel),
		},
	}

	for _, step := range steps {
		if info, err := os.Stat(step.output); err == nil && info.Size() > 0 {
			if vErr := utils.ValidateGvcf(step.output, verbose, true); vErr == nil {
				color.Green("[%s] Output %s already exists and is valid. Skipping...", step.name, step.output)
				continue
			}
			color.Yellow("[%s] Output %s exists but is invalid. Re-creating...", step.name, step.output)
		}

		_ = removeIfExists(step.output)
		_ = removeIfExists(step.output + idxExt)

		color.Cyan("[%s] Starting ...\n", step.name)
		var err error
		if step.run != nil {
			err = step.run()
		} else {
			err = runGatk(step.cmd, verbose)
		}
		if err != nil {
			return nil, fmt.Errorf("step %s failed: %w", step.name, err)
		}
	}

	return []string{filteredSNP, filteredINDEL}, nil
}

func createBootstrapKnownSitesOld(refFasta, input, gatkLogLevel string, verbose bool) ([]string, error) {
	base := strings.TrimSuffix(input, filepath.Ext(input))
	rawVCF := base + ".raw.vcf.gz"
	snpVCF := base + ".raw.SNP.vcf.gz"
	indelVCF := base + ".raw.INDEL.vcf.gz"
	snpColumns := base + ".raw.SNP.columns.vcf.gz"
	indelColumns := base + ".raw.INDEL.columns.vcf.gz"
	filteredSNP := base + ".raw.SNP.hard_filtered.vcf.gz"
	filteredINDEL := base + ".raw.INDEL.hard_filtered.vcf.gz"

	if err := removeIfExists(rawVCF, snpVCF, indelVCF, snpColumns, indelColumns, filteredSNP, filteredINDEL); err != nil {
		return nil, err
	}

	cmds := []string{
		fmt.Sprintf("gatk HaplotypeCaller -R %s -I %s -O %s", shQuote(refFasta), shQuote(input), shQuote(rawVCF)),
		fmt.Sprintf("gatk SelectVariants -V %s --select-type-to-include SNP -O %s --verbosity %s", shQuote(rawVCF), shQuote(snpVCF), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants -V %s --select-type-to-include INDEL -O %s --verbosity %s", shQuote(rawVCF), shQuote(indelVCF), gatkLogLevel),
		fmt.Sprintf(`gatk VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O %s --verbosity %s`, shQuote(snpVCF), shQuote(snpColumns), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants --exclude-filtered -V %s -O %s --verbosity %s", shQuote(snpColumns), shQuote(filteredSNP), gatkLogLevel),
		fmt.Sprintf(`gatk VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O %s --verbosity %s`, shQuote(indelVCF), shQuote(indelColumns), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants --exclude-filtered -V %s -O %s --verbosity %s", shQuote(indelColumns), shQuote(filteredINDEL), gatkLogLevel),
	}

	for _, cmd := range cmds {
		fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", cmd)
		if err := runBash(cmd, verbose); err != nil {
			return nil, err
		}
	}

	return []string{filteredSNP, filteredINDEL}, nil
}

// cleanupSampleOutputs removes all files in bamDir except:
//   - rgmd.cram and its index
//   - bqsr.cram and its index when present
//   - any .txt or .pdf files (logs, recal tables, plots)
//
// rgmdCramPath is the live path produced by the pipeline.
// bqsrCramPath is the final bqsr.cram path when available.
func cleanupSampleOutputs(bamDir, rgmdCramPath, bqsrCramPath string, requireBQSR bool) error {
	keep := make(map[string]struct{})

	addWithIndex := func(cramPath string) {
		keep[cramPath] = struct{}{}
		for _, idx := range indexCandidates(cramPath) {
			if _, err := os.Stat(idx); err == nil {
				keep[idx] = struct{}{}
			}
		}
	}

	if rgmdCramPath != "" {
		addWithIndex(rgmdCramPath)
	}
	if bqsrCramPath != "" {
		addWithIndex(bqsrCramPath)
	}

	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		fullPath := filepath.Join(bamDir, entry.Name())
		if _, ok := keep[fullPath]; ok {
			continue
		}
		// Always preserve log/report files.
		lower := strings.ToLower(entry.Name())
		if strings.HasSuffix(lower, ".txt") || strings.HasSuffix(lower, ".pdf") {
			continue
		}
		if err := os.Remove(fullPath); err != nil && !os.IsNotExist(err) {
			return fmt.Errorf("removing %s: %w", fullPath, err)
		}
	}

	// Scratch directories are skipped by the loop above because they are
	// directories, and the shard VCFs inside them are as large as the merged
	// call set, so they are removed explicitly once the final CRAMs exist.
	for _, dir := range scratchDirs(bamDir) {
		if err := os.RemoveAll(dir); err != nil && !os.IsNotExist(err) {
			return fmt.Errorf("removing %s: %w", dir, err)
		}
	}

	return nil
}

// Helper function to remove files if they exist.
// Small wrapper around os.Remove; could be inlined.
func removeIfExists(paths ...string) error {
	for _, path := range paths {
		if path == "" {
			continue
		}
		if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
			return err
		}
	}
	return nil
}

// Helper function to execute a bash command.
// Small wrapper around exec.Command; could be inlined.
func runBash(cmdStr string, verbose bool) error {
	cmd := exec.Command("bash", "-c", cmdStr)
	if verbose {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
	}
	return cmd.Run()
}

// Helper function to shell-quote a string.
// Small string manipulation; could be inlined.
func shQuote(value string) string {
	if value == "" {
		return "''"
	}
	return "'" + strings.ReplaceAll(value, `'`, `'\''`) + "'"
}
