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
	// ReasonAlignLongReads — no usable intermediates; align the sample's single
	// long-read FASTQ with pbmm2, then mark duplicates. The long-read twin of
	// ReasonAlignFromReads, kept as its own step so the plan printed for an
	// operator says which aligner a sample is about to get.
	ReasonAlignLongReads ProcessingReason = "align_long_reads"
	// ReasonAdoptExternalAlignment — a long-read sample arrived with an
	// externally produced alignment (MENINA_LONG_READS.aligned.cram and the
	// like). Rewrite its read group to this sample and normalise it to
	// <sample>.sorted.bam so the ordinary markdup → cram chain applies. The
	// external file is an input: never modified, never removed.
	ReasonAdoptExternalAlignment ProcessingReason = "adopt_external_alignment"
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
	// se is the single long-read FASTQ; empty for a short-read sample.
	se string
	// longRead marks a pbmm2 sample: no BQSR, ever.
	longRead bool
	// aligner is the aligner THIS sample's steps use — the run's --aligner for a
	// short-read sample, always "pbmm2" for a long-read one.
	//
	// It lives on the task rather than being passed to processSample because a
	// mixed run has two of them: a function holding both a parameter and a task
	// field is one typo away from running gatk MarkDuplicates on a PacBio BAM.
	aligner string
	// preset and refIndex are the pbmm2 alignment mode and the reference to
	// hand it (a preset-keyed .mmi, or the fasta where one could not be built).
	preset   string
	refIndex string
	// dupMarker is the duplicate marker for this sample, one of the
	// alignment.DupMarker* constants. Separate from aligner because the two are
	// different questions: every sample here is marked with GATK, including the
	// long-read ones, because pbmarkdup needs the native PacBio per-read tags
	// (zm, np, rq) that a FASTQ cannot carry and a pbmm2-from-FASTQ BAM
	// therefore lacks.
	dupMarker string
	// external is an externally produced alignment adopted as this sample's
	// input. It is never written into a SampleBamState slot: those five slots
	// describe artefacts the pipeline produced, and this is not one.
	external string
	state    SampleBamState
	plan     []ProcessingReason
}

type sampleResult struct {
	sample  string
	success bool
	err     error
}

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, preset string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
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

	if dicfErr := utils.EnsureGatkDict(resolvedFasta); dicfErr != nil {
		fmt.Printf("%v\n", dicfErr)
		return
	}

	// --aligner names the SHORT-read aligner here. Long-read samples are
	// detected by name and always use pbmm2, so pbmm2 is not a choice to make
	// for the run: alignPairedReads has no arm for it and would fail every
	// short-read sample individually, after the scan, hours in.
	if aligner == "pbmm2" {
		color.Red("--aligner pbmm2 is not valid in data-dir mode.\n")
		color.Red("Short-read samples need bwa-mem, bwa-mem2 or bowtie2; long-read samples are\n")
		color.Red("aligned with pbmm2 automatically, whatever --aligner says.\n")
		return
	}

	if presetErr := validatePreset(preset); presetErr != nil {
		color.Red("%v\n", presetErr)
		return
	}

	// BQSR is applied per sample, not per run: short-read samples are
	// recalibrated and long-read ones never are (see planInputs.longRead), so
	// the two kinds coexist under one --bqsr.
	if bqsr {
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
		tasks             []sampleTask
		alreadyComplete   int
		scanFailures      []string
		missingReads      []string
		missingLongReads  []string
		longReadSamples   []string
		ambiguousAdoption []string
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

		isLR := IsLongReadSample(sample)

		// Adoption is decided here, from names alone: this loop is serial, so
		// the header reads and the full validation belong in the worker, inside
		// ReasonAdoptExternalAlignment.
		var external string
		if isLR {
			candidates := adoptionCandidates(state)
			if len(candidates) > 1 {
				// Two movies, or a movie and a stale partial. Picking one
				// silently drops half a sample's data, so this is a question
				// for a person; nothing in the directory is touched.
				color.Red("[%s] More than one externally produced alignment in %s — not choosing between them:\n", sample, bamDir)
				for _, c := range candidates {
					color.Red("    %s\n", filepath.Base(c))
				}
				ambiguousAdoption = append(ambiguousAdoption, sample)
				continue
			}
			if len(candidates) == 1 {
				external = candidates[0]
				color.Blue("[%s] Long read sample with an existing alignment: %s\n", sample, filepath.Base(external))
			} else {
				color.Blue("[%s] Long read sample\n", sample)
			}
		}

		rgmdCramOK := state.RgmdCram.Present && state.RgmdCram.Valid && state.RgmdCram.IndexPresent && state.RgmdCram.IndexSize > 0
		bqsrCramOK := state.BqsrCram.Present && state.BqsrCram.Valid && state.BqsrCram.IndexPresent && state.BqsrCram.IndexSize > 0

		if isLR {
			// A long-read sample is finished at the rgmd.cram: BQSR is never
			// run on it, so requiring a bqsr.cram would mean no long-read
			// sample was ever complete.
			//
			// OtherFiles is not required to be empty either. An adopted
			// alignment stays in the directory permanently — it is an input,
			// and cleanup will not remove it — so demanding its absence would
			// re-plan and re-validate every adopted sample on every run, which
			// on a multi-hundred-gigabyte CRAM is not free.
			if rgmdCramOK && !state.SortedBam.Present && !state.RgmdBam.Present && !state.BqsrBam.Present {
				alreadyComplete++
				color.Green("[%s]- ✅ PASS (long read). Skipping ...\n", sample)
				continue
			}
		} else if rgmdCramOK && bqsrCramOK && !state.SortedBam.Present && !state.RgmdBam.Present && !state.BqsrBam.Present && len(state.OtherFiles) == 0 {
			alreadyComplete++
			color.Green("[%s]- ✅ PASS. Skipping ...\n", sample)
			continue
		}

		// Build the full ordered plan from the current state so processSample
		// can execute every step without re-scanning the filesystem.
		plan, needsReads := buildPlan(state, planInputs{
			bqsr:         bqsr && !isLR,
			gatkNeedsBam: gatkNeedsBam,
			longRead:     isLR,
			external:     external != "",
		})

		var fwd, rev, se string
		if needsReads {
			if isLR {
				longRead, readsErr := GetReadsLong(cleanReadsDir)
				if readsErr != nil {
					color.Red("[%s] %v\n", sample, readsErr)
					missingLongReads = append(missingLongReads, sample)
					continue
				}
				se = longRead
			} else {
				fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
				if readsErr != nil || len(fwdReads) != 1 || len(revReads) != 1 {
					color.Red("[%s] Forward and reverse reads not found in %s\n", sample, cleanReadsDir)
					missingReads = append(missingReads, sample)
					continue
				}
				fwd = fwdReads[0]
				rev = revReads[0]
			}
		}

		sampleAligner := aligner
		if isLR {
			sampleAligner = "pbmm2"
			longReadSamples = append(longReadSamples, sample)
		}

		color.Yellow("[%s] Sample queued (%s) — %d steps: %v\n\n", sample, sampleAligner, len(plan), plan)
		tasks = append(tasks, sampleTask{
			sample:        sample,
			bamDir:        bamDir,
			cleanReadsDir: cleanReadsDir,
			fwd:           fwd,
			rev:           rev,
			se:            se,
			longRead:      isLR,
			aligner:       sampleAligner,
			dupMarker:     alignment.DupMarkerGatk,
			preset:        preset,
			external:      external,
			state:         state,
			plan:          plan,
		})
	}

	color.Green("Samples complete: %d\n", alreadyComplete)
	color.Green("Samples queued:   %d (%d short-read, %d long-read)\n",
		len(tasks), len(tasks)-len(longReadSamples), len(longReadSamples))
	fmt.Printf("\n-------------------------------------- Queued Samples --------------------------------------\n")
	for _, task := range tasks {
		color.Yellow("%s (%s)\n", task.sample, task.aligner)
		fmt.Printf("\nTo DO:\n-------------------------------------------------\n\n")
		for _, reason := range task.plan {
			fmt.Printf("%s\n", reason)
		}
		if task.external != "" {
			// Named because it changes what the sample costs: the adopted file,
			// the bam it is normalised into, the rgmd.bam and the rgmd.cram are
			// all on disk at once until cleanup runs at the end.
			fmt.Printf("\nadopting: %s (kept — this pipeline did not create it)\n", task.external)
		}
		fmt.Printf("\n\n================================================================\n\n")
	}
	fmt.Printf("----------------------------------------------------------------------------------------------\n")

	if len(scanFailures) > 0 {
		color.Red("Scan failures:    %d\n", len(scanFailures))
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

	// pbmm2's tools are only required once the scan says the tree actually
	// holds long reads, which is why this is checked here rather than from the
	// --aligner flag: a mixed run asks for bwa-mem2 and still needs them.
	needsPbmm2 := false
	for _, task := range tasks {
		if containsReason(task.plan, ReasonAlignLongReads) {
			needsPbmm2 = true
			break
		}
	}
	if needsPbmm2 {
		// pbmm2 only: duplicate marking here is GATK's, including for long-read
		// samples. pbmarkdup cannot write a BAM aligned from a FASTQ — it needs
		// the native PacBio per-read tags a FASTQ cannot carry.
		if depErr := utils.CheckDeps([]string{"pbmm2"}); depErr != nil {
			color.Red("Long-read samples are queued but their tools are missing: %v\n", depErr)
			return
		}
		// Built once, here, while the run is still single-threaded. Building it
		// lazily inside the aligner meant every concurrent long-read sample
		// racing to write the same file.
		refIndex, idxErr := alignment.EnsurePbmm2Index(resolvedFasta, preset, verbose)
		if idxErr != nil {
			color.Red("Could not prepare the pbmm2 reference index: %v\n", idxErr)
			return
		}
		for i := range tasks {
			if tasks[i].longRead {
				tasks[i].refIndex = refIndex
			}
		}
	}

	results := make([]sampleResult, len(tasks))
	opts := newRuntimeOpts(threads, maxParallelJobs)

	color.Cyan("Processing %d samples using %d threads each (Max parallel jobs: %d)\n", len(tasks), threads, maxParallelJobs)
	color.Cyan("GATK: %d JVM slots, %s per sample step, %s per shard, %d interval shards, %d pair-HMM threads\n\n",
		cap(gatkSlots), opts.javaOpts, opts.shardJava, opts.shardCount, opts.pairHmm)

	// Short-read and long-read samples run in separate waves rather than
	// together. maxParallelJobs is derived from core count alone, and
	// newRuntimeOpts sizes the JVM heaps from the memory available *before* any
	// aligner index is resident. Overlapping the two kinds would hold a
	// bwa-mem2 index and a pbmm2 index in memory at the same time, several
	// copies each, against heaps sized as though neither existed. Long-read
	// samples are few in a normal estate, so the wall-clock cost is small.
	runWave := func(label string, want func(sampleTask) bool) {
		var wg sync.WaitGroup
		sem := make(chan struct{}, maxParallelJobs)
		queued := 0
		for i, task := range tasks {
			if !want(task) {
				continue
			}
			queued++
			wg.Add(1)
			go func(idx int, task sampleTask) {
				defer wg.Done()
				sem <- struct{}{}
				defer func() { <-sem }()

				fmt.Printf("\nProcessing sample: %s ....\n\n", task.sample)
				results[idx] = processSample(task, resolvedFasta, gatkLogLevel, verbose, quick, skipVer, bootstrap, knownSites, opts)
				if results[idx].success {
					color.Green("[%s] done\n", task.sample)
					return
				}
				color.Red("[%s] failed: %v\n", task.sample, results[idx].err)
			}(i, task)
		}
		if queued > 0 {
			color.Cyan("\n--------------------- %s wave: %d sample(s) ---------------------\n\n", label, queued)
		}
		wg.Wait()
	}

	runWave("Short read", func(t sampleTask) bool { return !t.longRead })
	runWave("Long read", func(t sampleTask) bool { return t.longRead })

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
		color.Yellow("Skipped missing reads:      %s\n", strings.Join(missingReads, ", "))
	}
	if len(missingLongReads) > 0 {
		color.Yellow("Skipped missing long reads: %s\n", strings.Join(missingLongReads, ", "))
	}
	if len(ambiguousAdoption) > 0 {
		color.Yellow("Skipped, more than one existing alignment: %s\n", strings.Join(ambiguousAdoption, ", "))
	}
	if len(scanFailures) > 0 {
		color.Red("Scan failures:              %s\n", strings.Join(scanFailures, ", "))
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

// planInputs is everything buildPlan decides from beyond the file state.
//
// A struct rather than a tail of bools: buildPlan(state, true, false) was
// already hard to read at its call sites, and long-read support would have made
// it four of them.
type planInputs struct {
	bqsr         bool
	gatkNeedsBam bool
	// longRead marks a pbmm2 sample. BQSR is never applied to one, so the
	// caller passes bqsr && !longRead and this only has to pick the aligner
	// step — but it is carried explicitly so the rule is stated in one place.
	longRead bool
	// external says an adoptable externally produced alignment was found, so
	// the sample needs normalising rather than aligning.
	external bool
}

// buildPlan inspects the current SampleBamState and returns the complete
// ordered list of steps still required to finish the sample, plus a bool
// indicating whether input FASTQs must be located before queuing.
// Steps are determined once from the snapshot; processSample executes them
// in order without any further filesystem re-scans between steps.
func buildPlan(state SampleBamState, in planInputs) ([]ProcessingReason, bool) {
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
			// Also the resume path for an adopted sample: <sample>.sorted.bam is
			// what adoption wrote, so a run that died after it picks up here
			// rather than normalising a hundred-gigabyte cram a second time.
			plan = append(plan, ReasonMarkDupSortedBam, ReasonConvertRgmdBam)
		case in.external:
			// Ranked below a usable sorted.bam and above aligning, and reached
			// whether or not a FASTQ exists — an adopted sample often has no
			// reads on disk at all. An unusable sorted.bam lands here too:
			// falling through to the aligner would be certain failure.
			plan = append(plan, ReasonAdoptExternalAlignment, ReasonMarkDupSortedBam, ReasonConvertRgmdBam)
		case in.longRead:
			plan = append(plan, ReasonAlignLongReads, ReasonConvertRgmdBam)
			needsReads = true
		default:
			plan = append(plan, ReasonAlignFromReads, ReasonConvertRgmdBam)
			needsReads = true
		}
		// After producing the cram we always need to index it.
		plan = append(plan, ReasonIndexRgmdCram)
	}

	// Never on a long-read sample: GATK's recalibration model is built for
	// short-read cycle/context covariates, and a PacBio sample carrying both a
	// .RGMD.cram and a .RGMD_bqsr.cram would additionally make
	// variants.selectAlignments ambiguous — both names contain "rgmd" — which
	// drops the sample from the call set without saying so.
	if in.bqsr && !in.longRead {
		if !isUsable(state.BqsrCram) {
			if isUsable(state.BqsrBam) {
				plan = append(plan, ReasonConvertBqsrBam)
			} else {
				// Where GATK cannot read an indexed cram, the recalibration
				// steps read the rgmd.bam instead, so the plan has to guarantee
				// an indexed one is in place by the time BQSR starts.
				if in.gatkNeedsBam {
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

// processSample executes one sample's plan.
//
// It takes no aligner argument on purpose. A mixed run has two of them, and the
// one a step must use is a property of the sample, not of the run — so it is
// read from task.aligner and nowhere else. A function holding both would be one
// typo away from running gatk MarkDuplicates on a PacBio BAM.
func processSample(task sampleTask, refFasta, gatkLogLevel string, verbose, quick, skipVer, bootstrap bool, knownSites []string, opts runtimeOpts) sampleResult {
	if err := os.MkdirAll(task.bamDir, 0o755); err != nil {
		return sampleResult{sample: task.sample, err: err}
	}

	sn := task.sample

	// Same pure function on the same reference as the one buildPlan was given,
	// so the steps in the plan and the file this picks below cannot disagree.
	gatkNeedsBam := gatkReadsNeedBam(refFasta)

	// Live paths: updated in-place as each step produces a new file, so later
	// steps always have the correct path without any filesystem re-scan.
	//
	// sortedBamPath is live for the same reason the others are: the adopt step
	// writes it, the align steps write it, and MarkDuplicates reads it. Reading
	// task.state.SortedBam.Path directly instead — as ReasonMarkDupSortedBam
	// used to — is a null path on the adoption route, where nothing was in the
	// sorted slot when the plan was built.
	sortedBamPath := task.state.SortedBam.Path
	rgmdBamPath := task.state.RgmdBam.Path
	rgmdCramPath := task.state.RgmdCram.Path
	bqsrBamPath := task.state.BqsrBam.Path
	bqsrCramPath := task.state.BqsrCram.Path

	// Everything this run wrote, so cleanup can tell its own intermediates from
	// files it must not touch. See cleanupSampleOutputs.
	var produced []string
	wrote := func(paths ...string) {
		produced = append(produced, paths...)
	}

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

			sortedBamPath = filepath.Join(task.bamDir, sn+".sorted.bam")
			color.Cyan("[%s] Aligning PE reads → %s using %s ...\n", sn, sortedBamPath, task.aligner)
			if err := alignPairedReads(task.fwd, task.rev, refFasta, sortedBamPath, sn, task.aligner, opts.threads, verbose); err != nil {
				color.Red("[%s] Alignment failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(sortedBamPath)
			color.Green("[%s] Alignment done: %s\n", sn, sortedBamPath)

			rgmdBamPath = alignment.RgmdBamPath(sortedBamPath)
			color.Cyan("[%s] MarkDuplicates: sorted.bam → rgmd.bam ...\n", sn)
			if err := alignment.MarkDuplicates(refFasta, sortedBamPath, verbose, task.dupMarker, gatkLogLevel, opts.javaOpts, opts.threads); err != nil {
				color.Red("[%s] MarkDuplicates failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(rgmdBamPath)
			color.Green("[%s] MarkDuplicates done: %s\n", sn, rgmdBamPath)

			// The next step (ReasonConvertRgmdBam) is always in the plan after
			// ReasonAlignFromReads; set rgmdCramPath so it has the right target.
			rgmdCramPath = strings.TrimSuffix(rgmdBamPath, filepath.Ext(rgmdBamPath)) + ".cram"

		// ------------------------------------------------------------------ //
		// Align the single long-read FASTQ with pbmm2 → sorted.bam → rgmd.bam //
		// ------------------------------------------------------------------ //
		case ReasonAlignLongReads:
			if task.se == "" {
				return sampleResult{sample: sn, err: fmt.Errorf("no usable alignment intermediates or long-read FASTQ found")}
			}
			if !skipVer {
				// ValidateFastqGz is gzip-based in both its quick and its full
				// form, so it can only speak for a compressed file. A plain
				// .fastq is a real layout, not an error, so say what was
				// skipped rather than inventing a check for it.
				if strings.HasSuffix(strings.ToLower(task.se), ".gz") {
					color.Cyan("[%s] Validating long reads: %s\n", sn, task.se)
					if err := utils.ValidateFastqGz(task.se, verbose, quick); err != nil {
						color.Red("[%s] Long read validation failed: %v\n", sn, err)
						return sampleResult{sample: sn, err: fmt.Errorf("long read validation failed: %w", err)}
					}
					color.Green("[%s] Long reads are valid\n", sn)
				} else {
					color.Yellow("[%s] %s is not gzipped — skipping the integrity check\n", sn, task.se)
				}
			}

			sortedBamPath = filepath.Join(task.bamDir, sn+".sorted.bam")
			color.Cyan("[%s] Aligning long reads → %s using pbmm2 (preset %s) ...\n", sn, sortedBamPath, task.preset)
			if err := alignLongReads(task.se, task.refIndex, sortedBamPath, sn, task.preset, opts.threads, verbose); err != nil {
				color.Red("[%s] Alignment failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(sortedBamPath)
			color.Green("[%s] Alignment done: %s\n", sn, sortedBamPath)

			rgmdBamPath = alignment.RgmdBamPath(sortedBamPath)
			color.Cyan("[%s] Marking duplicates: sorted.bam → rgmd.bam ...\n", sn)
			if err := alignment.MarkDuplicates(refFasta, sortedBamPath, verbose, task.dupMarker, gatkLogLevel, opts.javaOpts, opts.threads); err != nil {
				color.Red("[%s] Duplicate marking failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(rgmdBamPath)
			if err := checkMarkedDuplicates(rgmdBamPath); err != nil {
				color.Red("[%s] %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			color.Green("[%s] Duplicate marking done: %s\n", sn, rgmdBamPath)

			rgmdCramPath = strings.TrimSuffix(rgmdBamPath, filepath.Ext(rgmdBamPath)) + ".cram"

		// ------------------------------------------------------------------ //
		// Adopt an externally produced alignment → sorted.bam                 //
		// ------------------------------------------------------------------ //
		case ReasonAdoptExternalAlignment:
			if task.external == "" {
				return sampleResult{sample: sn, err: fmt.Errorf("adoption planned but no external alignment was recorded")}
			}
			color.Cyan("[%s] Adopting externally produced alignment: %s\n", sn, task.external)

			// Strict, and stricter than checkAlignmentContigs: adoption is new
			// trust in a file this pipeline never made, so a cram aligned to a
			// different build of a same-named assembly must not slip through on
			// a name-subset match.
			if err := verifyAdoptable(task.external, refFasta); err != nil {
				color.Red("[%s] %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			if !skipVer {
				color.Cyan("[%s] Validating %s ...\n", sn, filepath.Base(task.external))
				if err := utils.ValidateBam(task.external, refFasta, verbose, quick); err != nil {
					color.Red("[%s] Adopted alignment failed validation: %v\n", sn, err)
					return sampleResult{sample: sn, err: fmt.Errorf("adopted alignment %s failed validation: %w", task.external, err)}
				}
				color.Green("[%s] Adopted alignment is valid\n", sn)
			}

			sortedBamPath = filepath.Join(task.bamDir, sn+".sorted.bam")
			if err := normaliseAdoptedAlignment(task.external, sortedBamPath, refFasta, sn, opts.threads, verbose); err != nil {
				color.Red("[%s] Adoption failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(sortedBamPath)
			color.Green("[%s] Adopted as %s (%s left untouched)\n", sn, filepath.Base(sortedBamPath), filepath.Base(task.external))

		// ------------------------------------------------------------------ //
		// sorted.bam → rgmd.bam → rgmd.cram (MarkDuplicates path)            //
		// ------------------------------------------------------------------ //
		case ReasonMarkDupSortedBam:
			if sortedBamPath == "" {
				return sampleResult{sample: sn, err: fmt.Errorf("no sorted.bam to mark duplicates on")}
			}
			rgmdBamPath = alignment.RgmdBamPath(sortedBamPath)
			color.Cyan("[%s] MarkDuplicates: sorted.bam → rgmd.bam ...\n", sn)
			if err := alignment.MarkDuplicates(refFasta, sortedBamPath, verbose, task.dupMarker, gatkLogLevel, opts.javaOpts, opts.threads); err != nil {
				color.Red("[%s] MarkDuplicates failed: %v\n", sn, err)
				return sampleResult{sample: sn, err: err}
			}
			wrote(rgmdBamPath)
			if task.longRead {
				if err := checkMarkedDuplicates(rgmdBamPath); err != nil {
					color.Red("[%s] %v\n", sn, err)
					return sampleResult{sample: sn, err: err}
				}
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
			wrote(rgmdCramPath)
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
			wrote(rgmdBamPath)
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
			wrote(bqsrBamPath)
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
			wrote(bqsrCramPath)
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
			if err := cleanupSampleOutputs(task.bamDir, sn, produced, rgmdCramPath, bqsrCramPath); err != nil {
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
		fullPath := filepath.Join(bamDir, name)

		// classifyAlignmentFile (exports.go) holds the suffix rules, shared
		// with InspectBamDirMetadata so a scan and the pipeline cannot
		// disagree about what an RGMD cram is.
		slot := classifyAlignmentFile(name)
		if slot == slotIgnore {
			continue
		}

		if slot == slotOther {
			fileInfo, iErr := entry.Info()
			size := int64(0)
			if iErr == nil {
				size = fileInfo.Size()
			}
			state.OtherFiles = append(state.OtherFiles, FileInfo{
				Path:    fullPath,
				Size:    size,
				Present: true,
			})
			continue
		}

		if field := state.slotField(slot); field != nil {
			*field = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
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
	// Clear the target and every index spelling beside it. This used to remove
	// only a .bai for a BAM, which BamIndex has not written since it moved to
	// CSI: the stale .csi survived the re-alignment, and findIndex validates
	// and reuses an index it finds, against whatever data is now under that
	// name.
	if err := clearAlignment(sortedBam); err != nil {
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

// pipelineIntermediates lists the basenames this pipeline writes for one sample
// on the way to its durable CRAMs.
//
// Only these names, for this sample, are collectable by cleanup. It is a
// closed list rather than a pattern because it is the allow-list for deletion:
// anything not derivable from the sample's own name is, by definition, not
// something this pipeline produced.
func pipelineIntermediates(sample string) []string {
	stems := []string{
		sample + ".sorted.bam",
		sample + ".sorted.partial.bam",
		sample + ".RGMD.bam",
		sample + ".RGMD_bqsr.bam",
		sample + ".RGMD.cram.partial",
	}

	names := make([]string, 0, len(stems)*5)
	for _, stem := range stems {
		names = append(names, stem)
		// The index spellings the pipeline might have written beside it, and
		// the .partial an interrupted index or decode leaves behind.
		for _, idx := range indexCandidates(stem) {
			names = append(names, idx, idx+".partial")
		}
		names = append(names, stem+".partial")
	}
	return names
}

// cleanupSampleOutputs removes this sample's intermediates from bamDir, leaving
// the durable CRAMs — and anything the pipeline did not create.
//
// It deletes a file only when this run produced it, or when its name is one
// this pipeline writes for this sample (pipelineIntermediates). Everything else
// is left in place and reported.
//
// That polarity is the whole point, and it is not merely tidiness. A bams
// directory can hold data the pipeline did not create: a long-read sample may
// arrive with an externally produced MENINA_LONG_READS.aligned.cram, which is
// adopted as *input* and — where no long-read FASTQ was kept — is the only copy
// of that sample's reads in existence. The previous rule was "keep four paths,
// delete every other regular file", which removed it. Worse, it removed it on
// the *second* run: a finished sample is re-queued with a cleanup-only plan, so
// the deletion happened a run after the adoption.
//
// produced is what this process actually wrote, accumulated by processSample as
// each step succeeded. rgmdCramPath and bqsrCramPath are the durable outputs to
// keep; an empty rgmdCramPath means the sample never produced one, which is a
// refusal rather than a licence to delete.
func cleanupSampleOutputs(bamDir, sample string, produced []string, rgmdCramPath, bqsrCramPath string) error {
	// Refuse to clean a directory whose outcome was never established. Without
	// this, a bug that left rgmdCramPath empty turned cleanup into "delete
	// everything", which is the most destructive thing this pipeline can do.
	if rgmdCramPath == "" {
		return fmt.Errorf("refusing to clean %s: no rgmd.cram was produced for %s", bamDir, sample)
	}
	if _, err := os.Stat(rgmdCramPath); err != nil {
		return fmt.Errorf("refusing to clean %s: rgmd.cram %s is not there: %w", bamDir, rgmdCramPath, err)
	}

	keep := make(map[string]struct{})

	addWithIndex := func(cramPath string) {
		keep[cramPath] = struct{}{}
		for _, idx := range indexCandidates(cramPath) {
			if _, err := os.Stat(idx); err == nil {
				keep[idx] = struct{}{}
			}
		}
	}

	addWithIndex(rgmdCramPath)
	if bqsrCramPath != "" {
		addWithIndex(bqsrCramPath)
	}

	// Collectable: written by this run, or a name only this pipeline writes for
	// this sample. The second clause is what keeps an interrupted earlier run's
	// leftovers collectable, so a sample can reach the "already complete" state
	// again.
	collectable := make(map[string]struct{})
	for _, p := range produced {
		collectable[filepath.Base(p)] = struct{}{}
		for _, idx := range indexCandidates(p) {
			collectable[filepath.Base(idx)] = struct{}{}
		}
	}
	for _, name := range pipelineIntermediates(sample) {
		collectable[name] = struct{}{}
	}

	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return err
	}

	var left []string
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()
		fullPath := filepath.Join(bamDir, name)
		if _, ok := keep[fullPath]; ok {
			continue
		}
		// Always preserve log/report files.
		lower := strings.ToLower(name)
		if strings.HasSuffix(lower, ".txt") || strings.HasSuffix(lower, ".pdf") {
			continue
		}
		if _, ok := collectable[name]; !ok {
			left = append(left, name)
			continue
		}
		if err := os.Remove(fullPath); err != nil && !os.IsNotExist(err) {
			return fmt.Errorf("removing %s: %w", fullPath, err)
		}
	}

	// Said out loud because the safe choice is also the one that silently
	// accumulates disk: an operator who wants these gone has to know they are
	// there, and has to be the one to decide.
	if len(left) > 0 {
		sort.Strings(left)
		color.Yellow("[%s] left %d file(s) in %s that this run did not produce: %s\n",
			sample, len(left), bamDir, strings.Join(left, ", "))
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
		// Not os.Stdout: os/exec hands an *os.File straight to the child, whose
		// writes then reach the terminal without passing through anything a
		// progress bar can coordinate with. See utils.LogWriter.
		out := utils.LogWriter()
		defer out.Close()
		cmd.Stdout = out
		cmd.Stderr = out
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
