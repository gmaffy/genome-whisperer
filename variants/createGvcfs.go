package variants

import (
	"bufio"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/schollz/progressbar/v3"
)

// Options is everything the pipeline needs, for every stage and every input
// mode. Data-dir mode sets DataDir; config and inline modes set Bams and OutDir.
// Nothing below the stage entry points has to ask which mode it is running in.
type Options struct {
	DataDir string   // walk the standard data directory tree
	Bams    []string // explicit bam/cram paths (config or inline)
	OutDir  string   // output root; ignored in data-dir mode

	Species  string
	RefVer   string
	RefFasta string

	Caller       string // "gatk" or "deepvariant"
	Merger       string // "gatk" or "glnexus"
	DVVer        string // DeepVariant container version
	ModelType    string // DeepVariant model; empty means infer from the sample name
	GatkLogLevel string

	// Filtering. HardFilter carries the GATK thresholds; MinGQ is the extra
	// genotype-quality cutoff used by the DeepVariant profile, whose output has
	// none of the GATK annotations.
	HardFilter utils.HardFilterConfig
	MinGQ      int

	NoMerging    bool
	NoHardFilter bool

	Threads          int
	NoBqsr           bool // data-dir mode: use the rgmd bam/cram instead of bqsr
	Verbose          bool
	Quick            bool
	SkipVerification bool
}

// absRoots resolves DataDir, OutDir and RefFasta to absolute paths. Every stage
// calls it up front so the paths it prints, and the paths it passes to GATK,
// DeepVariant and GLnexus, do not depend on the current working directory.
func absRoots(opts Options) (Options, error) {
	var err error
	if opts.DataDir != "" {
		if opts.DataDir, err = filepath.Abs(opts.DataDir); err != nil {
			return opts, fmt.Errorf("resolving data directory: %w", err)
		}
	}
	if opts.OutDir != "" {
		if opts.OutDir, err = filepath.Abs(opts.OutDir); err != nil {
			return opts, fmt.Errorf("resolving output directory: %w", err)
		}
	}
	if opts.RefFasta != "" {
		if opts.RefFasta, err = filepath.Abs(opts.RefFasta); err != nil {
			return opts, fmt.Errorf("resolving reference fasta: %w", err)
		}
	}
	return opts, nil
}

// GvcfPath returns where a sample's gVCF for one chromosome belongs.
//
// Data-dir mode:  <sample>/reference_genomes/<refVer>/{gatk_gvcfs|dv_gvcfs}/<sample>.<chrom>.g.vcf.gz
// Otherwise:      <OutDir>/<chrom>/{gatk_gvcfs|dv_gvcfs}/<sample>.<chrom>.g.vcf.gz
//
// This is the only place either convention is written down. MergeGvcfs must call
// this rather than rebuilding paths of its own.
func GvcfPath(opts Options, s SampleWork, chrom string) string {
	gvcfDirName := "gatk_gvcfs"
	if strings.ToLower(opts.Caller) != "gatk" {
		gvcfDirName = "dv_gvcfs"
	}

	label := strings.ReplaceAll(chrom, ".", "_")
	base := filepath.Base(s.Cram)
	base = strings.TrimSuffix(base, filepath.Ext(base))
	gvcfName := fmt.Sprintf("%s.%s.g.vcf.gz", base, label)

	if opts.DataDir != "" {
		return filepath.Join(s.CramDir, opts.RefVer, gvcfDirName, gvcfName)
	}
	return filepath.Join(opts.OutDir, label, gvcfDirName, gvcfName)
}

// FindSampleAlignments resolves opts into the samples to process. In data-dir mode it
// walks the tree and picks each sample's alignment; otherwise it wraps each entry
// in opts.Bams. Samples it cannot use are returned in skipped rather than
// aborting the run.
func FindSampleAlignments(opts Options) ([]SampleWork, []string, error) {
	var samples []SampleWork
	var skipped []string

	numWorkers := opts.Threads
	if numWorkers <= 0 {
		numWorkers = runtime.NumCPU()
	}

	type evalResult struct {
		sampleWork *SampleWork
		skippedID  string
	}

	// ------------------------------------ config / inline ------------------------------------ //
	if opts.DataDir == "" {
		if numWorkers > len(opts.Bams) {
			numWorkers = len(opts.Bams)
		}
		if numWorkers <= 0 {
			numWorkers = 1
		}

		jobs := make(chan string, len(opts.Bams))
		for _, bam := range opts.Bams {
			jobs <- bam
		}
		close(jobs)

		results := make(chan evalResult, len(opts.Bams))
		inlineDesc := "Checking alignment files"
		if !opts.SkipVerification {
			if opts.Quick {
				inlineDesc = "Checking alignment files with quick validation"
			} else {
				inlineDesc = "Checking alignment files with deep validation"
			}
		}
		bar := progressbar.Default(int64(len(opts.Bams)), inlineDesc)

		var wg sync.WaitGroup
		for i := 0; i < numWorkers; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for bam := range jobs {
					_ = bar.Add(1)
					if _, err := os.Stat(bam); err != nil {
						results <- evalResult{skippedID: bam}
						continue
					}
					if !opts.SkipVerification {
						if valErr := utils.ValidateBam(bam, opts.RefFasta, opts.Verbose, opts.Quick); valErr != nil {
							results <- evalResult{skippedID: bam}
							continue
						}
					}
					base := filepath.Base(bam)
					results <- evalResult{
						sampleWork: &SampleWork{
							Sample: strings.TrimSuffix(base, filepath.Ext(base)),
							Cram:   bam,
						},
					}
				}
			}()
		}
		wg.Wait()
		close(results)
		_ = bar.Finish()

		for r := range results {
			if r.sampleWork != nil {
				samples = append(samples, *r.sampleWork)
			} else if r.skippedID != "" {
				skipped = append(skipped, r.skippedID)
			}
		}

		sort.Slice(samples, func(i, j int) bool {
			return samples[i].Sample < samples[j].Sample
		})
		sort.Strings(skipped)

		return samples, skipped, nil
	}

	// --------------------------------------- data-dir ---------------------------------------- //
	dataDirAbs, err := filepath.Abs(opts.DataDir)
	if err != nil {
		return nil, nil, fmt.Errorf("getting absolute path for data directory %s: %w", opts.DataDir, err)
	}

	// Lower-cased to match FindBams and JointVcfDir. On a case-sensitive
	// filesystem, mixing the two spellings means sample discovery and BAM
	// discovery look in different directories.
	pattern := filepath.Join(opts.DataDir, strings.ToLower(opts.Species), "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, nil, fmt.Errorf("finding samples in %s: %w", pattern, err)
	}

	seen := make(map[string]struct{}, len(matches))
	var sampleList []string
	for _, match := range matches {
		sample := filepath.Base(filepath.Dir(match))
		if _, ok := seen[sample]; ok {
			continue
		}
		seen[sample] = struct{}{}
		sampleList = append(sampleList, sample)
	}

	if numWorkers > len(sampleList) {
		numWorkers = len(sampleList)
	}
	if numWorkers <= 0 {
		numWorkers = 1
	}

	jobs := make(chan string, len(sampleList))
	for _, s := range sampleList {
		jobs <- s
	}
	close(jobs)

	results := make(chan evalResult, len(sampleList))
	barDesc := "Discovering alignments"
	if !opts.SkipVerification {
		if opts.Quick {
			barDesc = "Discovering alignments with quick validation"
		} else {
			barDesc = "Discovering alignments with deep validation"
		}
	}
	bar := progressbar.Default(int64(len(sampleList)), barDesc)

	var wg sync.WaitGroup
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for sample := range jobs {
				crams, bams, fErr := FindBams(dataDirAbs, opts.Species, sample, opts.RefVer, opts.NoBqsr)
				_ = bar.Add(1)
				if fErr != nil {
					results <- evalResult{skippedID: sample}
					continue
				}

				alignments := crams
				if len(alignments) == 0 {
					alignments = bams
				}
				switch len(alignments) {
				case 1:
					alignmentFile := alignments[0]
					if !opts.SkipVerification {
						if valErr := utils.ValidateBam(alignmentFile, opts.RefFasta, opts.Verbose, opts.Quick); valErr != nil {
							results <- evalResult{skippedID: sample}
							continue
						}
					}
					results <- evalResult{
						sampleWork: &SampleWork{
							Sample:  sample,
							Cram:    alignmentFile,
							CramDir: filepath.Dir(filepath.Dir(filepath.Dir(alignmentFile))),
						},
					}
				default:
					results <- evalResult{skippedID: sample}
				}
			}
		}()
	}

	wg.Wait()
	close(results)
	_ = bar.Finish()

	for r := range results {
		if r.sampleWork != nil {
			samples = append(samples, *r.sampleWork)
		} else if r.skippedID != "" {
			skipped = append(skipped, r.skippedID)
		}
	}

	sort.Slice(samples, func(i, j int) bool {
		return samples[i].Sample < samples[j].Sample
	})
	sort.Strings(skipped)

	return samples, skipped, nil
}

// gvcfJob is one sample × chromosome unit of gVCF work.
type gvcfJob struct {
	sample SampleWork
	chrom  string
	seqs   []SeqInfo
}

// gvcfJobResult is what discoverValidGvcfs reports for one gvcfJob: either an
// existing, valid gVCF to reuse (valid==true), or a job left for creation,
// tagged with whether a file was there at all (corrupt) so the create step's
// log line can say "creating" vs "re-creating corrupt".
type gvcfJobResult struct {
	job     gvcfJob
	path    string
	valid   bool
	corrupt bool
}

// listGvcfDirs lists, once each, every directory a data-dir mode job's gVCF
// could live in. Each sample owns its own directory
// (<CramDir>/<RefVer>/{gatk_gvcfs|dv_gvcfs}), so this costs one ReadDir per
// sample rather than one per sample × chromosome. A directory that does not
// exist yet (nothing created there so far) is simply recorded as empty.
func listGvcfDirs(opts Options, jobs []gvcfJob) map[string][]os.DirEntry {
	dirs := make(map[string][]os.DirEntry)
	for _, j := range jobs {
		dir := filepath.Dir(GvcfPath(opts, j.sample, j.chrom))
		if _, seen := dirs[dir]; seen {
			continue
		}
		entries, _ := os.ReadDir(dir)
		dirs[dir] = entries
	}
	return dirs
}

// findExistingGvcf looks for job's gVCF in dirEntries's listing of its
// directory, matching by the "<chrom>.g.vcf.gz" suffix rather than trusting a
// filename rebuilt from the sample's current alignment file. That tolerates
// the alignment having been renamed since the gVCF was created (an added
// "rgmd"/"bqsr" tag, a leftover ".cram" that an older run never stripped,
// etc.) — exactly the drift that made data-dir mode gVCFs go undiscovered.
//
// More than one match (e.g. a stale file in the old naming alongside a
// current one) is reported and resolved by taking the most recently modified.
func findExistingGvcf(dirEntries map[string][]os.DirEntry, dir, sample, chrom string) string {
	suffix := "." + strings.ReplaceAll(chrom, ".", "_") + ".g.vcf.gz"

	var matches []os.DirEntry
	for _, e := range dirEntries[dir] {
		if !e.IsDir() && strings.HasSuffix(e.Name(), suffix) {
			matches = append(matches, e)
		}
	}

	switch len(matches) {
	case 0:
		return ""
	case 1:
		return filepath.Join(dir, matches[0].Name())
	}

	sort.Slice(matches, func(a, b int) bool { return matches[a].Name() < matches[b].Name() })
	names := make([]string, len(matches))
	newest := matches[0]
	var newestMod time.Time
	for i, e := range matches {
		names[i] = e.Name()
		if info, iErr := e.Info(); iErr == nil && info.ModTime().After(newestMod) {
			newestMod = info.ModTime()
			newest = e
		}
	}
	color.Yellow("[%s] multiple gVCFs match %s in %s: %v — using most recently modified\n", sample, chrom, dir, names)
	return filepath.Join(dir, newest.Name())
}

// discoverValidGvcfs checks every job against what is already on disk. A
// missing gVCF is left for creation; an existing one is validated (unless
// --skip-verification) and either reused or, if corrupt, removed along with
// its index and left for re-creation. It runs with the same worker count as
// gVCF creation itself, since ValidateGvcf can be as slow as decompressing the
// file.
//
// Data-dir mode locates the existing file, if any, by scanning the sample's
// gVCF directory (see findExistingGvcf) rather than assuming GvcfPath's
// convention matches what is on disk. Inline/config mode keeps the exact
// GvcfPath + stat check: its gVCF directory holds every sample's output
// together, so a suffix match alone cannot tell samples apart.
func discoverValidGvcfs(opts Options, jobs []gvcfJob, maxParallel int) (valid map[string][]string, remaining []gvcfJobResult, reused int) {
	valid = make(map[string][]string)

	var dirEntries map[string][]os.DirEntry
	if opts.DataDir != "" {
		dirEntries = listGvcfDirs(opts, jobs)
	}

	jobCh := make(chan gvcfJob, len(jobs))
	for _, j := range jobs {
		jobCh <- j
	}
	close(jobCh)

	resultsCh := make(chan gvcfJobResult, len(jobs))

	desc := "Discovering gVCFs"
	if !opts.SkipVerification {
		if opts.Quick {
			desc = "Discovering valid gVCFs (quick validation)"
		} else {
			desc = "Discovering valid gVCFs (deep validation)"
		}
	}
	bar := progressbar.Default(int64(len(jobs)), desc)

	var wg sync.WaitGroup
	for w := 0; w < maxParallel; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobCh {
				canonicalPath := GvcfPath(opts, j.sample, j.chrom)
				existingPath := ""

				if opts.DataDir != "" {
					existingPath = findExistingGvcf(dirEntries, filepath.Dir(canonicalPath), j.sample.Sample, j.chrom)
				} else if _, sErr := os.Stat(canonicalPath); sErr == nil {
					existingPath = canonicalPath
				}
				_ = bar.Add(1)

				if existingPath == "" {
					resultsCh <- gvcfJobResult{job: j, path: canonicalPath}
					continue
				}
				if opts.SkipVerification {
					resultsCh <- gvcfJobResult{job: j, path: existingPath, valid: true}
					continue
				}
				if vErr := utils.ValidateGvcf(existingPath, opts.Verbose, opts.Quick); vErr != nil {
					// Remove the partial file and its index so a failed retry
					// cannot leave something that looks complete. Re-creation
					// always writes to the canonical path, so a repaired
					// sample stops carrying whatever legacy name it had.
					os.Remove(existingPath)
					os.Remove(existingPath + ".tbi")
					resultsCh <- gvcfJobResult{job: j, path: canonicalPath, corrupt: true}
					continue
				}
				resultsCh <- gvcfJobResult{job: j, path: existingPath, valid: true}
			}
		}()
	}
	wg.Wait()
	close(resultsCh)
	_ = bar.Finish()

	for r := range resultsCh {
		if r.valid {
			valid[r.job.chrom] = append(valid[r.job.chrom], r.path)
			reused++
			continue
		}
		remaining = append(remaining, r)
	}
	for chrom := range valid {
		sort.Strings(valid[chrom])
	}

	return valid, remaining, reused
}

// CreateGvcfs creates a gVCF for every sample × chromosome and returns them
// grouped by chromosome label, which is the shape MergeGvcfs needs. A chromosome
// is complete when len(result[chrom]) equals the number of samples.
//
// An existing valid gVCF is reused, a corrupt one is re-created, and a sample
// that fails is reported and skipped so the rest of the run continues.
func CreateGvcfs(opts Options) (map[string][]string, error) {

	// ==================================== Validate inputs ===================================== //

	if opts.Species == "" {
		return nil, fmt.Errorf("species name must not be empty")
	}
	if opts.RefVer == "" {
		return nil, fmt.Errorf("reference version must not be empty")
	}
	if opts.RefFasta == "" {
		return nil, fmt.Errorf("reference fasta must not be empty")
	}

	fastaInfo, err := os.Stat(opts.RefFasta)
	if err != nil {
		return nil, fmt.Errorf("accessing reference fasta %s: %w", opts.RefFasta, err)
	}
	if !fastaInfo.Mode().IsRegular() {
		return nil, fmt.Errorf("reference fasta %s is not a regular file", opts.RefFasta)
	}

	dictFilePath := opts.RefFasta[:len(opts.RefFasta)-len(filepath.Ext(opts.RefFasta))] + ".dict"
	if _, dErr := os.Stat(dictFilePath); dErr != nil {
		return nil, fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	caller := strings.ToLower(opts.Caller)
	if caller != "gatk" && caller != "deepvariant" {
		return nil, fmt.Errorf("caller must be gatk or deepvariant, got %q", opts.Caller)
	}

	if opts.DataDir == "" && len(opts.Bams) == 0 {
		return nil, fmt.Errorf("provide either a data directory or at least one bam file")
	}
	if opts.DataDir != "" && len(opts.Bams) > 0 {
		return nil, fmt.Errorf("provide either a data directory or bam files, not both")
	}
	if opts.DataDir == "" && opts.OutDir == "" {
		return nil, fmt.Errorf("an output directory is required when not using a data directory")
	}
	if opts.Threads <= 0 {
		opts.Threads = 4
	}

	// Make the roots absolute so every path logged or handed to GATK is a full
	// path, whatever the caller passed and whatever the working directory is.
	// opts is a value, so this only affects this call.
	if opts, err = absRoots(opts); err != nil {
		return nil, err
	}

	color.Green("All file paths valid\n....................................................\n\n")

	// ================================== Chroms and contigs ==================================== //

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return nil, fmt.Errorf("getting chromosomes and contigs: %w", err)
	}

	// =================================== Discover samples ===================================== //

	color.Cyan("================================== Discovering samples =================================\n\n")

	samples, skipped, err := FindSampleAlignments(opts)
	if err != nil {
		return nil, err
	}
	color.Green("\nAlignment Discovery Summary:\n")
	color.Green("  - Usable sample alignments found: %d\n", len(samples))
	if len(skipped) > 0 {
		color.Red("  - Samples with no valid alignments: %d %v\n", len(skipped), skipped)
		color.Red("  - Cannot proceed. Run AlignReads first.\n")
		color.Cyan("%s\n\n", strings.Repeat("=", 90))
		return nil, fmt.Errorf("Sample(s) : %v have invalid alignment files. Cannot proceed. Run AlignReads first.", skipped)
	}

	color.Green("  - Samples with no valid alignments: 0\n")
	color.Green("All samples have valid alignment files. Proceeding...\n")
	color.Cyan("%s\n\n", strings.Repeat("=", 80))

	// ===================================== Build job list ===================================== //

	var jobs []gvcfJob
	for _, s := range samples {
		for _, c := range chroms {
			jobs = append(jobs, gvcfJob{sample: s, chrom: c.ID, seqs: []SeqInfo{c}})
		}
		if len(contigs) > 0 {
			jobs = append(jobs, gvcfJob{sample: s, chrom: "contigs", seqs: contigs})
		}
	}

	totalCores := runtime.NumCPU()
	maxParallel := totalCores / opts.Threads
	if maxParallel < 1 {
		maxParallel = 1
	}
	color.Cyan("Machine has %d cores. %d threads per job -> max %d parallel jobs\n\n",
		totalCores, opts.Threads, maxParallel)

	// ================================ Discover valid gVCFs ==================================== //

	color.Cyan("================================ Discovering valid gVCFs ===============================\n\n")

	gvcfs, remaining, reused := discoverValidGvcfs(opts, jobs, maxParallel)

	color.Green("\ngVCF Discovery Summary:\n")
	color.Green("  - Valid gVCFs found:  %d\n", reused)
	color.Yellow("  - gVCFs to create:    %d\n", len(remaining))
	color.Cyan("%s\n\n", strings.Repeat("=", 90))

	// ======================================== Run jobs ======================================== //

	jobCh := make(chan gvcfJobResult, len(remaining))
	for _, r := range remaining {
		jobCh <- r
	}
	close(jobCh)

	var (
		mu          sync.Mutex
		wg          sync.WaitGroup
		failedTasks []FailedTask
		created     int
	)

	// One tick per job regardless of outcome, so the bar always reaches 100%.
	// AddDetail keeps a rolling window of the most recent worker activity
	// rendered above the bar (one row per worker, capped so it never outgrows
	// a typical terminal); the description underneath it carries a running
	// done/remaining/failed tally.
	detailRows := maxParallel
	if detailRows > 10 {
		detailRows = 10
	}
	bar := progressbar.NewOptions64(
		int64(len(remaining)),
		progressbar.OptionSetDescription("Creating gVCFs"),
		progressbar.OptionSetWriter(os.Stderr),
		progressbar.OptionSetWidth(10),
		progressbar.OptionShowTotalBytes(true),
		progressbar.OptionThrottle(65*time.Millisecond),
		progressbar.OptionShowCount(),
		progressbar.OptionShowIts(),
		progressbar.OptionOnCompletion(func() { fmt.Fprint(os.Stderr, "\n") }),
		progressbar.OptionSpinnerType(14),
		progressbar.OptionFullWidth(),
		progressbar.OptionSetRenderBlankState(true),
		progressbar.OptionSetMaxDetailRow(detailRows),
	)

	// reportProgress refreshes the description with how many jobs are done,
	// left and failed. Called with the current failedTasks count still held
	// under mu, so it never races the tally it reports.
	reportProgress := func(failed int) {
		st := bar.State()
		msg := fmt.Sprintf("Creating gVCFs (%d done, %d remaining", st.CurrentNum, st.Max-st.CurrentNum)
		if failed > 0 {
			msg += fmt.Sprintf(", %d failed", failed)
		}
		bar.Describe(msg + ")")
	}

	for w := 1; w <= maxParallel; w++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			for r := range jobCh {
				j := r.job
				theGVCF := r.path

				if mkErr := os.MkdirAll(filepath.Dir(theGVCF), 0755); mkErr != nil {
					_ = bar.Add(1)
					color.Red("[Worker %d] [%s] creating %s: %v\n\n", workerID, j.sample.label(), filepath.Dir(theGVCF), mkErr)
					mu.Lock()
					failedTasks = append(failedTasks, FailedTask{Sample: j.sample.Sample, Chrom: j.chrom, Reason: mkErr})
					failed := len(failedTasks)
					mu.Unlock()
					reportProgress(failed)
					continue
				}

				verb := "creating"
				if r.corrupt {
					verb = "re-creating corrupt"
				}

				// Honour --model-type when it is set. Only fall back to the data
				// directory's "...lr" long-read naming convention when it is not,
				// so the flag is never silently overridden.
				modelType := opts.ModelType
				if modelType == "" {
					modelType = "WGS"
					if strings.HasSuffix(strings.ToLower(j.sample.Sample), "lr") {
						modelType = "PACBIO"
					}
				}

				if caller == "gatk" {
					_ = bar.AddDetail(fmt.Sprintf("[Worker %d] %s %s %s gvcf with GATK", workerID, verb, j.sample.label(), j.chrom))
				} else {
					_ = bar.AddDetail(fmt.Sprintf("[Worker %d] %s %s %s gvcf with DeepVariant (%s)", workerID, verb, j.sample.label(), j.chrom, modelType))
				}

				var cErr error
				if caller == "gatk" {
					_, cErr = CreateGvcfGATK(j.sample.Cram, opts.RefFasta, j.seqs, theGVCF, opts.GatkLogLevel, opts.Verbose)
				} else {
					_, cErr = CreateGvcfDV(j.sample.Cram, opts.RefFasta, j.seqs, theGVCF, opts.DVVer, modelType, opts.Verbose, opts.Threads)
				}
				_ = bar.Add(1)

				if cErr != nil {
					_ = bar.AddDetail(fmt.Sprintf("[Worker %d] %s %s: FAILED", workerID, j.sample.label(), j.chrom))
					color.Red("[Worker %d] [%s] creating gVCF for %s FAILED: %v\n\n", workerID, j.sample.label(), j.chrom, cErr)
					mu.Lock()
					failedTasks = append(failedTasks, FailedTask{Sample: j.sample.Sample, Chrom: j.chrom, Reason: cErr})
					failed := len(failedTasks)
					mu.Unlock()
					reportProgress(failed)
					continue
				}

				_ = bar.AddDetail(fmt.Sprintf("[Worker %d] %s %s: done", workerID, j.sample.label(), j.chrom))
				mu.Lock()
				gvcfs[j.chrom] = append(gvcfs[j.chrom], theGVCF)
				created++
				failed := len(failedTasks)
				mu.Unlock()
				reportProgress(failed)
			}
		}(w)
	}
	wg.Wait()
	_ = bar.Finish()

	// ======================================== Summary ========================================= //

	// Sort each chromosome's gVCFs so the sample map and any downstream comparison
	// are deterministic; workers append in whatever order they finish.
	for chrom := range gvcfs {
		sort.Strings(gvcfs[chrom])
	}

	fmt.Printf("\n%s CREATE GVCFS SUMMARY %s\n", strings.Repeat("=", 25), strings.Repeat("=", 25))
	fmt.Printf("Samples:            %d\n", len(samples))
	fmt.Printf("Chromosome groups:  %d\n", len(gvcfs))
	fmt.Printf("gVCFs created:      %d\n", created)
	fmt.Printf("gVCFs reused:       %d\n", reused)

	if len(skipped) > 0 {
		color.Yellow("Samples skipped:    %d %v\n", len(skipped), skipped)
	}
	if len(failedTasks) > 0 {
		color.Red("Failed tasks:       %d\n", len(failedTasks))
		for _, f := range failedTasks {
			fmt.Printf("  - Sample: %s, Chrom: %s, Reason: %v\n", f.Sample, f.Chrom, f.Reason)
		}
	}

	// Report which chromosome groups are short of a full sample set, so the caller
	// can decide whether to merge them.
	for _, c := range chroms {
		if len(gvcfs[c.ID]) != len(samples) {
			color.Yellow("[%s] incomplete: %d/%d samples\n", c.ID, len(gvcfs[c.ID]), len(samples))
		}
	}
	if len(contigs) > 0 && len(gvcfs["contigs"]) != len(samples) {
		color.Yellow("[contigs] incomplete: %d/%d samples\n", len(gvcfs["contigs"]), len(samples))
	}
	fmt.Printf("%s\n\n", strings.Repeat("=", 72))

	if len(gvcfs) == 0 {
		return nil, fmt.Errorf("no gVCFs were produced (%d task(s) failed)", len(failedTasks))
	}
	return gvcfs, nil
}

// ---------------------------------------------------------------------------
// Absorbed from the retired RunVariantCaller.go / RunVariantCallerDir.go
// ---------------------------------------------------------------------------

type FailedTask struct {
	Sample string
	Chrom  string
	Reason error
}

type SampleWork struct {
	Sample  string
	Cram    string
	CramDir string // parent of the "bams" directory; used to derive gvcf output path
}

// projectDir returns the directory between the species and the sample in
// data-dir mode's <DataDir>/<species>/<project>/<sample>/... layout — in this
// dataset it is always a year (e.g. "2024"). It is empty in inline/config
// mode, where SampleWork has no CramDir and there is no such directory.
func (s SampleWork) projectDir() string {
	if s.CramDir == "" {
		return ""
	}
	return filepath.Base(filepath.Dir(filepath.Dir(s.CramDir)))
}

// label formats the sample for a log line, tagging on its project/year
// directory when one exists so "which run is this sample from" doesn't
// require looking anything up.
func (s SampleWork) label() string {
	if year := s.projectDir(); year != "" {
		return fmt.Sprintf("%s (%s)", s.Sample, year)
	}
	return s.Sample
}

func FindBams(dataDirAbs string, species string, sample string, refVer string, noBqsr bool) ([]string, []string, error) {
	cramPat := fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*.cram", dataDirAbs, strings.ToLower(species), sample, refVer)
	bamPat := fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*.bam", dataDirAbs, strings.ToLower(species), sample, refVer)
	cramMatches, err := filepath.Glob(cramPat)
	if err != nil {
		return nil, nil, fmt.Errorf("glob error: %w", err)
	}
	bamMatches, err := filepath.Glob(bamPat)
	if err != nil {
		return nil, nil, fmt.Errorf("glob error: %w", err)
	}

	if len(cramMatches) == 0 && len(bamMatches) == 0 {
		return nil, nil, fmt.Errorf("no cram files or bams found in %s", cramPat)
	}
	//fmt.Printf("Found %d cram files\n", len(cramMatches))
	//fmt.Printf("Found %d bam files\n", len(bamMatches))
	var validCrams []string
	var validBams []string
	if len(cramMatches) > 0 {

		for _, cramMatch := range cramMatches {
			if noBqsr {
				cramLower := strings.ToLower(cramMatch)
				if strings.Contains(strings.ToLower(cramLower), "rgmd") {
					validCrams = append(validCrams, cramMatch)
				}
			} else {
				cramLower := strings.ToLower(cramMatch)
				if !strings.HasSuffix(strings.ToLower(sample), "lr") && strings.Contains(strings.ToLower(cramLower), "bqsr") {
					validCrams = append(validCrams, cramMatch)
				} else if strings.HasSuffix(strings.ToLower(sample), "lr") && strings.Contains(strings.ToLower(cramLower), "rgmd") {
					validCrams = append(validCrams, cramMatch)
				}
			}
		}
	} else if len(bamMatches) > 0 {
		for _, bamMatch := range bamMatches {
			if noBqsr {
				bamLower := strings.ToLower(bamMatch)
				if strings.Contains(strings.ToLower(bamLower), "rgmd") {
					validBams = append(validBams, bamMatch)
				}
			} else {
				bamLower := strings.ToLower(bamMatch)
				if !strings.HasSuffix(strings.ToLower(sample), "lr") && strings.Contains(strings.ToLower(bamLower), "bqsr") {
					validBams = append(validBams, bamMatch)
				} else if strings.HasSuffix(strings.ToLower(sample), "lr") && strings.Contains(strings.ToLower(bamLower), "rgmd") {
					validBams = append(validBams, bamMatch)
				}
			}
		}

	}

	return validCrams, validBams, nil

}

type SeqInfo struct {
	ID  string
	Len int
}

func getChromsAndContigs(dictFilePath string) ([]SeqInfo, []SeqInfo, error) {
	dictFile, err := os.Open(dictFilePath)
	if err != nil {
		return nil, nil, fmt.Errorf("opening reference dict file %s: %w", dictFilePath, err)
	}
	defer dictFile.Close()

	scanner := bufio.NewScanner(dictFile)
	var seqs []SeqInfo

	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "@SQ") {
			continue
		}
		parts := strings.Split(line, "\t")
		seqID := strings.TrimPrefix(parts[1], "SN:")
		seqLenStr := strings.TrimPrefix(parts[2], "LN:")
		seqLen, aErr := strconv.Atoi(seqLenStr)
		if aErr != nil {
			return nil, nil, fmt.Errorf("parsing LN field for sequence %q in dict file: %w", seqID, aErr)
		}
		seqs = append(seqs, SeqInfo{ID: seqID, Len: seqLen})
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, fmt.Errorf("scanning dict file: %w", err)
	}

	// Sort by length descending so the biggest sequences become "chroms".
	sort.Slice(seqs, func(i, j int) bool { return seqs[i].Len > seqs[j].Len })

	var chroms, contigs []SeqInfo
	if len(seqs) > 21 {
		chroms = append(chroms, seqs[:21]...)
		contigs = append(contigs, seqs[21:]...)
	} else {
		chroms = append(chroms, seqs...)
	}

	// Always promote MT and Pltd into the chroms group.
	for i := 0; i < len(contigs); {
		if contigs[i].ID == "MT" || contigs[i].ID == "Pltd" {
			chroms = append(chroms, contigs[i])
			contigs = append(contigs[:i], contigs[i+1:]...)
		} else {
			i++
		}
	}

	return chroms, contigs, nil
}

func CreateGvcfGATK(bam string, refFile string, chroms []SeqInfo, theGVCF string, gatkLogLevel string, verbose bool) (string, error) {
	var regionArg string
	if len(chroms) == 1 {
		regionArg = chroms[0].ID
	} else {
		f, err := os.CreateTemp("", "gatk_intervals_contigs_*.list")
		if err != nil {
			return "", fmt.Errorf("creating GATK interval list: %w", err)
		}
		defer os.Remove(f.Name())
		defer f.Close()
		for _, c := range chroms {
			fmt.Fprintln(f, c.ID)
		}
		regionArg = f.Name()
	}

	hapCmdStr := fmt.Sprintf(
		`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`,
		refFile, bam, regionArg, theGVCF, gatkLogLevel,
	)
	fmt.Printf("\n%s\n\n", hapCmdStr)
	return theGVCF, utils.RunCmd(hapCmdStr, verbose)
}

func CreateGvcfDV(bam string, refFile string, chroms []SeqInfo, theGVCF string, dvVer string, modelType string, verbose bool, threadsPerJob int) (string, error) {
	var regions string
	bamDir := filepath.Dir(bam)
	bamName := filepath.Base(bam)
	refDir := filepath.Dir(refFile)
	refName := filepath.Base(refFile)
	gvcfName := filepath.Base(theGVCF)
	gvcfDir := filepath.Dir(theGVCF)
	vcfName := strings.Replace(gvcfName, ".g.vcf.gz", ".vcf.gz", 1)
	if len(chroms) == 1 {
		regions = chroms[0].ID
	} else {
		f, err := os.CreateTemp(gvcfDir, "deepvariant_intervals_*.bed")
		if err != nil {
			return "", fmt.Errorf("creating DeepVariant interval list: %w", err)
		}
		defer os.Remove(f.Name())
		defer f.Close()
		for _, c := range chroms {
			fmt.Fprintf(f, "%s\t0\t%d\n", c.ID, c.Len)
		}
		regions = "/output/" + filepath.Base(f.Name()) // BED lives in gvcfDir, mounted as /output
	}

	safeRegion := strings.NewReplacer(string(os.PathSeparator), "_", ".", "_").Replace(regions)
	intermediateName := fmt.Sprintf("tmp_%s_%s", strings.TrimSuffix(bamName, filepath.Ext(bamName)), safeRegion)
	intermediatePath := filepath.Join(gvcfDir, intermediateName)

	// Remove intermediate directory if it exists to ensure a clean start for re-runs
	if _, err := os.Stat(intermediatePath); err == nil {
		slog.Info("Removing existing intermediate directory", "path", intermediatePath)
		if err := os.RemoveAll(intermediatePath); err != nil {
			return "", fmt.Errorf("removing existing intermediate directory %s: %w", intermediatePath, err)
		}
	}

	dvCmdStr := fmt.Sprintf(
		`docker run -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
			`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
			`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
			`--num_shards=%d --intermediate_results_dir /output/%s`,
		bamDir, refDir, gvcfDir, dvVer,
		modelType, refName, bamName,
		regions, vcfName, gvcfName,
		threadsPerJob, intermediateName,
	)

	//fmt.Printf("\n%s\n\n", dvCmdStr)
	if err := utils.RunCmd(dvCmdStr, verbose); err != nil {
		return "", err
	}

	// Clean up this job's DeepVariant intermediate results. Only this job's
	// directory is removed (not a tmp_* glob) so concurrent jobs writing to the
	// same gvcf directory are not disturbed.
	if rErr := os.RemoveAll(intermediatePath); rErr != nil {
		slog.Warn("removing DeepVariant intermediate directory", "path", intermediatePath, "err", rErr)
	}

	return theGVCF, nil
}
