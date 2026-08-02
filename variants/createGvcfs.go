package variants

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
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

// FindGvcfSamples resolves opts into the samples to process. In data-dir mode it
// walks the tree and picks each sample's alignment; otherwise it wraps each entry
// in opts.Bams. Samples it cannot use are returned in skipped rather than
// aborting the run.
func FindGvcfSamples(opts Options) ([]SampleWork, []string, error) {
	var samples []SampleWork
	var skipped []string

	// ------------------------------------ config / inline ------------------------------------ //
	if opts.DataDir == "" {
		for _, bam := range opts.Bams {
			if _, err := os.Stat(bam); err != nil {
				color.Red("[%s] not a valid file path\n", bam)
				skipped = append(skipped, bam)
				continue
			}
			base := filepath.Base(bam)
			samples = append(samples, SampleWork{
				Sample: strings.TrimSuffix(base, filepath.Ext(base)),
				Cram:   bam,
			})
		}
		return samples, skipped, nil
	}

	// --------------------------------------- data-dir ---------------------------------------- //
	dataDirAbs, err := filepath.Abs(opts.DataDir)
	if err != nil {
		return nil, nil, fmt.Errorf("getting absolute path for data directory %s: %w", opts.DataDir, err)
	}

	pattern := filepath.Join(opts.DataDir, opts.Species, "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, nil, fmt.Errorf("finding samples in %s: %w", pattern, err)
	}

	seen := make(map[string]struct{}, len(matches))
	for _, match := range matches {
		sample := filepath.Base(filepath.Dir(match))
		if _, ok := seen[sample]; ok {
			continue
		}
		seen[sample] = struct{}{}

		crams, bams, fErr := FindBams(dataDirAbs, opts.Species, sample, opts.RefVer, opts.NoBqsr)
		if fErr != nil {
			color.Red("[%s] no alignment found: %v\n", sample, fErr)
			skipped = append(skipped, sample)
			continue
		}

		alignments := crams
		if len(alignments) == 0 {
			alignments = bams
		}
		switch len(alignments) {
		case 0:
			color.Red("[%s] no alignment found\n", sample)
			skipped = append(skipped, sample)
		case 1:
			color.Green("[%s] alignment found: %s\n", sample, alignments[0])
			samples = append(samples, SampleWork{
				Sample:  sample,
				Cram:    alignments[0],
				CramDir: filepath.Dir(filepath.Dir(filepath.Dir(alignments[0]))),
			})
		default:
			color.Red("[%s] multiple alignments found — please remove the extra file(s): %v\n", sample, alignments)
			skipped = append(skipped, sample)
		}
	}

	return samples, skipped, nil
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

	color.Green("All file paths valid\n....................................................\n\n")

	// ================================== Chroms and contigs ==================================== //

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return nil, fmt.Errorf("getting chromosomes and contigs: %w", err)
	}

	// =================================== Discover samples ===================================== //

	color.Cyan("================================== Discovering samples =================================\n\n")

	samples, skipped, err := FindGvcfSamples(opts)
	if err != nil {
		return nil, err
	}
	if len(samples) == 0 {
		return nil, fmt.Errorf("no usable samples found (skipped %d)", len(skipped))
	}
	color.Cyan("\nFound %d sample(s), skipped %d\n%s\n\n", len(samples), len(skipped), strings.Repeat("=", 34))

	// ===================================== Build job list ===================================== //

	type job struct {
		sample SampleWork
		chrom  string
		seqs   []SeqInfo
	}

	var jobs []job
	for _, s := range samples {
		for _, c := range chroms {
			jobs = append(jobs, job{sample: s, chrom: c.ID, seqs: []SeqInfo{c}})
		}
		if len(contigs) > 0 {
			jobs = append(jobs, job{sample: s, chrom: "contigs", seqs: contigs})
		}
	}

	// ======================================== Run jobs ======================================== //

	totalCores := runtime.NumCPU()
	maxParallel := totalCores / opts.Threads
	if maxParallel < 1 {
		maxParallel = 1
	}
	color.Cyan("Machine has %d cores. %d threads per job -> max %d parallel jobs\n\n",
		totalCores, opts.Threads, maxParallel)

	jobCh := make(chan job, len(jobs))
	for _, j := range jobs {
		jobCh <- j
	}
	close(jobCh)

	var (
		mu          sync.Mutex
		wg          sync.WaitGroup
		gvcfs       = make(map[string][]string)
		failedTasks []FailedTask
		reused      int
		created     int
	)

	for w := 1; w <= maxParallel; w++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			for j := range jobCh {
				theGVCF := GvcfPath(opts, j.sample, j.chrom)

				if mkErr := os.MkdirAll(filepath.Dir(theGVCF), 0755); mkErr != nil {
					color.Red("[Worker %d] [%s] creating %s: %v\n\n", workerID, j.sample.Sample, filepath.Dir(theGVCF), mkErr)
					mu.Lock()
					failedTasks = append(failedTasks, FailedTask{Sample: j.sample.Sample, Chrom: j.chrom, Reason: mkErr})
					mu.Unlock()
					continue
				}

				// ------------------------- reuse an existing gVCF ------------------------- //
				if _, sErr := os.Stat(theGVCF); sErr == nil {
					if opts.SkipVerification {
						color.Yellow("[Worker %d] [%s] %s exists, skipping integrity check\n\n", workerID, j.sample.Sample, j.chrom)
						mu.Lock()
						gvcfs[j.chrom] = append(gvcfs[j.chrom], theGVCF)
						reused++
						mu.Unlock()
						continue
					}
					if vErr := utils.ValidateGvcf(theGVCF, opts.Verbose, opts.Quick); vErr == nil {
						color.Green("[Worker %d] [%s] gVCF for %s is valid, reusing\n\n", workerID, j.sample.Sample, j.chrom)
						mu.Lock()
						gvcfs[j.chrom] = append(gvcfs[j.chrom], theGVCF)
						reused++
						mu.Unlock()
						continue
					}
					// Remove the partial file and its index so a failed retry cannot
					// leave something that looks complete.
					color.Yellow("[Worker %d] [%s] gVCF for %s is corrupt, re-creating\n\n", workerID, j.sample.Sample, j.chrom)
					os.Remove(theGVCF)
					os.Remove(theGVCF + ".tbi")
				}

				// ------------------------------ create the gVCF ---------------------------- //
				var cErr error
				if caller == "gatk" {
					color.Cyan("[Worker %d] [%s] creating %s with GATK ...\n\n", workerID, j.sample.Sample, j.chrom)
					_, cErr = CreateGvcfGATK(j.sample.Cram, opts.RefFasta, j.seqs, theGVCF, opts.GatkLogLevel, opts.Verbose)
				} else {
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
					color.Cyan("[Worker %d] [%s] creating %s with DeepVariant (%s) ...\n\n", workerID, j.sample.Sample, j.chrom, modelType)
					_, cErr = CreateGvcfDV(j.sample.Cram, opts.RefFasta, j.seqs, theGVCF, opts.DVVer, modelType, opts.Verbose, opts.Threads)
				}

				if cErr != nil {
					color.Red("[Worker %d] [%s] creating gVCF for %s FAILED: %v\n\n", workerID, j.sample.Sample, j.chrom, cErr)
					mu.Lock()
					failedTasks = append(failedTasks, FailedTask{Sample: j.sample.Sample, Chrom: j.chrom, Reason: cErr})
					mu.Unlock()
					continue
				}

				color.Green("[Worker %d] [%s] gVCF for %s created\n\n", workerID, j.sample.Sample, j.chrom)
				mu.Lock()
				gvcfs[j.chrom] = append(gvcfs[j.chrom], theGVCF)
				created++
				mu.Unlock()
			}
		}(w)
	}
	wg.Wait()

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
