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

	// Lower-cased to match FindBams and JointVcfDir. On a case-sensitive
	// filesystem, mixing the two spellings means sample discovery and BAM
	// discovery look in different directories.
	pattern := filepath.Join(opts.DataDir, strings.ToLower(opts.Species), "*", "*", "reference_genomes")
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
				// recreate only selects the verb in the log line below, so one message
				// says whether this is a first run or a repair.
				recreate := false

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
					os.Remove(theGVCF)
					os.Remove(theGVCF + ".tbi")
					recreate = true
				}

				// ------------------------------ create the gVCF ---------------------------- //
				verb := "creating"
				if recreate {
					verb = "re-creating corrupt"
				}

				var cErr error
				if caller == "gatk" {
					color.Cyan("[Worker %d] [%s] %s %s gvcf with GATK ... %s\n\n",
						workerID, j.sample.Sample, verb, j.chrom, theGVCF)
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
					color.Cyan("[Worker %d] [%s] %s %s gvcf with DeepVariant (%s) ... %s\n\n",
						workerID, j.sample.Sample, verb, j.chrom, modelType, theGVCF)
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
