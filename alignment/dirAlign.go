package alignment

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
)

type sampleWork struct {
	sample  string
	cram    string
	cramDir string // filepath.Dir of the cram, used to derive gvcf output path
}

// resumeStage records where in the pipeline a sample should resume.
type resumeStage int

const (
	resumeFromAlign            resumeStage = iota // no intermediate files found – start from scratch
	resumeFromMarkDupAndBQSR                      // sorted.bam exists – run MarkDuplicates → BQSR → convert
	resumeFromBQSR                                // rgmd.cram exists – run BQSR directly
	resumeFromCramRgmdThenBQSR                    // rgmd.bam exists – convert to cram, then BQSR
	resumeFromBamToCram                           // bqsr.bam exists – convert BAM→CRAM only
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

	return fwdReads, revReads, nil
}

func validateFastqGz(fastq string, verbose bool, quick bool) error {
	if quick {
		// Just check the gzip magic bytes and EOF integrity
		valStr := fmt.Sprintf("gzip -t %s", fastq)
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
		if verbose {
			return utils.RunBashCmdVerbose(valStr)
		}
		return utils.RunBashCmd(valStr)
	}

	// Full validation: gzip integrity + FASTQ format check (4-line records, quality scores)
	valStr := fmt.Sprintf(
		`bash -c 'gzip -cd %s | awk "NR%%4==1 && !/^@/ { print \"Bad header at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } NR%%4==3 && !/^\+/ { print \"Bad separator at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } END { if(NR%%4!=0) { print \"Truncated: \", NR, \"lines\" > \"/dev/stderr\"; exit 1 } }" '`,
		fastq,
	)
	fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
	if verbose {
		return utils.RunBashCmdVerbose(valStr)
	}
	return utils.RunBashCmd(valStr)
}

// GenomeRef holds a discovered reference genome on disk.
type GenomeRef struct {
	RefVer    string // e.g. "GRCh38"
	FastaPath string // path to .fa / .fasta / .fna
	DictPath  string // path to .dict
}

// GetValidGenomesFromDisk scans genomesDir for species → reference assemblies.
// Expected layout: <genomesDir>/<species>/<refVer>/assembly/  containing a fasta + dict file.
// Returns map[SPECIES_UPPER] → []GenomeRef.
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

			var fastaPath, dictPath string
			for _, f := range entries {
				if f.IsDir() {
					continue
				}
				name := f.Name()
				switch {
				case strings.HasSuffix(name, ".fa") ||
					strings.HasSuffix(name, ".fasta") ||
					strings.HasSuffix(name, ".fna"):
					fastaPath = filepath.Join(assemblyDir, name)
				case strings.HasSuffix(name, ".dict"):
					dictPath = filepath.Join(assemblyDir, name)
				}
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

	return result, nil
}

// resolveRefFasta returns refFasta as-is when provided, otherwise auto-discovers
// it from genomesDir using species and refVer. Returns the fasta path or an error.
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

func RunAlignReadsDir1(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
	// ============================================= Validate paths ================================================ //
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		fmt.Printf("Error getting absolute path for data directory: %s\n", dataDir)
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

	// ---- Resolve reference fasta (explicit path or auto-discover) ---- //
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

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	// ---- Validate BQSR parameters up-front ---- //
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

	color.Green("All file paths valid\n....................................................\n\n")
	fmt.Println("data dir abs:", dataDirAbs)

	// ================================== Discover samples ===================================== //
	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		panic(err)
	}

	color.Green("SAMPLES FOUND:\n---------------------------------------------------------------------\n\n ")
	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, match)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n", len(samples), color.GreenString(species))

	// ======================================= Resolve valid samples up-front ======================================= //
	fmt.Println("Looking for existing bam/cram files ....")

	var (
		validSamples []sampleWork
		missingBams  []string
		multipleBams []string
	)

	for _, sample := range samples {
		sampleName := filepath.Base(filepath.Dir(sample))
		sampleBaseDir := filepath.Dir(sample)
		verDir := filepath.Join(sampleBaseDir, "reference_genomes", refVer, "bams")

		// Helper: check if a file exists in verDir.
		fileExists := func(name string) (string, bool) {
			p := filepath.Join(verDir, name)
			if _, err := os.Stat(p); err == nil {
				return p, true
			}
			return "", false
		}

		bqsrCram := sampleName + ".bqsr.cram"
		bqsrBam := sampleName + ".bqsr.bam"
		rgmdCram := sampleName + ".rgmd.cram"
		rgmdBamFile := sampleName + ".rgmd.bam"
		sortedBamFile := sampleName + ".sorted.bam"

		// ---- Stage 1: bqsr.cram already present ---- //
		if cramPath, ok := fileExists(bqsrCram); ok {
			color.Green("[%s] bqsr.cram FOUND 😊: %s\n", sampleName, cramPath)
			if err := variants.ValidateBam(cramPath, resolvedFasta, verbose, quick); err != nil {
				color.Red("[%s] bqsr.cram not valid: %v — will re-process\n", sampleName, err)
				missingBams = append(missingBams, sample)
			} else {
				validSamples = append(validSamples, sampleWork{
					sample:  sample,
					cram:    cramPath,
					cramDir: filepath.Dir(filepath.Dir(filepath.Dir(cramPath))),
				})
			}
			continue
		}

		// ---- Stage 2: bqsr.bam present → convert BAM→CRAM only ---- //
		if bamPath, ok := fileExists(bqsrBam); ok {
			color.Yellow("[%s] bqsr.bam found – will convert BAM→CRAM only.\n", sampleName)
			missingBams = append(missingBams, sample)
			// Store intermediate path; the goroutine will pick it up via samplePair.
			// We encode this by appending a special marker; handled below in the work list.
			_ = bamPath // used when building samplePair further down
			continue
		}

		// ---- Stage 3: rgmd.cram present → run BQSR on it ---- //
		if cramPath, ok := fileExists(rgmdCram); ok {
			color.Yellow("[%s] rgmd.cram found – will run BQSR.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = cramPath
			continue
		}

		// ---- Stage 4: rgmd.bam present → convert to CRAM, then BQSR ---- //
		if bamPath, ok := fileExists(rgmdBamFile); ok {
			color.Yellow("[%s] rgmd.bam found – will convert to CRAM then run BQSR.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = bamPath
			continue
		}

		// ---- Stage 5: sorted.bam present → run MarkDuplicates + BQSR + convert ---- //
		if bamPath, ok := fileExists(sortedBamFile); ok {
			color.Yellow("[%s] sorted.bam found – will run MarkDuplicates → BQSR → convert.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = bamPath
			continue
		}

		// ---- No intermediate files found – needs full alignment ---- //
		color.Red("[%s] No bqsr.cram or intermediate files found 😒 – needs full alignment.\n", sampleName)
		missingBams = append(missingBams, sample)
	}

	color.Green("\nAlready aligned (valid): %d\n", len(validSamples))
	color.Red("Need alignment (missing): %d\n", len(missingBams))
	color.Yellow("Skipped (multiple crам files): %d\n", len(multipleBams))
	fmt.Printf("\n==============================================================\n\n")

	// ======================================= Build work list for missing samples ================================== //
	type samplePair struct {
		sample          string
		fwd             string
		rev             string
		bamDir          string
		resume          resumeStage // where to start in the pipeline
		intermediateBam string      // pre-existing BAM/CRAM to resume from (if any)
	}

	var missingCleanReads []string
	var multipleCleanReads []string
	var validCleanReads []samplePair

	if len(missingBams) > 0 {
		color.Green("Preparing to align %d sample(s):\n---------------------------------------------------------------------\n\n", len(missingBams))
		for _, cleanReadsDir := range missingBams {
			sample := filepath.Base(filepath.Dir(cleanReadsDir))
			fmt.Println(sample)

			if strings.HasSuffix(sample, "LR") {
				// Long-read samples: pbmm2 path — not handled in this function yet
				color.Yellow("[%s] LR sample detected (pbmm2). Skipping in dir-mode for now.\n", sample)
				missingCleanReads = append(missingCleanReads, sample)
				continue
			}

			// Ensure output bams dir exists: <sampleDir>/reference_genomes/<refVer>/bams/
			sampleBaseDir := filepath.Dir(cleanReadsDir)
			verDir := filepath.Join(sampleBaseDir, "reference_genomes", refVer, "bams")
			if mkErr := os.MkdirAll(verDir, 0755); mkErr != nil {
				color.Red("[%s] Error creating bams dir %s: %v\n", sample, verDir, mkErr)
				missingCleanReads = append(missingCleanReads, sample)
				continue
			}

			// ---- Determine resume stage by re-checking intermediate files ---- //
			fileInVerDir := func(name string) (string, bool) {
				p := filepath.Join(verDir, name)
				if _, err := os.Stat(p); err == nil {
					return p, true
				}
				return "", false
			}

			bqsrBamPath, hasBqsrBam := fileInVerDir(sample + ".bqsr.bam")
			rgmdCramPath, hasRgmdCram := fileInVerDir(sample + ".rgmd.cram")
			rgmdBamPath, hasRgmdBam := fileInVerDir(sample + ".rgmd.bam")
			sortedBamPath, hasSortedBam := fileInVerDir(sample + ".sorted.bam")

			switch {
			case hasBqsrBam:
				color.Yellow("[%s] Resuming: converting bqsr.bam → bqsr.cram\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromBamToCram,
					intermediateBam: bqsrBamPath,
				})

			case hasRgmdCram:
				color.Yellow("[%s] Resuming: running BQSR on rgmd.cram\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromBQSR,
					intermediateBam: rgmdCramPath,
				})

			case hasRgmdBam:
				color.Yellow("[%s] Resuming: converting rgmd.bam → cram, then BQSR\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromCramRgmdThenBQSR,
					intermediateBam: rgmdBamPath,
				})

			case hasSortedBam:
				color.Yellow("[%s] Resuming: sorted.bam found – MarkDuplicates → BQSR → convert\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromMarkDupAndBQSR,
					intermediateBam: sortedBamPath,
				})

			default:
				// No intermediate files – need reads for full alignment
				fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
				fmt.Printf("[%s] fwd: %v  rev: %v\n", sample, fwdReads, revReads)
				if readsErr != nil {
					color.Red("[%s] Error reading clean_reads dir: %v\n", sample, readsErr)
					missingCleanReads = append(missingCleanReads, sample)
					continue
				}

				switch {
				case len(fwdReads) == 0 && len(revReads) == 0:
					color.Red("[%s] No PE reads found in %s\n", sample, cleanReadsDir)
					missingCleanReads = append(missingCleanReads, sample)
				case len(fwdReads) == 1 && len(revReads) == 1:
					if skipVer {
						color.Yellow("[%s] Skipping fastq verification.\n", sample)
						validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir, resume: resumeFromAlign})
					} else {
						color.Blue("[%s] Validating fastq files...\n", sample)
						fwdErr := validateFastqGz(fwdReads[0], verbose, quick)
						revErr := validateFastqGz(revReads[0], verbose, quick)
						if fwdErr != nil || revErr != nil {
							color.Red("[%s] Fastq validation failed (fwd=%v rev=%v)\n", sample, fwdErr, revErr)
							missingCleanReads = append(missingCleanReads, sample)
						} else {
							color.Green("[%s] PE reads valid.\n", sample)
							validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir, resume: resumeFromAlign})
						}
					}
				default:
					color.Yellow("[%s] Multiple PE reads found — cannot determine which to use.\n", sample)
					multipleCleanReads = append(multipleCleanReads, cleanReadsDir)
					missingCleanReads = append(missingCleanReads, sample)
				}
			}
		}
	}

	color.Green("\nReady to align: %d\n", len(validCleanReads))
	color.Red("Missing/invalid reads: %d\n", len(missingCleanReads))
	color.Yellow("Skipped (multiple read files): %d\n", len(multipleCleanReads))
	fmt.Printf("\n==============================================================\n\n")

	if len(validCleanReads) == 0 {
		color.Yellow("No samples to align. Exiting.\n")
		return
	}

	// ======================================= Parallel alignment ================================================== //
	totalCores := runtime.NumCPU()
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}
	color.Green("Aligning %d sample(s) — up to %d jobs in parallel, %d threads each.\n\n",
		len(validCleanReads), maxParallelJobs, threads)

	var (
		mu               sync.Mutex
		failedAlignments []string
		wg               sync.WaitGroup
		sem              = make(chan struct{}, maxParallelJobs)
	)

	for _, rp := range validCleanReads {
		wg.Add(1)
		sem <- struct{}{}

		go func(readPair samplePair) {
			defer wg.Done()
			defer func() { <-sem }()

			sampleName := readPair.sample

			addFailure := func(reason string) {
				color.Red("[%s] FAILED: %s\n", sampleName, reason)
				mu.Lock()
				failedAlignments = append(failedAlignments, sampleName)
				mu.Unlock()
			}

			// ---- Paths ---- //
			sortedBam := filepath.Join(readPair.bamDir, sampleName+".sorted.bam")
			rgmdBam := filepath.Join(readPair.bamDir, sampleName+".rgmd.bam")
			rgmdCram := filepath.Join(readPair.bamDir, sampleName+".rgmd.cram")
			rgmdMetrics := filepath.Join(readPair.bamDir, sampleName+".rgmd.metrics.txt")
			bqsrCram := filepath.Join(readPair.bamDir, sampleName+".bqsr.cram")

			// bamToCram converts src BAM/CRAM to dst CRAM using samtools.
			bamToCram := func(src, dst string) error {
				cmd := fmt.Sprintf(`samtools view -C -T %s -o %s %s && samtools index %s`, resolvedFasta, dst, src, dst)
				color.Green("[%s] Converting to CRAM: %s\n", sampleName, cmd)
				if verbose {
					return utils.RunBashCmdVerbose(cmd)
				}
				return utils.RunBashCmd(cmd)
			}

			// runBQSR runs bootstrap/known-sites recalibration on inputBam and returns the output path.
			runBQSR := func(inputBam string) (string, error) {
				sampleKnownSites := make([]string, len(knownSites))
				copy(sampleKnownSites, knownSites)

				if len(sampleKnownSites) == 0 && bootstrap {
					color.Green("[%s] Creating known variants (bootstrap)...\n", sampleName)
					newKS, ksErr := CreateKnownVariants(resolvedFasta, inputBam, "", gatkLogLevel, verbose)
					if ksErr != nil {
						return "", fmt.Errorf("CreateKnownVariants: %v", ksErr)
					}
					sampleKnownSites = newKS
					color.Green("[%s] Known variants created: %v\n", sampleName, sampleKnownSites)
				}

				color.Green("[%s] Running Recalibrate...\n", sampleName)
				out, recErr := Recalibrate(resolvedFasta, inputBam, sampleKnownSites, "", verbose)
				if recErr != nil {
					return "", fmt.Errorf("Recalibrate: %v", recErr)
				}
				return out, nil
			}

			switch readPair.resume {

			// ------------------------------------------------------------------ //
			// resumeFromBamToCram: bqsr.bam exists → convert to bqsr.cram only   //
			// ------------------------------------------------------------------ //
			case resumeFromBamToCram:
				color.Green("[%s] Converting bqsr.bam → bqsr.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, bqsrCram); err != nil {
					addFailure(fmt.Sprintf("bqsr BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] Done. Final CRAM: %s\n", sampleName, bqsrCram)
				return

			// ------------------------------------------------------------------ //
			// resumeFromBQSR: rgmd.cram exists → run BQSR directly              //
			// ------------------------------------------------------------------ //
			case resumeFromBQSR:
				if !bqsr {
					color.Green("[%s] BQSR disabled — nothing to do (rgmd.cram already present).\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(readPair.intermediateBam)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromCramRgmdThenBQSR: rgmd.bam → cram it, then BQSR         //
			// ------------------------------------------------------------------ //
			case resumeFromCramRgmdThenBQSR:
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] rgmd.cram ready: %s\n", sampleName, rgmdCram)

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete.\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromMarkDupAndBQSR: sorted.bam → MarkDuplicates → BQSR → CRAM //
			// ------------------------------------------------------------------ //
			case resumeFromMarkDupAndBQSR:
				color.Green("[%s] Marking duplicates on sorted.bam...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, readPair.intermediateBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates complete.\n", sampleName)

				color.Green("[%s] Indexing rgmd.bam...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}

				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete. Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromAlign (default): full alignment from clean reads          //
			// ------------------------------------------------------------------ //
			default:
				// ----------------------------------------- BWA-MEM2 / BWA-MEM ----------------------------------------- //
				alignCmd := fmt.Sprintf(
					`bwa-mem2 mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
					threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				if aligner == "bwa-mem" {
					alignCmd = fmt.Sprintf(
						`bwa mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
						threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				}

				color.Green("[%s] Aligning with %s...\n", sampleName, aligner)
				fmt.Printf("[%s] %s\n", sampleName, alignCmd)
				var memErr error
				if verbose {
					memErr = utils.RunBashCmdVerbose(alignCmd)
				} else {
					memErr = utils.RunBashCmd(alignCmd)
				}
				if memErr != nil {
					addFailure(fmt.Sprintf("alignment: %v", memErr))
					return
				}
				color.Green("[%s] Alignment completed.\n", sampleName)

				// ----------------------------------------- Mark Duplicates -------------------------------------------- //
				color.Green("[%s] Marking duplicates...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates completed.\n", sampleName)

				// ----------------------------------------- BAM Index -------------------------------------------------- //
				color.Green("[%s] Indexing BAM...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				fmt.Printf("[%s] %s\n", sampleName, indexCmd)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}
				color.Green("[%s] BAM index completed.\n", sampleName)

				// ----------------------------------------- Convert to CRAM -------------------------------------------- //
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				// ============================================= BQSR ==================================================== //
				if !bqsr {
					color.Green("[%s] Pipeline complete (no BQSR). Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}

				color.Green("[%s] Starting BQSR...\n", sampleName)
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Final: %s\n", sampleName, bqsrOut)
			}

		}(rp)
	}

	wg.Wait()

	// ======================================= Final summary ======================================================= //
	fmt.Printf("\n==============================================================\n")
	color.Green("Already aligned (skipped): %d\n", len(validSamples))
	color.Green("Successfully processed:    %d\n", len(validCleanReads)-len(failedAlignments))
	if len(failedAlignments) > 0 {
		color.Red("Failed:                    %d\n", len(failedAlignments))
		for _, s := range failedAlignments {
			color.Red("  - %s\n", s)
		}
	}
	if len(missingCleanReads) > 0 {
		color.Yellow("Skipped (no valid reads):  %d\n", len(missingCleanReads))
		for _, s := range missingCleanReads {
			color.Yellow("  - %s\n", s)
		}
	}
	fmt.Printf("==============================================================\n\n")
}

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
	// ============================================= Validate paths ================================================ //
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		fmt.Printf("Error getting absolute path for data directory: %s\n", dataDir)
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

	// ---- Resolve reference fasta (explicit path or auto-discover) ---- //
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

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	// ---- Validate BQSR parameters up-front ---- //
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

	color.Green("All file paths valid\n....................................................\n\n")
	fmt.Println("data dir abs:", dataDirAbs)

	// ================================== Discover samples ===================================== //
	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		panic(err)
	}

	color.Green("SAMPLES FOUND:\n---------------------------------------------------------------------\n\n ")
	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, match)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n", len(samples), color.GreenString(species))

	// ============================================= Scan alignment files =========================================== //

	bamsStatResults, _ := ScanAlignments(dataDir, species, refVer, genomesDir, resolvedFasta, verbose, quick)

	for _, s := range bamsStatResults {
		bqsrCramOK := isUsable(s.BqsrCram) && hasIndex(s.BqsrCram)
		rgmdCramOK := isUsable(s.RgmdCram) && hasIndex(s.RgmdCram)

		// --- If already perfect ---
		if bqsrCramOK && rgmdCramOK &&
			!s.SortedBam.Present &&
			!s.RgmdBam.Present &&
			!s.BqsrBam.Present &&
			len(s.OtherFiles) == 0 {

			color.Green("[%s] PASS ✅ - Skipping'\n", s.Sample)
		} else if !rgmdCramOK {
			if isUsable(s.RgmdBam) {
				color.Green("[%s] Converting rgmd.bam → rgmd.cram\n", s.Sample)
				cramPath := strings.Replace(s.RgmdBam.Path, ".cram", ".bam", 1)
				err = BamToCram(s.RgmdBam.Path, cramPath, resolvedFasta, verbose)

			} else if isUsable(s.SortedBam) {
				//steps = append(steps, "Run markdup: sorted.bam → rgmd.bam")
				rgmdBam := strings.Replace(s.SortedBam.Path, ".sorted.bam", "rgmd.bam", 1)
				rgmdMetrics := strings.Replace(s.SortedBam.Path, ".sorted.bam", "rgmd.metrics.txt", 1)
				err = MarkDuplicates(resolvedFasta, s.SortedBam.Path, rgmdBam, rgmdMetrics, verbose, aligner, gatkLogLevel)
				if err != nil {
					color.Red("[%s] MarkDuplicates failed ❌\n", s.Sample)
				} else {
					rgmdCram := strings.Replace(rgmdBam, ".bam", ".cram", 1)
					err2 := BamToCram(rgmdBam, rgmdCram, resolvedFasta, verbose)
					if err2 != nil {
						color.Red("[%s] failed to convert rgmdBam to cram ❌\n", s.Sample)
					} else {
						color.Green("[%s] RGMD bam converted to cram ✅\n", s.Sample)
					}
				}
				//steps = append(steps, "Convert rgmd.bam → rgmd.cram")
			} else {
				//steps = append(steps, "Missing sorted.bam → must rerun alignment (bwa-mem)")
				color.Red("[%s] Now you have to start from scratch", s.Sample)
			}
		}

		// --- Ensure rgmd.cram index ---
		if isUsable(s.RgmdCram) && !hasIndex(s.RgmdCram) {
			//steps = append(steps, "Index rgmd.cram")
			err = BamIndex(s.RgmdBam.Path, verbose)
			if err != nil {
				color.Red("[%s] failed to index RGMD cram ❌\n", s.Sample)
			}
		}

		// --- Ensure bqsr.cram ---
		if !bqsrCramOK {
			if isUsable(s.RgmdCram) || isUsable(s.RgmdBam) {
				steps = append(steps, "Run BQSR → bqsr.bam")
				steps = append(steps, "Convert bqsr.bam → bqsr.cram")
			} else {
				steps = append(steps, "Cannot run BQSR: missing rgmd input")
			}
		}

		// --- Ensure bqsr.cram index ---
		if isUsable(s.BqsrCram) && !hasIndex(s.BqsrCram) {
			steps = append(steps, "Index bqsr.cram")
		}

		// --- Cleanup phase ---
		if bqsrCramOK && rgmdCramOK {
			if s.SortedBam.Present {
				steps = append(steps, "Delete sorted.bam")
			}
			if s.RgmdBam.Present {
				steps = append(steps, "Delete rgmd.bam")
			}
			if s.BqsrBam.Present {
				steps = append(steps, "Delete bqsr.bam")
			}
			if len(s.OtherFiles) > 0 {
				steps = append(steps, "Remove unexpected extra files")
			}
		}

		return Action{
			Sample:  s.Sample,
			Steps:   steps,
			Perfect: false,
		}

	}

	color.Green("\nReady to align: %d\n", len(validCleanReads))
	color.Red("Missing/invalid reads: %d\n", len(missingCleanReads))
	color.Yellow("Skipped (multiple read files): %d\n", len(multipleCleanReads))
	fmt.Printf("\n==============================================================\n\n")

	if len(validCleanReads) == 0 {
		color.Yellow("No samples to align. Exiting.\n")
		return
	}

	// ======================================= Parallel alignment ================================================== //
	totalCores := runtime.NumCPU()
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}
	color.Green("Aligning %d sample(s) — up to %d jobs in parallel, %d threads each.\n\n",
		len(validCleanReads), maxParallelJobs, threads)

	var (
		mu               sync.Mutex
		failedAlignments []string
		wg               sync.WaitGroup
		sem              = make(chan struct{}, maxParallelJobs)
	)

	for _, rp := range validCleanReads {
		wg.Add(1)
		sem <- struct{}{}

		go func(readPair samplePair) {
			defer wg.Done()
			defer func() { <-sem }()

			sampleName := readPair.sample

			addFailure := func(reason string) {
				color.Red("[%s] FAILED: %s\n", sampleName, reason)
				mu.Lock()
				failedAlignments = append(failedAlignments, sampleName)
				mu.Unlock()
			}

			// ---- Paths ---- //
			sortedBam := filepath.Join(readPair.bamDir, sampleName+".sorted.bam")
			rgmdBam := filepath.Join(readPair.bamDir, sampleName+".rgmd.bam")
			rgmdCram := filepath.Join(readPair.bamDir, sampleName+".rgmd.cram")
			rgmdMetrics := filepath.Join(readPair.bamDir, sampleName+".rgmd.metrics.txt")
			bqsrCram := filepath.Join(readPair.bamDir, sampleName+".bqsr.cram")

			// bamToCram converts src BAM/CRAM to dst CRAM using samtools.
			bamToCram := func(src, dst string) error {
				cmd := fmt.Sprintf(`samtools view -C -T %s -o %s %s && samtools index %s`, resolvedFasta, dst, src, dst)
				color.Green("[%s] Converting to CRAM: %s\n", sampleName, cmd)
				if verbose {
					return utils.RunBashCmdVerbose(cmd)
				}
				return utils.RunBashCmd(cmd)
			}

			// runBQSR runs bootstrap/known-sites recalibration on inputBam and returns the output path.
			runBQSR := func(inputBam string) (string, error) {
				sampleKnownSites := make([]string, len(knownSites))
				copy(sampleKnownSites, knownSites)

				if len(sampleKnownSites) == 0 && bootstrap {
					color.Green("[%s] Creating known variants (bootstrap)...\n", sampleName)
					newKS, ksErr := CreateKnownVariants(resolvedFasta, inputBam, "", gatkLogLevel, verbose)
					if ksErr != nil {
						return "", fmt.Errorf("CreateKnownVariants: %v", ksErr)
					}
					sampleKnownSites = newKS
					color.Green("[%s] Known variants created: %v\n", sampleName, sampleKnownSites)
				}

				color.Green("[%s] Running Recalibrate...\n", sampleName)
				out, recErr := Recalibrate(resolvedFasta, inputBam, sampleKnownSites, "", verbose)
				if recErr != nil {
					return "", fmt.Errorf("Recalibrate: %v", recErr)
				}
				return out, nil
			}

			switch readPair.resume {

			// ------------------------------------------------------------------ //
			// resumeFromBamToCram: bqsr.bam exists → convert to bqsr.cram only   //
			// ------------------------------------------------------------------ //
			case resumeFromBamToCram:
				color.Green("[%s] Converting bqsr.bam → bqsr.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, bqsrCram); err != nil {
					addFailure(fmt.Sprintf("bqsr BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] Done. Final CRAM: %s\n", sampleName, bqsrCram)
				return

			// ------------------------------------------------------------------ //
			// resumeFromBQSR: rgmd.cram exists → run BQSR directly              //
			// ------------------------------------------------------------------ //
			case resumeFromBQSR:
				if !bqsr {
					color.Green("[%s] BQSR disabled — nothing to do (rgmd.cram already present).\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(readPair.intermediateBam)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromCramRgmdThenBQSR: rgmd.bam → cram it, then BQSR         //
			// ------------------------------------------------------------------ //
			case resumeFromCramRgmdThenBQSR:
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] rgmd.cram ready: %s\n", sampleName, rgmdCram)

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete.\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromMarkDupAndBQSR: sorted.bam → MarkDuplicates → BQSR → CRAM //
			// ------------------------------------------------------------------ //
			case resumeFromMarkDupAndBQSR:
				color.Green("[%s] Marking duplicates on sorted.bam...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, readPair.intermediateBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates complete.\n", sampleName)

				color.Green("[%s] Indexing rgmd.bam...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}

				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete. Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromAlign (default): full alignment from clean reads          //
			// ------------------------------------------------------------------ //
			default:
				// ----------------------------------------- BWA-MEM2 / BWA-MEM ----------------------------------------- //
				alignCmd := fmt.Sprintf(
					`bwa-mem2 mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
					threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				if aligner == "bwa-mem" {
					alignCmd = fmt.Sprintf(
						`bwa mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
						threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				}

				color.Green("[%s] Aligning with %s...\n", sampleName, aligner)
				fmt.Printf("[%s] %s\n", sampleName, alignCmd)
				var memErr error
				if verbose {
					memErr = utils.RunBashCmdVerbose(alignCmd)
				} else {
					memErr = utils.RunBashCmd(alignCmd)
				}
				if memErr != nil {
					addFailure(fmt.Sprintf("alignment: %v", memErr))
					return
				}
				color.Green("[%s] Alignment completed.\n", sampleName)

				// ----------------------------------------- Mark Duplicates -------------------------------------------- //
				color.Green("[%s] Marking duplicates...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates completed.\n", sampleName)

				// ----------------------------------------- BAM Index -------------------------------------------------- //
				color.Green("[%s] Indexing BAM...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				fmt.Printf("[%s] %s\n", sampleName, indexCmd)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}
				color.Green("[%s] BAM index completed.\n", sampleName)

				// ----------------------------------------- Convert to CRAM -------------------------------------------- //
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				// ============================================= BQSR ==================================================== //
				if !bqsr {
					color.Green("[%s] Pipeline complete (no BQSR). Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}

				color.Green("[%s] Starting BQSR...\n", sampleName)
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Final: %s\n", sampleName, bqsrOut)
			}

		}(rp)
	}

	wg.Wait()

	// ======================================= Final summary ======================================================= //
	fmt.Printf("\n==============================================================\n")
	color.Green("Already aligned (skipped): %d\n", len(validSamples))
	color.Green("Successfully processed:    %d\n", len(validCleanReads)-len(failedAlignments))
	if len(failedAlignments) > 0 {
		color.Red("Failed:                    %d\n", len(failedAlignments))
		for _, s := range failedAlignments {
			color.Red("  - %s\n", s)
		}
	}
	if len(missingCleanReads) > 0 {
		color.Yellow("Skipped (no valid reads):  %d\n", len(missingCleanReads))
		for _, s := range missingCleanReads {
			color.Yellow("  - %s\n", s)
		}
	}
	fmt.Printf("==============================================================\n\n")
}
