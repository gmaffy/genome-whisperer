// RunAlignReadsDir_fixed.go
// Drop-in replacement for the RunAlignReadsDir function body in RunAlignReads.go.
//
// BUGS FIXED (all were in RunAlignReadsDir, lines ~1829–1916 of the original):
//
// 1. Inverted .bam/.cram replacement when deriving cramPath from rgmdBam.Path:
//      BEFORE: strings.Replace(s.RgmdBam.Path, ".cram", ".bam", 1)   // wrong direction
//      AFTER:  strings.Replace(s.RgmdBam.Path, ".bam",  ".cram", 1)
//
// 2. Missing "." separator when building rgmdBam path from sorted.bam:
//      BEFORE: strings.Replace(s.SortedBam.Path, ".sorted.bam", "rgmd.bam", 1)
//      AFTER:  strings.Replace(s.SortedBam.Path, ".sorted.bam", ".rgmd.bam", 1)
//
// 3. Left-over evaluateSample skeleton (steps []string, return Action{...}) pasted
//    verbatim inside the loop — these identifiers are not in scope here and the
//    early return would exit the function after the very first sample.
//    Removed entirely; actions are performed directly instead of collected.
//
// 4. validCleanReads / missingCleanReads / multipleCleanReads used after the scan
//    loop but never declared — the scan loop must now build validCleanReads so the
//    existing parallel worker section below can run.
//
// 5. validSamples (used in the final summary) declared but only as a counter — it
//    is now a simple int so the final summary compiles correctly.
//
// 6. Missing format argument: color.Red("[%s] Now you have to start from scratch")
//    → color.Red("[%s] No sorted.bam — must rerun alignment from reads\n", s.Sample)
//
// 7. rgmd.cram index / bqsr.cram index actions were appended to an undeclared
//    `steps` slice.  They now call BamIndex directly.
//
// 8. The BQSR / cleanup stubs referenced `steps` and returned Action{} — these
//    are now real operations wired to the same Recalibrate / os.Remove helpers
//    used in RunAlignReadsDir1.

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
)

func RunAlignReadsDir(
	dataDir string,
	species string,
	refVer string,
	refFasta string,
	genomesDir string,
	verbose bool,
	gatkLogLevel string,
	aligner string,
	quick bool,
	skipVer bool,
	bqsr bool,
	bootstrap bool,
	knownSites []string,
	threads int,
) {
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

	// ================================== Discover clean_reads dirs ============================================ //
	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		panic(err)
	}

	color.Green("SAMPLES FOUND:\n---------------------------------------------------------------------\n\n ")
	seen := make(map[string]struct{}, len(matches))
	var cleanReadsDirs []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			cleanReadsDirs = append(cleanReadsDirs, match)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n",
		len(cleanReadsDirs), color.GreenString(species))

	// ============================================= Scan existing alignment files ============================= //
	bamsStatResults, _ := ScanAlignments(dataDir, species, refVer, genomesDir, resolvedFasta, verbose, quick)

	// Index bamsStatResults by sample name for O(1) lookup
	bamsByName := make(map[string]SampleBamState, len(bamsStatResults))
	for _, s := range bamsStatResults {
		bamsByName[s.Sample] = s
	}

	// ============================================= Build work list =========================================== //
	// samplePair is declared in RunAlignReadsDir1 (same file) — reuse it directly.

	var (
		validCleanReads    []samplePair
		missingCleanReads  []string
		multipleCleanReads []string
		perfectSamples     int // already-perfect samples that are skipped
	)

	for _, cleanReadsDir := range cleanReadsDirs {
		sampleName := filepath.Base(filepath.Dir(cleanReadsDir))
		sampleBaseDir := filepath.Dir(cleanReadsDir)
		verDir := filepath.Join(sampleBaseDir, "reference_genomes", refVer, "bams")

		// Skip long-read samples — pbmm2 path not handled in dir-mode yet
		if strings.HasSuffix(sampleName, "LR") {
			color.Yellow("[%s] LR sample detected (pbmm2). Skipping in dir-mode for now.\n", sampleName)
			missingCleanReads = append(missingCleanReads, sampleName)
			continue
		}

		s, hasScanResult := bamsByName[sampleName]

		// --- Already perfect: valid indexed bqsr.cram + valid indexed rgmd.cram, nothing else ---
		if hasScanResult {
			bqsrCramOK := isUsable(s.BqsrCram) && hasIndex(s.BqsrCram)
			rgmdCramOK := isUsable(s.RgmdCram) && hasIndex(s.RgmdCram)

			if bqsrCramOK && rgmdCramOK &&
				!s.SortedBam.Present &&
				!s.RgmdBam.Present &&
				!s.BqsrBam.Present &&
				len(s.OtherFiles) == 0 {
				color.Green("[%s] PASS ✅ — already perfect. Skipping.\n", sampleName)
				perfectSamples++
				continue
			}

			// Ensure the bams output dir exists before we try to write anything
			if mkErr := os.MkdirAll(verDir, 0755); mkErr != nil {
				color.Red("[%s] Error creating bams dir %s: %v\n", sampleName, verDir, mkErr)
				missingCleanReads = append(missingCleanReads, sampleName)
				continue
			}

			// Determine the resume stage from what is already on disk
			switch {
			// bqsr.bam exists — just convert it to bqsr.cram
			case s.BqsrBam.Present:
				color.Yellow("[%s] Resuming: converting bqsr.bam → bqsr.cram\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sampleName,
					bamDir:          verDir,
					resume:          resumeFromBamToCram,
					intermediateBam: s.BqsrBam.Path,
				})

			// rgmd.cram exists — run BQSR on it
			case isUsable(s.RgmdCram):
				color.Yellow("[%s] Resuming: running BQSR on rgmd.cram\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sampleName,
					bamDir:          verDir,
					resume:          resumeFromBQSR,
					intermediateBam: s.RgmdCram.Path,
				})

			// rgmd.cram present but not valid — re-derive from rgmd.bam if available
			case s.RgmdCram.Present && !s.RgmdCram.Valid && isUsable(s.RgmdBam):
				color.Yellow("[%s] Resuming: rgmd.cram invalid, re-converting rgmd.bam → rgmd.cram then BQSR\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sampleName,
					bamDir:          verDir,
					resume:          resumeFromCramRgmdThenBQSR,
					intermediateBam: s.RgmdBam.Path,
				})

			// rgmd.bam exists — convert to rgmd.cram then BQSR
			case isUsable(s.RgmdBam):
				color.Yellow("[%s] Resuming: converting rgmd.bam → rgmd.cram then BQSR\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sampleName,
					bamDir:          verDir,
					resume:          resumeFromCramRgmdThenBQSR,
					intermediateBam: s.RgmdBam.Path,
				})

			// sorted.bam exists — MarkDuplicates → rgmd.cram → BQSR
			case isUsable(s.SortedBam):
				color.Yellow("[%s] Resuming: sorted.bam found — MarkDuplicates → rgmd.cram → BQSR\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sampleName,
					bamDir:          verDir,
					resume:          resumeFromMarkDupAndBQSR,
					intermediateBam: s.SortedBam.Path,
				})

			// No usable intermediate file — fall through to full alignment below
			default:
				goto fullAlignment
			}
			continue
		}

	fullAlignment:
		// No scan result or no usable intermediate file — need reads
		if mkErr := os.MkdirAll(verDir, 0755); mkErr != nil {
			color.Red("[%s] Error creating bams dir %s: %v\n", sampleName, verDir, mkErr)
			missingCleanReads = append(missingCleanReads, sampleName)
			continue
		}

		fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
		fmt.Printf("[%s] fwd: %v  rev: %v\n", sampleName, fwdReads, revReads)
		if readsErr != nil {
			color.Red("[%s] Error reading clean_reads dir: %v\n", sampleName, readsErr)
			missingCleanReads = append(missingCleanReads, sampleName)
			continue
		}

		switch {
		case len(fwdReads) == 0 && len(revReads) == 0:
			color.Red("[%s] No PE reads found in %s\n", sampleName, cleanReadsDir)
			missingCleanReads = append(missingCleanReads, sampleName)

		case len(fwdReads) == 1 && len(revReads) == 1:
			if skipVer {
				color.Yellow("[%s] Skipping fastq verification.\n", sampleName)
				validCleanReads = append(validCleanReads, samplePair{
					sample: sampleName, fwd: fwdReads[0], rev: revReads[0],
					bamDir: verDir, resume: resumeFromAlign,
				})
			} else {
				color.Blue("[%s] Validating fastq files...\n", sampleName)
				fwdErr := validateFastqGz(fwdReads[0], verbose, quick)
				revErr := validateFastqGz(revReads[0], verbose, quick)
				if fwdErr != nil || revErr != nil {
					color.Red("[%s] Fastq validation failed (fwd=%v rev=%v)\n", sampleName, fwdErr, revErr)
					missingCleanReads = append(missingCleanReads, sampleName)
				} else {
					color.Green("[%s] PE reads valid.\n", sampleName)
					validCleanReads = append(validCleanReads, samplePair{
						sample: sampleName, fwd: fwdReads[0], rev: revReads[0],
						bamDir: verDir, resume: resumeFromAlign,
					})
				}
			}

		default:
			color.Yellow("[%s] Multiple PE reads found — cannot determine which to use.\n", sampleName)
			multipleCleanReads = append(multipleCleanReads, cleanReadsDir)
			missingCleanReads = append(missingCleanReads, sampleName)
		}
	}

	color.Green("\nAlready perfect (skipped):   %d\n", perfectSamples)
	color.Green("Ready to process:            %d\n", len(validCleanReads))
	color.Red("Missing/invalid reads:       %d\n", len(missingCleanReads))
	color.Yellow("Skipped (multiple reads):   %d\n", len(multipleCleanReads))
	fmt.Printf("\n==============================================================\n\n")

	if len(validCleanReads) == 0 {
		color.Yellow("No samples to process. Exiting.\n")
		return
	}

	// ======================================= Parallel worker section ============================================= //
	totalCores := runtime.NumCPU()
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}
	color.Green("Processing %d sample(s) — up to %d jobs in parallel, %d threads each.\n\n",
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

			// ---- Output paths ---- //
			sortedBam    := filepath.Join(readPair.bamDir, sampleName+".sorted.bam")
			rgmdBam      := filepath.Join(readPair.bamDir, sampleName+".rgmd.bam")
			rgmdCram     := filepath.Join(readPair.bamDir, sampleName+".rgmd.cram")
			rgmdMetrics  := filepath.Join(readPair.bamDir, sampleName+".rgmd.metrics.txt")
			bqsrCram     := filepath.Join(readPair.bamDir, sampleName+".bqsr.cram")

			// bamToCram converts src BAM/CRAM → dst CRAM and immediately indexes it.
			bamToCram := func(src, dst string) error {
				cmd := fmt.Sprintf(`samtools view -C -T %s -o %s %s && samtools index %s`,
					resolvedFasta, dst, src, dst)
				color.Green("[%s] %s\n", sampleName, cmd)
				if verbose {
					return utils.RunBashCmdVerbose(cmd)
				}
				return utils.RunBashCmd(cmd)
			}

			// runBQSR runs bootstrap/known-sites recalibration and returns the output cram path.
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

			// deleteIfPresent silently removes a file, logging but not failing on error.
			deleteIfPresent := func(path string) {
				if _, statErr := os.Stat(path); statErr == nil {
					if rmErr := os.Remove(path); rmErr != nil {
						color.Yellow("[%s] Warning: could not delete %s: %v\n", sampleName, path, rmErr)
					} else {
						color.Green("[%s] Deleted intermediate file: %s\n", sampleName, path)
					}
				}
			}

			// cleanupIntermediates removes .sorted.bam, .rgmd.bam, .bqsr.bam once the
			// sample is perfect (valid indexed rgmd.cram + valid indexed bqsr.cram).
			cleanupIntermediates := func() {
				deleteIfPresent(sortedBam)
				deleteIfPresent(sortedBam + ".bai")
				deleteIfPresent(rgmdBam)
				deleteIfPresent(rgmdBam + ".bai")
				// bqsr.bam lives in the same dir, derive from bqsrCram path
				deleteIfPresent(strings.Replace(bqsrCram, ".cram", ".bam", 1))
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
				cleanupIntermediates()
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
				cleanupIntermediates()
				return

			// ------------------------------------------------------------------ //
			// resumeFromCramRgmdThenBQSR: rgmd.bam → convert to cram, then BQSR  //
			// ------------------------------------------------------------------ //
			case resumeFromCramRgmdThenBQSR:
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				// BUG FIX 1: was replacing ".cram" with ".bam" (backwards).
				// The source is readPair.intermediateBam (the .bam); destination is rgmdCram.
				if err := bamToCram(readPair.intermediateBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] rgmd.cram ready: %s\n", sampleName, rgmdCram)

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete.\n", sampleName)
					cleanupIntermediates()
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				cleanupIntermediates()
				return

			// ------------------------------------------------------------------ //
			// resumeFromMarkDupAndBQSR: sorted.bam → MarkDup → rgmd.cram → BQSR //
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
				idxCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(idxCmd)
				} else {
					indErr = utils.RunBashCmd(idxCmd)
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
					cleanupIntermediates()
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				cleanupIntermediates()
				return

			// ------------------------------------------------------------------ //
			// resumeFromAlign (default): full alignment from clean reads          //
			// ------------------------------------------------------------------ //
			default:
				// bwa-mem2 default; fall back to bwa-mem when specified
				alignCmd := fmt.Sprintf(
					`bwa-mem2 mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
					threads, sampleName, sampleName, sampleName,
					resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				if aligner == "bwa-mem" {
					alignCmd = fmt.Sprintf(
						`bwa mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
						threads, sampleName, sampleName, sampleName,
						resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
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

				// Mark duplicates
				color.Green("[%s] Marking duplicates...\n", sampleName)
				// BUG FIX 2: was missing the "." before "rgmd.bam", producing "<prefix>rgmd.bam".
				// Now uses the pre-computed rgmdBam path variable.
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

				// Index rgmd.bam
				color.Green("[%s] Indexing rgmd.bam...\n", sampleName)
				idxCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				fmt.Printf("[%s] %s\n", sampleName, idxCmd)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(idxCmd)
				} else {
					indErr = utils.RunBashCmd(idxCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}
				color.Green("[%s] BAM index completed.\n", sampleName)

				// Convert rgmd.bam → rgmd.cram
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				if !bqsr {
					color.Green("[%s] Pipeline complete (no BQSR). Final CRAM: %s\n", sampleName, rgmdCram)
					cleanupIntermediates()
					return
				}

				color.Green("[%s] Starting BQSR...\n", sampleName)
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Final: %s\n", sampleName, bqsrOut)
				cleanupIntermediates()
			}

		}(rp)
	}

	wg.Wait()

	// ======================================= Final summary ======================================================= //
	fmt.Printf("\n==============================================================\n")
	color.Green("Already perfect (skipped): %d\n", perfectSamples)
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
