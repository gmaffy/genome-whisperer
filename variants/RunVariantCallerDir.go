package variants

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

func runCmd(cmdStr string, verbose bool) error {
	return utils.RunCmd(cmdStr, verbose)
}

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

func FindBamOrVcfs(dataDirAbs string, species string, sample string, refVer string, bamOrVcf string, chrom string) ([]string, error) {
	var pattern string
	switch bamOrVcf {
	case "bams":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*bqsr.cram",
			dataDirAbs, strings.ToLower(species), sample, refVer)
	case "gvcfs":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/gvcfs/*%s.g.vcf.gz",
			dataDirAbs, strings.ToLower(species), sample, refVer, chrom)
	}

	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, fmt.Errorf("glob error: %w", err)
	}
	if len(matches) == 0 {
		return nil, fmt.Errorf("no %s files found in %s", bamOrVcf, sample)
	}
	if len(matches) > 1 {
		return matches, fmt.Errorf("multiple %s files found in %s", bamOrVcf, sample)
	}
	return matches, nil
}

func CreateGvcf(bam string, refFile string, chroms []SeqInfo, theGVCF string, gatkLogLevel string, caller string, dvVer string, modelType string, verbose bool) (string, error) {
	var regionArg string

	switch strings.ToLower(caller) {
	case "gatk":
		if len(chroms) == 1 {
			regionArg = chroms[0].ID
		} else {
			// Use a unique temp file per call to avoid races when multiple
			// goroutines run CreateGvcf concurrently.
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
		return theGVCF, runCmd(hapCmdStr, verbose)

	case "deepvariant":
		gvcfDir := filepath.Dir(theGVCF)
		gvcfName := filepath.Base(theGVCF)
		vcfName := strings.Replace(gvcfName, ".g.vcf.gz", ".vcf.gz", 1)

		var regions string
		if len(chroms) == 1 {
			regions = chroms[0].ID
		} else {
			// Write a BED file into gvcfDir (which is mounted as /output) so
			// the container can read it via its /output mount.
			f, err := os.CreateTemp(gvcfDir, "deepvariant_intervals_*.bed")
			if err != nil {
				return "", fmt.Errorf("creating DeepVariant interval list: %w", err)
			}
			defer os.Remove(f.Name())
			defer f.Close()
			for _, c := range chroms {
				fmt.Fprintf(f, "%s\t0\t%d\n", c.ID, c.Len)
			}
			regions = f.Name() // absolute host path inside gvcfDir → translatable to /output/...
		}

		bamDir := filepath.Dir(bam)
		bamName := filepath.Base(bam)
		refDir := filepath.Dir(refFile)
		refName := filepath.Base(refFile)

		safeRegion := strings.NewReplacer(string(os.PathSeparator), "_", ".", "_", "/", "_").Replace(regions)
		intermediateName := fmt.Sprintf("tmp_%s_%s", strings.TrimSuffix(bamName, filepath.Ext(bamName)), safeRegion)

		dvCmdStr := fmt.Sprintf(
			`docker run --rm -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
				`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
				`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
				`--intermediate_results_dir /output/%s`,
			bamDir, refDir, gvcfDir, dvVer,
			modelType, refName, bamName,
			regions, vcfName, gvcfName,
			intermediateName,
		)
		fmt.Printf("\n%s\n\n", dvCmdStr)
		return theGVCF, runCmd(dvCmdStr, verbose)

	default:
		return "", fmt.Errorf("unknown caller %q: must be \"gatk\" or \"deepvariant\"", caller)
	}
}

func VariantCallingDir(dataDir string, species string, refVer string, genomesDir string, refFasta string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool, gatkLogLevel string, quick bool, threadsPerJob int) error {

	// ── validate paths ────────────────────────────────────────────────────────
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		return fmt.Errorf("accessing data directory %s: %w", dataDir, err)
	}
	if !dInfo.IsDir() {
		return fmt.Errorf("data directory %s is not a directory", dataDir)
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		return fmt.Errorf("getting absolute path for data directory %s: %w", dataDir, err)
	}

	if species == "" {
		return fmt.Errorf("species name must not be empty")
	}
	if refVer == "" {
		return fmt.Errorf("reference version must not be empty")
	}

	resolvedFasta, resolveErr := utils.ResolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		return fmt.Errorf("could not resolve reference fasta: %w", resolveErr)
	}
	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil {
		return fmt.Errorf("accessing reference fasta file %s: %w", resolvedFasta, err)
	}
	if !fastaInfo.Mode().IsRegular() {
		return fmt.Errorf("reference fasta %s is not a regular file", resolvedFasta)
	}

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		return fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return fmt.Errorf("getting chromosomes and IDs: %w", err)
	}

	color.Cyan("All file paths valid\n…\n\n")

	// ── discover samples ──────────────────────────────────────────────────────
	color.Cyan("================================== Discovering samples =================================\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return fmt.Errorf("finding samples in %s: %w", pattern, err)
	}

	color.Green("SAMPLES FOUND:\n%s\n\n", strings.Repeat("-", 69))
	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, s)
			fmt.Println(s)
		}
	}
	color.Cyan("\nFound %d sample(s) for %s\n%s\n\n", len(samples), species, strings.Repeat("=", 34))

	// ── resolve valid samples ─────────────────────────────────────────────────
	fmt.Println("Looking for BAM files …")

	var (
		mu           sync.Mutex
		missingBams  []string
		multipleBams []string
		validSamples []SampleWork
		wg1          sync.WaitGroup
	)

	for _, sample := range samples {
		wg1.Add(1)
		go func(sample string) {
			defer wg1.Done()

			pat := fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*bqsr.cram",
				dataDirAbs, strings.ToLower(species), sample, refVer)
			cramFiles, cErr := filepath.Glob(pat)
			if cErr != nil {
				color.Red("Error finding BAM files: %v\n", cErr)
				mu.Lock()
				missingBams = append(missingBams, sample)
				mu.Unlock()
				return
			}

			switch len(cramFiles) {
			case 0:
				color.Red("[%s] bqsr.cram MISSING ❌\n", sample)
				mu.Lock()
				missingBams = append(missingBams, sample)
				mu.Unlock()
			case 1:
				color.Green("[%s] bqsr.cram FOUND ✅: %s\n", sample, cramFiles[0])
				color.Cyan("Validating BAM file …\n")
				if valErr := utils.ValidateBam(cramFiles[0], resolvedFasta, verbose, quick); valErr != nil {
					color.Red("[%s] bqsr.cram not valid: %v\n", sample, valErr)
					mu.Lock()
					missingBams = append(missingBams, sample)
					mu.Unlock()
				} else {
					color.Green("[%s] bqsr.cram valid ✅\n", sample)
					mu.Lock()
					validSamples = append(validSamples, SampleWork{
						Sample:  sample,
						Cram:    cramFiles[0],
						CramDir: filepath.Dir(filepath.Dir(filepath.Dir(cramFiles[0]))),
					})
					mu.Unlock()
				}
			default:
				color.Red("[%s] Multiple bqsr.cram files found — skipping ❌\n", sample)
				mu.Lock()
				multipleBams = append(multipleBams, sample)
				mu.Unlock()
			}
		}(sample)
	}
	wg1.Wait()

	color.Green("\nValid:    %d\n", len(validSamples))
	color.Red("Missing:  %d\n", len(missingBams))
	color.Yellow("Multiple: %d\n", len(multipleBams))
	fmt.Printf("\n%s\n\n", strings.Repeat("=", 89))

	if len(validSamples) == 0 {
		return fmt.Errorf("no valid samples found (missing: %d, multiple: %d)", len(missingBams), len(multipleBams))
	}

	color.Cyan("========= Creating gVCFs for %d valid samples =========\n\n", len(validSamples))

	// ── resource-aware concurrent gVCF creation ───────────────────────────────
	var (
		failedTasks []FailedTask
		failedMu    sync.Mutex
	)

	totalCores := runtime.NumCPU()
	if threadsPerJob <= 0 {
		threadsPerJob = 4
	}
	maxParallelJobs := totalCores / threadsPerJob
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
	}
	color.Cyan("Machine has %d cores. %d threads per job → max %d parallel jobs\n\n",
		totalCores, threadsPerJob, maxParallelJobs)

	type workItem struct {
		sw       SampleWork
		chroms   []SeqInfo
		label    string
		isContig bool
	}

	totalWork := len(validSamples) * (len(chroms) + 1)
	workCh := make(chan workItem, totalWork)
	var wg2 sync.WaitGroup

	for i := range maxParallelJobs {
		wg2.Add(1)
		go func(workerID int) {
			defer wg2.Done()
			for item := range workCh {
				sw := item.sw
				label := item.label

				cramName := filepath.Base(sw.Cram)
				var gvcfName string
				if item.isContig {
					gvcfName = strings.Replace(cramName, ".cram", ".contigs.g.vcf.gz", 1)
					gvcfName = strings.Replace(gvcfName, ".bam", ".contigs.g.vcf.gz", 1)
				} else {
					gvcfName = strings.Replace(cramName, ".cram", fmt.Sprintf(".%s.g.vcf.gz", label), 1)
					gvcfName = strings.Replace(gvcfName, ".bam", fmt.Sprintf(".%s.g.vcf.gz", label), 1)
				}
				gvcfPath := filepath.Join(sw.CramDir, refVer, "gvcfs", gvcfName)

				// Create the gvcfs directory if it doesn't exist
				if err := os.MkdirAll(filepath.Dir(gvcfPath), 0755); err != nil {
					color.Red("[Worker %d] [%s] Error creating directory for %s: %v\n\n", workerID, sw.Sample, label, err)
					failedMu.Lock()
					failedTasks = append(failedTasks, FailedTask{Sample: sw.Sample, Chrom: label, Reason: err})
					failedMu.Unlock()
					continue
				}

				var vcfFiles []string
				if item.isContig {
					if _, statErr := os.Stat(gvcfPath); statErr == nil {
						vcfFiles = []string{gvcfPath}
					}
				} else {
					vcfFiles, _ = FindBamOrVcfs(dataDirAbs, species, sw.Sample, refVer, "gvcfs", label)
				}

				switch len(vcfFiles) {
				case 0:
					color.Cyan("[Worker %d] [%s] gVCF not found for %s — creating …\n\n", workerID, sw.Sample, label)
					if _, err := CreateGvcf(sw.Cram, resolvedFasta, item.chroms, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
						color.Red("[Worker %d] [%s] Error creating gVCF for %s: %v\n\n", workerID, sw.Sample, label, err)
						failedMu.Lock()
						failedTasks = append(failedTasks, FailedTask{Sample: sw.Sample, Chrom: label, Reason: err})
						failedMu.Unlock()
					} else {
						color.Green("[Worker %d] [%s] gVCF for %s created successfully\n\n", workerID, sw.Sample, label)
					}

				case 1:
					vcf := vcfFiles[0]
					color.Green("[Worker %d] [%s] gVCF for %s exists: %s\n\n", workerID, sw.Sample, label, vcf)
					fmt.Printf("[Worker %d] [%s] checking integrity of %s …\n", workerID, sw.Sample, color.BlueString(vcf))
					if vErr := utils.ValidateGvcf(vcf, verbose, quick); vErr != nil {
						color.Red("[Worker %d] [%s] gVCF %s corrupted (%v) — re-creating\n", workerID, sw.Sample, color.BlueString(vcf), vErr)
						if _, err := CreateGvcf(sw.Cram, resolvedFasta, item.chroms, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
							color.Red("[Worker %d] [%s] Error re-creating gVCF %s: %v\n", workerID, sw.Sample, color.BlueString(vcf), err)
							failedMu.Lock()
							failedTasks = append(failedTasks, FailedTask{Sample: sw.Sample, Chrom: label, Reason: err})
							failedMu.Unlock()
						} else {
							color.Green("[Worker %d] [%s] gVCF %s re-created successfully\n", workerID, sw.Sample, color.BlueString(vcf))
						}
					} else {
						color.Green("[Worker %d] [%s] gVCF %s is valid ✅\n\n", workerID, sw.Sample, vcf)
					}

				default:
					color.Red("[Worker %d] [%s] Multiple gVCF files found for %s — please remove extras.\n\n", workerID, sw.Sample, label)
					failedMu.Lock()
					failedTasks = append(failedTasks, FailedTask{
						Sample: sw.Sample, Chrom: label,
						Reason: fmt.Errorf("multiple gVCF files found"),
					})
					failedMu.Unlock()
				}
			}
		}(i)
	}

	for _, sw := range validSamples {
		for _, chrom := range chroms {
			workCh <- workItem{sw: sw, chroms: []SeqInfo{chrom}, label: chrom.ID}
		}
	}
	if len(contigs) > 0 {
		for _, sw := range validSamples {
			workCh <- workItem{sw: sw, chroms: contigs, label: "contigs", isContig: true}
		}
	}
	close(workCh)
	wg2.Wait()

	// ── remove deepvariant intermediate directories ───────────────────────────
	if strings.ToLower(caller) == "deepvariant" {
		color.Cyan("Cleaning up DeepVariant intermediate results ...")
		for _, sw := range validSamples {
			// DeepVariant intermediate results are in sw.CramDir/refVer/gvcfs/tmp_*
			gvcfDir := filepath.Join(sw.CramDir, refVer, "gvcfs")
			pattern := filepath.Join(gvcfDir, "tmp_*")
			matches, _ := filepath.Glob(pattern)
			for _, m := range matches {
				os.RemoveAll(m)
			}
		}
	}

	// ── summary ───────────────────────────────────────────────────────────────
	fmt.Printf("\n%s FINAL SUMMARY %s\n", strings.Repeat("=", 29), strings.Repeat("=", 29))
	fmt.Printf("Machine cores:     %d\n", totalCores)
	fmt.Printf("Threads per job:   %d\n", threadsPerJob)
	fmt.Printf("Max parallel jobs: %d\n", maxParallelJobs)
	fmt.Printf("Samples processed: %d\n", len(validSamples))
	fmt.Printf("Missing BAMs:      %d %v\n", len(missingBams), missingBams)

	if len(failedTasks) > 0 {
		color.Red("Failed gVCF tasks: %d\n", len(failedTasks))
		for _, f := range failedTasks {
			fmt.Printf("  - Sample: %s, Chrom/Label: %s, Reason: %v\n", f.Sample, f.Chrom, f.Reason)
		}
	} else {
		color.Green("All gVCF tasks completed successfully!\n")
	}
	fmt.Printf("%s\n\n", strings.Repeat("=", 69))

	if !noMerging {
		if len(missingBams) == 0 && len(failedTasks) == 0 {
			color.Cyan("\nNo missing BAMs and no failed samples. Ready for merging.\n")

			var gvcfs []string
			for _, sw := range validSamples {
				gvcfDir := filepath.Join(sw.CramDir, refVer, "gvcfs")
				matches, _ := filepath.Glob(filepath.Join(gvcfDir, "*.g.vcf.gz"))
				gvcfs = append(gvcfs, matches...)
			}

			if len(gvcfs) > 0 {
				color.Cyan("Merging %d gVCFs ...\n", len(gvcfs))
				// We need a common output directory for merged results.
				// Using dataDir/species/refVer/VCFs
				outDir := filepath.Join(dataDirAbs, species, refVer, "VCFs")
				os.MkdirAll(outDir, 0755)

				// For now, we can call VariantCalling but it might be overkill as it re-runs gVCF creation.
				// Ideally we'd call the merge functions directly.
				// This is a placeholder for actual merging logic.
				color.Yellow("Direct merging of discovered gVCFs is not yet fully implemented in VariantCallingDir.")
			}
		} else {
			color.Yellow("Skipping merge due to missing BAMs (%d) or failed tasks (%d)\n", len(missingBams), len(failedTasks))
		}
	}

	if len(failedTasks) > 0 {
		return fmt.Errorf("%d gVCF task(s) failed", len(failedTasks))
	}
	return nil
}
