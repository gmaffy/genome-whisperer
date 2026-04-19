package alignment

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/variants"
)

// -----------------------------------------------------------------------------
// Types
// -----------------------------------------------------------------------------

type FileInfo struct {
	Path        string
	Size        int64
	Present     bool
	Valid       bool
	ValidateErr error
}

type SampleBamStatus struct {
	Sample string

	SortedBam   FileInfo
	SortedIndex FileInfo
	RgmdBam     FileInfo
	RgmdIndex   FileInfo
	BqsrBam     FileInfo
	BqsrIndex   FileInfo

	OtherFiles []FileInfo
}

type ScanResult struct {
	Samples []SampleBamStatus
	Errors  []ScanError
}

type ScanError struct {
	Sample string
	Err    error
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------

func ScanAlignmentStages(dataDir, species, refVer, genomesDir string, refFasta string, verbose, quick bool) ScanResult {

	color.Cyan("\n=================== STARTING ALIGNMENT SCAN ===================\n")
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		color.Red("Error accessing data directory: %s\n", dataDir)
		log.Fatal(err)
	}
	if !dInfo.IsDir() {
		color.Red("Data directory %s is not a directory\n", dataDir)
		log.Fatal(err)
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		color.Red("Error getting absolute path for data directory: %s\n", dataDir)
		log.Fatal(err)
	}

	if species == "" || refVer == "" {
		color.Red("Please provide species name and reference version name")
		log.Fatal("Missing species or refVer")
	}

	resolvedFasta, resolveErr := resolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		color.Red("Could not resolve reference fasta: %v\n", resolveErr)
		log.Fatal(resolveErr)
	}

	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil || !fastaInfo.Mode().IsRegular() {
		color.Red("Reference fasta file: %s is not a valid regular file\n", resolvedFasta)
		log.Fatal(err)
	}

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		color.Red("Reference dict file: %s does not exist\n", dictFilePath)
		log.Fatal(dicfErr)
	}

	// ---- 1. Discover sample bam-dirs ---- //
	pattern := filepath.Join(dataDirAbs, species, "*", "*", "reference_genomes", refVer, "bams")
	color.Yellow("Scanning with Glob Pattern: %s\n", pattern)

	matches, _ := filepath.Glob(pattern)
	color.Green("Found %d matching 'bams' directories.\n", len(matches))

	type entry struct {
		sample string
		bamDir string
	}
	entries := make([]entry, 0, len(matches))
	seen := make(map[string]struct{}, len(matches))

	for _, bamDir := range matches {
		cleanBamDir := filepath.Clean(bamDir)
		// Extract sample name safely: <dataDirAbs>/<species>/<group>/<sample>/reference_genomes/<refVer>/bams
		sample := filepath.Base(filepath.Dir(filepath.Dir(filepath.Dir(cleanBamDir))))

		color.HiBlack("  -> Extracted Sample Name: [%s] from path: %s\n", sample, cleanBamDir)

		if _, dup := seen[sample]; dup {
			continue
		}
		seen[sample] = struct{}{}
		entries = append(entries, entry{sample: sample, bamDir: cleanBamDir})
	}

	// ---- 2. Scan + validate all samples in parallel ---- //
	results := make([]SampleBamStatus, len(entries))
	var (
		errMu         sync.Mutex
		collectedErrs []ScanError
		wg            sync.WaitGroup
	)

	for i, e := range entries {
		wg.Add(1)
		go func(idx int, sample, bamDir string) {
			defer wg.Done()
			status, err := scanOneSample(sample, bamDir, resolvedFasta, verbose, quick)

			// Summary output for the sample
			bqsrStatus := "MISSING"
			if status.BqsrBam.Present {
				if status.BqsrBam.Valid {
					bqsrStatus = color.GreenString("VALID")
				} else {
					bqsrStatus = color.RedString("INVALID (Err: %v)", status.BqsrBam.ValidateErr)
				}
			}

			// Using safe formatting based on your previous print structure
			color.Cyan("[%s] Scan Complete | BqsrBam Status: %s\n", sample, bqsrStatus)

			results[idx] = status
			if err != nil {
				errMu.Lock()
				collectedErrs = append(collectedErrs, ScanError{Sample: sample, Err: err})
				errMu.Unlock()
			}
		}(i, e.sample, e.bamDir)
	}
	wg.Wait()

	color.Cyan("=================== ALIGNMENT SCAN FINISHED ===================\n\n")

	return ScanResult{
		Samples: results,
		Errors:  collectedErrs,
	}
}

// -----------------------------------------------------------------------------
// Per-sample scanner
// -----------------------------------------------------------------------------

var knownSuffixes = []string{
	".sorted.bam", ".sorted.bai", ".rgmd.bam", ".rgmd.bai", ".rgmd.cram",
	".rgmd.cram.crai", ".RGMD.bam", ".RGMD.bai", ".RGMD.cram", ".RGMD.cram.crai",
	".bqsr.bam", ".bqsr.bai", ".bqsr.cram", ".bqsr.cram.crai",
}

func scanOneSample(sample, bamDir, refFasta string, verbose, quick bool) (SampleBamStatus, error) {
	status := SampleBamStatus{Sample: sample}
	printPrefix := color.YellowString("[%s SCAN]", sample)

	fmt.Printf("%s Opening directory: %s\n", printPrefix, bamDir)
	dirEntries, err := os.ReadDir(bamDir)
	if err != nil {
		color.Red("%s Could not read directory: %v\n", printPrefix, err)
		return status, err
	}

	type bamToValidate struct {
		fi   *FileInfo
		path string
		name string
	}
	var toValidate []bamToValidate

	// Helper to check multiple suffixes safely
	hasSuffixAny := func(s string, suffixes ...string) bool {
		for _, suf := range suffixes {
			if strings.HasSuffix(s, suf) {
				return true
			}
		}
		return false
	}

	for _, de := range dirEntries {
		if de.IsDir() {
			continue
		}
		name := de.Name()
		fullPath := filepath.Join(bamDir, name)

		info, infoErr := de.Info()
		size := int64(0)
		if infoErr == nil {
			size = info.Size()
		}

		lowerName := strings.ToLower(name)
		matched := false

		// 1. Sorted Files
		if hasSuffixAny(lowerName, "sorted.bam") {
			status.SortedBam = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil}
			if status.SortedBam.Present {
				status.SortedBam.Valid = false
				toValidate = append(toValidate, bamToValidate{fi: &status.SortedBam, path: fullPath, name: name})
			}
			fmt.Printf("%s \t-> MATCHED: %s (Sorted BAM)\n", printPrefix, name)
			matched = true
		} else if hasSuffixAny(lowerName, "sorted.bai") {
			status.SortedIndex = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil, Valid: infoErr == nil}
			fmt.Printf("%s \t-> MATCHED: %s (Sorted Index)\n", printPrefix, name)
			matched = true
		} else

		// 2. RGMD Files
		if hasSuffixAny(lowerName, "rgmd.bam", "rgmd.cram") {
			status.RgmdBam = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil}
			if status.RgmdBam.Present {
				status.RgmdBam.Valid = false
				toValidate = append(toValidate, bamToValidate{fi: &status.RgmdBam, path: fullPath, name: name})
			}
			fmt.Printf("%s \t-> MATCHED: %s (RGMD BAM/CRAM)\n", printPrefix, name)
			matched = true
		} else if hasSuffixAny(lowerName, "rgmd.bai", "rgmd.crai", "rgmd.cram.crai") {
			status.RgmdIndex = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil, Valid: infoErr == nil}
			fmt.Printf("%s \t-> MATCHED: %s (RGMD Index)\n", printPrefix, name)
			matched = true
		} else

		// 3. BQSR Files
		if hasSuffixAny(lowerName, "bqsr.bam", "bqsr.cram") {
			status.BqsrBam = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil}
			if status.BqsrBam.Present {
				status.BqsrBam.Valid = false
				toValidate = append(toValidate, bamToValidate{fi: &status.BqsrBam, path: fullPath, name: name})
			}
			fmt.Printf("%s \t-> MATCHED: %s (BQSR BAM/CRAM)\n", printPrefix, name)
			matched = true
		} else if hasSuffixAny(lowerName, "bqsr.bai", "bqsr.crai", "bqsr.cram.crai") {
			status.BqsrIndex = FileInfo{Path: fullPath, Size: size, Present: infoErr == nil, Valid: infoErr == nil}
			fmt.Printf("%s \t-> MATCHED: %s (BQSR Index)\n", printPrefix, name)
			matched = true
		}

		// 4. Other Unrecognized Files
		if !matched {
			color.HiBlack("%s \t-> IGNORED: %s (Unrecognized file format)\n", printPrefix, name)
			status.OtherFiles = append(status.OtherFiles, FileInfo{
				Path:    fullPath,
				Size:    size,
				Present: infoErr == nil,
				Valid:   infoErr == nil,
			})
		}
	}

	// ---- Validate BAM/CRAM files concurrently ---- //
	if len(toValidate) > 0 {
		var valWg sync.WaitGroup
		for _, v := range toValidate {
			valWg.Add(1)
			go func(fi *FileInfo, path, name string) {
				defer valWg.Done()

				fmt.Printf("%s Validating %s with Samtools...\n", printPrefix, name)

				if valErr := variants.ValidateBam(path, refFasta, verbose, quick); valErr != nil {
					color.Red("%s ❌ VALIDATION FAILED for %s: %v\n", printPrefix, name, valErr)
					fi.Valid = false
					fi.ValidateErr = valErr
				} else {
					color.Green("%s ✅ VALIDATION PASSED for %s\n", printPrefix, name)
					fi.Valid = true
				}
			}(v.fi, v.path, v.name)
		}
		valWg.Wait()
	}

	return status, nil
}

func isKnownSuffix(name string) bool {
	for _, s := range knownSuffixes {
		if strings.HasSuffix(name, s) {
			return true
		}
	}
	return false
}
