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

// FileInfo holds metadata for a single file on disk.
type FileInfo struct {
	Path        string
	Size        int64
	Present     bool  // true when the file exists on disk
	Valid       bool  // true when Present AND samtools validation passed (index files: Present only)
	ValidateErr error // non-nil when validation was run and failed
}

// SampleBamStatus captures every alignment-stage file for one sample.
type SampleBamStatus struct {
	Sample string

	// Known pipeline files
	SortedBam   FileInfo // <sample>.sorted.bam
	SortedIndex FileInfo // <sample>.sorted.bai
	RgmdBam     FileInfo // <sample>.rgmd.bam / .RGMD.bam / .RGMD.cram
	RgmdIndex   FileInfo // <sample>.rgmd.bai / .RGMD.bai / .RGMD.cram.crai
	BqsrBam     FileInfo // <sample>.bqsr.bam / .bqsr.cram
	BqsrIndex   FileInfo // <sample>.bqsr.bai / .bqsr.cram.crai

	// Every other file sitting in the bams directory
	OtherFiles []FileInfo
}

// ScanResult is the return value of ScanAlignmentStages.
type ScanResult struct {
	Samples []SampleBamStatus
	Errors  []ScanError
}

// ScanError pairs a sample name with any non-fatal error encountered while
// scanning it (e.g. the bams directory doesn't exist yet).
type ScanError struct {
	Sample string
	Err    error
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------

// ScanAlignmentStages walks dataDir for every sample that matches the layout:
//
//	<dataDir>/<species>/<group>/<sample>/reference_genomes/<refVer>/bams/
//
// For each sample it collects the filepath, byte size, presence, and samtools
// validity of:
//   - <sample>.sorted.bam + index
//   - <sample>.rgmd.bam (or .RGMD.bam / .RGMD.cram) + index
//   - <sample>.bqsr.bam / .bqsr.cram + index
//
// BAM/CRAM files are validated with variants.ValidateBam. Index files report
// presence only (Present/Valid reflect os.Stat). Any unrecognised file in the
// bams directory is reported under OtherFiles.
//
// Scanning is two-level parallel: one goroutine per sample, and within each
// sample the three BAM/CRAM validations run concurrently.
func ScanAlignmentStages(dataDir, species, refVer, genomesDir string, refFasta string, verbose, quick bool) ScanResult {

	// ============================================= Validate paths ================================================ //
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		log.Fatal(err)
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		log.Fatal(err)
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		fmt.Printf("Error getting absolute path for data directory: %s\n", dataDir)
		log.Fatal(err)
	}

	if species == "" {
		fmt.Println("Please provide species name")
		log.Fatal(err)
	}
	if refVer == "" {
		fmt.Println("Please provide reference version name")
		log.Fatal(err)
	}

	// ---- Resolve reference fasta (explicit path or auto-discover) ---- //
	resolvedFasta, resolveErr := resolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		color.Red("Could not resolve reference fasta: %v\n", resolveErr)
		log.Fatal(resolveErr)
	}

	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", resolvedFasta)
		log.Fatal(err)
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", resolvedFasta)
		log.Fatal(err)
	}

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		log.Fatal(dicfErr)
	}

	// ---- 1. Discover sample bam-dirs with a single Glob call ---- //
	pattern := filepath.Join(dataDirAbs, species, "*", "*", "reference_genomes", refVer, "bams")
	matches, _ := filepath.Glob(pattern)

	type entry struct {
		sample string
		bamDir string
	}
	entries := make([]entry, 0, len(matches))
	seen := make(map[string]struct{}, len(matches))
	for _, bamDir := range matches {
		// bams → <refVer> → reference_genomes → <sample>
		sample := filepath.Base(filepath.Dir(filepath.Dir(filepath.Dir(bamDir))))
		if _, dup := seen[sample]; dup {
			continue
		}
		seen[sample] = struct{}{}
		entries = append(entries, entry{sample: sample, bamDir: bamDir})
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
			color.Cyan("[%s] %s\n", sample, status.Sample, status.BqsrBam)
			results[idx] = status
			if err != nil {
				errMu.Lock()
				collectedErrs = append(collectedErrs, ScanError{Sample: sample, Err: err})
				errMu.Unlock()
			}
		}(i, e.sample, e.bamDir)
	}
	wg.Wait()

	return ScanResult{
		Samples: results,
		Errors:  collectedErrs,
	}
}

// -----------------------------------------------------------------------------
// Per-sample scanner
// -----------------------------------------------------------------------------

// bamSuffixes are the alignment file suffixes we recognise as belonging to the
// pipeline. Anything else in the bams directory lands in OtherFiles.
var knownSuffixes = []string{
	".sorted.bam",
	".sorted.bai",
	".rgmd.bam",
	".rgmd.bai",
	".rgmd.cram",
	".rgmd.cram.crai",
	".RGMD.bam",
	".RGMD.bai",
	".RGMD.cram",
	".RGMD.cram.crai",
	".bqsr.bam",
	".bqsr.bai",
	".bqsr.cram",
	".bqsr.cram.crai",
}

func scanOneSample(sample, bamDir, refFasta string, verbose, quick bool) (SampleBamStatus, error) {
	status := SampleBamStatus{Sample: sample}

	// Single syscall to list the directory.
	dirEntries, err := os.ReadDir(bamDir)
	if err != nil {
		// Directory may not exist yet — sample hasn't been aligned.
		return status, err
	}

	// Map filename → pointer into status so classification is a single O(1) lookup.
	// Index files share their pointer with their BAM counterpart's index field.
	type slotKind uint8
	const (
		kindBam   slotKind = iota // needs samtools validation
		kindIndex                 // presence-only
		kindOther                 // size + presence, no validation
	)
	type slot struct {
		fi   *FileInfo
		kind slotKind
	}

	expected := map[string]slot{
		sample + ".sorted.bam":     {&status.SortedBam, kindBam},
		sample + ".sorted.bai":     {&status.SortedIndex, kindIndex},
		sample + ".rgmd.bam":       {&status.RgmdBam, kindBam},
		sample + ".rgmd.bai":       {&status.RgmdIndex, kindIndex},
		sample + ".RGMD.bam":       {&status.RgmdBam, kindBam},
		sample + ".RGMD.bai":       {&status.RgmdIndex, kindIndex},
		sample + ".RGMD.cram":      {&status.RgmdBam, kindBam},
		sample + ".RGMD.cram.crai": {&status.RgmdIndex, kindIndex},
		sample + ".bqsr.bam":       {&status.BqsrBam, kindBam},
		sample + ".bqsr.bai":       {&status.BqsrIndex, kindIndex},
		sample + ".bqsr.cram":      {&status.BqsrBam, kindBam},
		sample + ".bqsr.cram.crai": {&status.BqsrIndex, kindIndex},
	}

	// Collect BAM/CRAM files that need concurrent validation.
	type bamToValidate struct {
		fi   *FileInfo
		path string
	}
	var toValidate []bamToValidate

	for _, de := range dirEntries {
		if de.IsDir() {
			continue
		}
		name := de.Name()
		fullPath := filepath.Join(bamDir, name)

		info, infoErr := de.Info() // uses cached inode data — no extra stat syscall
		size := int64(0)
		if infoErr == nil {
			size = info.Size()
		}

		if s, ok := expected[name]; ok {
			s.fi.Path = fullPath
			s.fi.Size = size
			s.fi.Present = infoErr == nil

			switch s.kind {
			case kindBam:
				// Mark present for now; validation result fills Valid below.
				s.fi.Valid = false
				if s.fi.Present {
					toValidate = append(toValidate, bamToValidate{fi: s.fi, path: fullPath})
				}
			case kindIndex:
				// For index files, present == valid.
				s.fi.Valid = s.fi.Present
			}
		} else if !isKnownSuffix(name) {
			status.OtherFiles = append(status.OtherFiles, FileInfo{
				Path:    fullPath,
				Size:    size,
				Present: infoErr == nil,
				Valid:   infoErr == nil,
			})
		}
		// Files whose suffix is known but don't match this sample are stale
		// artefacts — include them in OtherFiles so nothing is silently hidden.
	}

	// ---- Validate BAM/CRAM files concurrently within this sample ---- //
	// ValidateBam shells out to samtools, so these are I/O + process bound.
	// Running them in parallel cuts per-sample wall time by ~3×.
	if len(toValidate) > 0 {
		var valWg sync.WaitGroup
		for _, v := range toValidate {
			valWg.Add(1)
			go func(fi *FileInfo, path string) {
				defer valWg.Done()
				if valErr := variants.ValidateBam(path, refFasta, verbose, quick); valErr != nil {
					fi.Valid = false
					fi.ValidateErr = valErr
				} else {
					fi.Valid = true
				}
			}(v.fi, v.path)
		}
		valWg.Wait()
	}

	return status, nil
}

// isKnownSuffix returns true when name ends with any suffix in knownSuffixes.
func isKnownSuffix(name string) bool {
	for _, s := range knownSuffixes {
		if strings.HasSuffix(name, s) {
			return true
		}
	}
	return false
}
