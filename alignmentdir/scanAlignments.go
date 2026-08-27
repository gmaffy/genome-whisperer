package alignmentdir

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// ============================= 1. Scan all dirs to see the stage of each sample ==================================== //

type FileInfo struct {
	Path         string
	Size         int64
	Present      bool
	Valid        bool
	ValidateErr  error
	IndexPresent bool
	IndexSize    int64
}

type SampleBamState struct {
	Sample    string
	SortedBam FileInfo
	RgmdBam   FileInfo
	RgmdCram  FileInfo
	BqsrBam   FileInfo
	BqsrCram  FileInfo

	OtherFiles []FileInfo
}

type ScanResult struct {
	Samples []SampleBamState
	Errors  []ScanError
}

type ScanError struct {
	Sample string
	Err    error
}

// --------------------------- Function to give details of files in bam directory ----------------------------------- //

func ScanAlignments(dataDir, species, refVer, genomesDir string, refFasta string, verbose, quick bool) ([]SampleBamState, error) {
	// Returns a list of sample structs. Each sample structs tell weather certain bams/crams are present or not
	// Imported alignment files are listed in SampleBamState
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

	dictFilePath := utils.DictPath(resolvedFasta)
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		color.Red("Reference dict file: %s does not exist\n", dictFilePath)
		log.Fatal(dicfErr)
	}

	// ---- 1. Discover sample bam-dirs ---- //
	pattern := filepath.Join(dataDirAbs, species, "*", "*", "reference_genomes", refVer, "bams")
	color.Yellow("Scanning with Glob Pattern: %s\n", pattern)

	matches, _ := filepath.Glob(pattern)
	color.Green("Found %d matching 'bams' directories.\n", len(matches))

	if len(matches) == 0 {
		color.Red("No matching 'bams' directories found.\n")
		return []SampleBamState{}, fmt.Errorf("no matching 'bams' directories found")
	}

	var scanResult []SampleBamState
	for _, bamDir := range matches {
		cleanBamDir := filepath.Clean(bamDir)
		// Extract sample name safely: <dataDirAbs>/<species>/<group>/<sample>/reference_genomes/<refVer>/bams
		sampleName := filepath.Base(filepath.Dir(filepath.Dir(filepath.Dir(cleanBamDir))))
		color.Green("Sample: %s: %s\n", sampleName, bamDir)

		sampleState := SampleBamState{Sample: sampleName}
		dirEntries, err := os.ReadDir(bamDir)

		if err != nil {
			color.Red("Error reading directory: %s\n", bamDir)
			return []SampleBamState{}, err
		}

		for _, de := range dirEntries {
			if de.IsDir() {
				color.Green("Found a dir: %s\n", de.Name())
			}
			name := de.Name()
			bamPath := filepath.Join(bamDir, name)

			info, infoErr := de.Info()
			bamSize := int64(0)
			if infoErr == nil {
				bamSize = info.Size()
			}
			if strings.HasSuffix(name, "sorted.bam") {
				sampleState.SortedBam = FileInfo{Path: bamPath, Size: bamSize, Present: true}

				valErr := utils.ValidateBam(bamPath, resolvedFasta, verbose, quick)
				if valErr == nil {
					sampleState.SortedBam.Valid = true
					//color.Green("[%s] has VALID sorted bam: %s\n", sampleName, bamPath)
				} else {
					sampleState.SortedBam.Valid = false
					//color.Yellow("[%s] has INVALID sorted bam: %s\n", sampleName, bamPath)
				}
				sampleState.SortedBam.ValidateErr = valErr
				setIndexInfo(&sampleState.SortedBam, verbose)

			} else if strings.HasSuffix(name, "rgmd.bam") || strings.HasSuffix(name, "RGMD.bam") {
				sampleState.RgmdBam = FileInfo{Path: bamPath, Size: bamSize, Present: true}

				valErr := utils.ValidateBam(bamPath, resolvedFasta, verbose, quick)
				if valErr == nil {
					sampleState.RgmdBam.Valid = true
					//color.Green("[%s] has VALID rgmd bam: %s\n", sampleName, bamPath)
				} else {
					sampleState.RgmdBam.Valid = false
					//color.Yellow("[%s] has INVALID rgmd bam: %s\n", sampleName, bamPath)
				}
				sampleState.RgmdBam.ValidateErr = valErr
				setIndexInfo(&sampleState.RgmdBam, verbose)

			} else if strings.HasSuffix(name, "rgmd.cram") || strings.HasSuffix(name, "RGMD.cram") {
				sampleState.RgmdCram = FileInfo{Path: bamPath, Size: bamSize, Present: true}

				valErr := utils.ValidateBam(bamPath, resolvedFasta, verbose, quick)
				if valErr == nil {
					sampleState.RgmdCram.Valid = true
					//color.Green("[%s] has VALID rgmd cram: %s\n", sampleName, bamPath)
				}
				sampleState.RgmdCram.ValidateErr = valErr
				setIndexInfo(&sampleState.RgmdCram, verbose)

			} else if strings.HasSuffix(name, "bqsr.bam") || strings.HasSuffix(name, "BQSR.bam") {
				sampleState.BqsrBam = FileInfo{Path: bamPath, Size: bamSize, Present: true}

				varErr := utils.ValidateBam(bamPath, resolvedFasta, verbose, quick)
				if varErr == nil {
					sampleState.BqsrBam.Valid = true
					color.Green("[%s] has VALID bqsr bam: %s\n", sampleName, bamPath)
				}

				setIndexInfo(&sampleState.BqsrBam, verbose)

			} else if strings.HasSuffix(name, "bqsr.cram") || strings.HasSuffix(name, "BQSR.cram") {
				sampleState.BqsrCram = FileInfo{Path: bamPath, Size: bamSize, Present: true}

				varErr := utils.ValidateBam(bamPath, resolvedFasta, verbose, quick)
				if varErr == nil {
					sampleState.BqsrCram.Valid = true
					color.Green("[%s] has VALID bqsr cram: %s\n", sampleName, bamPath)
				}
				setIndexInfo(&sampleState.BqsrCram, verbose)
			} else {
				if !strings.HasSuffix(name, ".bai") && !strings.HasSuffix(name, ".csi") && !strings.HasSuffix(name, ".crai") && !strings.HasSuffix(name, ".pdf") {
					color.Blue("[%s] Random artifact found: %s\n", sampleName, bamPath)
					sampleState.OtherFiles = append(sampleState.OtherFiles, FileInfo{Path: bamPath, Size: bamSize, Present: true})

				}

			}
		}
		scanResult = append(scanResult, sampleState)

	}

	return scanResult, nil
}

// ------------------------------- For each sample list steps needed to get perfect dir ----------------------------- //

type Step struct {
	Program string   // e.g. "samtools", "gatk", "bwa-mem"
	Inputs  []string // arguments, input files, flags, etc.
}

type Action struct {
	Sample  string
	Steps   []Step
	Perfect bool
}

func isUsable(f FileInfo) bool {
	return f.Present && f.Valid
}

func hasIndex(f FileInfo) bool {
	return f.IndexPresent && f.IndexSize > 0
}

func GetSampleActions(s SampleBamState) Action {
	var steps []Step

	// --- Check terminal goal ---
	bqsrCramOK := isUsable(s.BqsrCram) && hasIndex(s.BqsrCram)
	rgmdCramOK := isUsable(s.RgmdCram) && hasIndex(s.RgmdCram)

	// --- If already perfect ---
	if bqsrCramOK && rgmdCramOK &&
		!s.SortedBam.Present &&
		!s.RgmdBam.Present &&
		!s.BqsrBam.Present &&
		len(s.OtherFiles) == 0 {

		return Action{
			Sample: s.Sample,
			Steps: []Step{
				{
					Program: "Skip",
					Inputs:  []string{"✅PASS - All required cram files present"},
				},
			},
			Perfect: true,
		}

	}

	// --- Ensure rgmd.cram ---
	if !rgmdCramOK {
		if isUsable(s.RgmdBam) {
			steps = append(steps, Step{
				Program: "bam_to_cram",
				Inputs:  []string{s.RgmdBam.Path},
			})
		} else if isUsable(s.SortedBam) {
			//steps = append(steps, "Run markdup: sorted.bam → rgmd.bam")

			steps = append(steps, Step{
				Program: "markdup",
				Inputs:  []string{s.SortedBam.Path},
			})

			//steps = append(steps, "Convert rgmd.bam → rgmd.cram")
			steps = append(steps, Step{
				Program: "bam_to_cram",
				Inputs:  []string{s.RgmdBam.Path},
			})
		} else {
			//steps = append(steps, "Missing sorted.bam → must rerun alignment (bwa-mem)")
			steps = append(steps, Step{
				Program: "scratch",
				Inputs:  []string{s.Sample},
			})
		}
	}

	// --- Ensure rgmd.cram index ---
	if isUsable(s.RgmdCram) && !hasIndex(s.RgmdCram) {
		//steps = append(steps, "Index rgmd.cram")
		steps = append(steps, Step{
			Program: "indexBam",
			Inputs:  []string{s.RgmdBam.Path},
		})
	}

	// --- Ensure bqsr.cram ---
	if !bqsrCramOK {
		if isUsable(s.BqsrBam) {
			steps = append(steps, Step{
				Program: "bam_to_cram",
				Inputs:  []string{s.BqsrBam.Path},
			})

		} else if isUsable(s.RgmdBam) {
			//steps = append(steps, "Run BQSR → bqsr.bam")
			//steps = append(steps, "Convert bqsr.bam → bqsr.cram")
			if !isUsable(s.RgmdCram) {
				steps = append(steps, Step{
					Program: "bam_to_cram",
					Inputs:  []string{s.RgmdBam.Path},
				})
			}
			steps = append(steps, Step{
				Program: "bqsr",
				Inputs:  []string{s.RgmdBam.Path},
			})
		} else if isUsable(s.RgmdCram) {
			steps = append(steps, Step{
				Program: "bqsr",
				Inputs:  []string{s.RgmdBam.Path},
			})
			//steps = append(steps, "Cannot run BQSR: missing rgmd input")
		} else {
			steps = append(steps, Step{
				Program: "scratch",
				Inputs:  []string{s.Sample},
			})
		}
	}

	// --- Ensure bqsr.cram index ---
	if isUsable(s.BqsrCram) && !hasIndex(s.BqsrCram) {
		//steps = append(steps, "Index bqsr.cram")
		steps = append(steps, Step{
			Program: "indexBam",
			Inputs:  []string{s.BqsrBam.Path},
		})
	}

	// --- Cleanup phase ---
	if bqsrCramOK && rgmdCramOK {
		if s.SortedBam.Present {
			//steps = append(steps, "Delete sorted.bam")
			steps = append(steps, Step{
				Program: "rm",
				Inputs:  []string{s.SortedBam.Path},
			})
		}
		if s.RgmdBam.Present {
			//steps = append(steps, "Delete rgmd.bam")
			steps = append(steps, Step{
				Program: "rm",
				Inputs:  []string{s.RgmdBam.Path},
			})
		}
		if s.BqsrBam.Present {
			//steps = append(steps, "Delete bqsr.bam")
			steps = append(steps, Step{
				Program: "rm",
				Inputs:  []string{s.BqsrBam.Path},
			})
		}
		if len(s.OtherFiles) > 0 {
			//steps = append(steps, "Remove unexpected extra files")
			for _, file := range s.OtherFiles {
				steps = append(steps, Step{
					Program: "rm",
					Inputs:  []string{file.Path},
				})
			}
		}
	}

	return Action{
		Sample:  s.Sample,
		Steps:   steps,
		Perfect: false,
	}
}

//func EvaluateAll(samples []SampleBamState) []Action {
//	results := make([]Action, 0, len(samples))
//	for _, s := range samples {
//		results = append(results, evaluateSample(s))
//	}
//	return results
//}
