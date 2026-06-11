package variants

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

type sampleRefDir struct {
	Sample string
	RefDir string
}

func listSampleRefDirs(dataDirAbs string, speciesLower string) ([]sampleRefDir, error) {
	pattern := filepath.Join(dataDirAbs, speciesLower, "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, err
	}

	samples := make([]sampleRefDir, 0, len(matches))
	for _, match := range matches {
		sample := filepath.Base(filepath.Dir(match))
		samples = append(samples, sampleRefDir{
			Sample: sample,
			RefDir: match,
		})
	}
	return samples, nil
}

func listCramFiles(bamDir string) ([]string, error) {
	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return nil, err
	}
	crams := make([]string, 0)
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()
		if strings.HasSuffix(name, "bqsr.cram") {
			crams = append(crams, filepath.Join(bamDir, name))
		}
	}
	return crams, nil
}

func buildGvcfSet(gvcfDir string) (map[string]struct{}, error) {
	entries, err := os.ReadDir(gvcfDir)
	if err != nil {
		return nil, err
	}
	gvcfs := make(map[string]struct{}, len(entries))
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		gvcfs[entry.Name()] = struct{}{}
	}
	return gvcfs, nil
}

func validateGvcf(vcf string, verbose bool) error {
	valStr := fmt.Sprintf(`gunzip -t %s && bcftools view -h %s > /dev/null && bcftools view %s > /dev/null`, vcf, vcf, vcf)
	if verbose {
		return utils.RunBashCmdVerbose(valStr)
	}
	return utils.RunBashCmd(valStr)
}

// VariantCallingDirOptimised reduces repeated filesystem scans and skips unnecessary gVCF creation.
func VariantCallingDirOptimised(dataDir string, species string, refVer string, refFasta string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool, gatkLogLevel string) {
	// ============================================= Check paths ================================================ //
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
	if refFasta == "" {
		fmt.Println("Please provide reference name")
		return
	}

	fastaInfo, err := os.Stat(refFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", refFasta)
		return
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", refFasta)
		return
	}

	dictFilePath := refFasta[:len(refFasta)-len(filepath.Ext(refFasta))] + ".dict"
	_, dicfErr := os.Stat(dictFilePath)
	if dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		fmt.Printf("Error getting chromosomes and IDs: %v\n", err)
		return
	}

	color.Green("All file paths valid\n....................................................\n\n")

	// ================================== Checking Samples in dir Structure ===================================== //
	color.Green("Checking Samples in dir structure ...\n\n")
	sampleRefs, err := listSampleRefDirs(dataDirAbs, strings.ToLower(species))
	if err != nil {
		fmt.Printf("Error scanning sample directories: %v\n", err)
		return
	}

	fmt.Printf("there are %v Samples in the data directory for %s\n==================================\n\n", len(sampleRefs), species)

	// -------------------------------------- Check for bams ------------------------------------------------------- //
	fmt.Println("Looking for bams ....")

	missingBams := []string{}
	multipleBams := []string{}

	for _, sampleRef := range sampleRefs {
		bamDir := filepath.Join(sampleRef.RefDir, refVer, "bams")
		cramFiles, cramErr := listCramFiles(bamDir)
		if cramErr != nil {
			missingBams = append(missingBams, sampleRef.Sample)
			continue
		}

		if len(cramFiles) == 0 {
			missingBams = append(missingBams, sampleRef.Sample)
			continue
		}
		if len(cramFiles) > 1 {
			multipleBams = append(multipleBams, sampleRef.Sample)
			continue
		}

		cram := cramFiles[0]
		cramName := filepath.Base(cram)

		gvcfDir := filepath.Join(sampleRef.RefDir, refVer, "gvcfs")
		if err := os.MkdirAll(gvcfDir, 0755); err != nil {
			fmt.Printf("Error creating gvcf directory: %s: %v\n", gvcfDir, err)
			continue
		}

		gvcfSet, gErr := buildGvcfSet(gvcfDir)
		if gErr != nil {
			fmt.Printf("Error reading gvcf directory: %s: %v\n", gvcfDir, gErr)
			continue
		}

		fmt.Printf("Cram file found in %s: %v\n", sampleRef.Sample, cram)
		for _, chrom := range chroms {
			gvcfName := strings.Replace(cramName, ".cram", fmt.Sprintf("%s.g.vcf.gz", chrom.ID), 1)
			gvcfPath := filepath.Join(gvcfDir, gvcfName)

			if _, ok := gvcfSet[gvcfName]; ok {
				color.Green("CHECKING VALIDITY OF THE gVCF file %s", color.BlueString(gvcfPath))
				if vErr := validateGvcf(gvcfPath, verbose); vErr != nil {
					fmt.Printf("Error validating gvcf file: %v\n", vErr)
					fmt.Println("Re-creating gvcf file")
					if _, err4 := CreateGvcf(cram, refFasta, []SeqInfo{chrom}, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err4 != nil {
						fmt.Printf("Error creating gvcf file - %s: %v\n", gvcfPath, err4)
					}
				} else {
					fmt.Println("gVCF is valid")
				}
				continue
			}

			fmt.Printf("VCF file not found in %s: %s\n", sampleRef.Sample, chrom.ID)
			if _, err4 := CreateGvcf(cram, refFasta, []SeqInfo{chrom}, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err4 != nil {
				fmt.Printf("Error creating gvcf file - %s: %v\n", gvcfPath, err4)
			} else {
				gvcfSet[gvcfName] = struct{}{}
			}
		}

		contigGvcfName := strings.Replace(cramName, ".cram", ".contigs.g.vcf.gz", 1)
		contigGvcfPath := filepath.Join(gvcfDir, contigGvcfName)
		if _, ok := gvcfSet[contigGvcfName]; ok {
			color.Green("CHECKING VALIDITY OF THE gVCF file %s", color.BlueString(contigGvcfPath))
			if vErr := validateGvcf(contigGvcfPath, verbose); vErr != nil {
				fmt.Printf("Error validating gvcf file: %v\n", vErr)
				fmt.Println("Re-creating contig gvcf file")
				if _, err4 := CreateGvcf(cram, refFasta, contigs, contigGvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err4 != nil {
					fmt.Printf("Error creating gvcf file - %s: %v\n", contigGvcfPath, err4)
				}
			} else {
				fmt.Println("gVCF is valid")
			}
		} else {
			if _, err4 := CreateGvcf(cram, refFasta, contigs, contigGvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err4 != nil {
				fmt.Printf("Error creating gvcf file - %s: %v\n", contigGvcfPath, err4)
			}
		}
	}

	if len(multipleBams) > 0 {
		fmt.Printf("There are %v samples with multiple bams for %s\n==================================\n", len(multipleBams), species)
	}
	fmt.Printf("There are %v missing bams for %s\n==================================\n", len(missingBams), species)
}
