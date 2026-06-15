package variants

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/gmaffy/genome-whisperer/utils"
)

func VariantCallingConfig(configFile string, species string, refVer string, gatkLogLevel string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool, noHardFilter bool, hfcfg utils.HardFilterConfig, threads int) (string, error) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return "", fmt.Errorf("error reading config: %v", err)
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Bams:", cfg.Bams)
	fmt.Println("Output:", cfg.OutputDir)

	refFile := cfg.Reference
	if _, rErr := os.Stat(refFile); rErr != nil {
		fmt.Printf("Reference file %s does not exist\n", refFile)
		return "", fmt.Errorf("reference file %s does not exist", refFile)
	}

	bams := cfg.Bams
	if len(bams) == 0 {
		fmt.Println("You must provide at least one BAM file")
		return "", fmt.Errorf("you must provide at least one BAM file")
	}
	for _, b := range bams {
		if _, err := os.Stat(b); err != nil {
			fmt.Printf("BAM file %s is not a valid file path: %v\n", b, err)
			return "", fmt.Errorf("BAM file %s is not a valid file path: %v", b, err)
		}
	}

	outputDir := cfg.OutputDir
	outInfo, outErr := os.Stat(outputDir)
	if outErr != nil {
		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory %s does not exist. Attempting to create it.\n", outputDir)
			if createErr := os.MkdirAll(outputDir, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", outputDir, createErr)
				return "", fmt.Errorf("failed to create output directory %s: %v", outputDir, createErr)
			}
			fmt.Printf("Output directory %s created successfully.\n", outputDir)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", outputDir, outErr)
			return "", fmt.Errorf("error accessing output directory %s: %v", outputDir, outErr)
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output path %s is not a directory\n", outputDir)
		return "", fmt.Errorf("output path %s is not a directory", outputDir)
	}

	logFilePath := filepath.Join(outputDir, "variant_calling.log")
	finalVcf, fErr := VariantCalling(bams, refFile, species, refVer, outputDir, caller, merger, noMerging, noHardFilter, hfcfg, verbose, gatkLogLevel, threads, logFilePath, dvVer, modelType)
	fmt.Printf("Error calling variants: %v\nNo multisample VCF: %s\n", fErr, finalVcf)

	return finalVcf, fErr

}
