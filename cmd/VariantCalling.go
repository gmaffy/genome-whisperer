/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"
)

var VariantCallingCmd = &cobra.Command{
	Use:   "VariantCalling",
	Short: "Creates a multi-sample VCF file from bam files using GATK best practices OR Deepvariant and GLNEXUS",
	Long: `VariantCalling
        - Calls and hard filters SNPs and Indels using GATK best practices from bams generated from short reads
        - Calls and hard filters SNPs and Indels using Deepvariant from bams generated from long reads
        - Can use glenexus or GATK to merge gvcfs`,
	Run: func(cmd *cobra.Command, args []string) {

		fmt.Println("VariantCalling called")
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		threads, jErr := cmd.Flags().GetInt("threads")
		if jErr != nil {
			log.Fatalf("Error getting threads flag: %v", jErr)
		}

		refFile, refErr := cmd.Flags().GetString("reference")
		if refErr != nil {
			log.Fatalf("Error getting reference flag: %v", refErr)
		}

		speciesName, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			log.Fatalf("Error getting species flag: %v", sErr)
		}

		verbosity, vErr := cmd.Flags().GetString("verbosity")
		if vErr != nil {
			log.Fatalf("Error getting species flag: %v", vErr)
		}

		outDir, outErr := cmd.Flags().GetString("out")
		if outErr != nil {
			log.Fatalf("Error getting output directory flag: %v", outErr)
		}

		caller, callerErr := cmd.Flags().GetString("caller")
		if callerErr != nil {
			log.Fatalf("Error getting caller flag: %v", callerErr)
		}

		merger, mergerErr := cmd.Flags().GetString("merger")
		if mergerErr != nil {
			log.Fatalf("Error getting merger flag: %v", mergerErr)
		}

		dvVer, dvVerErr := cmd.Flags().GetString("deepvariant-version")
		if dvVerErr != nil {
			log.Fatalf("Error getting deepvariant version flag: %v", dvVerErr)
		}

		modelType, mtErr := cmd.Flags().GetString("model-type")
		if mtErr != nil {
			log.Fatalf("Error getting model type flag: %v", mtErr)
		}

		if speciesName == "" {
			fmt.Println("Please provide species name with flag --species ")
			return
		}

		//-------------------------------------------- Check dependencies ------------------------------------------ //
		switch caller {
		case "gatk":
			gatkErr := utils.CheckDeps([]string{"gatk"})
			if gatkErr != nil {
				fmt.Println("Dependency check failed ... ", gatkErr)
			}
			if merger == "glnexus" {
				glnErr := utils.CheckDeps([]string{"glnexus_cli", "bcftools", "bgzip"})
				if glnErr != nil {
					fmt.Println("Dependency check failed ... ", glnErr)
				}
			}

		case "deepvariant":
			docErr := utils.CheckDeps([]string{"docker"})
			if docErr != nil {
				fmt.Println("Dependency check failed ... ", docErr)
				return
			}
			if merger == "glnexus" {
				glnErr := utils.CheckDeps([]string{"glnexus_cli", "bcftools", "bgzip"})
				if glnErr != nil {
					fmt.Println("Dependency check failed ... ", glnErr)
					return
				}
			} else {
				fmt.Println("Use glnexus as merger if deepvariant is your chosen variant caller. ... ")
				return
			}
		}

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			_, err := os.Stat(configFile)
			if err != nil {
				fmt.Printf("Config file %s does not exist", configFile)
				return

			}

			verbose, verboseErr := cmd.Flags().GetBool("verbose")
			if verboseErr != nil {
				verbose = false
			}

			variants.VariantCallingConfig(configFile, speciesName, threads, verbosity, caller, merger, dvVer, modelType, verbose)

		} else {
			fmt.Printf("Running without config flag\n")
			bams, bamsErr := cmd.Flags().GetStringSlice("bam")
			if bamsErr != nil {
				log.Fatalf("Error getting bam flag: %v", bamsErr)
			}

			_, rErr := os.Stat(refFile)
			if rErr != nil {
				fmt.Printf("Reference file: %s does not exist\n\n", refFile)
				return
			}

			fmt.Printf("bams: %v\n", bams)
			if len(bams) == 0 {
				fmt.Println("You must provide at least one bam file")
				return
			} else {
				for i, _ := range bams {
					_, err := os.Stat(bams[i])
					if err != nil {
						fmt.Printf("Bam file: %s is not a valid file path", bams[i])
						log.Fatal(err)
					}
				}
			}
			outInfo, outErr := os.Stat(outDir)

			if outErr != nil {

				if os.IsNotExist(outErr) {
					fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", outDir)
					if createErr := os.MkdirAll(outDir, 0755); createErr != nil {
						fmt.Printf("Failed to create output directory %s: %v\n", outDir, createErr)
						return
					}
					fmt.Printf("Output directory %s created successfully.\n", outDir)
				} else {
					fmt.Printf("Error accessing output directory %s: %v\n", outDir, outErr)
					return
				}
			} else if !outInfo.IsDir() {
				fmt.Printf("Output Directory %s file path is not a directory\n", outDir)
				return
			}
			logFilePath := filepath.Join(outDir, "variant_calling.log")
			fmt.Printf("Bams: %v\n", bams)
			fmt.Printf("threads per sample: %v\n", threads)
			fmt.Printf("Reference: %v\n", refFile)

			verbose, verboseErr := cmd.Flags().GetBool("verbose")
			if verboseErr != nil {
				verbose = false
			}

			_, err := variants.VariantCalling(refFile, bams, outDir, speciesName, threads, verbosity, caller, merger, logFilePath, dvVer, modelType, verbose)
			if err != nil {
				return
			}
		}

	},
}

func init() {
	rootCmd.AddCommand(VariantCallingCmd)

	VariantCallingCmd.Flags().StringSliceP("bam", "b", []string{}, "Recalibrated bam file (Can specify multiple)")
	VariantCallingCmd.Flags().StringP("out", "o", "", "Recalibrated bam file")
	VariantCallingCmd.Flags().StringP("species", "s", "", "Species name")
	VariantCallingCmd.Flags().IntP("threads", "t", 4, "Number of threads per sample")
	VariantCallingCmd.Flags().String("verbosity", "ERROR", "Verbosity level for GATK")
	VariantCallingCmd.Flags().BoolP("verbose", "v", false, "Enable verbose output")
	VariantCallingCmd.Flags().StringP("config", "c", "", "Config file")
	VariantCallingCmd.Flags().StringP("reference", "r", "", "Reference file")
	VariantCallingCmd.Flags().String("caller", "gatk", "Variant caller to use. Options: gatk or deepvariant")
	VariantCallingCmd.Flags().StringP("merger", "m", "gatk", "GVCF merger to use. Options: gatk or glnexus")
	VariantCallingCmd.Flags().String("deepvariant-version", "1.9.0", "DeepVariant version")
	VariantCallingCmd.Flags().String("model-type", "WGS", "DeepVariant Model Type: WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA")
}
