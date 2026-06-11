/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"log"
	"os"
	"path/filepath"

	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"

	"github.com/fatih/color"
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

		color.Green("VariantCalling called")
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

		refVer, refVerErr := cmd.Flags().GetString("ref-version")
		if refVerErr != nil {
			log.Fatalf("Error getting reference version flag: %v", refVerErr)
		}

		speciesName, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			log.Fatalf("Error getting species flag: %v", sErr)
		}

		verbosity, vErr := cmd.Flags().GetString("verbosity")
		if vErr != nil {
			log.Fatalf("Error getting verbosity flag: %v", vErr)
		}

		outDir, outErr := cmd.Flags().GetString("out")
		if outErr != nil {
			log.Fatalf("Error getting output directory flag: %v", outErr)
		}

		caller, callerErr := cmd.Flags().GetString("caller")
		if callerErr != nil {
			log.Fatalf("Error getting caller flag: %v", callerErr)
		}

		nomerging, nErr := cmd.Flags().GetBool("no-merging")
		if nErr != nil {
			log.Fatalf("Error getting no-merging flag: %v", nErr)
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

		dataDir, dErr := cmd.Flags().GetString("data-dir")
		if dErr != nil {
			log.Fatalf("Error getting data-dir flag: %v", dErr)
		}

		verbose, verboseErr := cmd.Flags().GetBool("verbose")
		if verboseErr != nil {
			verbose = false
		}

		quick, qErr := cmd.Flags().GetBool("quick")
		if qErr != nil {
			log.Fatalf("Error getting quick flag: %v", qErr)
		}

		if speciesName == "" {
			fmt.Println("Please provide species name with flag --species ")
			return
		}
		genomesDir, gerr := cmd.Flags().GetString("genomes-dir")
		if gerr != nil {
			log.Fatalf("Error getting genomes-dir flag: %v", gerr)
		}

		caller = strings.ToLower(caller)
		merger = strings.ToLower(merger)

		//-------------------------------------------- Check dependencies ------------------------------------------ //
		switch caller {
		case "gatk":
			gatkErr := utils.CheckDeps([]string{"gatk"})
			if gatkErr != nil {
				fmt.Println("Dependency check failed ... ", gatkErr)
				return
			}
			if merger == "glnexus" {
				glnErr := utils.CheckDeps([]string{"glnexus_cli", "bcftools", "bgzip"})
				if glnErr != nil {
					fmt.Println("Dependency check failed ... ", glnErr)
					return
				}
			} else if merger != "gatk" {
				fmt.Println("merger must be either glnexus or gatk")
				return
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
		default:
			fmt.Println("caller must be either gatk or deepvariant")
			return
		}

		//-------------------------------------------- Run variant calling ------------------------------------------ //

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			_, err := os.Stat(configFile)
			if err != nil {
				fmt.Printf("Config file %s does not exist", configFile)
				return

			}

			variants.VariantCallingConfig(configFile, speciesName, threads, verbosity, caller, merger, dvVer, modelType, verbose, nomerging)

		} else if dataDir != "" {
			fmt.Printf("Running Variant Calling using Plennegy data directory structure: %s\n", color.BlueString(dataDir))
			_, err := os.Stat(dataDir)
			if err != nil {
				fmt.Printf("Data directory %s does not exist", dataDir)
				return
			}
			variants.VariantCallingDir(dataDir, speciesName, refVer, genomesDir, refFile, caller, merger, dvVer, modelType, verbose, nomerging, verbosity, quick, threads)
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
				for i := range bams {
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

			_, err := variants.VariantCalling(refFile, bams, outDir, speciesName, threads, verbosity, caller, merger, logFilePath, dvVer, modelType, verbose, nomerging)
			if err != nil {
				return
			}
		}

	},
}

func init() {
	rootCmd.AddCommand(VariantCallingCmd)

	VariantCallingCmd.Flags().StringSliceP("bam", "b", []string{}, "Recalibrated bam file (Can specify multiple)")
	VariantCallingCmd.Flags().StringP("out", "o", "", "Output directory")
	VariantCallingCmd.Flags().StringP("species", "s", "", "Species name")
	VariantCallingCmd.Flags().IntP("threads", "t", 4, "Number of threads per sample")
	VariantCallingCmd.Flags().String("verbosity", "WARNING", "Verbosity level for GATK")
	VariantCallingCmd.Flags().BoolP("verbose", "v", false, "Enable verbose output")
	VariantCallingCmd.Flags().StringP("config", "c", "", "Config file")
	VariantCallingCmd.Flags().StringP("reference", "r", "", "Reference file")
	VariantCallingCmd.Flags().String("ref-version", "", "Reference genome version")
	VariantCallingCmd.Flags().String("caller", "gatk", "Variant caller to use. Options: gatk or DeepVariant")
	VariantCallingCmd.Flags().StringP("merger", "m", "gatk", "GVCF merger to use. Options: gatk or glnexus")
	VariantCallingCmd.Flags().Bool("no-merging", false, "do not merge gvcfs.")
	VariantCallingCmd.Flags().String("deepvariant-version", "1.9.0", "DeepVariant version")
	VariantCallingCmd.Flags().String("model-type", "WGS", "DeepVariant Model Type: WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA")
	VariantCallingCmd.Flags().StringP("data-dir", "d", "", "Main data directory")
	VariantCallingCmd.Flags().Bool("quick", false, "Quick verification")
	VariantCallingCmd.Flags().StringP("genomes-dir", "g", "", "genomes-directory")
}
