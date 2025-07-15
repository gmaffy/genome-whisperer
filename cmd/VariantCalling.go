/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"
)

// VariantCallingCmd represents the VariantCalling command
var VariantCallingCmd = &cobra.Command{
	Use:   "VariantCalling",
	Short: "Creates a multi-sample VCF file from bam files using GATK best practices OR Deepvariant and GLNEXUS",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("VariantCalling called")
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		jobs, jErr := cmd.Flags().GetInt("jobs")
		if jErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", jErr)
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

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			_, err := os.Stat(configFile)
			if err != nil {
				fmt.Printf("Config file %s does not exist", configFile)
				return

			}
			variants.VariantCallingConfig(configFile, speciesName, jobs, verbosity, caller, merger, dvVer, modelType)

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
			fmt.Printf("Jobs: %v\n", jobs)
			fmt.Printf("Reference: %v\n", refFile)
			variants.VariantCalling(refFile, bams, outDir, speciesName, jobs, verbosity, caller, merger, logFilePath, dvVer, modelType)
		}

	},
}

func init() {
	rootCmd.AddCommand(VariantCallingCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// VariantCallingCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// VariantCallingCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	VariantCallingCmd.Flags().StringSliceP("bam", "b", []string{}, "Recalibrated bam file (Can specify multiple)")
	VariantCallingCmd.Flags().StringP("out", "o", "", "Recalibrated bam file")
	VariantCallingCmd.Flags().StringP("species", "s", "", "Species name")
	VariantCallingCmd.Flags().IntP("jobs", "j", 4, "Jobs per run")
	VariantCallingCmd.Flags().String("verbosity", "INFO", "Jobs per run")
	VariantCallingCmd.Flags().StringP("config", "c", "", "Config file")
	VariantCallingCmd.Flags().StringP("reference", "r", "", "Reference file")
	VariantCallingCmd.Flags().String("caller", "gatk", "Variant caller to use. Options: gatk or deepvariant")
	VariantCallingCmd.Flags().StringP("merger", "m", "gatk", "GVCF merger to use. Options: gatk or glnexus")
	VariantCallingCmd.Flags().String("deepvariant-version", "1.9.0", "DeepVariant version")
	VariantCallingCmd.Flags().String("model-type", "WGS", "DeepVariant Model Type: WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA")
}
