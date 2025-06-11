/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"os"

	"github.com/spf13/cobra"
)

// variantCallingCmd represents the variantCalling command
var variantCallingCmd = &cobra.Command{
	Use:   "variantCalling",
	Short: "Creates a multi-sample VCF file from bam files using GATK best practices",
	Long: `Runs the following pipeline:

1. gatk HaplotypeCaller on bam files to create gVCF files
2. gatk GenomicsDBImport to consolidate gVCFs into a single data store
3. gatk GenotypeVcfs to create a multi-sample VCF
4. Hard filter VCFS`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Printf("Checking dependencies ...\n\n")

		if err := utils.CheckDeps(); err != nil {
			log.Fatalf("Dependency check failed: %v", err)
		}

		fmt.Printf("Dependencies OK\n\n----------------------------------------------------------\n\n")

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
			variants.VariantCallingConfig(configFile, speciesName, jobs, verbosity)

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
			fmt.Printf("Bams: %v\n", bams)
			fmt.Printf("Jobs: %v\n", jobs)
			fmt.Printf("Reference: %v\n", refFile)
			variants.VariantCalling(refFile, bams, outDir, speciesName, jobs, verbosity)
		}
	},
}

func init() {
	rootCmd.AddCommand(variantCallingCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// variantCallingCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	variantCallingCmd.Flags().StringSliceP("bam", "b", []string{}, "Recalibrated bam file (Can specify multiple)")
	variantCallingCmd.Flags().StringP("out", "o", "", "Recalibrated bam file")
	variantCallingCmd.Flags().StringP("species", "s", "", "Species name")
	variantCallingCmd.Flags().IntP("jobs", "j", 4, "Jobs per run")
	variantCallingCmd.Flags().String("verbosity", "WARNING", "Jobs per run")
}
