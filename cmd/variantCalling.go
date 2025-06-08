/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
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
		fmt.Println("variantCalling called")
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

		outDir, outErr := cmd.Flags().GetString("out")
		if outErr != nil {
			log.Fatalf("Error getting output directory flag: %v", outErr)
		}

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)

		} else {
			fmt.Printf("Running without config flag\n")
			bams, bamsErr := cmd.Flags().GetStringSlice("bam")
			if bamsErr != nil {
				log.Fatalf("Error getting bam flag: %v", bamsErr)
			}

			_, rErr := os.Stat(refFile)
			if rErr != nil {
				fmt.Printf("Reference file: %s does not exist", refFile)
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
			fmt.Printf("Bams: %v\n", bams)
			fmt.Printf("Jobs: %v\n", jobs)
			fmt.Printf("Reference: %v\n", refFile)
			variants.VariantCalling(refFile, bams, outDir, speciesName)
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
	variantCallingCmd.Flags().StringSliceP("bam", "b", []string{}, "Recalibrated bam file")
	variantCallingCmd.Flags().StringP("out", "o", "", "Recalibrated bam file")
	variantCallingCmd.Flags().StringP("species", "s", "", "Species name")
	variantCallingCmd.Flags().IntP("jobs", "j", 4, "Jobs per run")
}
