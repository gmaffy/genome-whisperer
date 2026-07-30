/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/spf13/cobra"
)

// PrepareFastaFileCmd represents the PrepareFastaFile command
var PrepareFastaFileCmd = &cobra.Command{
	Use:   "PrepareFastaFile",
	Short: "Indexes fasta file with samtools, creates Gatk's fasta dictionary and indexes with aligner",
	Long:  `genome-whisperer PrepareFastaFile -r <reference fasta file> -a <bwa, bwa-mem2, bowtie2 or pbmm2 (default: bwa-mem2)>`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("PrepareFastaFile called")
		refFile, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}

		aligner, aErr := cmd.Flags().GetString("aligner")
		if aErr != nil {
			log.Fatalf("Error getting aligner flag: %v", aErr)
		}

		verbose, aErr := cmd.Flags().GetBool("verbose")
		if aErr != nil {
			log.Fatalf("Error getting verbose flag: %v", aErr)
		}

		err := utils.PrepareFasta(refFile, aligner, verbose)
		if err != nil {

			log.Fatalf("Error: %v", err)
		}
	},
}

func init() {
	rootCmd.AddCommand(PrepareFastaFileCmd)
	PrepareFastaFileCmd.Flags().SortFlags = false
	// --reference, --aligner and --verbose are persistent flags on the root command.
}
