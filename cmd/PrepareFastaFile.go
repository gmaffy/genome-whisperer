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

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// PrepareFastaFileCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	PrepareFastaFileCmd.Flags().StringP("reference", "r", "", "Reference fasta file ...")
	PrepareFastaFileCmd.Flags().StringP("aligner", "a", "bwa-mem2", "bwa, bwa-mem2, bowtie2 or pbmm2 ...")
	PrepareFastaFileCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
}
