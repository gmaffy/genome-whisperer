/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// variantCallingCmd represents the variantCalling command
var variantCallingCmd = &cobra.Command{
	Use:   "variantCalling",
	Short: "Creates a multi-sample vcf from bam files using the GATK best practices",
	Long: `Runs the following pipeline

1. gatk HaplotypeCaller on bam files to form gVCF (if there is more than one bam file) or VCF file(s) for each bam file 
2. gatk GenomeDBImport on all gVCFs to create a datastore
3. gatk GenotypeVcfs to create a multi-sample VCF
4. Hard filter vcfs

`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("variantCalling called")
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
	// variantCallingCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
