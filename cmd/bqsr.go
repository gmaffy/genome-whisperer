/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// bqsrCmd represents the bqsr command
var bqsrCmd = &cobra.Command{
	Use:   "bqsr",
	Short: "recalibrates a bam file using GATK BQSR pipeline",
	Long: `Runs the following commands on bam files (with duplicates marked)

1. gatk BaseRecalibrator
2. gatk ApplyBQSR

If no known-sites file is provided, a bootstrap method of generating one is run`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("bqsr called")
	},
}

func init() {
	rootCmd.AddCommand(bqsrCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// bqsrCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// bqsrCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	//alignSrMemCmd.Flags().StringP("config", "c", "", "config file path")
	//alignSrMemCmd.Flags().StringP("reference", "r", "", "Reference genome")
	bqsrCmd.Flags().StringSliceP("bam", "b", []string{}, "path to bam file (can specify multiple)")
	bqsrCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")

}
