/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/variants"
	"os"

	"github.com/spf13/cobra"
)

// HardFilterSNPsCmd represents the HardFilterSNPs command
var HardFilterSNPsCmd = &cobra.Command{
	Use:   "HardFilterSNPs",
	Short: "Hard Filter SNPs",
	Long:  `Hard Filter SNPs`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("HardFilterSNPs called")
		variant, vErr := cmd.Flags().GetString("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		gatkLogLevel, gatErr := cmd.Flags().GetString("gatk-log-level")
		if gatErr != nil {
			fmt.Println("Error getting gatk-log-level flag")
		}

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			fmt.Println("Error getting verbose flag")
		}

		_, err := os.Stat(variant)
		if err != nil {
			fmt.Printf("Variant file %s is not a valid file: %v\n", variant, err)
			return
		}
		_, vaErr := variants.HardFilterSNPs(variant, gatkLogLevel, verbose)
		if vaErr != nil {
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(HardFilterSNPsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// HardFilterSNPsCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// HardFilterSNPsCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	HardFilterSNPsCmd.Flags().StringP("variant", "V", "", "INDEL VCF file")
	HardFilterSNPsCmd.Flags().String("gatk-log-level", "INFO", "GATK log level")
	HardFilterSNPsCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
	err := HardFilterSNPsCmd.MarkFlagRequired("variant")
	if err != nil {
		return
	}
}
