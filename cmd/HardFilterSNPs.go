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
	HardFilterSNPsCmd.Flags().SortFlags = false
	HardFilterSNPsCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")
	if err := HardFilterSNPsCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
}
