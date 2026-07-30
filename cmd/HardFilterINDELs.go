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

// HardFilterINDELsCmd represents the HardFilterINDELs command
var HardFilterINDELsCmd = &cobra.Command{
	Use:   "HardFilterINDELs",
	Short: "Hard Filter Indels",
	Long:  `Hard Filter Indels`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("HardFilterINDELs called")
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
		_, vaErr := variants.HardFilterINDELs(variant, gatkLogLevel, verbose)
		if vaErr != nil {
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(HardFilterINDELsCmd)
	HardFilterINDELsCmd.Flags().SortFlags = false
	HardFilterINDELsCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")
	if err := HardFilterINDELsCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
}
