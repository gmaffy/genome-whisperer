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

// hardFilterIndelsCmd represents the hardFilterIndels command
var hardFilterIndelsCmd = &cobra.Command{
	Use:   "hardFilterIndels",
	Short: "Hard Filter Indels",
	Long:  `Hard Filter Indels`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("hardFilterIndels called")
		variant, vErr := cmd.Flags().GetString("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		_, err := os.Stat(variant)
		if err != nil {
			fmt.Printf("Variant file %s is not a valid file: %v\n", variant, err)
			return
		}
		variants.HardFilterINDELs(variant)
	},
}

func init() {
	rootCmd.AddCommand(hardFilterIndelsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// hardFilterIndelsCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	hardFilterIndelsCmd.Flags().StringP("variant", "V", "", "INDEL VCF file")
}
