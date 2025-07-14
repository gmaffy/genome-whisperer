/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// alignLrPbCmd represents the alignLrPb command
var alignLrPbCmd = &cobra.Command{
	Use:   "alignLrPb",
	Short: "Align PacBio long read paired-end reads to reference genome using pbmm2.",
	Long:  `Still to come`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("alignLrPb called")
	},
}

func init() {
	rootCmd.AddCommand(alignLrPbCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// alignLrPbCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// alignLrPbCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
