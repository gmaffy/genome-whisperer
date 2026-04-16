/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/spf13/cobra"
)

// ScanAlignmentsCmd represents the ScanAlignments command
var ScanAlignmentsCmd = &cobra.Command{
	Use:   "ScanAlignments",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("ScanAlignments called")
		referencePath, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}

		dataDir, dErr := cmd.Flags().GetString("data-dir")
		if dErr != nil {
			log.Fatalf("Error getting data-dir flag: %v", dErr)
		}
		speciesName, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			log.Fatalf("Error getting species flag: %v", sErr)
		}
		refVer, vErr := cmd.Flags().GetString("ref-version")
		if vErr != nil {
			log.Fatalf("Error getting reference version flag: %v", vErr)
		}

		quick, qErr := cmd.Flags().GetBool("quick")
		if qErr != nil {
			log.Fatalf("Error getting quick flag: %v", qErr)
		}

		genomesDir, gErr := cmd.Flags().GetString("genomes-dir")
		if gErr != nil {
			log.Fatalf("Error getting genomes-dir flag: %v", gErr)
		}

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vErr)
		}
		result := alignment.ScanAlignmentStages(dataDir, speciesName, refVer, genomesDir, referencePath, verbose, quick)
		for _, s := range result.Samples {
			fmt.Printf("%s bqsr present=%v valid=%v err=%v\n",
				s.Sample, s.BqsrBam.Present, s.BqsrBam.Valid, s.BqsrBam.ValidateErr)
		}
	},
}

func init() {
	rootCmd.AddCommand(ScanAlignmentsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// ScanAlignmentsCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	ScanAlignmentsCmd.Flags().StringP("reference", "r", "", "Path to reference genome")
	ScanAlignmentsCmd.Flags().StringP("data-dir", "d", "", "Main data directory")
	ScanAlignmentsCmd.Flags().StringP("species", "S", "", "Species name")
	ScanAlignmentsCmd.Flags().String("ref-version", "", "Reference genome version")
	ScanAlignmentsCmd.Flags().Bool("quick", false, "Quick verification")
	ScanAlignmentsCmd.Flags().StringP("genomes-dir", "g", "", "Skip verification")
	ScanAlignmentsCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
}
