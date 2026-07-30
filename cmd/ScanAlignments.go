/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/alignmentdir"
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
		result, _ := alignmentdir.ScanAlignments(dataDir, speciesName, refVer, genomesDir, referencePath, verbose, quick)

		fmt.Println(result)

		//for _, bamState := range result {
		//	if bamState.BqsrCram.Valid && bamState.RgmdCram.Valid && bamState.BqsrCram.IndexPresent && bamState.RgmdCram.IndexPresent {
		//		if !bamState.SortedBam.Present || !bamState.RgmdBam.Present || !bamState.BqsrBam.Present || len(bamState.OtherFiles)  == 0 {
		//			color.Green("[%s] PASS - CLEAN")
		//		} else {
		//			color.Cyan("[%s] PASS - DIRTY")
		//		}
		//	} else if bamState.BqsrCram.Valid && !bamState.BqsrCram.IndexPresent {
		//		color.Yellow(["[%s] resume from indexing bqsr"])
		//	}
		//}
	},
}

func init() {
	rootCmd.AddCommand(ScanAlignmentsCmd)
	ScanAlignmentsCmd.Flags().SortFlags = false
	// Every flag this command reads is a persistent flag on the root command.
}
