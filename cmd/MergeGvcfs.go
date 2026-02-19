/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/variants"
	"github.com/spf13/cobra"
)

// MergeGvcfsCmd represents the MergeGvcfs command
var MergeGvcfsCmd = &cobra.Command{
	Use:   "MergeGvcfs",
	Short: "Merge gVCFs into a multi-sample VCF",
	Long:  `Merge gvcfs from config file, inline or from standard data directory structure`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("MergeGvcfs called")

		merger, mErr := cmd.Flags().GetString("merger")
		if mErr != nil {
			log.Fatalf("Error getting config flag: %v", mErr)
		}
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}
		gvcfs, gErr := cmd.Flags().GetStringSlice("gvcf")
		if gErr != nil {
			log.Fatalf("Error getting gvcf flag: %v", gErr)
		}

		dataDir, dErr := cmd.Flags().GetString("data-dir")
		if dErr != nil {
			log.Fatalf("Error getting data-dir flag: %v", dErr)
		}

		species, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			log.Fatalf("Error getting species flag: %v", sErr)
		}
		refName, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}
		outDir, oErr := cmd.Flags().GetString("out-dir")
		if oErr != nil {
			log.Fatalf("Error getting out-dir flag: %v", oErr)
		}

		fmt.Printf("config: %v, gVcfs: %v, dataDir: %v, species: %v, refName: %v, outDir: %v, merger: %v \n", configFile, gvcfs, dataDir, species, refName, outDir, merger)
		variants.MergeGvcfs(configFile, gvcfs, dataDir, species, refName, outDir, merger)
	},
}

func init() {
	rootCmd.AddCommand(MergeGvcfsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	MergeGvcfsCmd.Flags().StringP("config", "c", "", "path to config file")
	MergeGvcfsCmd.Flags().StringSliceP("gvcf", "g", []string{}, "path to gvcf file (must be more than 1)")
	MergeGvcfsCmd.Flags().StringP("data-dir", "d", "", "path to data directory")
	MergeGvcfsCmd.Flags().StringP("species", "s", "", "species name")
	MergeGvcfsCmd.Flags().StringP("reference", "r", "", "reference genome name")
	MergeGvcfsCmd.Flags().StringP("out-dir", "o", "", "Output directory, (only pass if not using standard data directory structure)")
	MergeGvcfsCmd.Flags().StringP("merger", "m", "glnexus", "GVCF merger to use. Options: glnexus or gatk")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// MergeGvcfsCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
