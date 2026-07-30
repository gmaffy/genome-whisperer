/*
Copyright © 2026 Godwin Mafireyi (mafireyi@gmail.com)
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/utils"
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
			log.Fatalf("Error getting merger flag: %v", mErr)
		}
		caller, callErr := cmd.Flags().GetString("caller")
		if callErr != nil {
			log.Fatalf("Error getting caller flag: %v", callErr)
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
		refFasta, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}
		refVer, vErr := cmd.Flags().GetString("ref-version")
		if vErr != nil {
			log.Fatalf("Error getting reference version flag: %v", vErr)
		}
		outDir, oErr := cmd.Flags().GetString("out-dir")
		if oErr != nil {
			log.Fatalf("Error getting out-dir flag: %v", oErr)
		}

		logFile, lErr := cmd.Flags().GetString("log")
		if lErr != nil {
			log.Fatalf("Error getting log file flag: %v", lErr)
		}

		verbose, vbErr := cmd.Flags().GetBool("verbose")
		if vbErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vbErr)
		}

		quick, qErr := cmd.Flags().GetBool("quick")
		if qErr != nil {
			log.Fatalf("Error getting quick flag: %v", qErr)
		}

		skipVerification, svErr := cmd.Flags().GetBool("skip-verification")
		if svErr != nil {
			log.Fatalf("Error getting skip-verification flag: %v", svErr)
		}

		//---------------------------------------------- Check dependencies ------------------------------------------------------ //
		switch merger {
		case "glnexus":
			glnErr := utils.CheckDeps([]string{"glnexus_cli", "bcftools", "bgzip"})
			if glnErr != nil {
				fmt.Println("Dependency check failed ... ", glnErr)
				return
			}
		case "gatk":
			gatkErr := utils.CheckDeps([]string{"gatk"})
			if gatkErr != nil {
				fmt.Println("Dependency check failed ... ", gatkErr)
				return
			}
		default:
			fmt.Println("merger must be either glnexus or gatk")
			return
		}

		//--------------------------------------------- Run merge gvcfs ---------------------------------------------------------- //
		fmt.Printf("config: %v, gVcfs: %v, dataDir: %v, species: %v, refVer: %v, refFasta: %v, outDir: %v, caller: %v, merger: %v, log file: %s \n", configFile, gvcfs, dataDir, species, refVer, refFasta, outDir, caller, merger, logFile)
		variants.MergeGvcfs(configFile, gvcfs, dataDir, species, refVer, refFasta, outDir, caller, merger, verbose, quick, skipVerification)

	},
}

func init() {
	rootCmd.AddCommand(MergeGvcfsCmd)
	MergeGvcfsCmd.Flags().SortFlags = false
	MergeGvcfsCmd.Flags().StringSlice("gvcf", []string{}, "Path to a gVCF file (repeatable)")
}
