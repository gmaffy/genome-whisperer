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
		refFasta, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}
		refVer, vErr := cmd.Flags().GetString("ref-ver")
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

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vErr)
		}

		quick, qErr := cmd.Flags().GetBool("quick")
		if qErr != nil {
			log.Fatalf("Error getting quick flag: %v", qErr)
		}

		skipVerification, sErr := cmd.Flags().GetBool("skip-verification")
		if sErr != nil {
			log.Fatalf("Error getting skip-verification flag: %v", sErr)
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
		fmt.Printf("config: %v, gVcfs: %v, dataDir: %v, species: %v, refVer: %v, refFasta: %v, outDir: %v, merger: %v, log file: %s \n", configFile, gvcfs, dataDir, species, refVer, refFasta, outDir, merger, logFile)
		variants.MergeGvcfs(configFile, gvcfs, dataDir, species, refVer, refFasta, outDir, "gatk", merger, verbose, quick, skipVerification)

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
	MergeGvcfsCmd.Flags().String("ref-ver", "", "reference version")
	MergeGvcfsCmd.Flags().StringP("out-dir", "o", "", "Output directory, (only pass if not using standard data directory structure)")
	MergeGvcfsCmd.Flags().StringP("merger", "m", "glnexus", "GVCF merger to use. Options: glnexus or gatk")
	MergeGvcfsCmd.Flags().StringP("log", "l", "", "log file path")
	MergeGvcfsCmd.Flags().BoolP("verbose", "v", false, "verbose output")
	MergeGvcfsCmd.Flags().BoolP("quick", "q", false, "quick verification")
	MergeGvcfsCmd.Flags().Bool("skip-verification", false, "skip verification")

}
