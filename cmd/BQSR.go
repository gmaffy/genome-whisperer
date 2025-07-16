/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"log"
	"os"

	"github.com/spf13/cobra"
)

// BQSRCmd represents the BQSR command
var BQSRCmd = &cobra.Command{
	Use:   "BQSR",
	Short: "Recalibrates a BAM file using GATK BQSR pipeline.",
	Long: `Runs the following commands on bam files (with duplicates marked)

    1. gatk BaseRecalibrator
    2. gatk ApplyBQSR

If no known-sites file is provided, a bootstrap method of generating one is run`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("BQSR called")
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		refFile, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting ther refFile flag: %v", rErr)
		}

		logFile, logErr := cmd.Flags().GetString("log")
		if logErr != nil {
			log.Fatalf("Error getting log flag: %v", cErr)
		}

		bootstrap, bErr := cmd.Flags().GetBool("bootstrap")
		if bErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bErr)
		}

		jobs, jErr := cmd.Flags().GetInt("jobs")
		if jErr != nil {
			log.Fatalf("Error getting jobs flag: %v", jErr)
		}

		knownSites, ksErr := cmd.Flags().GetStringSlice("known-sites")
		if ksErr != nil {
			log.Fatalf("Error getting known-sites flag: %v", ksErr)
		}

		verbosity, vebErr := cmd.Flags().GetString("verbosity")
		if vebErr != nil {
			log.Fatalf("Error getting known-sites flag: %v", vebErr)
		}

		_, lErr := os.Stat(logFile)
		if lErr != nil {
			fmt.Printf("Log file: %s does not exist", logFile)
			fmt.Println("Provide a valid log file path using --log flag")
			return
		}

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			alignment.BQSRconfig(configFile, bootstrap, jobs, logFile, verbosity)

		} else {
			fmt.Printf("Running without config flag\n")

			bams, bamsErr := cmd.Flags().GetStringSlice("bam")
			if bamsErr != nil {
				log.Fatalf("Error getting bam flag: %v", bamsErr)
			}

			_, rErr := os.Stat(refFile)
			if rErr != nil {
				fmt.Printf("Reference file: %s does not exist", refFile)
				return
			}

			fmt.Printf("bams: %v\n", bams)
			if len(bams) == 0 {
				fmt.Println("You must provide at least one bam file")
				return
			} else {
				for i, _ := range bams {
					_, err := os.Stat(bams[i])
					if err != nil {
						fmt.Printf("Bam file: %s is not a valid file path", bams[i])
						log.Fatal(err)
					}
				}
			}

			if len(knownSites) == 0 && bootstrap == false {
				fmt.Println("Either pass a known-sites file or enable bootstrap method")
				return
			} else if len(knownSites) == 0 && bootstrap == true {
				fmt.Println("Running with bootstrap method")
				err := alignment.BootstrapBqsr(refFile, bams, jobs, logFile, verbosity)
				if err != nil {
					return
				}
			} else if len(knownSites) > 0 {
				fmt.Println("Running with known-sites flag")
				// ------------------------ Checking Known sites file paths ----------------------------------------- //
				for j, _ := range knownSites {
					_, err := os.Stat(knownSites[j])
					if err != nil {
						fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
						log.Fatal(err)
					}
				}

				// --------------------------- Running dbSnpBQSR ---------------------------------------------------- //
				err := alignment.DbSnpBqsr(refFile, bams, knownSites, jobs, logFile)
				if err != nil {
					return
				}

			} else {
				fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
				return
			}

		}

	},
}

func init() {
	rootCmd.AddCommand(BQSRCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// BQSRCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// BQSRCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	BQSRCmd.Flags().StringSliceP("bam", "b", []string{}, "path to bam file (can specify multiple)")
	BQSRCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")
	BQSRCmd.Flags().Bool("bootstrap", false, "Bootstrap method")
	BQSRCmd.Flags().IntP("jobs", "j", 4, "Number of jobs per run")
	BQSRCmd.Flags().String("log", "", "log file path")
	BQSRCmd.Flags().String("config", "", "config file path")
	BQSRCmd.Flags().StringP("reference", "r", "", "path to reference genome")
	BQSRCmd.Flags().String("verbosity", "ERROR", "GATK verbosity level (DEBUG, INFO, WARNING, ERROR)")
}
