/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/spf13/cobra"
	"log"
	"os"
)

// bqsrCmd represents the bqsr command
var bqsrCmd = &cobra.Command{
	Use:   "bqsr",
	Short: "recalibrates a bam file using GATK BQSR pipeline",
	Long: `Runs the following commands on bam files (with duplicates marked)

1. gatk BaseRecalibrator
2. gatk ApplyBQSR

If no known-sites file is provided, a bootstrap method of generating one is run`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("bqsr called")

		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		logFile, logErr := cmd.Flags().GetString("config")
		if logErr != nil {
			log.Fatalf("Error getting log flag: %v", cErr)
		}

		bootstrap, bErr := cmd.Flags().GetBool("bootstrap")
		if bErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bErr)
		}

		jobs, jErr := cmd.Flags().GetInt("jobs")
		if jErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", jErr)
		}

		knownSites, ksErr := cmd.Flags().GetStringSlice("known-sites")
		if ksErr != nil {
			log.Fatalf("Error getting known-sites flag: %v", ksErr)
		}

		_, lErr := os.Stat(logFile)
		if lErr != nil {
			fmt.Printf("Log file: %s does not exist", logFile)
			fmt.Println("Provide a valid log file path using --log flag")
			return
		}

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			alignment.BQSRconfig(configFile, bootstrap, jobs, logFile)

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
				err := alignment.BootstrapBqsr(refFile, bams, jobs, logFile)
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
	rootCmd.AddCommand(bqsrCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// bqsrCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// bqsrCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	bqsrCmd.Flags().StringSliceP("bam", "b", []string{}, "path to bam file (can specify multiple)")
	bqsrCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")
	bqsrCmd.Flags().Bool("bootstrap", false, "Bootstrap method")
	bqsrCmd.Flags().IntP("jobs", "j", 4, "Number of jobs per run")
	bqsrCmd.Flags().String("log", "", "log file path")

}
