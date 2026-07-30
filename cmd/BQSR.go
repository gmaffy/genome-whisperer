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
			log.Fatalf("Error getting log flag: %v", logErr)
		}

		bootstrap, bErr := cmd.Flags().GetBool("bootstrap")
		if bErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bErr)
		}

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vErr)
		}

		threads, tErr := cmd.Flags().GetInt("threads")
		if tErr != nil {
			log.Fatalf("Error getting threads flag: %v", tErr)
		}

		knownSites, ksErr := cmd.Flags().GetStringSlice("known-sites")
		if ksErr != nil {
			log.Fatalf("Error getting known-sites flag: %v", ksErr)
		}

		gatkLogLevel, gtErr := cmd.Flags().GetString("gatk-log-level")
		if gtErr != nil {
			log.Fatalf("Error getting gatk log level flag: %v", gtErr)
		}

		_, lErr := os.Stat(logFile)
		if lErr != nil {
			fmt.Printf("Log file: %s does not exist", logFile)
			fmt.Println("Provide a valid log file path using --log flag")
			return
		}

		if configFile != "" {
			fmt.Printf("Running with config file to %s\n", configFile)
			alignment.BQSRconfig(configFile, bootstrap, threads, logFile, gatkLogLevel, verbose)

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
				for i := range bams {
					_, err := os.Stat(bams[i])
					if err != nil {
						fmt.Printf("Bam file: %s is not a valid file path", bams[i])
						log.Fatal(err)
					}
				}
			}

			if len(knownSites) == 0 && !bootstrap {
				fmt.Println("Either pass a known-sites file or enable bootstrap method")
				return
			} else if len(knownSites) == 0 && bootstrap {
				fmt.Println("Running with bootstrap method")
				err := alignment.BootstrapBqsr(refFile, bams, threads, logFile, gatkLogLevel, verbose)
				if err != nil {
					return
				}
			} else if len(knownSites) > 0 {
				fmt.Println("Running with known-sites flag")
				// ------------------------ Checking Known sites file paths ----------------------------------------- //
				for j := range knownSites {
					_, err := os.Stat(knownSites[j])
					if err != nil {
						fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
						log.Fatal(err)
					}
				}

				// --------------------------- Running dbSnpBQSR ---------------------------------------------------- //
				err := alignment.DbSnpBqsr(refFile, bams, knownSites, threads, logFile)
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
	BQSRCmd.Flags().SortFlags = false
	// Every flag this command reads is a persistent flag on the root command.
}
