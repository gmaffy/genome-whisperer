/*
Copyright © 2025 God <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"os"

	"github.com/spf13/cobra"
)

var alignSrMemCmd = &cobra.Command{
	Use:   "alignSrMem",
	Short: "Short read alignment",
	Long:  `Aligns short PE reads to reference genome using bwa mem.`,
	Run: func(cmd *cobra.Command, args []string) {

		fmt.Printf("Checking dependencies ...\n\n")

		if err := utils.CheckDeps(); err != nil {
			log.Fatalf("Dependency check failed: %v", err)
		}

		fmt.Printf("Dependencies OK\n\n----------------------------------------------------------\n\n")

		// Get all flag values and handle errors
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		referencePath, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			log.Fatalf("Error getting reference flag: %v", rErr)
		}

		forwardPath, fErr := cmd.Flags().GetString("forward")
		if fErr != nil {
			log.Fatalf("Error getting forward reads flag: %v", fErr)
		}

		reversePath, revErr := cmd.Flags().GetString("reverse")
		if cErr != nil {
			log.Fatalf("Error getting reverse read flag: %v", revErr)
		}

		sampleName, sErr := cmd.Flags().GetString("sample")
		if cErr != nil {
			log.Fatalf("Error getting sample name flag: %v", sErr)
		}

		libName, lErr := cmd.Flags().GetString("library")
		if cErr != nil {
			log.Fatalf("Error getting library flag: %v", lErr)
		}

		outDir, oErr := cmd.Flags().GetString("output_dir")
		if cErr != nil {
			log.Fatalf("Error getting output dir flag: %v", oErr)
		}

		threads, tErr := cmd.Flags().GetInt("threads")
		if tErr != nil {
			log.Fatalf("Error getting threads flag: %v", tErr)
		}

		bootstrap, bErr := cmd.Flags().GetBool("bootstrap")
		if tErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bErr)
		}

		bqsr, bqsrErr := cmd.Flags().GetBool("bqsr")
		if tErr != nil {
			log.Fatalf("Error getting bqsr flag: %v", bqsrErr)
		}

		knownSites, ksErr := cmd.Flags().GetStringSlice("known-sites")
		if ksErr != nil {
			log.Fatalf("Error getting bqsr flag: %v", bqsrErr)
		}

		if configFile != "" {
			fmt.Println("Reading config file ...")
			_, confErr := os.Stat(configFile)
			if confErr != nil {
				log.Fatalf("Error reading config file: %v", confErr)
			}
			alignment.AlignShortReadsConfig(configFile, threads, bqsr, bootstrap)
		} else {
			fmt.Println("inline ...")
			_, refErr := os.Stat(referencePath)
			_, fwdErr := os.Stat(forwardPath)
			_, revErr := os.Stat(reversePath)
			outInfo, outErr := os.Stat(outDir)
			if refErr != nil {
				fmt.Printf("Reference genome path: %s, is not valid\n", referencePath)
				return
			}

			if fwdErr != nil {
				fmt.Printf("Forward reads path %s, is not valid\n", forwardPath)
				return
			}

			if revErr != nil {
				fmt.Printf("Reverse reads path %s, is not valid\n", reversePath)
				return
			}

			if outErr != nil {
				fmt.Printf("Output directory: %s is not a valid path\n", outDir)
				return
			}
			if !outInfo.IsDir() {
				fmt.Printf("Output Directory %s file path is not a directory", outDir)
				return
			}
			if sampleName == "" {
				fmt.Println("Please provide sample name is flag -s ")
				return
			}
			if libName == "" {
				fmt.Println("Please provide library name is flag -l ")
				return
			}
			fmt.Printf("All paths PASSED...\n ")
			// ----------------------------------------------- Check Paths if bqsr ------------------------------------------ //
			if bqsr {
				fmt.Println("Skipping BQSR")
				if len(knownSites) == 0 && bootstrap == false {
					fmt.Println("Either pass a known-sites file or enable bootstrap method")
					return
				} else if len(knownSites) > 0 {
					fmt.Println("Running with known-sites flag")
					// ---------------------------- Checking Known sites file paths ----------------------------------------- //
					for j, _ := range knownSites {
						_, err := os.Stat(knownSites[j])
						if err != nil {
							fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
							log.Fatal(err)
						}
					}
					if bootstrap == true {
						fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
						return
					}
				}
			}
			alignment.AlignShortReadsMem(referencePath, forwardPath, reversePath, sampleName, libName, outDir, threads)
		}
	},
}

func init() {
	rootCmd.AddCommand(alignSrMemCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// alignSrMemCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// alignSrMemCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	// --------------------------------------------- Config File ---------------------------------------------------- //
	//alignSrMemCmd.Flags().StringP("config", "c", "", "config file path")
	//alignSrMemCmd.Flags().StringP("reference", "r", "", "Reference genome")
	alignSrMemCmd.Flags().StringP("forward", "1", "", "Path to forward reads")
	alignSrMemCmd.Flags().StringP("reverse", "2", "", "Path to reverse reads")
	alignSrMemCmd.Flags().StringP("sample", "s", "", "Sample name")
	alignSrMemCmd.Flags().StringP("library", "l", "", "Library name")
	alignSrMemCmd.Flags().StringP("output_dir", "o", "", "output directory")
	alignSrMemCmd.Flags().IntP("threads", "t", 8, "number of threads")
	alignSrMemCmd.Flags().Bool("bqsr", false, "verbose")
	alignSrMemCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")
	alignSrMemCmd.Flags().Bool("bootstrap", false, "Bootstrap method")

	//alignSrMemCmd.Flags().Int32P("jobs", "j", 4, "Library name")
}
