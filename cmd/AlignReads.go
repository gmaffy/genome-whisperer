/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/spf13/cobra"
	"log"
	"os"
)

// AlignReadsCmd represents the AlignReads command
var AlignReadsCmd = &cobra.Command{
	Use:   "AlignReads",
	Short: "Reads alignment using bwa, bowtie2 or pbmm2.",
	Long: `AlignReads
        - Aligns short paired reads to reference using bwa mem or bowtie2
        - Aligns long reads to reference using pbmm2
        - Marks duplicates using picard tools
        - Recalibrates bam using GATK's BQSR pipeline`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("AlignReads called")
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
		if revErr != nil {
			log.Fatalf("Error getting reverse read flag: %v", revErr)
		}

		sePath, seErr := cmd.Flags().GetString("se")
		if seErr != nil {
			log.Fatalf("Error getting single end read flag: %v", seErr)
		}

		aligner, alErr := cmd.Flags().GetString("aligner")
		if alErr != nil {
			log.Fatalf("Error getting aligner flag: %v", seErr)
		}

		sampleName, sErr := cmd.Flags().GetString("sample")
		if sErr != nil {
			log.Fatalf("Error getting sample name flag: %v", sErr)
		}

		libName, lErr := cmd.Flags().GetString("library")
		if lErr != nil {
			log.Fatalf("Error getting library flag: %v", lErr)
		}

		preset, preErr := cmd.Flags().GetString("preset")
		if preErr != nil {
			log.Fatalf("Error getting output dir flag: %v", preErr)
		}

		outDir, oErr := cmd.Flags().GetString("output_dir")
		if oErr != nil {
			log.Fatalf("Error getting output dir flag: %v", oErr)
		}

		threads, tErr := cmd.Flags().GetInt("threads")
		if tErr != nil {
			log.Fatalf("Error getting threads flag: %v", tErr)
		}

		bootstrap, bErr := cmd.Flags().GetBool("bootstrap")
		if bErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bErr)
		}

		bqsr, bqsrErr := cmd.Flags().GetBool("bqsr")
		if bqsrErr != nil {
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

			logFilePath := outDir + "/alignment.log"

			alignment.RunAlignReadsConfig(configFile, threads, bqsr, bootstrap, aligner, logFilePath)
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
				fmt.Println("Running BQSR")
				if aligner == "pbmm2" {
					fmt.Println("We do not support BQSR for pbmm2 aligner. Please use bwa-mem or bowtie2 aligner or disable BQSR")
					return
				}
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
			logFilePath := outDir + "/alignment.log"
			err := alignment.RunAlignReads(referencePath, forwardPath, reversePath, sePath, sampleName, libName, outDir, threads, aligner, knownSites, bqsr, bootstrap, logFilePath, preset)
			if err != nil {
				return
			}
		}

	},
}

func init() {
	rootCmd.AddCommand(AlignReadsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// AlignReadsCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// AlignReadsCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	AlignReadsCmd.Flags().StringP("forward", "1", "", "Path to forward reads")
	AlignReadsCmd.Flags().StringP("reverse", "2", "", "Path to reverse reads")
	AlignReadsCmd.Flags().String("se", "", "Path to reverse reads")
	AlignReadsCmd.Flags().StringP("sample", "s", "", "Sample name")
	AlignReadsCmd.Flags().StringP("library", "l", "", "Library name")
	AlignReadsCmd.Flags().StringP("output_dir", "o", "", "output directory")
	AlignReadsCmd.Flags().IntP("threads", "t", 8, "number of threads")
	AlignReadsCmd.Flags().Bool("bqsr", false, "perform BQSR")
	AlignReadsCmd.Flags().String("aligner", "bwa-mem", "bwa-mem, bowtie2 or pbmm2")
	AlignReadsCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")
	AlignReadsCmd.Flags().Bool("bootstrap", false, "Bootstrap method")
	AlignReadsCmd.Flags().StringP("reference", "r", "", "Path to reference genome")
	AlignReadsCmd.Flags().StringP("config", "c", "", "Path to reference genome index")
	AlignReadsCmd.Flags().String("preset", "HIFI", "pbmm2 preset. Options: SUBREAD, CSS, HIFI, ISOSEQ and UNROLLED")

}
