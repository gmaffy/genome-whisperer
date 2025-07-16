/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
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

		verbose, vebErr := cmd.Flags().GetBool("verbose")
		if vebErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vebErr)
		}

		knownSites, ksErr := cmd.Flags().GetStringSlice("known-sites")
		if ksErr != nil {
			log.Fatalf("Error getting bqsr flag: %v", bqsrErr)
		}

		gatkLogLevel, glErr := cmd.Flags().GetString("gatk-log-evel")
		if glErr != nil {
			log.Fatalf("Error getting bqsr flag: %v", glErr)
		}

		if configFile != "" {
			fmt.Println("Reading config file ...")
			_, confErr := os.Stat(configFile)
			if confErr != nil {
				log.Fatalf("Error reading config file: %v", confErr)
			}

			//logFilePath := outDir + "/alignment.log"

			bams, err := alignment.RunAlignReadsConfig(configFile, threads, bqsr, bootstrap, aligner, preset, gatkLogLevel, verbose)
			if err != nil {
				fmt.Println(err)
				return
			} else {
				fmt.Printf("%v BAMs generated", len(bams))
			}
		} else {
			fmt.Println("inline ...")
			_, refErr := os.Stat(referencePath)
			_, fwdErr := os.Stat(forwardPath)
			_, revErr := os.Stat(reversePath)

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

			outInfo, outErr := os.Stat(outDir)

			if outErr != nil {

				if os.IsNotExist(outErr) {
					fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", outDir)
					if createErr := os.MkdirAll(outDir, 0755); createErr != nil {
						fmt.Printf("Failed to create output directory %s: %v\n", outDir, createErr)
						return
					}
					fmt.Printf("Output directory %s created successfully.\n", outDir)
				} else {
					fmt.Printf("Error accessing output directory %s: %v\n", outDir, outErr)
					return
				}
			} else if !outInfo.IsDir() {
				fmt.Printf("Output Directory %s file path is not a directory\n", outDir)
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
			if aligner == "pbmm2" {
				depErr := utils.CheckDeps([]string{"pbmm2", "samtools", "gatk"})
				if depErr != nil {
					log.Fatalf("Dependency check failed: %v", depErr)
				}

			} else if aligner == "bwa-mem" {
				depErr := utils.CheckDeps([]string{"bwa", "samtools", "gatk"})
				if depErr != nil {
					log.Fatalf("Dependency check failed: %v", depErr)
				}
			} else if aligner == "bowtie2" {
				depErr := utils.CheckDeps([]string{"bowtie2", "samtools", "gatk"})
				if depErr != nil {
					log.Fatalf("Dependency check failed: %v", depErr)
				}
			}
			bam, err := alignment.RunAlignReads(referencePath, forwardPath, reversePath, sePath, sampleName, libName, outDir, threads, aligner, knownSites, bqsr, bootstrap, logFilePath, preset, gatkLogLevel, verbose)
			if err != nil {
				return
			} else {
				fmt.Printf("BAM file generated: %s\n", bam)
				return
			}
		}

	},
}

func init() {
	AlignReadsCmd.Flags().StringP("forward", "1", "", "Path to forward reads")
	AlignReadsCmd.Flags().StringP("reverse", "2", "", "Path to reverse reads")
	AlignReadsCmd.Flags().String("se", "", "Path to reverse reads")
	AlignReadsCmd.Flags().StringP("sample", "s", "", "Sample name")
	AlignReadsCmd.Flags().StringP("library", "l", "", "Library name")
	AlignReadsCmd.Flags().StringP("output_dir", "o", "", "output directory")
	AlignReadsCmd.Flags().IntP("threads", "t", 8, "number of threads")
	AlignReadsCmd.Flags().Bool("bqsr", false, "perform BQSR")
	AlignReadsCmd.Flags().Bool("verbose", false, "verbose mode")
	AlignReadsCmd.Flags().String("aligner", "bwa-mem", "bwa-mem, bowtie2 or pbmm2")
	AlignReadsCmd.Flags().StringSliceP("known-sites", "k", []string{}, "Path to known sites vcf (can specify multiple)")
	AlignReadsCmd.Flags().Bool("bootstrap", false, "Bootstrap method")
	AlignReadsCmd.Flags().StringP("reference", "r", "", "Path to reference genome")
	AlignReadsCmd.Flags().StringP("config", "c", "", "Path to reference genome index")
	AlignReadsCmd.Flags().String("preset", "HIFI", "pbmm2 preset. Options: SUBREAD, CSS, HIFI, ISOSEQ and UNROLLED")
	AlignReadsCmd.Flags().String("gatk-log-evel", "ERROR", "GATK log level. Options: ERROR, INFO, DEBUG, TRACE")

}
