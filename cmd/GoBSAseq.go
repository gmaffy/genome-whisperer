/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmai.com>
*/
package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/gmaffy/genome-whisperer/bsaseq"

	"github.com/spf13/cobra"
)

// GoBSAseqCmd represents the GoBSAseq command
var GoBSAseqCmd = &cobra.Command{
	Use:   "GoBSAseq",
	Short: "Performs BSAseq analysis to detect QTLs",
	Long: `goBSAseq can detect QTLs given reads, bams or vcfs with the following samples:
	1. Two bulks only: provide --high-bulk (-A) and --low-bulk (-B) without parents
	2. Two bulks with one/two parents: provide --high-parent (-H), --low-parent (-L), --high-bulk (-A) and --low-bulk (-B)
	3. One bulk with one/two parents: provide --high-parent (-H), --low-parent (-L) and one of --high-bulk / --low-bulk`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("GoBSAseq called")
		fmt.Printf("Dependencies OK\n\n----------------------------------------------------------\n\n")

		// Get all flag values and handle errors
		variantFile, varErr := cmd.Flags().GetString("variant")
		if varErr != nil {
			log.Fatalf("Error getting variant flag: %v", varErr)
		}

		outputDir, outErr := cmd.Flags().GetString("out-dir")
		if outErr != nil {
			log.Fatalf("Error getting out-dir flag: %v", outErr)
		}

		highParent, hpErr := cmd.Flags().GetString("high-parent")
		if hpErr != nil {
			log.Fatalf("Error getting high-parent flag: %v", hpErr)
		}

		lowParent, lpErr := cmd.Flags().GetString("low-parent")
		if lpErr != nil {
			log.Fatalf("Error getting low-parent flag: %v", lpErr)
		}

		highBulk, hbErr := cmd.Flags().GetString("high-bulk")
		if hbErr != nil {
			log.Fatalf("Error getting high-bulk flag: %v", hbErr)
		}

		lowBulk, lbErr := cmd.Flags().GetString("low-bulk")
		if lbErr != nil {
			log.Fatalf("Error getting low-bulk flag: %v", lbErr)
		}

		minHighParentDepth, mhpdErr := cmd.Flags().GetInt("min-high-parent-depth")
		if mhpdErr != nil {
			log.Fatalf("Error getting min-high-parent-depth flag: %v", mhpdErr)
		}

		minLowParentDepth, mlpdErr := cmd.Flags().GetInt("min-low-parent-depth")
		if mlpdErr != nil {
			log.Fatalf("Error getting min-low-parent-depth flag: %v", mlpdErr)
		}

		minHighBulkDepth, mhbErr := cmd.Flags().GetInt("min-high-bulk-depth")
		if mhbErr != nil {
			log.Fatalf("Error getting min-high-bulk-depth flag: %v", mhbErr)
		}

		minLowBulkDepth, mlbErr := cmd.Flags().GetInt("min-low-bulk-depth")
		if mlbErr != nil {
			log.Fatalf("Error getting min-low-bulk-depth flag: %v", mlbErr)
		}

		highBulkSize, hbsErr := cmd.Flags().GetInt("high-bulk-size")
		if hbsErr != nil {
			log.Fatalf("Error getting high-bulk-size flag: %v", hbsErr)
		}

		lowBulkSize, lbsErr := cmd.Flags().GetInt("low-bulk-size")
		if lbsErr != nil {
			log.Fatalf("Error getting low-bulk-size flag: %v", lbsErr)
		}

		windowSize, wsErr := cmd.Flags().GetInt("window-size")
		if wsErr != nil {
			log.Fatalf("Error getting window-size flag: %v", wsErr)
		}

		stepSize, ssErr := cmd.Flags().GetInt("window-step")
		if ssErr != nil {
			log.Fatalf("Error getting window-step flag: %v", ssErr)
		}

		smoothing, smErr := cmd.Flags().GetBool("smooth")
		if smErr != nil {
			log.Fatalf("Error getting smooth flag: %v", smErr)
		}

		interactive, intErr := cmd.Flags().GetBool("interactive")
		if intErr != nil {
			log.Fatalf("Error getting interactive flag: %v", intErr)
		}

		bootstrap, bootErr := cmd.Flags().GetBool("bootstrap")
		if bootErr != nil {
			log.Fatalf("Error getting bootstrap flag: %v", bootErr)
		}

		popStructure, popErr := cmd.Flags().GetString("pop-structure")
		if popErr != nil {
			log.Fatalf("Error getting pop-structure flag: %v", popErr)
		}

		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		rep, repErr := cmd.Flags().GetInt("rep")
		if repErr != nil {
			log.Fatalf("Error getting rep flag: %v", repErr)
		}

		threads, tErr := cmd.Flags().GetInt("threads")
		if tErr != nil {
			log.Fatalf("Error getting threads flag: %v", tErr)
		}

		species, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			log.Fatalf("Error getting species flag: %v", sErr)
		}

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			log.Fatalf("Error getting verbose flag: %v", vErr)
		}

		bqsr, bqsrErr := cmd.Flags().GetBool("bqsr")
		if bqsrErr != nil {
			log.Fatalf("Error getting bqsr flag: %v", bqsrErr)
		}

		aligner, aErr := cmd.Flags().GetString("aligner")
		if aErr != nil {
			log.Fatalf("Error getting aligner flag: %v", aErr)
		}

		caller, callerErr := cmd.Flags().GetString("caller")
		if callerErr != nil {
			log.Fatalf("Error getting caller flag: %v", callerErr)
		}

		refVer, refvErr := cmd.Flags().GetString("ref-version")
		if refvErr != nil {
			log.Fatalf("Error getting ref-version flag: %v", refvErr)
		}

		noMerging, noMergingErr := cmd.Flags().GetBool("no-merging")
		if noMergingErr != nil {
			log.Fatalf("Error getting no-merging flag: %v", noMergingErr)
		}

		merger, mergerErr := cmd.Flags().GetString("merger")
		if mergerErr != nil {
			log.Fatalf("Error getting merger flag: %v", mergerErr)
		}

		preset, preErr := cmd.Flags().GetString("preset")
		if preErr != nil {
			log.Fatalf("Error getting preset flag: %v", preErr)
		}

		dvVer, dvErr := cmd.Flags().GetString("deepvariant-version")
		if dvErr != nil {
			log.Fatalf("Error getting deepvariant-version flag: %v", dvErr)
		}

		modelType, modelErr := cmd.Flags().GetString("model-type")
		if modelErr != nil {
			log.Fatalf("Error getting model-type flag: %v", modelErr)
		}

		gatkLogLevel, gatkErr := cmd.Flags().GetString("gatk-log-level")
		if gatkErr != nil {
			log.Fatalf("Error getting gatk-log-level flag: %v", gatkErr)
		}

		alignmentFmt, afErr := cmd.Flags().GetString("alignment-fmt")
		if afErr != nil {
			log.Fatalf("Error getting alignment-fmt flag: %v", afErr)
		}

		if interactive {
			fmt.Println("Running in interactive mode")
			bsaseq.InteractiveRun(variantFile, popStructure, rep, outputDir)

		} else {

			if variantFile == "" {

				fmt.Println("Reading config file ...")
				_, confErr := os.Stat(configFile)
				if confErr != nil {
					log.Fatalf("Error reading config file: %v", confErr)
				}

				if species == "" {
					log.Fatal("Please provide a species name")
				}
				fmt.Println("Running from config file")
				bsaseq.RunBsaSeqFromConfig(configFile, threads, species, minHighParentDepth, minLowParentDepth, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, bootstrap, bqsr, caller, merger, noMerging, aligner, preset, dvVer, modelType, gatkLogLevel, verbose, alignmentFmt, refVer)
			} else {
				if outputDir == "" {
					fmt.Println("No output directory provided. Supply one with --out-dir")
					return
				}
				outInfo, outErr := os.Stat(outputDir)

				if outErr != nil {

					if os.IsNotExist(outErr) {
						fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", outputDir)
						if createErr := os.MkdirAll(outputDir, 0755); createErr != nil {
							fmt.Printf("Failed to create output directory %s: %v\n", outputDir, createErr)
							return
						}
						fmt.Printf("Output directory %s created successfully.\n", outputDir)
					} else {
						fmt.Printf("Error accessing output directory %s: %v\n", outputDir, outErr)
						return
					}
				} else if !outInfo.IsDir() {
					fmt.Printf("Output Directory %s file path is not a directory\n", outputDir)
					return
				}
				fmt.Println("Running from vcf file")
				if highParent == "" && lowParent == "" && highBulk != "" && lowBulk != "" {
					fmt.Println("Running 2 bulks only analysis")
					bsaseq.TwoBulkOnlyRun(variantFile, highBulk, lowBulk, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputDir)
				} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk != "" {
					fmt.Println("Running 2 bulks 2 parents analysis")
					bsaseq.TwoBulkTwoParentsRun(variantFile, highParent, lowParent, highBulk, lowBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputDir)

				} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk == "" {
					fmt.Println("Running 1 high bulk, 2 parent analysis")
					outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_high_bsaseq_stats.tsv"
					bsaseq.OneBulkTwoParentsRun(variantFile, highParent, lowParent, highBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, highBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName, outputDir)

				} else if highParent != "" && lowParent != "" && highBulk == "" && lowBulk != "" {
					fmt.Println("Running 1 low bulk, 2 parent analysis")
					outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_low_bsaseq_stats.tsv"
					bsaseq.OneBulkTwoParentsRun(variantFile, highParent, lowParent, lowBulk, minHighParentDepth, minLowParentDepth, minLowBulkDepth, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName, outputDir)

				} else {
					log.Fatal("Invalid parameters. Valid combinations are:\n" +
						"1. Two bulks only: --high-bulk (-A) and --low-bulk (-B) without parents\n" +
						"2. Two bulks with two parents: --high-parent (-H), --low-parent (-L), --high-bulk (-A) and --low-bulk (-B)\n" +
						"3. One bulk with two parents: --high-parent (-H), --low-parent (-L) and one of --high-bulk / --low-bulk")
				}
			}
		}
	},
}

func init() {
	rootCmd.AddCommand(GoBSAseqCmd)
	GoBSAseqCmd.Flags().SortFlags = false
	GoBSAseqCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")

	// ----------------------------------- Parents and bulks ----------------------------------- //
	GoBSAseqCmd.Flags().StringP("high-parent", "H", "", "Name of the high (resistant) parent")
	GoBSAseqCmd.Flags().StringP("low-parent", "L", "", "Name of the low (susceptible) parent")
	GoBSAseqCmd.Flags().StringP("high-bulk", "A", "", "Name of the high (resistant) bulk")
	GoBSAseqCmd.Flags().StringP("low-bulk", "B", "", "Name of the low (susceptible) bulk")

	// ------------------------------------------ Depths ---------------------------------------- //
	GoBSAseqCmd.Flags().Int("min-high-parent-depth", 5, "Minimum depth for the high parent")
	GoBSAseqCmd.Flags().Int("min-low-parent-depth", 5, "Minimum depth for the low parent")
	GoBSAseqCmd.Flags().Int("min-high-bulk-depth", 40, "Minimum depth for the high bulk")
	GoBSAseqCmd.Flags().Int("min-low-bulk-depth", 40, "Minimum depth for the low bulk")

	// ------------------------------------------- Sizes ---------------------------------------- //
	GoBSAseqCmd.Flags().IntP("high-bulk-size", "n", 20, "Number of individuals in the high bulk")
	GoBSAseqCmd.Flags().IntP("low-bulk-size", "m", 20, "Number of individuals in the low bulk")

	// ---------------------------------------- Basic params ------------------------------------ //
	GoBSAseqCmd.Flags().StringP("pop-structure", "p", "F2", "Population structure: F2, BC or RIL")
	GoBSAseqCmd.Flags().Int("rep", 10000, "Replications for threshold calculation")

	// ------------------------------------------ Plotting -------------------------------------- //
	GoBSAseqCmd.Flags().IntP("window-size", "w", 2000000, "Plot window size")
	GoBSAseqCmd.Flags().Int("window-step", 10000, "Plot window step size")

	// ------------------------------------------ Toggles --------------------------------------- //
	GoBSAseqCmd.Flags().Bool("smooth", false, "Smooth the plot")
	GoBSAseqCmd.Flags().BoolP("interactive", "i", false, "Run in interactive mode")
	GoBSAseqCmd.Flags().Bool("bqsr", false, "Run base quality score recalibration")
}
