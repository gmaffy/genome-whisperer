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
	1. Two bulks only: provide high_bulk (-A) and low_bulk (-B) without parents
	2. Two bulks with one/two parents: provide high_parent (-H), low_parent (-L), high_bulk (-A), and low_bulk (-B)
	3. One bulk with one/two parents: provide high_parent (-H), low_parent (-L), and bulk (-X)`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("GoBSAseq called")
		fmt.Printf("Dependencies OK\n\n----------------------------------------------------------\n\n")

		// Get all flag values and handle errors
		variantFile, varErr := cmd.Flags().GetString("vcf")
		if varErr != nil {
			log.Fatalf("Error getting vcf flag: %v", varErr)
		}

		outputDir, outErr := cmd.Flags().GetString("out")
		if outErr != nil {
			log.Fatalf("Error getting vcf flag: %v", outErr)
		}

		highParent, hpErr := cmd.Flags().GetString("high_parent")
		if hpErr != nil {
			log.Fatalf("Error getting high_parent flag: %v", hpErr)
		}

		lowParent, lpErr := cmd.Flags().GetString("low_parent")
		if lpErr != nil {
			log.Fatalf("Error getting low_parent flag: %v", lpErr)
		}

		highBulk, hbErr := cmd.Flags().GetString("high_bulk")
		if hbErr != nil {
			log.Fatalf("Error getting high_bulk flag: %v", hbErr)
		}

		lowBulk, lbErr := cmd.Flags().GetString("low_bulk")
		if lbErr != nil {
			log.Fatalf("Error getting low_bulk flag: %v", lbErr)
		}

		minHighParentDepth, mhpdErr := cmd.Flags().GetInt("min_high_parent_depth")
		if mhpdErr != nil {
			log.Fatalf("Error getting min_high_parent_depth flag: %v", mhpdErr)
		}

		minLowParentDepth, mlpdErr := cmd.Flags().GetInt("min_low_parent_depth")
		if mlpdErr != nil {
			log.Fatalf("Error getting min_low_parent_depth flag: %v", mlpdErr)
		}

		minHighBulkDepth, mhbErr := cmd.Flags().GetInt("min_high_bulk_depth")
		if mhbErr != nil {
			log.Fatalf("Error getting min_high_bulk_depth flag: %v", mhbErr)
		}

		minLowBulkDepth, mlbErr := cmd.Flags().GetInt("min_low_bulk_depth")
		if mlbErr != nil {
			log.Fatalf("Error getting min_low_bulk_depth flag: %v", mlbErr)
		}

		highBulkSize, hbsErr := cmd.Flags().GetInt("high_bulk_size")
		if hbsErr != nil {
			log.Fatalf("Error getting high_bulk_size flag: %v", hbsErr)
		}

		lowBulkSize, lbsErr := cmd.Flags().GetInt("low_bulk_size")
		if lbsErr != nil {
			log.Fatalf("Error getting low_bulk_size flag: %v", lbsErr)
		}

		windowSize, wsErr := cmd.Flags().GetInt("window_size")
		if wsErr != nil {
			log.Fatalf("Error getting window_size flag: %v", wsErr)
		}

		stepSize, ssErr := cmd.Flags().GetInt("window_step")
		if ssErr != nil {
			log.Fatalf("Error getting window_step flag: %v", ssErr)
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
			log.Fatalf("Error getting bootstrap flag: %v", intErr)
		}

		popStructure, popErr := cmd.Flags().GetString("pop_structure")
		if popErr != nil {
			log.Fatalf("Error getting pop_structure flag: %v", popErr)
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

		refVer, refvErr := cmd.Flags().GetString("refVer")
		if refvErr != nil {
			log.Fatalf("Error getting refVer flag: %v", refvErr)
		}

		//noMerging, noMergingErr := cmd.Flags().GetBool("noMerging")
		//if noMergingErr != nil {
		//	log.Fatalf("Error getting noMerging flag: %v", noMergingErr)
		//}
		noMerging := false

		merger, mergerErr := cmd.Flags().GetString("merger")
		if mergerErr != nil {
			log.Fatalf("Error getting merger flag: %v", mergerErr)

		}

		preset, preErr := cmd.Flags().GetString("preset")
		if preErr != nil {
			log.Fatalf("Error getting output dir flag: %v", preErr)

		}
		dvVer, dvErr := cmd.Flags().GetString("dvVer")
		if dvErr != nil {
			log.Fatalf("Error getting dvVer flag: %v", dvErr)

		}
		modelType, modelErr := cmd.Flags().GetString("modelType")
		if modelErr != nil {
			log.Fatalf("Error getting modelType flag: %v", modelErr)

		}

		gatkLogLevel, gatkErr := cmd.Flags().GetString("gatkLogLevel")
		if gatkErr != nil {
			log.Fatalf("Error getting gatkLogLevel flag: %v", gatkErr)
		}

		alignmentFmt, afErr := cmd.Flags().GetString("alignment-fmt")
		if afErr != nil {
			log.Fatalf("Error getting alignment fmt flag: %v", afErr)
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
					fmt.Println("No output directory provided. Supply path to output directory with -o flag")
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
						"1. Two bulks only: provide high_bulk (-A) and low_bulk (-B) without parents\n" +
						"2. Two bulks with two parents: provide high_parent (-H), low_parent (-L), high_bulk (-A), and low_bulk (-B)\n" +
						"3. One bulk with two parents: provide high_parent (-H), low_parent (-L), and bulk (-X)")
				}
			}
		}
	},
}

func init() {
	rootCmd.AddCommand(GoBSAseqCmd)

	GoBSAseqCmd.Flags().StringP("vcf", "V", "", "Path to target vcf file or vcf table")
	GoBSAseqCmd.Flags().StringP("out", "o", "", "Output directory")
	// ------------------------------------------ PARENTS & BULKS --------------------------------------------------- //
	GoBSAseqCmd.Flags().StringP("high_parent", "H", "", "Name of high parent/resistant parent")
	GoBSAseqCmd.Flags().StringP("low_parent", "L", "", "Name of low parent/susceptible parent")
	GoBSAseqCmd.Flags().StringP("high_bulk", "A", "", "high bulk name (Resistance bulk)")
	GoBSAseqCmd.Flags().StringP("low_bulk", "B", "", "low bulk name (Susceptible bulk)")
	//GoBSAseqCmd.Flags().StringP("bulk", "X", "", "bulk name (if one bulk is used)")

	// --------------------------------------------------- DEPTHS --------------------------------------------------- //
	GoBSAseqCmd.Flags().IntP("min_high_parent_depth", "D", 5, "minimum depth for high parent")
	GoBSAseqCmd.Flags().IntP("min_low_parent_depth", "d", 5, "minimum depth for low parent")
	//GoBSAseqCmd.Flags().IntP("min_bulk_depth", "x", 40, "minimum depth for bulk (if one bulk is used)")
	GoBSAseqCmd.Flags().IntP("min_high_bulk_depth", "a", 40, "minimum depth for high bulk")
	GoBSAseqCmd.Flags().IntP("min_low_bulk_depth", "b", 40, "minimum depth for low bulk")

	// ---------------------------------------------------- SIZES --------------------------------------------------- //

	GoBSAseqCmd.Flags().IntP("bulk_size", "N", 20, "number of individuals in the bulk, (if one bulk is used)")
	GoBSAseqCmd.Flags().IntP("high_bulk_size", "n", 20, "number of individuals in the high bulk")
	GoBSAseqCmd.Flags().IntP("low_bulk_size", "m", 20, "number of individuals in the low bulk")

	// -------------------------------------------- BASIC PARAMS ---------------------------------------------------- //
	GoBSAseqCmd.Flags().StringP("pop_structure", "p", "F2", "F2, BC or RIL")
	GoBSAseqCmd.Flags().Int("rep", 10000, "Replications for threshold calculations ..")

	// -------------------------- Alignment and variant calling ----------------------------------------------------- //
	GoBSAseqCmd.Flags().String("aligner", "bwa-mem2", "bwa-mem, bwa-mem2 bowtie2")
	GoBSAseqCmd.Flags().String("preset", "HIFI", "Dpbmm2 preset. Options: SUBREAD, CSS, HIFI, ISOSEQ and UNROLLED")
	GoBSAseqCmd.Flags().String("caller", "gatk", "Variant caller to use. Options: gatk or deepvariant")
	GoBSAseqCmd.Flags().String("merger", "gatk", "GVCF merger to use. Options: gatk or glnexus")
	GoBSAseqCmd.Flags().String("modelType", "WGS", "DeepVariant Model Type: WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA")
	GoBSAseqCmd.Flags().String("dvVer", "1.9.0", "DeepVariant Version")
	GoBSAseqCmd.Flags().String("gatkLogLevel", "INFO", "GATK log level")

	// ------------------------------------------ PLOTTING PARAMS --------------------------------------------------- //
	GoBSAseqCmd.Flags().IntP("window_size", "w", 2000000, "window size for plotting")
	GoBSAseqCmd.Flags().IntP("window_step", "t", 10000, "step size for plotting")

	// ------------------------------------------------ TOGGLES ----------------------------------------------------- //
	GoBSAseqCmd.Flags().BoolP("smooth", "s", false, "smooth your plot")
	GoBSAseqCmd.Flags().BoolP("interactive", "i", false, "interactive")
	GoBSAseqCmd.Flags().Bool("bootstrap", false, "BSQR bootstrap")
	GoBSAseqCmd.Flags().Bool("bqsr", false, "enable base quality score recalibration")

	//----------------------------------------------- if config ----------------------------------------------------- //
	GoBSAseqCmd.Flags().StringP("config", "c", "", "path to config file")
	GoBSAseqCmd.Flags().Int("threads", 8, "number of threads")
	GoBSAseqCmd.Flags().String("species", "", "number of threads")
	GoBSAseqCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
	GoBSAseqCmd.Flags().String("alignment-fmt", "bam", "bam or cram")
	GoBSAseqCmd.Flags().Bool("noMerging", false, "no merging of gvcf files")
	GoBSAseqCmd.Flags().String("refVer", "", "Ref version")

}
