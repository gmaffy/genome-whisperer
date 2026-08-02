/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"

	"github.com/fatih/color"
	"github.com/spf13/cobra"
)

var VariantCallingCmd = &cobra.Command{
	Use:   "VariantCalling",
	Short: "Creates a multi-sample VCF file from bam files using GATK best practices OR Deepvariant and GLNEXUS",
	Long: `VariantCalling
        - Calls and hard filters SNPs and Indels using GATK best practices from bams generated from short reads
        - Calls and hard filters SNPs and Indels using Deepvariant from bams generated from long reads
        - Can use glnexus or GATK to merge gvcfs
        - Takes input from the standard data directory (--data-dir), a config file
          (--config) or explicit bam paths (--bam)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		color.Green("VariantCalling called")

		opts, err := variantOptions(cmd)
		if err != nil {
			return err
		}

		// A config file supplies whatever was not given on the command line.
		configFile, err := cmd.Flags().GetString("config")
		if err != nil {
			return err
		}
		if configFile != "" {
			cfg, cErr := utils.ReadConfig(configFile)
			if cErr != nil {
				return fmt.Errorf("reading config %s: %w", configFile, cErr)
			}
			if opts.RefFasta == "" {
				opts.RefFasta = cfg.Reference
			}
			if opts.OutDir == "" {
				opts.OutDir = cfg.OutputDir
			}
			if opts.Species == "" {
				opts.Species = cfg.Species
			}
			if len(opts.Bams) == 0 {
				opts.Bams = cfg.Bams
			}
		}

		// ------------------------------------ Stage toggles ------------------------------------- //
		if opts.NoHardFilter, err = cmd.Flags().GetBool("no-hard-filter"); err != nil {
			return err
		}
		// VariantCalling consumes BQSR-recalibrated alignments by default.
		bqsr, err := cmd.Flags().GetBool("bqsr")
		if err != nil {
			return err
		}
		opts.NoBqsr = !bqsr

		// ---------------------------------- Filter thresholds ----------------------------------- //
		lightFilter, _ := cmd.Flags().GetBool("light-filtering")
		minQDSNP, _ := cmd.Flags().GetFloat64("min-qd-snp")
		minQualSNP, _ := cmd.Flags().GetFloat64("min-qual-snp")
		maxSORSNP, _ := cmd.Flags().GetFloat64("max-sor-snp")
		maxFSSNP, _ := cmd.Flags().GetFloat64("max-fs-snp")
		minMQSNP, _ := cmd.Flags().GetFloat64("min-mq-snp")
		minMQRankSNP, _ := cmd.Flags().GetFloat64("min-mqranksum-snp")
		minReadPosRankSNP, _ := cmd.Flags().GetFloat64("min-readposranksum-snp")

		minQDIndel, _ := cmd.Flags().GetFloat64("min-qd-indel")
		minQualIndel, _ := cmd.Flags().GetFloat64("min-qual-indel")
		maxFSIndel, _ := cmd.Flags().GetFloat64("max-fs-indel")
		maxSORIndel, _ := cmd.Flags().GetFloat64("max-sor-indel")
		minReadPosRankIndel, _ := cmd.Flags().GetFloat64("min-readposranksum-indel")

		opts.HardFilter = utils.HardFilterConfig{
			LightFilter:              lightFilter,
			SNP_QD_Min:               minQDSNP,
			SNP_QUAL_Min:             minQualSNP,
			SNP_SOR_Max:              maxSORSNP,
			SNP_FS_Max:               maxFSSNP,
			SNP_MQ_Min:               minMQSNP,
			SNP_MQRankSum_Min:        minMQRankSNP,
			SNP_ReadPosRankSum_Min:   minReadPosRankSNP,
			INDEL_QD_Min:             minQDIndel,
			INDEL_QUAL_Min:           minQualIndel,
			INDEL_FS_Max:             maxFSIndel,
			INDEL_ReadPosRankSum_Min: minReadPosRankIndel,
			INDEL_SOR_Max:            maxSORIndel,
		}
		if opts.MinGQ, err = cmd.Flags().GetInt("min-gq"); err != nil {
			return err
		}

		// ---------------------------------- Check dependencies ---------------------------------- //
		switch opts.Caller {
		case "gatk":
			if depErr := utils.CheckDeps([]string{"gatk", "bcftools", "tabix"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
			if opts.Merger == "glnexus" {
				if depErr := utils.CheckDeps([]string{"glnexus_cli", "bgzip"}); depErr != nil {
					return fmt.Errorf("dependency check failed: %w", depErr)
				}
			} else if opts.Merger != "gatk" {
				return fmt.Errorf("merger must be gatk or glnexus, got %q", opts.Merger)
			}
		case "deepvariant":
			if depErr := utils.CheckDeps([]string{"docker", "bcftools", "tabix"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
			if opts.Merger != "glnexus" {
				return fmt.Errorf("use glnexus as the merger when the caller is deepvariant")
			}
			if depErr := utils.CheckDeps([]string{"glnexus_cli", "bgzip"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
		default:
			return fmt.Errorf("caller must be gatk or deepvariant, got %q", opts.Caller)
		}

		// ------------------------------------- Run pipeline ------------------------------------- //
		finalVCF, err := variants.RunPipeline(opts)
		if err != nil {
			return err
		}
		if finalVCF != "" {
			fmt.Printf("Final VCF: %s\n", finalVCF)
		}
		return nil
	},
}

func init() {
	rootCmd.AddCommand(VariantCallingCmd)
	VariantCallingCmd.Flags().SortFlags = false
	// VariantCalling consumes BQSR-recalibrated alignments by default;
	// --bqsr=false falls back to the RGMD bam/cram.
	VariantCallingCmd.Flags().Bool("bqsr", true, "Run base quality score recalibration")
	VariantCallingCmd.Flags().Bool("no-hard-filter", false, "Do not apply hard filters")

	// ----------------------------------------- Filtering -------------------------------------- //
	VariantCallingCmd.Flags().Bool("light-filtering", false, "Apply the relaxed (light) filter thresholds")
	VariantCallingCmd.Flags().Int("min-gq", 20, "DeepVariant: minimum genotype quality in at least one sample")

	VariantCallingCmd.Flags().Float64("min-qd-snp", 2.0, "SNPs: minimum QualByDepth")
	VariantCallingCmd.Flags().Float64("min-qual-snp", 30.0, "SNPs: minimum QUAL")
	VariantCallingCmd.Flags().Float64("max-sor-snp", 3.0, "SNPs: maximum StrandOddsRatio")
	VariantCallingCmd.Flags().Float64("max-fs-snp", 60.0, "SNPs: maximum FisherStrand")
	VariantCallingCmd.Flags().Float64("min-mq-snp", 40.0, "SNPs: minimum RMSMappingQuality")
	VariantCallingCmd.Flags().Float64("min-mqranksum-snp", -12.5, "SNPs: minimum MappingQualityRankSumTest")
	VariantCallingCmd.Flags().Float64("min-readposranksum-snp", -8.0, "SNPs: minimum ReadPosRankSumTest")

	VariantCallingCmd.Flags().Float64("min-qd-indel", 2.0, "INDELs: minimum QualByDepth")
	VariantCallingCmd.Flags().Float64("min-qual-indel", 30.0, "INDELs: minimum QUAL")
	VariantCallingCmd.Flags().Float64("max-fs-indel", 200.0, "INDELs: maximum FisherStrand")
	VariantCallingCmd.Flags().Float64("max-sor-indel", 10.0, "INDELs: maximum StrandOddsRatio")
	VariantCallingCmd.Flags().Float64("min-readposranksum-indel", -20.0, "INDELs: minimum ReadPosRankSumTest")
}
