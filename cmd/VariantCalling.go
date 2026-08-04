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
	Long: `VariantCalling runs the full pipeline: gVCF creation, joint genotyping,
concatenation and hard filtering.

  - Calls and hard filters SNPs and Indels using GATK best practices from bams generated from short reads
  - Calls and hard filters SNPs and Indels using Deepvariant from bams generated from long reads
  - Can use glnexus or GATK to merge gvcfs

INPUT PATHS
  Pick exactly one. Supplying more than one is an error.

  1. Data directory  --data-dir
     Discovers every sample under
       <data-dir>/<species>/<project>/<sample>/reference_genomes/<ref-version>/bams/
     and writes gVCFs beside each sample, in
       .../reference_genomes/<ref-version>/{gatk_gvcfs|dv_gvcfs}/
     VCFs go to
       <data-dir>/<species>/<ref-version>/VCFs/<caller>_<merger>/

       required: --data-dir --species --ref-version
                 and --reference or --genomes-dir
       --out-dir is ignored in this mode.

  2. Config file  --config
     Reads Reference, OutputDir and the bam: entries from the file. Anything also
     given on the command line wins. Generate a template with:
       genome-whisperer CreateTemplate

       required: --config --species --ref-version
                 plus Reference and OutputDir in the file (or --reference / --out-dir)

  3. Inline  --bam
     Every input on the command line. Repeat --bam once per sample.

       required: --bam --reference --out-dir --species --ref-version

  In modes 2 and 3 everything is written under --out-dir:
    <out-dir>/<chrom>/{gatk_gvcfs|dv_gvcfs}/   gVCFs
    <out-dir>/VCFs/<caller>_<merger>/          per-chromosome and final VCFs

CALLER / MERGER
  --caller gatk (default) may use --merger gatk or glnexus.
  --caller deepvariant requires --merger glnexus.

  Joint VCFs are kept in a subdirectory named for the combination that produced
  them: gatk_gatk, gatk_glnexus or dv_glnexus. Running a second combination
  therefore never overwrites the first, and both results stay side by side.

STAGES
  --no-merging      stop after gVCF creation
  --no-hard-filter  stop after merging
  Re-running skips any gVCF or joint VCF that already exists and is valid, so an
  interrupted run resumes where it stopped.`,
	Example: `  # 1. Standard data directory
  genome-whisperer VariantCalling \
      --data-dir /data/plennegy --species cotton --ref-version AD1.1 \
      --genomes-dir /data/genomes --threads 8

  # 2. Config file
  genome-whisperer VariantCalling \
      --config config.txt --species cotton --ref-version AD1.1

  # 3. Inline bams
  genome-whisperer VariantCalling \
      --bam s1.bqsr.cram --bam s2.bqsr.cram \
      --reference /data/genomes/cotton/AD1.1/ref.fa \
      --out-dir ./vcfs --species cotton --ref-version AD1.1

  # Long reads with DeepVariant (GLnexus is required as the merger)
  genome-whisperer VariantCalling \
      --data-dir /data/plennegy --species cotton --ref-version AD1.1 \
      --caller deepvariant --merger glnexus --model-type PACBIO

  # Use the RGMD alignment instead of the BQSR one
  genome-whisperer VariantCalling --data-dir /data/plennegy \
      --species cotton --ref-version AD1.1 --bqsr=false

  # Stop early
  genome-whisperer VariantCalling --data-dir /data/plennegy \
      --species cotton --ref-version AD1.1 --no-merging

  # Relax filtering, or tune one threshold
  genome-whisperer VariantCalling --config config.txt --species cotton \
      --ref-version AD1.1 --light-filtering
  genome-whisperer VariantCalling --config config.txt --species cotton \
      --ref-version AD1.1 --max-fs-snp 40 --min-qual-snp 50`,

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
