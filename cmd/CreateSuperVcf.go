/*
Copyright © 2026 Godwin Mafireyi (mafireyi@gmail.com)
*/
package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"

	"github.com/fatih/color"
	"github.com/spf13/cobra"
)

// CreateSuperVcfCmd runs the whole variant calling pipeline and then annotates
// the result. It lives in cmd rather than in variants because composing the two
// packages is the CLI layer's job: keeping it here means variants stays about
// producing VCFs and does not grow a dependency on snpEff, gene descriptions and
// PRG parsing.
var CreateSuperVcfCmd = &cobra.Command{
	Use:   "CreateSuperVcf",
	Short: "Run variant calling end to end and annotate the result",
	Long: `CreateSuperVcf
        - Runs the full variant calling pipeline: gVCF creation, merging and hard filtering
        - Annotates the final VCF with snpEff
        - Adds gene descriptions (--gene-description-tsv) and PRG BLAST hits (--prg) when given
        - Takes input from the standard data directory (--data-dir), a config file
          (--config) or explicit bam paths (--bam)

Hard filtering uses the default thresholds for the chosen caller. Use
VariantCalling followed by VariantAnnotation if you need to tune them.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		opts, err := variantOptions(cmd)
		if err != nil {
			return err
		}

		db, err := cmd.Flags().GetString("database")
		if err != nil {
			return err
		}
		if db == "" {
			return fmt.Errorf("a snpEff database is required: pass --database")
		}

		bsaseq, err := cmd.Flags().GetBool("bsaseq")
		if err != nil {
			return err
		}
		desc, err := cmd.Flags().GetString("gene-description-tsv")
		if err != nil {
			return err
		}
		prg, err := cmd.Flags().GetString("prg")
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

		// Default thresholds; this command does not expose the individual knobs.
		opts.HardFilter = utils.HardFilterConfig{
			SNP_QD_Min: 2.0, SNP_QUAL_Min: 30.0, SNP_SOR_Max: 3.0, SNP_FS_Max: 60.0,
			SNP_MQ_Min: 40.0, SNP_MQRankSum_Min: -12.5, SNP_ReadPosRankSum_Min: -8.0,
			INDEL_QD_Min: 2.0, INDEL_QUAL_Min: 30.0, INDEL_FS_Max: 200.0,
			INDEL_ReadPosRankSum_Min: -20.0, INDEL_SOR_Max: 10.0,
		}
		opts.MinGQ = 20

		// ---------------------------------- Check dependencies ---------------------------------- //
		deps := []string{"bcftools", "tabix", "snpEff", "java"}
		switch opts.Caller {
		case "gatk":
			deps = append(deps, "gatk")
		case "deepvariant":
			deps = append(deps, "docker")
			if opts.Merger != "glnexus" {
				return fmt.Errorf("use glnexus as the merger when the caller is deepvariant")
			}
		default:
			return fmt.Errorf("caller must be gatk or deepvariant, got %q", opts.Caller)
		}
		if opts.Merger == "glnexus" {
			deps = append(deps, "glnexus_cli", "bgzip")
		}
		if depErr := utils.CheckDeps(deps); depErr != nil {
			return fmt.Errorf("dependency check failed: %w", depErr)
		}

		// -------------------------------- Variant calling -------------------------------------- //
		finalVCF, err := variants.RunPipeline(opts)
		if err != nil {
			return err
		}
		if finalVCF == "" {
			return fmt.Errorf("no VCF was produced, nothing to annotate")
		}

		// ----------------------------------- Annotation ---------------------------------------- //
		color.Cyan("\n=============================== ANNOTATION ===============================\n\n")

		aErr, superVcfs := annotation.CreateSuperVcf([]string{finalVCF}, db, bsaseq, desc, prg)
		if aErr != nil {
			return fmt.Errorf("annotating %s: %w", finalVCF, aErr)
		}

		color.Green("Super VCF(s) created: %v\n", superVcfs)
		return nil
	},
}

func init() {
	rootCmd.AddCommand(CreateSuperVcfCmd)
	CreateSuperVcfCmd.Flags().SortFlags = false
}
