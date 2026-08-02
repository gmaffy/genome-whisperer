/*
Copyright © 2026 Godwin Mafireyi (mafireyi@gmail.com)
*/
package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"github.com/spf13/cobra"
)

// MergeGvcfsCmd represents the MergeGvcfs command
var MergeGvcfsCmd = &cobra.Command{
	Use:   "MergeGvcfs",
	Short: "Merge gVCFs into a multi-sample VCF",
	Long: `MergeGvcfs
        - Joint genotypes each chromosome with GATK or GLnexus
        - Concatenates the per-chromosome VCFs into one multi-sample VCF
        - Finds gVCFs in the standard data directory (--data-dir), or takes them
          from a config file (--config) or explicit paths (--gvcf)
        - Reuses a joint VCF that is valid and already holds the expected samples`,
	RunE: func(cmd *cobra.Command, args []string) error {
		opts, err := variantOptions(cmd)
		if err != nil {
			return err
		}

		gvcfs, err := cmd.Flags().GetStringSlice("gvcf")
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
			if len(gvcfs) == 0 {
				gvcfs = cfg.GVCFs
			}
		}

		// ---------------------------------- Check dependencies ---------------------------------- //
		switch opts.Merger {
		case "glnexus":
			if depErr := utils.CheckDeps([]string{"glnexus_cli", "bcftools", "bgzip", "tabix", "gatk"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
		case "gatk":
			if depErr := utils.CheckDeps([]string{"gatk", "bcftools", "tabix"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
		default:
			return fmt.Errorf("merger must be gatk or glnexus, got %q", opts.Merger)
		}

		// Explicit gVCF paths are merged as a single group. Without a chromosome to
		// key them on there is nothing to split them by, and a caller passing paths
		// by hand is merging one set.
		var grouped map[string][]string
		if len(gvcfs) > 0 {
			grouped = map[string][]string{"all": gvcfs}
		}

		finalVCF, err := variants.MergeGvcfs(opts, grouped)
		if err != nil {
			return err
		}
		fmt.Printf("Multi-sample VCF: %s\n", finalVCF)
		return nil
	},
}

func init() {
	rootCmd.AddCommand(MergeGvcfsCmd)
	MergeGvcfsCmd.Flags().SortFlags = false
	MergeGvcfsCmd.Flags().StringSlice("gvcf", []string{}, "Path to a gVCF file (repeatable)")
}
