package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/warehouse"
	"github.com/spf13/cobra"
)

// ScanWarehouseCmd inventories the configured data roots without modifying them.
var ScanWarehouseCmd = &cobra.Command{
	Use:   "ScanWarehouse",
	Short: "Inventory sequencing storage and write an NDJSON warehouse report",
	RunE: func(cmd *cobra.Command, args []string) error {
		mounts, err := cmd.Flags().GetStringSlice("mounts")
		if err != nil {
			return fmt.Errorf("reading mounts: %w", err)
		}
		report, err := cmd.Flags().GetString("report")
		if err != nil {
			return fmt.Errorf("reading report: %w", err)
		}
		genomesDir, err := cmd.Flags().GetString("genomes-dir")
		if err != nil {
			return fmt.Errorf("reading genomes-dir: %w", err)
		}
		validate, err := cmd.Flags().GetString("validate")
		if err != nil {
			return fmt.Errorf("reading validate: %w", err)
		}
		quick, err := cmd.Flags().GetBool("quick")
		if err != nil {
			return fmt.Errorf("reading quick: %w", err)
		}
		workers, err := cmd.Flags().GetInt("scan-workers")
		if err != nil {
			return fmt.Errorf("reading scan-workers: %w", err)
		}
		emitBatch, err := cmd.Flags().GetBool("emit-batch-seqids")
		if err != nil {
			return fmt.Errorf("reading emit-batch-seqids: %w", err)
		}
		verbose, err := cmd.Flags().GetBool("verbose")
		if err != nil {
			return fmt.Errorf("reading verbose: %w", err)
		}
		if len(mounts) == 0 {
			return fmt.Errorf("at least one --mounts path is required")
		}
		if report == "" {
			return fmt.Errorf("--report is required")
		}
		return warehouse.Scan(warehouse.Options{
			Roots:           mounts,
			GenomesDir:      genomesDir,
			ReportPath:      report,
			Validate:        validate,
			Quick:           quick,
			Workers:         workers,
			EmitBatchSeqIDs: emitBatch,
			Verbose:         verbose,
		})
	},
}

func init() {
	rootCmd.AddCommand(ScanWarehouseCmd)
	ScanWarehouseCmd.Flags().SortFlags = false
	ScanWarehouseCmd.Flags().StringSlice("mounts", nil, "Data roots or mounts to scan (repeatable)")
	ScanWarehouseCmd.Flags().String("report", "", "Output NDJSON report path")
	ScanWarehouseCmd.Flags().String("validate", warehouse.ValidateNone, "Validation scope: none, alignments or all")
	ScanWarehouseCmd.Flags().Int("scan-workers", warehouse.DefaultWorkers, "Concurrent directory inspections")
	ScanWarehouseCmd.Flags().Bool("emit-batch-seqids", false, "Include all sequence IDs in batched contig groups")
}
