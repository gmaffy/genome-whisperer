package cmd

import (
	"os"

	"github.com/gmaffy/genome-whisperer/alignmentdir"
	"github.com/spf13/cobra"
)

// ScanAlignmentsCmd represents the ScanAlignments command
var ScanAlignmentsCmd = &cobra.Command{
	Use:   "ScanAlignments",
	Short: "Report the state of alignment files in the data directory",
	Long: `ScanAlignments walks the standard data directory
  <data-dir>/<species>/<project>/<sample>/reference_genomes/<ref-version>/bams/
and reports, for every sample, which alignment files are present, whether each
is valid and indexed, which other files are lying around, and the steps
AlignReads would still run. It uses the same file rules and the same planning
as AlignReads, so the two never disagree. Nothing is changed.

Each sample gets one status:
  COMPLETE   finished, nothing to do
  CLEANUP    outputs done, only intermediates to remove
  RESUME     an alignment exists, later steps remain
  ALIGN      nothing usable, must align from reads
  AMBIGUOUS  several external alignments (long reads), needs a person
  ERROR      the bams directory could not be read

Plans assume AlignReads --bqsr (the default here); pass --bqsr=false to see the
plan for a run without recalibration. --quick checks headers only; the default
reads each file fully.

The reference comes from --reference or from --genomes-dir.`,
	Example: `  genome-whisperer ScanAlignments --data-dir /data --species squash --ref-version v1.1 \
      --genomes-dir /genomes

  # Header-only checks (faster)
  genome-whisperer ScanAlignments --data-dir /data --species squash --ref-version v1.1 \
      --genomes-dir /genomes --quick`,
	RunE: func(cmd *cobra.Command, args []string) error {
		referencePath, _ := cmd.Flags().GetString("reference")
		dataDir, _ := cmd.Flags().GetString("data-dir")
		species, _ := cmd.Flags().GetString("species")
		refVer, _ := cmd.Flags().GetString("ref-version")
		genomesDir, _ := cmd.Flags().GetString("genomes-dir")
		quick, _ := cmd.Flags().GetBool("quick")
		verbose, _ := cmd.Flags().GetBool("verbose")
		bqsr, _ := cmd.Flags().GetBool("bqsr")

		report, err := alignmentdir.ScanAlignments(dataDir, species, refVer, genomesDir, referencePath, bqsr, verbose, quick)
		if err != nil {
			return err
		}
		cmd.SilenceUsage = true
		alignmentdir.WriteScanReport(os.Stdout, report)
		return nil
	},
}

func init() {
	rootCmd.AddCommand(ScanAlignmentsCmd)
	ScanAlignmentsCmd.Flags().SortFlags = false
	ScanAlignmentsCmd.Flags().Bool("bqsr", true, "Run base quality score recalibration")
}
