/*
Copyright © 2026 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"strings"

	"github.com/gmaffy/genome-whisperer/markerdesign"
	"github.com/spf13/cobra"
)

// MarkerDesignPipelineCmd represents the MarkerDesignPipeline command
var MarkerDesignPipelineCmd = &cobra.Command{
	Use:   "MarkerDesignPipeline",
	Short: "Design markers in QTLs from a super VCF table",
	Long: `For every QTL (± --flank-size) in --qtls: slice the super VCF table, name sites from
the marker DB, simplify genotypes to one-letter codes, keep sites where resistant and
susceptible lines differ (RES_NOT_SUS, SUS_NOT_RES), and add 50 bp flanks wobbled with
every variant in --gvcf. Writes per-QTL TSV and coloured Excel files to --out-dir.`,
	Run: func(cmd *cobra.Command, args []string) {
		superVcf, _ := cmd.Flags().GetString("variant")
		qtls, _ := cmd.Flags().GetString("qtls")
		flankSize, _ := cmd.Flags().GetInt("flank-size")
		resLinesStr, _ := cmd.Flags().GetString("res-lines")
		susLinesStr, _ := cmd.Flags().GetString("sus-lines")
		reference, _ := cmd.Flags().GetString("reference")
		gvcf, _ := cmd.Flags().GetString("gvcf")
		markerDb, _ := cmd.Flags().GetString("marker-db")
		crop, _ := cmd.Flags().GetString("crop")
		trait, _ := cmd.Flags().GetString("trait")
		outDir, _ := cmd.Flags().GetString("out-dir")
		threads, _ := cmd.Flags().GetInt("threads")

		if reference == "" {
			fmt.Println("Error: --reference is required")
			cmd.Help()
			return
		}
		if outDir == "" {
			outDir = "."
		}

		err := markerdesign.MarkerDesignPipeline(superVcf, qtls, flankSize,
			strings.Split(resLinesStr, ","), strings.Split(susLinesStr, ","),
			reference, gvcf, markerDb, crop, trait, outDir, threads)
		if err != nil {
			fmt.Printf("Error: %v\n", err)
		}
	},
}

func init() {
	rootCmd.AddCommand(MarkerDesignPipelineCmd)
	MarkerDesignPipelineCmd.Flags().SortFlags = false
	MarkerDesignPipelineCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file (super VCF table)")
	MarkerDesignPipelineCmd.Flags().String("qtls", "", "GoBSAseq QTL Excel file (Chrom, Start, Stop)")
	MarkerDesignPipelineCmd.Flags().Int("flank-size", 250000, "Bases added either side of each QTL")
	MarkerDesignPipelineCmd.Flags().String("res-lines", "", "Comma-separated resistant line names")
	MarkerDesignPipelineCmd.Flags().StringP("sus-lines", "u", "", "Comma-separated susceptible line names")
	MarkerDesignPipelineCmd.Flags().String("gvcf", "", "Indexed joint VCF whose variants wobble the flanks")
	MarkerDesignPipelineCmd.Flags().String("marker-db", "", "Marker DB Excel file (chrom, pos, marker_name, typ)")
	MarkerDesignPipelineCmd.Flags().String("crop", "", "Crop prefix for new marker names, e.g. Cmos")
	MarkerDesignPipelineCmd.Flags().String("trait", "", "Trait suffix for new marker names, e.g. pc")
	for _, f := range []string{"variant", "qtls", "res-lines", "sus-lines", "gvcf", "marker-db", "crop", "trait"} {
		if err := MarkerDesignPipelineCmd.MarkFlagRequired(f); err != nil {
			return
		}
	}
}
