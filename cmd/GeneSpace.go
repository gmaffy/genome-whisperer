/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/genespace"
	"github.com/spf13/cobra"
)

// GeneSpaceCmd represents the GeneSpace command
var GeneSpaceCmd = &cobra.Command{
	Use:   "GeneSpace",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		gff, _ := cmd.Flags().GetString("gff")
		vcfTable, _ := cmd.Flags().GetString("vcf-table")
		chrom, _ := cmd.Flags().GetString("chrom")
		start, _ := cmd.Flags().GetInt("start")
		stop, _ := cmd.Flags().GetInt("stop")
		resLines, _ := cmd.Flags().GetString("res-lines")
		species, _ := cmd.Flags().GetString("species")
		descFile, _ := cmd.Flags().GetString("desc-file")
		prgFile, _ := cmd.Flags().GetString("prg-file")

		if gff == "" || vcfTable == "" || chrom == "" || resLines == "" {
			fmt.Println("Error: gff, vcf-table, chrom, and res-lines are required")
			cmd.Help()
			return
		}

		fmt.Println("Running GeneSpace analysis ...")
		_, err := genespace.GeneSpace(gff, vcfTable, chrom, start, stop, resLines, species, descFile, prgFile)
		if err != nil {
			fmt.Printf("Error: %v\n", err)
		}
	},
}

func init() {
	rootCmd.AddCommand(GeneSpaceCmd)

	GeneSpaceCmd.Flags().StringP("gff", "g", "", "Path to GFF file (required)")
	GeneSpaceCmd.Flags().StringP("vcf-table", "v", "", "Path to VCF table (required)")
	GeneSpaceCmd.Flags().StringP("chrom", "c", "", "Chromosome name (required)")
	GeneSpaceCmd.Flags().IntP("start", "s", 0, "Start position")
	GeneSpaceCmd.Flags().IntP("stop", "e", 2000000000, "Stop position")
	GeneSpaceCmd.Flags().StringP("res-lines", "r", "", "Comma-separated resistant lines (required)")
	GeneSpaceCmd.Flags().String("species", "", "Species name")
	GeneSpaceCmd.Flags().String("desc-file", "", "Gene description TSV file")
	GeneSpaceCmd.Flags().String("prg-file", "", "PRG blast results file")
}
