/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"strings"

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
		vcfTable, _ := cmd.Flags().GetString("variant")
		chrom, _ := cmd.Flags().GetString("chrom")
		start, _ := cmd.Flags().GetInt("start")
		stop, _ := cmd.Flags().GetInt("stop")
		resLinesStr, _ := cmd.Flags().GetString("res-lines")
		susLinesStr, _ := cmd.Flags().GetString("sus-lines")
		descFile, _ := cmd.Flags().GetString("gene-description-tsv")
		prgFile, _ := cmd.Flags().GetString("prg")
		outDir, _ := cmd.Flags().GetString("out-dir")

		if gff == "" || vcfTable == "" || chrom == "" || resLinesStr == "" || susLinesStr == "" {
			fmt.Println("Error: --gff, --variant, --chrom, --res-lines and --sus-lines are required")
			cmd.Help()
			return
		}

		resLines := strings.Split(resLinesStr, ",")
		susLines := strings.Split(susLinesStr, ",")

		_, err := genespace.GeneSpace(gff, vcfTable, chrom, start, stop, resLines, susLines, descFile, prgFile, outDir)
		if err != nil {
			fmt.Printf("Error: %v\n", err)
		}
	},
}

func init() {
	rootCmd.AddCommand(GeneSpaceCmd)
	GeneSpaceCmd.Flags().SortFlags = false
	GeneSpaceCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")
	GeneSpaceCmd.Flags().String("chrom", "", "Chromosome name")
	GeneSpaceCmd.Flags().Int("start", 0, "Start position")
	GeneSpaceCmd.Flags().IntP("stop", "e", 2000000000, "Stop position")
	GeneSpaceCmd.Flags().String("res-lines", "", "Comma-separated resistant line names")
	GeneSpaceCmd.Flags().StringP("sus-lines", "u", "", "Comma-separated susceptible line names")
	if err := GeneSpaceCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
	if err := GeneSpaceCmd.MarkFlagRequired("chrom"); err != nil {
		return
	}
	if err := GeneSpaceCmd.MarkFlagRequired("res-lines"); err != nil {
		return
	}
	if err := GeneSpaceCmd.MarkFlagRequired("sus-lines"); err != nil {
		return
	}
}
