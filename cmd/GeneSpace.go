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
	Short: "Rank candidate genes in a QTL window",
	Long: `GeneSpace ranks the genes in a chromosome window as candidates for a trait.

Each gene is scored on two independent axes and given a tier:

  Variant score  (0-4)  the most damaging variant whose genotypes follow the
                        resistant/susceptible split, tested under both a
                        dominant and a recessive pattern
  Function score (0-3)  PRG hit quality and gene-description terms from a
                        built-in trait profile

  A lead candidate · B strong · C plausible · D positional only

Input is a superVCF table (from CreateSuperVcf / VariantAnnotation or GoBSAseq).
Output is GeneSpace_<chrom>_<start>-<stop>.xlsx in --out-dir.

Example:
  genome-whisperer GeneSpace --gff genes.gff3 -V cohort_SUPER_VCF.tsv \
    --chrom Chr04 --start 1200000 --stop 2400000 \
    --res-lines R1,R2,R3 --sus-lines S1,S2 \
    --prg SPECIES_PRG_TOP_HITS.txt --gene-description-tsv desc.tsv -o genespace`,
	RunE: func(cmd *cobra.Command, args []string) error {
		gff, _ := cmd.Flags().GetString("gff")
		resLinesStr, _ := cmd.Flags().GetString("res-lines")
		susLinesStr, _ := cmd.Flags().GetString("sus-lines")
		opts := genespace.Options{GFF: gff}
		opts.VariantTable, _ = cmd.Flags().GetString("variant")
		opts.Chrom, _ = cmd.Flags().GetString("chrom")
		opts.Start, _ = cmd.Flags().GetInt("start")
		opts.Stop, _ = cmd.Flags().GetInt("stop")
		opts.DescFile, _ = cmd.Flags().GetString("gene-description-tsv")
		opts.PrgFile, _ = cmd.Flags().GetString("prg")
		opts.OutDir, _ = cmd.Flags().GetString("out-dir")
		opts.Trait, _ = cmd.Flags().GetString("trait")
		opts.MinConcordance, _ = cmd.Flags().GetFloat64("min-concordance")
		opts.PromoterLength, _ = cmd.Flags().GetInt("promoter-length")
		opts.ResLines = strings.Split(resLinesStr, ",")
		opts.SusLines = strings.Split(susLinesStr, ",")

		if gff == "" {
			return fmt.Errorf("--gff is required")
		}
		cmd.SilenceUsage = true
		_, _, err := genespace.GeneSpace(opts)
		return err
	},
}

func init() {
	rootCmd.AddCommand(GeneSpaceCmd)
	GeneSpaceCmd.Flags().SortFlags = false
	GeneSpaceCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file (superVCF table)")
	GeneSpaceCmd.Flags().String("chrom", "", "Chromosome name")
	GeneSpaceCmd.Flags().Int("start", 0, "Window start position")
	GeneSpaceCmd.Flags().IntP("stop", "e", 2000000000, "Window stop position")
	GeneSpaceCmd.Flags().String("res-lines", "", "Comma-separated resistant line names")
	GeneSpaceCmd.Flags().StringP("sus-lines", "u", "", "Comma-separated susceptible line names")
	GeneSpaceCmd.Flags().String("trait", genespace.DefaultProfile,
		"Trait profile: "+strings.Join(genespace.ProfileNames(), ", "))
	GeneSpaceCmd.Flags().Float64("min-concordance", 0.9,
		"Minimum genotype concordance (0-1] for a variant to count as following the R/S split")
	GeneSpaceCmd.Flags().Int("promoter-length", 2000,
		"Bases upstream of a gene whose variants count toward it (0 disables)")
	for _, f := range []string{"variant", "chrom", "res-lines", "sus-lines"} {
		_ = GeneSpaceCmd.MarkFlagRequired(f)
	}
}
