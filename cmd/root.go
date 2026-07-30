/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "genome-whisperer",
	Short: "Performs 1) Read Alignments 2) Variant Calling 3) Variant filtration 4) Variant Annotation 5) Pangenome assembly and 6) BSAseq",
	Long: `Tools Description:
    1. AlignReads
        - Aligns short paired reads to reference using bwa mem or bowtie2
        - Aligns long reads to reference using pbmm2
        - Marks duplicates using picard tools
        - Recalibrates bam using GATK's BQSR pipeline

    2. VariantCalling
        - Calls and hard filters SNPs and Indels using GATK best practices from bams generated from short reads
        - Calls and hard filters SNPs and Indels using Deepvariant from bams generated from long reads
        - Can use glenexus or GATK to merge gvcfs

    3. VariantAnnotation
        - Annotates VCFs using snpEff
        - Create snpeFF databases using species annotation files in database is not present

    4. Pangenome Assembly (GoPan)
        - Assembles a pangenome using the iterative mapping and assembly approach from short reads and a reference genome

    5. BSAseq (GoBSAseq)
        - Performs BSAseq Analysis.
        - Can take bam files, reads or VCF as input as well as a reference genome
        - Can use bulks only or bulks and parents as input.
        - Can use 1 or 2 bulks`,
	// Uncomment the following line if your bare application
	// has an action associated with it:
	// Run: func(cmd *cobra.Command, args []string) { },
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	fmt.Println("genome-whisperer ver 1.0.0")
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

func init() {
	// Flags used by more than one subcommand are defined once here as persistent
	// flags. Cobra makes them available to every subcommand and lists them under
	// "Global Flags:" in --help, so a shared flag has exactly one name, shorthand,
	// default and description across the whole CLI.
	//
	// Two flags that are shared by name cannot live here and stay in their own
	// commands instead:
	//   --variant  is a list in the annotation commands and a single file in the
	//              filter commands, and a persistent flag has one type.
	//   --bqsr     defaults to false for AlignReads/GoBSAseq but true for
	//              VariantCalling, and a persistent flag has one default.
	pf := rootCmd.PersistentFlags()
	pf.SortFlags = false

	// ---------------------------------- Reference and identity ------------------------------- //
	pf.StringP("reference", "r", "", "Path to the reference genome FASTA")
	pf.String("ref-version", "", "Reference genome version")
	pf.StringP("species", "s", "", "Species name")
	pf.StringP("genomes-dir", "g", "", "Directory containing prepared reference genomes")

	// ------------------------------------- Input selection ----------------------------------- //
	pf.StringP("config", "c", "", "Path to config file")
	pf.StringP("data-dir", "d", "", "Root of the standard data directory tree")
	pf.StringSliceP("bam", "b", []string{}, "Path to a BAM/CRAM file (repeatable)")

	// ----------------------------------------- Output ----------------------------------------- //
	pf.StringP("out-dir", "o", "", "Output directory")
	pf.StringP("log", "l", "", "Path to the run log file")
	pf.String("alignment-fmt", "bam", "Alignment output format: bam or cram")

	// ----------------------------------------- Runtime ---------------------------------------- //
	pf.IntP("threads", "t", 4, "Threads per job")
	pf.BoolP("verbose", "v", false, "Stream external tool output")
	pf.String("gatk-log-level", "INFO", "GATK log level: DEBUG, INFO, WARNING or ERROR")

	// --------------------------------------- Verification ------------------------------------- //
	pf.Bool("quick", false, "Use quick (header-only) integrity checks")
	pf.Bool("skip-verification", false, "Skip input integrity checks entirely")

	// ---------------------------------------- Alignment --------------------------------------- //
	pf.String("aligner", "bwa-mem2", "Aligner: bwa-mem, bwa-mem2, bowtie2 or pbmm2")
	pf.String("preset", "HIFI", "pbmm2 preset: SUBREAD, CSS, HIFI, ISOSEQ or UNROLLED")
	pf.Bool("bootstrap", false, "Bootstrap the BQSR known-sites set")
	pf.StringSliceP("known-sites", "k", []string{}, "Path to a known-sites VCF (repeatable)")

	// --------------------------------- Variant calling / merging ------------------------------ //
	pf.String("caller", "gatk", "Variant caller: gatk or deepvariant")
	pf.String("merger", "gatk", "gVCF merger: gatk or glnexus")
	pf.String("deepvariant-version", "1.10.0", "DeepVariant version")
	pf.String("model-type", "WGS", "DeepVariant model: WGS, WES, PACBIO, ONT_R104 or HYBRID_PACBIO_ILLUMINA")
	pf.Bool("no-merging", false, "Do not merge gVCFs")

	// ---------------------------------------- Annotation -------------------------------------- //
	pf.String("gff", "", "Path to the GFF3 annotation")
	pf.String("database", "", "snpEff database name")
	pf.Bool("bsaseq", false, "Emit BSAseq columns")
	pf.String("gene-description-tsv", "", "Gene description TSV (gene_id, gene_description)")
	pf.String("prg", "", "PRG BLAST results file")
}
