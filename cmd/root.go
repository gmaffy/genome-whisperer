/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
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
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

func init() {
	// Here you will define your flags and configuration settings.
	// Cobra supports persistent flags, which, if defined here,
	// will be global for your application.

	// rootCmd.PersistentFlags().StringVar(&cfgFile, "config", "", "config file (default is $HOME/.genome-whisperer.yaml)")

	// Cobra also supports local flags, which will only run
	// when this action is called directly.
	// rootCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
