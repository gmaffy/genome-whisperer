package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

// CreateTemplateCmd represents the CreateTemplate command
var CreateTemplateCmd = &cobra.Command{
	Use:   "CreateTemplate",
	Short: "Write a config.txt template for the config-file input mode",
	Long: `CreateTemplate writes config.txt in the current directory. Fill it in and
pass it to other commands with --config.

The file lists the output directory, the reference and its annotation (gff,
proteins, cds), and the inputs: read pairs or single-end reads with sample and
library names, bams, BSA-seq bams with their role, VCFs and gVCFs. Repeat a line
once per file. For BSA-seq the library name is HIGH_PARENT, LOW_PARENT,
HIGH_BULK or LOW_BULK.

Flags given on the command line win over values in the file.

An existing config.txt in the current directory is overwritten.`,
	Example: `  # Write config.txt in the current directory, then fill it in
  genome-whisperer CreateTemplate`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("CreateTemplate called")
		templateFile, err := os.Create("config.txt")
		if err != nil {
			fmt.Println("Error creating template file:", err)
			return
		}
		defer templateFile.Close()
		templateContent := `#===================================================================================================================== #
OutputDir: <path to output directory>
# =============================================== Reference ========================================================== #
Reference: <path to reference fasta file>
gff: <path to gff file>
proteins: <path to protein fasta file>
cds: <path to cds fasta file>
# =============================================== Read alignment ===================================================== #

ReadPair: <path to forward reads> <path to reverse reads> <sample name> <library name (for BSAseq this is HIGH_PARENT, LOW_PARENT, HIGH_BULK OR LOW_BULK)>
SingeEndReads: <path to reads> <sample name> <library name>

# =============================================== Bam files ========================================================== #
bam: <path to bam/cram file>
BSAseqBam: <path to bam/cram file> <HIGH_PARENT, LOW_PARENT, HIGH_BULK, LOW_BULK>

# =============================================== VCF files ========================================================== #
VCF: <path to vcf file>
gvcf: <path to gvcf file>

# ==================================================================================================================== #`
		_, err = templateFile.WriteString(templateContent)
		if err != nil {
			fmt.Println("Error writing to template file:", err)
			return
		}
		fmt.Println("Template file created successfully.")

	},
}

func init() {
	rootCmd.AddCommand(CreateTemplateCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// CreateTemplateCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// CreateTemplateCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
