/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

// CreateTemplateCmd represents the CreateTemplate command
var CreateTemplateCmd = &cobra.Command{
	Use:   "CreateTemplate",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
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
SingeEndReadsToAlign: <path to reads> <sample name> <library name>

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
