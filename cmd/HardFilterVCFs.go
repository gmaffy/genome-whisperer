/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

// HardFilterVCFsCmd represents the HardFilterVCFs command
var HardFilterVCFsCmd = &cobra.Command{
	Use:   "HardFilterVCFs",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("HardFilterVCFs called")
		variant, vErr := cmd.Flags().GetString("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		gatkLogLevel, gatErr := cmd.Flags().GetString("gatk-log-level")
		if gatErr != nil {
			fmt.Println("Error getting gatk-log-level flag")
		}

		verbose, vErr := cmd.Flags().GetBool("verbose")
		if vErr != nil {
			fmt.Println("Error getting verbose flag")
		}

		_, err := os.Stat(variant)
		if err != nil {
			fmt.Printf("Variant file %s is not a valid file: %v\n", variant, err)
			return
		}
		snps, vaErr := variants.HardFilterSNPs(variant, gatkLogLevel, verbose)
		if vaErr != nil {
			return
		}

		indels, inErr := variants.HardFilterINDELs(variant, gatkLogLevel, verbose)
		if inErr != nil {
			return
		}

		mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s --verbosity %s`, snps, indels, strings.Replace(snps, ".SNP.vcf.gz", ".hard_filtered.vcf.gz", 1), gatkLogLevel)

		var mErr error
		if verbose {
			mErr = utils.RunBashCmdVerbose(mergeCmdStr)
		} else {
			mErr = utils.RunBashCmd(mergeCmdStr)
		}

		if mErr != nil {
			fmt.Printf("Error running merge command: %v\n", mErr)
			return
		}

	},
}

func init() {
	rootCmd.AddCommand(HardFilterVCFsCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// HardFilterVCFsCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// HardFilterVCFsCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	HardFilterVCFsCmd.Flags().StringP("variant", "V", "", "INDEL VCF file")
	HardFilterVCFsCmd.Flags().String("gatk-log-level", "INFO", "GATK log level")
	HardFilterVCFsCmd.Flags().BoolP("verbose", "v", false, "Verbose output")
	err := HardFilterVCFsCmd.MarkFlagRequired("variant")
	if err != nil {
		return
	}
}
