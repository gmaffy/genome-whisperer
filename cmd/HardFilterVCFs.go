/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"

	"github.com/spf13/cobra"
)

// HardFilterVCFsCmd represents the HardFilterVCFs command
var HardFilterVCFsCmd = &cobra.Command{
	Use:   "HardFilterVCFs",
	Short: "Hard filter all variants (snps and indels)",
	Long:  `Filters both SNPs and INDELs using GATK recommended parameters and merges them to hard filtered VCF ... `,
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
	HardFilterVCFsCmd.Flags().SortFlags = false
	HardFilterVCFsCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")
	if err := HardFilterVCFsCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
}
