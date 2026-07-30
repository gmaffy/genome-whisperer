/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/variants"
	"os"

	"github.com/spf13/cobra"
)

// GetVariantTypeCmd represents the GetVariantType command
var GetVariantTypeCmd = &cobra.Command{
	Use:   "GetVariantType",
	Short: "Filters a vcf file to only include variants of a given type, e.g. SNPs, INDELs, MNPs, etc.",
	Long:  `Filters a vcf file to only include variants of a given type, e.g. SNPs, INDELs, MNPs, etc.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("GetVariantType called")

		variant, vErr := cmd.Flags().GetString("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		variantType, tErr := cmd.Flags().GetString("variant-type")
		if tErr != nil {
			fmt.Println("Error getting variant type flag")
		}

		gatkLogLevel, gatErr := cmd.Flags().GetString("gatk-log-level")
		if gatErr != nil {
			fmt.Println("Error getting gatk log level flag")
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
		_, gErr := variants.GetVariantType(variant, variantType, gatkLogLevel, verbose)
		if gErr != nil {
			fmt.Printf("Error getting variant type: %v\n", gErr)
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(GetVariantTypeCmd)
	GetVariantTypeCmd.Flags().SortFlags = false
	GetVariantTypeCmd.Flags().StringP("variant", "V", "", "Path to a VCF/variant file")
	GetVariantTypeCmd.Flags().String("variant-type", "", "Variant type to select: SNP, INDEL, MNP, etc.")
	if err := GetVariantTypeCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
	if err := GetVariantTypeCmd.MarkFlagRequired("variant-type"); err != nil {
		return
	}
}
