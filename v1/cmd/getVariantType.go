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

// getVariantTypeCmd represents the getVariantType command
var getVariantTypeCmd = &cobra.Command{
	Use:   "getVariantType",
	Short: "Filters a vcf file to only include variants of a given type, e.g. SNPs, INDELs, MNPs, etc.",
	Long:  `Filters a vcf file to only include variants of a given type, e.g. SNPs, INDELs, MNPs, etc.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("getVariantType called")
		variant, vErr := cmd.Flags().GetString("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		variantType, tErr := cmd.Flags().GetString("type")
		if tErr != nil {
			fmt.Println("Error getting variant type flag")
		}

		_, err := os.Stat(variant)
		if err != nil {
			fmt.Printf("Variant file %s is not a valid file: %v\n", variant, err)
			return
		}
		gErr := variants.GetVariantType(variant, variantType)
		if gErr != nil {
			fmt.Printf("Error getting variant type: %v\n", gErr)
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(getVariantTypeCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// getVariantTypeCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	getVariantTypeCmd.Flags().StringP("type", "t", "", "SNP, INDEL, MNP, etc.")
	getVariantTypeCmd.Flags().StringP("variant", "V", "", "VCF file")
}
