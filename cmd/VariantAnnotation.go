/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/gmaffy/genome-whisperer/annotation"

	"github.com/spf13/cobra"
)

// VariantAnnotationCmd represents the VariantAnnotation command
var VariantAnnotationCmd = &cobra.Command{
	Use:   "VariantAnnotation",
	Short: "Annotation of variants using snpEff",
	Long:  `Annotation of variants using snpEff`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("variantAnnotation called")
		vcfs, vErr := cmd.Flags().GetStringSlice("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}

		db, dErr := cmd.Flags().GetString("database")
		if dErr != nil {
			fmt.Println("Error getting database flag")
		}

		bsaseq, bErr := cmd.Flags().GetBool("bsaseq")
		if bErr != nil {
			fmt.Println("Error getting bsaseq flag")
		}

		desc, dErr := cmd.Flags().GetString("description")
		if dErr != nil {
			fmt.Println("Error getting description flag")
		}

		prg, pErr := cmd.Flags().GetString("prg")
		if pErr != nil {
			fmt.Println("Error getting prg flag")
		}

		if db == "" {
			fmt.Println("Please provide database name with flag --database ")
			return
		}

		if len(vcfs) == 0 {
			fmt.Println("Please provide at least one vcf file")
			return
		} else {
			for i := range vcfs {
				_, err := os.Stat(vcfs[i])
				if err != nil {
					fmt.Printf("Vcf file: %s is not a valid file path", vcfs[i])
					log.Fatal(err)
				}
			}
		}

		err, snpEffTsvFiles := annotation.CreateSuperVcf(vcfs, db, bsaseq, desc, prg)
		if err != nil {
			fmt.Println(err)
			return
		}
		fmt.Printf("SUPER VCFS created: %s\n---------------------\n\n", snpEffTsvFiles)
	},
}

func init() {
	rootCmd.AddCommand(VariantAnnotationCmd)

	VariantAnnotationCmd.Flags().StringSliceP("variant", "V", []string{}, "Variant file ...")
	VariantAnnotationCmd.Flags().StringP("database", "d", "", "database name")
	VariantAnnotationCmd.Flags().Bool("bsaseq", false, "output bsaseq columns")
	VariantAnnotationCmd.Flags().StringP("gene-description-file", "g", "", "Gene description file")
	VariantAnnotationCmd.Flags().StringP("prg", "p", "", "PRG blast file")
}
