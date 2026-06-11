/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/annotation"
	"log"
	"os"

	"github.com/spf13/cobra"
)

// variantAnnotationCmd represents the variantAnnotation command
var variantAnnotationCmd = &cobra.Command{
	Use:   "variantAnnotation",
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
		fmt.Println(bsaseq)
		fmt.Println(vcfs)
		fmt.Println(db)
		if db == "" {
			fmt.Println("Please provide database name with flag --database ")
			return
		}
		if len(vcfs) == 0 {
			fmt.Println("Please provide at least one vcf file")
			return
		} else {
			for i, _ := range vcfs {
				_, err := os.Stat(vcfs[i])
				if err != nil {
					fmt.Printf("Vcf file: %s is not a valid file path", vcfs[i])
					log.Fatal(err)
				}
			}
		}

		err := annotation.RunSnpEff(vcfs, db, bsaseq)
		if err != nil {
			fmt.Println(err)
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(variantAnnotationCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// variantAnnotationCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// variantAnnotationCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	variantAnnotationCmd.Flags().StringSliceP("variant", "V", []string{}, "Variant file ...")
	variantAnnotationCmd.Flags().StringP("database", "d", "", "Species name")
	variantAnnotationCmd.Flags().Bool("bsaseq", false, "output bsaseq columns")
}
