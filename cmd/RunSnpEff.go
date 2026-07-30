/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/spf13/cobra"
)

// RunSnpEffCmd represents the RunSnpEff command
var RunSnpEffCmd = &cobra.Command{
	Use:   "RunSnpEff",
	Short: "annotation of variants using snpEff",
	Long:  `Annotation of variants using snpEff. Outputs a vcf file and tsv file with annotations.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("RunSnpEff called")
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
			for i := range vcfs {
				_, err := os.Stat(vcfs[i])
				if err != nil {
					fmt.Printf("Vcf file: %s is not a valid file path", vcfs[i])
					log.Fatal(err)
				}
			}
		}

		err, snpEffTsvFiles, snpEffVcfFiles := annotation.RunSnpEff(vcfs, db, bsaseq)
		if err != nil {
			fmt.Println(err)
			return
		}
		fmt.Printf("snpEff vcfs created: %s\n---------------------\nTsv files created: %s\n------------------\n", snpEffVcfFiles, snpEffTsvFiles)
	},
}

func init() {
	rootCmd.AddCommand(RunSnpEffCmd)
	RunSnpEffCmd.Flags().SortFlags = false
	RunSnpEffCmd.Flags().StringSliceP("variant", "V", []string{}, "Path to a VCF/variant file (repeatable)")
	if err := RunSnpEffCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
}
