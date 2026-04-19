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

// AddPrgCmd represents the AddPrg command
var AddPrgCmd = &cobra.Command{
	Use:   "AddPrg",
	Short: "Adds protein hits stats  to the PRG",
	Long:  `Takes blast output and adds protein hits stats to the PRG`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("AddPrg called")
		fmt.Println("Adding gene descriptions to vcf files ................")
		vcfs, vErr := cmd.Flags().GetStringSlice("variant")
		if vErr != nil {
			fmt.Println("Error getting variant flag")
		}
		descFile, dErr := cmd.Flags().GetString("gene-description-tsv")
		if dErr != nil {
			fmt.Println("Error getting gene-description flag")
		}
		bsaseq, bErr := cmd.Flags().GetBool("bsaseq")
		if bErr != nil {
			fmt.Println("Error getting bsaseq flag")
		}
		//fmt.Println(bsaseq)
		//fmt.Println(vcfs)
		//fmt.Println(descFile)
		if descFile == "" {
			fmt.Println("Please provide gene-description-tsv name with flag -g ")
			return
		}
		if len(vcfs) == 0 {
			fmt.Println("Please provide at least one vcf file")
			return
		}
		for i := range vcfs {
			_, err := os.Stat(vcfs[i])
			if err != nil {
				fmt.Printf("Vcf file: %s is not a valid file path", vcfs[i])
				log.Fatal(err)
			}
		}

		err, _ := annotation.AddPrg(vcfs, descFile, bsaseq)
		if err != nil {
			fmt.Println(err)
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(AddPrgCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// AddPrgCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// AddPrgCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	AddPrgCmd.Flags().StringSliceP("variant", "V", []string{}, "snpEff tsv file ...")
	AddPrgCmd.Flags().StringP("gene-description-tsv", "g", "", "gene description tsv file ...")
	AddPrgCmd.Flags().Bool("bsaseq", false, "bsaseq output")
}
