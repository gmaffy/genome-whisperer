/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"os"

	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"

	"github.com/spf13/cobra"
)

// CreateSnpEffDBCmd represents the CreateSnpEffDB command
var CreateSnpEffDBCmd = &cobra.Command{
	Use:   "CreateSnpEffDB",
	Short: "Creates a snpEff database from a reference genome, protein fasta, cds fasta and gff3 file.",
	Long:  `Creates a snpEff database from a reference genome, protein fasta, cds fasta and gff3 file.`,

	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("CreateSnpEffDB called")
		err := utils.CheckDeps([]string{"gatk", "snpEff", "java"})
		if err != nil {
			return
		}
		refFile, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			fmt.Println("Error getting reference flag")
			return
		}

		protein, pErr := cmd.Flags().GetString("protein")
		if pErr != nil {
			fmt.Println("Error getting protein flag")
			return
		}

		cds, cErr := cmd.Flags().GetString("cds")
		if cErr != nil {
			fmt.Println("Error getting cds flag")
			return
		}

		gff, gErr := cmd.Flags().GetString("gff")
		if gErr != nil {
			fmt.Println("Error getting gff flag")
			return
		}

		species, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			fmt.Println("Error getting species flag")
			return

		}

		version, vErr := cmd.Flags().GetString("annotation-version")
		if vErr != nil {
			fmt.Println("Error getting annotation-version flag")
			return
		}

		config, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			fmt.Println("Error getting config flag")
			return
		}

		if config != "" {
			_, err := os.Stat(config)
			if err != nil {
				fmt.Printf("config file: %s is not a valid file path", config)
				return
			}
			er := annotation.CreateCustomDbFromConfig(config, species, version)
			if er != nil {
				fmt.Printf("Error creating custom db from config file: %v", er)
				return
			}

		} else {
			fmt.Println("Creating custom db from command line arguments")
			_, err := os.Stat(refFile)
			if err != nil {
				fmt.Printf("Reference file: %s is not a valid file path", refFile)
				return
			}

			_, err = os.Stat(protein)
			if err != nil {
				fmt.Printf("Protein file: %s is not a valid file path", protein)
				return
			}

			_, err = os.Stat(cds)
			if err != nil {
				fmt.Printf("CDS file: %s is not a valid file path", cds)
				return
			}

			_, err = os.Stat(gff)
			if err != nil {
				fmt.Printf("GFF file: %s is not a valid file path", gff)
				return
			}

			if species == "" {
				fmt.Println("Please provide species name")
				return
			}

			if version == "" {
				fmt.Println("Please provide the annotation version with --annotation-version")
				return
			}

			fmt.Println("All arguments passed are valid")

			err1 := annotation.CreateCustomDb(refFile, protein, cds, species, gff, version)
			if err1 != nil {
				fmt.Printf("Error creating custom db: %v", err1)
				return
			}

		}
	},
}

func init() {
	rootCmd.AddCommand(CreateSnpEffDBCmd)
	CreateSnpEffDBCmd.Flags().SortFlags = false
	CreateSnpEffDBCmd.Flags().String("protein", "", "Path to the protein FASTA")
	CreateSnpEffDBCmd.Flags().String("cds", "", "Path to the CDS FASTA")
	CreateSnpEffDBCmd.Flags().String("annotation-version", "", "Reference annotation version")
}
