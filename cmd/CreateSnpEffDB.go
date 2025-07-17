/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"
	"os"

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

		version, vErr := cmd.Flags().GetString("version")
		if vErr != nil {
			fmt.Println("Error getting version flag")
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
				fmt.Println("Please provide version")
				return
			}
			err1 := annotation.CreateCustomDb(refFile, protein, cds, species, gff, version)
			if err1 != nil {
				return
			}

		}
	},
}

func init() {
	rootCmd.AddCommand(CreateSnpEffDBCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// CreateSnpEffDBCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// CreateSnpEffDBCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
	CreateSnpEffDBCmd.Flags().String("protein", "", "Path to protein fasta file ...")
	CreateSnpEffDBCmd.Flags().String("cds", "", "Path to cds fasta file ...")
	CreateSnpEffDBCmd.Flags().String("gff", "", "Path to gff3 file ...")
	CreateSnpEffDBCmd.Flags().StringP("species", "s", "", "Species name (no spaces or special characters) ...")
	CreateSnpEffDBCmd.Flags().StringP("version", "v", "", "Reference annotation version ...")
	CreateSnpEffDBCmd.Flags().StringP("reference", "r", "", "Path to reference genome fasta file ...")
	CreateSnpEffDBCmd.Flags().StringP("config", "c", "", "Path to config file ...")

	// Check if -c flag is provided via persistent flags
	cFlag, _ := CreateSnpEffDBCmd.Flags().GetString("config")
	if cFlag == "" {
		// Mark flags as required only if -c is false
		requiredFlags := []string{"protein", "cds", "species", "gff", "reference", "version"}
		for _, flag := range requiredFlags {
			err := CreateSnpEffDBCmd.MarkFlagRequired(flag)
			if err != nil {
				return
			}
		}
	} else {
		requiredFlags := []string{"version", "species"}
		for _, flag := range requiredFlags {
			err := CreateSnpEffDBCmd.MarkFlagRequired(flag)
			if err != nil {
				return
			}
		}

	}
}
