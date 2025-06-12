/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"
	"os"

	"github.com/spf13/cobra"
)

// createSnpEffDbCmd represents the createSnpEffDb command
var createSnpEffDbCmd = &cobra.Command{
	Use:   "createSnpEffDb",
	Short: "Creates a snpEff database from a reference genome, protein fasta, cds fasta and gff3 file.",
	Long:  `Creates a snpEff database from a reference genome, protein fasta, cds fasta and gff3 file.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("createSnpEffDb called")
		err := utils.CheckDeps()
		if err != nil {
			return
		}
		ref, rErr := cmd.Flags().GetString("reference")
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
			err1 := annotation.CreateCustomDb(ref, protein, cds, species, gff, version)
			if err1 != nil {
				return
			}

		}

	},
}

func init() {
	rootCmd.AddCommand(createSnpEffDbCmd)

	// Add all flags
	createSnpEffDbCmd.Flags().String("protein", "", "Path to protein fasta file ...")
	createSnpEffDbCmd.Flags().String("cds", "", "Path to cds fasta file ...")
	createSnpEffDbCmd.Flags().String("gff", "", "Path to gff3 file ...")
	createSnpEffDbCmd.Flags().String("species", "", "Species name (no spaces or special characters) ...")
	createSnpEffDbCmd.Flags().String("version", "", "Reference annotation version ...")

	// Check if -c flag is provided via persistent flags
	cFlag, _ := rootCmd.PersistentFlags().GetString("config")
	if cFlag != "" {
		// Mark flags as required only if -c is false
		requiredFlags := []string{"protein", "cds", "species", "gff", "reference", "version"}
		for _, flag := range requiredFlags {
			err := createSnpEffDbCmd.MarkFlagRequired(flag)
			if err != nil {
				return
			}
		}
	} else {
		requiredFlags := []string{"version", "species"}
		for _, flag := range requiredFlags {
			err := createSnpEffDbCmd.MarkFlagRequired(flag)
			if err != nil {
				return
			}
		}

	}
}
