/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"

	"github.com/spf13/cobra"
)

// createSnpEffDbCmd represents the createSnpEffDb command
var createSnpEffDbCmd = &cobra.Command{
	Use:   "createSnpEffDb",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("createSnpEffDb called")
		err := utils.CheckDeps()
		if err != nil {
			return
		}
		ref, rErr := cmd.Flags().GetString("reference")
		if rErr != nil {
			fmt.Println("Error getting reference flag")
		}

		protein, pErr := cmd.Flags().GetString("protein")
		if pErr != nil {
			fmt.Println("Error getting protein flag")
		}

		cds, cErr := cmd.Flags().GetString("cds")
		if cErr != nil {
			fmt.Println("Error getting cds flag")
		}

		gff, gErr := cmd.Flags().GetString("gff")
		if gErr != nil {
			fmt.Println("Error getting reference flag")
		}

		species, sErr := cmd.Flags().GetString("species")
		if sErr != nil {
			fmt.Println("Error getting reference flag")
		}

		version, vErr := cmd.Flags().GetString("version")
		if vErr != nil {
			fmt.Println("Error getting reference flag")
		}

		err1 := annotation.CreateCustomDb(ref, protein, cds, species, gff, version)
		if err1 != nil {
			return
		}
	},
}

func init() {
	rootCmd.AddCommand(createSnpEffDbCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// createSnpEffDbCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	createSnpEffDbCmd.Flags().String("protein", "", "Path to protein fasta file ...")
	createSnpEffDbCmd.Flags().String("cds", "", "Path to cds fasta file ...")
	createSnpEffDbCmd.Flags().String("gff", "", "Path to gff3 file ...")
	createSnpEffDbCmd.Flags().String("species", "", "Species name (no spaces or special characters) ...")
	createSnpEffDbCmd.Flags().String("version", "", "Reference annotation version ...")
	err := createSnpEffDbCmd.MarkFlagRequired("protein")
	if err != nil {
		return
	}
	err = createSnpEffDbCmd.MarkFlagRequired("cds")
	if err != nil {
		return
	}

	err = createSnpEffDbCmd.MarkFlagRequired("cds")
	if err != nil {
		return
	}

	err = createSnpEffDbCmd.MarkFlagRequired("species")
	if err != nil {

	}

	err = createSnpEffDbCmd.MarkFlagRequired("gff")
	if err != nil {
		return
	}

	err = createSnpEffDbCmd.MarkFlagRequired("reference")
	if err != nil {
		return
	}
}
