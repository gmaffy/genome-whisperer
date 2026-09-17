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

	// RunE, not Run: a returned error makes cobra exit non-zero. With Run the
	// command printed build failures and still exited 0, so a caller checking $?
	// carried on as though the database existed.
	RunE: func(cmd *cobra.Command, args []string) error {
		if err := utils.CheckDeps([]string{"gatk", "snpEff", "java"}); err != nil {
			return fmt.Errorf("dependency check failed: %w", err)
		}

		flags := cmd.Flags()
		refFile, err := flags.GetString("reference")
		if err != nil {
			return err
		}
		protein, err := flags.GetString("protein")
		if err != nil {
			return err
		}
		cds, err := flags.GetString("cds")
		if err != nil {
			return err
		}
		gff, err := flags.GetString("gff")
		if err != nil {
			return err
		}
		species, err := flags.GetString("species")
		if err != nil {
			return err
		}
		version, err := flags.GetString("annotation-version")
		if err != nil {
			return err
		}
		config, err := flags.GetString("config")
		if err != nil {
			return err
		}

		if config != "" {
			if _, sErr := os.Stat(config); sErr != nil {
				return fmt.Errorf("config file %s is not readable: %w", config, sErr)
			}
			return annotation.CreateCustomDbFromConfig(config, species, version)
		}

		fmt.Println("Creating custom db from command line arguments")
		for _, f := range []struct{ flag, path string }{
			{"--reference", refFile},
			{"--protein", protein},
			{"--cds", cds},
			{"--gff", gff},
		} {
			if f.path == "" {
				return fmt.Errorf("%s is required", f.flag)
			}
			if _, sErr := os.Stat(f.path); sErr != nil {
				return fmt.Errorf("%s file %s is not readable: %w", f.flag, f.path, sErr)
			}
		}
		if species == "" {
			return fmt.Errorf("please provide a species name with --species")
		}
		if version == "" {
			return fmt.Errorf("please provide the annotation version with --annotation-version")
		}

		fmt.Println("All arguments passed are valid")
		return annotation.CreateCustomDb(refFile, protein, cds, species, gff, version)
	},
}

func init() {
	rootCmd.AddCommand(CreateSnpEffDBCmd)
	CreateSnpEffDBCmd.Flags().SortFlags = false
	CreateSnpEffDBCmd.Flags().String("cds", "", "Path to the CDS FASTA")
	// --protein and --annotation-version are persistent flags on the root
	// command: GetProtPrgTopHits needs both too.
}
