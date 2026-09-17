/*
Copyright © 2026 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"fmt"
	"os"

	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/spf13/cobra"
)

// AddPrgCmd represents the AddPrg command
var AddPrgCmd = &cobra.Command{
	Use:   "AddPrg",
	Short: "Adds protein hits stats  to the PRG",
	Long:  `Takes blast output and adds protein hits stats to the PRG`,
	RunE: func(cmd *cobra.Command, args []string) error {
		fmt.Println("Adding PRG blast hits to vcf files ................")
		vcfs, vErr := cmd.Flags().GetStringSlice("variant")
		if vErr != nil {
			return fmt.Errorf("reading variant: %w", vErr)
		}
		// The second parameter of annotation.AddPrg is the PRG blast file, so
		// this is --prg. It used to be handed --gene-description-tsv, which made
		// the command unable to do the one thing it is named for; gene
		// descriptions are AddGeneDescriptions' job.
		prgFile, pErr := cmd.Flags().GetString("prg")
		if pErr != nil {
			return fmt.Errorf("reading prg: %w", pErr)
		}
		bsaseq, bErr := cmd.Flags().GetBool("bsaseq")
		if bErr != nil {
			return fmt.Errorf("reading bsaseq: %w", bErr)
		}
		if prgFile == "" {
			return fmt.Errorf("provide a PRG blast file with --prg. Create one with: genome-whisperer GetProtPrgTopHits")
		}
		if _, err := os.Stat(prgFile); err != nil {
			return fmt.Errorf("PRG blast file %s is not readable: %w", prgFile, err)
		}
		if len(vcfs) == 0 {
			return fmt.Errorf("provide at least one vcf file with --variant")
		}
		for i := range vcfs {
			if _, err := os.Stat(vcfs[i]); err != nil {
				return fmt.Errorf("vcf file %s is not a valid file path: %w", vcfs[i], err)
			}
		}

		if err, _ := annotation.AddPrg(vcfs, prgFile, bsaseq); err != nil {
			return err
		}
		return nil
	},
}

func init() {
	rootCmd.AddCommand(AddPrgCmd)
	AddPrgCmd.Flags().SortFlags = false
	AddPrgCmd.Flags().StringSliceP("variant", "V", []string{}, "Path to a VCF/variant file (repeatable)")
	if err := AddPrgCmd.MarkFlagRequired("variant"); err != nil {
		return
	}
}
