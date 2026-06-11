/*
Copyright © 2025 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/pangenome"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"

	"github.com/spf13/cobra"
)

// goPanCmd represents the goPan command
var goPanCmd = &cobra.Command{
	Use:   "goPan -c <config file> -a <masurca, megahit or mac>",
	Short: "Creates a pangenome using the iterative mapping approach",
	Long:  `Creates a pangenome using the iterative mapping approach. Inputs are short reads and reference genome.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("goPan called")
		fmt.Printf("Checking dependencies ...\n\n")

		if err := utils.CheckDeps(); err != nil {
			log.Fatalf("Dependency check failed: %v", err)
		}

		fmt.Printf("Dependencies OK\n\n----------------------------------------------------------\n\n")
		configFile, cErr := cmd.Flags().GetString("config")
		if cErr != nil {
			log.Fatalf("Error getting config flag: %v", cErr)
		}

		assembler, aErr := cmd.Flags().GetString("assembler")
		if aErr != nil {
			log.Fatalf("Error getting assembler flag: %v", aErr)
		}

		fmt.Printf("Running with the following parameters:\nConfig file: %s\nAssembler: %s\n ...\n\n", configFile, assembler)
		pangenome.GoPan(configFile, assembler)
	},
}

func init() {
	rootCmd.AddCommand(goPanCmd)
	goPanCmd.Flags().StringP("assembler", "a", "mac", "masurca, megahit or mac")
	goPanCmd.Flags().StringP("config", "c", "", "Config file")
	err := goPanCmd.MarkFlagRequired("config")
	if err != nil {
		return
	}

}
