/*
Copyright © 2025 NAME HERE <EMAIL ADDRESS>
*/
package cmd

import (
	"os"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "genome-whisperer",
	Short: "A toolkit for genome analysis",
	Long: `A multi-purpose genome analysis tool for performing:
1.	Read alignment: ( bwa, bowtie2 and pbmm2)
2.	Variant calling & Variant Filtration : (GATK best practices)
3.	Variant annotation: ( SnpEff)
4.	BSAseq
5.	Pangenome analysis 
6.	Other utils
`,
	// Uncomment the following line if your bare application
	// has an action associated with it:
	// Run: func(cmd *cobra.Command, args []string) { },
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

var cfgFile string
var refFile string

func init() {
	rootCmd.PersistentFlags().StringVarP(&cfgFile, "config", "c", "", "path to config file ")
	rootCmd.PersistentFlags().StringVarP(&refFile, "reference", "r", "", "path to reference genome fasta file ")
	// Here you will define your flags and configuration settings.
	// Cobra supports persistent flags, which, if defined here,
	// will be global for your application.

	// rootCmd.PersistentFlags().StringVar(&cfgFile, "config", "", "config file (default is $HOME/.genome-whisperer-go.yaml)")

	// Cobra also supports local flags, which will only run
	// when this action is called directly.
	//rootCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
