package cmd

import (
	"fmt"

	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/spf13/cobra"
)

// GetProtPrgTopHitsCmd blasts a protein FASTA against the PRG database and
// writes the single best hit per protein beside the protein file.
var GetProtPrgTopHitsCmd = &cobra.Command{
	Use:   "GetProtPrgTopHits",
	Short: "BLASTs a protein FASTA against the PRG database and writes the top hit per protein",
	Long: `genome-whisperer GetProtPrgTopHits --protein <protein fasta> --prg-db <PRG blast db> --species <species> --annotation-version <version>

Runs blastp and reduces the result to one row per query protein — the hit with
the highest bitscore, ties going to the first one blastp reported. The output is
written next to the protein file as <SPECIES><ANNOTATION-VERSION>_PRG_TOP_HITS.txt,
which is the file --prg takes in AddPrg, CreateSuperVcf, VariantAnnotation and
GeneSpace.

The raw BLAST table is never written: blastp streams into the reducer, so a run
that would have produced a 127 MB table costs nothing but the answer.

--protein and --species may come from --config instead, from its "proteins:",
"Species:" and "Version:" keys. Flags win where both are given.

The PRG database is not built for you. If blastp cannot find one beside the path
given, build it with:

    makeblastdb -in <prgdb fasta> -dbtype prot -out <prgdb fasta>`,
	RunE: func(cmd *cobra.Command, args []string) error {
		protein, err := cmd.Flags().GetString("protein")
		if err != nil {
			return fmt.Errorf("reading protein: %w", err)
		}
		prgDb, err := cmd.Flags().GetString("prg-db")
		if err != nil {
			return fmt.Errorf("reading prg-db: %w", err)
		}
		species, err := cmd.Flags().GetString("species")
		if err != nil {
			return fmt.Errorf("reading species: %w", err)
		}
		annotationVersion, err := cmd.Flags().GetString("annotation-version")
		if err != nil {
			return fmt.Errorf("reading annotation-version: %w", err)
		}
		outDir, err := cmd.Flags().GetString("out-dir")
		if err != nil {
			return fmt.Errorf("reading out-dir: %w", err)
		}
		evalue, err := cmd.Flags().GetString("evalue")
		if err != nil {
			return fmt.Errorf("reading evalue: %w", err)
		}
		maxTargetSeqs, err := cmd.Flags().GetInt("max-target-seqs")
		if err != nil {
			return fmt.Errorf("reading max-target-seqs: %w", err)
		}
		threads, err := cmd.Flags().GetInt("threads")
		if err != nil {
			return fmt.Errorf("reading threads: %w", err)
		}
		verbose, err := cmd.Flags().GetBool("verbose")
		if err != nil {
			return fmt.Errorf("reading verbose: %w", err)
		}
		config, err := cmd.Flags().GetString("config")
		if err != nil {
			return fmt.Errorf("reading config: %w", err)
		}

		// A config file is a second source for the three identity values, not an
		// override: anything given on the command line stands.
		if config != "" {
			cfg, cErr := utils.ReadConfig(config)
			if cErr != nil {
				return fmt.Errorf("reading config %s: %w", config, cErr)
			}
			if protein == "" {
				protein = cfg.Proteins
			}
			if species == "" {
				species = cfg.Species
			}
			if annotationVersion == "" {
				annotationVersion = cfg.Version
			}
		}

		if protein == "" {
			return fmt.Errorf("provide a protein FASTA with --protein (or a \"proteins:\" line in --config)")
		}
		if prgDb == "" {
			return fmt.Errorf("provide the PRG protein BLAST database with --prg-db")
		}
		if species == "" {
			return fmt.Errorf("provide a species with --species (or a \"Species:\" line in --config); it names the output file")
		}

		out, err := annotation.BlastProtsToPrgDb(annotation.PrgBlastOptions{
			Proteins:          protein,
			PrgDb:             prgDb,
			Species:           species,
			AnnotationVersion: annotationVersion,
			OutDir:            outDir,
			Evalue:            evalue,
			Threads:           threads,
			MaxTargetSeqs:     maxTargetSeqs,
			Verbose:           verbose,
		})
		if err != nil {
			return err
		}

		fmt.Printf("PRG top hits written to %s\nPass it to the annotation commands with --prg %s\n", out, out)
		return nil
	},
}

func init() {
	rootCmd.AddCommand(GetProtPrgTopHitsCmd)
	GetProtPrgTopHitsCmd.Flags().SortFlags = false
	// No shorthands: -p and -e are already taken by GoBSAseq.
	GetProtPrgTopHitsCmd.Flags().String("prg-db", "", "blastp -db prefix for the PRG protein database")
	GetProtPrgTopHitsCmd.Flags().String("evalue", annotation.DefaultPrgEvalue, "blastp e-value cutoff")
	GetProtPrgTopHitsCmd.Flags().Int("max-target-seqs", 0, "Cap subjects per query (0 = uncapped; a cap is applied before the final ranking and can hide the best hit)")
	// --protein, --species, --annotation-version, --out-dir, --threads,
	// --verbose and --config are persistent flags on the root command.
}
