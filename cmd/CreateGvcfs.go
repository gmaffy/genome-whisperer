/*
Copyright © 2026 Godwin Mafireyi (mafireyi@gmail.com)
*/
package cmd

import (
	"fmt"
	"log"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"github.com/spf13/cobra"
)

// variantOptions reads the flags every pipeline command shares into a
// variants.Options. CreateGvcfs, MergeGvcfs and CreateSuperVcf all use it so the
// three stages cannot disagree about how a flag is interpreted.
//
// Filter thresholds are not read here: they are local to VariantCalling, which is
// the only command that exposes them. Everything else filters with the defaults.
func variantOptions(cmd *cobra.Command) (variants.Options, error) {
	get := cmd.Flags()

	dataDir, err := get.GetString("data-dir")
	if err != nil {
		return variants.Options{}, err
	}
	bams, err := get.GetStringSlice("bam")
	if err != nil {
		return variants.Options{}, err
	}
	outDir, err := get.GetString("out-dir")
	if err != nil {
		return variants.Options{}, err
	}
	species, err := get.GetString("species")
	if err != nil {
		return variants.Options{}, err
	}
	refVer, err := get.GetString("ref-version")
	if err != nil {
		return variants.Options{}, err
	}
	refFasta, err := get.GetString("reference")
	if err != nil {
		return variants.Options{}, err
	}
	genomesDir, err := get.GetString("genomes-dir")
	if err != nil {
		return variants.Options{}, err
	}
	caller, err := get.GetString("caller")
	if err != nil {
		return variants.Options{}, err
	}
	merger, err := get.GetString("merger")
	if err != nil {
		return variants.Options{}, err
	}
	dvVer, err := get.GetString("deepvariant-version")
	if err != nil {
		return variants.Options{}, err
	}
	modelType, err := get.GetString("model-type")
	if err != nil {
		return variants.Options{}, err
	}
	gatkLogLevel, err := get.GetString("gatk-log-level")
	if err != nil {
		return variants.Options{}, err
	}
	threads, err := get.GetInt("threads")
	if err != nil {
		return variants.Options{}, err
	}
	verbose, err := get.GetBool("verbose")
	if err != nil {
		return variants.Options{}, err
	}
	quick, err := get.GetBool("quick")
	if err != nil {
		return variants.Options{}, err
	}
	skipVerification, err := get.GetBool("skip-verification")
	if err != nil {
		return variants.Options{}, err
	}
	noMerging, err := get.GetBool("no-merging")
	if err != nil {
		return variants.Options{}, err
	}

	// --model-type only has a meaning for DeepVariant. Leaving it empty lets
	// CreateGvcfs fall back to the data directory's long-read naming convention,
	// so only pass it through when the user actually set it.
	if !get.Changed("model-type") {
		modelType = ""
	}

	// Resolve the reference from --genomes-dir when one is given. Otherwise leave
	// it as passed: it may still be empty here and be filled in from a config
	// file by the caller, and the stages validate it before use.
	resolvedFasta := refFasta
	if genomesDir != "" {
		resolved, rErr := utils.ResolveRefFasta(refFasta, genomesDir, species, refVer)
		if rErr != nil {
			return variants.Options{}, fmt.Errorf("could not resolve reference fasta: %w", rErr)
		}
		resolvedFasta = resolved
	}

	return variants.Options{
		DataDir:          dataDir,
		Bams:             bams,
		OutDir:           outDir,
		Species:          species,
		RefVer:           refVer,
		RefFasta:         resolvedFasta,
		Caller:           caller,
		Merger:           merger,
		DVVer:            dvVer,
		ModelType:        modelType,
		GatkLogLevel:     gatkLogLevel,
		Threads:          threads,
		Verbose:          verbose,
		Quick:            quick,
		SkipVerification: skipVerification,
		NoMerging:        noMerging,
	}, nil
}

// CreateGvcfsCmd represents the CreateGvcfs command
var CreateGvcfsCmd = &cobra.Command{
	Use:   "CreateGvcfs",
	Short: "Create per-sample, per-chromosome gVCFs from bam/cram files",
	Long: `CreateGvcfs
        - Calls variants per chromosome with GATK HaplotypeCaller or DeepVariant
        - Takes bams from the standard data directory (--data-dir), a config file
          (--config) or explicit paths (--bam)
        - Reuses gVCFs that already exist and are valid, and re-creates corrupt ones`,
	RunE: func(cmd *cobra.Command, args []string) error {
		opts, err := variantOptions(cmd)
		if err != nil {
			return err
		}

		// Consumes BQSR-recalibrated alignments by default, like VariantCalling.
		bqsr, err := cmd.Flags().GetBool("bqsr")
		if err != nil {
			return err
		}
		opts.NoBqsr = !bqsr

		// ------------------------------- Check dependencies ------------------------------- //
		switch opts.Caller {
		case "gatk":
			if depErr := utils.CheckDeps([]string{"gatk", "bcftools"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
		case "deepvariant":
			if depErr := utils.CheckDeps([]string{"docker", "bcftools"}); depErr != nil {
				return fmt.Errorf("dependency check failed: %w", depErr)
			}
		default:
			return fmt.Errorf("caller must be gatk or deepvariant, got %q", opts.Caller)
		}

		gvcfs, err := variants.CreateGvcfs(opts)
		if err != nil {
			return err
		}

		for chrom, paths := range gvcfs {
			log.Printf("%s: %d gVCF(s)", chrom, len(paths))
		}
		return nil
	},
}

func init() {
	rootCmd.AddCommand(CreateGvcfsCmd)
	CreateGvcfsCmd.Flags().SortFlags = false
	// Same default as VariantCalling: consume BQSR-recalibrated alignments, and
	// --bqsr=false falls back to the RGMD bam/cram.
	CreateGvcfsCmd.Flags().Bool("bqsr", true, "Run base quality score recalibration")
}
