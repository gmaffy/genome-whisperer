package variants

import (
	"fmt"
	"strings"

	"github.com/fatih/color"
)

// RunPipeline runs the whole variant calling pipeline: gVCF creation, joint
// genotyping and concatenation, then hard filtering. It returns the path of the
// final VCF.
//
// All three input modes reach this function through the same Options struct, so
// data-dir, config and inline runs execute exactly the same code. That is the
// point of the rewrite: previously the data-dir path stopped after gVCF creation
// and never merged at all.
//
// Named RunPipeline rather than runVariantCaller.go's VariantCalling so both can
// coexist while bsaseq and the original engine are migrated.
func RunPipeline(opts Options) (string, error) {

	// ================================= Stage 1: create gVCFs ================================== //

	color.Cyan("\n%s STAGE 1/3  CREATE GVCFS %s\n\n", strings.Repeat("=", 24), strings.Repeat("=", 24))

	gvcfs, err := CreateGvcfs(opts)
	if err != nil {
		return "", fmt.Errorf("creating gVCFs: %w", err)
	}

	if opts.NoMerging {
		color.Yellow("--no-merging set: stopping after gVCF creation\n")
		return "", nil
	}

	// ================================== Stage 2: merge gVCFs =================================== //

	color.Cyan("\n%s STAGE 2/3  MERGE GVCFS %s\n\n", strings.Repeat("=", 24), strings.Repeat("=", 24))

	jointVCF, err := MergeGvcfs(opts, gvcfs)
	if err != nil {
		return "", fmt.Errorf("merging gVCFs: %w", err)
	}

	if opts.NoHardFilter {
		color.Yellow("--no-hard-filter set: stopping after merging\n")
		color.Green("Final VCF: %s\n", jointVCF)
		return jointVCF, nil
	}

	// ==================================== Stage 3: filter ===================================== //

	color.Cyan("\n%s STAGE 3/3  HARD FILTER %s\n\n", strings.Repeat("=", 24), strings.Repeat("=", 24))

	filteredVCF, err := FilterVcf(opts, jointVCF)
	if err != nil {
		return "", fmt.Errorf("hard filtering %s: %w", jointVCF, err)
	}

	color.Green("Final VCF: %s\n", filteredVCF)
	return filteredVCF, nil
}
