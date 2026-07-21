package variants

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

//func glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper string) string {
//	lines := []string{
//		"cd ~/", "rm -r GLnexus.DB",
//		fmt.Sprintf("glnexus_cli --config gatk %s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz > %s/%s.%s.joint.bcf", speciesDir, refVer, sID, outDir, speciesUpper, sID),
//		fmt.Sprintf("bcftools view %s/%s.%s.joint.bcf | bgzip -@ 4 -c > %s/%s.%s.joint.vcf.gz", outDir, speciesUpper, sID, outDir, speciesUpper, sID),
//	}
//	return strings.Join(lines, "\n")
//}

func glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper string) string {
	lines := []string{
		"set -euo pipefail",
		"cd ~/",
		"rm -rf GLnexus.DB",
		fmt.Sprintf("glnexus_cli --config gatk %s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz > %s/%s.%s.joint.bcf", speciesDir, refVer, sID, outDir, speciesUpper, sID),
		fmt.Sprintf("bcftools view %s/%s.%s.joint.bcf | bgzip -@ 4 -c > %s/%s.%s.joint.vcf.gz", outDir, speciesUpper, sID, outDir, speciesUpper, sID),
	}
	return strings.Join(lines, " && \\\n")
}

func MergeGvcfs(config string, gvcfs []string, dataDir string, species string, refVer string, refFasta string, outDir string, caller string, merger string, verbose bool, quick bool, skipVerification bool) {
	fmt.Println("Merging GVCFs ...")
	if config != "" {
		fmt.Println("Using config file: ", config)
		// return
	}

	if len(gvcfs) > 0 {
		fmt.Println("Using gvcfs: ", gvcfs)
		// return
	}

	println("Merging GVCFs using plennegy data structure...")

	// ============================================= Validate paths ================================================ //

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		fmt.Printf("Error getting absolute path for data directory: %s\n", dataDir)
		return
	}

	if species == "" {
		fmt.Println("Please provide species name")
		return
	}
	if refVer == "" {
		fmt.Println("Please provide reference version name")
		return
	}
	if refFasta == "" {
		fmt.Println("Please provide reference fasta path")
		return
	}

	fastaInfo, err := os.Stat(refFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", refFasta)
		return
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", refFasta)
		return
	}

	dictFilePath := refFasta[:len(refFasta)-len(filepath.Ext(refFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	if outDir == "" {
		outDir = filepath.Join(dataDirAbs, species, "VCFs", refVer)
		fmt.Printf("No output directory provided. Using: %s\n", outDir)
	}
	if mkErr := os.MkdirAll(outDir, 0755); mkErr != nil {
		fmt.Printf("Failed to create output directory %s: %v\n", outDir, mkErr)
		return
	}

	color.Green("All file paths valid\n....................................................\n\n")

	// ============================================= Get chroms ================================================== //

	chroms, _, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		fmt.Printf("Error getting chromosomes: %v\n", err)
		return
	}

	// ============================================= Discover samples ============================================= //

	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		panic(err)
	}

	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, s)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n", len(samples), species)

	if len(samples) == 0 {
		color.Red("No samples found. Exiting.")
		return
	}

	// ============================================= Check & validate chrom gVCFs ================================ //
	//
	// For every (sample × chrom) pair, mirror VariantCallingDir's logic:
	//   0 files  → record as missing
	//   1 file   → validate; record as missing if corrupt
	//   2+ files → record as ambiguous (multiple)

	color.Green("Checking gVCF presence and integrity for all samples × chroms ...\n\n")

	type missingEntry struct {
		sample string
		chrom  string
		reason string // "missing", "corrupted", "multiple"
	}

	var missingGvcfs []missingEntry

	// Build the per-(sample,chrom) gvcf path the same way VariantCallingDir does, using FindBamOrVcfs.
	for _, sample := range samples {
		for _, chrom := range chroms {
			var vcfFiles []string
			if caller == "gatk" {
				vcfFiles, _ = FindGVCFs(dataDirAbs, species, sample, refVer, "gatk_gvcfs", chrom.ID)
			} else {
				vcfFiles, _ = FindGVCFs(dataDirAbs, species, sample, refVer, "dv_gvcfs", chrom.ID)
			}

			switch len(vcfFiles) {
			case 0:
				color.Red("[%s] gVCF MISSING for chrom %s\n", sample, chrom.ID)
				missingGvcfs = append(missingGvcfs, missingEntry{sample: sample, chrom: chrom.ID, reason: "missing"})

			case 1:
				vcf := vcfFiles[0]
				color.Green("[%s] gVCF found for chrom %s: %s\n", sample, chrom.ID, vcf)
				if skipVerification {
					color.Yellow("[%s] skipping integrity check for %s\n", sample, color.BlueString(vcf))
				} else {
					fmt.Printf("[%s] checking integrity of gVCF file: %s ..........\n", sample, color.BlueString(vcf))
					if vvvErr := utils.ValidateGvcf(vcf, verbose, quick); vvvErr != nil {
						color.Red("[%s] gVCF %s corrupted: %v\n", sample, color.BlueString(vcf), vvvErr)
						missingGvcfs = append(missingGvcfs, missingEntry{sample: sample, chrom: chrom.ID, reason: "corrupted"})
					} else {
						color.Green("[%s] gVCF %s is valid!!\n\n", sample, vcf)
					}
				}

			default:
				color.Red("[%s] Multiple gVCF files found for chrom %s — please remove the extra file(s): %v\n\n", sample, chrom.ID, vcfFiles)
				missingGvcfs = append(missingGvcfs, missingEntry{sample: sample, chrom: chrom.ID, reason: "multiple"})
			}
		}
	}

	// ============================================= Bail if anything is missing ================================= //

	if len(missingGvcfs) > 0 {
		color.Red("\n\nCannot proceed with merging. The following gVCFs are missing or invalid:\n")
		color.Red("%-30s %-20s %s\n", "SAMPLE", "CHROM", "REASON")
		color.Red("%s\n", strings.Repeat("-", 60))
		for _, m := range missingGvcfs {
			color.Red("%-30s %-20s %s\n", m.sample, m.chrom, m.reason)
		}
		fmt.Printf("\nTotal issues: %d\n", len(missingGvcfs))
		return
	}

	color.Green("\nAll gVCFs present and valid. Proceeding with merge ...\n\n")

	// ============================================= Merge ======================================================= //
	//TODO:
	// a) verify missing gvcfs for each sample, continue for samples that are complete but skip if r
	//1. Check VCFs/<ver>/ directiory for latest ALL.vcf.gz file (from timestamped directory)
	//2. Use bcftools to check the sample names (or vcfgo)
	//3. check sample names of current scan and check if there are new samples
	//4. If new samples create another time stamp
	//5. If gatk caller merge with GATK
	//6. If deepvariant merge with glnexus

	switch merger {
	case "gatk":
		fmt.Println("Using GATK MergeVcfs")
		// TODO: implement GATK merging

	case "glnexus":
		fmt.Println("Using GLnexus merge")

		for _, chrom := range chroms {
			sID := chrom.ID

			speciesDir := filepath.Join(dataDirAbs, strings.ToLower(species))
			speciesUpper := strings.ToUpper(species)
			glnexusCmdStr := glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper)

			var glErr error
			if verbose {
				glErr = utils.RunBashCmdVerbose(glnexusCmdStr)
			} else {
				glErr = utils.RunBashCmd(glnexusCmdStr)
			}

			if glErr != nil {
				color.Red("GLnexus FAILED for chrom %s: %v\n", sID, glErr)
				return
			}

			color.Green("GLnexus completed for chrom %s\n", sID)
			fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
		}

	default:
		fmt.Println("Please provide a valid merger: either gatk or glnexus")
	}
}
