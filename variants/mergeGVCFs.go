package variants

import (
	"bufio"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

func gvcfSampleName(gvcf string) ([]string, error) {
	in, cleanup, err := openVCF(gvcf)
	if err != nil {
		return nil, fmt.Errorf("open %q: %w", gvcf, err)
	}
	defer cleanup()

	scanner := bufio.NewScanner(in)
	// header lines can be long when there are many samples; grow the buffer
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) <= 9 {
				return []string{}, fmt.Errorf("No sample columns") // sites-only VCF, no sample columns
			}
			return fields[9:], nil
		}
		if !strings.HasPrefix(line, "#") {
			break // hit variant records without ever seeing #CHROM
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanning %s: %w", gvcf, err)
	}
	return nil, fmt.Errorf("%s: no #CHROM header line found", gvcf)

}

func allGvcfSampleNames(gvcfs []string) ([]string, error) {
	var allNames []string
	for _, gvcf := range gvcfs {
		names, err := gvcfSampleName(gvcf)
		if err != nil {
			return nil, err
		}
		for _, name := range names {
			allNames = append(allNames, name)
		}

	}
	return allNames, nil
}

func MergeGvcfsGATKDir(gvcfs []string, refFile string, chrom string, species string, refVer string, verbose bool, logLevel string, outDir string, quick bool) (string, error) {
	chromDirName := "contigs"
	if chrom != "contigs" {
		chromDirName = strings.ReplaceAll(chrom, ".", "_")
	}
	chromDir := filepath.Join(outDir, chromDirName)
	theDB := filepath.Join(chromDir, chrom+"DB")
	tmpDir := filepath.Join(chromDir, "tmp")
	tmpDir2 := filepath.Join(chromDir, "tmp2")
	vcfDir := filepath.Join(chromDir, "VCFs")

	for _, dir := range []string{tmpDir, tmpDir2, vcfDir} {
		if _, err := os.Stat(dir); os.IsNotExist(err) {
			if cErr := os.MkdirAll(dir, 0755); cErr != nil {
				return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
			}
		}
	}

	jointVCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.vcf.gz")

	_, vcfErr := os.Stat(jointVCF)
	if vcfErr == nil {
		fmt.Printf("merged VCF %s already exists. Validating ...\n", vcfDir)
		vErr := utils.ValidateGvcf(jointVCF, verbose, quick)
		if vErr != nil {
			color.Yellow("merged VCF %s is corrupted: %v", jointVCF, vErr)
			color.Yellow("re-merging ...")
		}

	}

	// ── GenomicsDBImport ────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGenomicsDBImport, "ALL", chrom) {
		if rErr := os.RemoveAll(theDB); rErr != nil {
			return "", fmt.Errorf("removing %s: %w", theDB, rErr)
		}

		if len(gvcfs) == 0 {
			return "", fmt.Errorf("no gVCFs to merge for %s", chrom)
		}
		// Write the sample map alongside the gVCFs. Their directory is caller-specific
		// (gatk_gvcfs / dv_gvcfs), so derive it from the gVCF paths rather than assuming.
		sampleMap, err := CreateSampleMap(gvcfs, filepath.Dir(gvcfs[0]))
		if err != nil {
			return "", fmt.Errorf("creating sample map: %w", err)
		}

		gDBCmd := fmt.Sprintf(
			`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport --sample-name-map %s `+
				`--genomicsdb-workspace-path %s --tmp-dir %s -L %s `+
				`--genomicsdb-shared-posixfs-optimizations true --batch-size 50 `+
				`--bypass-feature-reader --verbosity %s`,
			sampleMap, theDB, tmpDir, chrom, logLevel,
		)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", gDBCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(gDBCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("gatk GenomicsDBImport: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenomicsDBImport, chrom))
	}

	// ── GenotypeGVCFs ───────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGenotypeGVCFs, "ALL", chrom) {
		genoCmd := fmt.Sprintf(
			`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
			refFile, theDB, jointVCF, tmpDir2, logLevel,
		)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", genoCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(genoCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("gatk GenotypeGVCFs: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenotypeGVCFs, chrom))
	}

	return jointVCF, nil
}

func MergeGvcfsGlnexusDir(gvcfs []string, chrom string, species string, refVer string, caller string, verbose bool, outDir string, logged []utils.LogEntry, jlog *slog.Logger) (string, error) {

	chromDirName := "contigs"
	if chrom != "contigs" {
		chromDirName = strings.ReplaceAll(chrom, ".", "_")
	}
	chromDir := filepath.Join(outDir, chromDirName)

	vcfDir := filepath.Join(chromDir, "VCFs")
	for _, dir := range []string{vcfDir} {
		if _, err := os.Stat(dir); os.IsNotExist(err) {
			if cErr := os.MkdirAll(dir, 0755); cErr != nil {
				return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
			}
			//fmt.Printf("Created directory: %s\n", dir)
		}
	}
	var glnexusPreset string
	if caller == "gatk" {
		glnexusPreset = "gatk"
	} else {
		glnexusPreset = "DeepVariant"
	}
	jointVCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.vcf.gz")
	jointBCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.bcf")
	gvcfFiles := strings.Join(gvcfs, " ")

	if !utils.StageHasCompleted(logged, stageGLNEXUS, "ALL", chrom) {
		glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s %s > %s`, glnexusPreset, filepath.Join(vcfDir, "GLnexus.DB"), gvcfFiles, jointBCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", glnexusCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(glnexusCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("glnexus_cli: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUS, chrom))
	}

	// ── bcftools view | bgzip ──────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGLNEXUSBCFTOOLS, "ALL", chrom) {
		bcfCmd := fmt.Sprintf("bcftools view %s | bgzip -@ 4 -c > %s", jointBCF, jointVCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", bcfCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(bcfCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("bcftools view: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSBCFTOOLS, chrom))
	}

	// ── tabix index ──────────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGLNEXUSIndex, "ALL", chrom) {
		indexCmd := fmt.Sprintf(`tabix -p vcf %s`, jointVCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", indexCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(indexCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("indexing joint VCF with tabix: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSIndex, chrom))
	}

	return jointVCF, nil
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

	color.Green("Checking gVCF presence and integrity for all samples × chroms ...\n\n")

	type missingEntry struct {
		sample string
		chrom  string
		reason string // "missing", "corrupted", "multiple"
	}

	badChromsMap := make(map[string][]missingEntry)
	for _, chrom := range chroms {
		var missingGvcfs []missingEntry
		for _, sample := range samples {
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
		if len(missingGvcfs) > 0 {
			badChromsMap[chrom.ID] = missingGvcfs
		}
	}

	// ============================================= Bail if anything is missing ================================= //

	if len(badChromsMap) > 0 {
		color.Cyan("========================== Skipping the following chromosomes with missing sample gvcfs: %v\n=============================\n ")
		for c, entries := range badChromsMap {
			var missingSamples []string
			for _, entry := range entries {
				missingSamples = append(missingSamples, entry.sample)
			}
			color.Yellow("Samples with missing/corrupted gvcfs for chromosome [%s]:\n ", c)
			color.Yellow("%s\n---------------------------------------------------------------------\n\n", missingSamples)

		}
	}

	// ==================================== Get

	for _, chrom := range chroms {
		if entries, ok := badChromsMap[chrom.ID]; ok {
			color.Red("Skipping chromosome %s due to missing/corrupted gVCFs for samples: %v\n", chrom.ID, entries)
			continue
		}
		switch merger {
		case "gatk":

		}

	}

	//if len(missingGvcfs) > 0 {
	//	color.Red("\n\nCannot proceed with merging. The following gVCFs are missing or invalid:\n")
	//	color.Red("%-30s %-20s %s\n", "SAMPLE", "CHROM", "REASON")
	//	color.Red("%s\n", strings.Repeat("-", 60))
	//	for _, m := range missingGvcfs {
	//		color.Red("%-30s %-20s %s\n", m.sample, m.chrom, m.reason)
	//	}
	//	fmt.Printf("\nTotal issues: %d\n", len(missingGvcfs))
	//	return
	//}
	//
	//color.Green("\nAll gVCFs present and valid. Proceeding with merge ...\n\n")

	// ============================================= Merge ======================================================= //
	//TODO:
	// a) verify missing gvcfs for each sample, continue for samples that are complete but skip if r
	//1. Check VCFs/<ver>/ directiory for latest ALL.vcf.gz file (from timestamped directory)
	//2. Use bcftools to check the sample names (or vcfgo)
	//3. check sample names of current scan and check if there are new samples
	//4. If new samples create another time stamp
	//5. If gatk caller merge with GATK
	//6. If deepvariant merge with glnexus

	// switch merger {
	// case "gatk":
	// 	fmt.Println("Using GATK MergeVcfs")
	// 	// TODO: implement GATK merging

	// case "glnexus":
	// 	fmt.Println("Using GLnexus merge")

	// 	for _, chrom := range chroms {
	// 		sID := chrom.ID

	// 		speciesDir := filepath.Join(dataDirAbs, strings.ToLower(species))
	// 		speciesUpper := strings.ToUpper(species)
	// 		glnexusCmdStr := glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper)

	// 		var glErr error
	// 		if verbose {
	// 			glErr = utils.RunBashCmdVerbose(glnexusCmdStr)
	// 		} else {
	// 			glErr = utils.RunBashCmd(glnexusCmdStr)
	// 		}

	// 		if glErr != nil {
	// 			color.Red("GLnexus FAILED for chrom %s: %v\n", sID, glErr)
	// 			return
	// 		}

	// 		color.Green("GLnexus completed for chrom %s\n", sID)
	// 		fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
	// 	}

	// default:
	// 	fmt.Println("Please provide a valid merger: either gatk or glnexus")
	// }
}
