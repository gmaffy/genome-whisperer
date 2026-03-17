package annotation

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"regexp"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/schollz/progressbar/v3"

	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

func checkSnpEffDB(db string) error {
	getDbCmdStr := fmt.Sprintf(`snpEff databases | grep "%s" | awk '{print $1}'`, db)
	getDbCmd := exec.Command("bash", "-c", getDbCmdStr)

	var stdoutBuf bytes.Buffer
	getDbCmd.Stdout = &stdoutBuf
	getDbCmd.Stderr = os.Stderr

	getDbErr := getDbCmd.Run()
	if getDbErr != nil {
		fmt.Printf("could not get database: %s\n", getDbErr)
		if exitError, ok := getDbErr.(*exec.ExitError); ok {

			if exitError.ExitCode() == 1 {
				fmt.Printf("No database matching '%s' found by grep.\n", db)
				return getDbErr
			}
		}
		return getDbErr

	}
	output := stdoutBuf.String()
	cleanedOutput := strings.TrimSpace(output)

	if cleanedOutput == "" {
		fmt.Printf("No database matching '%s' found in snpEff databases output.\n", db)
		return fmt.Errorf("No database matching '%s' found in snpEff databases output.\n", db)
	}

	matches := strings.Split(cleanedOutput, "\n")

	for _, match := range matches {
		trimmedMatch := strings.TrimSpace(match)
		if trimmedMatch == db {
			fmt.Printf("Database '%s' found in snpEff databases output.\n", db)
			return nil
		}

	}
	fmt.Printf("Here are some matches from installed databases:\n\n")
	for _, match := range matches {
		fmt.Printf("%s\n", match)
	}
	fmt.Printf("\n")
	return fmt.Errorf("Database '%s' not found in snpEff databases output.\n", db)

}

func RunSnpEff(vcfs []string, db string, bsaseq bool) (error, []string, []string) {
	// --------------------------------------------------- Get Config file ---------------------------------------- //
	snEffPath, err := exec.LookPath("snpEff")
	if err != nil {
		return fmt.Errorf("snpEff not found: %w", err), []string{}, []string{}
	}

	scriptsDir := filepath.Dir(snEffPath)
	snpEffDir := filepath.Dir(scriptsDir)
	configPath := filepath.Join(snpEffDir, "snpEff.config")

	_, rErr := os.Stat(configPath)
	if rErr != nil {
		fmt.Printf("Tried to find snpEff config file at: %s. It does not exist\n\n", configPath)
		return rErr, []string{}, []string{}
	}
	//fmt.Printf("snpEff config file found at %s\n", configPath)

	// --------------------------------------------- Check if database exists --------------------------------------- //
	fmt.Printf("Checking if database %s is installed ...\n\n", db)
	dbErr := checkSnpEffDB(db)
	if dbErr != nil {
		fmt.Printf("No exact match for %s found in snpEff databases\n\n", db)
		dbDir := filepath.Join(snpEffDir, "data")
		dbs, dataErr := os.ReadDir(dbDir)
		if dataErr != nil {
			fmt.Printf("You have not installed any custom databases\n\n")
			fmt.Printf("Please install a custom database or look for installed databases by running: snpEff databases\n\n")
		} else {
			fmt.Printf("Here are the databases you have installed yourself \n\n")
			for _, dir := range dbs {

				if dir.IsDir() {
					fmt.Printf("%s\n", dir.Name())
				}
			}
			fmt.Printf("\n")
			fmt.Printf("Please install a custom database or look for installed databases by running: snpEff databases\n\n")
		}
		return dbErr, []string{}, []string{}
	}
	fmt.Printf("DATABASE FOUND: %s\n", db)

	// ------------------------------------------------ Run SnpEff -------------------------------------------------- //
	var tsvFiles []string
	var vcfFiles []string
	for _, vcf := range vcfs {
		var snpEffVcf string
		var snpEffTsv string

		if strings.HasSuffix(vcf, ".vcf.gz") {
			snpEffVcf = strings.TrimSuffix(vcf, ".vcf.gz") + ".snpEff.vcf"
			snpEffTsv = strings.TrimSuffix(vcf, ".vcf.gz") + ".snpEff.tsv"
		} else if strings.HasSuffix(vcf, ".vcf") {
			snpEffVcf = strings.TrimSuffix(vcf, ".vcf") + ".snpEff.vcf"
			snpEffTsv = strings.TrimSuffix(vcf, ".vcf") + ".snpEff.tsv"
		} else {
			fmt.Println("vcf file must be in vcf or vcf.gz format")
			return fmt.Errorf("vcf file must be in vcf or vcf.gz format"), []string{}, []string{}
		}

		snpEffCmdStr := fmt.Sprintf(`snpEff -c %s -v -o gatk %s %s > %s`, configPath, db, vcf, snpEffVcf)
		fmt.Println(snpEffCmdStr)
		err := utils.RunBashCmdVerbose(snpEffCmdStr)
		if err != nil {
			return err, []string{}, []string{}
		}

		fmt.Println("Converting vcf to table")
		var vtCmdStr string
		if bsaseq {
			vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`, snpEffVcf, snpEffTsv)
		} else {
			vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`, snpEffVcf, snpEffTsv)
		}
		fmt.Println(vtCmdStr)
		terr := utils.RunBashCmdVerbose(vtCmdStr)
		if terr != nil {
			return terr, []string{}, []string{}
		}

		sErr, splitSnpEffTsv := splitEffColumns(snpEffTsv)
		if sErr != nil {
			return sErr, []string{}, []string{}
		}
		tsvFiles = append(tsvFiles, splitSnpEffTsv)
		vcfFiles = append(vcfFiles, snpEffVcf)

	}

	return nil, tsvFiles, vcfFiles
}

type SnpEffEffect struct {
	Effect            string
	Impact            string
	FunctionalClass   string
	CodonChange       string
	AminoAcidChange   string
	GeneName          string
	TranscriptBiotype string
	GeneCoding        string
	TranscriptID      string
	Errors            string
	Warnings          string
}

func parseEffColumn(effColValue string) SnpEffEffect {
	eff := SnpEffEffect{}

	matchEffect := regexp.MustCompile(`^([^(]+)`).FindStringSubmatch(effColValue)
	if len(matchEffect) > 1 {
		eff.Effect = strings.TrimSpace(matchEffect[1])
	}

	var innerContent string
	matchInside := regexp.MustCompile(`\(([^)]*)\)`).FindStringSubmatch(effColValue)
	if len(matchInside) > 1 {
		innerContent = matchInside[1]
	}

	parts := strings.Split(innerContent, "|")

	fieldPointers := []*string{
		&eff.Impact, &eff.FunctionalClass, &eff.CodonChange,
		&eff.AminoAcidChange, &eff.GeneName, &eff.TranscriptBiotype,
		&eff.GeneCoding, &eff.TranscriptID, &eff.Errors, &eff.Warnings,
	}

	for i, part := range parts {
		if i < len(fieldPointers) {
			*fieldPointers[i] = strings.TrimSpace(part)
		}
	}

	return eff
}

func splitEffColumns(effFile string) (error, string) {
	fmt.Printf("Splitting EFF column in file: %s \n\n ", effFile)
	inputFile, err := os.Open(effFile)
	if err != nil {
		return fmt.Errorf("failed to open input file %s: %w", effFile, err), ""
	}
	defer inputFile.Close()
	scanner := bufio.NewScanner(inputFile)

	if !scanner.Scan() {
		return fmt.Errorf("input file %s is empty or has no header", effFile), ""
	}
	headerLine := scanner.Text()
	originalHeaders := strings.Split(headerLine, "\t")

	effColIndex := -1
	var otherHeaders []string
	for i, h := range originalHeaders {
		if h == "EFF" {
			effColIndex = i
		} else {
			otherHeaders = append(otherHeaders, h)
		}
	}

	if effColIndex == -1 {
		return fmt.Errorf("EFF column not found in header: %s", headerLine), ""
	}

	snpEffNewHeaders := []string{
		"SNPEFF_EFFECT", "SNPEff_IMPACT", "SNPEFF_FUNCTIONAL_CLASS",
		"SNPEFF_CODON_CHANGE", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME",
		"SNPEFF_TRANSCRIPT_BIOTYPE", "SNPEFF_GENE_CODING",
		"SNPEFF_TRANSCRIPT_ID", "ERRORS", "WARNINGS",
	}

	newHeader := strings.Join(otherHeaders, "\t") + "\t" + strings.Join(snpEffNewHeaders, "\t")

	outputFileName := strings.Replace(effFile, ".tsv", "_EFF.tsv", 1)
	if !strings.HasSuffix(effFile, ".tsv") {
		outputFileName = effFile + "_EFF.tsv"
	}

	fmt.Printf("Writing to file: %s ...\n", outputFileName)
	outputFile, err := os.Create(outputFileName)
	if err != nil {
		return fmt.Errorf("failed to create output file %s: %w", outputFileName, err), ""
	}
	defer outputFile.Close()

	writer := bufio.NewWriter(outputFile)
	_, err = writer.WriteString(newHeader + "\n")
	if err != nil {
		return fmt.Errorf("failed to write header to output file: %w", err), ""
	}

	lineNum := 1
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()
		fields := strings.Split(line, "\t")

		if len(fields) != len(originalHeaders) {
			fmt.Fprintf(os.Stderr, "Warning: Skipping line %d due to inconsistent column count. Expected %d, got %d. Line: %s\n",
				lineNum, len(originalHeaders), len(fields), line)
			continue
		}

		effColValue := fields[effColIndex]
		snpEffData := parseEffColumn(effColValue)

		var outputFields []string
		for i, fieldVal := range fields {
			if i != effColIndex { // Add all original fields except the "EFF" column
				outputFields = append(outputFields, fieldVal)
			}
		}

		outputFields = append(outputFields,
			snpEffData.Effect,
			snpEffData.Impact,
			snpEffData.FunctionalClass,
			snpEffData.CodonChange,
			snpEffData.AminoAcidChange,
			snpEffData.GeneName,
			snpEffData.TranscriptBiotype,
			snpEffData.GeneCoding,
			snpEffData.TranscriptID,
			snpEffData.Errors,
			snpEffData.Warnings,
		)

		_, err = writer.WriteString(strings.Join(outputFields, "\t") + "\n")
		if err != nil {
			return fmt.Errorf("failed to write line %d to output file: %w", lineNum, err), ""
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading input file: %w", err), ""
	}

	err = writer.Flush()
	if err != nil {
		return err, ""
	}
	return nil, outputFileName
}

func AddDescriptions(variants []string, desc string, bsaseq bool) (error, []string) {
	fmt.Println("Variants: ", variants)
	fmt.Println("Gene Description File: ", desc)

	// ========================================== Read description file ============================================= //
	fmt.Printf("Reading description file: %s ...............................\n", desc)

	descFile, err := os.Open(desc)
	if err != nil {
		return fmt.Errorf("failed to open description file %s: %w", desc, err), nil
	}
	defer descFile.Close()

	scanner := bufio.NewScanner(descFile)

	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	geneDesc := make(map[string]string)

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")

		if len(fields) < 2 {
			return fmt.Errorf("description file %s invalid format", desc), nil
		}

		geneDesc[fields[0]] = fields[1]
	}

	fmt.Printf("Loaded %d gene descriptions\n-------------------------------------------------\n\n", len(geneDesc))

	// ============================================ Process variant files =========================================== //

	var descTsvFiles []string

	for _, variantFile := range variants {

		var snpEffTsv io.Reader
		var outputFileName string

		// -------------------------------------------- Input detection --------------------------------------------- //

		if strings.HasSuffix(variantFile, ".tsv") || strings.HasSuffix(variantFile, ".txt") {
			fmt.Printf("Processing tsv file: %s\n\n", variantFile)

			f, err := os.Open(variantFile)
			if err != nil {
				return fmt.Errorf("failed to open %s: %w", variantFile, err), nil
			}
			defer f.Close()

			snpEffTsv = f
			outputFileName = strings.TrimSuffix(strings.TrimSuffix(variantFile, ".tsv"), ".txt") + "_DESC.tsv"

		} else if strings.HasSuffix(variantFile, ".vcf") || strings.HasSuffix(variantFile, ".vcf.gz") {
			fmt.Printf("Processing vcf file: %s\n\n", variantFile)
			fmt.Println("Converting VCF to table")

			var variantTsvFile string

			if strings.HasSuffix(variantFile, ".vcf") {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf") + ".tsv"
			} else {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf.gz") + ".tsv"
			}

			var vtCmdStr string

			if bsaseq {
				vtCmdStr = fmt.Sprintf(
					`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`,
					variantFile, variantTsvFile)
			} else {
				vtCmdStr = fmt.Sprintf(
					`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`,
					variantFile, variantTsvFile)
			}

			if err := utils.RunBashCmdVerbose(vtCmdStr); err != nil {
				return err, nil
			}

			sErr, splitFile := splitEffColumns(variantTsvFile)
			if sErr != nil {
				return sErr, nil
			}

			f, err := os.Open(splitFile)
			if err != nil {
				return fmt.Errorf("failed to open %s: %w", splitFile, err), nil
			}
			defer f.Close()

			snpEffTsv = f
			outputFileName = strings.TrimSuffix(splitFile, ".tsv") + "_DESC.tsv"

		} else {
			return fmt.Errorf("invalid variant file format: %s", variantFile), nil
		}

		// ====================================================== Output ============================================ //

		out, err := os.Create(outputFileName)
		if err != nil {
			return err, nil
		}
		defer out.Close()

		writer := bufio.NewWriter(out)

		innerScanner := bufio.NewScanner(snpEffTsv)
		innerScanner.Buffer(buf, 10*1024*1024)

		// ----------------------------------------- Check variant file header ------------------------------------- //

		if !innerScanner.Scan() {
			return fmt.Errorf("file %s empty", variantFile), nil
		}

		headerLine := innerScanner.Text()
		headers := strings.Split(headerLine, "\t")

		geneIdx := -1
		transIdx := -1

		for i, h := range headers {
			if h == "SNPEFF_GENE_NAME" {
				geneIdx = i
			}
			if h == "SNPEFF_TRANSCRIPT_ID" {
				transIdx = i
			}
		}

		if geneIdx == -1 && transIdx == -1 {
			return fmt.Errorf("neither SNPEFF_GENE_NAME nor SNPEFF_TRANSCRIPT_ID found"), nil
		}

		// ----------------------------- Read all rows to compare gene/transcript notation -------------------------- //

		fmt.Println("Comparing gene/transcript names in variant file against ones in description file ...")

		var allRows [][]string
		genesInVariants := make(map[string]struct{})
		transInVariants := make(map[string]struct{})

		readBar := progressbar.Default(-1, "reading rows")
		for innerScanner.Scan() {
			fields := strings.Split(innerScanner.Text(), "\t")
			allRows = append(allRows, fields)

			if geneIdx >= 0 && geneIdx < len(fields) {
				genesInVariants[fields[geneIdx]] = struct{}{}
			}
			if transIdx >= 0 && transIdx < len(fields) {
				transInVariants[fields[transIdx]] = struct{}{}
			}
			err := readBar.Add(1)
			if err != nil {
				return err, nil
			}
		}
		err3 := readBar.Finish()
		if err3 != nil {
			return err3, nil
		}

		// Count overlaps with the description map
		geneOverlap := 0
		for g := range genesInVariants {
			if _, ok := geneDesc[g]; ok {
				geneOverlap++
			}
		}

		transOverlap := 0
		for t := range transInVariants {
			if _, ok := geneDesc[t]; ok {
				transOverlap++
			}
		}

		// Decide lookup strategy, mirroring the Python priority logic
		var useGeneIdx int
		switch {
		case geneOverlap > 0:
			fmt.Println("YOUR GENE DESCRIPTION FILE HAS GENE NAMES ... using SNPEFF_GENE_NAME ...")
			useGeneIdx = geneIdx
		case transOverlap > 0:
			fmt.Println("YOUR GENE DESCRIPTION FILE HAS TRANSCRIPT NAMES ... using SNPEFF_TRANSCRIPT_ID ...")
			useGeneIdx = transIdx
		default:
			fmt.Printf("YOUR GENE_DESCRIPTION FILE %s HAS GENES THAT ARE IN A DIFFERENT FORMAT TO THE FORMAT OF YOUR VCF GENES ...\n", desc)
			return fmt.Errorf("no overlap between variant gene/transcript IDs and description file keys"), nil
		}

		// ----------------------------------------------- Write output -------------------------------------------- //

		fmt.Printf("\nWriting output to %s\n\n", outputFileName)

		writer.WriteString(headerLine + "\tGENE_DESC\n")
		writeBar := progressbar.Default(int64(len(allRows)), "writing output")
		for _, fields := range allRows {
			description := "NA"

			if useGeneIdx >= 0 && useGeneIdx < len(fields) {
				if d, ok := geneDesc[fields[useGeneIdx]]; ok {
					description = d
				}
			}

			fields = append(fields, description)
			writer.WriteString(strings.Join(fields, "\t") + "\n")
			writeBar.Add(1)
		}
		writeBar.Finish()

		writer.Flush()

		descTsvFiles = append(descTsvFiles, outputFileName)

		fmt.Printf("Created: %v\n\n", outputFileName)
	}

	return nil, descTsvFiles
}

func AddPrg(variants []string, prgBlastFile string, bsaseq bool) (error, []string) {
	fmt.Println("Variants: ", variants)
	fmt.Println("PRG Blast File: ", prgBlastFile)

	// ========================================== Read PRG blast file ============================================= //
	fmt.Printf("Reading PRG blast file: %s ...............................\n", prgBlastFile)

	prgFile, err := os.Open(prgBlastFile)
	if err != nil {
		return fmt.Errorf("failed to open PRG blast file %s: %w", prgBlastFile, err), nil
	}
	defer prgFile.Close()

	scanner := bufio.NewScanner(prgFile)
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	// ----------------------------------------- Read PRG header ------------------------------------------------- //

	if !scanner.Scan() {
		return fmt.Errorf("PRG blast file %s is empty", prgBlastFile), nil
	}

	prgHeaders := strings.Split(scanner.Text(), "\t")

	qseqIdx := -1
	percIdentIdx := -1
	lengthIdx := -1
	qlenIdx := -1
	evalIdx := -1

	for i, h := range prgHeaders {
		switch h {
		case "QseqID":
			qseqIdx = i
		case "PercIdent":
			percIdentIdx = i
		case "Length":
			lengthIdx = i
		case "Qlen":
			qlenIdx = i
		case "Eval":
			evalIdx = i
		}
	}

	missingCols := []string{}
	if qseqIdx == -1 {
		missingCols = append(missingCols, "QseqID")
	}
	if percIdentIdx == -1 {
		missingCols = append(missingCols, "PercIdent")
	}
	if lengthIdx == -1 {
		missingCols = append(missingCols, "Length")
	}
	if qlenIdx == -1 {
		missingCols = append(missingCols, "Qlen")
	}
	if evalIdx == -1 {
		missingCols = append(missingCols, "Eval")
	}
	if len(missingCols) > 0 {
		return fmt.Errorf("PRG blast file missing required columns: %s", strings.Join(missingCols, ", ")), nil
	}

	// ----------------------------------------- Read PRG rows ---------------------------------------------------- //

	type prgEntry struct {
		percIdent    string
		qlenMatchLen string
		eval         string
	}

	prgMap := make(map[string]prgEntry)

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")

		if qseqIdx >= len(fields) || lengthIdx >= len(fields) || qlenIdx >= len(fields) {
			continue
		}

		qseqID := fields[qseqIdx]
		qlenMatchLen := fields[lengthIdx] + "/" + fields[qlenIdx]

		prgMap[qseqID] = prgEntry{
			percIdent:    fields[percIdentIdx],
			qlenMatchLen: qlenMatchLen,
			eval:         fields[evalIdx],
		}
	}

	fmt.Printf("Loaded %d PRG entries\n-------------------------------------------------\n\n", len(prgMap))

	// ============================================ Process variant files ========================================== //

	var prgTsvFiles []string

	for _, variantFile := range variants {

		var snpEffTsv io.Reader
		var outputFileName string

		// -------------------------------------------- Input detection ------------------------------------------- //

		if strings.HasSuffix(variantFile, ".tsv") || strings.HasSuffix(variantFile, ".txt") {
			fmt.Printf("Processing tsv file: %s\n\n", variantFile)

			f, err := os.Open(variantFile)
			if err != nil {
				return fmt.Errorf("failed to open %s: %w", variantFile, err), nil
			}
			defer f.Close()

			snpEffTsv = f
			outputFileName = strings.TrimSuffix(strings.TrimSuffix(variantFile, ".tsv"), ".txt") + "_PRG.tsv"

		} else if strings.HasSuffix(variantFile, ".vcf") || strings.HasSuffix(variantFile, ".vcf.gz") {
			fmt.Printf("Processing vcf file: %s\n\n", variantFile)
			fmt.Println("Converting VCF to table ...")

			var variantTsvFile string
			if strings.HasSuffix(variantFile, ".vcf") {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf") + ".tsv"
			} else {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf.gz") + ".tsv"
			}

			var vtCmdStr string
			if bsaseq {
				vtCmdStr = fmt.Sprintf(
					`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`,
					variantFile, variantTsvFile)
			} else {
				vtCmdStr = fmt.Sprintf(
					`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`,
					variantFile, variantTsvFile)
			}

			if err := utils.RunBashCmdVerbose(vtCmdStr); err != nil {
				return err, nil
			}

			sErr, splitFile := splitEffColumns(variantTsvFile)
			if sErr != nil {
				return sErr, nil
			}

			f, err := os.Open(splitFile)
			if err != nil {
				return fmt.Errorf("failed to open %s: %w", splitFile, err), nil
			}
			defer f.Close()

			snpEffTsv = f
			outputFileName = strings.TrimSuffix(splitFile, ".tsv") + "_PRG.tsv"

		} else {
			return fmt.Errorf("invalid variant file format: %s", variantFile), nil
		}

		// ----------------------------------------- Check variant file header ------------------------------------ //

		innerScanner := bufio.NewScanner(snpEffTsv)
		innerScanner.Buffer(buf, 10*1024*1024)

		if !innerScanner.Scan() {
			return fmt.Errorf("file %s is empty", variantFile), nil
		}

		headerLine := innerScanner.Text()
		headers := strings.Split(headerLine, "\t")

		geneIdx := -1
		transIdx := -1

		for i, h := range headers {
			if h == "SNPEFF_GENE_NAME" {
				geneIdx = i
			}
			if h == "SNPEFF_TRANSCRIPT_ID" {
				transIdx = i
			}
		}

		if geneIdx == -1 && transIdx == -1 {
			return fmt.Errorf("neither SNPEFF_GENE_NAME nor SNPEFF_TRANSCRIPT_ID found in %s", variantFile), nil
		}

		// ----------------------------- Read all rows and collect gene/transcript sets ---------------------------- //

		fmt.Println("Comparing gene/transcript names in variant file against PRG file ...")

		var allRows [][]string
		genesInVariants := make(map[string]struct{})
		transInVariants := make(map[string]struct{})

		readBar := progressbar.Default(-1, "reading rows")
		for innerScanner.Scan() {
			fields := strings.Split(innerScanner.Text(), "\t")
			allRows = append(allRows, fields)

			if geneIdx >= 0 && geneIdx < len(fields) {
				genesInVariants[fields[geneIdx]] = struct{}{}
			}
			if transIdx >= 0 && transIdx < len(fields) {
				transInVariants[fields[transIdx]] = struct{}{}
			}
			readBar.Add(1)
		}
		readBar.Finish()

		// ----------------------------- Check overlap and decide lookup column ------------------------------------ //

		geneOverlap := 0
		for g := range genesInVariants {
			if _, ok := prgMap[g]; ok {
				geneOverlap++
			}
		}

		transOverlap := 0
		for t := range transInVariants {
			if _, ok := prgMap[t]; ok {
				transOverlap++
			}
		}

		var useIdx int
		switch {
		case geneOverlap > 0:
			fmt.Println("PRG file uses GENE NAMES ... using SNPEFF_GENE_NAME ...")
			useIdx = geneIdx
		case transOverlap > 0:
			fmt.Println("PRG file uses TRANSCRIPT IDs ... using SNPEFF_TRANSCRIPT_ID ...")
			useIdx = transIdx
		default:
			fmt.Println("YOUR PRG FILE GENES ARE IN A DIFFERENT FORMAT TO THE FORMAT OF YOUR VCF GENES ...")
			fmt.Println("Re format your PRG FILE gene names to match those in GATK VCF file ...")
			return fmt.Errorf("no overlap between variant gene/transcript IDs and PRG QseqID"), nil
		}

		// ----------------------------------------------- Write output ------------------------------------------- //

		out, err := os.Create(outputFileName)
		if err != nil {
			return err, nil
		}
		defer out.Close()

		writer := bufio.NewWriter(out)

		writer.WriteString(headerLine + "\tPercIdent\tQLEN/MATCH_LEN\tEval\n")

		fmt.Printf("\nWriting output to %s\n\n", outputFileName)

		writeBar := progressbar.Default(int64(len(allRows)), "writing output")
		for _, fields := range allRows {
			percIdent := "0"
			qlenMatchLen := "0"
			eval := "0"

			if useIdx >= 0 && useIdx < len(fields) {
				if entry, ok := prgMap[fields[useIdx]]; ok {
					percIdent = entry.percIdent
					qlenMatchLen = entry.qlenMatchLen
					eval = entry.eval
				}
			}

			fields = append(fields, percIdent, qlenMatchLen, eval)
			writer.WriteString(strings.Join(fields, "\t") + "\n")
			writeBar.Add(1)
		}
		writeBar.Finish()

		writer.Flush()

		prgTsvFiles = append(prgTsvFiles, outputFileName)

		fmt.Printf("\nCreated: %s\n", outputFileName)
	}

	return nil, prgTsvFiles
}

func inputDetection(variantFile string, bsaseq bool) (error, string, string) {
	fmt.Println("Variants: ", variantFile)
	// -------------------------------------------- Input detection ------------------------------------------- //
	var snpEffTsv string
	var outputFileName string
	if strings.HasSuffix(variantFile, ".tsv") || strings.HasSuffix(variantFile, ".txt") {
		fmt.Printf("Processing tsv file: %s\n\n", variantFile)
		snpEffTsv = variantFile
		outputFileName = strings.TrimSuffix(strings.TrimSuffix(variantFile, ".tsv"), ".txt") + "_PRG.tsv"

	} else if strings.HasSuffix(variantFile, ".vcf") || strings.HasSuffix(variantFile, ".vcf.gz") {
		fmt.Printf("Processing vcf file: %s\n\n", variantFile)
		fmt.Println("Converting VCF to table ...")

		var variantTsvFile string
		if strings.HasSuffix(variantFile, ".vcf") {
			variantTsvFile = strings.TrimSuffix(variantFile, ".vcf") + ".tsv"
		} else {
			variantTsvFile = strings.TrimSuffix(variantFile, ".vcf.gz") + ".tsv"
		}

		var vtCmdStr string
		if bsaseq {
			vtCmdStr = fmt.Sprintf(
				`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`,
				variantFile, variantTsvFile)
		} else {
			vtCmdStr = fmt.Sprintf(
				`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`,
				variantFile, variantTsvFile)
		}

		if err := utils.RunBashCmdVerbose(vtCmdStr); err != nil {
			return err, "", ""
		}

		sErr, splitFile := splitEffColumns(variantTsvFile)
		if sErr != nil {
			return sErr,"", ""
		}


		snpEffTsv = splitFile
		outputFileName = strings.TrimSuffix(splitFile, ".tsv") + "_PRG.tsv"

	} else {
		return fmt.Errorf("invalid variant file format: %s", variantFile), "", ""
	}

	return nil, snpEffTsv, outputFileName
}

func CreateSuperVcf(variants []string, db string, bsaseq bool, desc string, prgBlastFile string) (error, []string) {
	fmt.Println("Variants: ", variants)
	fmt.Println("Gene Description File: ", desc)

	// ========================================== Read description file ============================================= //
	fmt.Printf("Reading description file: %s ...............................\n", desc)

	descFile, err := os.Open(desc)
	if err != nil {
		return fmt.Errorf("failed to open description file %s: %w", desc, err), nil
	}
	defer descFile.Close()

	descScanner := bufio.NewScanner(descFile)

	descBuf := make([]byte, 0, 1024*1024)
	descScanner.Buffer(descBuf, 10*1024*1024)

	geneDesc := make(map[string]string)

	for descScanner.Scan() {
		fields := strings.Split(descScanner.Text(), "\t")

		if len(fields) < 2 {
			return fmt.Errorf("description file %s invalid format", desc), nil
		}

		geneDesc[fields[0]] = fields[1]
	}

	fmt.Printf("Loaded %d gene descriptions\n-------------------------------------------------\n\n", len(geneDesc))

	// ========================================== Read PRG blast file ============================================= //
	fmt.Printf("Reading PRG blast file: %s ...............................\n", prgBlastFile)

	prgFile, err := os.Open(prgBlastFile)
	if err != nil {
		return fmt.Errorf("failed to open PRG blast file %s: %w", prgBlastFile, err), nil
	}
	defer prgFile.Close()

	scanner := bufio.NewScanner(prgFile)
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	// ----------------------------------------- Read PRG header ------------------------------------------------- //

	if !scanner.Scan() {
		return fmt.Errorf("PRG blast file %s is empty", prgBlastFile), nil
	}

	prgHeaders := strings.Split(scanner.Text(), "\t")

	qseqIdx := -1
	percIdentIdx := -1
	lengthIdx := -1
	qlenIdx := -1
	evalIdx := -1

	for i, h := range prgHeaders {
		switch h {
		case "QseqID":
			qseqIdx = i
		case "PercIdent":
			percIdentIdx = i
		case "Length":
			lengthIdx = i
		case "Qlen":
			qlenIdx = i
		case "Eval":
			evalIdx = i
		}
	}

	missingCols := []string{}
	if qseqIdx == -1 {
		missingCols = append(missingCols, "QseqID")
	}
	if percIdentIdx == -1 {
		missingCols = append(missingCols, "PercIdent")
	}
	if lengthIdx == -1 {
		missingCols = append(missingCols, "Length")
	}
	if qlenIdx == -1 {
		missingCols = append(missingCols, "Qlen")
	}
	if evalIdx == -1 {
		missingCols = append(missingCols, "Eval")
	}
	if len(missingCols) > 0 {
		return fmt.Errorf("PRG blast file missing required columns: %s", strings.Join(missingCols, ", ")), nil
	}

	// ----------------------------------------- Read PRG rows ---------------------------------------------------- //

	type prgEntry struct {
		percIdent    string
		qlenMatchLen string
		eval         string
	}

	prgMap := make(map[string]prgEntry)

	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")

		if qseqIdx >= len(fields) || lengthIdx >= len(fields) || qlenIdx >= len(fields) {
			continue
		}

		qseqID := fields[qseqIdx]
		qlenMatchLen := fields[lengthIdx] + "/" + fields[qlenIdx]

		prgMap[qseqID] = prgEntry{
			percIdent:    fields[percIdentIdx],
			qlenMatchLen: qlenMatchLen,
			eval:         fields[evalIdx],
		}
	}

	fmt.Printf("Loaded %d PRG entries\n-------------------------------------------------\n\n", len(prgMap))

	var prgTsvFiles []string

	for _, variantFile := range variants {

		var snpEffTsv io.Reader
		var outputFileName string

		// -------------------------------------------- Input detection ------------------------------------------- //

		snpEffTsvFile, outputFileName, err := inputDetection(variantFile, bsaseq)
		snpEffTsv = open(snpEffTsvFile)
		// ----------------------------------------- Check variant file header ------------------------------------ //

		innerScanner := bufio.NewScanner(snpEffTsv)
		innerScanner.Buffer(buf, 10*1024*1024)

		if !innerScanner.Scan() {
			return fmt.Errorf("file %s is empty", variantFile), nil
		}

		headerLine := innerScanner.Text()
		headers := strings.Split(headerLine, "\t")

		geneIdx := -1
		transIdx := -1

		for i, h := range headers {
			if h == "SNPEFF_GENE_NAME" {
				geneIdx = i
			}
			if h == "SNPEFF_TRANSCRIPT_ID" {
				transIdx = i
			}
		}

		if geneIdx == -1 && transIdx == -1 {
			return fmt.Errorf("neither SNPEFF_GENE_NAME nor SNPEFF_TRANSCRIPT_ID found in %s", variantFile), nil
		}

		// ----------------------------- Read all rows and collect gene/transcript sets ---------------------------- //

		fmt.Println("Comparing gene/transcript names in variant file against PRG file ...")

		var allRows [][]string
		genesInVariants := make(map[string]struct{})
		transInVariants := make(map[string]struct{})

		readBar := progressbar.Default(-1, "reading rows")
		for innerScanner.Scan() {
			fields := strings.Split(innerScanner.Text(), "\t")
			allRows = append(allRows, fields)

			if geneIdx >= 0 && geneIdx < len(fields) {
				genesInVariants[fields[geneIdx]] = struct{}{}
			}
			if transIdx >= 0 && transIdx < len(fields) {
				transInVariants[fields[transIdx]] = struct{}{}
			}
			readBar.Add(1)
		}
		readBar.Finish()

		// ----------------------------- Check overlap and decide lookup column ------------------------------------ //

		geneOverlap := 0
		for g := range genesInVariants {
			if _, ok := prgMap[g]; ok {
				geneOverlap++
			}
		}

		transOverlap := 0
		for t := range transInVariants {
			if _, ok := prgMap[t]; ok {
				transOverlap++
			}
		}

		var useIdx int
		switch {
		case geneOverlap > 0:
			fmt.Println("PRG file uses GENE NAMES ... using SNPEFF_GENE_NAME ...")
			useIdx = geneIdx
		case transOverlap > 0:
			fmt.Println("PRG file uses TRANSCRIPT IDs ... using SNPEFF_TRANSCRIPT_ID ...")
			useIdx = transIdx
		default:
			fmt.Println("YOUR PRG FILE GENES ARE IN A DIFFERENT FORMAT TO THE FORMAT OF YOUR VCF GENES ...")
			fmt.Println("Re format your PRG FILE gene names to match those in GATK VCF file ...")
			return fmt.Errorf("no overlap between variant gene/transcript IDs and PRG QseqID"), nil
		}

		// ----------------------------------------------- Write output ------------------------------------------- //

		out, err := os.Create(outputFileName)
		if err != nil {
			return err, nil
		}
		defer out.Close()

		writer := bufio.NewWriter(out)

		writer.WriteString(headerLine + "\tPercIdent\tQLEN/MATCH_LEN\tEval\n")

		fmt.Printf("\nWriting output to %s\n\n", outputFileName)

		writeBar := progressbar.Default(int64(len(allRows)), "writing output")
		for _, fields := range allRows {
			percIdent := "0"
			qlenMatchLen := "0"
			eval := "0"

			if useIdx >= 0 && useIdx < len(fields) {
				if entry, ok := prgMap[fields[useIdx]]; ok {
					percIdent = entry.percIdent
					qlenMatchLen = entry.qlenMatchLen
					eval = entry.eval
				}
			}

			fields = append(fields, percIdent, qlenMatchLen, eval)
			writer.WriteString(strings.Join(fields, "\t") + "\n")
			writeBar.Add(1)
		}
		writeBar.Finish()

		writer.Flush()

		prgTsvFiles = append(prgTsvFiles, outputFileName)

		fmt.Printf("\nCreated: %s\n", outputFileName)
	}

	return nil, prgTsvFiles



	return nil, descTsvFiles
}
