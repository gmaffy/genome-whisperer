package annotation

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"regexp"

	"github.com/gmaffy/genome-whisperer/utils"

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
	// ================================== Reading description file ===================================================//
	descFile, err := os.Open(desc)
	if err != nil {
		return fmt.Errorf("failed to open description file %s: %w", desc, err), []string{}
	}
	defer descFile.Close()

	scanner := bufio.NewScanner(descFile)
	geneDesc := make(map[string]string)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) != 2 {
			return fmt.Errorf("description file %s has invalid format. Expected 2 columns, got %d", desc, len(fields)), []string{}
		}
		gene, description := fields[0], fields[1]
		geneDesc[gene] = description
	}

	//================================ Adding descriptions to tsv files =============================================//

	fmt.Printf("Adding descriptions to variant files: %s\n", variants)
	var descTsvFiles []string
	for _, variantFile := range variants {
		fmt.Printf("Adding descriptions to variant file: %s\n", variantFile)
		var snpEffTsv io.Reader
		var outputFileName string
		if strings.HasSuffix(variantFile, ".tsv") || strings.HasSuffix(variantFile, ".txt") {
			var snpEffTsvFile *os.File
			snpEffTsvFile, err = os.Open(variantFile)
			if err != nil {
				return fmt.Errorf("failed to open input file %s: %w", variantFile, err), []string{}
			}
			defer snpEffTsvFile.Close()
			snpEffTsv = snpEffTsvFile
			outputFileName = strings.Replace(strings.Replace(variantFile, ".txt", "_DESC.tsv", 1), ".tsv", "_DESC.tsv", 1)

		} else if strings.HasSuffix(variantFile, ".vcf") || strings.HasSuffix(variantFile, "vcf.gz") {
			fmt.Println("Converting vcf to table")
			var vtCmdStr string
			var variantTsvFile string
			if strings.HasSuffix(variantFile, ".vcf") {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf") + ".tsv"
			} else {
				variantTsvFile = strings.TrimSuffix(variantFile, ".vcf.gz") + ".tsv"
			}
			if bsaseq {
				vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`, variantFile, variantTsvFile)
			} else {
				vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`, variantFile, variantTsvFile)
			}
			fmt.Println(vtCmdStr)
			terr := utils.RunBashCmdVerbose(vtCmdStr)
			if terr != nil {
				return terr, []string{}
			}

			sErr, splitSnpEffTsv := splitEffColumns(variantTsvFile)
			if sErr != nil {
				return sErr, []string{}
			}
			var snpEffTsvFile *os.File
			snpEffTsvFile, err = os.Open(splitSnpEffTsv)
			if err != nil {
				return fmt.Errorf("failed to open input file %s: %w", splitSnpEffTsv, err), []string{}
			}
			defer snpEffTsvFile.Close()
			snpEffTsv = snpEffTsvFile
			outputFileName = strings.Replace(strings.Replace(splitSnpEffTsv, ".txt", "_DESC.tsv", 1), ".tsv", "_DESC.tsv", 1)

		} else {
			return fmt.Errorf("variant file %s has invalid format. Expected .vcf or .tsv", variantFile), []string{}
		}

		outputFile, err := os.Create(outputFileName)
		if err != nil {
			return fmt.Errorf("failed to create output file %s: %w", outputFileName, err), []string{}
		}
		defer outputFile.Close()

		writer := bufio.NewWriter(outputFile)

		innerScanner := bufio.NewScanner(snpEffTsv)

		if !innerScanner.Scan() {
			return fmt.Errorf("input file %s is empty or has no header", variantFile), []string{}
		}
		headerLine := innerScanner.Text()
		originalHeaders := strings.Split(headerLine, "\t")

		transcriptIdx := -1
		for i, h := range originalHeaders {
			if h == "SNPEFF_TRANSCRIPT_ID" {
				transcriptIdx = i
			}
		}

		if transcriptIdx == -1 {
			return fmt.Errorf("SNPEFF_TRANSCRIPT_ID column not found in header: %s", headerLine), []string{}
		}

		newHeader := strings.Join(originalHeaders, "\t") + "\tDESCRIPTION"
		for innerScanner.Scan() {
			line := innerScanner.Text()
			fields := strings.Split(line, "\t")
			if fields[0] == "CHROM" {

				_, err2 := writer.WriteString(newHeader + "\n")
				if err2 != nil {
					return err2, []string{}
				}
				continue
			}
			transcriptID := fields[transcriptIdx]
			if description, ok := geneDesc[transcriptID]; ok {
				fields = append(fields, description)
			} else {
				fields = append(fields, "")
			}

			_, err2 := writer.WriteString(strings.Join(fields, "\t") + "\n")
			if err2 != nil {
				return err2, []string{}
			}
		}

		err2 := writer.Flush()
		if err2 != nil {
			return err2, []string{}
		}
		descTsvFiles = append(descTsvFiles, outputFileName)
	}
	return nil, descTsvFiles
}

//func AddPrg(vcfs []string, prg string) error {
//	fmt.Printf("Adding program to vcfs %s\n", vcfs)
//	for _, vcf := range vcfs {
//		fmt.Printf("Adding program to vcf: %s\n", vcf)
//		prgCmdStr := fmt.Sprintf(`bcftools annotate -h "##PRG=%s" %s > %s`, prg, vcf, vcf)
//		fmt.Println(prgCmdStr)
//	}
//}

func CreateSuperVcf(vcfs []string, db string, bsaseq bool, desc string, prg string) (error, []string) {
	fmt.Printf("Running snpEff for vcfs %s\n", vcfs)
	err, snpEffTsvFiles, snpEffvcfFiles := RunSnpEff(vcfs, db, bsaseq)
	if err != nil {
		return err, []string{}
	}
	fmt.Printf("SnpEff done for vcfs %s complete: tsvFiles: %s, annotated vcf files: %s \n", vcfs, snpEffTsvFiles, snpEffvcfFiles)

	err, descTsvFiles := AddDescriptions(snpEffTsvFiles, desc, bsaseq)
	if err != nil {
		return err, []string{}
	}

	return nil, descTsvFiles
}
