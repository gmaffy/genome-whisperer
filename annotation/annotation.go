package annotation

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"regexp"

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
			// grep returns 1 if no lines are selected
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

func RunSnpEff(vcfs []string, db string, bsaseq bool) error {
	// --------------------------------------------------- Get Config file ---------------------------------------- //
	snEffPath, err := exec.LookPath("snpEff")
	if err != nil {
		return fmt.Errorf("snpEff not found: %w", err)
	}

	scriptsDir := filepath.Dir(snEffPath)
	snpEffDir := filepath.Dir(scriptsDir)
	configPath := filepath.Join(snpEffDir, "snpEff.config")

	_, rErr := os.Stat(configPath)
	if rErr != nil {
		fmt.Printf("Tried to find snpEff config file at: %s. It does not exist\n\n", configPath)
		return rErr
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
		return dbErr
	}
	fmt.Printf("DATABASE FOUND: %s\n", db)
	// ------------------------------------------------ Run SnpEff -------------------------------------------------- //
	for _, vcf := range vcfs {
		var snpEffVcf string
		var snpEffTxt string
		if strings.HasSuffix(vcf, ".vcf.gz") {
			snpEffVcf = strings.TrimSuffix(vcf, ".vcf.gz") + ".snpEff.vcf"
			snpEffTxt = strings.TrimSuffix(vcf, ".vcf.gz") + ".snpEff.txt"
		} else if strings.HasSuffix(vcf, ".vcf") {
			snpEffVcf = strings.TrimSuffix(vcf, ".vcf") + ".snpEff.vcf"
			snpEffTxt = strings.TrimSuffix(vcf, ".vcf") + ".snpEff.txt"
		} else {
			fmt.Println("vcf file must be in vcf or vcf.gz format")
			return fmt.Errorf("vcf file must be in vcf or vcf.gz format")
		}

		snpEffCmdStr := fmt.Sprintf(`snpEff -c %s -v -o gatk %s %s > %s`, configPath, db, vcf, snpEffVcf)
		fmt.Println(snpEffCmdStr)
		err := utils.RunBashCmdVerbose(snpEffCmdStr)
		if err != nil {
			return err
		}

		fmt.Println("Converting vcf to table")
		var vtCmdStr string
		if bsaseq {
			vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -GF AD -GF DP -GF GQ -F EFF -O %s`, snpEffVcf, snpEffTxt)
		} else {
			vtCmdStr = fmt.Sprintf(`gatk VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -F QUAL -F TYPE -GF GT -F EFF -O %s`, snpEffVcf, snpEffTxt)
		}
		fmt.Println(vtCmdStr)
		terr := utils.RunBashCmdVerbose(vtCmdStr)
		if terr != nil {
			return terr
		}

		sErr := splitEffColumns(snpEffTxt)
		if sErr != nil {
			return sErr
		}

	}

	return nil
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

func splitEffColumns(effFile string) error {
	fmt.Printf("Splitting EFF column in file: %s \n\n ", effFile)
	inputFile, err := os.Open(effFile)
	if err != nil {
		return fmt.Errorf("failed to open input file %s: %w", effFile, err)
	}
	defer inputFile.Close()
	scanner := bufio.NewScanner(inputFile)

	if !scanner.Scan() {
		return fmt.Errorf("input file %s is empty or has no header", effFile)
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
		return fmt.Errorf("EFF column not found in header: %s", headerLine)
	}

	snpEffNewHeaders := []string{
		"SNPEFF_EFFECT", "SNPEff_IMPACT", "SNPEFF_FUNCTIONAL_CLASS",
		"SNPEFF_CODON_CHANGE", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME",
		"SNPEFF_TRANSCRIPT_BIOTYPE", "SNPEFF_GENE_CODING",
		"SNPEFF_TRANSCRIPT_ID", "ERRORS", "WARNINGS",
	}

	newHeader := strings.Join(otherHeaders, "\t") + "\t" + strings.Join(snpEffNewHeaders, "\t")

	outputFileName := strings.Replace(effFile, ".txt", "_EFF.tsv", 1)
	if !strings.HasSuffix(effFile, ".txt") {
		outputFileName = effFile + "_EFF.txt"
	}

	fmt.Printf("Writing to file: %s ...\n", outputFileName)
	outputFile, err := os.Create(outputFileName)
	if err != nil {
		return fmt.Errorf("failed to create output file %s: %w", outputFileName, err)
	}
	defer outputFile.Close()

	writer := bufio.NewWriter(outputFile)
	_, err = writer.WriteString(newHeader + "\n")
	if err != nil {
		return fmt.Errorf("failed to write header to output file: %w", err)
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
			return fmt.Errorf("failed to write line %d to output file: %w", lineNum, err)
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading input file: %w", err)
	}

	return writer.Flush()
}

func CreateCustomDb(ref, prot, cds, species, gff, version string) error {
	// --------------------------------------------------- Get Config file ---------------------------------------- //
	snEffPath, err := exec.LookPath("snpEff")
	if err != nil {
		return fmt.Errorf("snpEff not found: %w", err)
	}

	scriptsDir := filepath.Dir(snEffPath)
	snpEffDir := filepath.Dir(scriptsDir)
	configPath := filepath.Join(snpEffDir, "snpEff.config")

	_, rErr := os.Stat(configPath)
	if rErr != nil {
		fmt.Printf("Tried to find snpEff config file at: %s. It does not exist\n\n", configPath)
		return rErr
	}

	fmt.Printf("Editing snpEff config file ...\n\n")

	cmdStr1 := fmt.Sprintf(`sed -i "164 a # %s genome, version %s%s"  %s`, species, species, version, configPath)
	fmt.Printf(cmdStr1)
	err1 := utils.RunBashCmdVerbose(cmdStr1)
	if err1 != nil {
		return err1
	}

	cmdStr2 := fmt.Sprintf(`sed -i "165 a %s%s.genome : %s" %s`, species, version, species, configPath)
	fmt.Printf(cmdStr2)
	err2 := utils.RunBashCmdVerbose(cmdStr2)
	if err2 != nil {
		return err2
	}

	fmt.Printf("Creating directories ...\n\n")
	dataDir := filepath.Join(snpEffDir, "data")
	db := fmt.Sprintf("%s%s", species, version)
	dbDir := filepath.Join(dataDir, db)
	dErr := os.MkdirAll(dataDir, 0755)
	if dErr != nil {
		return dErr
	}
	dErr1 := os.MkdirAll(dbDir, 0755)
	if dErr1 != nil {
		return dErr1
	}

	fmt.Printf(" mkdir %s", dbDir)
	fmt.Printf("Copying reference, cds and protein files ...\n\n")

	fmt.Printf("Copying %s to %s ...\n\n", ref, filepath.Join(dataDir, "sequences.fa"))
	crErr := CopyFile(ref, filepath.Join(dataDir, "sequences.fa"))
	if crErr != nil {
		fmt.Printf("failed to copy  %s to %s", ref, filepath.Join(dataDir, "sequences.fa"))
		return crErr
	}

	fmt.Printf("Copying %s to %s ...\n\n", prot, filepath.Join(dataDir, "protein.fa"))

	cpErr := CopyFile(prot, filepath.Join(dataDir, "protein.fa"))
	if cpErr != nil {
		fmt.Printf("failed to copy %s to %s", prot, filepath.Join(dataDir, "protein.fa"))
		return cpErr
	}

	fmt.Printf("Copying %s to %s ...\n\n", cds, filepath.Join(dataDir, "cds.fa"))

	ccErr := CopyFile(cds, filepath.Join(dataDir, "cds.fa"))
	if ccErr != nil {
		fmt.Printf("failed to copy  sequences.fa: %v", ccErr)
		return ccErr
	}

	cmdStr3 := fmt.Sprintf(`gffread %s -T -o %s`, gff, filepath.Join(dataDir, "genes.gtf"))
	fmt.Printf(cmdStr3)
	err3 := utils.RunBashCmdVerbose(cmdStr3)
	if err3 != nil {
		return err3
	}

	fmt.Printf("Building database ...")

	cmdStr4 := fmt.Sprintf(`snpEff build -gtf22 -v %s`, db)
	fmt.Printf(cmdStr4)
	err4 := utils.RunBashCmdVerbose(cmdStr4)
	if err4 != nil {
		return err4
	}

	return nil
}

func CopyFile(src, dst string) error {

	sourceFile, sErr := os.Open(src)
	if sErr != nil {
		return fmt.Errorf("couldn't open reference file %s: %w", src, sErr)
	}
	defer sourceFile.Close()

	dstFile, dstErr := os.Open(dst)
	if dstFile != nil {
		return fmt.Errorf("couldn't open sequences.fa: %w", dstErr)
	}
	defer dstFile.Close()

	_, err := io.Copy(dstFile, sourceFile)
	if err != nil {
		return fmt.Errorf("failed to copy file contents: %w", err)
	}

	return nil

}
