package annotation

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

func check_err(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
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

	// ------------------------------------ Check if db is not already installed ------------------------------------ //

	//------------------------------------------ Editing config file ------------------------------------------------ //
	fmt.Printf("Editing snpEff config file ...\n\n")
	firstLine := fmt.Sprintf("# %s genome, version %s%s", species, species, version)
	secondLine := fmt.Sprintf("%s%s.genome : %s", species, version, species)

	fmt.Printf("Adding %s to snpEff config file ...\n\n", firstLine)
	fmt.Printf("Adding %s to snpEff config file ...\n\n", secondLine)

	f, err := os.OpenFile(configPath, os.O_APPEND|os.O_WRONLY, 0644)
	check_err(err)
	defer f.Close()
	_, err = f.WriteString(firstLine + "\n" + secondLine + "\n")
	check_err(err)

	//---------------------------------- Copying annotation files to snpEff data dir -------------------------------- //

	fmt.Printf("Creating directories ...\n\n")

	dataDir := filepath.Join(snpEffDir, "data")
	fmt.Printf(" mkdir -p %s", dataDir)
	err = os.MkdirAll(dataDir, 0755)
	check_err(err)

	db := fmt.Sprintf("%s%s", species, version)
	dbDir := filepath.Join(dataDir, db)
	fmt.Printf(" mkdir -p %s", dbDir)
	err = os.MkdirAll(dbDir, 0755)
	check_err(err)

	fmt.Printf("Copying reference, cds and protein files ...\n\n")

	var seqs string
	var pro string
	var cd string

	if strings.HasSuffix(ref, ".gz") {
		seqs = filepath.Join(dbDir, "sequences.fa.gz")
	} else {
		seqs = filepath.Join(dbDir, "sequences.fa")
	}

	if strings.HasSuffix(prot, ".gz") {
		pro = filepath.Join(dbDir, "protein.fa.gz")
	} else {
		pro = filepath.Join(dbDir, "protein.fa")
	}

	if strings.HasSuffix(cds, ".gz") {
		cd = filepath.Join(dbDir, "cds.fa.gz")
	} else {
		cd = filepath.Join(dbDir, "cds.fa")
	}

	fmt.Printf("Copying %s to %s ...\n\n", ref, seqs)
	crErr := CopyFile(ref, seqs)
	check_err(crErr)

	fmt.Printf("Copying %s to %s ...\n\n", prot, pro)
	err = CopyFile(prot, pro)
	check_err(err)

	fmt.Printf("Copying %s to %s ...\n\n", cds, cd)
	ccErr := CopyFile(cds, cd)
	check_err(ccErr)

	var gffreadCmd string
	if strings.HasSuffix(gff, ".gz") {
		// For gzipped files, pipe through gunzip
		gffreadCmd = fmt.Sprintf(`gunzip -c %s | gffread - -T -o %s`, gff, filepath.Join(dbDir, "genes.gtf"))
	} else {
		// For uncompressed files, use direct gffread
		gffreadCmd = fmt.Sprintf(`gffread %s -T -o %s`, gff, filepath.Join(dbDir, "genes.gtf"))
	}

	fmt.Printf("Running: %s\n", gffreadCmd)
	err3 := utils.RunBashCmdVerbose(gffreadCmd)
	if err3 != nil {
		return fmt.Errorf("failed to convert GFF to GTF: %w", err3)
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
		return fmt.Errorf("couldn't open source file %s: %w", src, sErr)
	}
	defer sourceFile.Close()

	dstFile, dErr := os.Create(dst)
	if dErr != nil {
		return fmt.Errorf("couldn't create destination file %s: %w", dst, dErr)
	}
	defer dstFile.Close()

	_, err := io.Copy(dstFile, sourceFile)
	if err != nil {
		return fmt.Errorf("failed to copy file contents: %w", err)
	}

	return nil
}
