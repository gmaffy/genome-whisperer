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

func checkErr(err error) {
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
	fmt.Printf("Checking if database %s is installed ...\n\n", species)
	db := fmt.Sprintf("%s%s", species, version)
	dbErr := checkSnpEffDB(db)
	if dbErr == nil {
		fmt.Printf("Database %s already. Database creation stopped.\n\n", db)
		return dbErr
	}

	//------------------------------------------ Editing config file ------------------------------------------------ //
	fmt.Printf("Editing snpEff config file ...\n\n")
	firstLine := fmt.Sprintf("# %s genome, version %s%s", species, species, version)
	secondLine := fmt.Sprintf("%s%s.genome : %s", species, version, species)

	fmt.Printf("Adding %s to snpEff config file ...\n\n", firstLine)
	fmt.Printf("Adding %s to snpEff config file ...\n\n", secondLine)

	f, err := os.OpenFile(configPath, os.O_APPEND|os.O_WRONLY, 0644)
	checkErr(err)
	defer f.Close()
	_, err = f.WriteString(firstLine + "\n" + secondLine + "\n")
	checkErr(err)

	//---------------------------------- Copying annotation files to snpEff data dir -------------------------------- //

	fmt.Printf("Creating directories ...\n\n")

	dataDir := filepath.Join(snpEffDir, "data")
	fmt.Printf(" mkdir -p %s", dataDir)
	err = os.MkdirAll(dataDir, 0755)
	checkErr(err)

	dbDir := filepath.Join(dataDir, db)
	fmt.Printf(" mkdir -p %s", dbDir)
	err = os.MkdirAll(dbDir, 0755)
	checkErr(err)

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
	checkErr(crErr)

	fmt.Printf("Copying %s to %s ...\n\n", prot, pro)
	err = CopyFile(prot, pro)
	checkErr(err)

	fmt.Printf("Copying %s to %s ...\n\n", cds, cd)
	ccErr := CopyFile(cds, cd)
	checkErr(ccErr)

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

func CreateCustomDbFromConfig(configFile, species, version string) error {
	fmt.Println("Reading config file ...")
	cfg, err1 := utils.ReadConfig(configFile)
	if err1 != nil {
		fmt.Printf("Error reading config: %v\n", err1)
		return err1
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Proteins:", cfg.Proteins)
	fmt.Println("CDS:", cfg.CDS)
	fmt.Println("GFF:", cfg.GFF)

	refFile := cfg.Reference
	protein := cfg.Proteins
	cds := cfg.CDS
	gff := cfg.GFF
	_, err := os.Stat(refFile)
	if err != nil {
		fmt.Printf("Reference file: %s is not a valid file path", refFile)
		return err
	}

	_, err = os.Stat(protein)
	if err != nil {
		fmt.Printf("Protein file: %s is not a valid file path", protein)
		return err
	}

	_, err = os.Stat(cds)
	if err != nil {
		fmt.Printf("CDS file: %s is not a valid file path", cds)
		return err
	}

	_, err = os.Stat(gff)
	if err != nil {
		fmt.Printf("GFF file: %s is not a valid file path", gff)
		return err
	}

	if species == "" {
		fmt.Println("Please provide species name")
		return err
	}

	if version == "" {
		fmt.Println("Please provide version")
		return err
	}

	err = CreateCustomDb(refFile, protein, cds, species, gff, version)
	if err != nil {
		return err
	}

	return nil
}
