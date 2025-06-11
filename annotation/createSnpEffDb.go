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
	fmt.Println(cmdStr1)
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

	var seqs string
	var pro string
	var cd string

	if strings.HasSuffix(ref, ".gz") {
		seqs = filepath.Join(dataDir, "sequences.fa.gz")
	} else {
		seqs = filepath.Join(dataDir, "sequences.fa")
	}

	if strings.HasSuffix(prot, ".gz") {
		pro = filepath.Join(dataDir, "protein.fa.gz")
	} else {
		pro = filepath.Join(dataDir, "protein.fa")
	}

	if strings.HasSuffix(cds, ".gz") {
		cd = filepath.Join(dataDir, "cds.fa.gz")
	} else {
		cd = filepath.Join(dataDir, "cds.fa")
	}

	fmt.Printf("Copying %s to %s ...\n\n", ref, seqs)
	crErr := CopyFile(ref, seqs)
	if crErr != nil {
		fmt.Printf("failed to copy  %s to %s", ref, seqs)
		return crErr
	}

	fmt.Printf("Copying %s to %s ...\n\n", prot, pro)

	cpErr := CopyFile(prot, pro)
	if cpErr != nil {
		fmt.Printf("failed to copy %s to %s", prot, pro)
		return cpErr
	}

	fmt.Printf("Copying %s to %s ...\n\n", cds, cd)

	ccErr := CopyFile(cds, cd)
	if ccErr != nil {
		fmt.Printf("failed to copy  %s: %v", cds, ccErr)
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
