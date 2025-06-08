package utils

import (
	"bufio"
	"os"
	"os/exec"
	"strings"
)

type Config struct {
	Reference   string
	GFF         string
	Proteins    string
	CDS         string
	Species     string
	OutputDir   string
	BaseName    string
	Bams        []string
	ReadPairs   [][]string
	VCF         string
	Version     string
	VCFs        []string
	SelectChrom string
	SelectStart string
	SelectStop  string
	SelectVCF   string
	KnownSites  []string
	SnpEff      string

	Java8      string
	Threads    string
	InputDir   string
	DataType   string
	GVCFsDir   string
	CallerName string
	GVCFs      []string
}

func ReadConfig(configPath string) (Config, error) {
	configFile, err := os.Open(configPath)
	if err != nil {
		return Config{}, err
	}
	defer configFile.Close()
	var cfg Config

	scanner := bufio.NewScanner(configFile)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if line == "" {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])

		switch key {
		case "Reference":
			cfg.Reference = value
		case "gff":
			cfg.GFF = value
		case "proteins":
			cfg.Proteins = value
		case "cds":
			cfg.CDS = value
		case "Species":
			cfg.Species = value
		case "OutputDir":
			cfg.OutputDir = value
		case "BaseName":
			cfg.BaseName = value
		case "bam":
			cfg.Bams = append(cfg.Bams, value)

		case "ReadPair":
			pairs := strings.Fields(value)
			cfg.ReadPairs = append(cfg.ReadPairs, pairs)
		case "VCF":
			cfg.VCF = value
			cfg.VCFs = append(cfg.VCFs, value)
		case "gvcf":
			cfg.GVCFs = append(cfg.GVCFs, value)
		case "select_chrom":
			cfg.SelectChrom = value
		case "select_start":
			cfg.SelectStart = value
		case "select_stop":
			cfg.SelectStop = value
		case "select_vcf":
			cfg.SelectVCF = value
		case "known-sites":
			cfg.KnownSites = append(cfg.KnownSites, value)
		case "snpEff":
			cfg.SnpEff = value

		case "java_8":
			cfg.Java8 = value
		case "threads":
			cfg.Threads = value
		case "InputDir":
			cfg.InputDir = value
		case "DATA_TYPE":
			cfg.DataType = value
		case "GVCFS_Dir":
			cfg.GVCFsDir = value
		case "Caller_name":
			cfg.CallerName = value
		}
	}

	if err := scanner.Err(); err != nil {
		return cfg, err
	}

	return cfg, nil

}

func RunBashCmdVerbose(cmdStr string) error {
	cmd := exec.Command("bash", "-c", cmdStr)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	if err != nil {
		return err
	}
	return nil
}
