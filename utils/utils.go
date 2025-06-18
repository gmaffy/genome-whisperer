package utils

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
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

type LogEntry struct {
	Timestamp  string
	Tool       string
	Program    string
	Sample     string
	Chromosome string
	Status     string
	Cmd        string
}

func CheckDeps() error {
	deps := []string{"gatk", "samtools", "bwa", "java", "snpEff", "gffread", "masurca", "MAC2.0", "megahit", "seqtk", "bowtie2", "bedtools"}

	for _, dep := range deps {
		if _, err := exec.LookPath(dep); err != nil {
			fmt.Printf("%s not found!\n\n", dep)
			return fmt.Errorf("%s not found: %w", dep, err)
		}
		fmt.Printf("%s OK\n", dep)
	}

	for _, prog := range deps {
		path, _ := exec.LookPath(prog)
		fmt.Printf("Using %s at %s\n", prog, path)
	}

	return nil
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
		fmt.Println("CMD error:", err)
		return err
	}
	return nil
}

//type LogEntry struct {
//	Timestamp  string
//	Chromosome string
//	Program    string
//	Bam        string
//	Status     string
//	Cmd        string
//}
//
//type ChromosomeSamplePair struct {
//	Chromosome string
//	Sample     string
//}

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

func ParseLogFile(logFilePath string) []LogEntry {
	var data []LogEntry
	file, err := os.Open(logFilePath)
	if err != nil {
		fmt.Printf("Log file '%s' not found, starting fresh or assuming no previous runs.\n", logFilePath)
		return data
	}
	defer file.Close()
	fmt.Println("Parsing log file ...")

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		var entry map[string]interface{}
		lineBytes := scanner.Bytes()
		if len(lineBytes) == 0 { // Skip empty lines
			continue
		}
		if err := json.Unmarshal(lineBytes, &entry); err != nil {
			fmt.Printf("Warning: Skipping malformed log line: %v - Line: %s\n", err, string(lineBytes))
			continue // skip malformed line
		}

		timestamp, timeOk := entry["time"]
		tool, toolOk := entry["msg"]
		chromVal, chromOk := entry["CHROM"]
		statusVal, statusOk := entry["STATUS"]
		sampleVal, sampleOk := entry["SAMPLE"]
		programVal, programOk := entry["PROGRAM"]
		if chromOk && statusOk && sampleOk && programOk && timeOk && toolOk {
			r := LogEntry{
				Timestamp:  timestamp.(string),
				Tool:       tool.(string),
				Program:    programVal.(string),
				Sample:     sampleVal.(string),
				Chromosome: chromVal.(string),
				Status:     statusVal.(string),
			}
			data = append(data, r)
		}

	}

	return data
}

func StageHasCompleted(logEntries []LogEntry, prog string, sample string, chrom string) bool {
	for _, entry := range logEntries {
		if entry.Program == prog && entry.Sample == sample && entry.Chromosome == chrom && entry.Status == "COMPLETED" {
			return true
		}
	}
	return false
}

//
//func GetHapFinished(logEntries []LogEntry) []ChromosomeSamplePair {
//	var hapFinished []ChromosomeSamplePair
//	for _, entry := range logEntries {
//		if entry.Program == "HaplotypeCaller" && entry.Status == "FINISHED" {
//			pair := ChromosomeSamplePair{Chromosome: entry.Chromosome, Sample: entry.Bam}
//			hapFinished = append(hapFinished, pair)
//		}
//	}
//	return hapFinished
//}
//
//func GetCompletedStages(logEntries []LogEntry) map[string]map[string]bool {
//	cs := make(map[string]map[string]bool)
//	for _, entry := range logEntries {
//
//		if entry.Program == "GenomicsDBImport" && entry.Status == "FINISHED" {
//			if _, ok := cs["GenomicsDBImport"]; !ok {
//				cs["GenomicsDBImport"] = make(map[string]bool)
//			}
//			cs["GenomicsDBImport"][entry.Chromosome] = true
//		} else if entry.Program == "GenotypeGVCFs" && entry.Status == "FINISHED" {
//			if _, ok := cs["GenotypeGVCFs"]; !ok {
//				cs["GenotypeGVCFs"] = make(map[string]bool)
//			}
//			cs["GenotypeGVCFs"][entry.Chromosome] = true
//		} else if entry.Program == "SELECT_SNPS" && entry.Status == "FINISHED" {
//			if _, ok := cs["SELECT_SNPS"]; !ok {
//				cs["SELECT_SNPS"] = make(map[string]bool)
//			}
//			cs["SELECT_SNPS"][entry.Chromosome] = true
//
//		} else if entry.Program == "SELECT_INDELS" && entry.Status == "FINISHED" {
//			if _, ok := cs["SELECT_INDELS"]; !ok {
//				cs["SELECT_INDELS"] = make(map[string]bool)
//			}
//			cs["SELECT_INDELS"][entry.Chromosome] = true
//		} else if entry.Program == "HardFilteringSNPS" && entry.Status == "FINISHED" {
//			if _, ok := cs["HardFilteringSNPS"]; !ok {
//				cs["HardFilteringSNPS"] = make(map[string]bool)
//			}
//			cs["HardFilteringSNPS"][entry.Chromosome] = true
//		} else if entry.Program == "HardFilteringINDELS" && entry.Status == "FINISHED" {
//			if _, ok := cs["HardFilteringINDELS"]; !ok {
//				cs["HardFilteringINDELS"] = make(map[string]bool)
//			}
//			cs["HardFilteringINDELS"][entry.Chromosome] = true
//		} else if entry.Program == "MergeVcfs" && entry.Status == "FINISHED" {
//			if _, ok := cs["MergeVcfs"]; !ok {
//				cs["MergeVcfs"] = make(map[string]bool)
//			}
//			cs["MergeVcfs"][entry.Chromosome] = true
//		}
//	}
//	return cs
//}

//func hasCompletedStatus(logFile string, targetChrom int) (bool, error) {
//	file, err := os.Open(logFile)
//	if err != nil {
//		return false, err
//	}
//	defer file.Close()
//
//	scanner := bufio.NewScanner(file)
//	for scanner.Scan() {
//		var entry map[string]interface{}
//		if err := json.Unmarshal(scanner.Bytes(), &entry); err != nil {
//			continue // skip malformed line
//		}
//
//		chromVal, chromOk := entry["CHROM"]
//		statusVal, statusOk := entry["STATUS"]
//
//		// Check CHROM and STATUS
//		if chromOk && statusOk {
//			if intVal, ok := toInt(chromVal); ok && intVal == targetChrom {
//				if statusStr, ok := statusVal.(string); ok && statusStr == "COMPLETED" {
//					return true, nil
//				}
//			}
//		}
//	}
//
//	return false, scanner.Err()
//}
