package alignment

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"os/exec"
	"runtime"
	"strings"
	"sync"
)

type Config struct {
	Reference   string
	GFF         string
	Proteins    string
	CDS         string
	Species     string
	OutputDir   string
	BaseName    string
	RawBams     []string
	RGMDBams    []string
	ReadPairs   [][]string
	VCF         string
	Version     string
	VCFs        []string
	SelectChrom string
	SelectStart string
	SelectStop  string
	SelectVCF   string
	BQSRBams    []string
	DbSNPs      []string
	SnpEff      string

	Java8      string
	Threads    string
	InputDir   string
	DataType   string
	GVCFsDir   string
	CallerName string
	GVCFs      []string
}

func AlignShortReadsMem(referencePath string, forwardPath string, reversePath string, sampleName string, libName string, outputDir string, threads int) {
	fmt.Println("Reading ...")
	fmt.Printf("ReferencePath: %s\nFwd: %s\nRev: %s\nSample Name: %s\nlib Name: %s\nOuDir: %s\n", referencePath, forwardPath, reversePath, sampleName, libName, outputDir)

	lineDir := fmt.Sprintf("%s/%s", outputDir, sampleName)
	bErr := os.MkdirAll(lineDir, 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
	}
	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	sortedBam := fmt.Sprintf("%s/%s.sorted.bam", lineDir, sampleName)
	rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sampleName)
	rgmdMetrics := fmt.Sprintf("%s/%s.RGMD.metrics.txt", lineDir, sampleName)
	rgmdIndex := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sampleName)

	// Construct the one-liner shell command
	cmdStr := fmt.Sprintf(`bwa mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
	fmt.Println(cmdStr)
	cmd := exec.Command("bash", "-c", cmdStr)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	if err != nil {
		return
	}
	fmt.Printf("Marking duplicates ....")
	mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -I %s -O %s -M %s`, sortedBam, rgmdBam, rgmdMetrics)
	mDupCmd := exec.Command("bash", "-c", mDupCmdStr)
	mDupCmd.Stdout = os.Stdout
	mDupCmd.Stderr = os.Stderr

	mErr := mDupCmd.Run()
	if mErr != nil {
		return
	}

	fmt.Printf("Index Bam ....")
	indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
	indexCmd := exec.Command("bash", "-c", indexCmdStr)
	indexCmd.Stdout = os.Stdout
	indexCmd.Stderr = os.Stderr

	iErr := indexCmd.Run()
	if iErr != nil {
		return
	}

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
		case "rawBam":
			cfg.RawBams = append(cfg.RawBams, value)
		case "rgmdBam":
			cfg.RGMDBams = append(cfg.RGMDBams, value)
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
		case "bqsrBam":
			cfg.BQSRBams = append(cfg.BQSRBams, value)
		case "dbSNP":
			cfg.DbSNPs = append(cfg.DbSNPs, value)
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

func AlignShortReadsConfig(configPath string, threadsPerSample int) {
	fmt.Println("Reading config file ...")
	cfg, err := ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Species:", cfg.Species)
	fmt.Println("Read Pairs:", cfg.ReadPairs)
	fmt.Println("Threads", cfg.Threads)
	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	//num, nErr := strconv.Atoi(cfg.Threads)
	//if nErr != nil {
	//	fmt.Printf("Threads should be an integer: %v\n", nErr)
	//	return
	//}
	//threadsPerSample := threads
	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threadsPerSample = totalCores
	}

	fmt.Printf("Running up to %d jobs in parallel with %d threads each\n", maxParallelJobs, threadsPerSample)

	// Check paths
	for _, pair := range cfg.ReadPairs {
		if len(pair) < 4 {
			fmt.Printf("This read pair is wrongly formated %s\n", pair)
			fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			continue
		}
		ref := cfg.Reference
		out := cfg.OutputDir
		fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]
		_, refErr := os.Stat(ref)
		_, fwdErr := os.Stat(fwd)
		_, revErr := os.Stat(rev)
		outInfo, outErr := os.Stat(out)
		if refErr != nil {
			fmt.Printf("Reference genome path: %s, is not valid\n", ref)
			return
		}

		if fwdErr != nil {
			fmt.Printf("Forward reads path %s, is not valid\n", fwd)
			return
		}

		if revErr != nil {
			fmt.Printf("Reverse reads path %s, is not valid\n", rev)
			return
		}

		if outErr != nil {
			fmt.Printf("Output directory: %s is not a valid path\n", out)
			return
		}
		if !outInfo.IsDir() {
			fmt.Printf("Output Directory %s file path is not a directory", out)
			return
		}
		if sn == "" {
			fmt.Println("Please provide sample name is flag -s ")
			return
		}
		if lb == "" {
			fmt.Println("Please provide library name is flag -l ")
			return
		}
	}

	// ============================================== Run Alignments ================================================ //
	var wg sync.WaitGroup
	sem := make(chan struct{}, maxParallelJobs)
	for _, pair := range cfg.ReadPairs {

		wg.Add(1)
		sem <- struct{}{}
		go func(pair []string) {
			defer wg.Done()
			defer func() { <-sem }() // release slot

			ref := cfg.Reference
			out := cfg.OutputDir
			fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]

			AlignShortReadsMem(ref, fwd, rev, sn, lb, out, threadsPerSample)
		}(pair)
		//fmt.Println("Read Pair:", pair)
	}

}
