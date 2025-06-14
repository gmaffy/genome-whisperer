package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"os"
	"os/exec"
	"runtime"
	"sync"
)

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

func AlignShortReadsConfig(configPath string, threadsPerSample int) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
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

	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threadsPerSample = totalCores
	}

	fmt.Printf("Running up to %d jobs in parallel with %d threads each\n", maxParallelJobs, threadsPerSample)

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

	}
	wg.Wait()

}

func AlignShortReadsBt(ref string, fwd string, rev string, sn string, lb string, sampleDir string, threads int) error {
	fmt.Println("Aligning reads using bowtie2 ...")

	cmdStr := fmt.Sprintf(`bowtie2 -I 0 -X 1000 -x %s -1 %s -2 %s --end-to-end --sensitive --threads %v  --rg-id %s.1 --rg PL:BGISEQ --rg SM:%s --rg LB:%s | samtools sort -o %s/%s.sorted.bam`, ref, fwd, rev, threads, sn, sn, lb, sampleDir, sn)
	fmt.Println(cmdStr)
	bErr := utils.RunBashCmdVerbose(cmdStr)
	if bErr != nil {
		fmt.Printf("Error running bowtie2: %v\n", bErr)
		return bErr
	}

	fmt.Printf("Indexing Bam ....")
	mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s/%s.sorted.bam -O %s/%s.sorted.bai`, sampleDir, sn, sampleDir, sn)
	mErr := utils.RunBashCmdVerbose(mDupCmdStr)
	if mErr != nil {
		fmt.Printf("Error indexing bam file: %v\n", mErr)
		return mErr
	}
	return nil
}
