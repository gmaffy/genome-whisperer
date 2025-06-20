package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sync"
)

func AlignShortReadsMem(referencePath string, forwardPath string, reversePath string, sampleName string, libName string, outputDir string, threads int) error {

	fmt.Println("Reading ...")
	// ----------------------------------------- Output Paths ------------------------------------------------------- //
	lineDir := fmt.Sprintf("%s/%s", outputDir, sampleName)
	rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sampleName)
	rgmdMetrics := fmt.Sprintf("%s/%s.RGMD.metrics.txt", lineDir, sampleName)
	rgmdIndex := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sampleName)
	sortedBam := fmt.Sprintf("%s/%s.sorted.bam", lineDir, sampleName)

	// ----------------------------------------- Run bwa mem -------------------------------------------------------- //
	bErr := os.MkdirAll(lineDir, 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
	}
	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	cmdStr := fmt.Sprintf(`bwa mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
	fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

	err := utils.RunBashCmdVerbose(cmdStr)
	if err != nil {
		return err
	}

	// ------------------------------------------- Mark Duplicates -------------------------------------------------- //
	fmt.Printf("Marking duplicates ....")
	mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -I %s -O %s -M %s`, sortedBam, rgmdBam, rgmdMetrics)
	fmt.Printf("%s\n-----------------------------------------------\n\n", mDupCmdStr)

	err = utils.RunBashCmdVerbose(mDupCmdStr)
	if err != nil {
		return err
	}

	// ------------------------------------------------- Index Bam -------------------------------------------------- //
	fmt.Printf("Index Bam ....")
	indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
	fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

	err = utils.RunBashCmdVerbose(indexCmdStr)
	if err != nil {
		return err
	}
	return nil

}

func AlignShortReadsConfig(configPath string, threadsPerSample int, knownSites []string, bqsr bool, bootstrap bool) {

	// ---------------------------------------- Check Paths --------------------------------------------------------- //
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}

	ref := cfg.Reference
	_, refErr := os.Stat(ref)
	if refErr != nil {
		fmt.Printf("Reference genome path: %s, is not valid\n", ref)
		return
	}

	out := cfg.OutputDir
	outInfo, outErr := os.Stat(out)

	if outErr != nil {

		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", out)
			if createErr := os.MkdirAll(out, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", out, createErr)
				return
			}
			fmt.Printf("Output directory %s created successfully.\n", out)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", out, outErr)
			return
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output Directory %s file path is not a directory\n", out)
		return
	}

	i := 0
	for _, pair := range cfg.ReadPairs {
		if len(pair) < 4 {
			fmt.Printf("This read pair is wrongly formated %s\n", pair)
			fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			continue
		}

		fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]

		_, fwdErr := os.Stat(fwd)
		_, revErr := os.Stat(rev)

		if fwdErr != nil {
			fmt.Printf("Forward reads path %s, is not valid\n", fwd)
			return
		}

		if revErr != nil {
			fmt.Printf("Reverse reads path %s, is not valid\n", rev)
			return
		}

		if sn == "" {
			fmt.Println("Please provide sample name ")
			fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			return
		}
		if lb == "" {
			fmt.Println("Please provide library name  ")
			fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			return
		}
		i++
	}

	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Output directory:", cfg.OutputDir)
	fmt.Printf("Running short read alignment for %v read pairs", i)

	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(out, "alignMem.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	//mw := io.MultiWriter(logFile, os.Stdout)
	//log.SetOutput(mw)
	fmt.Println("Log file created.")

	jsonHandler := slog.NewJSONHandler(logFile, nil)

	jlog := slog.New(jsonHandler)

	logged := utils.ParseLogFile(logFilePath)

	// ============================================== Run Alignments ================================================ //

	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threadsPerSample = totalCores
	}

	fmt.Printf("Running up to %d jobs in parallel with %d threads each\n", maxParallelJobs, threadsPerSample)

	var wg sync.WaitGroup
	sem := make(chan struct{}, maxParallelJobs)

	// ----------------------------------------------- Check Paths if bqsr ------------------------------------------ //
	if bqsr {
		fmt.Println("Skipping BQSR")
		if len(knownSites) == 0 && bootstrap == false {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return
		} else if len(knownSites) > 0 {
			fmt.Println("Running with known-sites flag")
			// ---------------------------- Checking Known sites file paths ----------------------------------------- //
			for j, _ := range knownSites {
				_, err := os.Stat(knownSites[j])
				if err != nil {
					fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
					log.Fatal(err)
				}
			}
			if bootstrap == true {
				fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
				return
			}
		}
	}
	//------------------------------------------ Run alignment ------------------------------------------------------ //
	var bams []string
	for _, pair := range cfg.ReadPairs {

		wg.Add(1)
		sem <- struct{}{}
		go func(pair []string) {
			defer wg.Done()
			defer func() { <-sem }()

			fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]
			lineDir := fmt.Sprintf("%s/%s", out, sn)
			rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sn)
			rgmdBai := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sn)

			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "STATUS", "STARTED")

			_, rgmdBamErr := os.Stat(rgmdBam)
			_, rgmdBaiErr := os.Stat(rgmdBai)

			isDone := utils.StageHasCompleted(logged, "BWA_MEM", sn, "ALL")
			if isDone && rgmdBamErr == nil && rgmdBaiErr == nil {
				msg := fmt.Sprintf("Bwa and MarkDuplicates already completed for %s. Skipping.\n\n------------------------------\n\n", sn)
				slog.Info(msg)

			} else {
				alErr := AlignShortReadsMem(ref, fwd, rev, sn, lb, out, threadsPerSample)
				if alErr != nil {
					jlog.Error("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED- %v", alErr))
					slog.Error("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED- %v", alErr))

					return
				}
				jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sn, "STATUS", "COMPLETED")
			}
			bams = append(bams, rgmdBam)

		}(pair)

	}
	wg.Wait()
	bqsrLogFile := filepath.Join(out, "bqsr.log")
	if bqsr {
		fmt.Println("Running BQSR")
		if len(knownSites) == 0 && bootstrap == true {
			fmt.Println("Running with bootstrap method")
			err := BootstrapBqsr(ref, bams, maxParallelJobs, bqsrLogFile)
			if err != nil {
				fmt.Println("BQSR failed", err)
				return
			}
		} else if len(knownSites) > 0 && bootstrap == false {
			err := DbSnpBqsr(ref, bams, knownSites, maxParallelJobs, bqsrLogFile)
			if err != nil {
				fmt.Println("BQSR failed", err)
				return
			}

		}
	}

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
