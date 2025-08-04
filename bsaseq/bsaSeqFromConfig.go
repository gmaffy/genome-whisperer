package bsaseq

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sync"

	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
)

func RunBsaSeqFromConfig(
	configFile string,
	threads int,
	species string,
	minHighParentDepth int,
	minLowParentDepth int,
	minHighBulkDepth int,
	minLowBulkDepth int,
	highBulkSize int,
	lowBulkSize int,
	windowSize int,
	stepSize int,
	smoothing bool,
	popStructure string,
	rep int,
	bootstrap bool,
	bqsr bool,
	caller string,
	merger string,
	aligner string,
	preset string,
	dvVer string,
	modelType string,
	gatkLogLevel string,
	verbose bool,
	alignmentFmt string) {

	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}

	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	// ------------------------------------------ Check Paths ------------------------------------------------------- //

	//knownSites := cfg.KnownSites
	refFile := cfg.Reference
	_, err = os.Stat(refFile)
	if err != nil {
		fmt.Printf("Reference file: %s is not a valid file path", refFile)
		return
	}

	outDir := cfg.OutputDir
	fmt.Printf("Output Directory: %s\n", outDir)
	outInfo, outErr := os.Stat(outDir)
	if outErr != nil {

		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", outDir)
			if createErr := os.MkdirAll(outDir, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", outDir, createErr)
				return
			}
			fmt.Printf("Output directory %s created successfully.\n", outDir)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", outDir, outErr)
			return
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output Directory %s file path is not a directory\n", outDir)
		return
	}

	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}

	configBams := cfg.BSAseqBams
	readPairs := cfg.ReadPairs
	//seReads := cfg.SeReads

	libMap := map[string]bool{"HIGH_BULK": true, "LOW_BULK": true, "HIGH_PARENT": true, "LOW_PARENT": true}

	libSampleMap := make(map[string]string)

	// --------------------------------------------- Log file ------------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(outDir, "bsaseq.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	logged := utils.ParseLogFile(logFilePath)

	// ----------------------------------------------- Check Paths if bqsr ------------------------------------------ //
	knownSites := cfg.KnownSites
	if bqsr {
		fmt.Println("----------------------------------------------- Check Paths if bqsr ------------------------------------------")
		if aligner == "pbmm2" {
			fmt.Println("We do not support BQSR for pbmm2 aligner. Please use bwa-mem or bowtie2 aligner or disable BQSR")
			return
		}
		if len(knownSites) == 0 && !bootstrap {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return
		} else if len(knownSites) > 0 {
			fmt.Println("Running with known-sites flag")
			// ---------------------------- Checking Known sites file paths ----------------------------------------- //
			for j := range knownSites {
				_, err := os.Stat(knownSites[j])
				if err != nil {
					fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
					return
				}
			}
			if bootstrap{
				fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
				return
			}
		}
	}

	var allBams []string
	finalVcf := ""

	jlog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
	slog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

	// ------------------------------------ Get parents and bulks & bam files names -------------------------------------------------- //
	if len(readPairs) == 0 && len(configBams) == 0 {
		fmt.Println("You must provide at least one read pair or bam file")
		return
	} else if len(readPairs) > 0 && len(configBams) > 0 {
		fmt.Println("You must provide either read pairs or bam files not both!!")
		return
	} else if len(readPairs) > 0 {
		fmt.Println("Working with read pairs")
		i := 0
		for _, pair := range readPairs {
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
				return
			}
			if lb == "" {
				fmt.Println("Please provide library name ")
				return
			}
			if _, ok := libMap[lb]; !ok {
				fmt.Printf("Library name %s is not valid\n", lb)
				fmt.Println("Valid library names are: HIGH_BULK, LOW_BULK, HIGH_PARENT, LOW_PARENT")
				return
			} else {
				libSampleMap[lb] = sn
				fmt.Printf("%s : %s \n", lb, sn)
				i++
			}
			if i > 4 {
				fmt.Println("You can only provide 4 read pairs")
				return
			}
			bam := ""

			if bqsr {
				lineDir := fmt.Sprintf("%s/%s", outDir, sn)
				bam = fmt.Sprintf("%s/%s.RGMD_bqsr.bam", lineDir, sn)
			} else {
				lineDir := fmt.Sprintf("%s/%s", outDir, sn)
				bam = fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sn)
			}
			allBams = append(allBams, bam)
		}

		fmt.Printf("BSAseq inputs fetched .....\n\n")
		fmt.Printf("Total read pairs: %d\n", len(readPairs))
		fmt.Printf("Bams to be created : %s\n\n", allBams)
		
		

		fmt.Printf("Aligning reads to ref ...........\n\n")

		fmt.Printf("Running up to %d jobs in parallel with %d threads each\n", maxParallelJobs, threads)
		

		if utils.StageHasCompleted(logged, "BSASEQ_ALIGNMENT", "ALL", "ALL") {
			fmt.Printf("Alignment already completed. Skipping ........ \n ")
			
		} else {
			jlog.Info("BSASEQ", "PROGRAM", "BSASEQ_ALIGNMENT", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("BSASEQ", "PROGRAM", "BSASEQ_ALIGNMENT", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

			
			var wg sync.WaitGroup
			sem := make(chan struct{}, maxParallelJobs)
			for _, pair := range cfg.ReadPairs {
				if len(pair) < 4 {
					fmt.Printf("This read pair is wrongly formated %s\n", pair)
					fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name>")
					continue
				}

				wg.Add(1)
				sem <- struct{}{}
				go func(pair []string) {
					defer wg.Done()
					defer func() { <-sem }()

					fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]
					
					jlog.Info("BSASEQ", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
					slog.Info("BSASEQ", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", "STARTED")

					isDone := utils.StageHasCompleted(logged, "PE_ALIGNMENT", sn, "ALL")
					if isDone {
						msg := fmt.Sprintf("%s and MarkDuplicates already completed for %s. Skipping.\n\n-------------------------------------------------------\n\n", aligner, sn)
						slog.Info(msg)
						

					} else {
						_, alErr := alignment.RunAlignReads(cfg.Reference, fwd, rev, "", sn, lb, cfg.OutputDir, threads, aligner, knownSites, bqsr, bootstrap, logFilePath, preset, gatkLogLevel, verbose, alignmentFmt)
						if alErr != nil {
							jlog.Error("BSASEQ", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", alErr))
							slog.Error("BSASEQ", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED - %v", alErr))
							return
						}
						jlog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
						slog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", "COMPLETED")
						
					}
				}(pair)

			}
			wg.Wait()

			jlog.Info("BSASEQ", "PROGRAM", "BSASEQ_ALIGNMENT", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("BSASEQ", "PROGRAM", "BSASEQ_ALIGNMENT", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			fmt.Printf("Alignment completed. Bams %s ...........\n\n", allBams)
			

		}

	} else {
		fmt.Println("Working with bam files")

		for _, bamPair := range configBams {
			if len(bamPair) != 2 {
				fmt.Printf("The BSAseq bam is wrongly formated - %s\n", bamPair)
				fmt.Println("Supply bam files in this format: BSAseqBam: <path to bam file> <bulk or parent>: (HIGH_PARENT, LOW_PARENT, HIGH_BULK, LOW_BULK)> ")
				continue
			}
			bam, bamType := bamPair[0], bamPair[1]
			if _, ok := libMap[bamType]; !ok {
				fmt.Printf("Bam types %s is not valid\n", bamType)
				fmt.Println("Valid bam type names are: HIGH_BULK, LOW_BULK, HIGH_PARENT, LOW_PARENT")
				return
			} else {
				libSampleMap[bamType] = bam
				fmt.Printf("Library name %s for sample %s is valid\n", bamType, bam)
				allBams = append(allBams, bam)
			}
		}

	}

	// ---------------------------------------- Variant Calling ------------------------------------------------- //
	if utils.StageHasCompleted(logged, "BSASEQ_VARIANT_CALLING", "ALL", "ALL") {
		fmt.Printf("Variant calling already completed. Skipping ........ \n ")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

		hardFilteredVCF, err := variants.VariantCalling(refFile, allBams, outDir, species, threads, gatkLogLevel, caller, merger, logFilePath, dvVer, modelType, verbose)
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", err))
			return
		}
		jlog.Error("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Error("BSASEQ", "PROGRAM", "BSASEQ_VARIANT_CALLING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		finalVcf = hardFilteredVCF

	}

	fmt.Println("VARIANT CALLING DONE STARING BSAseq")
	fmt.Println(libSampleMap)

	highBulk := libSampleMap["HIGH_BULK"]
	lowBulk := libSampleMap["LOW_BULK"]
	highParent := libSampleMap["HIGH_PARENT"]
	lowParent := libSampleMap["LOW_PARENT"]
	fmt.Println("highParent: ", highParent, "lowParent: ", lowParent, "highBulk: ", highBulk, "lowBulk: ", lowBulk)
	if highParent == "" && lowParent == "" && highBulk != "" && lowBulk != "" {
		fmt.Println("Running 2 bulks only analysis")
		TwoBulkOnlyRun(finalVcf, highBulk, lowBulk, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outDir)
	} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk != "" {
		fmt.Println("Running 2 bulks 2 parents analysis")
		TwoBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, lowBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outDir)

	} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk == "" {
		fmt.Println("Running 1 high bulk, 2 parent analysis")
		outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_high_bsaseq_stats.tsv"
		OneBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, highBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName, outDir)

	} else if highParent != "" && lowParent != "" && highBulk == "" && lowBulk != "" {
		fmt.Println("Running 1 low bulk, 2 parent analysis")
		outputName := highParent + "_samp_" + lowParent + "_samp_" + lowBulk + "_samp_low_bsaseq_stats.tsv"
		OneBulkTwoParentsRun(finalVcf, highParent, lowParent, lowBulk, minHighParentDepth, minLowParentDepth, minLowBulkDepth, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName, outDir)

	} else {
		log.Fatal("Invalid parameters. Valid combinations are:\n" +
			"1. Two bulks only: provide high_bulk (-A) and low_bulk (-B) without parents\n" +
			"2. Two bulks with two parents: provide high_parent (-H), low_parent (-L), high_bulk (-A), and low_bulk (-B)\n" +
			"3. One bulk with two parents: provide high_parent (-H), low_parent (-L), and bulk (-X)")
	}

}
