package bsaseq

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"os"
	"path/filepath"
	"strings"

	"runtime"
	"sync"
)

func RunBsaSeqFromConfig(
	configFile string,
	threads int,
	species string,
	bootstrap bool,
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
	rep int) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Species:", cfg.Species)
	fmt.Println("Read Pairs:", cfg.ReadPairs)
	fmt.Println("Threads", threads)
	fmt.Println("Output directory:", cfg.OutputDir)
	fmt.Println("Bams: ", cfg.Bams)
	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	knownSites := cfg.KnownSites

	refFile := cfg.Reference
	_, err = os.Stat(refFile)
	if err != nil {
		fmt.Printf("Reference file: %s is not a valid file path", refFile)
		return
	}

	outDir := cfg.OutputDir
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

	configBams := cfg.Bams
	readPairs := cfg.ReadPairs

	libMap := map[string]bool{"HIGH_BULK": true, "LOW_BULK": true, "HIGH_PARENT": true, "LOW_PARENT": true}

	libSampleMap := make(map[string]string)

	var rgmdBams []string
	var bqsrBams []string
	if len(readPairs) == 0 && len(configBams) == 0 {
		fmt.Println("You must provide at least one read pair or bam file")
		return
	} else if len(readPairs) > 0 && len(configBams) > 0 {
		fmt.Println("You must provide either read pairs or bam files not both!!")
		return
	} else if len(readPairs) > 0 {
		fmt.Println("Working with read pairs")
		//------------------------------------ Align reads to ref --------------------------------------------------- //
		if len(knownSites) == 0 && bootstrap == false {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return
		} else if len(knownSites) == 0 && bootstrap == true {
			fmt.Println("We will run BQSR using the bootstrap method")

		} else if len(knownSites) > 0 {
			fmt.Println("We will run BQSR using known-sites")
			// ------------------------ Checking Known sites file paths ----------------------------------------- //
			for j, _ := range knownSites {
				_, err := os.Stat(knownSites[j])
				if err != nil {
					fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
					log.Fatal(err)
				}
			}

		} else {
			fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
			return
		}

		fmt.Println("Aligning reads to reference")
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
				fmt.Printf("Library name %s is valid\n", lb)
				libSampleMap[lb] = sn
			}
		}

		var wg sync.WaitGroup
		sem := make(chan struct{}, maxParallelJobs)
		for _, pair := range cfg.ReadPairs {

			wg.Add(1)
			sem <- struct{}{}
			go func(pair []string) {
				defer wg.Done()
				defer func() { <-sem }() // release slot

				fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]
				lineDir := fmt.Sprintf("%s/%s", outDir, sn)
				rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sn)
				bqsrBam := strings.TrimSuffix(rgmdBam, ".bam") + "_bqsr.bam"
				rgmdBams = append(rgmdBams, rgmdBam)
				bqsrBams = append(bqsrBams, bqsrBam)

				alignment.AlignShortReadsMem(refFile, fwd, rev, sn, lb, outDir, threads)
			}(pair)

		}
		wg.Wait()

		// --------------------------------------- Recalibrating Bams ----------------------------------------------- //
		if len(knownSites) == 0 && bootstrap == false {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return
		} else if len(knownSites) == 0 && bootstrap == true {
			fmt.Println("Running with bootstrap method")
			alignment.BootstrapBqsr(refFile, rgmdBams, maxParallelJobs)
		} else if len(knownSites) > 0 {
			fmt.Println("Running with known-sites flag")
			// ---------------------------------- Running dbSnpBQSR ------------------------------------------------- //
			alignment.DbSnpBqsr(refFile, rgmdBams, knownSites, maxParallelJobs)

		} else {
			fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
			return
		}

		// ---------------------------------------- Variant Calling ------------------------------------------------- //
		variants.VariantCalling(refFile, bqsrBams, outDir, species, maxParallelJobs, "INFO")
		finalVcf := filepath.Join(outDir, species+"_"+".joint_hard_filtered.vcf.gz")

		// --------------------------------------------- BSAseq ----------------------------------------------------- //
		highParent, lowParent, highBulk, lowBulk := "", "", "", ""
		for _, lb := range libSampleMap {
			if lb == "HIGH_PARENT" {
				highParent = libSampleMap[lb]
			} else if lb == "LOW_PARENT" {
				lowParent = libSampleMap[lb]
			} else if lb == "HIGH_BULK" {
				highBulk = libSampleMap[lb]
			} else if lb == "LOW_BULK" {
				lowBulk = libSampleMap[lb]
			} else {
				fmt.Println("Something went wrong")
			}
		}
		if highParent == "" && lowParent == "" && highBulk != "" && lowBulk != "" {
			fmt.Println("Running 2 bulks only analysis")
			TwoBulkOnlyRun(finalVcf, highBulk, lowBulk, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)
		} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk != "" {
			fmt.Println("Running 2 bulks 2 parents analysis")
			TwoBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, lowBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)

		} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk == "" {
			fmt.Println("Running 1 high bulk, 2 parent analysis")
			outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_high_bsaseq_stats.tsv"
			OneBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, highBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

		} else if highParent != "" && lowParent != "" && highBulk == "" && lowBulk != "" {
			fmt.Println("Running 1 low bulk, 2 parent analysis")
			outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_low_bsaseq_stats.tsv"
			OneBulkTwoParentsRun(finalVcf, highParent, lowParent, lowBulk, minHighParentDepth, minLowParentDepth, minLowBulkDepth, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

		} else {
			log.Fatal("Invalid parameters. Valid combinations are:\n" +
				"1. Two bulks only: provide high_bulk (-A) and low_bulk (-B) without parents\n" +
				"2. Two bulks with two parents: provide high_parent (-H), low_parent (-L), high_bulk (-A), and low_bulk (-B)\n" +
				"3. One bulk with two parents: provide high_parent (-H), low_parent (-L), and bulk (-X)")
		}

	} else {
		fmt.Println("Starting from bam files")
		// ---------------------------------------- Variant Calling ------------------------------------------------- //
		variants.VariantCalling(refFile, configBams, outDir, species, maxParallelJobs, "INFO")
		finalVcf := filepath.Join(outDir, species+"_"+".joint_hard_filtered.vcf.gz")

		// --------------------------------------------- BSAseq ----------------------------------------------------- //
		highParent, lowParent, highBulk, lowBulk := "", "", "", ""
		for _, lb := range libSampleMap {
			if lb == "HIGH_PARENT" {
				highParent = libSampleMap[lb]
			} else if lb == "LOW_PARENT" {
				lowParent = libSampleMap[lb]
			} else if lb == "HIGH_BULK" {
				highBulk = libSampleMap[lb]
			} else if lb == "LOW_BULK" {
				lowBulk = libSampleMap[lb]
			} else {
				fmt.Println("Something went wrong")
			}
		}
		if highParent == "" && lowParent == "" && highBulk != "" && lowBulk != "" {
			fmt.Println("Running 2 bulks only analysis")
			TwoBulkOnlyRun(finalVcf, highBulk, lowBulk, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)
		} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk != "" {
			fmt.Println("Running 2 bulks 2 parents analysis")
			TwoBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, lowBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, minLowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)

		} else if highParent != "" && lowParent != "" && highBulk != "" && lowBulk == "" {
			fmt.Println("Running 1 high bulk, 2 parent analysis")
			outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_high_bsaseq_stats.tsv"
			OneBulkTwoParentsRun(finalVcf, highParent, lowParent, highBulk, minHighParentDepth, minLowParentDepth, minHighBulkDepth, highBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

		} else if highParent != "" && lowParent != "" && highBulk == "" && lowBulk != "" {
			fmt.Println("Running 1 low bulk, 2 parent analysis")
			outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_low_bsaseq_stats.tsv"
			OneBulkTwoParentsRun(finalVcf, highParent, lowParent, lowBulk, minHighParentDepth, minLowParentDepth, minLowBulkDepth, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

		} else {
			log.Fatal("Invalid parameters. Valid combinations are:\n" +
				"1. Two bulks only: provide high_bulk (-A) and low_bulk (-B) without parents\n" +
				"2. Two bulks with two parents: provide high_parent (-H), low_parent (-L), high_bulk (-A), and low_bulk (-B)\n" +
				"3. One bulk with two parents: provide high_parent (-H), low_parent (-L), and bulk (-X)")
		}

	}

}
