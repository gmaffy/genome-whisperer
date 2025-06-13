package bsaseq

import (
	"fmt"
	"github.com/go-gota/gota/dataframe"
	"github.com/go-gota/gota/series"

	"os"
	"os/exec"

	"path/filepath"
	"sort"

	"strings"
)

func InteractiveRun(vcfFile string, popStructure string, rep int) {
	// ======================================== Create Results dir =================================================== #
	resultsDir := createResultsDir()
	// =========================================== Convert VCF to Tsv ================================================ #
	fmt.Printf("Converting table to dataframe ...\n\n")
	var tsvFile string
	if strings.ToLower(filepath.Ext(vcfFile)) == ".vcf" || strings.ToLower(filepath.Ext(vcfFile)) == ".gz" {
		fmt.Printf("Working with VCF file ...\n\n")
		fmt.Printf("Converting VCF file to table format ...\n\n")
		tsvFile = filepath.Join(resultsDir, "raw_table.tsv")
		cmd := exec.Command("gatk", "VariantsToTable", "-V", vcfFile, "-F",
			"CHROM", "-F", "POS", "-F", "REF", "-F", "ALT", "-F", "QUAL", "-F", "TYPE", "-GF", "GT", "-GF", "AD", "-GF",
			"DP", "-GF", "GQ", "-O", tsvFile)

		varTabErr := cmd.Run()
		if varTabErr != nil {
			fmt.Println("error running gatk SelectVariants: ", varTabErr)
			return
		}

	} else {
		fmt.Println("Working with tsv file ...")
		tsvFile = vcfFile
	}

	// ============================================== Creating DF ==================================================== #
	fmt.Printf("Reading data from table file  %s\n", tsvFile)

	vcfTableFileOpen, vcfErr := os.Open(tsvFile)
	if vcfErr != nil {
		fmt.Printf("Error opening vcf file %s: %s\n", tsvFile, vcfErr)
		panic(vcfErr)
		return
	}
	defer func(vcfTableFileOpen *os.File) {
		voErr := vcfTableFileOpen.Close()
		if voErr != nil {
			fmt.Printf("Error closing vcf file %s: %s\n", vcfFile, voErr)
			panic(voErr)
			return
		}
	}(vcfTableFileOpen)

	vcfDF := dataframe.ReadCSV(vcfTableFileOpen, dataframe.WithDelimiter('\t'))

	fmt.Printf("Dataframe created. Here: %s \n\n", vcfDF)

	// ======================================== Replacing | with / for GT columns ==================================== #

	fmt.Printf("Replacing | with / for GT columns ...\n\n")
	sampleParametersDic := map[int]string{0: "None"}
	sampleID := 0
	var ids []int
	ids = append(ids, 0)
	for _, colName := range vcfDF.Names() {
		if strings.HasSuffix(colName, ".GT") {
			sampleID++
			ids = append(ids, sampleID)
			sampleParametersDic[sampleID] = strings.TrimSuffix(colName, ".GT")
			updatedGTs := vcfDF.Col(colName).Records()
			for i := range updatedGTs {
				updatedGTs[i] = strings.ReplaceAll(updatedGTs[i], "|", "/")
			}
			vcfDF = vcfDF.Mutate(series.New(updatedGTs, series.String, colName))
		}
	}

	// ============================================ SELECT SAMPLES =================================================== #

	fmt.Printf("================================ SAMPLE SELECTION ==========================================\n\n")

	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	sort.Ints(ids)
	for id := range ids {
		fmt.Printf("%v : %v\n", id, sampleParametersDic[id])
	}

	fmt.Printf("Enter number corresponding to the sample ...\n\n")

	var highParentChoice int
	var lowParentChoice int
	var highBulkChoice int
	var lowBulkChoice int

	var highParentDepth int
	var lowParentDepth int
	var highBulkDepth int
	var lowBulkDepth int

	var highBulkSize int
	var lowBulkSize int

	var windowSize int
	var stepSize int

	var smoothing bool

	fmt.Printf("------------------------------------- PARENT CHOICES --------------------------------------\n\n")

	smoothing = true

	fmt.Println("Enter HIGH PARENT number:")
	_, highParErr := fmt.Scan(&highParentChoice)
	if highParErr != nil {
		fmt.Printf("HIGH PARENT number should be numerical and part of the list above: %s\n", highParErr)
		return
	}

	//if highParentChoice == 0 {
	//	fmt.Printf("HIGH PARENT number should not be None\n")
	//	return
	//}

	highParent := sampleParametersDic[highParentChoice]
	fmt.Printf("HIGH Parent is: %s \n\n", highParent)

	fmt.Println("Enter LOW PARENT number:")
	_, lowParErr := fmt.Scan(&lowParentChoice)
	if lowParErr != nil {
		fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", lowParErr)
		return
	}

	//if lowParentChoice == 0 {
	//	fmt.Printf("LOW PARENT number should not be None\n")
	//	return
	//}

	if lowParentChoice == highParentChoice && lowParentChoice != 0 {
		fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
		return
	}

	lowParent := sampleParametersDic[lowParentChoice]
	fmt.Printf("LOW parent is: %s \n\n", lowParent)

	fmt.Printf("------------------------------------- BULK CHOICES ----------------------------------------\n\n")
	fmt.Println("Enter HIGH BULK number:")
	_, highBulkErr := fmt.Scan(&highBulkChoice)
	if highBulkErr != nil {
		fmt.Printf("HIGH BULK number should be numerical and part of the list above: %s\n", highBulkErr)
		return
	}

	if highBulkChoice == highParentChoice || highBulkChoice == lowParentChoice {
		fmt.Println("Your HIGH bulk cannot be the same as any of the parents")
		return
	}

	highBulk := sampleParametersDic[highBulkChoice]
	fmt.Printf("HIGH bulk is: %s \n\n", highBulk)

	fmt.Println("Enter LOW BULK number:")
	_, lowBulkErr := fmt.Scan(&lowBulkChoice)
	if lowBulkErr != nil {
		fmt.Printf("LOW BULK number should be numerical and part of the list above: %s\n", lowBulkErr)
		return
	}

	if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
		fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
		return
	}
	lowBulk := sampleParametersDic[lowBulkChoice]
	fmt.Printf("LOW bulk is: %s \n\n", lowBulk)

	fmt.Printf("----------------------------------- SAMPLE MIN DEPTH --------------------------------------\n\n")

	fmt.Printf("------------------------------------- PARENT DEPTHS ---------------------------------------\n\n")

	if highParentChoice != 0 {
		fmt.Printf("Enter minimum depth for HIGH PARENT %s (integer): \n", highParent)
		_, highParDepErr := fmt.Scan(&highParentDepth)
		if highParDepErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("HIGH parent (%s) min depth is: %v \n\n", highParent, highParentDepth)
	}

	if lowBulkDepth != 0 {
		fmt.Printf("Enter minimum depth for LOW PARENT %s (integer): \n", lowParent)
		_, lowParDepErr := fmt.Scan(&lowParentDepth)
		if lowParDepErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW parent (%s) min depth: %v \n\n", lowParent, lowParentDepth)

	}

	fmt.Printf("--------------------------------------- BULK DEPTHS ---------------------------------------\n\n")

	if lowBulkChoice != 0 && highBulkChoice != 0 && highParentChoice == 0 && lowParentChoice == 0 {
		fmt.Println("Runnng bulks only")
		fmt.Printf("Enter minimum depth for HIGH BULK %s (integer): \n", highBulk)

		_, highBulkDepErr := fmt.Scan(&highBulkDepth)
		if highBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("HIGH bulk (%s) minimum depth is: %v \n\n", highBulk, highBulkDepth)

		fmt.Printf("Enter minimum depth for LOW BULK %s (integer): \n", lowBulk)
		_, lowBulkDepErr := fmt.Scan(&lowBulkDepth)
		if lowBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW bulk (%s) minimum depth is: %v \n", lowBulk, lowBulkDepth)

		fmt.Printf("\n\n--------------------------- BULK SIZES -----------------------------------\n\n")

		fmt.Printf("Enter BULK SIZE for HIGH BULK: %s (integer): \n", highBulk)
		_, highBulkSizeErr := fmt.Scan(&highBulkSize)
		if highBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("(%s) Bulk size is: %v \n\n", highBulk, highBulkSize)

		fmt.Printf("Enter BULK SIZE for LOW BULK %s (integer): \n", lowBulk)
		_, lowBulkSizeErr := fmt.Scan(&lowBulkSize)
		if lowBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW BULK (%s) DEPTH: %v \n\n", lowBulk, lowBulkSize)

		fmt.Printf("\n\n---------------------------------- PLOTTING PARAMETERS -------------------------------------\n\n")

		fmt.Printf("Enter WINDOW SIZE for PLOTTING: (integer, e.g 2000000): \n")
		_, winSizeErr := fmt.Scan(&windowSize)
		if winSizeErr != nil {
			fmt.Println("Window size value must be numerical")
			return
		}

		fmt.Printf("Enter STEP SIZE for PLOTTING: (integer, e.g 10000): \n")
		_, stepSizeErr := fmt.Scan(&stepSize)
		if stepSizeErr != nil {
			fmt.Println("Step size value must be numerical")
			return
		}

		fmt.Println("LETS RUN ....")
		TwoBulkOnlyRun(vcfFile, highBulk, lowBulk, highBulkDepth, lowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)
	} else if lowBulkChoice == 0 && highBulkChoice != 0 && lowParentChoice != 0 && highParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
		fmt.Printf("Enter minimum depth for HIGH BULK %s (integer): \n", highBulk)
		_, highBulkDepErr := fmt.Scan(&highBulkDepth)
		if highBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("HIGH bulk (%s) minimum depth is: %v \n\n", highBulk, highBulkDepth)

		fmt.Printf("\n\n---------------------------------- BULK SIZES -------------------------------------\n\n")

		fmt.Printf("Enter BULK SIZE for HIGH BULK: %s (integer): \n", highBulk)
		_, highBulkSizeErr := fmt.Scan(&highBulkSize)
		if highBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("(%s) Bulk size is: %v \n\n", highBulk, highBulkSize)
		fmt.Printf("Filtering VCF file ....\n\n")

		fmt.Printf("\n\n---------------------------------- PLOTTING PARAMETERS -------------------------------------\n\n")

		fmt.Printf("Enter WINDOW SIZE for PLOTTING: (integer, e.g 2000000): \n")
		_, winSizeErr := fmt.Scan(&windowSize)
		if winSizeErr != nil {
			fmt.Println("Window size value must be numerical")
			return
		}

		fmt.Printf("Enter STEP SIZE for PLOTTING: (integer, e.g 10000): \n")
		_, stepSizeErr := fmt.Scan(&stepSize)
		if stepSizeErr != nil {
			fmt.Println("Step size value must be numerical")
			return
		}
		outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_high_bsaseq_stats.tsv"
		OneBulkTwoParentsRun(vcfFile, highParent, lowParent, highBulk, highParentDepth, lowParentDepth, highBulkDepth, highBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

	} else if highBulkChoice == 0 && highParentChoice != 0 && lowParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (LOW bulk)...")
		fmt.Printf("Enter minimum depth for LOW BULK %s (integer): \n ", lowBulk)
		_, lowBulkDepErr := fmt.Scan(&lowBulkDepth)
		if lowBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW bulk (%s) minimum depth is: %v \n\n", lowBulk, lowBulkDepth)

		fmt.Printf("\n\n------------------------------------- BULK SIZES ----------------------------------\n\n")

		fmt.Printf("Enter BULK SIZE for LOW BULK %s (integer): \n", lowBulk)
		_, lowBulkSizeErr := fmt.Scan(&lowBulkSize)
		if lowBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW BULK (%s) DEPTH: %v \n\n", lowBulk, lowBulkSize)
		fmt.Printf("Filtering VCF file ....\n\n")

		fmt.Printf("\n\n---------------------------------- PLOTTING PARAMETERS -------------------------------------\n\n")

		fmt.Printf("Enter WINDOW SIZE for PLOTTING: (integer, e.g 2000000): \n")
		_, winSizeErr := fmt.Scan(&windowSize)
		if winSizeErr != nil {
			fmt.Println("Window size value must be numerical")
			return
		}

		fmt.Printf("Enter STEP SIZE for PLOTTING: (integer, e.g 10000): \n")
		_, stepSizeErr := fmt.Scan(&stepSize)
		if stepSizeErr != nil {
			fmt.Println("Step size value must be numerical")
			return
		}
		outputName := highParent + "_samp_" + lowParent + "_samp_" + highBulk + "_samp_low_bsaseq_stats.tsv"
		OneBulkTwoParentsRun(vcfFile, highParent, lowParent, lowBulk, highParentDepth, lowParentDepth, lowBulkDepth, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep, outputName)

	} else {
		fmt.Println("Working with two bulks")
		fmt.Printf("Enter minimum depth for HIGH BULK %s (integer): \n", highBulk)

		_, highBulkDepErr := fmt.Scan(&highBulkDepth)
		if highBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("HIGH bulk (%s) minimum depth is: %v \n\n", highBulk, highBulkDepth)

		fmt.Printf("Enter minimum depth for LOW BULK %s (integer): \n", lowBulk)
		_, lowBulkDepErr := fmt.Scan(&lowBulkDepth)
		if lowBulkDepErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW bulk (%s) minimum depth is: %v \n", lowBulk, lowBulkDepth)

		fmt.Printf("\n\n--------------------------- BULK SIZES -----------------------------------\n\n")

		fmt.Printf("Enter BULK SIZE for HIGH BULK: %s (integer): \n", highBulk)
		_, highBulkSizeErr := fmt.Scan(&highBulkSize)
		if highBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
			return
		}

		fmt.Printf("(%s) Bulk size is: %v \n\n", highBulk, highBulkSize)

		fmt.Printf("Enter BULK SIZE for LOW BULK %s (integer): \n", lowBulk)
		_, lowBulkSizeErr := fmt.Scan(&lowBulkSize)
		if lowBulkSizeErr != nil {
			fmt.Println("Depth value must be numerical")
		}

		fmt.Printf("LOW BULK (%s) DEPTH: %v \n\n", lowBulk, lowBulkSize)

		fmt.Printf("\n\n---------------------------------- PLOTTING PARAMETERS -------------------------------------\n\n")

		fmt.Printf("Enter WINDOW SIZE for PLOTTING: (integer, e.g 2000000): \n")
		_, winSizeErr := fmt.Scan(&windowSize)
		if winSizeErr != nil {
			fmt.Println("Window size value must be numerical")
			return
		}

		fmt.Printf("Enter STEP SIZE for PLOTTING: (integer, e.g 10000): \n")
		_, stepSizeErr := fmt.Scan(&stepSize)
		if stepSizeErr != nil {
			fmt.Println("Step size value must be numerical")
			return
		}

		fmt.Println("LETS RUN ....")

		TwoBulkTwoParentsRun(vcfFile, highParent, lowParent, highBulk, lowBulk, highParentDepth, lowParentDepth, highBulkDepth, lowBulkDepth, highBulkSize, lowBulkSize, windowSize, stepSize, smoothing, popStructure, rep)

	}
}
