package bsaseq

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/go-echarts/go-echarts/v2/components"
	"golang.org/x/sync/errgroup"
	"log"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

// ========================================= Filtering ============================================================== //

func oneBulkTwoParTsvFilter(tsvFile string, highPar string, highParDP int, lowPar string, lowParDP int, bulk string, bulkDP int, winSize int, stepSize int, resultsDir string) ([]OneBulkTwoParentsRecord, error) {
	// ----------------------------------------------- Read to struct ------------------------------------------------//
	var filteredRecords []OneBulkTwoParentsRecord
	fmt.Printf("Reading VCF file %s ...\n\n", tsvFile)
	filterStart := time.Now()
	vcfStruct, err := readTsvToStructOneBulkTwoPar(tsvFile, highPar, lowPar, bulk)

	if err != nil {
		fmt.Println("Error reading tsv file: ", err)
		return filteredRecords, err
	}
	fmt.Printf("#Variants in original VCF: %d\n\n", len(vcfStruct))
	fmt.Printf("Removing short contigs ...\n\n")

	// --------------------------------------------- Remove short contigs --------------------------------------------//
	bsaStruct := removeShortContigsOneBulkTwoPar(vcfStruct, winSize, stepSize)

	bsaCount := len(bsaStruct)
	fmt.Printf("#Variants after removing short contigs: %d\n\n", bsaCount)

	fmt.Printf("Filtering variants ...\n\n")
	// ---------------------------------------------- Filter with Params -------------------------------------------- //

	for _, rec := range bsaStruct {
		hGT := rec.HighParGT
		hAD := rec.HighParAD
		hDP := rec.HighParDP
		hParts := strings.Split(hGT, "/")
		hADComms := strings.Count(hAD, ",")

		lGT := rec.LowParGT
		lAD := rec.LowParAD
		lDP := rec.LowParDP
		lParts := strings.Split(lGT, "/")
		lADComms := strings.Count(lAD, ",")

		BGT := rec.BulkGT
		BAD := rec.BulkAD
		BDP := rec.BulkDP
		BADComms := strings.Count(BAD, ",")

		if hGT != "./." && lGT != "./." && BGT != "./." {
			if len(hParts) == 2 && hParts[0] == hParts[1] && len(lParts) == 2 && lParts[0] == lParts[1] {
				if hGT != lGT {
					if hDP >= highParDP && lDP >= lowParDP && BDP >= bulkDP {

						if BADComms == 1 && hADComms == 1 && lADComms == 1 {
							filteredRecords = append(filteredRecords, rec)
						}
					}

				}

			}
		}
	}
	fmt.Println("#Variants after filtering: ", len(filteredRecords))

	filterEnd := time.Now()
	filterElapsed := filterEnd.Sub(filterStart)
	fmt.Printf("Filtering took %s\n", filterElapsed)

	// ---------------------------------------------- Write to file ------------------------------------------------- //

	writeOneBulkTwoPar(filteredRecords, highPar, lowPar, bulk, filepath.Join(resultsDir, "filtered.tsv"))
	fmt.Println("Filtered tsv file saved at: ", filepath.Join(resultsDir, "filtered.tsv"))
	return filteredRecords, nil

}

// ========================================= Statistics ============================================================= //
func oneBulkTwoParStats(filteredRecords []OneBulkTwoParentsRecord, bulkSize int, popStructure string, rep int) ([]OneBulkTwoParentsRecord, error) {
	fmt.Println("Calculating SNPIndex & Thresholds ... ")
	statsStart := time.Now()
	fmt.Println("simulating AF. ")
	smAF := simulateAF(popStructure, float64(bulkSize), rep)
	fmt.Println("simulated AF. ")
	var statsRecords = make([]OneBulkTwoParentsRecord, len(filteredRecords))
	fmt.Println("Calculating SNPIndex & Thresholds ... ")
	var g errgroup.Group

	for i := range filteredRecords {
		i := i
		g.Go(func() error {
			statsRecords[i] = calculateStatsRecordOne(filteredRecords[i], rep, smAF)
			return nil
		})
	}

	if err := g.Wait(); err != nil {
		fmt.Println("Stats calc Error: ", err)
		return nil, err
	}

	statsEnd := time.Now()
	statsElapsed := statsEnd.Sub(statsStart)
	fmt.Printf("Statistics and Thresholds took ... %s\n", statsElapsed)
	return statsRecords, nil

}

// ========================================= Plotting =============================================================== //
func slidingWindowAnalysisOne(snps []OneBulkTwoParentsRecord, variantType string, windowSize int, stepSize int) []OneBulkTwoParentsRecord {
	var results []OneBulkTwoParentsRecord
	if len(snps) == 0 {
		return results
	}

	fmt.Printf("Getting chromosome maps for variant type: %s ...\n\n", variantType)
	chromMap := make(map[string][]OneBulkTwoParentsRecord)
	for _, snp := range snps {
		if snp.Type == variantType {
			chromMap[snp.Chrom] = append(chromMap[snp.Chrom], snp)
		}
	}

	fmt.Printf("Performing sliding window analysis ...\n\n")

	for chrom, chromSNPs := range chromMap {
		fmt.Printf("Processing chromosome: %s ...\n\n", chrom)
		//fmt.Printf("Getting max position for %s ...\n\n", chrom)
		sort.Slice(chromSNPs, func(i, j int) bool {
			return chromSNPs[i].Pos < chromSNPs[j].Pos
		})

		maxPos := chromSNPs[len(chromSNPs)-1].Pos

		for start := 0.0; start <= maxPos; start += float64(stepSize) {
			end := start + float64(windowSize)
			var window []OneBulkTwoParentsRecord
			for _, snp := range chromSNPs {
				if snp.Pos >= start && snp.Pos < end {
					window = append(window, snp)
				}
			}

			if len(window) == 0 {
				continue
			}

			var (
				sumI, sump99, sump95, sump99L, sump95L, sumPos float64
			)
			for _, snp := range window {
				sumI += snp.SI
				sump99 += snp.P99
				sump95 += snp.P95
				sump99L += snp.Mp99
				sump95L += snp.Mp95
				sumPos += snp.Pos

			}
			windowCenter := (start + start + float64(windowSize)) / 2.0
			n := float64(len(window))
			results = append(results, OneBulkTwoParentsRecord{
				Chrom: chrom,
				Pos:   windowCenter,
				Type:  variantType,
				SI:    sumI / n,
				P99:   sump99 / n,
				P95:   sump95 / n,
				Mp99:  sump99L / n,
				Mp95:  sump95L / n,
			})
		}
	}
	fmt.Printf("Sliding window analysis done for all chromosomes... now to writing to file!\n\n")
	return results
}

func plottingChartsOne(SmoothedSnps []OneBulkTwoParentsRecord, outputHTML string, outputCSV string, smoothing bool) error {
	fmt.Printf("Creating charts ...\n\n")
	fmt.Printf("Getting chromosome maps ...\n\n")
	chromMap := make(map[string][]OneBulkTwoParentsRecord)
	for _, snp := range SmoothedSnps {
		chromMap[snp.Chrom] = append(chromMap[snp.Chrom], snp)
	}

	fmt.Printf("Sorting chromosomes ...\n")
	chroms := make([]string, 0, len(chromMap))
	for chrom := range chromMap {
		chroms = append(chroms, chrom)
	}
	sort.Strings(chroms)

	page := components.NewPage()

	qtlFile, _ := os.Create(outputCSV)

	defer qtlFile.Close()

	_, err := fmt.Fprintf(qtlFile, "CHROM\tQTLstart\tQTLend\tPEAK\tSTATISTIC\tTHRESHOLD\n")
	if err != nil {
		return err
	}

	for _, chrom := range chroms {
		chromSNPs := chromMap[chrom]
		var x []int
		var Iy, Ip99y, Ip95y, Ip99Ly, Ip95Ly []float64
		for _, snp := range chromSNPs {
			x = append(x, int(snp.Pos))
			Iy = append(Iy, snp.SI)
			Ip99y = append(Ip99y, snp.P99)
			Ip95y = append(Ip95y, snp.P95)
			Ip99Ly = append(Ip99Ly, snp.Mp99)
			Ip95Ly = append(Ip95Ly, snp.Mp95)
		}

		Chart := createLineChart(x, Iy, Ip99y, Ip95y, Ip99Ly, Ip95Ly, chrom+" SNP-index", "SNP-index", smoothing)
		_, i9peakY, i9start, i9end, i9Found := detectQtlPeaks(x, Iy, Ip99y)
		_, i5peakY, i5start, i5end, i5Found := detectQtlPeaks(x, Iy, Ip95y)
		_, i9LpeakY, i9Lstart, i9Lend, i9LFound := detectQtlValleys(x, Iy, Ip99Ly)
		_, i5LpeakY, i5Lstart, i5Lend, i5LFound := detectQtlValleys(x, Iy, Ip95Ly)
		page.AddCharts(Chart)

		if i9Found {
			_, hi9Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tSNP-index\t99per (peak)\n", chrom, i9start, i9end, i9peakY)
			if hi9Err != nil {
				return hi9Err
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (95perc confidene)", chrom, hi5start, hi5end, hi5peakY)
		}

		if i5Found {
			_, hi5Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tSNP-index\t95per (peak)\n", chrom, i5start, i5end, i5peakY)
			if hi5Err != nil {
				return hi5Err
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (95perc confidene)", chrom, hi5start, hi5end, hi5peakY)
		}

		if i9LFound {
			_, hi9lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tSNP-index\t99per (valley)\n", chrom, i9Lstart, i9Lend, i9LpeakY)
			if hi9lErr != nil {
				return hi9lErr
			}

		}
		if i5LFound {
			_, hi5lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tHigh-SNP-index\t95per (valley)\n", chrom, i5Lstart, i5Lend, i5LpeakY)
			if hi5lErr != nil {
				return hi5lErr
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (99perc valley)", chrom, hi9Lstart, hi9Lend, hi9LpeakY)
		}

	}

	f, err := os.Create(outputHTML)
	if err != nil {
		return err
	}
	return page.Render(f)
}

func OneBulkTwoParentsRun(
	vcfFile string,
	highParent string,
	lowParent string,
	bulk string,

	minHighParentDepth int,
	minLowParentDepth int,
	minBulkDepth int,

	bulkSize int,

	windowSize int,
	stepSize int,
	smoothing bool,

	popStructure string,
	rep int,
	outputName string,
	outputDir string,
) {
	// --------------------------------------------- Log file ------------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(outputDir, "bsaseq.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)
	logged := utils.ParseLogFile(logFilePath)

	jlog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
	slog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

	// ======================================== Create Results dir =================================================== #
	resultsDir, err := createResultsDir(outputDir)
	if err != nil {
		fmt.Println("Error creating results dir: ", err)
		return
	}

	statsFile := filepath.Join(resultsDir, outputName)
	// ======================================== Filter vcf file ====================================================== #
	fmt.Printf("================================== Filtering Start ======================================\n\n")
	var filteredRecords []OneBulkTwoParentsRecord
	var tsvFile string
	if strings.ToLower(filepath.Ext(vcfFile)) == ".vcf" || strings.ToLower(filepath.Ext(vcfFile)) == ".gz" {
		fmt.Printf("Working with VCF file ...\n\n")
		tsvFile = filepath.Join(resultsDir, "rawVariants.tsv")
		cmd3 := exec.Command("gatk", "VariantsToTable", "-V", vcfFile, "-F",
			"CHROM", "-F", "POS", "-F", "REF", "-F", "ALT", "-F", "QUAL", "-F", "TYPE", "-GF", "GT", "-GF", "AD", "-GF",
			"DP", "-GF", "GQ", "-O", tsvFile)

		varTabErr := cmd3.Run()
		if varTabErr != nil {
			log.Fatalf("error running gatk VariantsToTable: %s", varTabErr)
		}

	} else {
		fmt.Println("Working with tsv file ...")
		tsvFile = vcfFile

	}
	fmt.Printf("Filtering DF ...\n\n")

	if utils.StageHasCompleted(logged, "FILTERING", "ALL", "ALL") {
		fmt.Println("FILTERING has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		filteredRecords, err = oneBulkTwoParTsvFilter(vcfFile, highParent, minHighParentDepth, lowParent, minLowParentDepth, bulk, minBulkDepth, windowSize, stepSize, resultsDir)
		if err != nil {
			fmt.Println("Error filtering tsv file: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("ERROR: %s", err))
			slog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("ERROR: %s", err))
			return

		}
		jlog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		fmt.Printf("# Variants after filtering: %d\n\n", len(filteredRecords))
	}

	fmt.Printf("================================== Filtering End ======================================\n\n")

	// ============================================= STATISTICS ====================================================== #

	fmt.Printf("====================================== BSAseq Statistics Start ========================================== \n\n")
	var statsRecords []OneBulkTwoParentsRecord
	if utils.StageHasCompleted(logged, "STATS", "ALL", "ALL") {
		fmt.Println("STATS has already completed. Skipping.")

	} else {
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		statsRecords, err = oneBulkTwoParStats(filteredRecords, bulkSize, popStructure, rep)
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}

		fmt.Printf("Writing stats file to %s... \n\n", statsFile)
		err = writeOneBulkTwoPar(statsRecords, highParent, lowParent, bulk, statsFile)
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

	}

	fmt.Printf("====================================== BSAseq Statistics End ========================================== \n\n")
	// =============================================== PLOTTING ====================================================== #
	fmt.Printf("====================================== BSAseq Plotting Start ========================================== \n\n")
	fmt.Printf("Performing sliding window analysis for SNP ...\n\n")
	var slidingRecordsSNPs []OneBulkTwoParentsRecord
	var slidingRecordsINDELs []OneBulkTwoParentsRecord
	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "SNP", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slidingRecordsSNPs = slidingWindowAnalysisOne(statsRecords, "SNP", windowSize, stepSize)

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, "sliding_window_stats.tsv"))
		err = writeOneBulkTwoPar(slidingRecordsSNPs, highParent, lowParent, bulk, filepath.Join(resultsDir, "sliding_window_stats.tsv"))
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

	}

	if utils.StageHasCompleted(logged, "PLOTTING", "SNP", "ALL") {
		fmt.Println("Plotting is already done. Skipping ...")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		err = plottingChartsOne(slidingRecordsSNPs, filepath.Join(resultsDir, "goBSAseq_plot_snp.html"), filepath.Join(resultsDir, "goBSAseq_plot_snp.tsv"), smoothing)
		if err != nil {
			//fmt.Println("Error plotting charts: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

	}

	fmt.Printf("====================================== BSAseq Plotting Start ========================================== \n\n")
	fmt.Printf("Performing sliding window analysis for INDEL ...\n\n")
	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "INDEL", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slidingRecordsINDELs = slidingWindowAnalysisOne(statsRecords, "INDEL", windowSize, stepSize)

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, "sliding_window_stats_INDEL.tsv"))
		err = writeOneBulkTwoPar(slidingRecordsINDELs, highParent, lowParent, bulk, filepath.Join(resultsDir, "sliding_window_stats_INDEL.tsv"))
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	}

	if utils.StageHasCompleted(logged, "PLOTTING", "INDEL", "ALL") {
		fmt.Println("Plotting is already done. Skipping ...")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		err = plottingChartsOne(slidingRecordsINDELs, filepath.Join(resultsDir, "goBSAseq_plot_INDEL.html"), filepath.Join(resultsDir, "goBSAseq_plot_INDEL.tsv"), smoothing)
		if err != nil {
			fmt.Println("Error plotting charts: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		fmt.Printf("====================================== BSAseq Plotting End ========================================== \n\n")
	}

}
