package bsaseq

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
	"golang.org/x/sync/errgroup"
	"io"
	"log"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

// =================================================== Filtering ==================================================== //
func twoBulkTwoParTsvFilter(tsvFile string, highPar string, highParDP int, lowPar string, lowParDP int, highBulk string, highBulkDP int, lowBulk string, lowBulkDP int, winSize int, stepSize int, resultsDir string, logFile io.Writer) ([]TwoBulkTwoParentsRecord, error) {

	var filteredRecords []TwoBulkTwoParentsRecord
	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)
	// ----------------------------------------------- Read to struct ------------------------------------------------//
	fmt.Printf("Reading VCF file %s ...\n\n", tsvFile)
	filterStart := time.Now()

	jlog.Info("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	vcfStruct, err := readTsvToStructTwoBulkTwoPar(tsvFile, highPar, lowPar, highBulk, lowBulk)
	if err != nil {
		jlog.Error("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		slog.Error("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		return filteredRecords, err
	}
	jlog.Info("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("BSASEQ", "PROGRAM", "READ_TSV", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

	fmt.Printf("#Variants in original VCF: %d\n\n", len(vcfStruct))

	// --------------------------------------------- Remove short contigs --------------------------------------------//
	fmt.Printf("Removing short contigs ...\n\n")
	jlog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	bsaStruct, err := removeShortContigsTwoBulkTwoPar(vcfStruct, winSize, stepSize)
	if err != nil {
		jlog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		slog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		return nil, err

	}
	jlog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("BSASEQ", "PROGRAM", "REMOVE_SHORT_CONTIGS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

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

		hBGT := rec.HighBulkGT
		hBAD := rec.HighBulkAD
		hBDP := rec.HighBulkDP
		hBADComms := strings.Count(hBAD, ",")

		lBGT := rec.LowBulkGT
		lBAD := rec.LowBulkAD
		lBDP := rec.LowBulkDP
		lBADComms := strings.Count(lBAD, ",")

		if hGT != "./." && lGT != "./." && hBGT != "./." && lBGT != "./." {
			if len(hParts) == 2 && hParts[0] == hParts[1] && len(lParts) == 2 && lParts[0] == lParts[1] {
				if hGT != lGT {
					if hDP >= highParDP && lDP >= lowParDP && hBDP >= highBulkDP && lBDP >= lowBulkDP {

						if hBADComms == 1 && lBADComms == 1 && hADComms == 1 && lADComms == 1 {
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

	err = writeTwoBulkTwoPar(filteredRecords, highPar, lowPar, highBulk, lowBulk, filepath.Join(resultsDir, "filtered.tsv"))
	if err != nil {
		jlog.Info("BSASEQ", "PROGRAM", "WRITE_TO_TSV", "SAMPLE", "filtered.tsv", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		slog.Info("BSASEQ", "PROGRAM", "WRITE_TO_TSV", "SAMPLE", "filtered.tsv", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		return nil, err

	}
	fmt.Println("Filtered tsv file saved at: ", filepath.Join(resultsDir, "filtered.tsv"))
	return filteredRecords, nil

}

func twoBulksOnlyFilter(tsvFile string, highBulk string, highBulkDP int, lowBulk string, lowBulkDP int, winSize int, stepSize int, resultsDir string) ([]TwoBulkTwoParentsRecord, error) {
	// ----------------------------------------------- Read to struct ------------------------------------------------//
	fmt.Printf("Reading VCF file %s ...\n\n", tsvFile)
	filterStart := time.Now()
	vcfStruct, err := readTsvToStructTwoBulkOnly(tsvFile, highBulk, lowBulk)
	if err != nil {
		fmt.Println("Error reading VCF file: ", err)
		return nil, err
	}
	fmt.Printf("#Variants in original VCF: %d\n\n", len(vcfStruct))
	fmt.Printf("Removing short contigs ...\n\n")

	// --------------------------------------------- Remove short contigs --------------------------------------------//
	bsaStruct, err := removeShortContigsTwoBulkOnly(vcfStruct, winSize, stepSize)
	if err != nil {
		fmt.Println("Error removing short contigs: ", err)
		return nil, err
	}

	bsaCount := len(bsaStruct)
	fmt.Printf("#Variants after removing short contigs: %d\n\n", bsaCount)

	fmt.Printf("Filtering variants ...\n\n")
	// ---------------------------------------------- Filter with Params -------------------------------------------- //
	var filteredRecords []TwoBulkTwoParentsRecord
	for _, rec := range bsaStruct {

		hBGT := rec.HighBulkGT
		hBAD := rec.HighBulkAD
		hBDP := rec.HighBulkDP
		hBADComms := strings.Count(hBAD, ",")

		lBGT := rec.LowBulkGT
		lBAD := rec.LowBulkAD
		lBDP := rec.LowBulkDP
		lBADComms := strings.Count(lBAD, ",")

		if hBGT != "./." && lBGT != "./." {
			if hBDP >= highBulkDP && lBDP >= lowBulkDP {
				if hBADComms == 1 && lBADComms == 1 {
					filteredRecords = append(filteredRecords, rec)
				}
			}

		}
	}
	fmt.Println("#Variants after filtering: ", len(filteredRecords))

	filterEnd := time.Now()
	filterElapsed := filterEnd.Sub(filterStart)
	fmt.Printf("Filtering took %s\n", filterElapsed)

	// ---------------------------------------------- Write to file ------------------------------------------------- //

	err = writeTwoBulkOnly(filteredRecords, highBulk, lowBulk, filepath.Join(resultsDir, "filtered.tsv"))
	if err != nil {
		fmt.Println("Error writing to file: ", err)
		return nil, err
	}
	fmt.Println("Filtered tsv file saved at: ", filepath.Join(resultsDir, "filtered.tsv"))
	return filteredRecords, nil

}

// =================================================== Statistics =================================================== //

func twoBulkTwoParStats(filteredRecords []TwoBulkTwoParentsRecord, highBulkSize int, lowBulkSize int, popStructure string, rep int) ([]TwoBulkTwoParentsRecord, error) {
	fmt.Println("Calculating BSA-seq Statistics SNPIndex, Delta SNPIndex, G-Statistics & Thresholds ... ")
	statsStart := time.Now()
	highSmAF := simulateAF(popStructure, float64(highBulkSize), rep)
	lowSmAF := simulateAF(popStructure, float64(lowBulkSize), rep)

	var statsRecords = make([]TwoBulkTwoParentsRecord, len(filteredRecords))
	var g errgroup.Group

	for i := range filteredRecords {
		i := i
		g.Go(func() error {
			statsRecords[i] = calculateStatsRecord(filteredRecords[i], rep, highSmAF, lowSmAF)
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

func twoBulkOnlyStats(filteredRecords []TwoBulkTwoParentsRecord, highBulkSize int, lowBulkSize int, popStructure string, rep int) ([]TwoBulkTwoParentsRecord, error) {
	fmt.Println("Calculating BSA-seq Statistics SNPIndex, Delta SNPIndex, G-Statistics & Thresholds ... ")
	statsStart := time.Now()
	highSmAF := simulateAF(popStructure, float64(highBulkSize), rep)
	lowSmAF := simulateAF(popStructure, float64(lowBulkSize), rep)

	var statsRecords = make([]TwoBulkTwoParentsRecord, len(filteredRecords))
	var g errgroup.Group

	for i := range filteredRecords {
		i := i
		g.Go(func() error {
			statsRecords[i] = calculateStatsRecordBulksOnly(filteredRecords[i], rep, highSmAF, lowSmAF)
			return nil
		})
	}

	if err := g.Wait(); err != nil {
		fmt.Println("Stats calc Error: ", err)
	}

	statsEnd := time.Now()
	statsElapsed := statsEnd.Sub(statsStart)
	fmt.Printf("Statistics and Thresholds took ... %s\n", statsElapsed)
	return statsRecords, nil

}

// =============================================== Plotting ========================================================= //
func slidingWindowAnalysis(snps []TwoBulkTwoParentsRecord, variantType string, windowSize int, stepSize int) []TwoBulkTwoParentsRecord {
	var results []TwoBulkTwoParentsRecord
	if len(snps) == 0 {
		return results
	}

	fmt.Printf("Getting chromosome maps for variant type: %s ...\n\n", variantType)
	chromMap := make(map[string][]TwoBulkTwoParentsRecord)
	for _, snp := range snps {
		if snp.Type == variantType {
			chromMap[snp.Chrom] = append(chromMap[snp.Chrom], snp)
		}
	}

	fmt.Printf("Performing sliding window analysis ...\n\n")

	for chrom, chromSNPs := range chromMap {

		sort.Slice(chromSNPs, func(i, j int) bool {
			return chromSNPs[i].Pos < chromSNPs[j].Pos
		})

		maxPos := chromSNPs[len(chromSNPs)-1].Pos

		for start := 0.0; start <= maxPos; start += float64(stepSize) {
			end := start + float64(windowSize)
			var window []TwoBulkTwoParentsRecord
			for _, snp := range chromSNPs {
				if snp.Pos >= start && snp.Pos < end {
					window = append(window, snp)
				}
			}

			if len(window) == 0 {
				continue
			}

			var (
				sumDS, sumDSp99, sumDSp95, sumDSp99L, sumDSp95L, sumGS, sumGSp99, sumGSp95, sumHI, sumHp99, sumHp95, sumHp99L, sumHp95L, sumLI, sumLp99, sumLp95, sumLp99L, sumLp95L, sumPos float64
			)
			for _, snp := range window {
				sumDS += snp.DSI
				sumDSp99 += snp.DsiP99
				sumDSp95 += snp.DsiP95
				sumDSp99L += snp.DsiMp99
				sumDSp95L += snp.DsiMp95

				sumGS += snp.GS
				sumGSp99 += snp.GsP99
				sumGSp95 += snp.GsP95

				sumHI += snp.HighSI
				sumHp99 += snp.HighP99
				sumHp95 += snp.HighP95
				sumHp99L += snp.HighMp99
				sumHp95L += snp.HighMp95

				sumLI += snp.LowSI
				sumLp99 += snp.LowP99
				sumLp95 += snp.LowP95
				sumLp99L += snp.LowMp99
				sumLp95L += snp.LowMp95
				sumPos += snp.Pos

			}
			windowCenter := (start + start + float64(windowSize)) / 2.0
			n := float64(len(window))
			results = append(results, TwoBulkTwoParentsRecord{
				Chrom:    chrom,
				Pos:      windowCenter,
				Type:     variantType,
				HighSI:   sumHI / n,
				HighP99:  sumHp99 / n,
				HighP95:  sumHp95 / n,
				HighMp99: sumHp99L / n,
				HighMp95: sumHp95L / n,

				LowSI:   sumLI / n,
				LowP99:  sumLp99 / n,
				LowP95:  sumLp95 / n,
				LowMp99: sumLp99L / n,
				LowMp95: sumLp95L / n,

				DSI:     sumDS / n,
				DsiP99:  sumDSp99 / n,
				DsiP95:  sumDSp95 / n,
				DsiMp99: sumDSp99L / n,
				DsiMp95: sumDSp95L / n,

				GS:    sumGS / n,
				GsP99: sumGSp99 / n,
				GsP95: sumGSp95 / n,
			})
		}
	}
	fmt.Printf("Sliding window analysis done for all chromosomes... now to writing to file!\n\n")
	return results
}

func createLineChart(x []int, y []float64, y9 []float64, y5 []float64, y9L []float64, y5L []float64, title, ylabel string, smoothing bool) *charts.Line {
	line := charts.NewLine()
	line.SetGlobalOptions(
		charts.WithInitializationOpts(opts.Initialization{Theme: types.ThemeWesteros}),
		charts.WithTitleOpts(opts.Title{Title: title}),
		charts.WithYAxisOpts(opts.YAxis{Name: ylabel}),
		charts.WithXAxisOpts(opts.XAxis{Name: "Position (bp)"}),
	)
	var yData []opts.LineData
	for _, v := range y {
		yData = append(yData, opts.LineData{Value: v})
	}

	var y9Data []opts.LineData
	for _, v9 := range y9 {
		y9Data = append(y9Data, opts.LineData{Value: v9})
	}

	var y5Data []opts.LineData
	for _, v5 := range y5 {
		y5Data = append(y5Data, opts.LineData{Value: v5})
	}

	var y9LData []opts.LineData
	for _, v9L := range y9L {
		y9LData = append(y9LData, opts.LineData{Value: v9L})
	}

	var y5LData []opts.LineData
	for _, v5L := range y5L {
		y5LData = append(y5LData, opts.LineData{Value: v5L})
	}
	line.SetXAxis(x).AddSeries(title, yData).
		AddSeries("p99", y9Data).
		AddSeries("p95", y5Data).
		AddSeries("p99L", y9LData).
		AddSeries("p95L", y5LData).
		SetSeriesOptions(charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}))
	return line
}

func plottingCharts(SmoothedSnps []TwoBulkTwoParentsRecord, outputHTML string, outputCSV string, smoothing bool) error {
	fmt.Printf("Creating charts and Detecting QTLs ...\n\n")

	// ------------------------------ QTL File -------------------------------------------//

	chromMap := make(map[string][]TwoBulkTwoParentsRecord)
	for _, snp := range SmoothedSnps {
		chromMap[snp.Chrom] = append(chromMap[snp.Chrom], snp)
	}

	chroms := make([]string, 0, len(chromMap))
	for chrom := range chromMap {
		chroms = append(chroms, chrom)
	}
	sort.Strings(chroms)

	page := components.NewPage()
	// -------------- QTL file -----------------//
	qtlFile, _ := os.Create(outputCSV)

	defer qtlFile.Close()

	_, err := fmt.Fprintf(qtlFile, "CHROM\tQTLstart\tQTLend\tPEAK\tSTATISTIC\tTHRESHOLD\n")
	if err != nil {
		return err
	}

	//var qtl_data []TwoBulkQtlRecord
	for _, chrom := range chroms {
		chromSNPs := chromMap[chrom]
		var x []int
		var HIy, HIp99y, HIp95y, HIp99Ly, HIp95Ly, LIy, LIp99y, LIp95y, LIp99Ly, LIp95Ly, DSy, DSp99y, DSp95y, DSp99Ly, DSp95Ly, GSy, GSp99y, GSp95y []float64
		for _, snp := range chromSNPs {
			x = append(x, int(snp.Pos))
			HIy = append(HIy, snp.HighSI)
			HIp99y = append(HIp99y, snp.HighP99)
			HIp95y = append(HIp95y, snp.HighP95)
			HIp99Ly = append(HIp99Ly, snp.HighMp99)
			HIp95Ly = append(HIp95Ly, snp.HighMp95)

			LIy = append(LIy, snp.LowSI)
			LIp99y = append(LIp99y, snp.LowP99)
			LIp95y = append(LIp95y, snp.LowP95)
			LIp99Ly = append(LIp99Ly, snp.LowMp99)
			LIp95Ly = append(LIp95Ly, snp.LowMp95)

			DSy = append(DSy, snp.DSI)
			DSp99y = append(DSp99y, snp.DsiP99)
			DSp95y = append(DSp95y, snp.DsiP95)
			DSp99Ly = append(DSp99Ly, snp.DsiMp99)
			DSp95Ly = append(DSp95Ly, snp.DsiMp95)

			GSy = append(GSy, snp.GS)
			GSp99y = append(GSp99y, snp.GsP99)
			GSp95y = append(GSp95y, snp.GsP95)

		}

		//fmt.Printf("Creating charts for chromosome: %s ...\n", chrom)
		deltaChart := createLineChart(x, DSy, DSp99y, DSp95y, DSp99Ly, DSp95Ly, chrom+" ΔSNP-index", "ΔSNP-index", smoothing)
		hiChart := createLineChart(x, HIy, HIp99y, HIp95y, HIp99Ly, HIp95Ly, chrom+" High-SNP-index", "High-SNP-index", smoothing)
		liChart := createLineChart(x, LIy, LIp99y, LIp95y, LIp99Ly, LIp95Ly, chrom+" Low-SNP-index", "Low-SNP-index", smoothing)
		gChart := createLineChart(x, GSy, GSp99y, GSp95y, nil, nil, chrom+" G-Statistics", "G-Statistics", smoothing)

		_, ds9peakY, ds9start, ds9end, ds9Found := detectQtlPeaks(x, DSy, DSp99y)
		_, ds5peakY, ds5start, ds5end, ds5Found := detectQtlPeaks(x, DSy, DSp95y)
		_, ds9LpeakY, ds9Lstart, ds9Lend, ds9LFound := detectQtlValleys(x, DSy, DSp99Ly)
		_, ds5LpeakY, ds5Lstart, ds5Lend, ds5LFound := detectQtlValleys(x, DSy, DSp95Ly)

		_, gs9peakY, gs9start, gs9end, gs9Found := detectQtlPeaks(x, GSy, GSp99y)
		_, gs5peakY, gs5start, gs5end, gs5Found := detectQtlPeaks(x, GSy, GSp95y)

		_, hi9peakY, hi9start, hi9end, hi9Found := detectQtlPeaks(x, HIy, HIp99y)
		_, hi5peakY, hi5start, hi5end, hi5Found := detectQtlPeaks(x, HIy, HIp95y)
		_, hi9LpeakY, hi9Lstart, hi9Lend, hi9LFound := detectQtlValleys(x, HIy, HIp99Ly)
		_, hi5LpeakY, hi5Lstart, hi5Lend, hi5LFound := detectQtlValleys(x, HIy, HIp95Ly)

		_, li9peakY, li9start, li9end, li9Found := detectQtlPeaks(x, LIy, LIp99y)
		_, li5peakY, li5start, li5end, li5Found := detectQtlPeaks(x, LIy, LIp95y)
		_, li9LpeakY, li9Lstart, li9Lend, li9LFound := detectQtlValleys(x, LIy, LIp99Ly)
		_, li5LpeakY, li5Lstart, li5Lend, li5LFound := detectQtlValleys(x, LIy, LIp95Ly)

		if ds9Found {
			_, dErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tΔSNP-index\t99per (peak)\n", chrom, ds9start, ds9end, ds9peakY)
			if dErr != nil {
				return dErr
			}
			//fmt.Printf("%s, %v, %v, %v, ΔSNP-index, 99per (peak)", chrom, ds9start, ds9end, ds9peakY)
		}

		if ds5Found {
			_, d5Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tΔSNP-index\t95per (peak)\n", chrom, ds5start, ds5end, ds5peakY)
			if d5Err != nil {
				return d5Err
			}
			//fmt.Printf("%s, %v, %v, %v ΔSNP-index (95perc confidene)", chrom, ds5start5, ds5end5, ds5peakY)
		}

		if ds9LFound {
			_, ds9lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tΔSNP-index\t99per (valley)\n", chrom, ds9Lstart, ds9Lend, ds9LpeakY)
			if ds9lErr != nil {
				return ds9lErr
			}
			//fmt.Printf("%s, %v, %v, %v, ΔSNP-index, 99per (valley)", chrom, ds5Lstart, ds5Lend, ds5LpeakY)
		}

		if ds5LFound {
			_, d5lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tΔSNP-index\t95per (valley)\n", chrom, ds5Lstart, ds5Lend, ds5LpeakY)
			if d5lErr != nil {
				return d5lErr
			}
			//fmt.Printf("%s, %v, %v, %v ΔSNP-index (95perc confidene)", chrom, ds5start5, ds5end5, ds5peakY)
		}

		if gs9Found {
			_, gs9Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tG-statistic\t99per (peak)\n", chrom, gs9start, gs9end, gs9peakY)
			if gs9Err != nil {
				return gs9Err
			}
			//fmt.Printf("%s, %v, %v, %v G-statistic (99perc peak)", chrom, gs9start, gs9end, gs9peakY)
		}

		if gs5Found {
			_, d5Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tG-statistic\t95per (peak)\n", chrom, gs5start, gs5end, gs5peakY)
			if d5Err != nil {
				return d5Err
			}
			//fmt.Printf("%s, %v, %v, %v G-statistic (95perc confidene)", chrom, gs5start, gs5end, gs5peakY)
		}

		if hi9Found {
			_, hi9Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tHigh-SNP-index\t99per (peak)\n", chrom, hi9start, hi9end, hi9peakY)
			if hi9Err != nil {
				return hi9Err
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (99perc confidene)", chrom, hi9start, hi9end, hi9peakY)
		}
		if hi5Found {
			_, hi5Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tHigh-SNP-index\t95per (peak)\n", chrom, hi5start, hi5end, hi5peakY)
			if hi5Err != nil {
				return hi5Err
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (95perc confidene)", chrom, hi5start, hi5end, hi5peakY)
		}

		if hi9LFound {
			_, hi9lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tHi-SNP-index\t99per (valley)\n", chrom, hi9Lstart, hi9Lend, hi9LpeakY)
			if hi9lErr != nil {
				return hi9lErr
			}

		}
		if hi5LFound {
			_, hi5lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tHigh-SNP-index\t95per (valley)\n", chrom, hi5Lstart, hi5Lend, hi5LpeakY)
			if hi5lErr != nil {
				return hi5lErr
			}
			//fmt.Printf("%s, %v, %v, %v High-SNP-index (99perc valley)", chrom, hi9Lstart, hi9Lend, hi9LpeakY)
		}
		if li9Found {
			_, li9Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tLow-SNP-index\t99per (peak)\n", chrom, li9start, li9end, li9peakY)
			if li9Err != nil {
				return li9Err
			}
			//fmt.Printf("%s, %v, %v, %v Low-SNP-index (99perc confidene)", chrom, li9start, li9end, li9peakY)
		}
		if li5Found {
			_, li5Err := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tLow-SNP-index\t95per (peak)\n", chrom, li5start, li5end, li5peakY)
			if li5Err != nil {
				return li5Err
			}
			//fmt.Printf("%s, %v, %v, %v Low-SNP-index (95perc confidene)", chrom, li5start, li5end, li5peakY)
		}

		if li9LFound {
			_, li9lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tLow-SNP-index\t99per (valley)\n", chrom, li9Lstart, li9Lend, li9LpeakY)
			if li9lErr != nil {
				return li9lErr
			}
			//fmt.Printf("%s, %v, %v, %v Low-SNP-index (99perc confidene)", chrom, li9Lstart, li9Lend, li9LpeakY)
		}
		if li5LFound {
			_, li5lErr := fmt.Fprintf(qtlFile, "%s\t%v\t%v\t%v\tLow-SNP-index\t95per (valley)\n", chrom, li5Lstart, li5Lend, li5LpeakY)
			if li5lErr != nil {
				return li5lErr
			}
			//fmt.Printf("%s, %v, %v, %v Low-SNP-index (99perc confidene)", chrom, li5Lstart, li5Lend, li5LpeakY)
		}

		page.SetLayout(components.PageFlexLayout)
		page.AddCharts(hiChart, deltaChart, liChart, gChart)

	}

	f, err := os.Create(outputHTML)
	if err != nil {
		return err
	}
	return page.Render(f)
}

func TwoBulkTwoParentsRun(
	vcfFile string,
	highParent string,
	lowParent string,
	highBulk string,
	lowBulk string,
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
	outDir string,
) {

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

	jlog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
	slog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

	// ======================================== Create Results dir ================================================== //
	resultsDir, err := createResultsDir(outDir)
	if err != nil {
		jlog.Error("BSASEQ", "PROGRAM", "RESULTS_DIR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		slog.Error("BSASEQ", "PROGRAM", "RESULTS_DIR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		fmt.Println("Error creating results directory: ", err)
		return
	}

	// ========================================= Get tsv table ====================================================== //
	fmt.Printf("================================== Filtering Start ======================================\n\n")
	var tsvFile string
	var filteredRecords []TwoBulkTwoParentsRecord
	if strings.ToLower(filepath.Ext(vcfFile)) == ".vcf" || strings.ToLower(filepath.Ext(vcfFile)) == ".gz" {
		fmt.Printf("Working with VCF file ...\n\n")

		fmt.Println("Running gatk VariantsToTable ...")
		tsvFile = filepath.Join(resultsDir, "rawVariants.tsv")
		cmd3 := exec.Command("gatk", "VariantsToTable", "-V", vcfFile, "-F", "CHROM", "-F", "POS", "-F",
			"REF", "-F", "ALT", "-F", "QUAL", "-F", "TYPE", "-GF", "GT", "-GF", "AD", "-GF", "DP", "-GF", "GQ", "-O", tsvFile)

		varTabErr := cmd3.Run()
		if varTabErr != nil {
			log.Fatalf("error running gatk SelectVariants: %s", varTabErr)
		}

	} else {
		fmt.Println("Working with tsv file ...")
		tsvFile = vcfFile
	}

	if utils.StageHasCompleted(logged, "FILTERING", "ALL", "ALL") {
		fmt.Println("FILTERING has already completed. Skipping.")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

		fmt.Printf("Filtering DF ...\n\n")
		filteredRecords, err = twoBulkTwoParTsvFilter(tsvFile, highParent, minHighParentDepth, lowParent, minLowParentDepth,
			highBulk, minHighBulkDepth, lowBulk, minLowBulkDepth, windowSize, stepSize, resultsDir, logFile)

		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
		}
		jlog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		fmt.Printf("# Variants after filtering: %d\n\n", len(filteredRecords))
	}

	fmt.Printf("================================== Filtering End ======================================\n\n")

	// ============================================= STATISTICS ===================================================== //

	fmt.Printf("====================================== BSAseq Statistics Start ========================================== \n\n")
	statsFile := filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_stats.tsv")
	var statsRecords []TwoBulkTwoParentsRecord
	if utils.StageHasCompleted(logged, "STATS", "ALL", "ALL") {
		fmt.Println("STATS has already completed. Skipping.")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

		statsRecords, err = twoBulkTwoParStats(filteredRecords, highBulkSize, lowBulkSize, popStructure, rep)

		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}

		fmt.Printf("Writing stats file to %s... \n\n", statsFile)
		err = writeTwoBulkTwoPar(statsRecords, highParent, lowParent, highBulk, lowBulk, statsFile)
		if err != nil {
			jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

	}

	fmt.Printf("====================================== BSAseq Statistics End ========================================== \n\n")

	// ============================================= PLOTTING ======================================================= //
	fmt.Printf("====================================== BSAseq Plotting Start ========================================== \n\n")
	fmt.Printf("Performing sliding window analysis for SNP ...\n\n")
	var slidingRecordsSNPs []TwoBulkTwoParentsRecord
	var slidingRecordsINDELs []TwoBulkTwoParentsRecord

	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "SNP", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slidingRecordsSNPs = slidingWindowAnalysis(statsRecords, "SNP", windowSize, stepSize)
		if len(slidingRecordsSNPs) == 0 {
			jlog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_snp.tsv"))

		err = writeTwoBulkTwoPar(slidingRecordsSNPs, highParent, lowParent, highBulk, lowBulk, filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_snp.tsv"))
		if err != nil {
			fmt.Println("Error writing stats file to file")
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

	}

	if utils.StageHasCompleted(logged, "PLOTTING", "SNP", "ALL") {
		fmt.Println("Plotting is SNPs already done. Skipping ...")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "STARTED")

		err = plottingCharts(slidingRecordsSNPs, filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_plot_snp.html"), filepath.Join(resultsDir, "goBSAseq_snp.tsv"), smoothing)
		if err != nil {
			fmt.Println("Error plotting charts: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		fmt.Printf("====================================== BSAseq Plotting End ========================================== \n\n")

	}

	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "INDEL", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slidingRecordsINDELs = slidingWindowAnalysis(statsRecords, "INDEL", windowSize, stepSize)
		if len(slidingRecordsINDELs) == 0 {
			jlog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_indel.tsv"))

		err = writeTwoBulkTwoPar(slidingRecordsINDELs, highParent, lowParent, highBulk, lowBulk, filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_indel.tsv"))
		if err != nil {
			fmt.Println("Error writing stats file to file")
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

	}

	if utils.StageHasCompleted(logged, "PLOTTING", "INDEL", "INDEL") {
		fmt.Println("Plotting SNPs is already done. Skipping ...")
		return
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

		err = plottingCharts(slidingRecordsINDELs, filepath.Join(resultsDir, highParent+"_samp_"+lowParent+"_samp_"+highBulk+"_samp_"+lowBulk+"_both_bsaseq_plot.html"), filepath.Join(resultsDir, "goBSAseq.tsv"), smoothing)
		if err != nil {
			fmt.Println("Error plotting charts: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", statsFile, "CHROMOSOME", "ALL", "STATUS", "FAILED")
			slog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", statsFile, "CHROMOSOME", "ALL", "STATUS", "FAILED")
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", statsFile, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", statsFile, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		fmt.Printf("====================================== BSAseq Plotting End ========================================== \n\n")

	}

}

func TwoBulkOnlyRun(
	vcfFile string,
	highBulk string,
	lowBulk string,
	minHighBulkDepth int,
	minLowBulkDepth int,
	highBulkSize int,
	lowBulkSize int,
	windowSize int,
	stepSize int,
	smoothing bool,
	popStructure string,
	rep int,
	outDir string,
) {
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

	jlog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
	slog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

	// ======================================== Create Results dir =================================================== #
	resultsDir, err := createResultsDir(outDir)
	if err != nil {
		fmt.Println("Error creating results directory: ", err)
		return
	}

	// ======================================== Filter vcf file ====================================================== #
	fmt.Printf("================================== Filtering Start ======================================\n\n")
	var filteredRecords []TwoBulkTwoParentsRecord
	var statsRecords []TwoBulkTwoParentsRecord
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
		fmt.Printf("Filtering DF ...\n\n")
		tsvFile = vcfFile

	}

	if utils.StageHasCompleted(logged, "FILTERING", "ALL", "ALL") {
		fmt.Println("FILTERING has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

		filteredRecords, err = twoBulksOnlyFilter(tsvFile, highBulk, minHighBulkDepth, lowBulk, minLowBulkDepth, windowSize, stepSize, resultsDir)
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "FILTERING", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		fmt.Printf("# Variants after filtering: %d\n\n", len(filteredRecords))
		jlog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BSASEQ", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		fmt.Printf("================================== Filtering End ======================================\n\n")

	}

	// ============================================= STATISTICS ====================================================== #

	fmt.Printf("====================================== BSAseq Statistics Start ========================================== \n\n")

	if utils.StageHasCompleted(logged, "STATS", "ALL", "ALL") {
		fmt.Println("STATS has already completed. Skipping.")

	} else {
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		statsRecords, err = twoBulkOnlyStats(filteredRecords, highBulkSize, lowBulkSize, popStructure, rep)
		if err != nil {
			jlog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		statsFile := filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_stats.tsv")

		fmt.Printf("Writing stats file to %s... \n\n", statsFile)
		err = writeTwoBulkOnly(statsRecords, highBulk, lowBulk, statsFile)
		if err != nil {
			fmt.Println("Error writing stats file to file")
			jlog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "STATS", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

		fmt.Printf("====================================== BSAseq Statistics End ========================================== \n\n")

	}

	// ============================================= PLOTTING ====================================================== #
	var slidingRecordsSNPs []TwoBulkTwoParentsRecord
	var slidingRecordsINDELs []TwoBulkTwoParentsRecord
	fmt.Printf("====================================== BSAseq Plotting Start ========================================== \n\n")
	fmt.Printf("Performing sliding window analysis for SNP ...\n\n")

	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "SNP", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

		slidingRecordsSNPs = slidingWindowAnalysis(statsRecords, "SNP", windowSize, stepSize)

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats.tsv"))
		err = writeTwoBulkOnly(slidingRecordsSNPs, highBulk, lowBulk, filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_SNP.tsv"))
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

		err = plottingCharts(slidingRecordsSNPs, filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_plot_SNP.html"), filepath.Join(resultsDir, "goBSAseq_SNP.tsv"), smoothing)
		if err != nil {
			fmt.Println("Error plotting charts: ", err)
			jlog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BSASEQ", "PROGRAM", "PLOTTING", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "SNP", "CHROMOSOME", "ALL", "STATUS", "COMPLETED", "CMD", "ALL")

	}

	if utils.StageHasCompleted(logged, "SLIDING_WINDOW", "INDEL", "ALL") {
		fmt.Println("SLIDING_WINDOW has already completed. Skipping.")
	} else {
		jlog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
		slog.Info("BSASEQ", "PROGRAM", "SLIDING_WINDOW", "SAMPLE", "INDEL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

		fmt.Printf("Performing sliding window analysis for INDELs ...\n\n")
		slidingRecordsIndels := slidingWindowAnalysis(statsRecords, "INDEL", windowSize, stepSize)

		fmt.Printf("Writing sliding window analysis file to %s... \n\n", filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_INDELS.tsv"))
		err = writeTwoBulkOnly(slidingRecordsIndels, highBulk, lowBulk, filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_sliding_stats_INDELS.tsv"))
		if err != nil {
			fmt.Println("Error writing stats file to file")
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

		err = plottingCharts(slidingRecordsINDELs, filepath.Join(resultsDir, highBulk+"_samp_"+lowBulk+"_both_bsaseq_plot_indel.html"), filepath.Join(resultsDir, "goBSAseq_indel.tsv"), smoothing)
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
