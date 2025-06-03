package bsaseq

import (
	"encoding/csv"
	"fmt"
	"golang.org/x/exp/rand"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

var thresholdCacheOne = make(map[string]map[string]float64)
var thresholdMuOne sync.Mutex

type OneBulkTwoParentsRecord struct {
	Chrom     string
	Pos       float64
	Ref       string
	Alt       string
	Type      string
	HighParGT string
	LowParGT  string
	BulkGT    string
	HighParDP int
	LowParDP  int
	BulkDP    int
	HighParAD string
	LowParAD  string
	BulkAD    string
	SI        float64
	P99       float64
	P95       float64
	Mp99      float64
	Mp95      float64
}

func readTsvToStructOneBulkTwoPar(tsvFile string, highPar string, lowPar string, bulk string) []OneBulkTwoParentsRecord {
	file, err := os.Open(tsvFile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.TrimLeadingSpace = true

	records, err := reader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}
	if len(records) < 2 {
		log.Fatalf("Expected at least 2 rows, got %d", len(records))
	}

	header := records[0]
	colIndex := make(map[string]int)

	requiredCols := []string{highPar + ".GT", lowPar + ".GT", bulk + ".GT", highPar + ".DP", lowPar + ".DP", bulk + ".DP", highPar + ".AD", lowPar + ".AD", bulk + ".AD"}
	for _, col := range requiredCols {
		found := false
		for i, headerCol := range header {
			if headerCol == col {
				found = true
				colIndex[headerCol] = i
				break
			}
		}
		if !found {
			log.Fatalf("required column %s not found in header", col)
		}
	}

	for i, col := range header {
		colIndex[col] = i
	}

	var data []OneBulkTwoParentsRecord
	for _, row := range records[1:] {
		pos, _ := strconv.ParseFloat(row[colIndex["POS"]], 64)
		hPdp, _ := strconv.Atoi(row[colIndex[highPar+".DP"]])
		lPdp, _ := strconv.Atoi(row[colIndex[lowPar+".DP"]])
		bDP, _ := strconv.Atoi(row[colIndex[bulk+".DP"]])

		r := OneBulkTwoParentsRecord{
			Chrom: row[colIndex["CHROM"]],
			Pos:   pos,
			Ref:   row[colIndex["REF"]],
			Alt:   row[colIndex["ALT"]],
			Type:  row[colIndex["TYPE"]],

			HighParGT: strings.ReplaceAll(row[colIndex[highPar+".GT"]], "|", "/"),
			HighParDP: hPdp,
			HighParAD: row[colIndex[highPar+".AD"]],

			LowParGT: strings.ReplaceAll(row[colIndex[lowPar+".GT"]], "|", "/"),
			LowParDP: lPdp,
			LowParAD: row[colIndex[lowPar+".AD"]],

			BulkGT: strings.ReplaceAll(row[colIndex[bulk+".GT"]], "|", "/"),
			BulkDP: bDP,
			BulkAD: row[colIndex[bulk+".AD"]],
		}
		data = append(data, r)
	}
	return data
}

func removeShortContigsOneBulkTwoPar(variants []OneBulkTwoParentsRecord, winSize int, stepSize int) []OneBulkTwoParentsRecord {
	chromMax := make(map[string]float64)
	for _, v := range variants {
		if v.Pos > chromMax[v.Chrom] {
			chromMax[v.Chrom] = v.Pos
		}
	}

	var goodChroms = make(map[string]bool)
	var chroms []string
	for chrom, maxPos := range chromMax {
		if int(maxPos) > winSize+stepSize {
			goodChroms[chrom] = true
			chroms = append(chroms, chrom)
		}
	}

	sort.Strings(chroms)

	fmt.Printf("Chromosomes suitable for BSA-seq analysis: \n%s\n\n", chroms)

	var bsaSeqRecords []OneBulkTwoParentsRecord
	for _, v := range variants {
		if goodChroms[v.Chrom] {
			bsaSeqRecords = append(bsaSeqRecords, v)
		}
	}

	return bsaSeqRecords
}

func writeOneBulkTwoPar(variants []OneBulkTwoParentsRecord, highPar string, lowPar string, bulk string, outputFile string) {
	file, err := os.Create(outputFile)
	if err != nil {
		log.Fatalf("Failed to create output CSV: %v", err)
	}
	defer func(file *os.File) {
		err := file.Close()
		if err != nil {
			panic(err)
		}
	}(file)

	writer := csv.NewWriter(file)
	writer.Comma = '\t'

	defer writer.Flush()

	header := []string{
		"CHROM", "POS", "REF", "ALT", "TYPE",
		lowPar + ".GT", lowPar + ".AD", lowPar + ".DP",
		highPar + ".GT", highPar + ".AD", highPar + ".DP",
		bulk + ".GT", bulk + ".AD", bulk + ".DP",
		bulk + "_SNP_INDEX", bulk + "_p99", bulk + "_p95", bulk + "_m_p99", bulk + "_m_p95",
	}
	if hErr := writer.Write(header); hErr != nil {
		log.Fatalf("Failed to write header: %v", hErr)
	}

	for _, rec := range variants {
		record := []string{
			rec.Chrom, strconv.FormatInt(int64(rec.Pos), 10), rec.Ref, rec.Alt, rec.Type,
			rec.LowParGT, rec.LowParAD, strconv.FormatInt(int64(rec.LowParDP), 10),
			rec.HighParGT, rec.HighParAD, strconv.FormatInt(int64(rec.HighParDP), 10),
			rec.BulkGT, rec.BulkAD, strconv.FormatInt(int64(rec.BulkDP), 10),

			strconv.FormatFloat(rec.SI, 'f', 6, 64),
			strconv.FormatFloat(rec.P99, 'f', 6, 64), strconv.FormatFloat(rec.P95, 'f', 6, 64),
			strconv.FormatFloat(rec.Mp99, 'f', 6, 64), strconv.FormatFloat(rec.Mp95, 'f', 6, 64),
		}
		if wErr := writer.Write(record); wErr != nil {
			log.Fatalf("Error: %s\n", wErr)
		}
	}

}

func oneBulkThresholdsCached(bulkDP float64, smAF float64, rep int) map[string]float64 {
	key := fmt.Sprintf("%f", bulkDP)
	thresholdMuOne.Lock()
	if thresholds, ok := thresholdCacheOne[key]; ok {
		thresholdMuOne.Unlock()
		return thresholds
	}
	thresholdMuOne.Unlock()

	thresholds := oneBulkThresholds(int(bulkDP), smAF, rep)
	thresholdMuOne.Lock()
	thresholdCacheOne[key] = thresholds
	thresholdMuOne.Unlock()
	return thresholds
}

func oneBulkThresholds(bulkDP int, smAF float64, rep int) map[string]float64 {
	src := rand.NewSource(uint64(time.Now().UnixNano()))
	binomial := distuv.Binomial{N: float64(bulkDP), P: smAF, Src: src}
	smBulkAFArray := make([]float64, rep)

	for i := 0; i < rep; i++ {
		smBulkAltDP := binomial.Rand()
		smBulkAFArray[i] = math.Round((smBulkAltDP/float64(bulkDP))*1e6) / 1e6
	}
	sort.Float64s(smBulkAFArray)

	BulkAf99 := stat.Quantile(0.995, stat.Empirical, smBulkAFArray, nil)
	BulkAf95 := stat.Quantile(0.95, stat.Empirical, smBulkAFArray, nil)
	BulkAfL99 := stat.Quantile(0.005, stat.Empirical, smBulkAFArray, nil)
	BulkAfL95 := stat.Quantile(0.05, stat.Empirical, smBulkAFArray, nil)

	thresholds := map[string]float64{"bulk99": math.Round(BulkAf99*1e6) / 1e6, "bulk95": math.Round(BulkAf95*1e6) / 1e6,
		"bulkL99": math.Round(BulkAfL99*1e6) / 1e6, "bulkL95": math.Round(BulkAfL95*1e6) / 1e6}

	return thresholds

}

func calculateStatsRecordOne(rec OneBulkTwoParentsRecord, rep int, smAF float64) OneBulkTwoParentsRecord {
	lGT := rec.LowParGT
	ref := rec.Ref
	BDP := rec.BulkDP
	BAD := rec.BulkAD

	lowParGT := strings.Split(lGT, "/")[0]
	bulkRefAltAD := strings.Split(BAD, ",")

	if ref == lowParGT {
		bulkAltAD, _ := strconv.ParseFloat(bulkRefAltAD[1], 64)
		rec.SI = math.Round(bulkAltAD/float64(BDP)*1e6) / 1e6
		threshMap := oneBulkThresholdsCached(float64(BDP), smAF, rep)
		rec.P99 = threshMap["bulk99"]
		rec.P95 = threshMap["bulk95"]
		rec.Mp99 = threshMap["bulkL99"]
		rec.Mp95 = threshMap["bulkL95"]
	} else {
		bulkAltAD, _ := strconv.ParseFloat(bulkRefAltAD[0], 64)
		rec.SI = math.Round(bulkAltAD/float64(BDP)*1e6) / 1e6
		threshMap := oneBulkThresholdsCached(float64(BDP), smAF, rep)
		rec.P99 = threshMap["bulk99"]
		rec.P95 = threshMap["bulk95"]
		rec.Mp99 = threshMap["bulkL99"]
		rec.Mp95 = threshMap["bulkL95"]
	}

	return rec
}
