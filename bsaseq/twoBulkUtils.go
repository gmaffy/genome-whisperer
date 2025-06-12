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

var thresholdCache = make(map[string]map[string]float64)
var thresholdMu sync.Mutex

type TwoBulkTwoParentsRecord struct {
	Chrom      string
	Pos        float64
	Ref        string
	Alt        string
	Type       string
	HighParGT  string
	LowParGT   string
	HighBulkGT string
	LowBulkGT  string
	HighParDP  int
	LowParDP   int
	HighBulkDP int
	LowBulkDP  int
	HighParAD  string
	LowParAD   string
	HighBulkAD string
	LowBulkAD  string
	HighSI     float64
	LowSI      float64
	DSI        float64
	GS         float64
	HighP99    float64
	HighP95    float64
	HighMp99   float64
	HighMp95   float64
	LowP99     float64
	LowP95     float64
	LowMp99    float64
	LowMp95    float64
	DsiP99     float64
	DsiP95     float64
	DsiMp99    float64
	DsiMp95    float64
	GsP99      float64
	GsP95      float64
}

type TwoBulkOnlyRecord struct {
	Chrom      string
	Pos        float64
	Ref        string
	Alt        string
	Type       string
	HighBulkGT string
	LowBulkGT  string
	HighBulkDP int
	LowBulkDP  int
	HighBulkAD string
	LowBulkAD  string
	HighSI     float64
	LowSI      float64
	DSI        float64
	GS         float64
	HighP99    float64
	HighP95    float64
	HighMp99   float64
	HighMp95   float64
	LowP99     float64
	LowP95     float64
	LowMp99    float64
	LowMp95    float64
	DsiP99     float64
	DsiP95     float64
	DsiMp99    float64
	DsiMp95    float64
	GsP99      float64
	GsP95      float64
}
type TwoBulkQtlRecord struct {
	Chrom string
	HIp99 string

	HIp95  string
	HIp99L string
	HIp95L string
	Lip99  string
	Lip95  string
	Lip99L string
	Lip95L string
	DSp99  string
	DSp95  string
	DSp99L string
	DSp95L string
	GSp99  string
	GSp95  string
}

func readTsvToStructTwoBulkTwoPar(tsvFile string, highPar string, lowPar string, highBulk string, lowBulk string) []TwoBulkTwoParentsRecord {
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

	requiredCols := []string{highPar + ".GT", lowPar + ".GT", highBulk + ".GT", lowBulk + ".GT", highPar + ".DP", lowPar + ".DP", highBulk + ".DP", lowBulk + ".DP", highPar + ".AD", lowPar + ".AD", highBulk + ".AD", lowBulk + ".AD"}
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

	var data []TwoBulkTwoParentsRecord
	for _, row := range records[1:] {
		pos, _ := strconv.ParseFloat(row[colIndex["POS"]], 64)
		hPdp, _ := strconv.Atoi(row[colIndex[highPar+".DP"]])
		lPdp, _ := strconv.Atoi(row[colIndex[lowPar+".DP"]])
		hBdp, _ := strconv.Atoi(row[colIndex[highBulk+".DP"]])
		lBdp, _ := strconv.Atoi(row[colIndex[lowBulk+".DP"]])

		r := TwoBulkTwoParentsRecord{
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

			HighBulkGT: strings.ReplaceAll(row[colIndex[highBulk+".GT"]], "|", "/"),
			HighBulkDP: hBdp,
			HighBulkAD: row[colIndex[highBulk+".AD"]],

			LowBulkGT: strings.ReplaceAll(row[colIndex[lowBulk+".GT"]], "|", "/"),
			LowBulkDP: lBdp,
			LowBulkAD: row[colIndex[lowBulk+".AD"]],
		}
		data = append(data, r)
	}
	return data
}

func readTsvToStructTwoBulkOnly(tsvFile string, highBulk string, lowBulk string) []TwoBulkOnlyRecord {
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

	requiredCols := []string{highBulk + ".GT", lowBulk + ".GT", highBulk + ".DP", lowBulk + ".DP", highBulk + ".AD", lowBulk + ".AD"}
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

	var data []TwoBulkOnlyRecord
	for _, row := range records[1:] {
		pos, _ := strconv.ParseFloat(row[colIndex["POS"]], 64)

		hBdp, _ := strconv.Atoi(row[colIndex[highBulk+".DP"]])
		lBdp, _ := strconv.Atoi(row[colIndex[lowBulk+".DP"]])

		r := TwoBulkOnlyRecord{
			Chrom: row[colIndex["CHROM"]],
			Pos:   pos,
			Ref:   row[colIndex["REF"]],
			Alt:   row[colIndex["ALT"]],
			Type:  row[colIndex["TYPE"]],

			HighBulkGT: strings.ReplaceAll(row[colIndex[highBulk+".GT"]], "|", "/"),
			HighBulkDP: hBdp,
			HighBulkAD: row[colIndex[highBulk+".AD"]],

			LowBulkGT: strings.ReplaceAll(row[colIndex[lowBulk+".GT"]], "|", "/"),
			LowBulkDP: lBdp,
			LowBulkAD: row[colIndex[lowBulk+".AD"]],
		}
		data = append(data, r)
	}
	return data
}

func gStatistic(a1, b1, a2, b2 int) float64 {
	n1 := a1 + b1
	n2 := a2 + b2
	if n1 == 0 || n2 == 0 {
		return 0.0
	}

	totalRef := a1 + a2
	totalDepth := n1 + n2
	p := float64(totalRef) / float64(totalDepth)

	g := 0.0

	if a1 > 0 {
		g += float64(a1) * math.Log(float64(a1)/(float64(n1)*p))
	}
	if b1 > 0 {
		g += float64(b1) * math.Log(float64(b1)/(float64(n1)*(1-p)))
	}
	if a2 > 0 {
		g += float64(a2) * math.Log(float64(a2)/(float64(n2)*p))
	}
	if b2 > 0 {
		g += float64(b2) * math.Log(float64(b2)/(float64(n2)*(1-p)))
	}

	return 2 * g
}

func twoBulkThresholdsCached(highBulkDP float64, lowBulkDP float64, resSmAF float64, susSmAF float64, rep int) map[string]float64 {
	key := fmt.Sprintf("%f_%f", highBulkDP, lowBulkDP)
	thresholdMu.Lock()
	if thresholds, ok := thresholdCache[key]; ok {
		thresholdMu.Unlock()
		return thresholds
	}
	thresholdMu.Unlock()

	thresholds := twoBulkThresholds(int(highBulkDP), int(lowBulkDP), resSmAF, susSmAF, rep)
	thresholdMu.Lock()
	thresholdCache[key] = thresholds
	thresholdMu.Unlock()
	return thresholds
}

func twoBulkThresholds(highBulkDP int, lowBulkDP int, highSmAF float64, lowSmAF float64, rep int) map[string]float64 {
	src := rand.NewSource(uint64(time.Now().UnixNano()))
	highBinomial := distuv.Binomial{N: float64(highBulkDP), P: highSmAF, Src: src}
	lowBinomial := distuv.Binomial{N: float64(lowBulkDP), P: lowSmAF, Src: src}

	smhighBulkAFArray := make([]float64, rep)
	smlowBulkAFArray := make([]float64, rep)
	smDsiArray := make([]float64, rep)
	smGsArray := make([]float64, rep)
	for i := 0; i < rep; i++ {
		smResAltDP := highBinomial.Rand()
		smSusAltDP := lowBinomial.Rand()
		smResRefDP := float64(highBulkDP) - smResAltDP
		smSusRefDP := float64(lowBulkDP) - smSusAltDP
		smhighBulkAFArray[i] = math.Round((smResAltDP/float64(highBulkDP))*1e6) / 1e6
		smlowBulkAFArray[i] = math.Round((smSusAltDP/float64(lowBulkDP))*1e6) / 1e6
		smDsiArray[i] = math.Round((smhighBulkAFArray[i]-smlowBulkAFArray[i])*1e6) / 1e6
		smGsArray[i] = math.Round(gStatistic(int(smResAltDP), int(smResRefDP), int(smSusAltDP), int(smSusRefDP))*1e6) / 1e6

	}
	sort.Float64s(smlowBulkAFArray)
	sort.Float64s(smhighBulkAFArray)
	sort.Float64s(smDsiArray)
	sort.Float64s(smGsArray)

	highBulkAf99 := stat.Quantile(0.995, stat.Empirical, smhighBulkAFArray, nil)
	highBulkAf95 := stat.Quantile(0.95, stat.Empirical, smhighBulkAFArray, nil)
	highBulkAfL99 := stat.Quantile(0.005, stat.Empirical, smhighBulkAFArray, nil)
	highBulkAfL95 := stat.Quantile(0.05, stat.Empirical, smhighBulkAFArray, nil)

	lowBulkAf99 := stat.Quantile(0.995, stat.Empirical, smlowBulkAFArray, nil)
	lowBulkAf95 := stat.Quantile(0.95, stat.Empirical, smlowBulkAFArray, nil)
	lowBulkAfL99 := stat.Quantile(0.005, stat.Empirical, smlowBulkAFArray, nil)
	lowBulkAfL95 := stat.Quantile(0.05, stat.Empirical, smlowBulkAFArray, nil)

	dsiAf99 := stat.Quantile(0.995, stat.Empirical, smDsiArray, nil)
	dsiAf95 := stat.Quantile(0.95, stat.Empirical, smDsiArray, nil)
	dsiAfL99 := stat.Quantile(0.005, stat.Empirical, smDsiArray, nil)
	dsiAfL95 := stat.Quantile(0.05, stat.Empirical, smDsiArray, nil)

	gs99 := stat.Quantile(0.995, stat.Empirical, smGsArray, nil)
	gs95 := stat.Quantile(0.95, stat.Empirical, smGsArray, nil)

	thresholds := map[string]float64{"res99": math.Round(highBulkAf99*1e6) / 1e6, "res95": math.Round(highBulkAf95*1e6) / 1e6,
		"resL99": math.Round(highBulkAfL99*1e6) / 1e6, "resL95": math.Round(highBulkAfL95*1e6) / 1e6,
		"sus99": math.Round(lowBulkAf99*1e6) / 1e6, "sus95": math.Round(lowBulkAf95*1e6) / 1e6,
		"susL99": math.Round(lowBulkAfL99*1e6) / 1e6, "susL95": math.Round(lowBulkAfL95*1e6) / 1e6,
		"dsi99": math.Round(dsiAf99*1e6) / 1e6, "dsi95": math.Round(dsiAf95*1e6) / 1e6, "dsiL99": math.Round(dsiAfL99*1e6) / 1e6,
		"dsiL95": math.Round(dsiAfL95*1e6) / 1e6, "gs99": math.Round(gs99*1e6) / 1e6, "gs95": math.Round(gs95*1e6) / 1e6}

	return thresholds

}

func removeShortContigsTwoBulkTwoPar(variants []TwoBulkTwoParentsRecord, winSize int, stepSize int) []TwoBulkTwoParentsRecord {
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

	var bsaSeqRecords []TwoBulkTwoParentsRecord
	for _, v := range variants {
		if goodChroms[v.Chrom] {
			bsaSeqRecords = append(bsaSeqRecords, v)
		}
	}

	return bsaSeqRecords
}

func removeShortContigsTwoBulkOnly(variants []TwoBulkOnlyRecord, winSize int, stepSize int) []TwoBulkOnlyRecord {
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

	var bsaSeqRecords []TwoBulkOnlyRecord
	for _, v := range variants {
		if goodChroms[v.Chrom] {
			bsaSeqRecords = append(bsaSeqRecords, v)
		}
	}

	return bsaSeqRecords
}

func writeTwoBulkTwoPar(variants []TwoBulkTwoParentsRecord, highPar string, lowPar string, highBulk string, lowBulk string, outputFile string) {
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

	// Write header
	header := []string{
		"CHROM", "POS", "REF", "ALT", "TYPE",
		lowPar + ".GT", lowPar + ".AD", lowPar + ".DP",
		highPar + ".GT", highPar + ".AD", highPar + ".DP",
		lowBulk + ".GT", lowBulk + ".AD", lowBulk + ".DP",
		highBulk + ".GT", highBulk + ".AD", highBulk + ".DP",
		highBulk + "_SNP_INDEX", lowBulk + "_SNP_INDEX", "DELTA_SNP_INDEX", "G_STATISTIC",
		highBulk + "_p99", highBulk + "_p95", highBulk + "_m_p99", highBulk + "_m_p95",
		lowBulk + "_p99", lowBulk + "_p95", lowBulk + "_m_p99", lowBulk + "_m_p95",
		"dsi_p99", "dsi_p95", "dsi_m_p99", "dsi_m_p95",
		"gs_p99", "gs_p95",
	}
	if hErr := writer.Write(header); hErr != nil {
		log.Fatalf("Failed to write header: %v", hErr)
	}

	// Write records
	for _, rec := range variants {
		record := []string{
			rec.Chrom, strconv.FormatInt(int64(rec.Pos), 10), rec.Ref, rec.Alt, rec.Type,
			rec.LowParGT, rec.LowParAD, strconv.FormatInt(int64(rec.LowParDP), 10),
			rec.HighParGT, rec.HighParAD, strconv.FormatInt(int64(rec.HighParDP), 10),
			rec.LowBulkGT, rec.LowBulkAD, strconv.FormatInt(int64(rec.LowBulkDP), 10),
			rec.HighBulkGT, rec.HighBulkAD, strconv.FormatInt(int64(rec.HighBulkDP), 10),

			strconv.FormatFloat(rec.HighSI, 'f', 6, 64), strconv.FormatFloat(rec.LowSI, 'f', 6, 64),
			strconv.FormatFloat(rec.DSI, 'f', 6, 64), strconv.FormatFloat(rec.GS, 'f', 6, 64),
			strconv.FormatFloat(rec.HighP99, 'f', 6, 64), strconv.FormatFloat(rec.HighP95, 'f', 6, 64),
			strconv.FormatFloat(rec.HighMp99, 'f', 6, 64), strconv.FormatFloat(rec.HighMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.LowP99, 'f', 6, 64), strconv.FormatFloat(rec.LowP95, 'f', 6, 64),
			strconv.FormatFloat(rec.LowMp99, 'f', 6, 64), strconv.FormatFloat(rec.LowMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.DsiP99, 'f', 6, 64), strconv.FormatFloat(rec.DsiP95, 'f', 6, 64),
			strconv.FormatFloat(rec.DsiMp99, 'f', 6, 64), strconv.FormatFloat(rec.DsiMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.GsP99, 'f', 6, 64), strconv.FormatFloat(rec.GsP95, 'f', 6, 64),
		}
		if wErr := writer.Write(record); wErr != nil {
			log.Fatalf("Error: %s\n", wErr)
		}
	}

}

func writeTwoBulkOnly(variants []TwoBulkOnlyRecord, highBulk string, lowBulk string, outputFile string) {
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

	// Write header
	header := []string{
		"CHROM", "POS", "REF", "ALT", "TYPE",
		lowBulk + ".GT", lowBulk + ".AD", lowBulk + ".DP",
		highBulk + ".GT", highBulk + ".AD", highBulk + ".DP",
		highBulk + "_SNP_INDEX", lowBulk + "_SNP_INDEX", "DELTA_SNP_INDEX", "G_STATISTIC",
		highBulk + "_p99", highBulk + "_p95", highBulk + "_m_p99", highBulk + "_m_p95",
		lowBulk + "_p99", lowBulk + "_p95", lowBulk + "_m_p99", lowBulk + "_m_p95",
		"dsi_p99", "dsi_p95", "dsi_m_p99", "dsi_m_p95",
		"gs_p99", "gs_p95",
	}
	if hErr := writer.Write(header); hErr != nil {
		log.Fatalf("Failed to write header: %v", hErr)
	}

	// Write records
	for _, rec := range variants {
		record := []string{
			rec.Chrom, strconv.FormatInt(int64(rec.Pos), 10), rec.Ref, rec.Alt, rec.Type,
			rec.LowBulkGT, rec.LowBulkAD, strconv.FormatInt(int64(rec.LowBulkDP), 10),
			rec.HighBulkGT, rec.HighBulkAD, strconv.FormatInt(int64(rec.HighBulkDP), 10),

			strconv.FormatFloat(rec.HighSI, 'f', 6, 64), strconv.FormatFloat(rec.LowSI, 'f', 6, 64),
			strconv.FormatFloat(rec.DSI, 'f', 6, 64), strconv.FormatFloat(rec.GS, 'f', 6, 64),
			strconv.FormatFloat(rec.HighP99, 'f', 6, 64), strconv.FormatFloat(rec.HighP95, 'f', 6, 64),
			strconv.FormatFloat(rec.HighMp99, 'f', 6, 64), strconv.FormatFloat(rec.HighMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.LowP99, 'f', 6, 64), strconv.FormatFloat(rec.LowP95, 'f', 6, 64),
			strconv.FormatFloat(rec.LowMp99, 'f', 6, 64), strconv.FormatFloat(rec.LowMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.DsiP99, 'f', 6, 64), strconv.FormatFloat(rec.DsiP95, 'f', 6, 64),
			strconv.FormatFloat(rec.DsiMp99, 'f', 6, 64), strconv.FormatFloat(rec.DsiMp95, 'f', 6, 64),
			strconv.FormatFloat(rec.GsP99, 'f', 6, 64), strconv.FormatFloat(rec.GsP95, 'f', 6, 64),
		}
		if wErr := writer.Write(record); wErr != nil {
			log.Fatalf("Error: %s\n", wErr)
		}
	}

}

func calculateStatsRecord(rec TwoBulkTwoParentsRecord, rep int, resSmAF float64, susSmAF float64) TwoBulkTwoParentsRecord {
	lGT := rec.LowParGT
	ref := rec.Ref
	hBDP := rec.HighBulkDP
	lBDP := rec.LowBulkDP
	hBAD := rec.HighBulkAD
	lBAD := rec.LowBulkAD

	lowParGT := strings.Split(lGT, "/")[0]
	highBulkRefAltAD := strings.Split(hBAD, ",")
	lowBulkRefAltAD := strings.Split(lBAD, ",")

	if ref == lowParGT {
		highBulkRefAD, _ := strconv.ParseFloat(highBulkRefAltAD[0], 64)
		highBulkAltAD, _ := strconv.ParseFloat(highBulkRefAltAD[1], 64)

		lowBulkRefAD, _ := strconv.ParseFloat(lowBulkRefAltAD[0], 64)
		lowBulkAltAD, _ := strconv.ParseFloat(lowBulkRefAltAD[1], 64)

		rec.HighSI = math.Round(highBulkAltAD/float64(hBDP)*1e6) / 1e6
		rec.LowSI = math.Round(lowBulkAltAD/float64(lBDP)*1e6) / 1e6
		rec.DSI = math.Round((rec.HighSI-rec.LowSI)*1e6) / 1e6
		rec.GS = math.Round(gStatistic(int(highBulkRefAD), int(highBulkAltAD), int(lowBulkRefAD), int(lowBulkAltAD))*1e6) / 1e6
		threshMap := twoBulkThresholdsCached(float64(hBDP), float64(lBDP), resSmAF, susSmAF, rep)
		rec.HighP99 = threshMap["res99"]
		rec.HighP95 = threshMap["res95"]
		rec.HighMp99 = threshMap["resL99"]
		rec.HighMp95 = threshMap["resL95"]
		rec.LowP99 = threshMap["sus99"]
		rec.LowP95 = threshMap["sus95"]
		rec.LowMp99 = threshMap["susL99"]
		rec.LowMp95 = threshMap["susL95"]
		rec.DsiP99 = threshMap["dsi99"]
		rec.DsiP95 = threshMap["dsi95"]
		rec.DsiMp99 = threshMap["dsiL99"]
		rec.DsiMp95 = threshMap["dsiL95"]
		rec.GsP99 = threshMap["gs99"]
		rec.GsP95 = threshMap["gs95"]
	} else {
		highBulkRefAD, _ := strconv.ParseFloat(highBulkRefAltAD[1], 64)
		highBulkAltAD, _ := strconv.ParseFloat(highBulkRefAltAD[0], 64)

		lowBulkRefAD, _ := strconv.ParseFloat(lowBulkRefAltAD[1], 64)
		lowBulkAltAD, _ := strconv.ParseFloat(lowBulkRefAltAD[0], 64)

		rec.HighSI = math.Round(highBulkAltAD/float64(hBDP)*1e6) / 1e6
		rec.LowSI = math.Round(lowBulkAltAD/float64(lBDP)*1e6) / 1e6
		rec.DSI = math.Round((rec.HighSI-rec.LowSI)*1e6) / 1e6
		rec.GS = math.Round(gStatistic(int(highBulkRefAD), int(highBulkAltAD), int(lowBulkRefAD), int(lowBulkAltAD))*1e6) / 1e6

		threshMap := twoBulkThresholdsCached(float64(hBDP), float64(lBDP), resSmAF, susSmAF, rep)
		rec.HighP99 = threshMap["res99"]
		rec.HighP95 = threshMap["res95"]
		rec.HighMp99 = threshMap["resL99"]
		rec.HighMp95 = threshMap["resL95"]
		rec.LowP99 = threshMap["sus99"]
		rec.LowP95 = threshMap["sus95"]
		rec.LowMp99 = threshMap["susL99"]
		rec.LowMp95 = threshMap["susL95"]
		rec.DsiP99 = threshMap["dsi99"]
		rec.DsiP95 = threshMap["dsi95"]
		rec.DsiMp99 = threshMap["dsiL99"]
		rec.DsiMp95 = threshMap["dsiL95"]
		rec.GsP99 = threshMap["gs99"]
		rec.GsP95 = threshMap["gs95"]
	}

	return rec
}

func calculateStatsRecordBulksOnly(rec TwoBulkOnlyRecord, rep int, resSmAF float64, susSmAF float64) TwoBulkOnlyRecord {

	hBDP := rec.HighBulkDP
	lBDP := rec.LowBulkDP
	hBAD := rec.HighBulkAD
	lBAD := rec.LowBulkAD

	highBulkRefAltAD := strings.Split(hBAD, ",")
	lowBulkRefAltAD := strings.Split(lBAD, ",")

	highBulkRefAD, _ := strconv.ParseFloat(highBulkRefAltAD[0], 64)
	highBulkAltAD, _ := strconv.ParseFloat(highBulkRefAltAD[1], 64)

	lowBulkRefAD, _ := strconv.ParseFloat(lowBulkRefAltAD[0], 64)
	lowBulkAltAD, _ := strconv.ParseFloat(lowBulkRefAltAD[1], 64)

	rec.HighSI = math.Round(highBulkAltAD/float64(hBDP)*1e6) / 1e6
	rec.LowSI = math.Round(lowBulkAltAD/float64(lBDP)*1e6) / 1e6
	rec.DSI = math.Round((rec.HighSI-rec.LowSI)*1e6) / 1e6
	rec.GS = math.Round(gStatistic(int(highBulkRefAD), int(highBulkAltAD), int(lowBulkRefAD), int(lowBulkAltAD))*1e6) / 1e6
	threshMap := twoBulkThresholdsCached(float64(hBDP), float64(lBDP), resSmAF, susSmAF, rep)
	rec.HighP99 = threshMap["res99"]
	rec.HighP95 = threshMap["res95"]
	rec.HighMp99 = threshMap["resL99"]
	rec.HighMp95 = threshMap["resL95"]
	rec.LowP99 = threshMap["sus99"]
	rec.LowP95 = threshMap["sus95"]
	rec.LowMp99 = threshMap["susL99"]
	rec.LowMp95 = threshMap["susL95"]
	rec.DsiP99 = threshMap["dsi99"]
	rec.DsiP95 = threshMap["dsi95"]
	rec.DsiMp99 = threshMap["dsiL99"]
	rec.DsiMp95 = threshMap["dsiL95"]
	rec.GsP99 = threshMap["gs99"]
	rec.GsP95 = threshMap["gs95"]
	return rec
}
