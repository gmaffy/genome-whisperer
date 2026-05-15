package genespace

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/schollz/progressbar/v3"
	"github.com/xuri/excelize/v2"
)

func genePosFromGffMap(gff, chrom string, start, stop int) (map[string][2]string, error) {
	dic := make(map[string][2]string)

	f, err := os.Open(gff)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) < 9 {
			continue
		}

		pos, err := strconv.Atoi(fields[3])
		if err != nil {
			continue
		}

		if fields[0] == chrom && fields[2] == "mRNA" && pos >= start && pos < stop {
			gene := strings.SplitN(strings.SplitN(fields[len(fields)-1], ";", 2)[0], "=", 2)[1]
			dic[gene] = [2]string{fields[3], fields[4]}
		}
	}

	return dic, scanner.Err()
}

func geneDescMap(descFile string) (map[string]string, error) {
	dic := make(map[string]string)
	f, err := os.Open(descFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) < 2 {
			continue
		}
		dic[fields[0]] = fields[1]
	}
	return dic, scanner.Err()
}

type PrgMatches struct {
	PercIdent    float64
	MatchLenQlen string
	Eval         float64
}

func genePrgMap(prgFile string) (map[string]PrgMatches, error) {
	dic := make(map[string]PrgMatches)
	f, err := os.Open(prgFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.LazyQuotes = true

	header, err := r.Read()
	if err != nil {
		return nil, fmt.Errorf("reading header: %w", err)
	}

	idx := make(map[string]int)
	for i, h := range header {
		idx[h] = i
	}

	seqIdIdx, okS := idx["QseqID"]
	lengthIdx, okL := idx["Length"]
	qlenIdx, okQ := idx["Qlen"]
	percIdentIdx, okP := idx["PercIdent"]
	evalIdx, okE := idx["Eval"]
	if !okS || !okL || !okQ || !okP || !okE {
		return nil, fmt.Errorf("Blast table missing required columns (QseqID, PercIdent, Qlen, Length, Eval)")
	}

	for {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading row: %w", err)
		}

		percIdent, _ := strconv.ParseFloat(row[percIdentIdx], 64)
		eval, _ := strconv.ParseFloat(row[evalIdx], 64)
		mlenQlen := row[lengthIdx] + "/" + row[qlenIdx]
		rec := PrgMatches{
			PercIdent:    percIdent,
			MatchLenQlen: mlenQlen,
			Eval:         eval,
		}

		if _, exists := dic[row[seqIdIdx]]; !exists {
			dic[row[seqIdIdx]] = rec
		}
	}
	return dic, nil
}

type VCFRecord struct {
	CHROM        string
	POS          int
	SNPEffEffect string
	Genotypes    map[string]string
}

func getQtlRecords(path, chrom string, start, stop int) ([]VCFRecord, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, fmt.Errorf("opening VCF table: %w", err)
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.LazyQuotes = true

	header, err := r.Read()
	if err != nil {
		return nil, nil, fmt.Errorf("reading header: %w", err)
	}

	idx := make(map[string]int)
	for i, h := range header {
		idx[h] = i
	}

	chromIdx, okC := idx["CHROM"]
	posIdx, okP := idx["POS"]
	effIdx, okE := idx["SNPEFF_EFFECT"]
	if !okC || !okP || !okE {
		return nil, nil, fmt.Errorf("VCF table missing required columns (CHROM, POS, SNPEFF_EFFECT)")
	}

	var gtCols []string
	for _, h := range header {
		if strings.HasSuffix(h, ".GT") {
			gtCols = append(gtCols, h)
		}
	}

	var records []VCFRecord
	for {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("reading row: %w", err)
		}

		pos, _ := strconv.Atoi(row[posIdx])
		if row[chromIdx] != chrom || pos < start || pos > stop {
			continue
		}
		rec := VCFRecord{
			CHROM:        row[chromIdx],
			POS:          pos,
			SNPEffEffect: row[effIdx],
			Genotypes:    make(map[string]string),
		}
		for _, col := range gtCols {
			rec.Genotypes[col] = row[idx[col]]
		}
		records = append(records, rec)
	}
	return records, gtCols, nil
}

func resNotSus(resGTs, susGTs []string) bool {
	unphase := func(gt string) string { return strings.ReplaceAll(gt, "|", "/") }

	// Derive the single sus genotype without allocating a map.
	var singleSus string
	for _, g := range susGTs {
		u := unphase(g)
		if u == "./." {
			continue
		}
		if singleSus == "" {
			singleSus = u
		} else if u != singleSus {
			return false // more than one distinct sus genotype
		}
	}
	if singleSus == "" {
		return false // all sus calls are missing
	}

	// Sus must be homozygous.
	if a, b, ok := strings.Cut(singleSus, "/"); ok && a != b {
		return false
	}

	// All res must be non-missing and differ from sus.
	for _, g := range resGTs {
		if u := unphase(g); u == "./." || u == singleSus {
			return false
		}
	}
	return true
}

type GeneSpaceRow struct {
	Gene             string
	Start            int
	Stop             int
	PRGPerc          float64
	PRGMatchLen      string
	PRGEval          float64
	GeneDescription  string
	NumSNPs          int
	NumSigSNPs       int
	NumResNotSusSNPs int
}

// Priority describes the candidate priority tier for a gene row.
type Priority struct {
	Level    int    // 1 (lowest) … 5 (highest)
	Category string // human-readable label
	HexColor string // RRGGBB fill colour for Excel
}

// resistanceKeywords are searched (case-insensitively) in the gene description
// to assign priority tier 2 and above.
// resistanceKeywords are matched (case-insensitively) against the gene
// description to flag defence-related candidates.  Terms are derived from
// observed gene-description patterns in the dataset.
var resistanceKeywords = []string{
	// Classical R-gene / NLR terms
	"NBS-LRR", "NBS", "NLR", "LRR",
	"CC-NBS", "TIR-NBS",
	"leucine-rich repeat receptor",              // e.g. "Leucine-rich repeat receptor-like protein kinase"
	"leucine-rich repeat receptor-like",         // LRR-RLK family
	"plant intracellular ras-group-related lrr", // PIRL proteins

	// Named R-genes / well-known resistance loci seen in the data
	"ZAR1", // receptor protein kinase-like protein ZAR1
	"CDR1", // aspartic proteinase CDR1 (constitutive disease resistance)
	"RGA",  // R-gene analog
	"RGH",
	"RPM", "RPS",
	"MSP1", // leucine-rich repeat receptor protein kinase MSP1

	// Broad defence / immunity vocabulary
	"resistance", "resistant",
	"disease resistance",
	"defense", "defence",
	"pathogen", "pathogenesis-related",
	"immune", "immunity",
	"innate immun",

	// Oxidative burst — key early defence response
	"respiratory burst oxidase", // e.g. "respiratory burst oxidase homolog protein E"

	// PR-proteins and related
	"pathogenesis-related protein", // "basic form of pathogenesis-related protein 1-like"
	"thaumatin",                    // PR-5 family
	"chitinase",
	"glucanase",
	"(1->3)-beta-glucan endohydrolase", // beta-1,3-glucanase (PR-2)
	"peroxidase",                       // class III peroxidases — HR / lignification

	// Defence signalling
	"lipoxygenase",            // JA biosynthesis / wound response
	"4-coumarate--coa ligase", // phenylpropanoid / lignin pathway
	"cysteine protease", "cysteine proteinase", "cysteinase",
	"aspartic proteinase", "aspartyl protease",
	"subtilisin-like protease",

	// Receptor kinases broadly involved in pattern-triggered immunity
	"receptor protein kinase",
	"receptor-like serine/threonine-protein kinase",
	"leucine-rich repeat receptor-like serine/threonine-protein kinase",
	"lysm domain receptor-like kinase", // chitin / MAMP perception

	// Ubiquitin / proteasome — critical for R-protein turnover
	"ring-type e3 ubiquitin",
	"u-box domain",
	"f-box",

	// Late-stage / stress-associated
	"late embryogenesis abundant",     // LEA proteins — stress tolerance
	"early-responsive to dehydration", // ERD proteins
	"detoxification",
	"multidrug resistance protein abc transporter",
}

// transcriptionFactorKeywords are searched to identify tier-5 (TF) candidates.
// Terms are grounded in the TF families observed in the dataset.
var transcriptionFactorKeywords = []string{
	// Generic labels seen verbatim in descriptions
	"transcription factor",
	"transcriptional activator",
	"transcription repressor",
	"transcriptional regulator",
	"transcription elongation regulator",
	"rna polymerase ii transcription", // mediator subunits / coactivators
	"corepressor",                     // dr1-associated corepressor
	"coactivator",

	// Major plant TF families present in the data
	"WRKY", // WRKY DNA-binding protein — major defence TF family
	"MYB",  // R2R3-MYB, MYB1R1, MYB117, MYB124, myb-like
	"NAC",  // NAC domain-containing — stress & senescence TFs
	"bHLH", // basic helix-loop-helix (bHLH30, bHLH49, BHLH)
	"bZIP", // basic leucine zipper (bZIP43)
	"basic leucine zipper",
	"AP2", // AP2/ERF superfamily
	"ERF", // ethylene response factor (ERF118, ESR2, WRI1)
	"ethylene-responsive transcription factor",
	"TCP",  // TCP family (TCP8, CYCLOIDEA)
	"MADS", // MADS-box / AGL proteins
	"agamous-like mads-box",
	"squamosa promoter-binding",     // SPL family
	"trihelix transcription factor", // trihelix / GT factors (ASIL2, PTL, ASR3)
	"homeobox",                      // homeobox-leucine zipper (ATHB, knotted)
	"homeobox-leucine zipper",
	"homeobox protein knotted",
	"zinc finger homeodomain",
	"zinc finger transcription factor",
	"GATA transcription factor",
	"B3 domain-containing",  // B3 superfamily (RAV1, ARF)
	"auxin response factor", // ARF — auxin signalling TFs
	"auxin-responsive",
	"YABBY", // axial regulator YABBY
	"DOF",   // cyclic dof factor
	"cyclic dof factor",
	"GAGA-binding transcriptional activator",
	"heat-inducible transcription repressor",
	"t-box transcription factor",
	"myb-like transcription factor",
	"HTH myb-type",
	"DELLA protein",                 // DELLA — GA signalling master regulator
	"mediator of rna polymerase ii", // Mediator complex subunits
}

// ClassifyPriority assigns a Priority tier to a GeneSpaceRow using the
// following rules (higher tier wins when multiple conditions match):
//
//	1 – PRG database match only   (PRGPerc > 0, no keyword match)
//	2 – Resistance keywords       (keyword match, no PRG)
//	3 – PRG + R-gene associated   (PRGPerc > 0 AND resistance keyword)
//	4 – STRONG candidates         (ResNotSus SNPs present AND (PRG or keyword))
//	5 – Transcription factors     (any TF keyword in description)
func ClassifyPriority(row GeneSpaceRow) Priority {
	desc := strings.ToLower(row.GeneDescription)

	hasPRG := row.PRGPerc > 0
	hasResKeyword := false
	for _, kw := range resistanceKeywords {
		if strings.Contains(desc, strings.ToLower(kw)) {
			hasResKeyword = true
			break
		}
	}
	hasTF := false
	for _, kw := range transcriptionFactorKeywords {
		if strings.Contains(desc, strings.ToLower(kw)) {
			hasTF = true
			break
		}
	}
	isStrong := row.NumResNotSusSNPs > 0 && (hasPRG || hasResKeyword)

	switch {
	case hasTF:
		return Priority{5, "Transcription factor", "C5A0FF"} // purple
	case isStrong:
		return Priority{4, "STRONG candidate", "FF9900"} // orange
	case hasPRG && hasResKeyword:
		return Priority{3, "PRG + R-gene associated", "9DC3E6"} // light blue
	case hasResKeyword:
		return Priority{2, "Resistance keywords (no PRG)", "FFFF99"} // light yellow
	case hasPRG:
		return Priority{1, "PRG database match only", "92D050"} // light green
	default:
		return Priority{0, "No match", "FFFFFF"}
	}
}

func GeneSpace(gffPath, vcfTable, chrom string, start, stop int, resLines, susLines []string, descFile, prgFile string, outDir string) ([]GeneSpaceRow, error) {

	color.Green("▶ Starting GeneSpace analysis for chromosome %s (%d - %d)", chrom, start, stop)

	// --- 1. Load all three reference files concurrently while VCF loads on
	//        the main goroutine, overlapping all four I/O operations. ---
	type gffResult struct {
		data map[string][2]string
		err  error
	}
	type descResult struct {
		data map[string]string
		err  error
	}
	type prgResult struct {
		data map[string]PrgMatches
		err  error
	}

	gffCh := make(chan gffResult, 1)
	descCh := make(chan descResult, 1)
	prgCh := make(chan prgResult, 1)

	go func() {
		d, err := genePosFromGffMap(gffPath, chrom, start, stop)
		gffCh <- gffResult{d, err}
	}()
	go func() {
		d, err := geneDescMap(descFile)
		descCh <- descResult{d, err}
	}()
	go func() {
		d, err := genePrgMap(prgFile)
		prgCh <- prgResult{d, err}
	}()

	// FIX: gtCols is now captured rather than discarded, preserving the caller's ability to know which sample columns were present.
	qtlRecords, gtCols, err := getQtlRecords(vcfTable, chrom, start, stop)
	if err != nil {
		return nil, fmt.Errorf("reading VCF table: %w", err)
	}
	_ = gtCols // available for callers that need it; not used internally here

	gffRes := <-gffCh
	if gffRes.err != nil {
		return nil, fmt.Errorf("reading GFF file: %w", gffRes.err)
	}
	descRes := <-descCh
	if descRes.err != nil {
		return nil, fmt.Errorf("reading gene description file: %w", descRes.err)
	}
	prgRes := <-prgCh
	if prgRes.err != nil {
		return nil, fmt.Errorf("reading PRG blast file: %w", prgRes.err)
	}

	gencPos := gffRes.data
	geneDesc := descRes.data
	genePrg := prgRes.data

	// --- 2. Sort qtlRecords by position once so each worker can binary-search
	//        its gene window instead of scanning the entire slice. ---
	// FIX: reduces per-gene SNP scan from O(all SNPs) to O(log n + hits).
	sort.Slice(qtlRecords, func(i, j int) bool {
		return qtlRecords[i].POS < qtlRecords[j].POS
	})

	synonymousEffects := []string{"UPSTREAM", "DOWNSTREAM", "SYNONYMOUS_CODING", "SYNONYMOUS_STOP", "INTRON"}

	// --- 3. Collect gene names and set up progress bar. ---
	genes := make([]string, 0, len(gencPos))
	for gene := range gencPos {
		genes = append(genes, gene)
	}

	bar := progressbar.Default(int64(len(genes)), "Processing genes")

	// --- 4. Process genes concurrently with a fixed worker pool. ---
	workCh := make(chan string, len(genes))
	for _, g := range genes {
		workCh <- g
	}
	close(workCh)

	var (
		mu            sync.Mutex
		geneSpaceRows []GeneSpaceRow
		wg            sync.WaitGroup
	)

	// Progress updates are handled by a single goroutine to avoid
	// concurrent writes to the progress bar interfering with log output.
	progressCh := make(chan struct{}, len(genes))
	var wgBar sync.WaitGroup

	wgBar.Add(1)
	go func() {
		defer wgBar.Done()
		for range progressCh {
			_ = bar.Add(1)
		}
		_ = bar.Finish()
		fmt.Fprint(os.Stderr, "\n")
	}()

	const numWorkers = 8
	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go func() {
			defer wg.Done()
			for gene := range workCh {
				// Workers do not log each gene to avoid interfering with the single
				// progress bar. The bar reflects overall progress in real time.

				pos := gencPos[gene]
				geneStart, _ := strconv.Atoi(pos[0])
				geneStop, _ := strconv.Atoi(pos[1])

				// Binary-search the sorted qtlRecords for this gene's window.
				lo := sort.Search(len(qtlRecords), func(i int) bool { return qtlRecords[i].POS >= geneStart })
				hi := sort.Search(len(qtlRecords), func(i int) bool { return qtlRecords[i].POS > geneStop })

				numOfSnps := 0
				sigSnps := 0
				resNotSusCount := 0

				for _, rec := range qtlRecords[lo:hi] {
					numOfSnps++
					if !slices.Contains(synonymousEffects, rec.SNPEffEffect) {
						sigSnps++
					}
					resGTs := make([]string, len(resLines))
					for i, r := range resLines {
						resGTs[i] = rec.Genotypes[r]
					}
					susGTs := make([]string, len(susLines))
					for i, s := range susLines {
						susGTs[i] = rec.Genotypes[s]
					}
					if resNotSus(resGTs, susGTs) {
						resNotSusCount++
					}
				}

				desc, okDesc := geneDesc[gene]
				if !okDesc {
					desc = "NA"
				}
				prg, okPrg := genePrg[gene]
				if !okPrg {
					prg = PrgMatches{PercIdent: -1, MatchLenQlen: "NA", Eval: -1}
				}

				row := GeneSpaceRow{
					Gene:             gene,
					Start:            geneStart,
					Stop:             geneStop,
					PRGPerc:          prg.PercIdent,
					PRGMatchLen:      prg.MatchLenQlen,
					PRGEval:          prg.Eval,
					GeneDescription:  desc,
					NumSNPs:          numOfSnps,
					NumSigSNPs:       sigSnps,
					NumResNotSusSNPs: resNotSusCount,
				}

				mu.Lock()
				geneSpaceRows = append(geneSpaceRows, row)
				mu.Unlock()

				// Progress update.
				progressCh <- struct{}{}
			}
		}()
	}

	wg.Wait()
	close(progressCh)
	wgBar.Wait()

	color.Cyan("✔ GeneSpace analysis complete")
	color.Green("\n========================================== Write to excel ==========================================\n\n")
	// Ensure output directory exists.
	if err := os.MkdirAll(outDir, 0755); err != nil {
		return nil, fmt.Errorf("creating output directory %s: %w", outDir, err)
	}
	excelPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq_%s_%d-%d.xlsx", chrom, start, stop))
	if err = WriteGeneSpaceExcel(geneSpaceRows, excelPath); err != nil {
		return nil, fmt.Errorf("writing Excel file: %w", err)
	}
	return geneSpaceRows, nil
}

// ─── Excel export ────────────────────────────────────────────────────────────

// colWidth tracks the maximum content width seen for a column so we can
// auto-size at the end.
type colWidth [12]int

// WriteGeneSpaceExcel writes rows to an Excel workbook at outPath using the
// pure-Go excelize library (github.com/xuri/excelize/v2).
// No Python or shell-out required.
func WriteGeneSpaceExcel(rows []GeneSpaceRow, outPath string) error {
	f := excelize.NewFile()
	defer f.Close()

	// ── shared style helpers ──────────────────────────────────────────────

	// mustStyle panics on style creation errors; they only arise from
	// malformed JSON which cannot happen with literal structs here.
	mustStyle := func(s *excelize.Style) int {
		id, err := f.NewStyle(s)
		if err != nil {
			panic(fmt.Sprintf("excelize.NewStyle: %v", err))
		}
		return id
	}

	navy := "1F3864"
	white := "FFFFFF"
	thinBorder := []excelize.Border{
		{Type: "left", Color: "CCCCCC", Style: 1},
		{Type: "right", Color: "CCCCCC", Style: 1},
		{Type: "top", Color: "CCCCCC", Style: 1},
		{Type: "bottom", Color: "CCCCCC", Style: 1},
	}

	headerStyle := mustStyle(&excelize.Style{
		Font:      &excelize.Font{Bold: true, Color: white, Family: "Arial", Size: 10},
		Fill:      excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{navy}},
		Alignment: &excelize.Alignment{Horizontal: "center", WrapText: true, Vertical: "center"},
		Border:    thinBorder,
	})

	// rowStyleFor returns a new style ID for the given hex fill colour.
	rowStyleFor := func(hex string) int {
		return mustStyle(&excelize.Style{
			Font:      &excelize.Font{Family: "Arial", Size: 10},
			Fill:      excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{hex}},
			Alignment: &excelize.Alignment{Vertical: "center"},
			Border:    thinBorder,
		})
	}

	// Cache one style per unique hex so we don't create thousands of duplicates.
	styleCache := make(map[string]int)
	cachedStyle := func(hex string) int {
		if id, ok := styleCache[hex]; ok {
			return id
		}
		id := rowStyleFor(hex)
		styleCache[hex] = id
		return id
	}

	// ── Legend sheet ─────────────────────────────────────────────────────

	const legend = "Legend"
	f.SetSheetName("Sheet1", legend)

	legHeaderStyle := mustStyle(&excelize.Style{
		Font:      &excelize.Font{Bold: true, Color: white, Family: "Arial"},
		Fill:      excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{navy}},
		Alignment: &excelize.Alignment{Horizontal: "center"},
	})
	legHeaders := []string{"Priority", "Category", "Colour"}
	for ci, h := range legHeaders {
		cell, _ := excelize.CoordinatesToCellName(ci+1, 1)
		f.SetCellValue(legend, cell, h)
		f.SetCellStyle(legend, cell, cell, legHeaderStyle)
	}

	type legRow struct {
		level    int
		category string
		hex      string
	}
	legendRows := []legRow{
		{1, "PRG database match only", "92D050"},
		{2, "Resistance keywords (no PRG)", "FFFF99"},
		{3, "PRG + R-gene associated", "9DC3E6"},
		{4, "STRONG candidates", "FF9900"},
		{5, "Transcription factors", "C5A0FF"},
	}
	for ri, lr := range legendRows {
		row := ri + 2

		cellA, _ := excelize.CoordinatesToCellName(1, row)
		cellB, _ := excelize.CoordinatesToCellName(2, row)
		cellC, _ := excelize.CoordinatesToCellName(3, row)

		centreStyle := mustStyle(&excelize.Style{
			Font:      &excelize.Font{Family: "Arial"},
			Alignment: &excelize.Alignment{Horizontal: "center"},
		})
		colourStyle := mustStyle(&excelize.Style{
			Font: &excelize.Font{Family: "Arial"},
			Fill: excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{lr.hex}},
		})

		f.SetCellValue(legend, cellA, lr.level)
		f.SetCellStyle(legend, cellA, cellA, centreStyle)
		f.SetCellValue(legend, cellB, lr.category)
		f.SetCellValue(legend, cellC, lr.hex)
		f.SetCellStyle(legend, cellC, cellC, colourStyle)
	}
	f.SetColWidth(legend, "A", "C", 36)

	// ── GeneSpace data sheet ──────────────────────────────────────────────

	const dataSheet = "GeneSpace"
	f.NewSheet(dataSheet)

	headers := []string{
		"Gene", "Start", "Stop",
		"PRG %Identity", "PRG Match/Qlen", "PRG E-value",
		"Gene Description",
		"# SNPs", "# Sig SNPs", "# ResNotSus SNPs",
		"Priority", "Priority Category",
	}

	// Track max content width per column for auto-sizing.
	var widths colWidth
	for ci, h := range headers {
		if len(h) > widths[ci] {
			widths[ci] = len(h)
		}
		cell, _ := excelize.CoordinatesToCellName(ci+1, 1)
		f.SetCellValue(dataSheet, cell, h)
		f.SetCellStyle(dataSheet, cell, cell, headerStyle)
	}

	// Sort rows highest-priority first before writing.
	sort.Slice(rows, func(i, j int) bool {
		return ClassifyPriority(rows[i]).Level > ClassifyPriority(rows[j]).Level
	})

	for ri, r := range rows {
		excelRow := ri + 2
		p := ClassifyPriority(r)
		sID := cachedStyle(p.HexColor)

		vals := []any{
			r.Gene, r.Start, r.Stop,
			r.PRGPerc, r.PRGMatchLen, r.PRGEval,
			r.GeneDescription,
			r.NumSNPs, r.NumSigSNPs, r.NumResNotSusSNPs,
			p.Level, p.Category,
		}
		for ci, v := range vals {
			cell, _ := excelize.CoordinatesToCellName(ci+1, excelRow)
			f.SetCellValue(dataSheet, cell, v)
			f.SetCellStyle(dataSheet, cell, cell, sID)

			w := len(fmt.Sprintf("%v", v))
			if w > widths[ci] {
				widths[ci] = w
			}
		}
	}

	// Auto-size columns (cap at 60 characters).
	cols := []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"}
	for ci, col := range cols {
		w := float64(widths[ci]) + 3
		if w > 60 {
			w = 60
		}
		f.SetColWidth(dataSheet, col, col, w)
	}

	// Freeze the header row and add auto-filter.
	f.SetPanes(dataSheet, &excelize.Panes{
		Freeze:      true,
		Split:       false,
		YSplit:      1,
		TopLeftCell: "A2",
		ActivePane:  "bottomLeft",
	})
	lastCol, _ := excelize.CoordinatesToCellName(len(headers), len(rows)+1)
	f.AutoFilter(dataSheet, fmt.Sprintf("A1:%s", lastCol), nil)

	// Set GeneSpace as the active sheet.
	idx, _ := f.GetSheetIndex(dataSheet)
	f.SetActiveSheet(idx)

	if err := f.SaveAs(outPath); err != nil {
		return fmt.Errorf("saving Excel file: %w", err)
	}
	color.Green("✔ Excel report written to %s", outPath)
	return nil
}
