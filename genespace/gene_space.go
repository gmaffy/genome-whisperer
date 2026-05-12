// ========================================================================== //
//                   Author: Godwin Mafireyi                                  //
//                   Email:  mafireyi@gmail.com                               //
//                   Tel:    +27788673138                                     //
//                  GeneSpa                                                   //
// ========================================================================== //

package genespace

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/fatih/color"
	"github.com/schollz/progressbar/v3"
	"github.com/xuri/excelize/v2"
)

// GeneRange holds the start and stop positions for a gene.
type GeneRange struct {
	Start string
	Stop  string
}

// GeneSpaceRow represents one row in the output Excel sheet.
type GeneSpaceRow struct {
	Gene             string
	Start            string
	Stop             string
	PRGPerc          string
	PRGMatchLen      string
	PRGEval          string
	GeneDescription  string
	NumSNPs          int
	NumSigSNPs       int
	NumResUniqueSNPs int
	NumResSigSNPs    int
}

// -------------------------------------------------------------------------- //
// geneStartStop parses a GFF file and returns genes within [start, stop) on
// the given chromosome.  Only 'mRNA' features are considered.
// -------------------------------------------------------------------------- //
func geneStartStop(gffPath, chrom string, start, stop int) (map[string]GeneRange, error) {
	result := make(map[string]GeneRange)

	f, err := os.Open(gffPath)
	if err != nil {
		return nil, fmt.Errorf("opening GFF: %w", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 9 {
			continue
		}
		if fields[0] != chrom || fields[2] != "mRNA" {
			continue
		}
		pos, err := strconv.Atoi(fields[3])
		if err != nil {
			continue
		}
		if pos < start || pos >= stop {
			continue
		}
		// Extract gene name from attributes column (ID=<name>;...)
		attrs := strings.Split(fields[8], ";")
		if len(attrs) == 0 {
			continue
		}
		parts := strings.SplitN(attrs[0], "=", 2)
		if len(parts) != 2 {
			continue
		}
		gene := parts[1]
		result[gene] = GeneRange{Start: fields[3], Stop: fields[4]}
	}
	return result, scanner.Err()
}

// -------------------------------------------------------------------------- //
// resIsUnique mirrors the Python res_is_unique logic.
// resGTs and susGTs are slices of genotype strings (e.g. "0/1", "1|1", "./.").
// -------------------------------------------------------------------------- //
func resIsUnique(resGTs, susGTs []string) bool {
	// Normalise phased to unphased
	unphase := func(gt string) string { return strings.ReplaceAll(gt, "|", "/") }

	unphSus := make([]string, len(susGTs))
	for i, g := range susGTs {
		unphSus[i] = unphase(g)
	}
	unphRes := make([]string, len(resGTs))
	for i, g := range resGTs {
		unphRes[i] = unphase(g)
	}

	// Non-missing susceptibles
	nonMissingSus := []string{}
	for _, s := range unphSus {
		if s != "./." {
			nonMissingSus = append(nonMissingSus, s)
		}
	}
	if len(nonMissingSus) == 0 {
		return false
	}

	// All susceptibles must share the same genotype
	susSet := make(map[string]struct{})
	for _, s := range nonMissingSus {
		susSet[s] = struct{}{}
	}
	if len(susSet) != 1 {
		return false
	}

	// That single genotype must be homozygous
	singleSus := nonMissingSus[0]
	alleles := strings.Split(singleSus, "/")
	if len(alleles) == 2 && alleles[0] != alleles[1] {
		return false // heterozygous
	}

	// res unique means res is not sus.
	// i.e sus is homozygous and different from res.
	// We already verified sus is homozygous.
	// Now check that ALL resistant lines are different from sus.
	// And resistant lines must not be missing (as per previous requirement,
	// though "res is not sus" might imply missing is okay as long as it's not sus,
	// but usually we want valid genotypes for comparisons).
	// Let's stick to: if any res is missing, we might not want to call it unique,
	// or if any res matches sus, it's definitely not unique.

	for _, r := range unphRes {
		if r == "./." {
			return false
		}
		if r == singleSus {
			return false
		}
	}

	return true
}

// -------------------------------------------------------------------------- //
// VCFRecord is a minimal, flat representation of one row in the VCF table.
// Column names that end with ".GT" are stored in the Genotypes map.
// -------------------------------------------------------------------------- //
type VCFRecord struct {
	CHROM        string
	POS          int
	SNPEffEffect string
	Genotypes    map[string]string // key: "<sample>.GT"
}

// -------------------------------------------------------------------------- //
// readVCFTable reads a tab-separated VCF table into a slice of VCFRecords.
// The file must have a header row; columns "CHROM", "POS", and
// "SNPEFF_EFFECT" are mandatory; any column ending in ".GT" is a genotype.
// -------------------------------------------------------------------------- //
func readVCFTable(path string) ([]VCFRecord, []string, error) {
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

	// Index key columns
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

	// Collect all GT sample column names
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

// -------------------------------------------------------------------------- //
// resistantLinesUnique filters records to those where the resistant lines
// carry a genotype that is uniquely different from all susceptibles.
// -------------------------------------------------------------------------- //
func resistantLinesUnique(records []VCFRecord, resLines []string, susLines []string) []VCFRecord {
	var out []VCFRecord
	for _, rec := range records {
		resGTs := make([]string, len(resLines))
		for i, r := range resLines {
			resGTs[i] = rec.Genotypes[r]
		}
		susGTs := make([]string, len(susLines))
		for i, s := range susLines {
			susGTs[i] = rec.Genotypes[s]
		}
		if resIsUnique(resGTs, susGTs) {
			out = append(out, rec)
		}
	}
	return out
}

// -------------------------------------------------------------------------- //
// Stubs for external lookups (gene_descriptions / add_prg).
// Replace with real implementations or database lookups as needed.
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// nonSigEffects mirrors the Python filter (rows whose effect IS in this list
// are counted as significant — note the original code uses isin for sig but
// ~isin for unique_sig, so both sets are captured below).
// -------------------------------------------------------------------------- //
var nonCodingEffects = map[string]bool{
	"UPSTREAM":          true,
	"DOWNSTREAM":        true,
	"SYNONYMOUS_CODING": true,
	"SYNONYMOUS_STOP":   true,
	"INTRON":            true,
}

// -------------------------------------------------------------------------- //
// geneSpace is the main analysis function, mirroring gene_space() in Python.
// -------------------------------------------------------------------------- //
func GeneSpace(gffPath, vcfTable, chrom string, start, stop int, resLines, susLines []string, descFile, prgFile string) ([]GeneSpaceRow, error) {
	color.Green("▶ Starting GeneSpace analysis for chromosome %s (%d - %d)", chrom, start, stop)

	fmt.Print("  Reading VCF table... ")
	records, allGTCols, err := readVCFTable(vcfTable)
	if err != nil {
		fmt.Println("FAILED")
		return nil, err
	}
	fmt.Printf("DONE (%d records, %d columns)\n", len(records), len(allGTCols))

	fmt.Print("  Filtering QTL region... ")
	var qtlRecords []VCFRecord
	for _, rec := range records {
		if rec.CHROM == chrom && rec.POS >= start && rec.POS < stop {
			qtlRecords = append(qtlRecords, rec)
		}
	}
	fmt.Printf("DONE (%d records)\n", len(qtlRecords))

	// OPTIMIZATION: Index records by position for faster lookup within gene ranges
	posToRecords := make(map[int][]VCFRecord)
	for _, rec := range qtlRecords {
		posToRecords[rec.POS] = append(posToRecords[rec.POS], rec)
	}
	// Sorted unique positions
	var sortedPos []int
	for p := range posToRecords {
		sortedPos = append(sortedPos, p)
	}
	sort.Ints(sortedPos)

	fmt.Print("  Identifying resistant-unique variants... ")
	uniqueRecords := resistantLinesUnique(qtlRecords, resLines, susLines)
	fmt.Printf("DONE (%d records)\n", len(uniqueRecords))

	posToUniqueRecords := make(map[int][]VCFRecord)
	for _, rec := range uniqueRecords {
		posToUniqueRecords[rec.POS] = append(posToUniqueRecords[rec.POS], rec)
	}
	var sortedUniquePos []int
	for p := range posToUniqueRecords {
		sortedUniquePos = append(sortedUniquePos, p)
	}
	sort.Ints(sortedUniquePos)

	fmt.Print("  Parsing GFF for genes... ")
	geneDic, err := geneStartStop(gffPath, chrom, start, stop)
	if err != nil {
		fmt.Println("FAILED")
		return nil, err
	}
	fmt.Printf("DONE (%d genes)\n", len(geneDic))

	// Load external data
	descMap := make(map[string]string)
	if descFile != "" {
		fmt.Printf("  Loading gene descriptions from %s... ", descFile)
		df, err := os.Open(descFile)
		if err == nil {
			scanner := bufio.NewScanner(df)
			for scanner.Scan() {
				flds := strings.Split(scanner.Text(), "\t")
				if len(flds) >= 2 {
					descMap[flds[0]] = flds[1]
				}
			}
			df.Close()
			fmt.Printf("DONE (%d loaded)\n", len(descMap))
		} else {
			fmt.Println("FAILED")
			color.Red("    Warning: could not open description file: %v", err)
		}
	}

	type prgEntry struct {
		percIdent    string
		qlenMatchLen string
		eval         string
	}
	prgMap := make(map[string]prgEntry)
	if prgFile != "" {
		fmt.Printf("  Loading PRG stats from %s... ", prgFile)
		pf, err := os.Open(prgFile)
		if err == nil {
			scanner := bufio.NewScanner(pf)
			if scanner.Scan() {
				header := strings.Split(scanner.Text(), "\t")
				idx := make(map[string]int)
				for i, h := range header {
					idx[h] = i
				}
				qseqIdx, ok1 := idx["QseqID"]
				percIdx, ok2 := idx["PercIdent"]
				lenIdx, ok3 := idx["Length"]
				qlenIdx, ok4 := idx["Qlen"]
				evalIdx, ok5 := idx["Eval"]
				if ok1 && ok2 && ok3 && ok4 && ok5 {
					for scanner.Scan() {
						flds := strings.Split(scanner.Text(), "\t")
						if len(flds) > max(qseqIdx, percIdx, lenIdx, qlenIdx, evalIdx) {
							// Use first occurrence for each gene
							geneID := flds[qseqIdx]
							if _, exists := prgMap[geneID]; !exists {
								prgMap[geneID] = prgEntry{
									percIdent:    flds[percIdx],
									qlenMatchLen: flds[lenIdx] + "/" + flds[qlenIdx],
									eval:         flds[evalIdx],
								}
							}
						}
					}
				}
			}
			pf.Close()
			fmt.Printf("DONE (%d loaded)\n", len(prgMap))
		} else {
			fmt.Println("FAILED")
			color.Red("    Warning: could not open PRG file: %v", err)
		}
	}

	var rows []GeneSpaceRow

	// Process genes
	fmt.Println("  Processing genes...")
	bar := progressbar.NewOptions(len(geneDic),
		progressbar.OptionEnableColorCodes(true),
		progressbar.OptionShowCount(),
		progressbar.OptionSetWidth(15),
		progressbar.OptionSetDescription("[cyan]Analyzing genes...[reset]"),
		progressbar.OptionSetTheme(progressbar.Theme{
			Saucer:        "[green]■[reset]",
			SaucerHead:    "[green]■[reset]",
			SaucerPadding: " ",
			BarStart:      "|",
			BarEnd:        "|",
		}))

	for gene, gr := range geneDic {
		bar.Describe(fmt.Sprintf("[cyan]Analyzing gene: [yellow]%s[reset]", gene))
		gStart, _ := strconv.Atoi(gr.Start)
		gStop, _ := strconv.Atoi(gr.Stop)

		// SNPs within gene
		var geneRecords []VCFRecord
		idx := sort.SearchInts(sortedPos, gStart)
		for i := idx; i < len(sortedPos); i++ {
			p := sortedPos[i]
			if p > gStop {
				break
			}
			geneRecords = append(geneRecords, posToRecords[p]...)
		}

		// Significant SNPs (non-coding effects)
		numSig := 0
		for _, rec := range geneRecords {
			if nonCodingEffects[rec.SNPEffEffect] {
				numSig++
			}
		}

		// Unique SNPs within gene
		var uniqueGene []VCFRecord
		idxU := sort.SearchInts(sortedUniquePos, gStart)
		for i := idxU; i < len(sortedUniquePos); i++ {
			p := sortedUniquePos[i]
			if p > gStop {
				break
			}
			uniqueGene = append(uniqueGene, posToUniqueRecords[p]...)
		}
		// Unique significant SNPs (coding / non-synonymous)
		resSig := 0
		for _, rec := range uniqueGene {
			if !nonCodingEffects[rec.SNPEffEffect] {
				resSig++
			}
		}

		desc := "N/A"
		if d, ok := descMap[gene]; ok {
			desc = d
		}
		stats := [3]string{"0.0", "NA", "NA"}
		if p, ok := prgMap[gene]; ok {
			stats = [3]string{p.percIdent, p.qlenMatchLen, p.eval}
		}

		rows = append(rows, GeneSpaceRow{
			Gene:             gene,
			Start:            gr.Start,
			Stop:             gr.Stop,
			PRGPerc:          stats[0],
			PRGMatchLen:      stats[1],
			PRGEval:          stats[2],
			GeneDescription:  desc,
			NumSNPs:          len(geneRecords),
			NumSigSNPs:       numSig,
			NumResUniqueSNPs: len(uniqueGene),
			NumResSigSNPs:    resSig,
		})
		bar.Add(1)
	}
	bar.Finish()
	fmt.Println()

	// Sort by Start position
	sort.Slice(rows, func(i, j int) bool {
		si, _ := strconv.Atoi(rows[i].Start)
		sj, _ := strconv.Atoi(rows[j].Start)
		return si < sj
	})

	// Write Excel output
	outFile := strings.Join([]string{chrom, strconv.Itoa(start), strconv.Itoa(stop), "_GENE_SPACE.xlsx"}, "_")
	if err := writeExcel(outFile, rows); err != nil {
		return nil, fmt.Errorf("writing Excel: %w", err)
	}
	color.Green("✔ Results written to %s", outFile)

	return rows, nil
}

// -------------------------------------------------------------------------- //
// writeExcel writes the results to an .xlsx file using excelize.
// -------------------------------------------------------------------------- //
func writeExcel(path string, rows []GeneSpaceRow) error {
	f := excelize.NewFile()
	sheet := "Sheet1"

	headers := []string{
		"GENE", "START", "STOP", "PRG_PERC", "PRG_MATCH/LEN", "PRG_EVAL",
		"GENE_DESCRIPTION", "#SNPs", "#SIG SNPs", "#RES_UNIQUE SNPs", "#RES_UNIQUE SIG SNPs",
	}
	for col, h := range headers {
		cell, _ := excelize.CoordinatesToCellName(col+1, 1)
		f.SetCellValue(sheet, cell, h)
	}

	for rowIdx, row := range rows {
		vals := []interface{}{
			row.Gene, row.Start, row.Stop,
			row.PRGPerc, row.PRGMatchLen, row.PRGEval,
			row.GeneDescription,
			row.NumSNPs, row.NumSigSNPs, row.NumResUniqueSNPs, row.NumResSigSNPs,
		}
		for col, v := range vals {
			cell, _ := excelize.CoordinatesToCellName(col+1, rowIdx+2)
			f.SetCellValue(sheet, cell, v)
		}
	}

	return f.SaveAs(path)
}

// -------------------------------------------------------------------------- //
// main
// -------------------------------------------------------------------------- //
//func main() {
//	if len(os.Args) < 8 {
//		fmt.Println("USAGE: gene_space <gff> <super_vcf> <chrom> <start> <stop> <res-lines> <sus-lines>")
//		os.Exit(1)
//	}
//
//	gff := os.Args[1]
//	superVCF := os.Args[2]
//	chrom := os.Args[3]
//	startStr := os.Args[4]
//	stopStr := os.Args[5]
//	resLinesStr := os.Args[6]
//	susLinesStr := os.Args[7]
//
//	start, err := strconv.Atoi(startStr)
//	if err != nil {
//		fmt.Fprintf(os.Stderr, "Invalid start: %v\n", err)
//		os.Exit(1)
//	}
//	stop, err := strconv.Atoi(stopStr)
//	if err != nil {
//		fmt.Fprintf(os.Stderr, "Invalid stop: %v\n", err)
//		os.Exit(1)
//	}
//
//	resLines := strings.Split(resLinesStr, ",")
//	susLines := strings.Split(susLinesStr, ",")
//
//	if _, err := GeneSpace(gff, superVCF, chrom, start, stop, resLines, susLines, "", ""); err != nil {
//		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
//		os.Exit(1)
//	}
//}
