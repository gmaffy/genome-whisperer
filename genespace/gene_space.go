// ========================================================================== //
//                   Author: Godwin Mafireyi                                  //
//                   Email:  mafireyi@gmail.com                               //
//                   Tel:    +27788673138                                     //
//                  GeneSpa                                      //
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

	// Resistant lines must not be missing
	for _, r := range unphRes {
		if r == "./." {
			return false
		}
	}

	// Resistant lines must not carry the susceptible genotype
	for _, r := range unphRes {
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
func resistantLinesUnique(records []VCFRecord, allGTCols []string, resLines []string) []VCFRecord {
	resSet := make(map[string]struct{})
	for _, r := range resLines {
		resSet[r] = struct{}{}
	}

	var susLines []string
	for _, col := range allGTCols {
		if _, ok := resSet[col]; !ok {
			susLines = append(susLines, col)
		}
	}

	fmt.Println("Resistant lines:", resLines)
	fmt.Println("Susceptible lines:", susLines)

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
func geneDesc(_ string, gene string) string {
	// TODO: implement gene_descriptions_old.gene_desc equivalent
	return "N/A"
}

func prgStats(_ string, gene string) [3]string {
	// TODO: implement add_prg_old.prot_prg_dic equivalent
	return [3]string{"0.0", "NA", "NA"}
}

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
func geneSpace(gffPath, vcfTable, chrom string, start, stop int, csResLines, species string) ([]GeneSpaceRow, error) {
	fmt.Println("Reading VCF table ...")
	records, allGTCols, err := readVCFTable(vcfTable)
	if err != nil {
		return nil, err
	}

	fmt.Println("Filtering VCF table ...")
	var qtlRecords []VCFRecord
	for _, rec := range records {
		if rec.CHROM == chrom && rec.POS >= start && rec.POS < stop {
			qtlRecords = append(qtlRecords, rec)
		}
	}
	fmt.Printf("QTL records: %d\n", len(qtlRecords))

	fmt.Println("Getting resistant lines unique data frame ...")
	resLines := strings.Split(csResLines, ",")
	uniqueRecords := resistantLinesUnique(qtlRecords, allGTCols, resLines)
	fmt.Printf("Unique records: %d\n", len(uniqueRecords))

	fmt.Println("Getting gene start and stop positions ...")
	geneDic, err := geneStartStop(gffPath, chrom, start, stop)
	if err != nil {
		return nil, err
	}

	var rows []GeneSpaceRow

	for gene, gr := range geneDic {
		fmt.Printf("Processing gene: %s\n", gene)
		gStart, _ := strconv.Atoi(gr.Start)
		gStop, _ := strconv.Atoi(gr.Stop)

		// SNPs within gene
		var geneRecords []VCFRecord
		for _, rec := range qtlRecords {
			if rec.POS >= gStart && rec.POS <= gStop {
				geneRecords = append(geneRecords, rec)
			}
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
		for _, rec := range uniqueRecords {
			if rec.POS >= gStart && rec.POS <= gStop {
				uniqueGene = append(uniqueGene, rec)
			}
		}
		// Unique significant SNPs (coding / non-synonymous)
		resSig := 0
		for _, rec := range uniqueGene {
			if !nonCodingEffects[rec.SNPEffEffect] {
				resSig++
			}
		}

		stats := prgStats(species, gene)
		desc := geneDesc(species, gene)

		fmt.Printf("  #SNPs: %d  #SigSNPs: %d  #ResUnique: %d  #ResUniqueSig: %d\n",
			len(geneRecords), numSig, len(uniqueGene), resSig)
		fmt.Printf("  Gene description: %s\n", desc)
		fmt.Printf("  %%PRG: %v\n", stats)

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
	}

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
	fmt.Printf("Output written to %s\n", outFile)

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
//		fmt.Println("USAGE: gene_space <gff> <super_vcf> <chrom> <start> <stop> <comma-separated-resistant-lines> <species>")
//		os.Exit(1)
//	}
//
//	gff := os.Args[1]
//	superVCF := os.Args[2]
//	chrom := os.Args[3]
//	startStr := os.Args[4]
//	stopStr := os.Args[5]
//	csResLines := os.Args[6]
//	species := os.Args[7]
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
//	if _, err := geneSpace(gff, superVCF, chrom, start, stop, csResLines, species); err != nil {
//		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
//		os.Exit(1)
//	}
//}
