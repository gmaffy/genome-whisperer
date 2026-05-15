package gene_space

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/fatih/color"
)

type GeneDetails struct {
	Gene  string
	Start int
	Stop  int
}

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
			gene := strings.SplitN(strings.SplitN(fields[8], ";", 2)[0], "=", 2)[1]
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

func genePrgMap(prgFile string) (map[string][]PrgMatches, error) {
	dic := make(map[string][]PrgMatches)
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
		return nil, fmt.Errorf("reading header: %w", err) //nil //, nil, fmt.Errorf("reading header: %w", err)
	}

	// Index key columns
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

	//var records []PrgMatches
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

		dic[row[seqIdIdx]] = append(dic[row[seqIdIdx]], rec)
	}
	return dic, nil

}

type VCFRecord struct {
	CHROM        string
	POS          int
	SNPEffEffect string
	Genotypes    map[string]string // key: "<sample>.GT"
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

func GeneSpace(gffPath, vcfTable, chrom string, start, stop int, resLines, susLines []string, descFile, prgFile string) error {

	color.Green("▶ Starting GeneSpace analysis for chromosome %s (%d - %d)", chrom, start, stop)
	qtlRecords, gtCols, err := getQtlRecords(vcfTable)
	if err != nil {
		return fmt.Errorf("reading VCF table: %w", err)
	}
	gencPos, err := genePosFromGffMap(gffPath, chrom, start, stop)
	if err != nil {
		return fmt.Errorf("reading GFF file: %w", err)
	}
	geneDesc, err := geneDescMap(descFile)
	if err != nil {
		return fmt.Errorf("reading gene description file: %w", err)
	}
	genePrg, err := genePrgMap(prgFile)
	if err != nil {
		return fmt.Errorf("reading PRG blast file: %w", err)
	}

	for gene, pos := range gencPos {
		color.Green("Working on gene %s (%s)", gene, pos[0])
		for _, rec := range qtlRecords {
			
		}
	}

	return nil

}
