// Package markerdesign designs markers inside QTLs: it slices a super VCF table
// to each QTL, simplifies genotypes to one-letter codes, keeps sites that split
// resistant from susceptible lines and adds wobbled flanking sequence.
// Port of utils/marker_design_pipeline_fast.py.
package markerdesign

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"slices"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/xuri/excelize/v2"
	"golang.org/x/sync/errgroup"
)

// Keys are the sorted, distinct one-letter codes called at a site.
var iupac = map[string]string{
	"AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K", "+-": "0",
}

const hets = "MRWSYK0"

var wobbleCodes = map[string]string{
	"AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K", "+-": "0",
	"CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N",
}

var trailingDigits = regexp.MustCompile(`(\d+)$`)

type qtl struct {
	chrom  string
	lo, hi int
}

type variant struct {
	pos      int
	ref, alt string
}

// setKey returns the distinct bytes of b, sorted, for the iupac/wobble lookups.
func setKey(b []byte) string {
	slices.Sort(b)
	out := b[:0]
	for i, c := range b {
		if i == 0 || c != b[i-1] {
			out = append(out, c)
		}
	}
	return string(out)
}

// openText opens a plain or gzip/bgzip-compressed text file.
func openText(path string) (io.Reader, func() error, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	if !strings.HasSuffix(path, ".gz") {
		return f, f.Close, nil
	}
	gz, err := gzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, nil, err
	}
	return gz, func() error { gz.Close(); return f.Close() }, nil
}

// readSheet returns the first sheet of an Excel file as one map per row, keyed by header.
func readSheet(path string) ([]map[string]string, error) {
	f, err := excelize.OpenFile(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	rows, err := f.GetRows(f.GetSheetName(0), excelize.Options{RawCellValue: true})
	if err != nil || len(rows) == 0 {
		return nil, err
	}
	var out []map[string]string
	for _, r := range rows[1:] {
		m := make(map[string]string, len(rows[0]))
		for i, h := range rows[0] {
			if i < len(r) {
				m[h] = r[i]
			}
		}
		out = append(out, m)
	}
	return out, nil
}

// excelInt parses an Excel numeric cell, which may be stored as "1234" or "1234.0".
func excelInt(s string) (int, bool) {
	v, err := strconv.ParseFloat(s, 64)
	return int(v), err == nil
}

// splitTSV splits a line on tabs, falling back to a CSV reader so quoted fields
// are unquoted the way pandas did.
func splitTSV(line string) ([]string, error) {
	if !strings.Contains(line, `"`) {
		return strings.Split(line, "\t"), nil
	}
	r := csv.NewReader(strings.NewReader(line))
	r.Comma, r.LazyQuotes, r.FieldsPerRecord = '\t', true, -1
	return r.Read()
}

// readSuperVcf streams the super VCF table and keeps only rows inside a QTL window.
func readSuperVcf(path string, qtls []qtl) ([]string, [][]string, error) {
	in, closeFn, err := openText(path)
	if err != nil {
		return nil, nil, err
	}
	defer closeFn()
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 1<<20), 1<<30)
	if !sc.Scan() {
		return nil, nil, fmt.Errorf("%s: empty file", path)
	}
	header, err := splitTSV(sc.Text())
	if err != nil {
		return nil, nil, err
	}
	if header[0] != "CHROM" || header[1] != "POS" {
		return nil, nil, fmt.Errorf("%s: first two columns must be CHROM and POS", path)
	}
	var rows [][]string
	for sc.Scan() {
		line := sc.Text()
		chrom, rest, _ := strings.Cut(line, "\t")
		posStr, _, _ := strings.Cut(rest, "\t")
		pos, err := strconv.Atoi(posStr)
		if err != nil {
			continue
		}
		for _, q := range qtls {
			if chrom == q.chrom && pos >= q.lo && pos <= q.hi {
				row, err := splitTSV(line)
				if err != nil {
					return nil, nil, err
				}
				for len(row) < len(header) {
					row = append(row, "")
				}
				rows = append(rows, row)
				break
			}
		}
	}
	return header, rows, sc.Err()
}

// readGenome loads only the wanted chromosomes from a FASTA.
func readGenome(path string, chroms map[string]bool) (map[string]string, error) {
	in, closeFn, err := openText(path)
	if err != nil {
		return nil, err
	}
	defer closeFn()
	genome := make(map[string]string)
	var name string
	var seq strings.Builder
	flush := func() {
		if chroms[name] {
			genome[name] = seq.String()
		}
		seq.Reset()
	}
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 1<<20), 1<<30)
	for sc.Scan() {
		line := sc.Bytes()
		if len(line) > 0 && line[0] == '>' {
			flush()
			name = strings.Fields(string(line[1:]) + " ")[0]
		} else if chroms[name] {
			seq.Write(line)
		}
	}
	flush()
	return genome, sc.Err()
}

// simplify replaces each genotype with a one-character code per allele; '*' is a
// spanning deletion. It appends REF_CODE and INDEL_MULTILEN.
func simplify(header []string, rows [][]string, gtIdx []int, refI, altI int) ([]string, [][]string) {
	out := make([][]string, len(rows))
	for r, row := range rows {
		ref := row[refI]
		var alleles []string
		for _, a := range append([]string{ref}, strings.Split(row[altI], ",")...) {
			if a != "*" {
				alleles = append(alleles, a)
			}
		}
		lengths := map[int]bool{}
		minLen := math.MaxInt
		for _, a := range alleles {
			lengths[len(a)] = true
			minLen = min(minLen, len(a))
		}
		shortest := 0
		for _, a := range alleles {
			if len(a) == minLen {
				shortest++
			}
		}
		codes := map[string]byte{"*": '-', ".": 'N'}
		for _, a := range alleles {
			switch {
			case len(a) > minLen:
				codes[a] = '+'
			case len(a) == 1 && (len(lengths) == 1 || shortest > 1):
				codes[a] = a[0] // SNP allele (incl. SNP alleles at MIXED sites)
			default:
				codes[a] = '-'
			}
		}

		newRow := append(make([]string, 0, len(row)+2), row...)
		for _, c := range gtIdx {
			var called []byte
			for a := range strings.SplitSeq(strings.ReplaceAll(row[c], "|", "/"), "/") {
				code, ok := codes[a]
				if !ok {
					code = 'N'
				}
				called = append(called, code)
			}
			key := setKey(called)
			switch {
			case strings.Contains(key, "N"):
				newRow[c] = "N"
			case len(key) == 1:
				newRow[c] = key
			default:
				newRow[c] = iupac[key] // base/indel hets have no code
				if newRow[c] == "" {
					newRow[c] = "N"
				}
			}
		}
		multilen := "False"
		if len(lengths) > 2 {
			multilen = "True"
		}
		out[r] = append(newRow, string(codes[ref]), multilen)
	}
	return append(append([]string{}, header...), "REF_CODE", "INDEL_MULTILEN"), out
}

// splits reports whether the non-missing calls in a are one homozygous code that
// b, fully called, never carries (b may be het). splits(sus, res) is RES_NOT_SUS,
// splits(res, sus) is SUS_NOT_RES.
func splits(a, b []string) bool {
	call := ""
	for _, g := range a {
		if g == "N" {
			continue
		}
		if call != "" && g != call {
			return false
		}
		call = g
	}
	if call == "" || strings.Contains(hets, call) {
		return false
	}
	for _, g := range b {
		if g == "N" || g == call {
			return false
		}
	}
	return true
}

// wobble returns seq[lo:hi] (0-based) with IUPAC codes at substitution sites, the
// number of variants hit, and whether an indel overlaps.
func wobble(seq string, lo, hi int, vars []variant, maxRefLen int, target variant) (string, int, string) {
	hi = min(hi, len(seq))
	lo = min(lo, hi)
	flank := []byte(seq[lo:hi])
	nHits, hasIndel := 0, false
	i := sort.Search(len(vars), func(k int) bool { return vars[k].pos > lo-maxRefLen }) // earliest variant whose REF could reach the flank
	j := sort.Search(len(vars), func(k int) bool { return vars[k].pos > hi })
	for _, v := range vars[i:max(i, j)] {
		vLo := v.pos - 1
		var alleles []string
		for _, a := range append([]string{v.ref}, strings.Split(v.alt, ",")...) {
			if a != "*" {
				alleles = append(alleles, a)
			}
		}
		if vLo+len(v.ref) <= lo || len(alleles) < 2 || v == target {
			continue
		}
		nHits++
		indel := false
		for _, a := range alleles {
			indel = indel || len(a) != len(v.ref)
		}
		if indel {
			hasIndel = true
			continue
		}
		for k := range len(v.ref) { // SNPs and MNPs: wobble every overlapping base
			if lo <= vLo+k && vLo+k < hi {
				bases := make([]byte, len(alleles))
				for n, a := range alleles {
					bases[n] = strings.ToUpper(a[k : k+1])[0]
				}
				code := wobbleCodes[setKey(bases)]
				if code == "" {
					code = "N"
				}
				flank[vLo+k-lo] = code[0]
			}
		}
	}
	indel := "NO"
	if hasIndel {
		indel = "YES"
	}
	return string(flank), nHits, indel
}

var primerCols = []string{"FWD_PRIMER", "REV_PRIMER", "FWD_WOBBLED", "#FWD_WOBBLES", "FWD_INDEL WOBBLES",
	"REV_WOBBLED", "#REV_WOBBLES", "REV_INDEL WOBBLES", "UN_WOBBLED_SEQ", "WOBBLED_SEQ"}

// addPrimers adds flanking sequence around each row's variant, wobbled with every variant in vars.
func addPrimers(header []string, rows [][]string, col map[string]int, seq string, vars []variant, flank int) ([]string, [][]string) {
	if len(rows) == 0 {
		return header, rows
	}
	maxRefLen := 0
	for _, v := range vars {
		maxRefLen = max(maxRefLen, len(v.ref))
	}
	out := make([][]string, len(rows))
	for r, row := range rows {
		pos, _ := strconv.Atoi(row[col["POS"]])
		ref, alt := row[col["REF"]], row[col["ALT"]]
		target := variant{pos, ref, alt}

		// Move the VCF padding base(s) into the left flank: TA/T -> T[A/-]
		alleles := append([]string{ref}, strings.Split(alt, ",")...)
		var real []string
		for _, a := range alleles {
			if a != "*" {
				real = append(real, a)
			}
		}
		n := 0 // common prefix length
		if len(real) > 1 {
			n = len(real[0])
			for _, a := range real[1:] {
				k := 0
				for k < n && k < len(a) && a[k] == real[0][k] {
					k++
				}
				n = k
			}
		}
		trimmed := make([]string, len(alleles))
		for k, a := range alleles {
			switch {
			case a == "*":
				trimmed[k] = a
			case a[n:] == "":
				trimmed[k] = "-"
			default:
				trimmed[k] = a[n:]
			}
		}
		snp := "[" + strings.Join(trimmed, "/") + "]"
		start, end := pos-1+n, pos-1+len(ref) // 0-based span of the trimmed REF allele

		clip := func(i int) int { return max(0, min(i, len(seq))) }
		fwdSeq := seq[clip(start-flank):clip(start)]
		revSeq := seq[clip(end):clip(end+flank)]
		fwd, fwdN, fwdIndel := wobble(seq, max(0, start-flank), start, vars, maxRefLen, target)
		rev, revN, revIndel := wobble(seq, end, end+flank, vars, maxRefLen, target)
		out[r] = append(append(make([]string, 0, len(row)+len(primerCols)), row...),
			fwdSeq, revSeq, fwd, strconv.Itoa(fwdN), fwdIndel, rev, strconv.Itoa(revN), revIndel,
			fwdSeq+snp+revSeq, fwd+snp+rev)
	}
	return append(append([]string{}, header...), primerCols...), out
}

func writeTSV(path string, header []string, rows [][]string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	w := bufio.NewWriterSize(f, 1<<20)
	for _, row := range append([][]string{header}, rows...) {
		for i, v := range row {
			if i > 0 {
				w.WriteByte('\t')
			}
			if strings.ContainsAny(v, "\"\t\n\r") {
				v = `"` + strings.ReplaceAll(v, `"`, `""`) + `"`
			}
			w.WriteString(v)
		}
		w.WriteByte('\n')
	}
	if err := w.Flush(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// writeColouredExcel writes the table with conditional formatting on the GT
// columns: het orange, REF lightblue, other calls yellow.
func writeColouredExcel(path string, header []string, rows [][]string) error {
	f := excelize.NewFile()
	defer f.Close()
	const sheet = "Sheet1"

	var gtCols []string
	refCol := ""
	isGT := make([]bool, len(header))
	for i, h := range header {
		name, _ := excelize.ColumnNumberToName(i + 1)
		if strings.HasSuffix(h, ".GT") {
			gtCols, isGT[i] = append(gtCols, name), true
		}
		if h == "REF_CODE" {
			refCol = name
		}
	}
	if len(rows) > 0 && len(gtCols) > 0 && refCol != "" {
		last := len(rows) + 1
		ranges := make([]string, len(gtCols))
		for i, c := range gtCols {
			ranges[i] = fmt.Sprintf("%s2:%s%d", c, c, last)
		}
		cell := gtCols[0] + "2"
		het := make([]string, len(hets))
		for i, h := range hets {
			het[i] = fmt.Sprintf(`%s="%c"`, cell, h)
		}
		var opts []excelize.ConditionalFormatOptions
		for _, rule := range [][2]string{
			{"OR(" + strings.Join(het, ",") + ")", "FFA500"},
			{fmt.Sprintf("EXACT(%s,$%s2)", cell, refCol), "ADD8E6"},
			{cell + `<>"N"`, "FFFF00"},
		} {
			style, err := f.NewConditionalStyle(&excelize.Style{Fill: excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{rule[1]}}})
			if err != nil {
				return err
			}
			opts = append(opts, excelize.ConditionalFormatOptions{Type: "formula", Criteria: rule[0], Format: &style})
		}
		// Must be set before the stream writer takes over the sheet
		if err := f.SetConditionalFormat(sheet, strings.Join(ranges, " "), opts); err != nil {
			return err
		}
	}

	// Columns where every value is a number are written as numbers, as pandas did
	numeric := make([]bool, len(header))
	for i := range header {
		numeric[i] = !isGT[i]
		for _, row := range rows {
			if row[i] == "" {
				continue
			}
			if v, err := strconv.ParseFloat(row[i], 64); err != nil || math.IsInf(v, 0) || math.IsNaN(v) {
				numeric[i] = false
				break
			}
		}
	}

	sw, err := f.NewStreamWriter(sheet)
	if err != nil {
		return err
	}
	cells := make([]any, len(header))
	for i, h := range header {
		cells[i] = h
	}
	if err := sw.SetRow("A1", cells); err != nil {
		return err
	}
	for r, row := range rows {
		cells := make([]any, len(header))
		for i, v := range row {
			switch {
			case v == "":
				cells[i] = nil
			case numeric[i]:
				cells[i], _ = strconv.ParseFloat(v, 64)
			case header[i] == "INDEL_MULTILEN":
				cells[i] = v == "True"
			default:
				cells[i] = v
			}
		}
		cellName, _ := excelize.CoordinatesToCellName(1, r+2)
		if err := sw.SetRow(cellName, cells); err != nil {
			return err
		}
	}
	if err := sw.Flush(); err != nil {
		return err
	}
	return f.SaveAs(path)
}

// MarkerDesignPipeline runs the marker design steps for every QTL in qtlsFile,
// writing per-QTL TSV and coloured Excel files to outDir.
func MarkerDesignPipeline(superVcf, qtlsFile string, flankSize int, resLines, susLines []string,
	genomeFasta, gvcf, markerDbExcel, crop, trait, outDir string, threads int) error {
	t0 := time.Now()
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		return err
	}

	fmt.Println("Read qtl file")
	qtlRows, err := readSheet(qtlsFile)
	if err != nil {
		return fmt.Errorf("reading QTL file: %w", err)
	}
	var qtls []qtl
	chroms := map[string]bool{}
	for _, r := range qtlRows {
		start, ok1 := excelInt(r["Start"])
		stop, ok2 := excelInt(r["Stop"])
		if r["Chrom"] == "" || !ok1 || !ok2 {
			continue
		}
		qtls = append(qtls, qtl{r["Chrom"], max(0, start-flankSize), stop + flankSize})
		chroms[r["Chrom"]] = true
	}

	dbRows, err := readSheet(markerDbExcel)
	if err != nil {
		return fmt.Errorf("reading marker DB: %w", err)
	}
	type dbEntry struct {
		chrom string
		pos   int
		name  string
	}
	var markerDb []dbEntry
	dbNames, dbTypes := map[string][]string{}, map[string][]string{}
	for _, r := range dbRows {
		pos, ok := excelInt(r["pos"])
		if r["chrom"] == "" || !ok {
			continue
		}
		markerDb = append(markerDb, dbEntry{r["chrom"], pos, r["marker_name"]})
		key := r["chrom"] + "\t" + strconv.Itoa(pos)
		dbNames[key] = append(dbNames[key], r["marker_name"])
		if !slices.Contains(dbTypes[key], r["typ"]) {
			dbTypes[key] = append(dbTypes[key], r["typ"])
		}
	}

	fmt.Println("Reading super vcf (QTL regions only) ...")
	header, rows, err := readSuperVcf(superVcf, qtls)
	if err != nil {
		return fmt.Errorf("reading super vcf: %w", err)
	}
	fmt.Printf("%d super vcf sites fall in QTLs\n", len(rows))

	fmt.Println("Reading genome ...")
	genome, err := readGenome(genomeFasta, chroms)
	if err != nil {
		return fmt.Errorf("reading genome: %w", err)
	}

	// Column order: MARKER_NAME, SOURCE, site columns, res GTs, sus GTs, other GTs, rest
	idx := map[string]int{}
	for i, h := range header {
		idx[h] = i
	}
	var resGts, susGts, otherGts []string
	for _, l := range resLines {
		resGts = append(resGts, strings.TrimSpace(l)+".GT")
	}
	for _, l := range susLines {
		susGts = append(susGts, strings.TrimSpace(l)+".GT")
	}
	order := append([]string{"CHROM", "POS", "REF", "ALT", "QUAL", "TYPE"}, append(resGts, susGts...)...)
	inOrder := map[string]bool{"MARKER_NAME": true, "SOURCE": true}
	for _, c := range order {
		if _, ok := idx[c]; !ok {
			return fmt.Errorf("super vcf has no %s column", c)
		}
		inOrder[c] = true
	}
	for _, h := range header {
		if strings.HasSuffix(h, ".GT") && !inOrder[h] {
			otherGts = append(otherGts, h)
		}
	}
	order = append(order, otherGts...)
	for _, c := range otherGts {
		inOrder[c] = true
	}
	for _, h := range header {
		if !inOrder[h] {
			order = append(order, h)
		}
	}
	header = append([]string{"MARKER_NAME", "SOURCE"}, order...)

	// Marker name/source from the DB at the same CHROM/POS, else a NEW name
	for r, row := range rows {
		chrom, posStr := row[0], row[1]
		pos, _ := strconv.Atoi(posStr)
		key := chrom + "\t" + strconv.Itoa(pos)
		name, source := strings.Join(dbNames[key], ";"), strings.Join(dbTypes[key], ";")
		if _, ok := dbNames[key]; !ok {
			num := chrom
			if m := trailingDigits.FindString(chrom); m != "" {
				n, _ := strconv.Atoi(m) // Cma_Chr01 -> 1
				num = strconv.Itoa(n)
			}
			name, source = fmt.Sprintf("%s_%s_%.3f_%s", crop, num, float64(pos)/1e6, trait), "NEW"
		}
		newRow := make([]string, 0, len(header))
		newRow = append(newRow, name, source)
		for _, c := range order {
			newRow = append(newRow, row[idx[c]])
		}
		rows[r] = newRow
	}

	col := map[string]int{}
	var gtIdx, resIdx, susIdx []int
	for i, h := range header {
		col[h] = i
		if strings.HasSuffix(h, ".GT") {
			gtIdx = append(gtIdx, i)
		}
	}
	for _, g := range resGts {
		resIdx = append(resIdx, col[g])
	}
	for _, g := range susGts {
		susIdx = append(susIdx, col[g])
	}
	pick := func(row []string, cols []int) []string {
		out := make([]string, len(cols))
		for i, c := range cols {
			out[i] = row[c]
		}
		return out
	}

	fmt.Println("extract QTLs from super vcf and analyse ...")
	var g errgroup.Group
	g.SetLimit(max(1, threads))
	for _, q := range qtls {
		g.Go(func() error {
			tag := fmt.Sprintf("[%s:%d-%d]", q.chrom, q.lo, q.hi)
			var qtlRows [][]string
			for _, row := range rows {
				if pos, _ := strconv.Atoi(row[col["POS"]]); row[col["CHROM"]] == q.chrom && pos >= q.lo && pos <= q.hi {
					qtlRows = append(qtlRows, row)
				}
			}
			base := filepath.Join(outDir, fmt.Sprintf("%s_%d_%d.SUPERVCF", q.chrom, q.lo, q.hi))
			if err := writeTSV(base+".tsv", header, qtlRows); err != nil {
				return err
			}
			fmt.Printf("%s QTL tsv %s.tsv created ...\n", tag, base)

			simpleHeader, simpleRows := simplify(header, qtlRows, gtIdx, col["REF"], col["ALT"])
			multilen := 0
			for _, row := range simpleRows {
				if row[len(row)-1] == "True" {
					multilen++
				}
			}
			if err := writeTSV(base+".FORMATED.tsv", simpleHeader, simpleRows); err != nil {
				return err
			}
			fmt.Printf("%s Reformated %s.FORMATED.tsv: %d sites, %d with >2 indel lengths\n", tag, base, len(simpleRows), multilen)
			if err := writeColouredExcel(base+".FORMATED.COLOURED.xlsx", simpleHeader, simpleRows); err != nil {
				return err
			}

			var resNotSusRows, susNotResRows [][]string
			for _, row := range simpleRows {
				sus, res := pick(row, susIdx), pick(row, resIdx)
				if splits(sus, res) {
					resNotSusRows = append(resNotSusRows, row)
				}
				if splits(res, sus) {
					susNotResRows = append(susNotResRows, row)
				}
			}

			fmt.Printf("%s Adding primers, wobbled with all variants in %s ...\n", tag, gvcf)
			region := fmt.Sprintf("%s:%d-%d", q.chrom, max(1, q.lo-1000), q.hi+1000) // margin for flanks and long deletions
			var stderr strings.Builder
			bcf := exec.Command("bcftools", "query", "-r", region, "-f", `%POS\t%REF\t%ALT\n`, gvcf)
			bcf.Stderr = &stderr
			query, err := bcf.Output()
			if err != nil {
				return fmt.Errorf("%s bcftools query: %w: %s", tag, err, stderr.String())
			}
			var vars []variant
			for line := range strings.SplitSeq(strings.TrimSpace(string(query)), "\n") {
				f := strings.Split(line, "\t")
				if len(f) != 3 || f[2] == "." {
					continue
				}
				pos, _ := strconv.Atoi(f[0])
				vars = append(vars, variant{pos, f[1], f[2]})
			}
			sort.SliceStable(vars, func(i, j int) bool { return vars[i].pos < vars[j].pos })

			seq := genome[q.chrom]
			for _, out := range []struct {
				suffix string
				rows   [][]string
			}{{".FORMATED.RES_NOT_SUS", resNotSusRows}, {".FORMATED.SUS_NOT_RES", susNotResRows}} {
				h, r := addPrimers(simpleHeader, out.rows, col, seq, vars, 50)
				if err := writeTSV(base+out.suffix+".tsv", h, r); err != nil {
					return err
				}
				if err := writeColouredExcel(base+out.suffix+".COLOURED.xlsx", h, r); err != nil {
					return err
				}
			}

			placedPos := map[int]bool{}
			for _, row := range qtlRows {
				pos, _ := strconv.Atoi(row[col["POS"]])
				placedPos[pos] = true
			}
			inQtl, placed := 0, 0
			var missed []string
			for _, m := range markerDb {
				if m.chrom != q.chrom || m.pos < q.lo || m.pos > q.hi {
					continue
				}
				inQtl++
				if placedPos[m.pos] {
					placed++
				} else {
					missed = append(missed, m.name)
				}
			}
			fmt.Printf("SUMMARY %s:%d-%d: %d markers already in QTL, %d placed, %d did not (%s)\n",
				q.chrom, q.lo, q.hi, inQtl, placed, len(missed), strings.Join(missed, ", "))
			return nil
		})
	}
	if err := g.Wait(); err != nil {
		return err
	}
	fmt.Printf("Finished in %s\n", time.Since(t0).Round(time.Second))
	return nil
}
