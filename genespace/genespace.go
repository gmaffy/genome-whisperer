// Package genespace ranks the genes in a QTL window as candidates for a trait.
//
// Each gene is scored on independent axes and the tier is derived from them:
//
//	Variant score  (0-4)  the most damaging variant whose genotypes follow the
//	                      resistant/susceptible split (dominant or recessive)
//	Function score (0-3)  PRG hit quality plus description terms from a
//	                      built-in trait profile
//
// Genetic evidence comes first; annotation adjusts the ranking within it.
package genespace

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/xuri/excelize/v2"
)

// PRG hit thresholds for a hit to count as a real PRG match. Weaker hits still
// earn a function score of 1.
const (
	PrgMinIdentity = 40.0  // percent identity
	PrgMinCoverage = 60.0  // alignment length as a percent of query length
	PrgMaxEvalue   = 1e-10 // e-value
)

// Options configures a GeneSpace run.
type Options struct {
	GFF            string
	VariantTable   string // superVCF table (genome-whisperer or GoBSAseq)
	Chrom          string
	Start, Stop    int
	ResLines       []string
	SusLines       []string
	DescFile       string // optional: id<TAB>description
	PrgFile        string // optional: PRG top-hits table
	OutDir         string
	Trait          string  // built-in profile name
	MinConcordance float64 // 0 < x <= 1
	PromoterLength int     // bp upstream of a gene attributed to it; 0 disables
}

// ─── Genes ───────────────────────────────────────────────────────────────────

// Gene is one gene in the window, with its isoforms collapsed.
type Gene struct {
	ID          string
	Start, Stop int
	Strand      string
	Transcripts []string
}

func parseAttrs(s string) map[string]string {
	m := make(map[string]string)
	for _, kv := range strings.Split(s, ";") {
		k, v, ok := strings.Cut(strings.TrimSpace(kv), "=")
		if ok && k != "" {
			m[k] = v
		}
	}
	return m
}

// loadGenes reads the GFF and returns the genes on chrom that overlap
// [start, stop], sorted by position. Transcripts (mRNA/transcript rows) are
// grouped by their Parent gene; if the GFF has no transcript rows, gene rows
// are used. Malformed lines are skipped.
func loadGenes(gff, chrom string, start, stop int) ([]Gene, error) {
	f, err := os.Open(gff)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	fromTx := make(map[string]*Gene)
	fromGene := make(map[string]*Gene)

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 9 || fields[0] != chrom {
			continue
		}
		s, err1 := strconv.Atoi(fields[3])
		e, err2 := strconv.Atoi(fields[4])
		if err1 != nil || err2 != nil {
			continue
		}
		attrs := parseAttrs(fields[8])

		switch fields[2] {
		case "mRNA", "transcript":
			id := attrs["ID"]
			parent, _, _ := strings.Cut(attrs["Parent"], ",")
			if parent == "" {
				parent = id
			}
			if parent == "" {
				continue
			}
			g, ok := fromTx[parent]
			if !ok {
				g = &Gene{ID: parent, Start: s, Stop: e, Strand: fields[6]}
				fromTx[parent] = g
			}
			g.Start = min(g.Start, s)
			g.Stop = max(g.Stop, e)
			if id != "" {
				g.Transcripts = append(g.Transcripts, id)
			}
		case "gene":
			if id := attrs["ID"]; id != "" {
				fromGene[id] = &Gene{ID: id, Start: s, Stop: e, Strand: fields[6]}
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	src := fromTx
	if len(src) == 0 {
		src = fromGene
	}
	var genes []Gene
	for _, g := range src {
		if g.Start <= stop && g.Stop >= start {
			genes = append(genes, *g)
		}
	}
	sort.Slice(genes, func(i, j int) bool {
		if genes[i].Start != genes[j].Start {
			return genes[i].Start < genes[j].Start
		}
		return genes[i].ID < genes[j].ID
	})
	return genes, nil
}

// ─── Descriptions and PRG hits ───────────────────────────────────────────────

func geneDescMap(descFile string) (map[string]string, error) {
	dic := make(map[string]string)
	if descFile == "" {
		return dic, nil
	}
	f, err := os.Open(descFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	for scanner.Scan() {
		id, desc, ok := strings.Cut(scanner.Text(), "\t")
		if !ok {
			continue
		}
		desc, _, _ = strings.Cut(desc, "\t")
		dic[id] = strings.TrimSpace(desc)
	}
	return dic, scanner.Err()
}

// PrgHit is the best PRGdb BLAST hit for one protein.
type PrgHit struct {
	PercIdent float64
	Coverage  float64 // Length / Qlen, percent
	Length    int
	Qlen      int
	Eval      float64
}

// Pass reports whether the hit clears all three PRG thresholds.
func (h PrgHit) Pass() bool {
	return h.PercIdent >= PrgMinIdentity && h.Coverage >= PrgMinCoverage && h.Eval <= PrgMaxEvalue
}

func genePrgMap(prgFile string) (map[string]PrgHit, error) {
	dic := make(map[string]PrgHit)
	if prgFile == "" {
		return dic, nil
	}
	f, err := os.Open(prgFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.LazyQuotes = true
	r.FieldsPerRecord = -1

	header, err := r.Read()
	if err != nil {
		return nil, fmt.Errorf("reading header: %w", err)
	}
	idx := make(map[string]int)
	for i, h := range header {
		idx[h] = i
	}
	cols := []string{"QseqID", "PercIdent", "Length", "Qlen", "Eval"}
	for _, c := range cols {
		if _, ok := idx[c]; !ok {
			return nil, fmt.Errorf("PRG table missing required columns (%s)", strings.Join(cols, ", "))
		}
	}
	maxIdx := 0
	for _, c := range cols {
		maxIdx = max(maxIdx, idx[c])
	}

	for {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading row: %w", err)
		}
		if len(row) <= maxIdx {
			continue
		}
		pi, _ := strconv.ParseFloat(row[idx["PercIdent"]], 64)
		ev, err := strconv.ParseFloat(row[idx["Eval"]], 64)
		if err != nil {
			continue
		}
		length, _ := strconv.Atoi(row[idx["Length"]])
		qlen, _ := strconv.Atoi(row[idx["Qlen"]])
		hit := PrgHit{PercIdent: pi, Length: length, Qlen: qlen, Eval: ev}
		if qlen > 0 {
			hit.Coverage = 100 * float64(length) / float64(qlen)
		}
		id := row[idx["QseqID"]]
		if old, ok := dic[id]; !ok || hit.Eval < old.Eval {
			dic[id] = hit
		}
	}
	return dic, nil
}

// ─── Variant table ───────────────────────────────────────────────────────────

// Variant is one row of the superVCF table.
type Variant struct {
	POS       int
	Effect    string
	Impact    string
	Genotypes map[string]string // "<sample>.GT" -> genotype
}

// readVariants reads the variants on chrom within [start, stop] from a superVCF
// table, sorted by position. Only CHROM, POS and the *.GT columns are required;
// SNPEFF_EFFECT and SNPEFF_IMPACT are used when present (matched
// case-insensitively, since genome-whisperer writes "SNPEff_IMPACT").
func readVariants(path, chrom string, start, stop int) ([]Variant, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, fmt.Errorf("opening variant table: %w", err)
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.LazyQuotes = true
	r.FieldsPerRecord = -1
	r.ReuseRecord = true

	header, err := r.Read()
	if err != nil {
		return nil, nil, fmt.Errorf("reading header: %w", err)
	}
	header = slices.Clone(header)

	col := func(name string) int {
		for i, h := range header {
			if strings.EqualFold(strings.TrimSpace(h), name) {
				return i
			}
		}
		return -1
	}
	chromIdx, posIdx := col("CHROM"), col("POS")
	if chromIdx < 0 || posIdx < 0 {
		return nil, nil, fmt.Errorf("variant table missing required columns (CHROM, POS)")
	}
	effIdx, impIdx := col("SNPEFF_EFFECT"), col("SNPEFF_IMPACT")

	var gtCols []string
	var gtIdx []int
	for i, h := range header {
		if strings.HasSuffix(h, ".GT") {
			gtCols = append(gtCols, h)
			gtIdx = append(gtIdx, i)
		}
	}
	if len(gtCols) == 0 {
		return nil, nil, fmt.Errorf("variant table has no <sample>.GT columns")
	}

	field := func(row []string, i int) string {
		if i < 0 || i >= len(row) {
			return ""
		}
		return row[i]
	}

	var out []Variant
	for {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("reading row: %w", err)
		}
		if field(row, chromIdx) != chrom {
			continue
		}
		pos, err := strconv.Atoi(field(row, posIdx))
		if err != nil || pos < start || pos > stop {
			continue
		}
		v := Variant{
			POS:       pos,
			Effect:    field(row, effIdx),
			Impact:    field(row, impIdx),
			Genotypes: make(map[string]string, len(gtCols)),
		}
		for k, c := range gtCols {
			v.Genotypes[c] = field(row, gtIdx[k])
		}
		out = append(out, v)
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i].POS < out[j].POS })
	return out, gtCols, nil
}

// sampleGTCols maps line names (e.g. "R1", or "R1.GT") to their genotype
// column headers ("R1.GT"). Names are trimmed and blanks skipped; a name with
// no matching column in the variant table is an error, so a typo can't
// silently yield zero concordant SNPs.
func sampleGTCols(lines, gtCols []string) ([]string, error) {
	var cols, missing []string
	for _, l := range lines {
		l = strings.TrimSpace(l)
		if l == "" {
			continue
		}
		col := strings.TrimSuffix(l, ".GT") + ".GT"
		if !slices.Contains(gtCols, col) {
			missing = append(missing, l)
			continue
		}
		cols = append(cols, col)
	}
	if len(missing) > 0 {
		return nil, fmt.Errorf("no genotype column for %s (available: %s)",
			strings.Join(missing, ", "), strings.Join(gtCols, ", "))
	}
	if len(cols) == 0 {
		return nil, fmt.Errorf("no line names given")
	}
	return cols, nil
}

// ─── Genotype concordance ────────────────────────────────────────────────────

// normGT parses a diploid genotype into sorted alleles. Phasing is dropped,
// "1/0" becomes "0/1", a haploid call "1" becomes "1/1", and any missing
// allele makes the whole call missing (ok == false).
func normGT(gt string) (a, b string, ok bool) {
	gt = strings.ReplaceAll(strings.TrimSpace(gt), "|", "/")
	parts := strings.Split(gt, "/")
	switch len(parts) {
	case 1:
		a, b = parts[0], parts[0]
	case 2:
		a, b = parts[0], parts[1]
	default:
		return "", "", false
	}
	if a == "" || a == "." || b == "" || b == "." {
		return "", "", false
	}
	if a > b {
		a, b = b, a
	}
	return a, b, true
}

// SNPTest is the result of testing one variant against the R/S split.
type SNPTest struct {
	Concordance float64 // 0-1, balanced between the R and S groups
	Mode        string  // "DOM", "REC", "DOM/REC" or "" when untestable
	Allele      string  // the allele that marks the resistant group
	Testable    bool
}

// testSNP tests a variant under both inheritance modes, for every allele seen
// (the reference allele included, since the reference genome may itself be
// resistant):
//
//	DOM: resistant lines carry the allele, susceptible lines don't.
//	REC: resistant lines are homozygous for it, susceptible lines aren't.
//
// Concordance is the mean of the R and S match fractions, so a large group
// can't hide a mismatch in a small one. Missing calls are ignored, but the
// variant is untestable unless both groups have a call and at least half of
// all lines are called.
func testSNP(resGTs, susGTs []string) SNPTest {
	type call struct{ a, b string }
	parse := func(gts []string) []call {
		var out []call
		for _, g := range gts {
			if a, b, ok := normGT(g); ok {
				out = append(out, call{a, b})
			}
		}
		return out
	}
	rc, sc := parse(resGTs), parse(susGTs)
	if len(rc) == 0 || len(sc) == 0 || 2*(len(rc)+len(sc)) < len(resGTs)+len(susGTs) {
		return SNPTest{}
	}

	alleleSet := make(map[string]struct{})
	for _, c := range append(slices.Clone(rc), sc...) {
		alleleSet[c.a] = struct{}{}
		alleleSet[c.b] = struct{}{}
	}
	alleles := make([]string, 0, len(alleleSet))
	for a := range alleleSet {
		alleles = append(alleles, a)
	}
	sort.Strings(alleles)

	frac := func(calls []call, pred func(call) bool) float64 {
		n := 0
		for _, c := range calls {
			if pred(c) {
				n++
			}
		}
		return float64(n) / float64(len(calls))
	}

	best := SNPTest{Testable: true}
	const eps = 1e-9
	for _, al := range alleles {
		carries := func(c call) bool { return c.a == al || c.b == al }
		hom := func(c call) bool { return c.a == al && c.b == al }

		dom := (frac(rc, carries) + frac(sc, func(c call) bool { return !carries(c) })) / 2
		rec := (frac(rc, hom) + frac(sc, func(c call) bool { return !hom(c) })) / 2

		score, mode := dom, "DOM"
		switch {
		case math.Abs(dom-rec) < eps:
			mode = "DOM/REC"
		case rec > dom:
			score, mode = rec, "REC"
		}
		if score > best.Concordance+eps {
			best.Concordance, best.Mode, best.Allele = score, mode, al
		}
	}
	return best
}

// ─── Scoring ─────────────────────────────────────────────────────────────────

var (
	highEffects = []string{
		"STOP_GAINED", "STOP_LOST", "START_LOST", "FRAME_SHIFT",
		"SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "EXON_DELETED", "RARE_AMINO_ACID",
	}
	moderateEffects = []string{
		"NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", "CODON_CHANGE",
		"CODON_INSERTION", "CODON_DELETION",
		"CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_CHANGE_PLUS_CODON_DELETION",
		"UTR_5_DELETED", "UTR_3_DELETED", "SPLICE_SITE_BRANCH_U12",
	}
	regulatoryEffects = []string{
		"UTR_5_PRIME", "UTR_3_PRIME", "SPLICE_SITE_REGION", "START_GAINED", "NON_SYNONYMOUS_STOP",
	}
	silentEffects = []string{"SYNONYMOUS_CODING", "SYNONYMOUS_STOP", "SYNONYMOUS_START"}
)

// impactScore scores one concordant variant:
//
//	4 HIGH · 3 MODERATE · 2 UTR / splice region / other LOW · 1 upstream · 0 other
//
// Synonymous changes score 0 even though snpEff calls them LOW. A variant that
// only reaches the gene through its promoter window scores at most 1. When the
// table has no impact column, the class is derived from the effect name.
func impactScore(impact, effect string, upstream bool) int {
	eff := strings.ToUpper(strings.TrimSpace(effect))
	if upstream {
		return 1
	}
	if slices.Contains(silentEffects, eff) {
		return 0
	}
	switch strings.ToUpper(strings.TrimSpace(impact)) {
	case "HIGH":
		return 4
	case "MODERATE":
		return 3
	case "LOW":
		return 2
	}
	switch {
	case slices.Contains(highEffects, eff):
		return 4
	case slices.Contains(moderateEffects, eff):
		return 3
	case slices.Contains(regulatoryEffects, eff):
		return 2
	case eff == "UPSTREAM":
		return 1
	}
	return 0
}

// Tier is a candidate class, A (best) to D.
type Tier struct {
	Rank  int // 1 = A … 4 = D, for sorting
	Code  string
	Label string
	Hex   string // row fill colour
}

var (
	TierA = Tier{1, "A", "Lead candidate", "FF9900"}
	TierB = Tier{2, "B", "Strong", "FFD966"}
	TierC = Tier{3, "C", "Plausible", "9DC3E6"}
	TierD = Tier{4, "D", "Positional only", "FFFFFF"}
)

// ClassifyTier combines the scores:
//
//	A  variant >= 3 and function >= 2
//	B  variant >= 2, or (any concordant variant and function >= 2)
//	C  function >= 2
//	D  everything else
func ClassifyTier(variantScore, functionScore, numConcordant int) Tier {
	switch {
	case variantScore >= 3 && functionScore >= 2:
		return TierA
	case variantScore >= 2 || (numConcordant > 0 && functionScore >= 2):
		return TierB
	case functionScore >= 2:
		return TierC
	default:
		return TierD
	}
}

// functionScore scores how plausible a gene is for the trait:
//
//	3 PRG pass and strong term · 2 PRG pass or strong term
//	1 weak term, weak PRG hit or TF · 0 nothing
func functionScore(prg *PrgHit, strong, weak, tf []string) int {
	prgPass := prg != nil && prg.Pass()
	switch {
	case prgPass && len(strong) > 0:
		return 3
	case prgPass || len(strong) > 0:
		return 2
	case prg != nil || len(weak) > 0 || len(tf) > 0:
		return 1
	}
	return 0
}

// GeneSpaceRow is one gene in the output.
type GeneSpaceRow struct {
	Gene        string
	Start, Stop int
	Strand      string

	NumVariants     int     // gene body + promoter window
	VariantsPerKb   float64 // gene body only
	NumConcordant   int
	BestConcordance float64 // over all testable variants
	BestMode        string  // mode of the top concordant variant
	TopEffect       string  // effect of the top concordant variant
	VariantScore    int

	PRG           *PrgHit // nil when no hit
	Description   string
	MatchedTerms  []string
	IsTF          bool
	FunctionScore int

	Tier Tier
}

// scoreGene attributes variants to g and scores it. variants must be sorted.
func scoreGene(g Gene, variants []Variant, resCols, susCols []string, minConc float64, promoter int,
	desc string, prg *PrgHit, prof Profile) GeneSpaceRow {

	lo, hi := g.Start, g.Stop
	switch g.Strand {
	case "+":
		lo -= promoter
	case "-":
		hi += promoter
	}
	i := sort.Search(len(variants), func(i int) bool { return variants[i].POS >= lo })
	j := sort.Search(len(variants), func(i int) bool { return variants[i].POS > hi })

	row := GeneSpaceRow{Gene: g.ID, Start: g.Start, Stop: g.Stop, Strand: g.Strand,
		Description: desc, PRG: prg}

	bodyVariants := 0
	bestConc := -1.0
	resGTs := make([]string, len(resCols))
	susGTs := make([]string, len(susCols))
	for _, v := range variants[i:j] {
		row.NumVariants++
		upstream := v.POS < g.Start || v.POS > g.Stop
		if !upstream {
			bodyVariants++
		}
		for k, c := range resCols {
			resGTs[k] = v.Genotypes[c]
		}
		for k, c := range susCols {
			susGTs[k] = v.Genotypes[c]
		}
		t := testSNP(resGTs, susGTs)
		if !t.Testable {
			continue
		}
		row.BestConcordance = max(row.BestConcordance, t.Concordance)
		if t.Concordance+1e-9 < minConc {
			continue
		}
		row.NumConcordant++
		s := impactScore(v.Impact, v.Effect, upstream)
		if s > row.VariantScore || (s == row.VariantScore && t.Concordance > bestConc) || row.NumConcordant == 1 {
			row.VariantScore = s
			bestConc = t.Concordance
			row.BestMode = t.Mode
			row.TopEffect = v.Effect
			if upstream {
				row.TopEffect = "UPSTREAM (promoter window)"
			}
		}
	}
	if kb := float64(g.Stop-g.Start+1) / 1000; kb > 0 {
		row.VariantsPerKb = math.Round(float64(bodyVariants)/kb*100) / 100
	}

	strong := matchTerms(prof.Strong, desc)
	weak := matchTerms(prof.Weak, desc)
	tf := matchTerms(prof.TF, desc)
	row.MatchedTerms = append(append(append([]string{}, strong...), weak...), tf...)
	row.IsTF = len(tf) > 0
	row.FunctionScore = functionScore(prg, strong, weak, tf)
	row.Tier = ClassifyTier(row.VariantScore, row.FunctionScore, row.NumConcordant)
	return row
}

// SortRows orders rows best-first: tier, variant score, function score, best
// concordance, number of concordant variants, then position.
func SortRows(rows []GeneSpaceRow) {
	sort.SliceStable(rows, func(i, j int) bool {
		a, b := rows[i], rows[j]
		switch {
		case a.Tier.Rank != b.Tier.Rank:
			return a.Tier.Rank < b.Tier.Rank
		case a.VariantScore != b.VariantScore:
			return a.VariantScore > b.VariantScore
		case a.FunctionScore != b.FunctionScore:
			return a.FunctionScore > b.FunctionScore
		case a.BestConcordance != b.BestConcordance:
			return a.BestConcordance > b.BestConcordance
		case a.NumConcordant != b.NumConcordant:
			return a.NumConcordant > b.NumConcordant
		case a.Start != b.Start:
			return a.Start < b.Start
		}
		return a.Gene < b.Gene
	})
}

// ─── Run ─────────────────────────────────────────────────────────────────────

func (o *Options) validate() error {
	switch {
	case o.GFF == "" || o.VariantTable == "" || o.Chrom == "":
		return fmt.Errorf("--gff, --variant and --chrom are required")
	case o.Start < 0 || o.Stop < o.Start:
		return fmt.Errorf("invalid window %d-%d", o.Start, o.Stop)
	case o.MinConcordance <= 0 || o.MinConcordance > 1:
		return fmt.Errorf("--min-concordance must be in (0, 1], got %g", o.MinConcordance)
	case o.PromoterLength < 0:
		return fmt.Errorf("--promoter-length must be >= 0, got %d", o.PromoterLength)
	}
	if o.OutDir == "" {
		o.OutDir = "."
	}
	if o.Trait == "" {
		o.Trait = DefaultProfile
	}
	return nil
}

// GeneSpace scores every gene in the window and writes the Excel report. It
// returns the rows (best first) and the path written.
func GeneSpace(opts Options) ([]GeneSpaceRow, string, error) {
	if err := opts.validate(); err != nil {
		return nil, "", err
	}
	prof, err := GetProfile(opts.Trait)
	if err != nil {
		return nil, "", err
	}
	if err := os.MkdirAll(opts.OutDir, 0755); err != nil {
		return nil, "", fmt.Errorf("creating output directory %s: %w", opts.OutDir, err)
	}
	if opts.DescFile == "" {
		color.Yellow("⚠ No --gene-description-tsv given: description terms will not be scored")
	}
	if opts.PrgFile == "" {
		color.Yellow("⚠ No --prg given: PRG hits will not be scored")
	}

	color.Green("▶ GeneSpace %s:%d-%d (trait profile: %s)", opts.Chrom, opts.Start, opts.Stop, prof.Name)

	// Load the reference files while the variant table loads.
	var (
		wg               sync.WaitGroup
		genes            []Gene
		descs            map[string]string
		prgs             map[string]PrgHit
		gErr, dErr, pErr error
	)
	wg.Add(3)
	go func() { defer wg.Done(); genes, gErr = loadGenes(opts.GFF, opts.Chrom, opts.Start, opts.Stop) }()
	go func() { defer wg.Done(); descs, dErr = geneDescMap(opts.DescFile) }()
	go func() { defer wg.Done(); prgs, pErr = genePrgMap(opts.PrgFile) }()

	// Read past the window edges so promoter windows of edge genes are complete.
	variants, gtCols, vErr := readVariants(opts.VariantTable, opts.Chrom,
		max(0, opts.Start-opts.PromoterLength), opts.Stop+opts.PromoterLength)
	wg.Wait()

	switch {
	case vErr != nil:
		return nil, "", fmt.Errorf("reading variant table: %w", vErr)
	case gErr != nil:
		return nil, "", fmt.Errorf("reading GFF file: %w", gErr)
	case dErr != nil:
		return nil, "", fmt.Errorf("reading gene description file: %w", dErr)
	case pErr != nil:
		return nil, "", fmt.Errorf("reading PRG file: %w", pErr)
	}

	resCols, err := sampleGTCols(opts.ResLines, gtCols)
	if err != nil {
		return nil, "", fmt.Errorf("--res-lines: %w", err)
	}
	susCols, err := sampleGTCols(opts.SusLines, gtCols)
	if err != nil {
		return nil, "", fmt.Errorf("--sus-lines: %w", err)
	}
	if len(genes) == 0 {
		return nil, "", fmt.Errorf("no genes found on %s between %d and %d in %s", opts.Chrom, opts.Start, opts.Stop, opts.GFF)
	}

	bar := utils.NewBar(int64(len(genes)), "Scoring genes")
	release := utils.AttachBar(bar)
	rows := make([]GeneSpaceRow, 0, len(genes))
	for _, g := range genes {
		ids := append([]string{g.ID}, g.Transcripts...)
		desc := ""
		for _, id := range ids {
			if d, ok := descs[id]; ok && d != "" {
				desc = d
				break
			}
		}
		var prg *PrgHit
		for _, id := range ids {
			if h, ok := prgs[id]; ok && (prg == nil || h.Eval < prg.Eval) {
				h := h
				prg = &h
			}
		}
		rows = append(rows, scoreGene(g, variants, resCols, susCols, opts.MinConcordance, opts.PromoterLength, desc, prg, prof))
		_ = bar.Add(1)
	}
	_ = bar.Finish()
	release()
	fmt.Fprintln(os.Stderr)

	SortRows(rows)

	counts := map[string]int{}
	for _, r := range rows {
		counts[r.Tier.Code]++
	}
	color.Cyan("✔ Scored %d genes, %d variants: A=%d B=%d C=%d D=%d",
		len(rows), len(variants), counts["A"], counts["B"], counts["C"], counts["D"])

	excelPath := filepath.Join(opts.OutDir, fmt.Sprintf("GeneSpace_%s_%d-%d.xlsx", opts.Chrom, opts.Start, opts.Stop))
	if err := WriteGeneSpaceExcel(rows, prof, excelPath); err != nil {
		return nil, "", fmt.Errorf("writing Excel file: %w", err)
	}
	return rows, excelPath, nil
}

// ─── Excel export ────────────────────────────────────────────────────────────

var excelHeaders = []string{
	"Gene", "Start", "Stop", "Strand",
	"# Variants", "Variants/kb", "# Concordant", "Best concordance", "Best mode", "Top concordant effect", "Variant score",
	"PRG %Identity", "PRG coverage %", "PRG E-value",
	"Gene description", "Matched terms", "TF", "Function score",
	"Tier", "Tier label",
}

func (r GeneSpaceRow) cells() []any {
	var pi, cov, ev any = "", "", ""
	if r.PRG != nil {
		pi = r.PRG.PercIdent
		cov = math.Round(r.PRG.Coverage*10) / 10
		ev = r.PRG.Eval
	}
	tf := "N"
	if r.IsTF {
		tf = "Y"
	}
	mode, eff := r.BestMode, r.TopEffect
	if r.NumConcordant == 0 {
		mode, eff = "-", ""
	}
	return []any{
		r.Gene, r.Start, r.Stop, r.Strand,
		r.NumVariants, r.VariantsPerKb, r.NumConcordant, math.Round(r.BestConcordance*1000) / 1000, mode, eff, r.VariantScore,
		pi, cov, ev,
		r.Description, strings.Join(r.MatchedTerms, "; "), tf, r.FunctionScore,
		r.Tier.Code, r.Tier.Label,
	}
}

// WriteGeneSpaceExcel writes the rows (already sorted) to a workbook with a
// GeneSpace data sheet and a Legend sheet.
func WriteGeneSpaceExcel(rows []GeneSpaceRow, prof Profile, outPath string) error {
	f := excelize.NewFile()
	defer f.Close()

	newStyle := func(s *excelize.Style) (int, error) { return f.NewStyle(s) }
	border := []excelize.Border{
		{Type: "left", Color: "CCCCCC", Style: 1},
		{Type: "right", Color: "CCCCCC", Style: 1},
		{Type: "top", Color: "CCCCCC", Style: 1},
		{Type: "bottom", Color: "CCCCCC", Style: 1},
	}
	headerStyle, err := newStyle(&excelize.Style{
		Font:      &excelize.Font{Bold: true, Color: "FFFFFF", Family: "Arial", Size: 10},
		Fill:      excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{"1F3864"}},
		Alignment: &excelize.Alignment{Horizontal: "center", Vertical: "center", WrapText: true},
		Border:    border,
	})
	if err != nil {
		return err
	}
	tierStyle := make(map[string]int)
	for _, t := range []Tier{TierA, TierB, TierC, TierD} {
		id, err := newStyle(&excelize.Style{
			Font:      &excelize.Font{Family: "Arial", Size: 10},
			Fill:      excelize.Fill{Type: "pattern", Pattern: 1, Color: []string{t.Hex}},
			Alignment: &excelize.Alignment{Vertical: "center"},
			Border:    border,
		})
		if err != nil {
			return err
		}
		tierStyle[t.Code] = id
	}

	// ── GeneSpace data sheet ──
	const data = "GeneSpace"
	if err := f.SetSheetName("Sheet1", data); err != nil {
		return err
	}
	widths := make([]int, len(excelHeaders))
	if err := f.SetSheetRow(data, "A1", &excelHeaders); err != nil {
		return err
	}
	lastHeader, _ := excelize.CoordinatesToCellName(len(excelHeaders), 1)
	_ = f.SetCellStyle(data, "A1", lastHeader, headerStyle)
	for i, h := range excelHeaders {
		widths[i] = len(h)
	}
	for ri, r := range rows {
		vals := r.cells()
		first, _ := excelize.CoordinatesToCellName(1, ri+2)
		last, _ := excelize.CoordinatesToCellName(len(vals), ri+2)
		if err := f.SetSheetRow(data, first, &vals); err != nil {
			return err
		}
		_ = f.SetCellStyle(data, first, last, tierStyle[r.Tier.Code])
		for ci, v := range vals {
			widths[ci] = max(widths[ci], len(fmt.Sprint(v)))
		}
	}
	for ci, w := range widths {
		name, _ := excelize.ColumnNumberToName(ci + 1)
		_ = f.SetColWidth(data, name, name, float64(min(w+3, 60)))
	}
	_ = f.SetPanes(data, &excelize.Panes{Freeze: true, YSplit: 1, TopLeftCell: "A2", ActivePane: "bottomLeft"})
	lastCell, _ := excelize.CoordinatesToCellName(len(excelHeaders), len(rows)+1)
	_ = f.AutoFilter(data, "A1:"+lastCell, nil)

	// ── Legend sheet ──
	const legend = "Legend"
	if _, err := f.NewSheet(legend); err != nil {
		return err
	}
	legendRows := [][]any{
		{"Tier", "Label", "Rule"},
		{TierA.Code, TierA.Label, "Variant score >= 3 and Function score >= 2"},
		{TierB.Code, TierB.Label, "Variant score >= 2, or a concordant variant and Function score >= 2"},
		{TierC.Code, TierC.Label, "Function score >= 2, no qualifying concordant variant"},
		{TierD.Code, TierD.Label, "In the window, nothing else"},
		{},
		{"Score", "Value", "Meaning"},
		{"Variant score", 4, "Concordant HIGH-impact variant (stop gained, frameshift, splice site, start lost)"},
		{"Variant score", 3, "Concordant MODERATE-impact variant (missense, in-frame indel)"},
		{"Variant score", 2, "Concordant UTR, splice-region or other LOW-impact (non-synonymous) variant"},
		{"Variant score", 1, "Concordant variant in the promoter window upstream of the gene"},
		{"Variant score", 0, "No concordant variant, or only intronic/synonymous/downstream ones"},
		{"Function score", 3, fmt.Sprintf("PRG pass (identity >= %g%%, coverage >= %g%%, E <= %g) and a strong term", PrgMinIdentity, PrgMinCoverage, PrgMaxEvalue)},
		{"Function score", 2, "PRG pass or a strong term"},
		{"Function score", 1, "Weak term, weaker PRG hit, or transcription factor"},
		{"Function score", 0, "Nothing"},
		{},
		{"Concordance", "", "Mean of the fraction of resistant and of susceptible lines that fit the pattern; missing calls ignored"},
		{"Best mode", "DOM", "Resistant lines carry the allele, susceptible lines do not"},
		{"Best mode", "REC", "Resistant lines are homozygous for the allele, susceptible lines are not"},
		{"Best mode", "DOM/REC", "Both patterns fit equally (e.g. R all 1/1, S all 0/0)"},
		{},
		{"Trait profile", prof.Name, prof.Description},
	}
	for ri, lr := range legendRows {
		cell, _ := excelize.CoordinatesToCellName(1, ri+1)
		if err := f.SetSheetRow(legend, cell, &lr); err != nil {
			return err
		}
	}
	_ = f.SetCellStyle(legend, "A1", "C1", headerStyle)
	_ = f.SetCellStyle(legend, "A7", "C7", headerStyle)
	for i, t := range []Tier{TierA, TierB, TierC, TierD} {
		c := fmt.Sprintf("A%d", i+2)
		e := fmt.Sprintf("C%d", i+2)
		_ = f.SetCellStyle(legend, c, e, tierStyle[t.Code])
	}
	_ = f.SetColWidth(legend, "A", "B", 18)
	_ = f.SetColWidth(legend, "C", "C", 100)

	idx, _ := f.GetSheetIndex(data)
	f.SetActiveSheet(idx)

	if err := f.SaveAs(outPath); err != nil {
		return fmt.Errorf("saving Excel file: %w", err)
	}
	color.Green("✔ Excel report written to %s", outPath)
	return nil
}
