package genespace

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/xuri/excelize/v2"
)

func TestNormGT(t *testing.T) {
	cases := []struct {
		in      string
		a, b    string
		missing bool
	}{
		{"0/1", "0", "1", false},
		{"1/0", "0", "1", false},
		{"1|0", "0", "1", false},
		{"1", "1", "1", false},
		{".", "", "", true},
		{"./.", "", "", true},
		{"./1", "", "", true},
		{"", "", "", true},
		{"0/1/1", "", "", true},
	}
	for _, c := range cases {
		a, b, ok := normGT(c.in)
		if ok == c.missing || a != c.a || b != c.b {
			t.Errorf("normGT(%q) = %q,%q,%v", c.in, a, b, ok)
		}
	}
}

func TestTestSNP(t *testing.T) {
	cases := []struct {
		name     string
		res, sus []string
		conc     float64
		mode     string
		allele   string
		testable bool
	}{
		{"both modes fit", []string{"1/1", "1/1"}, []string{"0/0", "0/0"}, 1, "DOM/REC", "1", true},
		{"dominant with het R", []string{"0/1", "1/1"}, []string{"0/0", "0/0"}, 1, "DOM", "1", true},
		{"recessive", []string{"1/1", "1/1"}, []string{"0/1", "0/0"}, 1, "REC", "1", true},
		{"ref allele resistant", []string{"0/0", "0/0"}, []string{"1/1", "1/1"}, 1, "DOM/REC", "0", true},
		{"unordered het", []string{"1/0", "0|1"}, []string{"0/0"}, 1, "DOM", "1", true},
		{"missing ignored", []string{"1/1", "./."}, []string{"0/0", "0/0"}, 1, "DOM/REC", "1", true},
		{"monomorphic", []string{"0/0", "0/0", "0/0"}, []string{"0/0"}, 0.5, "DOM/REC", "0", true},
		{"no sus calls", []string{"1/1"}, []string{"./."}, 0, "", "", false},
		{"low call rate", []string{"1/1", ".", "."}, []string{"0/0", "."}, 0, "", "", false},
	}
	for _, c := range cases {
		got := testSNP(c.res, c.sus)
		if got.Testable != c.testable || abs(got.Concordance-c.conc) > 1e-9 || got.Mode != c.mode || got.Allele != c.allele {
			t.Errorf("%s: got %+v", c.name, got)
		}
	}

	// One mismatching S line out of two must not be hidden by a large R group.
	got := testSNP([]string{"1/1", "1/1", "1/1", "1/1", "1/1", "1/1", "1/1", "1/1"}, []string{"0/0", "1/1"})
	if got.Concordance >= 0.9 {
		t.Errorf("unbalanced groups: concordance %v should be < 0.9", got.Concordance)
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func TestProfileTerms(t *testing.T) {
	p, err := GetProfile(DefaultProfile)
	if err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		desc             string
		strong, weak, tf bool
	}{
		{"MFS superfamily transporter", false, false, false},
		{"organic cation transporter", false, false, false},
		{"spinach-type protein", false, false, false},
		{"40S ribosomal protein S6 RPS6", false, false, false},
		{"WRKY transcription factor 33", false, false, true},
		{"ethylene-responsive transcription factor ERF118", false, false, true},
		{"MLO-like protein 1", true, false, false},
		{"eukaryotic translation initiation factor 4E-1", true, false, false},
		{"disease resistance protein RPP13", true, false, false},
		{"putative NBS-LRR protein", true, false, false},
		{"LRR receptor-like serine/threonine-protein kinase FLS2", true, true, false},
		{"basic form of pathogenesis-related protein 1", false, true, false},
	}
	for _, c := range cases {
		s, w, f := matchTerms(p.Strong, c.desc), matchTerms(p.Weak, c.desc), matchTerms(p.TF, c.desc)
		if (len(s) > 0) != c.strong || (len(w) > 0) != c.weak || (len(f) > 0) != c.tf {
			t.Errorf("%q: strong=%v weak=%v tf=%v", c.desc, s, w, f)
		}
	}
	if _, err := GetProfile("nope"); err == nil || !strings.Contains(err.Error(), DefaultProfile) {
		t.Errorf("unknown profile error should list available profiles, got %v", err)
	}
}

func TestImpactScore(t *testing.T) {
	cases := []struct {
		impact, effect string
		upstream       bool
		want           int
	}{
		{"HIGH", "STOP_GAINED", false, 4},
		{"MODERATE", "NON_SYNONYMOUS_CODING", false, 3},
		{"LOW", "SYNONYMOUS_CODING", false, 0},
		{"LOW", "SPLICE_SITE_REGION", false, 2},
		{"MODIFIER", "UTR_5_PRIME", false, 2},
		{"MODIFIER", "INTRON", false, 0},
		{"", "FRAME_SHIFT", false, 4},
		{"HIGH", "STOP_GAINED", true, 1},
	}
	for _, c := range cases {
		if got := impactScore(c.impact, c.effect, c.upstream); got != c.want {
			t.Errorf("impactScore(%q,%q,%v) = %d, want %d", c.impact, c.effect, c.upstream, got, c.want)
		}
	}
}

func TestClassifyTier(t *testing.T) {
	cases := []struct {
		v, f, n int
		want    string
	}{
		{4, 2, 1, "A"},
		{3, 3, 1, "A"},
		{4, 1, 1, "B"},
		{2, 0, 1, "B"},
		{0, 2, 1, "B"}, // intronic concordant + strong function
		{1, 2, 1, "B"},
		{0, 3, 0, "C"},
		{1, 1, 1, "D"},
		{0, 0, 0, "D"},
	}
	for _, c := range cases {
		if got := ClassifyTier(c.v, c.f, c.n).Code; got != c.want {
			t.Errorf("ClassifyTier(%d,%d,%d) = %s, want %s", c.v, c.f, c.n, got, c.want)
		}
	}
}

func TestFunctionScore(t *testing.T) {
	pass := &PrgHit{PercIdent: 60, Coverage: 90, Eval: 1e-50}
	weak := &PrgHit{PercIdent: 30, Coverage: 20, Eval: 1e-5}
	x := []string{"x"}
	cases := []struct {
		prg            *PrgHit
		strong, wk, tf []string
		want           int
	}{
		{pass, x, nil, nil, 3},
		{pass, nil, nil, nil, 2},
		{nil, x, nil, nil, 2},
		{weak, nil, nil, nil, 1},
		{nil, nil, x, nil, 1},
		{nil, nil, nil, x, 1},
		{nil, nil, nil, nil, 0},
	}
	for i, c := range cases {
		if got := functionScore(c.prg, c.strong, c.wk, c.tf); got != c.want {
			t.Errorf("case %d: got %d want %d", i, got, c.want)
		}
	}
}

func write(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestLoadGenes(t *testing.T) {
	dir := t.TempDir()
	gff := write(t, dir, "a.gff", strings.Join([]string{
		"##gff-version 3",
		"1\t.\tgene\t100\t400\t.\t+\t.\tID=g1",
		"1\t.\tmRNA\t100\t300\t.\t+\t.\tID=g1.1;Parent=g1",
		"1\t.\tmRNA\t150\t400\t.\t+\t.\tID=g1.2;Parent=g1",
		"1\t.\tmRNA\t900\t1000\t.\t-\t.\tbroken-attribute",
		"1\t.\tmRNA\t5000\t6000\t.\t-\t.\tID=g3.1;Parent=g3",
		"2\t.\tmRNA\t100\t200\t.\t+\t.\tID=g4.1;Parent=g4",
	}, "\n"))

	// Window 350-5500 overlaps g1 (starts before the window) and g3.
	genes, err := loadGenes(gff, "1", 350, 5500)
	if err != nil {
		t.Fatal(err)
	}
	if len(genes) != 2 || genes[0].ID != "g1" || genes[1].ID != "g3" {
		t.Fatalf("got %+v", genes)
	}
	g := genes[0]
	if g.Start != 100 || g.Stop != 400 || len(g.Transcripts) != 2 || g.Strand != "+" {
		t.Errorf("isoforms not collapsed: %+v", g)
	}
}

// TestGeneSpaceEndToEnd runs the whole pipeline on a small fixture and reads
// the workbook back.
func TestGeneSpaceEndToEnd(t *testing.T) {
	dir := t.TempDir()
	gff := write(t, dir, "a.gff", strings.Join([]string{
		"1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=nlr.1;Parent=nlr",       // NLR + missense → A
		"1\t.\tmRNA\t3000\t4000\t.\t-\t.\tID=mlo.1;Parent=mlo",       // S-gene, recessive stop → A
		"1\t.\tmRNA\t6000\t7000\t.\t+\t.\tID=prom.1;Parent=prom",     // promoter variant only → D
		"1\t.\tmRNA\t8000\t9000\t.\t+\t.\tID=tf.1;Parent=tf",         // TF, no variants → D
		"1\t.\tmRNA\t10000\t11000\t.\t+\t.\tID=super.1;Parent=super", // "superfamily", non-concordant → D
	}, "\n"))
	vt := write(t, dir, "v.tsv", strings.Join([]string{
		"CHROM\tPOS\tREF\tALT\tR1.GT\tR2.GT\tS1.GT\tS2.GT\tSNPEFF_EFFECT\tSNPEff_IMPACT",
		"1\t1500\tA\tG\t0/1\t1/1\t0/0\t0/0\tNON_SYNONYMOUS_CODING\tMODERATE",
		"1\t3500\tC\tT\t1/1\t1/1\t0/1\t0/0\tSTOP_GAINED\tHIGH",
		"1\t5500\tC\tT\t1/1\t1/1\t0/0\t0/0\tUPSTREAM\tMODIFIER",
		"1\t10500\tG\tA\t0/1\t0/0\t0/1\t0/0\tNON_SYNONYMOUS_CODING\tMODERATE",
		"2\t1500\tA\tG\t1/1\t1/1\t0/0\t0/0\tSTOP_GAINED\tHIGH",
	}, "\n"))
	desc := write(t, dir, "d.tsv", strings.Join([]string{
		"nlr\tdisease resistance protein RPP13",
		"mlo.1\tMLO-like protein 6",
		"tf\tWRKY transcription factor 33",
		"super\tMFS superfamily transporter",
	}, "\n"))
	prg := write(t, dir, "p.tsv", "QseqID\tQlen\tSseqID\tSlen\tPercIdent\tLength\tEval\nnlr.1\t900\t1\t900\t75\t880\t1e-100\n")

	rows, path, err := GeneSpace(Options{
		GFF: gff, VariantTable: vt, Chrom: "1", Start: 0, Stop: 20000,
		ResLines: []string{"R1", " R2"}, SusLines: []string{"S1.GT", "S2"},
		DescFile: desc, PrgFile: prg, OutDir: filepath.Join(dir, "out"),
		MinConcordance: 0.9, PromoterLength: 2000,
	})
	if err != nil {
		t.Fatal(err)
	}

	byGene := map[string]GeneSpaceRow{}
	for _, r := range rows {
		byGene[r.Gene] = r
	}
	check := func(gene, tier string, vs, fs int) {
		t.Helper()
		r := byGene[gene]
		if r.Tier.Code != tier || r.VariantScore != vs || r.FunctionScore != fs {
			t.Errorf("%s: tier=%s variant=%d function=%d, want %s/%d/%d (%+v)",
				gene, r.Tier.Code, r.VariantScore, r.FunctionScore, tier, vs, fs, r)
		}
	}
	check("nlr", "A", 3, 3)
	check("mlo", "A", 4, 2)
	check("prom", "D", 1, 0)
	check("tf", "D", 0, 1)
	check("super", "D", 0, 0)
	if byGene["mlo"].BestMode != "REC" {
		t.Errorf("mlo mode = %q, want REC", byGene["mlo"].BestMode)
	}
	if byGene["super"].NumVariants != 1 || byGene["super"].NumConcordant != 0 {
		t.Errorf("super: %+v", byGene["super"])
	}
	if rows[0].Gene != "mlo" || rows[1].Gene != "nlr" {
		t.Errorf("sort order: %s, %s", rows[0].Gene, rows[1].Gene)
	}

	f, err := excelize.OpenFile(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	xr, err := f.GetRows("GeneSpace")
	if err != nil {
		t.Fatal(err)
	}
	if len(xr) != 6 || xr[0][0] != "Gene" || xr[1][0] != "mlo" || xr[1][len(xr[1])-2] != "A" {
		t.Errorf("unexpected sheet contents: %v", xr[:2])
	}
	if filepath.Base(path) != "GeneSpace_1_0-20000.xlsx" {
		t.Errorf("output name %s", path)
	}
}

func TestGeneSpaceErrors(t *testing.T) {
	dir := t.TempDir()
	gff := write(t, dir, "a.gff", "1\t.\tmRNA\t1\t10\t.\t+\t.\tID=a.1;Parent=a\n")
	vt := write(t, dir, "v.tsv", "CHROM\tPOS\tR1.GT\tS1.GT\n1\t5\t1/1\t0/0\n")
	base := Options{GFF: gff, VariantTable: vt, Chrom: "1", Stop: 100,
		ResLines: []string{"R1"}, SusLines: []string{"S1"}, OutDir: dir, MinConcordance: 0.9}

	o := base
	o.SusLines = []string{"S9"}
	if _, _, err := GeneSpace(o); err == nil || !strings.Contains(err.Error(), "S9") {
		t.Errorf("unknown line: %v", err)
	}
	o = base
	o.Trait = "flowering"
	if _, _, err := GeneSpace(o); err == nil || !strings.Contains(err.Error(), "disease-resistance") {
		t.Errorf("unknown trait: %v", err)
	}
	// No PRG or description file is fine.
	if _, _, err := GeneSpace(base); err != nil {
		t.Errorf("optional inputs: %v", err)
	}
}
