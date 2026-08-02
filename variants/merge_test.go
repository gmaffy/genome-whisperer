package variants

import (
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// writeVCF writes a minimal but valid VCF. Uncompressed, which openVCF accepts,
// so no bgzip is needed.
func writeVCF(t *testing.T, path string, samples []string, records ...string) string {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		t.Fatal(err)
	}

	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##contig=<ID=A01,length=100000>\n")
	for _, k := range []string{"QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"} {
		b.WriteString("##INFO=<ID=" + k + ",Number=1,Type=Float,Description=\"" + k + "\">\n")
	}
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	b.WriteString("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(samples) > 0 {
		b.WriteString("\tFORMAT\t" + strings.Join(samples, "\t"))
	}
	b.WriteString("\n")
	for _, r := range records {
		b.WriteString(r + "\n")
	}

	if err := os.WriteFile(path, []byte(b.String()), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

// ---------------------------------------------------------------------------
// Sample names
// ---------------------------------------------------------------------------

func TestVcfSampleNames(t *testing.T) {
	p := writeVCF(t, filepath.Join(t.TempDir(), "x.vcf"), []string{"S1", "S2", "S3"})
	got, err := vcfSampleNames(p)
	if err != nil {
		t.Fatalf("vcfSampleNames: %v", err)
	}
	if strings.Join(got, ",") != "S1,S2,S3" {
		t.Errorf("got %v", got)
	}
}

// A sites-only VCF has no sample columns; treating that as "zero samples" would
// let writeSampleMap emit a blank entry.
func TestVcfSampleNamesSitesOnly(t *testing.T) {
	p := writeVCF(t, filepath.Join(t.TempDir(), "sites.vcf"), nil)
	if _, err := vcfSampleNames(p); err == nil {
		t.Error("expected an error for a VCF with no sample columns")
	}
}

// GATK and GLnexus emit sorted sample columns while the gVCFs arrive in discovery
// order. An order-sensitive comparison reported a mismatch on nearly every run and
// forced a pointless re-merge.
func TestSampleNamesMatchIgnoresOrder(t *testing.T) {
	cases := []struct {
		a, b []string
		want bool
	}{
		{[]string{"S1", "S2"}, []string{"S2", "S1"}, true},
		{[]string{"S1", "S2"}, []string{"S1", "S2"}, true},
		{[]string{"S1", "S2"}, []string{"S1", "S2", "S3"}, false},
		{[]string{"S1", "S2"}, []string{"S1", "S3"}, false},
		{nil, nil, true},
		{[]string{"S1"}, nil, false},
	}
	for _, c := range cases {
		if got := sampleNamesMatch(c.a, c.b); got != c.want {
			t.Errorf("sampleNamesMatch(%v, %v) = %v, want %v", c.a, c.b, got, c.want)
		}
	}
}

// The comparison must not reorder its inputs: the caller still uses the gVCF list
// afterwards to build the sample map.
func TestSampleNamesMatchDoesNotMutateInputs(t *testing.T) {
	a := []string{"S2", "S1"}
	b := []string{"S1", "S2"}
	sampleNamesMatch(a, b)
	if a[0] != "S2" || b[0] != "S1" {
		t.Errorf("inputs were reordered: a=%v b=%v", a, b)
	}
}

// ---------------------------------------------------------------------------
// Sample map
// ---------------------------------------------------------------------------

func TestWriteSampleMapUsesHeaderNames(t *testing.T) {
	dir := t.TempDir()
	// Filename and header sample name deliberately disagree: GenomicsDBImport keys
	// on the header, so that is what has to be written.
	g1 := writeVCF(t, filepath.Join(dir, "wrong_name_1.g.vcf"), []string{"SAMPLE_A"})
	g2 := writeVCF(t, filepath.Join(dir, "wrong_name_2.g.vcf"), []string{"SAMPLE_B"})

	out := filepath.Join(dir, "sample_map.txt")
	if err := writeSampleMap(out, []string{g1, g2}); err != nil {
		t.Fatalf("writeSampleMap: %v", err)
	}

	data, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d: %q", len(lines), data)
	}
	for i, want := range []struct{ name, path string }{{"SAMPLE_A", g1}, {"SAMPLE_B", g2}} {
		fields := strings.Split(lines[i], "\t")
		if len(fields) != 2 || fields[0] != want.name || fields[1] != want.path {
			t.Errorf("line %d = %q, want %q\\t%q", i, lines[i], want.name, want.path)
		}
	}
}

// Duplicate sample names make GenomicsDBImport fail minutes in, with an error that
// does not name the files. Catch it before the import starts.
func TestWriteSampleMapRejectsDuplicateSamples(t *testing.T) {
	dir := t.TempDir()
	g1 := writeVCF(t, filepath.Join(dir, "a.g.vcf"), []string{"SAME"})
	g2 := writeVCF(t, filepath.Join(dir, "b.g.vcf"), []string{"SAME"})

	err := writeSampleMap(filepath.Join(dir, "map.txt"), []string{g1, g2})
	if err == nil {
		t.Fatal("expected an error for duplicate sample names")
	}
	if !strings.Contains(err.Error(), "duplicate sample") {
		t.Errorf("unhelpful error: %v", err)
	}
}

// The old CreateSampleMap silently skipped anything whose name did not end in
// .g.vcf.gz, which could yield an empty sample map and a baffling GATK error.
func TestWriteSampleMapRejectsMultiSampleInput(t *testing.T) {
	dir := t.TempDir()
	multi := writeVCF(t, filepath.Join(dir, "joint.vcf"), []string{"S1", "S2"})

	err := writeSampleMap(filepath.Join(dir, "map.txt"), []string{multi})
	if err == nil {
		t.Fatal("expected an error: a per-sample gVCF must hold exactly one sample")
	}
	if !strings.Contains(err.Error(), "expected 1 sample") {
		t.Errorf("unhelpful error: %v", err)
	}
}

// ---------------------------------------------------------------------------
// GATK intervals — the contigs bug
// ---------------------------------------------------------------------------

func TestGatkIntervalArgsSingleSequence(t *testing.T) {
	region, extra, err := gatkIntervalArgs(t.TempDir(), []SeqInfo{{ID: "A01", Len: 100}})
	if err != nil {
		t.Fatalf("gatkIntervalArgs: %v", err)
	}
	if region != "A01" {
		t.Errorf("region = %q, want A01", region)
	}
	if extra != "" {
		t.Errorf("a single chromosome needs no extra arguments, got %q", extra)
	}
}

// The contigs group covers many sequences and "contigs" is not a sequence name.
// Passing the label made GenomicsDBImport fail with "contig contigs not present in
// the sequence dictionary", so contigs never merged under --merger gatk.
func TestGatkIntervalArgsWritesIntervalListForContigs(t *testing.T) {
	dir := t.TempDir()
	seqs := []SeqInfo{{ID: "ctg1", Len: 10}, {ID: "ctg2", Len: 20}, {ID: "ctg3", Len: 30}}

	region, extra, err := gatkIntervalArgs(dir, seqs)
	if err != nil {
		t.Fatalf("gatkIntervalArgs: %v", err)
	}
	if region == "contigs" {
		t.Fatal("the label must never be passed to -L; it is not a sequence name")
	}
	if filepath.Dir(region) != dir {
		t.Errorf("interval list should be written into %s, got %s", dir, region)
	}
	if !strings.Contains(extra, "--merge-input-intervals") {
		t.Errorf("many intervals should be merged, got extra = %q", extra)
	}

	data, err := os.ReadFile(region)
	if err != nil {
		t.Fatalf("interval list not written: %v", err)
	}
	got := strings.Split(strings.TrimSpace(string(data)), "\n")
	if len(got) != len(seqs) {
		t.Fatalf("expected %d intervals, got %d: %q", len(seqs), len(got), data)
	}
	for i, s := range seqs {
		if got[i] != s.ID {
			t.Errorf("interval %d = %q, want %q", i, got[i], s.ID)
		}
	}
}

func TestGatkIntervalArgsRejectsEmptySequenceList(t *testing.T) {
	if _, _, err := gatkIntervalArgs(t.TempDir(), nil); err == nil {
		t.Error("expected an error when no sequences are given")
	}
}

// ---------------------------------------------------------------------------
// Concatenation
// ---------------------------------------------------------------------------

// With one input ConcatenateVcfs copies instead of shelling out, so this path is
// testable without GATK. It must bring the index along.
func TestConcatenateVcfsSingleInputCopiesWithIndex(t *testing.T) {
	dir := t.TempDir()
	src := writeVCF(t, filepath.Join(dir, "A01.joint.vcf"), []string{"S1"})
	if err := os.WriteFile(src+".tbi", []byte("fake-index"), 0644); err != nil {
		t.Fatal(err)
	}
	dst := filepath.Join(dir, "COTTON.ad1.1.gatk_gatk.all.vcf.gz")

	got, err := ConcatenateVcfs([]string{src}, dst, false)
	if err != nil {
		t.Fatalf("ConcatenateVcfs: %v", err)
	}
	if got != dst {
		t.Errorf("returned %q, want %q", got, dst)
	}
	if _, err := os.Stat(dst); err != nil {
		t.Errorf("output not written: %v", err)
	}
	if _, err := os.Stat(dst + ".tbi"); err != nil {
		t.Errorf("index not copied alongside: %v", err)
	}
}

func TestConcatenateVcfsRejectsEmptyInput(t *testing.T) {
	if _, err := ConcatenateVcfs(nil, filepath.Join(t.TempDir(), "out.vcf.gz"), false); err == nil {
		t.Error("expected an error when there is nothing to concatenate")
	}
}

// ---------------------------------------------------------------------------
// Dictionary ordering of the merge list
// ---------------------------------------------------------------------------

// gatk MergeVcfs needs its inputs in dictionary order. The previous code appended
// them in whatever order workers finished, so the concatenated VCF could be
// mis-sorted. This reproduces the sort MergeGvcfs applies.
func TestLabelsSortIntoDictionaryOrderWithContigsLast(t *testing.T) {
	dict := writeDict(t, []SeqInfo{
		{ID: "A01", Len: 500},
		{ID: "A02", Len: 900}, // longer, but second in the dictionary
		{ID: "A03", Len: 100},
	})
	order, err := dictOrder(dict)
	if err != nil {
		t.Fatal(err)
	}

	// Deliberately shuffled, as workers would produce them.
	labels := []string{"contigs", "A03", "A01", "A02"}
	sort.Slice(labels, func(i, j int) bool {
		pi, oki := order[labels[i]]
		pj, okj := order[labels[j]]
		if !oki {
			pi = len(order)
		}
		if !okj {
			pj = len(order)
		}
		return pi < pj
	})

	want := "A01,A02,A03,contigs"
	if got := strings.Join(labels, ","); got != want {
		t.Errorf("sorted order = %s, want %s", got, want)
	}
}
