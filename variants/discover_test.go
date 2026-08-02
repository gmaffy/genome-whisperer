package variants

import (
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// touch creates an empty file and any missing parent directories.
func touch(t *testing.T, path string) string {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, nil, 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

// alignmentPath builds a path in the standard data directory layout:
//
//	<root>/<species>/<project>/<sample>/reference_genomes/<refVer>/bams/<file>
//
// The species component is lower-cased, matching FindBams.
func alignmentPath(root, species, project, sample, refVer, file string) string {
	return filepath.Join(root, strings.ToLower(species), project, sample,
		"reference_genomes", refVer, "bams", file)
}

// writeRef creates a reference fasta with the .dict beside it, which is where
// every stage looks for it.
func writeRef(t *testing.T, dir string, seqs []SeqInfo) string {
	t.Helper()
	fasta := filepath.Join(dir, "ref.fa")
	if err := os.MkdirAll(dir, 0755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(fasta, []byte(">A01\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}

	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\n")
	for _, s := range seqs {
		b.WriteString("@SQ\tSN:" + s.ID + "\tLN:" + itoa(s.Len) + "\n")
	}
	if err := os.WriteFile(filepath.Join(dir, "ref.dict"), []byte(b.String()), 0644); err != nil {
		t.Fatal(err)
	}
	return fasta
}

func sampleNames(samples []SampleWork) []string {
	var out []string
	for _, s := range samples {
		out = append(out, s.Sample)
	}
	sort.Strings(out)
	return out
}

// ---------------------------------------------------------------------------
// Inline / config mode
// ---------------------------------------------------------------------------

func TestFindGvcfSamplesInline(t *testing.T) {
	dir := t.TempDir()
	b1 := touch(t, filepath.Join(dir, "S1.bqsr.cram"))
	b2 := touch(t, filepath.Join(dir, "S2.bqsr.bam"))
	missing := filepath.Join(dir, "S3.bqsr.cram")

	samples, skipped, err := FindGvcfSamples(Options{Bams: []string{b1, b2, missing}})
	if err != nil {
		t.Fatalf("FindGvcfSamples: %v", err)
	}

	if got := strings.Join(sampleNames(samples), ","); got != "S1.bqsr,S2.bqsr" {
		t.Errorf("samples = %s", got)
	}
	if len(skipped) != 1 || skipped[0] != missing {
		t.Errorf("skipped = %v, want just the missing file", skipped)
	}
	// Inline mode has no per-sample home directory; GvcfPath falls back to OutDir.
	for _, s := range samples {
		if s.CramDir != "" {
			t.Errorf("[%s] CramDir should be empty in inline mode, got %q", s.Sample, s.CramDir)
		}
	}
}

// A missing bam must not abort the run: the rest of the cohort still gets called.
func TestFindGvcfSamplesInlineAllMissing(t *testing.T) {
	dir := t.TempDir()
	samples, skipped, err := FindGvcfSamples(Options{Bams: []string{
		filepath.Join(dir, "a.cram"), filepath.Join(dir, "b.cram"),
	}})
	if err != nil {
		t.Fatalf("FindGvcfSamples should not error on missing inputs: %v", err)
	}
	if len(samples) != 0 || len(skipped) != 2 {
		t.Errorf("got %d samples and %d skipped, want 0 and 2", len(samples), len(skipped))
	}
}

// ---------------------------------------------------------------------------
// Data directory mode
// ---------------------------------------------------------------------------

func TestFindGvcfSamplesDataDir(t *testing.T) {
	root := t.TempDir()
	const species, refVer = "cotton", "AD1.1"

	// Short-read sample: the bqsr alignment is the one to use.
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.bqsr.cram"))
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.rgmd.cram"))
	// Long-read sample ("...lr"): the rgmd alignment is the one to use, because
	// BQSR is not applied to long reads.
	touch(t, alignmentPath(root, species, "proj2", "S2lr", refVer, "S2lr.rgmd.cram"))
	// No alignment at all: must be skipped, not fatal.
	if err := os.MkdirAll(filepath.Join(root, species, "proj1", "S3", "reference_genomes"), 0755); err != nil {
		t.Fatal(err)
	}

	samples, skipped, err := FindGvcfSamples(Options{
		DataDir: root, Species: species, RefVer: refVer, Caller: "gatk",
	})
	if err != nil {
		t.Fatalf("FindGvcfSamples: %v", err)
	}

	if got := strings.Join(sampleNames(samples), ","); got != "S1,S2lr" {
		t.Errorf("samples = %s, want S1,S2lr", got)
	}
	if len(skipped) != 1 || skipped[0] != "S3" {
		t.Errorf("skipped = %v, want [S3]", skipped)
	}

	for _, s := range samples {
		// CramDir must be the reference_genomes directory: GvcfPath appends
		// <refVer>/<gatk_gvcfs|dv_gvcfs> to it.
		if filepath.Base(s.CramDir) != "reference_genomes" {
			t.Errorf("[%s] CramDir = %s, want .../reference_genomes", s.Sample, s.CramDir)
		}
		if !filepath.IsAbs(s.Cram) {
			t.Errorf("[%s] alignment path should be absolute, got %s", s.Sample, s.Cram)
		}
		switch s.Sample {
		case "S1":
			if !strings.Contains(s.Cram, "bqsr") {
				t.Errorf("[S1] should use the bqsr alignment, got %s", s.Cram)
			}
		case "S2lr":
			if !strings.Contains(s.Cram, "rgmd") {
				t.Errorf("[S2lr] should use the rgmd alignment, got %s", s.Cram)
			}
		}
	}
}

// --bqsr=false (NoBqsr) selects the rgmd alignment for short-read samples too.
func TestFindGvcfSamplesDataDirNoBqsr(t *testing.T) {
	root := t.TempDir()
	const species, refVer = "cotton", "AD1.1"
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.bqsr.cram"))
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.rgmd.cram"))

	samples, _, err := FindGvcfSamples(Options{
		DataDir: root, Species: species, RefVer: refVer, Caller: "gatk", NoBqsr: true,
	})
	if err != nil {
		t.Fatalf("FindGvcfSamples: %v", err)
	}
	if len(samples) != 1 {
		t.Fatalf("expected 1 sample, got %d", len(samples))
	}
	if !strings.Contains(samples[0].Cram, "rgmd") {
		t.Errorf("NoBqsr should select the rgmd alignment, got %s", samples[0].Cram)
	}
}

// Two candidate alignments are ambiguous. Picking one silently risks calling from
// a stale file, so the sample is reported and skipped.
func TestFindGvcfSamplesDataDirAmbiguousAlignment(t *testing.T) {
	root := t.TempDir()
	const species, refVer = "cotton", "AD1.1"
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.bqsr.cram"))
	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.old.bqsr.cram"))

	samples, skipped, err := FindGvcfSamples(Options{
		DataDir: root, Species: species, RefVer: refVer, Caller: "gatk",
	})
	if err != nil {
		t.Fatalf("FindGvcfSamples: %v", err)
	}
	if len(samples) != 0 {
		t.Errorf("an ambiguous sample must not be used, got %+v", samples)
	}
	if len(skipped) != 1 || skipped[0] != "S1" {
		t.Errorf("skipped = %v, want [S1]", skipped)
	}
}

// Species case must not change which directory is searched. Sample discovery used
// the raw --species while FindBams lower-cased it, so on a case-sensitive
// filesystem the two looked in different places.
func TestFindGvcfSamplesSpeciesCaseInsensitive(t *testing.T) {
	root := t.TempDir()
	const refVer = "AD1.1"
	touch(t, alignmentPath(root, "cotton", "proj1", "S1", refVer, "S1.bqsr.cram"))

	for _, species := range []string{"cotton", "Cotton", "COTTON"} {
		samples, _, err := FindGvcfSamples(Options{
			DataDir: root, Species: species, RefVer: refVer, Caller: "gatk",
		})
		if err != nil {
			t.Fatalf("--species %s: %v", species, err)
		}
		if len(samples) != 1 {
			t.Errorf("--species %s found %d samples, want 1", species, len(samples))
		}
	}
}

// ---------------------------------------------------------------------------
// Existing gVCF inventory
// ---------------------------------------------------------------------------

// FindExistingGvcfs must group by chromosome and look exactly where CreateGvcfs
// writes, since it builds paths with the same GvcfPath function. The old code
// globbed *.g.vcf.gz per sample, which mixed every chromosome into one flat list.
func TestFindExistingGvcfsGroupsByChromosome(t *testing.T) {
	root := t.TempDir()
	const species, refVer = "cotton", "AD1.1"

	opts := Options{
		DataDir: root, Species: species, RefVer: refVer, Caller: "gatk", Merger: "gatk",
		RefFasta: writeRef(t, filepath.Join(root, "genome"), []SeqInfo{
			{ID: "A01", Len: 1000}, {ID: "A02", Len: 900},
		}),
		// Avoids shelling out to bcftools; presence is what this test is about.
		SkipVerification: true,
	}

	touch(t, alignmentPath(root, species, "proj1", "S1", refVer, "S1.bqsr.cram"))
	touch(t, alignmentPath(root, species, "proj1", "S2", refVer, "S2.bqsr.cram"))

	samples, _, err := FindGvcfSamples(opts)
	if err != nil {
		t.Fatal(err)
	}
	if len(samples) != 2 {
		t.Fatalf("expected 2 samples, got %d", len(samples))
	}

	// Both samples have A01; only S1 has A02.
	for _, s := range samples {
		touch(t, GvcfPath(opts, s, "A01"))
	}
	for _, s := range samples {
		if s.Sample == "S1" {
			touch(t, GvcfPath(opts, s, "A02"))
		}
	}

	gvcfs, err := FindExistingGvcfs(opts)
	if err != nil {
		t.Fatalf("FindExistingGvcfs: %v", err)
	}

	if len(gvcfs["A01"]) != 2 {
		t.Errorf("A01 should have 2 gVCFs, got %d", len(gvcfs["A01"]))
	}
	if len(gvcfs["A02"]) != 1 {
		t.Errorf("A02 should have 1 gVCF, got %d", len(gvcfs["A02"]))
	}
	if _, ok := gvcfs["contigs"]; ok {
		t.Error("a 2-sequence reference has no contigs group")
	}

	// Sorted so the sample map and any sample-name comparison are deterministic.
	if !sort.StringsAreSorted(gvcfs["A01"]) {
		t.Errorf("A01 gVCFs are not sorted: %v", gvcfs["A01"])
	}

	// This is the completeness rule MergeGvcfs applies: the cohort size is the
	// largest sample set seen, and A02 is short of it.
	expected := 0
	for _, paths := range gvcfs {
		if len(paths) > expected {
			expected = len(paths)
		}
	}
	if expected != 2 {
		t.Fatalf("cohort size = %d, want 2", expected)
	}
	if len(gvcfs["A02"]) >= expected {
		t.Error("A02 should be detected as incomplete")
	}
}
