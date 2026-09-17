package variants

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// ExpectedUnits is the denominator every completeness figure divides by, so
// these tests pin both the rule and — more importantly — that it still agrees
// with the functions that name the files.
//
// The stakes are on record: a hand-built inventory of this same estate used the
// raw @SQ count as the denominator and reported 81,202 missing gVCFs for a
// sample that was complete.

// dotted mirrors pepo v4.1, whose sequence IDs all contain a dot.
func dotted(n int) []SeqInfo {
	seqs := make([]SeqInfo, 0, n)
	for i := 0; i < n; i++ {
		// Descending lengths keep file order and length order aligned.
		seqs = append(seqs, SeqInfo{ID: "Cp4.1LG" + pad2(i), Len: 50_000_000 - i*1_000_000})
	}
	return seqs
}

func pad2(i int) string {
	s := itoa(i)
	if len(s) < 2 {
		return "0" + s
	}
	return s
}

func TestExpectedUnitsSanitisesDottedSeqids(t *testing.T) {
	// 21 sequences: every one gets its own job, so there is no batch.
	units, err := ExpectedUnits(writeDict(t, dotted(21)))
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	if len(units) != 21 {
		t.Fatalf("got %d units, want 21", len(units))
	}
	for _, u := range units {
		if u.Kind != UnitChrom {
			t.Errorf("unit %q kind = %q, want %q", u.ID, u.Kind, UnitChrom)
		}
		if strings.Contains(u.Label, ".") {
			t.Errorf("label %q still contains a dot", u.Label)
		}
		if !strings.Contains(u.ID, ".") {
			t.Errorf("ID %q lost its dot; the raw seqid must survive", u.ID)
		}
	}
	if units[0].ID != "Cp4.1LG00" || units[0].Label != "Cp4_1LG00" {
		t.Errorf("first unit = {%q, %q}, want {Cp4.1LG00, Cp4_1LG00}", units[0].ID, units[0].Label)
	}
}

func TestExpectedUnitsAddsBatchOnlyAboveTwentyOne(t *testing.T) {
	// The boundary is the rule: at 21 every sequence is called individually.
	at, err := ExpectedUnits(writeDict(t, dotted(21)))
	if err != nil {
		t.Fatalf("ExpectedUnits(21): %v", err)
	}
	for _, u := range at {
		if u.Kind == UnitContigs {
			t.Fatalf("21 sequences produced a %q unit; there should be none", BatchedLabel)
		}
	}

	above, err := ExpectedUnits(writeDict(t, dotted(22)))
	if err != nil {
		t.Fatalf("ExpectedUnits(22): %v", err)
	}
	if len(above) != 22 {
		t.Fatalf("got %d units for 22 sequences, want 22 (21 chroms + 1 batch)", len(above))
	}
	last := above[len(above)-1]
	if last.Kind != UnitContigs || last.ID != BatchedLabel || last.Label != BatchedLabel {
		t.Fatalf("last unit = {%q, %q, %q}, want the %q batch", last.ID, last.Label, last.Kind, BatchedLabel)
	}
	if len(last.Members) != 1 {
		t.Errorf("batch holds %d members, want 1", len(last.Members))
	}
}

func TestExpectedUnitsPromotesOrganelles(t *testing.T) {
	// The pepper shape: many sequences, plus MT and Pltd short enough to fall
	// outside the top 21 and be promoted back in.
	seqs := append(dotted(30),
		SeqInfo{ID: "MT", Len: 500_000},
		SeqInfo{ID: "Pltd", Len: 150_000},
	)

	units, err := ExpectedUnits(writeDict(t, seqs))
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	// 21 longest + MT + Pltd called individually, everything else batched.
	if len(units) != 24 {
		t.Fatalf("got %d units, want 24", len(units))
	}

	var sawMT, sawPltd bool
	for _, u := range units {
		switch u.ID {
		case "MT":
			sawMT = true
		case "Pltd":
			sawPltd = true
		}
	}
	if !sawMT || !sawPltd {
		t.Errorf("MT promoted = %v, Pltd promoted = %v; both must be called individually", sawMT, sawPltd)
	}
}

func TestExpectedUnitsBatchLengthIsTheSum(t *testing.T) {
	seqs := dotted(24)
	units, err := ExpectedUnits(writeDict(t, seqs))
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	batch := units[len(units)-1]
	if batch.Kind != UnitContigs {
		t.Fatalf("last unit kind = %q, want %q", batch.Kind, UnitContigs)
	}

	want := 0
	for _, m := range batch.Members {
		want += m.Len
	}
	if batch.Length != want {
		t.Errorf("batch length = %d, want %d (the sum of its members)", batch.Length, want)
	}
}

// TestExpectedUnitsAgreesWithPathBuilders is the reason ExpectedUnits lives in
// this package. If anyone changes how a label is spelled in a filename, this
// fails rather than silently teaching an auditor to look for the wrong names.
func TestExpectedUnitsAgreesWithPathBuilders(t *testing.T) {
	seqs := append(dotted(30), SeqInfo{ID: "MT", Len: 500_000})
	units, err := ExpectedUnits(writeDict(t, seqs))
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	sample := SampleWork{
		Sample:  "S1",
		Cram:    "/data/pepo/2026/S1/reference_genomes/v4.1/bams/S1.RGMD_bqsr.cram",
		CramDir: "/data/pepo/2026/S1/reference_genomes",
	}
	opts := Options{DataDir: "/data", Species: "pepo", RefVer: "v4.1", Caller: "gatk", Merger: "gatk"}

	for _, u := range units {
		gvcf := filepath.Base(GvcfPath(opts, sample, u.ID))
		if want := "." + u.Label + ".g.vcf.gz"; !strings.HasSuffix(gvcf, want) {
			t.Errorf("GvcfPath(%q) = %q, want suffix %q", u.ID, gvcf, want)
		}

		joint := filepath.Base(JointVcfPath(opts, u.ID))
		if want := "." + u.Label + ".joint.vcf.gz"; !strings.HasSuffix(joint, want) {
			t.Errorf("JointVcfPath(%q) = %q, want suffix %q", u.ID, joint, want)
		}
	}
}

func TestGvcfDirNamesMatchGvcfPath(t *testing.T) {
	sample := SampleWork{
		Sample:  "S1",
		Cram:    "/data/pepo/2026/S1/reference_genomes/v4.1/bams/S1.RGMD.cram",
		CramDir: "/data/pepo/2026/S1/reference_genomes",
	}

	// Every name the standard claims must be one GvcfPath actually writes.
	got := map[string]bool{}
	for _, caller := range []string{"gatk", "deepvariant"} {
		opts := Options{DataDir: "/data", Species: "pepo", RefVer: "v4.1", Caller: caller}
		got[filepath.Base(filepath.Dir(GvcfPath(opts, sample, "Cp4.1LG01")))] = true
	}

	for _, name := range GvcfDirNames() {
		if !got[name] {
			t.Errorf("GvcfDirNames lists %q, but no caller makes GvcfPath produce it", name)
		}
	}
	if len(got) != len(GvcfDirNames()) {
		t.Errorf("GvcfPath produces %v, GvcfDirNames claims %v", got, GvcfDirNames())
	}
}

func TestCallerMergerTagsMatchJointVcfDir(t *testing.T) {
	// The three combinations the pipeline can be asked for.
	combos := []struct{ caller, merger string }{
		{"gatk", "gatk"},
		{"gatk", "glnexus"},
		{"deepvariant", "glnexus"},
	}

	got := map[string]bool{}
	for _, c := range combos {
		opts := Options{DataDir: "/data", Species: "pepo", RefVer: "v4.1", Caller: c.caller, Merger: c.merger}
		got[filepath.Base(JointVcfDir(opts))] = true
	}

	for _, tag := range CallerMergerTags() {
		if !got[tag] {
			t.Errorf("CallerMergerTags lists %q, but no combination makes JointVcfDir produce it", tag)
		}
	}
	if len(got) != len(CallerMergerTags()) {
		t.Errorf("JointVcfDir produces %v, CallerMergerTags claims %v", got, CallerMergerTags())
	}
}

func TestExpectedUnitsMissingDict(t *testing.T) {
	if _, err := ExpectedUnits(filepath.Join(t.TempDir(), "absent.dict")); err == nil {
		t.Fatal("ExpectedUnits on a missing dict returned no error; a scanner must be able to tell resolved from missing")
	}
}

// ---------------------------------------------------------------------------
// On-disk truth. These skip when the reference volume is not mounted, following
// the pattern in alignmentdir/scatter_test.go. They are the tests that pin the
// rule against what the pipeline actually wrote, rather than against a fixture
// built from the same assumptions.
// ---------------------------------------------------------------------------

func TestExpectedUnitsMatchesPepperOnDisk(t *testing.T) {
	const dict = "/mnt/z/genomes/pepper/UCD10Xv1.1/assembly/GCF_002878395.1_UCD10Xv1.1_genomic_corrected.dict"
	if _, err := os.Stat(dict); err != nil {
		t.Skip("pepper UCD10Xv1.1 dict not mounted")
	}

	units, err := ExpectedUnits(dict)
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	// 81,202 sequences in the dict; 24 units called. The gap is the whole
	// point of the rule.
	if len(units) != 24 {
		t.Fatalf("got %d units, want 24", len(units))
	}

	byID := map[string]Unit{}
	for _, u := range units {
		byID[u.ID] = u
	}
	for _, want := range []string{"chr1", "chr12", "MT", "Pltd", BatchedLabel} {
		if _, ok := byID[want]; !ok {
			t.Errorf("unit %q missing from the expected set", want)
		}
	}
	if batch := byID[BatchedLabel]; len(batch.Members) < 80_000 {
		t.Errorf("batch holds %d members, want the ~81k that are not called individually", len(batch.Members))
	}
}

func TestExpectedUnitsMatchesPepoOnDisk(t *testing.T) {
	const dict = "/mnt/z/genomes/pepo/v4.1/assembly/Cpepo_genome_v4.1.dict"
	if _, err := os.Stat(dict); err != nil {
		t.Skip("pepo v4.1 dict not mounted")
	}

	units, err := ExpectedUnits(dict)
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	// Exactly 21 sequences, so every one is called individually and there is
	// no batch. A scanner expecting a "contigs" gVCF here would report every
	// sample as one short.
	if len(units) != 21 {
		t.Fatalf("got %d units, want 21", len(units))
	}
	for _, u := range units {
		if u.Kind == UnitContigs {
			t.Fatalf("pepo v4.1 produced a %q unit; it has only 21 sequences", BatchedLabel)
		}
		// Every pepo seqid is dotted, which is what makes the two filename
		// spellings diverge between GATK and DeepVariant.
		if !strings.Contains(u.ID, ".") || strings.Contains(u.Label, ".") {
			t.Errorf("unit {ID:%q Label:%q}: want a dotted ID and an undotted label", u.ID, u.Label)
		}
	}
}

// TestPepperGvcfFilenamesResolveToUnits reads the gVCF filenames the pipeline
// actually wrote for one pepper sample and checks every one maps to a unit.
func TestPepperGvcfFilenamesResolveToUnits(t *testing.T) {
	const (
		dict    = "/mnt/z/genomes/pepper/UCD10Xv1.1/assembly/GCF_002878395.1_UCD10Xv1.1_genomic_corrected.dict"
		gvcfDir = "/mnt/u/DATA/pepper/2025/HP2831/reference_genomes/UCD10Xv1.1/gvcfs"
	)
	if _, err := os.Stat(dict); err != nil {
		t.Skip("pepper dict not mounted")
	}
	entries, err := os.ReadDir(gvcfDir)
	if err != nil {
		t.Skip("pepper sample gVCF directory not mounted")
	}

	units, err := ExpectedUnits(dict)
	if err != nil {
		t.Fatalf("ExpectedUnits: %v", err)
	}

	// Both spellings, because GATK sanitises dots and DeepVariant does not.
	spellings := map[string]string{}
	for _, u := range units {
		spellings[u.Label] = u.ID
		spellings[u.ID] = u.ID
	}

	matched := map[string]bool{}
	for _, e := range entries {
		name := e.Name()
		if !strings.HasSuffix(name, ".g.vcf.gz") {
			continue
		}
		stem := strings.TrimSuffix(name, ".g.vcf.gz")

		found := ""
		for spelling, id := range spellings {
			if strings.HasSuffix(stem, "."+spelling) {
				// Longest wins, so chr1 cannot claim a chr12 file.
				if len(spelling) > len(found) {
					found, matched[id] = spelling, true
				}
			}
		}
		if found == "" {
			t.Errorf("gVCF %q matched no expected unit", name)
		}
	}

	if len(matched) != len(units) {
		t.Errorf("matched %d of %d units on disk", len(matched), len(units))
	}
}
