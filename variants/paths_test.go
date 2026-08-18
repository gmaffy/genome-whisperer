package variants

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Path construction is where the original bugs lived: two output roots, two
// filenames for the concatenated VCF, and a "latest directory" scan that could
// return a chromosome directory. These tests pin every name the package writes.

func TestCallerMergerTag(t *testing.T) {
	cases := []struct {
		caller, merger, want string
	}{
		{"gatk", "gatk", "gatk_gatk"},
		{"gatk", "glnexus", "gatk_glnexus"},
		{"deepvariant", "glnexus", "dv_glnexus"},
		// Case from the command line must not leak into the path.
		{"GATK", "GLnexus", "gatk_glnexus"},
		{"DeepVariant", "GLNEXUS", "dv_glnexus"},
	}
	for _, c := range cases {
		got := callerMergerTag(Options{Caller: c.caller, Merger: c.merger})
		if got != c.want {
			t.Errorf("callerMergerTag(%q, %q) = %q, want %q", c.caller, c.merger, got, c.want)
		}
	}
}

func TestGvcfPathDataDirMode(t *testing.T) {
	sample := SampleWork{
		Sample:  "S1",
		Cram:    "/data/cotton/proj/S1/reference_genomes/AD1.1/bams/S1.bqsr.cram",
		CramDir: "/data/cotton/proj/S1/reference_genomes",
	}
	opts := Options{DataDir: "/data", Species: "cotton", RefVer: "AD1.1", Caller: "gatk"}

	got := GvcfPath(opts, sample, "A01")
	want := filepath.Join("/data/cotton/proj/S1/reference_genomes", "AD1.1", "gatk_gvcfs", "S1.bqsr.A01.g.vcf.gz")
	if got != want {
		t.Errorf("GvcfPath =\n  %s\nwant\n  %s", got, want)
	}

	// DeepVariant writes to the sibling directory.
	opts.Caller = "deepvariant"
	if got := GvcfPath(opts, sample, "A01"); !strings.Contains(got, "dv_gvcfs") {
		t.Errorf("deepvariant should use dv_gvcfs, got %s", got)
	}
}

func TestGvcfPathOutDirMode(t *testing.T) {
	sample := SampleWork{Sample: "S1", Cram: "/in/S1.bqsr.bam"}
	opts := Options{OutDir: "/out", Species: "cotton", RefVer: "AD1.1", Caller: "gatk"}

	got := GvcfPath(opts, sample, "A01")
	want := filepath.Join("/out", "A01", "gatk_gvcfs", "S1.bqsr.A01.g.vcf.gz")
	if got != want {
		t.Errorf("GvcfPath =\n  %s\nwant\n  %s", got, want)
	}
}

// A chromosome ID containing a dot must be sanitised consistently in both the
// directory and the filename, or the reuse check looks for a file that was never
// written under that name.
func TestGvcfPathSanitisesChromLabel(t *testing.T) {
	sample := SampleWork{Sample: "S1", Cram: "/in/S1.bqsr.cram"}
	opts := Options{OutDir: "/out", Species: "cotton", RefVer: "AD1.1", Caller: "gatk"}

	got := GvcfPath(opts, sample, "Chr1.1")
	if strings.Contains(filepath.Base(got), "Chr1.1") {
		t.Errorf("dot in chromosome ID should be replaced: %s", got)
	}
	if !strings.HasSuffix(got, "S1.bqsr.Chr1_1.g.vcf.gz") {
		t.Errorf("unexpected filename: %s", got)
	}
	if !strings.Contains(got, filepath.Join("out", "Chr1_1")) {
		t.Errorf("directory should use the sanitised label too: %s", got)
	}
}

// The base name comes from TrimSuffix(base, Ext), not a chained
// strings.Replace(".cram")/.Replace(".bam"), which silently produced a wrong path
// for an alignment matching neither.
func TestGvcfPathTrimsAlignmentExtension(t *testing.T) {
	opts := Options{OutDir: "/out", Species: "cotton", RefVer: "AD1.1", Caller: "gatk"}
	for _, aln := range []string{"/in/S1.bqsr.cram", "/in/S1.bqsr.bam", "/in/S1.bqsr.sam"} {
		got := filepath.Base(GvcfPath(opts, SampleWork{Sample: "S1", Cram: aln}, "A01"))
		if got != "S1.bqsr.A01.g.vcf.gz" {
			t.Errorf("%s -> %s, want S1.bqsr.A01.g.vcf.gz", aln, got)
		}
	}
}

func TestJointVcfDirIncludesCombination(t *testing.T) {
	dd := Options{DataDir: "/data", Species: "Cotton", RefVer: "AD1.1", Caller: "gatk", Merger: "glnexus"}
	want := filepath.Join("/data", "cotton", "MERGED_VCFs", "AD1.1", "gatk_glnexus")
	if got := JointVcfDir(dd); got != want {
		t.Errorf("JointVcfDir =\n  %s\nwant\n  %s", got, want)
	}

	od := Options{OutDir: "/out", Species: "Cotton", RefVer: "AD1.1", Caller: "deepvariant", Merger: "glnexus"}
	want = filepath.Join("/out", "VCFs", "dv_glnexus")
	if got := JointVcfDir(od); got != want {
		t.Errorf("JointVcfDir =\n  %s\nwant\n  %s", got, want)
	}
}

// mergeWorkDir and gatkScratchDir are rooted at opts.LocalWorkDir / the OS temp
// dir, not JointVcfDir, since GenomicsDBImport and glnexus_cli do small-file
// random I/O that performs badly on the NAS/NFS storage JointVcfDir usually
// lives on.
func TestMergeWorkDirAndGatkScratchDirAreLocal(t *testing.T) {
	opts := Options{DataDir: "/data", LocalWorkDir: "/home/user/.genome-whisperer/merge-work",
		Species: "Cotton", RefVer: "AD1.1", Caller: "gatk", Merger: "glnexus"}

	want := filepath.Join("/home/user/.genome-whisperer/merge-work", "cotton", "AD1.1", "gatk_glnexus", "work", "A01")
	if got := mergeWorkDir(opts, "A01"); got != want {
		t.Errorf("mergeWorkDir =\n  %s\nwant\n  %s", got, want)
	}
	if got := mergeWorkDir(opts, "A01"); strings.HasPrefix(got, JointVcfDir(opts)) {
		t.Errorf("mergeWorkDir must not live under JointVcfDir: %s", got)
	}

	wantTmp := filepath.Join(os.TempDir(), "genome-whisperer-merge", "cotton", "gatk_glnexus", "A01")
	if got := gatkScratchDir(opts, "A01"); got != wantTmp {
		t.Errorf("gatkScratchDir =\n  %s\nwant\n  %s", got, wantTmp)
	}
}

// Concurrent merge jobs run one chromosome group each on a shared LocalWorkDir
// (see mergeJobs), so every chromosome and every caller/merger combination must
// get its own scratch directory or two jobs would collide mid-import.
func TestMergeWorkDirAndGatkScratchDirArePerChromAndCombination(t *testing.T) {
	base := Options{LocalWorkDir: "/home/user/.genome-whisperer/merge-work", Species: "cotton", RefVer: "AD1.1"}

	seenWork := make(map[string]bool)
	seenScratch := make(map[string]bool)
	for _, c := range []struct{ caller, merger string }{
		{"gatk", "gatk"}, {"gatk", "glnexus"}, {"deepvariant", "glnexus"},
	} {
		opts := base
		opts.Caller, opts.Merger = c.caller, c.merger
		for _, chrom := range []string{"A01", "A02"} {
			if w := mergeWorkDir(opts, chrom); seenWork[w] {
				t.Errorf("duplicate mergeWorkDir: %s", w)
			} else {
				seenWork[w] = true
			}
			if s := gatkScratchDir(opts, chrom); seenScratch[s] {
				t.Errorf("duplicate gatkScratchDir: %s", s)
			} else {
				seenScratch[s] = true
			}
		}
	}
	if len(seenWork) != 6 || len(seenScratch) != 6 {
		t.Errorf("expected 6 distinct paths each, got %d work, %d scratch", len(seenWork), len(seenScratch))
	}
}

// Each combination must land in its own directory. Sharing one directory let a
// GATK-merged VCF be reused for a GLnexus run: the sample names matched, so the
// reuse check passed and the wrong file was accepted silently.
func TestEachCombinationGetsItsOwnPath(t *testing.T) {
	seen := make(map[string]string)
	for _, c := range []struct{ caller, merger string }{
		{"gatk", "gatk"}, {"gatk", "glnexus"}, {"deepvariant", "glnexus"},
	} {
		opts := Options{DataDir: "/data", Species: "cotton", RefVer: "AD1.1", Caller: c.caller, Merger: c.merger}
		for _, p := range []string{JointVcfPath(opts, "A01"), FinalVcfPath(opts)} {
			if prev, dup := seen[p]; dup {
				t.Errorf("%s/%s collides with %s at %s", c.caller, c.merger, prev, p)
			}
			seen[p] = c.caller + "/" + c.merger
		}
	}
	if len(seen) != 6 {
		t.Errorf("expected 6 distinct paths, got %d", len(seen))
	}
}

// The combination sits before the chromosome so the trailing
// ".<chrom>.joint.vcf.gz" and ".all.vcf.gz" stay intact for suffix and glob
// checks.
func TestVcfNameSuffixesArePreserved(t *testing.T) {
	opts := Options{DataDir: "/data", Species: "cotton", RefVer: "AD1.1", Caller: "gatk", Merger: "glnexus"}

	joint := filepath.Base(JointVcfPath(opts, "A01"))
	if joint != "COTTON.ad1.1.gatk_glnexus.A01.joint.vcf.gz" {
		t.Errorf("joint name = %s", joint)
	}
	if !strings.HasSuffix(joint, ".A01.joint.vcf.gz") {
		t.Errorf("per-chromosome tail must stay .<chrom>.joint.vcf.gz: %s", joint)
	}

	final := filepath.Base(FinalVcfPath(opts))
	if final != "COTTON.ad1.1.gatk_glnexus.all.vcf.gz" {
		t.Errorf("final name = %s", final)
	}
	if !strings.HasSuffix(final, ".all.vcf.gz") {
		t.Errorf("final tail must stay .all.vcf.gz: %s", final)
	}
}

func TestAbsRoots(t *testing.T) {
	got, err := absRoots(Options{OutDir: "vcfs", RefFasta: filepath.Join("ref", "ref.fa")})
	if err != nil {
		t.Fatalf("absRoots: %v", err)
	}
	if !filepath.IsAbs(got.OutDir) {
		t.Errorf("OutDir not absolute: %q", got.OutDir)
	}
	if !filepath.IsAbs(got.RefFasta) {
		t.Errorf("RefFasta not absolute: %q", got.RefFasta)
	}
	// An unset root must stay unset: the mode checks distinguish the two by
	// emptiness, and filepath.Abs("") would turn it into the working directory.
	if got.DataDir != "" {
		t.Errorf("empty DataDir became %q", got.DataDir)
	}
	// Unlike DataDir/OutDir, an unset LocalWorkDir gets a real default rather
	// than staying empty: it is never optional the way "no data-dir mode" is.
	if !filepath.IsAbs(got.LocalWorkDir) {
		t.Errorf("LocalWorkDir not defaulted to an absolute path: %q", got.LocalWorkDir)
	}
}

// An explicit --local-work-dir must still work when $HOME is not what the
// caller wants (e.g. it turns out to be network storage too).
func TestAbsRootsRespectsExplicitLocalWorkDir(t *testing.T) {
	got, err := absRoots(Options{LocalWorkDir: filepath.Join("scratch", "work")})
	if err != nil {
		t.Fatalf("absRoots: %v", err)
	}
	if !filepath.IsAbs(got.LocalWorkDir) || filepath.Base(got.LocalWorkDir) != "work" {
		t.Errorf("LocalWorkDir = %q, want an absolute path ending in .../scratch/work", got.LocalWorkDir)
	}
}

func TestAbsRootsIsIdempotent(t *testing.T) {
	once, err := absRoots(Options{DataDir: "data", OutDir: "out", RefFasta: "ref.fa"})
	if err != nil {
		t.Fatal(err)
	}
	twice, err := absRoots(once)
	if err != nil {
		t.Fatal(err)
	}
	// Options holds a slice, so compare the fields absRoots actually touches.
	if once.DataDir != twice.DataDir || once.OutDir != twice.OutDir || once.RefFasta != twice.RefFasta ||
		once.LocalWorkDir != twice.LocalWorkDir {
		t.Errorf("absRoots not idempotent:\n  %s %s %s %s\n  %s %s %s %s",
			once.DataDir, once.OutDir, once.RefFasta, once.LocalWorkDir,
			twice.DataDir, twice.OutDir, twice.RefFasta, twice.LocalWorkDir)
	}
}

// ---------------------------------------------------------------------------
// Reference dictionary
// ---------------------------------------------------------------------------

// writeDict builds a .dict whose @SQ order deliberately differs from length
// order, so a function that sorts by length cannot pass a dictionary-order test
// by accident.
func writeDict(t *testing.T, seqs []SeqInfo) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.dict")

	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\n")
	for _, s := range seqs {
		b.WriteString("@SQ\tSN:" + s.ID + "\tLN:" + itoa(s.Len) + "\tUR:file:/ref.fa\n")
	}
	if err := os.WriteFile(path, []byte(b.String()), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	var d []byte
	for i > 0 {
		d = append([]byte{byte('0' + i%10)}, d...)
		i /= 10
	}
	return string(d)
}

func TestDictOrderFollowsFileOrderNotLength(t *testing.T) {
	// Shortest sequence first in the file.
	dict := writeDict(t, []SeqInfo{
		{ID: "A01", Len: 10},
		{ID: "A02", Len: 900},
		{ID: "A03", Len: 500},
	})

	order, err := dictOrder(dict)
	if err != nil {
		t.Fatalf("dictOrder: %v", err)
	}
	for id, want := range map[string]int{"A01": 0, "A02": 1, "A03": 2} {
		if order[id] != want {
			t.Errorf("order[%s] = %d, want %d", id, order[id], want)
		}
	}
	if len(order) != 3 {
		t.Errorf("expected 3 entries, got %d", len(order))
	}
}

// getChromsAndContigs sorts by length, which is why dictOrder has to read the
// dict itself. This pins the difference so nobody swaps one for the other.
func TestGetChromsAndContigsSortsByLengthNotDictOrder(t *testing.T) {
	dict := writeDict(t, []SeqInfo{
		{ID: "small", Len: 10},
		{ID: "big", Len: 1000},
	})

	chroms, contigs, err := getChromsAndContigs(dict)
	if err != nil {
		t.Fatalf("getChromsAndContigs: %v", err)
	}
	if len(contigs) != 0 {
		t.Errorf("with 2 sequences everything should be a chromosome, got %d contigs", len(contigs))
	}
	if len(chroms) != 2 || chroms[0].ID != "big" {
		t.Fatalf("expected length-descending order, got %+v", chroms)
	}

	order, err := dictOrder(dict)
	if err != nil {
		t.Fatal(err)
	}
	if order["small"] != 0 {
		t.Errorf("dictOrder should keep file order (small first), got %d", order["small"])
	}
}

// The split keeps the 21 longest sequences as chromosomes and always promotes MT
// and Pltd. Undocumented in the original code, so pinned here.
func TestGetChromsAndContigsTop21PlusOrganelles(t *testing.T) {
	var seqs []SeqInfo
	for i := 0; i < 30; i++ {
		seqs = append(seqs, SeqInfo{ID: "seq" + itoa(i), Len: 1000 - i})
	}
	seqs = append(seqs, SeqInfo{ID: "MT", Len: 1}, SeqInfo{ID: "Pltd", Len: 2})

	chroms, contigs, err := getChromsAndContigs(writeDict(t, seqs))
	if err != nil {
		t.Fatalf("getChromsAndContigs: %v", err)
	}

	in := func(list []SeqInfo, id string) bool {
		for _, s := range list {
			if s.ID == id {
				return true
			}
		}
		return false
	}

	// 21 longest, plus MT and Pltd pulled back out of the contigs.
	if len(chroms) != 23 {
		t.Errorf("expected 21 longest + MT + Pltd = 23 chromosomes, got %d", len(chroms))
	}
	for _, id := range []string{"MT", "Pltd"} {
		if !in(chroms, id) {
			t.Errorf("%s should be promoted to chroms", id)
		}
		if in(contigs, id) {
			t.Errorf("%s should not remain in contigs", id)
		}
	}
	if !in(chroms, "seq0") {
		t.Error("longest sequence should be a chromosome")
	}
	if !in(contigs, "seq29") {
		t.Error("shortest non-organelle sequence should be a contig")
	}
	if len(chroms)+len(contigs) != len(seqs) {
		t.Errorf("every sequence must appear exactly once: %d + %d != %d",
			len(chroms), len(contigs), len(seqs))
	}
}

func TestGetChromsAndContigsMissingDict(t *testing.T) {
	if _, _, err := getChromsAndContigs(filepath.Join(t.TempDir(), "nope.dict")); err == nil {
		t.Error("expected an error for a missing dict file")
	}
}
