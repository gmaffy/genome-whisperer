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

// ---------------------------------------------------------------------------
// GenomicsDB workspace reuse
// ---------------------------------------------------------------------------

// completeWorkspace fakes what GenomicsDBImport leaves behind on success: the
// workspace with its metadata, and the marker holding the input fingerprint.
func completeWorkspace(t *testing.T, workDir, fingerprint string) (theDB, marker string) {
	t.Helper()
	theDB = filepath.Join(workDir, "db")
	if err := os.MkdirAll(theDB, 0755); err != nil {
		t.Fatal(err)
	}
	for _, meta := range []string{"callset.json", "vidmap.json"} {
		if err := os.WriteFile(filepath.Join(theDB, meta), []byte("{}"), 0644); err != nil {
			t.Fatal(err)
		}
	}
	marker = filepath.Join(workDir, gdbImportDoneFile)
	if err := os.WriteFile(marker, []byte(fingerprint), 0644); err != nil {
		t.Fatal(err)
	}
	return theDB, marker
}

func TestGdbReuseAcceptsCompleteWorkspace(t *testing.T) {
	theDB, marker := completeWorkspace(t, t.TempDir(), "abc123")

	reuse, reason := gdbReuse(theDB, marker, "abc123")
	if !reuse {
		t.Errorf("a complete workspace should be reused, got reason %q", reason)
	}
	if reason != "" {
		t.Errorf("nothing is being discarded, got reason %q", reason)
	}
}

// The marker is written only after GenomicsDBImport returns. An import killed
// part-way leaves a workspace holding an unknown subset of the cohort, which
// would silently genotype fewer samples than asked for.
func TestGdbReuseRejectsIncompleteImport(t *testing.T) {
	workDir := t.TempDir()
	theDB, marker := completeWorkspace(t, workDir, "abc123")
	if err := os.Remove(marker); err != nil {
		t.Fatal(err)
	}

	reuse, reason := gdbReuse(theDB, marker, "abc123")
	if reuse {
		t.Fatal("a workspace with no completion marker must be rebuilt")
	}
	if reason == "" {
		t.Error("an existing workspace being discarded should say why")
	}
}

// A workspace built before a sample was added, or before a gVCF was re-created,
// no longer matches the gVCFs on disk.
func TestGdbReuseRejectsStaleFingerprint(t *testing.T) {
	theDB, marker := completeWorkspace(t, t.TempDir(), "old-cohort")

	reuse, reason := gdbReuse(theDB, marker, "new-cohort")
	if reuse {
		t.Fatal("a workspace built from different inputs must be rebuilt")
	}
	if !strings.Contains(reason, "different") {
		t.Errorf("reason should name the mismatch, got %q", reason)
	}
}

func TestGdbReuseRejectsWorkspaceMissingMetadata(t *testing.T) {
	workDir := t.TempDir()
	theDB, marker := completeWorkspace(t, workDir, "abc123")
	if err := os.Remove(filepath.Join(theDB, "vidmap.json")); err != nil {
		t.Fatal(err)
	}

	reuse, reason := gdbReuse(theDB, marker, "abc123")
	if reuse {
		t.Fatal("a workspace missing its metadata must be rebuilt")
	}
	if !strings.Contains(reason, "vidmap.json") {
		t.Errorf("reason should name the missing file, got %q", reason)
	}
}

// A first run has nothing to discard and should say nothing about it.
func TestGdbReuseSilentWhenNoWorkspace(t *testing.T) {
	dir := t.TempDir()
	reuse, reason := gdbReuse(filepath.Join(dir, "db"), filepath.Join(dir, gdbImportDoneFile), "abc123")
	if reuse {
		t.Error("there is no workspace to reuse")
	}
	if reason != "" {
		t.Errorf("nothing existed, so nothing is being discarded; got %q", reason)
	}
}

func TestGdbFingerprintChangesWithInputs(t *testing.T) {
	dir := t.TempDir()
	a := writeVCF(t, filepath.Join(dir, "s1.g.vcf"), []string{"S1"})
	b := writeVCF(t, filepath.Join(dir, "s2.g.vcf"), []string{"S2"})
	seqs := []SeqInfo{{ID: "A01", Len: 100}}

	base, err := gdbFingerprint([]string{a}, seqs)
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}

	same, err := gdbFingerprint([]string{a}, seqs)
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}
	if same != base {
		t.Error("the same inputs must fingerprint the same, or a workspace is never reusable")
	}

	added, err := gdbFingerprint([]string{a, b}, seqs)
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}
	if added == base {
		t.Error("adding a sample must change the fingerprint")
	}

	widened, err := gdbFingerprint([]string{a}, append(seqs, SeqInfo{ID: "A02", Len: 50}))
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}
	if widened == base {
		t.Error("changing the imported intervals must change the fingerprint")
	}
}

// A gVCF re-created since the import holds different calls under the same path.
func TestGdbFingerprintChangesWhenGvcfIsRewritten(t *testing.T) {
	dir := t.TempDir()
	p := writeVCF(t, filepath.Join(dir, "s1.g.vcf"), []string{"S1"})
	seqs := []SeqInfo{{ID: "A01", Len: 100}}

	before, err := gdbFingerprint([]string{p}, seqs)
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}

	writeVCF(t, p, []string{"S1"}, "A01\t100\t.\tA\tT\t50\t.\tQD=2.0\tGT:GQ\t0/1:40")

	after, err := gdbFingerprint([]string{p}, seqs)
	if err != nil {
		t.Fatalf("gdbFingerprint: %v", err)
	}
	if after == before {
		t.Error("a re-created gVCF must invalidate the workspace built from the old one")
	}
}

func TestGdbFingerprintFailsOnMissingGvcf(t *testing.T) {
	_, err := gdbFingerprint([]string{filepath.Join(t.TempDir(), "gone.g.vcf")}, []SeqInfo{{ID: "A01"}})
	if err == nil {
		t.Error("expected an error when a gVCF cannot be stat'ed")
	}
}

// ---------------------------------------------------------------------------
// Merge parallelism
// ---------------------------------------------------------------------------

func TestMergeJobsHonoursExplicitSetting(t *testing.T) {
	jobs, _ := mergeJobs(Options{MergeJobs: 6, Threads: 4}, "gatk", 10)
	if jobs != 6 {
		t.Errorf("--merge-jobs 6 should give 6 jobs, got %d", jobs)
	}
}

// Starting more workers than there are chromosome groups just leaves them idle.
func TestMergeJobsNeverExceedsGroupCount(t *testing.T) {
	jobs, _ := mergeJobs(Options{MergeJobs: 16, Threads: 1}, "gatk", 3)
	if jobs != 3 {
		t.Errorf("3 groups cannot use more than 3 jobs, got %d", jobs)
	}
}

func TestMergeJobsAlwaysAtLeastOne(t *testing.T) {
	for _, opts := range []Options{
		{MergeJobs: 0, Threads: 1024}, // more threads per job than the machine has cores
		{MergeJobs: -1},
		{},
	} {
		if jobs, _ := mergeJobs(opts, "gatk", 5); jobs < 1 {
			t.Errorf("%+v gave %d jobs; the merge must always run", opts, jobs)
		}
	}
	if jobs, _ := mergeJobs(Options{}, "gatk", 0); jobs < 1 {
		t.Errorf("no groups still needs a valid worker count, got %d", jobs)
	}
}

// GLnexus is handed this figure, and glnexus_cli rejects a non-positive
// --mem-gbytes. It must stay usable even on a machine whose free memory cannot
// be read.
func TestMergeJobsMemoryShareIsAlwaysUsable(t *testing.T) {
	for _, groups := range []int{0, 1, 24} {
		_, memGB := mergeJobs(Options{Threads: 4}, "glnexus", groups)
		if memGB < glnexusJobMinMemGB {
			t.Errorf("%d groups gave %d GB per job, want at least %d", groups, memGB, glnexusJobMinMemGB)
		}
	}
}

// A job's heaps must fit inside the share it was given, or running the jobs
// concurrently oversubscribes the machine — the whole point of deriving them.
func TestGatkMergeHeapsFitWithinTheJobShare(t *testing.T) {
	for _, share := range []int{8, 9, 12, 16, 24, 29, 64} {
		importGB, genotypeGB := gatkMergeHeaps(share)
		if genotypeGB > share {
			t.Errorf("share %d GB: genotype heap %dg exceeds it", share, genotypeGB)
		}
		if importGB > share {
			t.Errorf("share %d GB: import heap %dg exceeds it", share, importGB)
		}
		// Headroom above the larger heap for metaspace, thread stacks and
		// GenomicsDB's off-heap buffers.
		if genotypeGB == share {
			t.Errorf("share %d GB: genotype heap %dg leaves no headroom", share, genotypeGB)
		}
	}
}

// The old hardcoded 8g/12g are the ceiling: they were sized for a large cohort,
// and a bigger heap on an idle machine buys nothing but longer GC pauses.
func TestGatkMergeHeapsCappedAtTheOldFixedSizes(t *testing.T) {
	importGB, genotypeGB := gatkMergeHeaps(1024)
	if importGB != 8 || genotypeGB != 12 {
		t.Errorf("got import %dg genotype %dg, want the former fixed 8g/12g", importGB, genotypeGB)
	}
}

// An oversubscribed --merge-jobs, or a host whose free memory cannot be read,
// must still produce heaps GATK will start with.
func TestGatkMergeHeapsFlooredForTinyShares(t *testing.T) {
	for _, share := range []int{0, 1, 2, 4} {
		importGB, genotypeGB := gatkMergeHeaps(share)
		if importGB < 2 || genotypeGB < 4 {
			t.Errorf("share %d GB gave import %dg genotype %dg; too small to run", share, importGB, genotypeGB)
		}
	}
}

// The 20-core / 32 GB machine this was tuned on: the old 16 GB budget allowed a
// single job, serialising 21 chromosomes.
func TestMergeJobsMemoryShareSizesUsableHeaps(t *testing.T) {
	jobs, memGB := mergeJobs(Options{Threads: 4}, "gatk", 21)
	importGB, genotypeGB := gatkMergeHeaps(memGB)

	if memGB < gatkMergeJobMemGB {
		t.Errorf("share %d GB is below the %d GB floor a job is budgeted", memGB, gatkMergeJobMemGB)
	}
	if jobs*genotypeGB < genotypeGB {
		t.Fatal("nonsensical job count")
	}
	if importGB < 2 || genotypeGB < 4 {
		t.Errorf("%d job(s) of %d GB gave unusable heaps: import %dg genotype %dg",
			jobs, memGB, importGB, genotypeGB)
	}
}
