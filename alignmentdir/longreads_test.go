package alignmentdir

import (
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// touchFile writes an empty file, creating its directory.
func touchFile(t *testing.T, path string) {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("mkdir %s: %v", filepath.Dir(path), err)
	}
	if err := os.WriteFile(path, []byte("x"), 0o644); err != nil {
		t.Fatalf("write %s: %v", path, err)
	}
}

// dirNames lists the regular files left in a directory.
func dirNames(t *testing.T, dir string) []string {
	t.Helper()
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("read %s: %v", dir, err)
	}
	var names []string
	for _, e := range entries {
		if !e.IsDir() {
			names = append(names, e.Name())
		}
	}
	sort.Strings(names)
	return names
}

// ---------------------------------------------------------------------------
// Cleanup
// ---------------------------------------------------------------------------

func TestCleanupKeepsWhatThisRunDidNotProduce(t *testing.T) {
	// The case this exists for: a long-read sample whose alignment was produced
	// somewhere else, adopted as input, and — where no FASTQ was kept — the
	// only copy of that sample's reads. The previous rule kept four paths and
	// deleted every other regular file, which removed it.
	dir := t.TempDir()

	adopted := filepath.Join(dir, "MENINA_LONG_READS.aligned.cram")
	adoptedIdx := adopted + ".crai"
	sortedBam := filepath.Join(dir, "MO971_LR.sorted.bam")
	rgmdBam := filepath.Join(dir, "MO971_LR.RGMD.bam")
	rgmdBamIdx := rgmdBam + ".csi"
	rgmdCram := filepath.Join(dir, "MO971_LR.RGMD.cram")
	rgmdCramIdx := rgmdCram + ".crai"
	otherSample := filepath.Join(dir, "OTHER_SAMPLE.RGMD.bam")
	recal := filepath.Join(dir, "MO971_LR.recal.txt")
	notes := filepath.Join(dir, "provenance.md")

	for _, p := range []string{adopted, adoptedIdx, sortedBam, rgmdBam, rgmdBamIdx, rgmdCram, rgmdCramIdx, otherSample, recal, notes} {
		touchFile(t, p)
	}

	produced := []string{sortedBam, rgmdBam, rgmdCram}
	if err := cleanupSampleOutputs(dir, "MO971_LR", produced, rgmdCram, ""); err != nil {
		t.Fatalf("cleanupSampleOutputs: %v", err)
	}

	got := dirNames(t, dir)
	want := []string{
		"MENINA_LONG_READS.aligned.cram",      // adopted input — never ours to delete
		"MENINA_LONG_READS.aligned.cram.crai", // its index goes with it
		"MO971_LR.RGMD.cram",                  // the durable artefact
		"MO971_LR.RGMD.cram.crai",
		"MO971_LR.recal.txt",    // .txt is always kept
		"OTHER_SAMPLE.RGMD.bam", // another sample's file: a name we write, but not for this sample
		"provenance.md",         // nothing this pipeline writes
	}
	sort.Strings(want)
	if strings.Join(got, ",") != strings.Join(want, ",") {
		t.Errorf("after cleanup the directory holds\n  %v\nwant\n  %v", got, want)
	}
}

func TestCleanupRemovesItsOwnIntermediatesFromAnEarlierRun(t *testing.T) {
	// The other half of the rule: an interrupted earlier run's leftovers carry
	// names only this pipeline writes for this sample, so they stay collectable
	// even though this process did not write them. Without that, a sample could
	// never reach the "already complete" state again.
	dir := t.TempDir()

	rgmdCram := filepath.Join(dir, "ON1708.RGMD.cram")
	for _, p := range []string{
		filepath.Join(dir, "ON1708.sorted.bam"),
		filepath.Join(dir, "ON1708.sorted.bam.csi"),
		filepath.Join(dir, "ON1708.sorted.partial.bam"),
		filepath.Join(dir, "ON1708.RGMD.bam"),
		filepath.Join(dir, "ON1708.RGMD_bqsr.bam"),
		rgmdCram,
	} {
		touchFile(t, p)
	}

	// produced is empty: this run did nothing but clean up.
	if err := cleanupSampleOutputs(dir, "ON1708", nil, rgmdCram, ""); err != nil {
		t.Fatalf("cleanupSampleOutputs: %v", err)
	}

	if got := dirNames(t, dir); strings.Join(got, ",") != "ON1708.RGMD.cram" {
		t.Errorf("after cleanup the directory holds %v, want just the rgmd.cram", got)
	}
}

func TestCleanupRefusesWithoutAnRgmdCram(t *testing.T) {
	// An empty rgmdCramPath used to mean an empty keep-set, which turned
	// cleanup into "delete everything in this directory".
	dir := t.TempDir()
	touchFile(t, filepath.Join(dir, "MENINA_LONG_READS.aligned.cram"))
	touchFile(t, filepath.Join(dir, "MO971_LR.sorted.bam"))

	if err := cleanupSampleOutputs(dir, "MO971_LR", nil, "", ""); err == nil {
		t.Fatal("cleanupSampleOutputs succeeded with no rgmd.cram, want a refusal")
	}
	if got := len(dirNames(t, dir)); got != 2 {
		t.Errorf("refused cleanup left %d files, want both still present", got)
	}
}

func TestCleanupRefusesWhenTheRgmdCramIsGone(t *testing.T) {
	// A path that was set but never written means a step failed silently; the
	// files still on disk are the only way back.
	dir := t.TempDir()
	touchFile(t, filepath.Join(dir, "MO971_LR.sorted.bam"))

	missing := filepath.Join(dir, "MO971_LR.RGMD.cram")
	if err := cleanupSampleOutputs(dir, "MO971_LR", nil, missing, ""); err == nil {
		t.Fatal("cleanupSampleOutputs succeeded with a missing rgmd.cram, want a refusal")
	}
	if got := len(dirNames(t, dir)); got != 1 {
		t.Errorf("refused cleanup left %d files, want the sorted.bam still present", got)
	}
}

// ---------------------------------------------------------------------------
// Long-read FASTQ discovery
// ---------------------------------------------------------------------------

func TestGetReadsLongIgnoresTheForwardReadHeuristics(t *testing.T) {
	// A PacBio movie name ends in a six-digit timestamp, so
	// m64011_200521_073741.fastq.gz ends in "1.fastq.gz" and isFwd calls it a
	// forward read. Identifying the long-read file by elimination — the one
	// ClassifyRead cannot place — would find nothing here, and every real HiFi
	// sample would be reported as having no reads.
	dir := t.TempDir()
	movie := filepath.Join(dir, "m64011_200521_073741.fastq.gz")
	touchFile(t, movie)

	if ClassifyRead(filepath.Base(movie)) != ReadRoleForward {
		t.Fatalf("precondition changed: %s no longer classifies as a forward read, "+
			"so this test no longer guards what it was written for", filepath.Base(movie))
	}

	got, err := GetReadsLong(dir)
	if err != nil {
		t.Fatalf("GetReadsLong: %v", err)
	}
	if got != movie {
		t.Errorf("GetReadsLong = %q, want %q", got, movie)
	}
}

func TestGetReadsLong(t *testing.T) {
	tests := []struct {
		name    string
		files   []string
		want    string
		wantErr bool
	}{
		{
			name:  "one fastq beside files that are not reads",
			files: []string{"MO971_LR.hifi_reads.fastq.gz", "notes.txt", "md5sums"},
			want:  "MO971_LR.hifi_reads.fastq.gz",
		},
		{
			name:  "uncompressed fastq",
			files: []string{"reads.fq"},
			want:  "reads.fq",
		},
		{
			// Two movies, or two subsets of one. Choosing silently would drop
			// half the sample.
			name:    "two fastqs",
			files:   []string{"m64011_200521_073741.fastq.gz", "m64011_200522_081530.fastq.gz"},
			wantErr: true,
		},
		{
			name:    "no fastq at all",
			files:   []string{"notes.txt"},
			wantErr: true,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			for _, f := range tc.files {
				touchFile(t, filepath.Join(dir, f))
			}

			got, err := GetReadsLong(dir)
			if tc.wantErr {
				if err == nil {
					t.Fatalf("GetReadsLong = %q, want an error", got)
				}
				return
			}
			if err != nil {
				t.Fatalf("GetReadsLong: %v", err)
			}
			if got != filepath.Join(dir, tc.want) {
				t.Errorf("GetReadsLong = %q, want %q", got, tc.want)
			}
		})
	}
}

// ---------------------------------------------------------------------------
// Adoption
// ---------------------------------------------------------------------------

func TestAdoptableAlignmentIsInvisibleToTheWarehouse(t *testing.T) {
	// The whole warehouse-compatibility argument rests on this: adoption is an
	// alignmentdir decision, and the estate report must go on calling these
	// files exactly what it called them before.
	for _, name := range []string{
		"MENINA_LONG_READS.aligned.cram",
		"S1_LONG_READS.cram",
		"MO971.pbmm2.bam",
	} {
		if !adoptableAlignment(name) {
			t.Errorf("adoptableAlignment(%q) = false, want true", name)
		}
		if role := AlignmentRole(name); role != RoleOther {
			t.Errorf("AlignmentRole(%q) = %q, want %q — adoption must not change what the scanner reports",
				name, role, RoleOther)
		}
	}
}

func TestAdoptableAlignmentRejectsOurOwnNames(t *testing.T) {
	for _, name := range []string{
		"S1.sorted.bam",
		"S1.RGMD.bam",
		"S1.RGMD.cram",
		"S1.RGMD_bqsr.cram",
		// Written by the old single-sample pbmm2 path: a leftover
		// intermediate, not somebody's input.
		"S1.RG.bam",
		// Not alignments at all.
		"S1.g.vcf.gz",
		"README",
		"S1.RGMD.cram.crai",
	} {
		if adoptableAlignment(name) {
			t.Errorf("adoptableAlignment(%q) = true, want false", name)
		}
	}
}

func TestAdoptionCandidates(t *testing.T) {
	external := FileInfo{Path: "/d/bams/MENINA_LONG_READS.aligned.cram", Present: true}
	secondMovie := FileInfo{Path: "/d/bams/MENINA_LONG_READS_2.aligned.cram", Present: true}
	junk := FileInfo{Path: "/d/bams/notes.md", Present: true}

	tests := []struct {
		name  string
		state SampleBamState
		want  int
	}{
		{
			// Precedence is decided by slot, not by directory order. This is
			// what stops a genuine .RGMD.cram from losing to an external file
			// that happens to sort after it.
			name:  "a usable rgmd.cram outranks anything external",
			state: SampleBamState{RgmdCram: indexed(), OtherFiles: []FileInfo{external}},
			want:  0,
		},
		{
			name:  "a usable rgmd.bam outranks anything external",
			state: SampleBamState{RgmdBam: usable(), OtherFiles: []FileInfo{external}},
			want:  0,
		},
		{
			// The resume path: <sample>.sorted.bam is what adoption wrote last
			// run, so re-normalising the cram would be hours of repeated work.
			name:  "a usable sorted.bam outranks anything external",
			state: SampleBamState{SortedBam: usable(), OtherFiles: []FileInfo{external}},
			want:  0,
		},
		{
			name:  "nothing of our own, one external alignment",
			state: SampleBamState{OtherFiles: []FileInfo{external, junk}},
			want:  1,
		},
		{
			// Reported, never chosen between.
			name:  "two movies",
			state: SampleBamState{OtherFiles: []FileInfo{external, secondMovie}},
			want:  2,
		},
		{
			name:  "nothing adoptable",
			state: SampleBamState{OtherFiles: []FileInfo{junk}},
			want:  0,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if got := adoptionCandidates(tc.state); len(got) != tc.want {
				t.Errorf("adoptionCandidates = %v (%d), want %d candidates", got, len(got), tc.want)
			}
		})
	}
}

// ---------------------------------------------------------------------------
// Planning
// ---------------------------------------------------------------------------

func TestBuildPlanLongReadNeverRunsBQSR(t *testing.T) {
	// Asked for with every flag that would turn BQSR on for a short-read
	// sample, including the reference shape that adds the extra bam steps.
	plan, _ := buildPlan(SampleBamState{}, planInputs{bqsr: true, gatkNeedsBam: true, longRead: true})

	for _, forbidden := range []ProcessingReason{
		ReasonRunBQSR,
		ReasonConvertBqsrBam,
		ReasonIndexBqsrCram,
		ReasonMaterializeRgmdBam,
	} {
		if containsReason(plan, forbidden) {
			t.Errorf("long-read plan %v contains %s", plan, forbidden)
		}
	}
}

func TestBuildPlanPicksTheAlignerStepByKind(t *testing.T) {
	longPlan, longNeedsReads := buildPlan(SampleBamState{}, planInputs{longRead: true})
	if !containsReason(longPlan, ReasonAlignLongReads) {
		t.Errorf("long-read plan %v is missing %s", longPlan, ReasonAlignLongReads)
	}
	if containsReason(longPlan, ReasonAlignFromReads) {
		t.Errorf("long-read plan %v uses the paired-read aligner", longPlan)
	}
	if !longNeedsReads {
		t.Error("long-read plan without an existing alignment must ask for reads")
	}

	shortPlan, shortNeedsReads := buildPlan(SampleBamState{}, planInputs{})
	if !containsReason(shortPlan, ReasonAlignFromReads) {
		t.Errorf("short-read plan %v is missing %s", shortPlan, ReasonAlignFromReads)
	}
	if containsReason(shortPlan, ReasonAlignLongReads) {
		t.Errorf("short-read plan %v uses the long-read aligner", shortPlan)
	}
	if !shortNeedsReads {
		t.Error("short-read plan without an existing alignment must ask for reads")
	}
}

func TestBuildPlanAdoptedExternalAlignmentSkipsTheAligner(t *testing.T) {
	plan, needsReads := buildPlan(SampleBamState{}, planInputs{longRead: true, external: true})

	want := []ProcessingReason{
		ReasonAdoptExternalAlignment,
		ReasonMarkDupSortedBam,
		ReasonConvertRgmdBam,
		ReasonIndexRgmdCram,
		ReasonCleanup,
	}
	if len(plan) != len(want) {
		t.Fatalf("plan = %v, want %v", plan, want)
	}
	for i := range want {
		if plan[i] != want[i] {
			t.Fatalf("plan = %v, want %v", plan, want)
		}
	}

	// An adopted sample frequently has no FASTQ on disk at all, so asking for
	// reads would strand it.
	if needsReads {
		t.Error("an adopted sample must not require reads to be located")
	}
}

func TestBuildPlanAdoptionOutranksAligningButNotOurOwnSortedBam(t *testing.T) {
	// A sorted.bam is what a previous run's adoption produced. Re-normalising
	// the external cram instead would repeat hours of work for the same file.
	plan, _ := buildPlan(SampleBamState{SortedBam: usable()}, planInputs{longRead: true, external: true})
	if containsReason(plan, ReasonAdoptExternalAlignment) {
		t.Errorf("plan %v re-adopts even though a usable sorted.bam is present", plan)
	}
	if !containsReason(plan, ReasonMarkDupSortedBam) {
		t.Errorf("plan %v should resume at duplicate marking", plan)
	}
}

func TestBuildPlanLongReadIsCompleteAtTheRgmdCram(t *testing.T) {
	// No bqsr.cram is ever produced for a long-read sample, so a finished one
	// has nothing left to do but the cleanup pass.
	plan, needsReads := buildPlan(SampleBamState{RgmdCram: indexed()}, planInputs{longRead: true})
	if len(plan) != 1 || plan[0] != ReasonCleanup {
		t.Errorf("plan = %v, want just %s", plan, ReasonCleanup)
	}
	if needsReads {
		t.Error("a finished sample must not require reads")
	}
}

// ---------------------------------------------------------------------------
// Preset validation
// ---------------------------------------------------------------------------

func TestValidatePreset(t *testing.T) {
	// The flag's help text long read "CSS", a transposition of CCS, so a value
	// copied straight out of --help was invalid. Checking up front turns that
	// into an immediate error rather than one raised after the short-read
	// samples have already run.
	for _, ok := range []string{"HIFI", "hifi", "CCS", "SUBREAD", "ISOSEQ", "UNROLLED"} {
		if err := validatePreset(ok); err != nil {
			t.Errorf("validatePreset(%q) = %v, want nil", ok, err)
		}
	}
	for _, bad := range []string{"CSS", "", "PACBIO", "HiFi2"} {
		if err := validatePreset(bad); err == nil {
			t.Errorf("validatePreset(%q) = nil, want an error", bad)
		}
	}
}

func TestAdoptionIgnoresOurOwnInterruptedWrites(t *testing.T) {
	// Found by running the pipeline: a run killed mid-alignment leaves
	// <sample>.sorted.partial.bam, which is not in any slot and so arrives here
	// in OtherFiles. Adopting it would be wrong twice over — it is this
	// pipeline's own file, and it is truncated — and it would mean a killed run
	// silently "adopted" its own wreckage instead of re-aligning.
	state := SampleBamState{
		Sample: "MO971_LR",
		OtherFiles: []FileInfo{
			{Path: "/d/bams/MO971_LR.sorted.partial.bam", Present: true},
			{Path: "/d/bams/MO971_LR.RGMD.cram.partial", Present: true},
			{Path: "/d/bams/MO971_LR.RG.bam", Present: true},
		},
	}
	if got := adoptionCandidates(state); len(got) != 0 {
		t.Errorf("adoptionCandidates = %v, want none — these are all our own leavings", got)
	}

	// A genuine external alignment sitting beside them is still found.
	state.OtherFiles = append(state.OtherFiles, FileInfo{Path: "/d/bams/MENINA_LONG_READS.aligned.cram", Present: true})
	got := adoptionCandidates(state)
	if len(got) != 1 || filepath.Base(got[0]) != "MENINA_LONG_READS.aligned.cram" {
		t.Errorf("adoptionCandidates = %v, want just the external alignment", got)
	}
}
