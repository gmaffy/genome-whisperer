package alignmentdir

import (
	"os"
	"path/filepath"
	"testing"
)

func usable() FileInfo { return FileInfo{Present: true, Valid: true} }
func indexed() FileInfo {
	return FileInfo{Present: true, Valid: true, IndexPresent: true, IndexSize: 1}
}

// stepIndex returns where a step sits in a plan, or -1.
func stepIndex(plan []ProcessingReason, want ProcessingReason) int {
	for i, step := range plan {
		if step == want {
			return i
		}
	}
	return -1
}

func TestGatkReadsNeedBam(t *testing.T) {
	tests := []struct {
		name    string
		entries []faiEntry
		want    bool
	}{
		{"contig over the BAI ceiling", []faiEntry{{"chr1", baiMaxPosition + 1}}, true},
		{"contig exactly at the ceiling", []faiEntry{{"chr1", baiMaxPosition}}, false},
		{"one big contig among small ones", []faiEntry{{"chr1", 1000}, {"chr2", 1_056_413_696}, {"chr3", 1000}}, true},
		{"an ordinary genome", []faiEntry{{"chr1", 30_000_000}, {"chr2", 25_000_000}}, false},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			refFasta, _ := writeFai(t, t.TempDir(), tc.entries)
			if got := gatkReadsNeedBam(refFasta); got != tc.want {
				t.Fatalf("gatkReadsNeedBam = %v, want %v", got, tc.want)
			}
		})
	}
}

func TestGatkReadsNeedBamWithoutFai(t *testing.T) {
	refFasta := filepath.Join(t.TempDir(), "ref.fa")
	if err := os.WriteFile(refFasta, []byte(">chr1\nACGT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	// No .fai to read: the scatter path reports that far better than this does,
	// so it must not decide the route by itself.
	if gatkReadsNeedBam(refFasta) {
		t.Fatal("a reference with no .fai should not be reported as needing a bam")
	}
}

// On a reference GATK can read as CRAM, the plan must be exactly what it always
// was — the new steps are for the big-contig case only.
func TestBuildPlanUnchangedWhenCramIsReadable(t *testing.T) {
	state := SampleBamState{RgmdCram: indexed()}
	plan, _ := buildPlan(state, true, false)

	if stepIndex(plan, ReasonMaterializeRgmdBam) != -1 {
		t.Fatalf("plan should not decode the cram on a readable reference: %v", plan)
	}
	if stepIndex(plan, ReasonIndexRgmdBam) != -1 {
		t.Fatalf("plan should not index a bam it does not need: %v", plan)
	}
	if stepIndex(plan, ReasonRunBQSR) == -1 {
		t.Fatalf("plan should still run BQSR: %v", plan)
	}
}

func TestBuildPlanDecodesCramWhenItIsTheOnlyCopy(t *testing.T) {
	state := SampleBamState{RgmdCram: indexed()}
	plan, needsReads := buildPlan(state, true, true)

	materialize := stepIndex(plan, ReasonMaterializeRgmdBam)
	runBQSR := stepIndex(plan, ReasonRunBQSR)
	if materialize == -1 {
		t.Fatalf("cram-only sample needs the cram decoded first: %v", plan)
	}
	if materialize > runBQSR {
		t.Fatalf("the bam must exist before BQSR reads it: %v", plan)
	}
	if needsReads {
		t.Fatal("decoding a cram does not need the FASTQs")
	}
}

func TestBuildPlanIndexesAnExistingRgmdBam(t *testing.T) {
	// The bam is right there — decoding the cram again would be wasted hours.
	state := SampleBamState{RgmdCram: indexed(), RgmdBam: usable()}
	plan, _ := buildPlan(state, true, true)

	if got := stepIndex(plan, ReasonMaterializeRgmdBam); got != -1 {
		t.Fatalf("should not decode the cram when an rgmd.bam exists: %v", plan)
	}
	indexBam := stepIndex(plan, ReasonIndexRgmdBam)
	if indexBam == -1 || indexBam > stepIndex(plan, ReasonRunBQSR) {
		t.Fatalf("an unindexed rgmd.bam must be indexed before BQSR: %v", plan)
	}
}

func TestBuildPlanLeavesAnIndexedRgmdBamAlone(t *testing.T) {
	state := SampleBamState{RgmdCram: indexed(), RgmdBam: indexed()}
	plan, _ := buildPlan(state, true, true)

	if stepIndex(plan, ReasonMaterializeRgmdBam) != -1 || stepIndex(plan, ReasonIndexRgmdBam) != -1 {
		t.Fatalf("an indexed rgmd.bam needs no preparation: %v", plan)
	}
}

func TestBuildPlanIndexesTheBamItIsAboutToWrite(t *testing.T) {
	// Nothing on disk: the sample aligns from reads, and the rgmd.bam that
	// MarkDuplicates writes is what BQSR will read, so it needs an index.
	plan, needsReads := buildPlan(SampleBamState{}, true, true)

	if !needsReads {
		t.Fatal("aligning from scratch needs the FASTQs")
	}
	if stepIndex(plan, ReasonMaterializeRgmdBam) != -1 {
		t.Fatalf("no cram to decode when aligning from scratch: %v", plan)
	}
	indexBam := stepIndex(plan, ReasonIndexRgmdBam)
	if indexBam == -1 || indexBam > stepIndex(plan, ReasonRunBQSR) {
		t.Fatalf("the fresh rgmd.bam must be indexed before BQSR: %v", plan)
	}
}

func TestBuildPlanAddsNothingWithoutBQSR(t *testing.T) {
	plan, _ := buildPlan(SampleBamState{RgmdCram: indexed()}, false, true)

	if stepIndex(plan, ReasonMaterializeRgmdBam) != -1 || stepIndex(plan, ReasonIndexRgmdBam) != -1 {
		t.Fatalf("no BQSR means no reason to produce a bam: %v", plan)
	}
}

const onionSQ = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:CM060928.1\tLN:1056413696\tM5:693c4fe25d0a29c007961af5e001d8e0\tUR:/mnt/u/gone.fna.gz\n" +
	"@SQ\tSN:CM060930.1\tLN:275542313\tM5:d0bd96eab210454c0f849cfa3dca5c0c\tUR:/mnt/u/gone.fna.gz\n"

const renamedSQ = "@HD\tVN:1.6\n" +
	"@SQ\tSN:chr1_0-1056413696\tLN:1056413696\tM5:693c4fe25d0a29c007961af5e001d8e0\n" +
	"@SQ\tSN:chr1_2112827392-2388369705\tLN:275542313\tM5:d0bd96eab210454c0f849cfa3dca5c0c\n"

func TestSameSequencesSeesThroughARename(t *testing.T) {
	if !sameSequences(parseSQ(onionSQ), parseSQ(renamedSQ)) {
		t.Fatal("identical M5 checksums and lengths mean the same sequences, renamed")
	}
}

func TestSameSequencesRejectsDifferentAssemblies(t *testing.T) {
	other := "@SQ\tSN:chr1\tLN:1056413696\tM5:00000000000000000000000000000000\n" +
		"@SQ\tSN:chr2\tLN:275542313\tM5:11111111111111111111111111111111\n"
	if sameSequences(parseSQ(onionSQ), parseSQ(other)) {
		t.Fatal("different checksums are a different assembly, not a rename")
	}
	// A dictionary with no M5 to compare cannot be called a rename either.
	noM5 := "@SQ\tSN:chr1\tLN:1056413696\n@SQ\tSN:chr2\tLN:275542313\n"
	if sameSequences(parseSQ(onionSQ), parseSQ(noM5)) {
		t.Fatal("a missing M5 is not evidence of sameness")
	}
}

func TestParseSQKeepsFileOrder(t *testing.T) {
	got := parseSQ(onionSQ)
	if len(got) != 2 {
		t.Fatalf("expected 2 @SQ records, got %d", len(got))
	}
	if got[0].name != "CM060928.1" || got[0].length != "1056413696" || got[0].m5 == "" {
		t.Fatalf("first record parsed wrong: %+v", got[0])
	}
	if got[1].name != "CM060930.1" {
		t.Fatalf("@SQ order not preserved: %+v", got)
	}
}
