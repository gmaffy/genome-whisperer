package alignmentdir

import (
	"bytes"
	"strings"
	"testing"
)

func file(path string, f FileInfo) FileInfo {
	f.Path = path
	return f
}

func TestAssessSampleStatus(t *testing.T) {
	stray := FileInfo{Path: "/b/S1.notes.tmp", Present: true}
	cases := []struct {
		name  string
		state SampleBamState
		bqsr  bool
		want  SampleStatus
	}{
		{"short read finished", SampleBamState{Sample: "S1", RgmdCram: indexed(), BqsrCram: indexed()}, true, StatusComplete},
		{"short read leftover bam", SampleBamState{Sample: "S1", RgmdCram: indexed(), BqsrCram: indexed(), RgmdBam: usable()}, true, StatusCleanup},
		{"short read stray file", SampleBamState{Sample: "S1", RgmdCram: indexed(), BqsrCram: indexed(), OtherFiles: []FileInfo{stray}}, true, StatusCleanup},
		{"short read needs bqsr", SampleBamState{Sample: "S1", RgmdCram: indexed()}, true, StatusResume},
		{"short read without bqsr", SampleBamState{Sample: "S1", RgmdCram: indexed()}, false, StatusCleanup},
		{"short read cram missing index", SampleBamState{Sample: "S1", RgmdCram: usable(), BqsrCram: indexed()}, true, StatusResume},
		{"short read nothing", SampleBamState{Sample: "S1"}, true, StatusAlign},
		{"long read finished", SampleBamState{Sample: "MENINA_LR", RgmdCram: indexed()}, true, StatusComplete},
		{"long read with external kept", SampleBamState{Sample: "MENINA_LR", RgmdCram: indexed(),
			OtherFiles: []FileInfo{file("/b/movie.aligned.cram", usable())}}, true, StatusComplete},
		{"long read adopt", SampleBamState{Sample: "MENINA_LR",
			OtherFiles: []FileInfo{file("/b/movie.aligned.cram", usable())}}, true, StatusResume},
		{"long read two externals", SampleBamState{Sample: "MENINA_LR", OtherFiles: []FileInfo{
			file("/b/movie1.aligned.cram", usable()), file("/b/movie2.aligned.cram", usable())}}, true, StatusAmbiguous},
		{"long read nothing", SampleBamState{Sample: "MENINA_LR"}, true, StatusAlign},
	}
	for _, c := range cases {
		r := SampleReport{SampleAssessment: assessSample(c.state, c.bqsr, false), State: c.state}
		if got := r.Status(); got != c.want {
			t.Errorf("%s: status %s, want %s (plan %v)", c.name, got, c.want, r.Plan)
		}
	}
}

// The scan must plan exactly what a run would: same function, same inputs.
func TestAssessSampleMatchesBuildPlan(t *testing.T) {
	state := SampleBamState{Sample: "S1", SortedBam: usable()}
	a := assessSample(state, true, true)
	want, wantReads := buildPlan(state, planInputs{bqsr: true, gatkNeedsBam: true})
	if strings.Join(reasonStrings(a.Plan), ",") != strings.Join(reasonStrings(want), ",") || a.NeedsReads != wantReads {
		t.Errorf("plan %v / %v, want %v / %v", a.Plan, a.NeedsReads, want, wantReads)
	}
}

func reasonStrings(p []ProcessingReason) []string {
	out := make([]string, len(p))
	for i, r := range p {
		out[i] = string(r)
	}
	return out
}

func TestWriteScanReport(t *testing.T) {
	stray := FileInfo{Path: "/b/S2.notes.tmp", Present: true, Size: 2048}
	samples := []SampleBamState{
		{Sample: "S1", RgmdCram: file("/b/S1.RGMD.cram", indexed()), BqsrCram: file("/b/S1.RGMD_bqsr.cram", indexed())},
		{Sample: "S2", RgmdCram: file("/b/S2.RGMD.cram", indexed()), BqsrCram: file("/b/S2.RGMD_bqsr.cram", indexed()), OtherFiles: []FileInfo{stray}},
		{Sample: "S3", RgmdCram: file("/b/S3.RGMD.cram", usable())},
		{Sample: "S4"},
	}
	report := ScanReport{Species: "squash", RefVer: "v1.1", BQSR: true}
	for _, s := range samples {
		report.Samples = append(report.Samples, SampleReport{
			SampleAssessment: assessSample(s, true, false), Project: "P1", State: s})
	}

	var buf bytes.Buffer
	WriteScanReport(&buf, report)
	out := buf.String()
	if testing.Verbose() {
		t.Log("\n" + out)
	}
	for _, want := range []string{
		"P1/S1", "COMPLETE",
		"S2.notes.tmp", "unrecognised files block completion", "Blocked by unrecognised files: P1/S2",
		"S3.RGMD.cram", "ok, NO INDEX", "next: index RGMD.cram → run BQSR",
		"ALIGN", "needs the sample's reads",
		"Samples: 4",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("report missing %q:\n%s", want, out)
		}
	}
}

func TestHumanSize(t *testing.T) {
	for in, want := range map[int64]string{0: "0 B", 1023: "1023 B", 2048: "2.0 KiB", 5 << 30: "5.0 GiB"} {
		if got := humanSize(in); got != want {
			t.Errorf("humanSize(%d) = %q, want %q", in, got, want)
		}
	}
}
