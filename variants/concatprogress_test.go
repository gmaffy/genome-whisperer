package variants

import "testing"

// The bar for the concatenation is driven entirely by this pattern, so the exact
// shape of a MergeVcfs progress line is load-bearing. These are real lines.
func TestMergeVcfsProgressParsing(t *testing.T) {
	tests := []struct {
		name   string
		line   string
		want   int64 // 0 with wantOK false means the line carries no count
		wantOK bool
	}{
		{
			name: "a real progress line",
			line: "INFO\t2026-09-03 10:35:54\tMergeVcfs\tProcessed     9,920,000 records.  Elapsed time: 00:04:52s.  Time for last 10,000:    0s.  Last read position: scaffold_17:256,545",
			want: 9_920_000, wantOK: true,
		},
		{
			name: "single digit, no separators",
			line: "INFO\tMergeVcfs\tProcessed     1 records.  Elapsed time: 00:00:01s.",
			want: 1, wantOK: true,
		},
		{
			// "Time for last 10,000" must not be mistaken for the running total.
			name: "the interval count is not the total",
			line: "INFO\tMergeVcfs\tProcessed        10,000 records.  Time for last 10,000:    0s.",
			want: 10_000, wantOK: true,
		},
		{
			name: "the completion line carries no count",
			line: "[Thu Sep 03 10:55:21 SAST 2026] picard.vcf.MergeVcfs done. Elapsed time: 4.38 minutes.",
			want: 0, wantOK: false,
		},
		{
			name: "an unrelated line is ignored",
			line: "Using GATK jar /home/godwin/tools/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar",
			want: 0, wantOK: false,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, ok := mergeVcfsRecordCount(tc.line)
			if got != tc.want || ok != tc.wantOK {
				t.Errorf("mergeVcfsRecordCount() = (%d, %v), want (%d, %v)", got, ok, tc.want, tc.wantOK)
			}
		})
	}
}

// The advance function must survive the lines that carry no count. It cannot be
// checked by watching the bar move, because a bar has no terminal to draw on
// under `go test` and every Add on an invisible bar is a no-op — which is why
// the decision it makes lives in mergeVcfsRecordCount, tested above.
func TestConcatProgressToleratesEveryLine(t *testing.T) {
	// No real VCFs, so totalVcfRecords fails and this is the spinning bar.
	bar, advance := concatProgress([]string{"/nonexistent/a.vcf.gz"})
	if bar == nil {
		t.Fatal("no bar returned")
	}
	for _, line := range []string{
		"INFO\tMergeVcfs\tProcessed     10,000 records.",
		"[Thu Sep 03 10:55:21 SAST 2026] picard.vcf.MergeVcfs done.",
		"Processed  records.",
		"",
	} {
		advance(line)
	}
}

// A missing or countless index must degrade to the spinning bar, never to a
// bar with a bogus total.
func TestTotalVcfRecordsReportsFailure(t *testing.T) {
	if _, ok := totalVcfRecords([]string{"/nonexistent/a.vcf.gz"}); ok {
		t.Error("totalVcfRecords claimed a total for a file that does not exist")
	}
	if _, ok := totalVcfRecords(nil); ok {
		t.Error("totalVcfRecords claimed a total for no inputs")
	}
}
