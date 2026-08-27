package alignment

import "testing"

func TestRgmdBamPath(t *testing.T) {
	tests := []struct {
		sortedBam string
		want      string
	}{
		// The ".sorted" tag is dropped so the durable artefacts keep the
		// <sample>.RGMD.* names the data directories are built around.
		{"/data/ON1708/bams/ON1708.sorted.bam", "/data/ON1708/bams/ON1708.RGMD.bam"},
		{"/data/ON1708/bams/ON1708.bam", "/data/ON1708/bams/ON1708.RGMD.bam"},
		{"/data/ON1708/bams/ON1708.sorted.cram", "/data/ON1708/bams/ON1708.RGMD.bam"},
	}
	for _, tc := range tests {
		if got := RgmdBamPath(tc.sortedBam); got != tc.want {
			t.Errorf("RgmdBamPath(%q) = %q, want %q", tc.sortedBam, got, tc.want)
		}
	}
}

func TestIndexPathPicksCSIForBam(t *testing.T) {
	// A BAI cannot address the contigs this pipeline works on; CRAM carries its
	// own .crai. See the note on IndexPath.
	if got := IndexPath("/x/ON1708.RGMD.bam"); got != "/x/ON1708.RGMD.bam.csi" {
		t.Errorf("bam index = %q, want .csi", got)
	}
	if got := IndexPath("/x/ON1708.RGMD.cram"); got != "/x/ON1708.RGMD.cram.crai" {
		t.Errorf("cram index = %q, want .crai", got)
	}
}
