package alignment

import (
	"strings"
	"testing"
)

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

func TestMarkDupCmdUsesPbmarkdupForPbmm2(t *testing.T) {
	// pbmm2 has no "markdup" subcommand — pbmarkdup is a separate binary, which
	// is the one cmd/AlignReads.go adds to the dependency list for this aligner.
	// This branch shipped for a while as "pbmm2 markdup", which no other kind of
	// test can see: the tool name lives inside a format string.
	cmd, rgmdBam := markDupCmd("/ref/genome.fa", "/data/bams/MO971_LR.sorted.bam", DupMarkerPbmarkdup, "INFO", "-Xmx8G", 4)

	if !strings.HasPrefix(cmd, "pbmarkdup ") {
		t.Errorf("pbmm2 markdup command = %q, want it to start with %q", cmd, "pbmarkdup ")
	}
	if strings.Contains(cmd, "pbmm2 markdup") {
		t.Errorf("command still shells the non-existent pbmm2 subcommand: %q", cmd)
	}
	if strings.Contains(cmd, "MarkDuplicates") {
		t.Errorf("pbmm2 path must not reach GATK MarkDuplicates: %q", cmd)
	}
	if want := "/data/bams/MO971_LR.RGMD.bam"; rgmdBam != want {
		t.Errorf("rgmdBam = %q, want %q", rgmdBam, want)
	}
	// The output must be the last argument: pbmarkdup's usage is
	// "pbmarkdup [options] <INFILE> <OUTFILE>".
	if !strings.HasSuffix(cmd, rgmdBam) {
		t.Errorf("command = %q, want it to end with the output %q", cmd, rgmdBam)
	}
}

func TestMarkDupCmdUsesGatkByDefault(t *testing.T) {
	// Anything that is not an explicit request for pbmarkdup gets GATK, which
	// is what every sample the directory pipeline produces uses — long-read
	// ones included, because pbmarkdup cannot write a BAM aligned from a FASTQ.
	for _, aligner := range []string{DupMarkerGatk, "bwa-mem2", "pbmm2", ""} {
		cmd, rgmdBam := markDupCmd("/ref/genome.fa", "/data/bams/ON1708.sorted.bam", aligner, "INFO", "-Xmx8G", 4)

		if !strings.Contains(cmd, "MarkDuplicates") {
			t.Errorf("%s: command = %q, want GATK MarkDuplicates", aligner, cmd)
		}
		// A whole-genome duplicate marking spills roughly the size of the read
		// set, so losing --TMP_DIR fills /tmp rather than the spill disk.
		if !strings.Contains(cmd, "--TMP_DIR ") {
			t.Errorf("%s: command = %q, want --TMP_DIR", aligner, cmd)
		}
		if !strings.Contains(cmd, "-M /data/bams/ON1708.RGMD.metrics.txt") {
			t.Errorf("%s: command = %q, want the metrics file beside the rgmd bam", aligner, cmd)
		}
		if want := "/data/bams/ON1708.RGMD.bam"; rgmdBam != want {
			t.Errorf("%s: rgmdBam = %q, want %q", aligner, rgmdBam, want)
		}
	}
}

func TestPbmm2AlignCmdStampsTheReadGroup(t *testing.T) {
	cmd := pbmm2AlignCmd(
		"/data/MO971_LR/clean_reads/m64011_200521_073741.fastq.gz",
		"/ref/genome.fa.HIFI.mmi",
		"/data/bams/MO971_LR.sorted.partial.bam",
		"MO971_LR", "MO971_LR_1", "HIFI", 8)

	// HaplotypeCaller names the gVCF's sample column from SM, so a missing or
	// wrong read group is a sample-identity bug, not a cosmetic one. pbmm2
	// documents --rg for FASTA/Q input, which is the only input this takes.
	for _, want := range []string{
		`--rg '@RG\tID:MO971_LR.1\tSM:MO971_LR\tLB:MO971_LR_1\tPL:PACBIO'`,
		"--preset HIFI",
		"/ref/genome.fa.HIFI.mmi",
		"m64011_200521_073741.fastq.gz",
		"-j 8",
	} {
		if !strings.Contains(cmd, want) {
			t.Errorf("pbmm2 command = %q,\nwant it to contain %q", cmd, want)
		}
	}

	// Piped into the repo's own sort rather than pbmm2 --sort: pbmm2 would
	// spill to a directory this pipeline does not control, and the GATK steps
	// running alongside depend on that space.
	if !strings.Contains(cmd, "| samtools sort") {
		t.Errorf("pbmm2 command = %q, want it piped into samtools sort", cmd)
	}
	if strings.Contains(cmd, "--sort ") {
		t.Errorf("pbmm2 command = %q, want pbmm2's own --sort not used", cmd)
	}
	if !strings.Contains(cmd, "-T ") {
		t.Errorf("pbmm2 command = %q, want the sort given an explicit -T spill prefix", cmd)
	}

	// The index is built once per run before any sample starts; building it
	// here is what made concurrent long-read samples race to write one file.
	if strings.Contains(cmd, "pbmm2 index") {
		t.Errorf("pbmm2 command = %q, want no lazy index build", cmd)
	}
}

func TestPbmm2AlignCmdNormalisesThePreset(t *testing.T) {
	// pbmm2 only accepts the upper-case spellings, and --preset comes straight
	// from a command line.
	cmd := pbmm2AlignCmd("/r.fq", "/ref.mmi", "/out.bam", "S1", "S1_1", "hifi", 4)
	if !strings.Contains(cmd, "--preset HIFI") {
		t.Errorf("pbmm2 command = %q, want the preset upper-cased", cmd)
	}
}
