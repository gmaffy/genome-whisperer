package warehouse

import (
	"testing"
)

func TestIsYear(t *testing.T) {
	valid := []string{"1999", "2020", "2025", "2026", "2099"}
	for _, y := range valid {
		if !IsYear(y) {
			t.Errorf("IsYear(%q) = false, want true", y)
		}
	}

	invalid := []string{
		"1899", "2100", "202", "20251", "PEPPER_NAMES.xlsx",
		"MO834", "MERGED_VCFs", "VCFs", "HP2831", "",
	}
	for _, y := range invalid {
		if IsYear(y) {
			t.Errorf("IsYear(%q) = true, want false", y)
		}
	}
}

func TestGvcfDirCaller(t *testing.T) {
	tests := []struct {
		dir        string
		wantCaller string
		wantOk     bool
	}{
		{"gatk_gvcfs", "gatk", true},
		{"dv_gvcfs", "dv", true},
		{"gvcfs", "", false},
		{"VCFs", "", false},
		{"bams", "", false},
		{"clean_reads", "", false},
	}

	for _, tt := range tests {
		caller, ok := GvcfDirCaller(tt.dir)
		if ok != tt.wantOk || caller != tt.wantCaller {
			t.Errorf("GvcfDirCaller(%q) = (%q, %v), want (%q, %v)",
				tt.dir, caller, ok, tt.wantCaller, tt.wantOk)
		}
	}
}

func TestSplitCallerMergerTag(t *testing.T) {
	tests := []struct {
		tag        string
		wantCaller string
		wantMerger string
		wantOk     bool
	}{
		{"gatk_gatk", "gatk", "gatk", true},
		{"gatk_glnexus", "gatk", "glnexus", true},
		{"dv_glnexus", "dv", "glnexus", true},
		{"deepvariant_glnexus", "", "", false},
		{"gatk", "", "", false},
		{"random_tag", "", "", false},
	}

	for _, tt := range tests {
		caller, merger, ok := SplitCallerMergerTag(tt.tag)
		if ok != tt.wantOk || caller != tt.wantCaller || merger != tt.wantMerger {
			t.Errorf("SplitCallerMergerTag(%q) = (%q, %q, %v), want (%q, %q, %v)",
				tt.tag, caller, merger, ok, tt.wantCaller, tt.wantMerger, tt.wantOk)
		}
	}
}

func TestIsPruned(t *testing.T) {
	pruned := []string{
		"", ".git", ".snapshot", "@Recycle", "@Recently-Snapshot",
		"lost+found", "#recycle", "work", "shards", "QC",
		"tmp_dv_12345", "tmp_intermediate",
	}
	for _, name := range pruned {
		if !IsPruned(name) {
			t.Errorf("IsPruned(%q) = false, want true", name)
		}
	}

	notPruned := []string{
		"2025", "pepper", "clean_reads", "reference_genomes",
		"UCD10Xv1.1", "bams", "gatk_gvcfs", "dv_gvcfs", "MERGED_VCFs",
	}
	for _, name := range notPruned {
		if IsPruned(name) {
			t.Errorf("IsPruned(%q) = true, want false", name)
		}
	}
}

func TestIsFastq(t *testing.T) {
	fastq := []string{
		"ForwardReads.fq.gz", "ReverseReads.fastq.gz",
		"reads.fq", "reads.fastq", "READ1.FQ.GZ", "read_2.FASTQ",
	}
	for _, name := range fastq {
		if !IsFastq(name) {
			t.Errorf("IsFastq(%q) = false, want true", name)
		}
	}

	notFastq := []string{
		"sample.bam", "sample.cram", "sample.g.vcf.gz", "sample.vcf.gz",
		"notes.txt", "data.csv", "metrics.pdf",
	}
	for _, name := range notFastq {
		if IsFastq(name) {
			t.Errorf("IsFastq(%q) = true, want false", name)
		}
	}
}
