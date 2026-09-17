package warehouse

import (
	"bufio"
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
)

func TestScanSyntheticEstate(t *testing.T) {
	tempDir := t.TempDir()

	// 1. Setup genomes directory with a minimal .dict file
	genomesDir := filepath.Join(tempDir, "genomes")
	assemblyDir := filepath.Join(genomesDir, "pepper", "UCD10Xv1.1", "assembly")
	if err := os.MkdirAll(assemblyDir, 0o755); err != nil {
		t.Fatal(err)
	}
	dictContent := "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:270000000\n"
	if err := os.WriteFile(filepath.Join(assemblyDir, "UCD10Xv1.1.dict"), []byte(dictContent), 0o644); err != nil {
		t.Fatal(err)
	}

	// 2. Setup DATA directory structure
	dataDir := filepath.Join(tempDir, "DATA")
	sampleDir := filepath.Join(dataDir, "pepper", "2025", "HP2831")
	cleanReadsDir := filepath.Join(sampleDir, DirCleanReads)
	refDir := filepath.Join(sampleDir, DirReferenceGenomes, "UCD10Xv1.1")
	bamsDir := filepath.Join(refDir, DirBams)
	gatkGvcfsDir := filepath.Join(refDir, "gatk_gvcfs")
	dvGvcfsDir := filepath.Join(refDir, "dv_gvcfs")
	legacyGvcfsDir := filepath.Join(refDir, DirLegacyGvcfs)
	mergedDir := filepath.Join(dataDir, "pepper", DirMergedVcfs, "UCD10Xv1.1", "gatk_gatk")

	dirs := []string{
		cleanReadsDir, bamsDir, gatkGvcfsDir, dvGvcfsDir, legacyGvcfsDir, mergedDir,
	}
	for _, d := range dirs {
		if err := os.MkdirAll(d, 0o755); err != nil {
			t.Fatal(err)
		}
	}

	// Helper to write dummy non-empty files
	writeFile := func(path string, content string) {
		if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
	}

	// Stray file at sample position
	writeFile(filepath.Join(dataDir, "pepper", "2025", "PEPPER_NAMES.xlsx"), "stray excel")

	// Clean reads + stray file in clean reads
	writeFile(filepath.Join(cleanReadsDir, "ForwardReads.fq.gz"), "fastq1")
	writeFile(filepath.Join(cleanReadsDir, "ReverseReads.fq.gz"), "fastq2")
	writeFile(filepath.Join(cleanReadsDir, "notes.txt"), "stray notes")

	// BAMs and indices
	writeFile(filepath.Join(bamsDir, "HP2831.RGMD.cram"), "cramcontent")
	writeFile(filepath.Join(bamsDir, "HP2831.RGMD.cram.crai"), "craicontent")

	// GATK gVCF
	writeFile(filepath.Join(gatkGvcfsDir, "HP2831.RGMD_bqsr.chr1.g.vcf.gz"), "gvcf1")
	writeFile(filepath.Join(gatkGvcfsDir, "HP2831.RGMD_bqsr.chr1.g.vcf.gz.tbi"), "tbi1")

	// DeepVariant gVCF + sidecars
	writeFile(filepath.Join(dvGvcfsDir, "HP2831.RGMD.chr1.g.vcf.gz"), "dv_gvcf")
	writeFile(filepath.Join(dvGvcfsDir, "HP2831.RGMD.chr1.g.vcf.gz.tbi"), "dv_tbi")
	writeFile(filepath.Join(dvGvcfsDir, "HP2831.RGMD.chr1.vcf.gz"), "dv_sidecar_vcf")
	writeFile(filepath.Join(dvGvcfsDir, "visual_report.html"), "<html>report</html>")

	// Legacy gVCF
	writeFile(filepath.Join(legacyGvcfsDir, "HP2831.chr1.g.vcf.gz"), "legacy_gvcf")

	// Joint VCFs
	writeFile(filepath.Join(mergedDir, "PEPPER.chr1.joint.vcf.gz"), "joint_vcf")
	writeFile(filepath.Join(mergedDir, "PEPPER.chr1.joint.vcf.gz.tbi"), "joint_tbi")
	writeFile(filepath.Join(mergedDir, "PEPPER.all.hard_filtered.vcf.gz"), "hard_filtered")
	writeFile(filepath.Join(mergedDir, "PEPPER.all.hard_filtered.vcf.gz.tbi"), "hard_filtered_tbi")

	// 3. Run Scan
	reportPath := filepath.Join(tempDir, "output_report.ndjson")
	opts := Options{
		Roots:      []string{dataDir},
		GenomesDir: genomesDir,
		ReportPath: reportPath,
		Validate:   ValidateNone,
		Workers:    2,
	}

	if err := Scan(opts); err != nil {
		t.Fatalf("Scan failed: %v", err)
	}

	// 4. Verify generated NDJSON
	f, err := os.Open(reportPath)
	if err != nil {
		t.Fatalf("Open report: %v", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	counts := map[string]int{}
	var lastLine string

	for scanner.Scan() {
		lastLine = scanner.Text()
		var row struct {
			Kind string `json:"kind"`
		}
		if err := json.Unmarshal([]byte(lastLine), &row); err != nil {
			t.Fatalf("unmarshal error on line: %s, err: %v", lastLine, err)
		}
		counts[row.Kind]++
	}

	if err := scanner.Err(); err != nil {
		t.Fatal(err)
	}

	// Ensure header was emitted
	if counts[KindReport] != 1 {
		t.Errorf("expected 1 report header, got %d", counts[KindReport])
	}
	// Volume, crop, reference, sample, sample_reference, joint_set
	if counts[KindVolume] != 1 {
		t.Errorf("expected 1 volume row, got %d", counts[KindVolume])
	}
	if counts[KindCrop] != 1 {
		t.Errorf("expected 1 crop row, got %d", counts[KindCrop])
	}
	if counts[KindReference] != 1 {
		t.Errorf("expected 1 reference row, got %d", counts[KindReference])
	}
	if counts[KindSample] != 1 {
		t.Errorf("expected 1 sample row, got %d", counts[KindSample])
	}
	if counts[KindSampleReference] != 1 {
		t.Errorf("expected 1 sample_reference row, got %d", counts[KindSampleReference])
	}
	if counts[KindJointSet] != 1 {
		t.Errorf("expected 1 joint_set row, got %d", counts[KindJointSet])
	}
	// Nonconformance rows: PEPPER_NAMES.xlsx, notes.txt, legacy gvcfs
	if counts[KindNonconformance] < 3 {
		t.Errorf("expected at least 3 nonconformance rows, got %d", counts[KindNonconformance])
	}

	// Summary must be last and complete
	var summary Summary
	if err := json.Unmarshal([]byte(lastLine), &summary); err != nil {
		t.Fatalf("unmarshal summary: %v", err)
	}
	if summary.Kind != KindSummary {
		t.Fatalf("last line kind = %q, want %q", summary.Kind, KindSummary)
	}
	if !summary.Complete {
		t.Errorf("summary Complete = false, want true")
	}
	if summary.FilesInventoried <= 0 || summary.BytesInventoried <= 0 {
		t.Errorf("files/bytes not counted: files=%d, bytes=%d", summary.FilesInventoried, summary.BytesInventoried)
	}
}

func TestReportPathInsideDataDirRefused(t *testing.T) {
	tempDir := t.TempDir()
	dataDir := filepath.Join(tempDir, "DATA")
	if err := os.MkdirAll(filepath.Join(dataDir, "pepper", "2025", "HP2831", DirCleanReads), 0o755); err != nil {
		t.Fatal(err)
	}

	reportPath := filepath.Join(dataDir, "report.ndjson")
	opts := Options{
		Roots:      []string{dataDir},
		ReportPath: reportPath,
	}

	if err := Scan(opts); err == nil {
		t.Fatalf("expected Scan to refuse report path inside data directory, but it succeeded")
	}
}

func TestDirectoryWithoutSampleMarkersIsJunk(t *testing.T) {
	tempDir := t.TempDir()
	dataDir := filepath.Join(tempDir, "DATA")
	junkDir := filepath.Join(dataDir, "beans", "2016", "reassembled_reads")
	if err := os.MkdirAll(junkDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(junkDir, "assembly.fa"), []byte("data"), 0o644); err != nil {
		t.Fatal(err)
	}
	reportPath := filepath.Join(tempDir, "report.ndjson")
	if err := Scan(Options{Roots: []string{dataDir}, ReportPath: reportPath, Workers: 1}); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(reportPath)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	counts := map[string]int{}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		var row struct {
			Kind string `json:"kind"`
		}
		if err := json.Unmarshal(scanner.Bytes(), &row); err != nil {
			t.Fatal(err)
		}
		counts[row.Kind]++
	}
	if counts[KindSample] != 0 || counts[KindJunkCandidate] != 1 {
		t.Fatalf("samples=%d junk=%d, want 0/1", counts[KindSample], counts[KindJunkCandidate])
	}
}
