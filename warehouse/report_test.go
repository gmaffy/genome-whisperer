package warehouse

import (
	"bufio"
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
)

func TestReportWriterRoundtrip(t *testing.T) {
	tmpDir := t.TempDir()
	reportPath := filepath.Join(tmpDir, "report.ndjson")

	w, err := NewWriter(reportPath)
	if err != nil {
		t.Fatalf("NewWriter failed: %v", err)
	}

	header := Header{
		Kind:            KindReport,
		SchemaVersion:   SchemaVersion,
		StandardVersion: StandardVersion,
		Tool:            "genome-whisperer",
		ToolVersion:     "1.0.0",
		Command:         "ScanWarehouse",
		RunID:           "test_run_123",
		StartedAt:       "2026-09-08T15:00:00Z",
		Host:            "nas",
		OS:              "linux",
		RootsRequested:  []string{"/mnt/u/DATA"},
		GenomesDir:      "/mnt/z/genomes",
		Validate:        ValidateNone,
		ScanWorkers:     16,
		UnitRule:        "top21-by-length+MT+Pltd+contigs",
	}
	if err := w.Write(KindReport, header); err != nil {
		t.Fatalf("Write(KindReport) failed: %v", err)
	}

	vol := Volume{
		Kind:           KindVolume,
		Volume:         "/mnt/u",
		FSDev:          72,
		Roots:          []string{"/mnt/u/DATA", "/mnt/u/data"},
		DataDir:        "/mnt/u/DATA",
		Spelling:       "DATA",
		Present:        true,
		CapacityKnown:  true,
		BytesTotal:     7585801011200,
		BytesFree:      471443374080,
		BytesAvailable: 471443374080,
		BlockSize:      4096,
		Crops:          []string{"pepper"},
	}
	if err := w.Write(KindVolume, vol); err != nil {
		t.Fatalf("Write(KindVolume) failed: %v", err)
	}

	crop := Crop{
		Kind:         KindCrop,
		Volume:       "/mnt/u",
		Crop:         "pepper",
		Spelling:     "pepper",
		Path:         "/mnt/u/DATA/pepper",
		Years:        []string{"2025"},
		SampleCount:  1,
		References:   []string{"UCD10Xv1.1"},
		JointLayouts: []string{JointStandard},
		HasMerged:    true,
	}
	if err := w.Write(KindCrop, crop); err != nil {
		t.Fatalf("Write(KindCrop) failed: %v", err)
	}

	ref := Reference{
		Kind:              KindReference,
		Volume:            "/mnt/u",
		Crop:              "pepper",
		Reference:         "UCD10Xv1.1",
		Spelling:          "UCD10Xv1.1",
		DictStatus:        DictResolved,
		DictPath:          "/mnt/z/genomes/pepper/UCD10Xv1.1/assembly/test.dict",
		GenomesRef:        "UCD10Xv1.1",
		Rule:              "top21-by-length+MT+Pltd+contigs",
		SeqCount:          81202,
		ExpectedUnitCount: 24,
		Units: []Unit{
			{ID: "chr1", Label: "chr1", UnitKind: "chrom", Length: 270000000, MemberCount: 1},
			{ID: "contigs", Label: "contigs", UnitKind: "contigs", Length: 150000000, MemberCount: 81179},
		},
		SampleReferenceCount: 1,
	}
	if err := w.Write(KindReference, ref); err != nil {
		t.Fatalf("Write(KindReference) failed: %v", err)
	}

	batch := ReferenceBatch{
		Kind:        KindReferenceBatch,
		Crop:        "pepper",
		Reference:   "UCD10Xv1.1",
		MemberCount: 2,
		TotalLength: 150000000,
		SeqIDs:      []string{"ctg1", "ctg2"},
	}
	if err := w.Write(KindReferenceBatch, batch); err != nil {
		t.Fatalf("Write(KindReferenceBatch) failed: %v", err)
	}

	sample := Sample{
		Kind:                KindSample,
		Volume:              "/mnt/u",
		Crop:                "pepper",
		Year:                "2025",
		Sample:              "HP2831",
		Path:                "/mnt/u/DATA/pepper/2025/HP2831",
		LongRead:            false,
		HasCleanReads:       true,
		HasReferenceGenomes: true,
		ReadLayout:          LayoutPaired,
		ReadsFwd:            1,
		ReadsRev:            1,
		ReadsSE:             0,
		ReadsBytes:          2000000,
		ReadFiles: []ReadFile{
			{Name: "ForwardReads.fq.gz", Role: "fwd", Bytes: 1000000, MTime: "2025-06-01T12:00:00Z"},
			{Name: "ReverseReads.fq.gz", Role: "rev", Bytes: 1000000, MTime: "2025-06-01T12:00:00Z"},
		},
		References: []string{"UCD10Xv1.1"},
	}
	if err := w.Write(KindSample, sample); err != nil {
		t.Fatalf("Write(KindSample) failed: %v", err)
	}

	sampleRef := SampleReference{
		Kind:                KindSampleReference,
		Volume:              "/mnt/u",
		Crop:                "pepper",
		Year:                "2025",
		Sample:              "HP2831",
		Reference:           "UCD10Xv1.1",
		Spelling:            "UCD10Xv1.1",
		Path:                "/mnt/u/DATA/pepper/2025/HP2831/reference_genomes/UCD10Xv1.1",
		AlignmentDirPresent: true,
		Alignments: []AlignmentFile{
			{
				Name:       "HP2831.RGMD.cram",
				Role:       "rgmd_cram",
				Bytes:      2000000,
				MTime:      "2025-06-02T10:00:00Z",
				Index:      "HP2831.RGMD.cram.crai",
				IndexBytes: 1000,
				Validated:  false,
			},
		},
		AlignmentBytes:  2000000,
		AlignmentRoles:  []string{"rgmd_cram"},
		InventoryStatus: InventoryStandard,
		GvcfSets: []GvcfSet{
			{
				Dir:               "gatk_gvcfs",
				Caller:            "gatk",
				ExpectedUnitCount: 24,
				PresentCount:      24,
				IndexedCount:      24,
				Bytes:             5000000,
				Stems:             []string{"HP2831.RGMD_bqsr"},
				LabelSpellings:    []string{SpellingSanitised},
			},
		},
	}
	if err := w.Write(KindSampleReference, sampleRef); err != nil {
		t.Fatalf("Write(KindSampleReference) failed: %v", err)
	}

	joint := JointSet{
		Kind:                 KindJointSet,
		Volume:               "/mnt/u",
		Crop:                 "pepper",
		Reference:            "UCD10Xv1.1",
		Tag:                  "gatk_gatk",
		Caller:               "gatk",
		Merger:               "gatk",
		Layout:               JointStandard,
		Dir:                  "/mnt/u/DATA/pepper/MERGED_VCFs/UCD10Xv1.1/gatk_gatk",
		Container:            ContainerVcfGz,
		ExpectedUnitCount:    24,
		PresentCount:         24,
		IndexedCount:         24,
		JointBytes:           10000000,
		HasAllVcf:            true,
		HasAllIndex:          true,
		HasHardFiltered:      true,
		HasHardFilteredIndex: true,
		HasSnpEffVcf:         true,
		HasSnpEffTsv:         true,
		HasEffTsv:            true,
		HasDescTsv:           false,
		HasSuperVcfTsv:       false,
	}
	if err := w.Write(KindJointSet, joint); err != nil {
		t.Fatalf("Write(KindJointSet) failed: %v", err)
	}

	nc := Nonconformance{
		Kind:     KindNonconformance,
		Rule:     RuleLegacyGvcfDir,
		Severity: SeverityWarning,
		Volume:   "/mnt/u",
		Crop:     "pepper",
		Year:     "2025",
		Sample:   "HP2831",
		Path:     "/mnt/u/DATA/pepper/2025/HP2831/reference_genomes/CM334/gvcfs",
		Found:    "gvcfs",
		Expected: "gatk_gvcfs or dv_gvcfs",
		Count:    24,
		Detail:   "legacy gvcfs dir",
	}
	if err := w.Write(KindNonconformance, nc); err != nil {
		t.Fatalf("Write(KindNonconformance) failed: %v", err)
	}

	w.CountFiles(50, 17000000)
	w.CountError()

	if err := w.Close(true); err != nil {
		t.Fatalf("Close failed: %v", err)
	}

	// Verify file was written and can be read line by line as valid JSON.
	f, err := os.Open(reportPath)
	if err != nil {
		t.Fatalf("Open report failed: %v", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var lineCount int
	kinds := map[string]int{}

	for scanner.Scan() {
		lineCount++
		var generic map[string]interface{}
		if err := json.Unmarshal(scanner.Bytes(), &generic); err != nil {
			t.Fatalf("line %d is not valid JSON: %v", lineCount, err)
		}
		kind, ok := generic["kind"].(string)
		if !ok || kind == "" {
			t.Fatalf("line %d missing 'kind'", lineCount)
		}
		kinds[kind]++
	}

	if err := scanner.Err(); err != nil {
		t.Fatalf("scanner error: %v", err)
	}

	// 9 written rows + 1 summary row = 10 rows
	if lineCount != 10 {
		t.Fatalf("expected 10 lines, got %d", lineCount)
	}
	if kinds[KindSummary] != 1 {
		t.Errorf("expected 1 summary row, got %d", kinds[KindSummary])
	}
	if kinds[KindReport] != 1 {
		t.Errorf("expected 1 report row, got %d", kinds[KindReport])
	}
}

func TestReportWriterAbandon(t *testing.T) {
	tmpDir := t.TempDir()
	reportPath := filepath.Join(tmpDir, "abandoned.ndjson")

	w, err := NewWriter(reportPath)
	if err != nil {
		t.Fatalf("NewWriter failed: %v", err)
	}

	if err := w.Write(KindReport, Header{Kind: KindReport}); err != nil {
		t.Fatalf("Write failed: %v", err)
	}

	w.Abandon()

	if _, err := os.Stat(reportPath); !os.IsNotExist(err) {
		t.Fatalf("expected destination file to not exist after Abandon, but stat err = %v", err)
	}

	// Also ensure no leftover tmp files in tmpDir
	entries, err := os.ReadDir(tmpDir)
	if err != nil {
		t.Fatalf("ReadDir failed: %v", err)
	}
	if len(entries) != 0 {
		t.Fatalf("expected 0 files in temp dir after Abandon, got %d", len(entries))
	}
}

func TestParseContractFixture(t *testing.T) {
	fixturePath := "/mnt/e/GitHub/plennegybfx/plennegybfx_refactor/plennegybfx/data_warehouse/tests/fixtures/warehouse_scan.ndjson"
	f, err := os.Open(fixturePath)
	if err != nil {
		t.Skipf("Contract fixture not found at %s: %v", fixturePath, err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	line := 0
	counts := map[string]int{}

	for scanner.Scan() {
		line++
		data := scanner.Bytes()
		var generic struct {
			Kind string `json:"kind"`
		}
		if err := json.Unmarshal(data, &generic); err != nil {
			t.Fatalf("line %d unmarshal kind: %v", line, err)
		}
		counts[generic.Kind]++

		switch generic.Kind {
		case KindReport:
			var h Header
			if err := json.Unmarshal(data, &h); err != nil {
				t.Fatalf("line %d unmarshal Header: %v", line, err)
			}
			if h.SchemaVersion != SchemaVersion {
				t.Errorf("header schema version %d != %d", h.SchemaVersion, SchemaVersion)
			}
		case KindVolume:
			var v Volume
			if err := json.Unmarshal(data, &v); err != nil {
				t.Fatalf("line %d unmarshal Volume: %v", line, err)
			}
		case KindCrop:
			var c Crop
			if err := json.Unmarshal(data, &c); err != nil {
				t.Fatalf("line %d unmarshal Crop: %v", line, err)
			}
		case KindReference:
			var r Reference
			if err := json.Unmarshal(data, &r); err != nil {
				t.Fatalf("line %d unmarshal Reference: %v", line, err)
			}
		case KindReferenceBatch:
			var rb ReferenceBatch
			if err := json.Unmarshal(data, &rb); err != nil {
				t.Fatalf("line %d unmarshal ReferenceBatch: %v", line, err)
			}
		case KindSample:
			var s Sample
			if err := json.Unmarshal(data, &s); err != nil {
				t.Fatalf("line %d unmarshal Sample: %v", line, err)
			}
		case KindSampleReference:
			var sr SampleReference
			if err := json.Unmarshal(data, &sr); err != nil {
				t.Fatalf("line %d unmarshal SampleReference: %v", line, err)
			}
		case KindJointSet:
			var js JointSet
			if err := json.Unmarshal(data, &js); err != nil {
				t.Fatalf("line %d unmarshal JointSet: %v", line, err)
			}
		case KindNonconformance:
			var nc Nonconformance
			if err := json.Unmarshal(data, &nc); err != nil {
				t.Fatalf("line %d unmarshal Nonconformance: %v", line, err)
			}
		case KindSummary:
			var s Summary
			if err := json.Unmarshal(data, &s); err != nil {
				t.Fatalf("line %d unmarshal Summary: %v", line, err)
			}
			if !s.Complete {
				t.Errorf("summary Complete = false")
			}
			for k, count := range counts {
				if s.RowCounts[k] != count {
					t.Errorf("summary count for %s = %d, counted %d", k, s.RowCounts[k], count)
				}
			}
		default:
			t.Fatalf("unexpected kind %q at line %d", generic.Kind, line)
		}
	}

	if err := scanner.Err(); err != nil {
		t.Fatalf("scanner error: %v", err)
	}

	if line != 26 {
		t.Errorf("expected 26 lines, read %d", line)
	}
}

