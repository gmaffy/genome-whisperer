package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/variants"
)

// =============================================================================
// Types
// =============================================================================

// AlignmentStage represents how far along a sample's alignment pipeline is.
type AlignmentStage int

const (
	StageUnknown  AlignmentStage = iota
	StageScratch                 // No usable files at all — must align from raw reads
	StageSorted                  // sorted.bam exists; needs markdup
	StageRgmd                    // rgmd bam/cram exists; needs BQSR
	StageBqsr                    // bqsr bam/cram exists; needs conversion/indexing
	StageComplete                // bqsr.cram + rgmd.cram, both valid + indexed, no intermediates
)

// FileInfo holds presence, validity, and index metadata for a single alignment file.
type FileInfo struct {
	Path         string
	Size         int64
	Present      bool
	Valid        bool
	ValidateErr  error
	IndexPresent bool
	IndexSize    int64
}

func (f FileInfo) IsReady() bool {
	return f.Present && f.Valid
}

func (f FileInfo) IsIndexed() bool {
	return f.IndexPresent && f.IndexSize > 0
}

func (f FileInfo) IsComplete() bool {
	return f.IsReady() && f.IsIndexed()
}

// SampleAlignState describes every alignment file found for one sample.
type SampleAlignState struct {
	Sample     string
	BamDir     string
	CleanReads string // path to clean_reads dir, populated if scratch alignment needed

	SortedBam FileInfo
	RgmdBam   FileInfo
	RgmdCram  FileInfo
	BqsrBam   FileInfo
	BqsrCram  FileInfo

	OtherFiles []FileInfo
}

// Stage inspects the state and returns the earliest incomplete pipeline stage.
func (s SampleAlignState) Stage() AlignmentStage {
	switch {
	case s.BqsrCram.IsComplete() && s.RgmdCram.IsComplete() &&
		!s.SortedBam.Present && !s.RgmdBam.Present && !s.BqsrBam.Present:
		return StageComplete

	case s.BqsrCram.IsReady() || s.BqsrBam.IsReady():
		return StageBqsr

	case s.RgmdCram.IsReady() || s.RgmdBam.IsReady():
		return StageRgmd

	case s.SortedBam.IsReady():
		return StageSorted

	default:
		return StageScratch
	}
}

// =============================================================================
// Scanning
// =============================================================================

// ScanAlignments discovers all bam directories for a given species + refVer
// and concurrently inspects each sample's alignment files.
// Validation (variants.ValidateBam) runs per-file within each sample's goroutine.
func ScanAlignments(
	dataDir, species, refVer string,
	resolvedFasta string,
	verbose, quick bool,
) ([]SampleAlignState, error) {

	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		return nil, fmt.Errorf("resolving data dir: %w", err)
	}

	// Glob: <dataDir>/<species>/*/*/reference_genomes/<refVer>/bams
	pattern := filepath.Join(dataDirAbs, species, "*", "*", "reference_genomes", refVer, "bams")
	color.Yellow("\nScanning: %s\n", pattern)

	bamDirs, err := filepath.Glob(pattern)
	if err != nil {
		return nil, fmt.Errorf("glob pattern error: %w", err)
	}
	if len(bamDirs) == 0 {
		return nil, fmt.Errorf("no 'bams' directories found under %s", pattern)
	}
	color.Green("Found %d sample bam-director(ies).\n\n", len(bamDirs))

	// Scan all samples concurrently.
	type result struct {
		state SampleAlignState
		err   error
	}
	results := make([]result, len(bamDirs))
	var wg sync.WaitGroup

	for i, bamDir := range bamDirs {
		wg.Add(1)
		go func(idx int, dir string) {
			defer wg.Done()
			state, err := scanOneSampleDir(dir, dataDirAbs, species, refVer, resolvedFasta, verbose, quick)
			results[idx] = result{state, err}
		}(i, bamDir)
	}
	wg.Wait()

	var states []SampleAlignState
	for _, r := range results {
		if r.err != nil {
			color.Red("Scan error for a sample: %v\n", r.err)
			continue
		}
		states = append(states, r.state)
	}

	return states, nil
}

// scanOneSampleDir reads a single bam directory and returns the filled SampleAlignState.
// Directory layout: <dataDirAbs>/<species>/<group>/<sample>/reference_genomes/<refVer>/bams
func scanOneSampleDir(
	bamDir, dataDirAbs, species, refVer string,
	resolvedFasta string,
	verbose, quick bool,
) (SampleAlignState, error) {

	cleanDir := filepath.Clean(bamDir)

	// Walk up the path to extract sample name from the fixed-depth layout.
	// bams → refVer → reference_genomes → sample → group → species → dataDirAbs
	sampleName := filepath.Base(filepath.Dir(filepath.Dir(filepath.Dir(cleanDir))))

	// Derive the clean_reads path (used if scratch alignment is needed later).
	// <dataDirAbs>/<species>/<group>/<sample>/clean_reads
	sampleRoot := filepath.Dir(filepath.Dir(filepath.Dir(cleanDir)))
	cleanReadsDir := filepath.Join(sampleRoot, "clean_reads")

	state := SampleAlignState{
		Sample:     sampleName,
		BamDir:     cleanDir,
		CleanReads: cleanReadsDir,
	}

	entries, err := os.ReadDir(cleanDir)
	if err != nil {
		return state, fmt.Errorf("[%s] reading bam dir: %w", sampleName, err)
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue // bam dirs should be flat
		}

		name := entry.Name()
		fullPath := filepath.Join(cleanDir, name)

		// Skip index files — they are picked up alongside their parent.
		if strings.HasSuffix(name, ".bai") || strings.HasSuffix(name, ".crai") {
			continue
		}

		size := int64(0)
		if info, err := entry.Info(); err == nil {
			size = info.Size()
		}

		switch {
		case hasSuffix(name, "sorted.bam"):
			state.SortedBam = inspectBam(fullPath, size, resolvedFasta, verbose, quick)

		case hasSuffix(name, "rgmd.bam", "RGMD.bam"):
			state.RgmdBam = inspectBam(fullPath, size, resolvedFasta, verbose, quick)

		case hasSuffix(name, "rgmd.cram", "RGMD.cram"):
			state.RgmdCram = inspectCram(fullPath, size, resolvedFasta, verbose, quick)

		case hasSuffix(name, "bqsr.bam", "BQSR.bam"):
			state.BqsrBam = inspectBam(fullPath, size, resolvedFasta, verbose, quick)

		case hasSuffix(name, "bqsr.cram", "BQSR.cram"):
			state.BqsrCram = inspectCram(fullPath, size, resolvedFasta, verbose, quick)

		default:
			if !strings.HasSuffix(name, ".pdf") {
				color.Blue("[%s] Unexpected file: %s\n", sampleName, name)
				state.OtherFiles = append(state.OtherFiles, FileInfo{
					Path: fullPath, Size: size, Present: true,
				})
			}
		}
	}

	printSampleScanSummary(state)
	return state, nil
}

// inspectBam validates a BAM file and looks for its .bai index.
func inspectBam(path string, size int64, refFasta string, verbose, quick bool) FileInfo {
	f := FileInfo{Path: path, Size: size, Present: true}
	f.ValidateErr = variants.ValidateBam(path, refFasta, verbose, quick)
	f.Valid = f.ValidateErr == nil

	baiPath := strings.TrimSuffix(path, ".bam") + ".bai"
	if fi, err := os.Stat(baiPath); err == nil && fi.Mode().IsRegular() {
		f.IndexPresent = true
		f.IndexSize = fi.Size()
	}
	return f
}

// inspectCram validates a CRAM file and looks for its .crai index
// (both the <name>.crai and <name>.cram.crai conventions).
func inspectCram(path string, size int64, refFasta string, verbose, quick bool) FileInfo {
	f := FileInfo{Path: path, Size: size, Present: true}
	f.ValidateErr = variants.ValidateBam(path, refFasta, verbose, quick)
	f.Valid = f.ValidateErr == nil

	for _, craiPath := range []string{
		strings.TrimSuffix(path, ".cram") + ".crai",
		path + ".crai",
	} {
		if fi, err := os.Stat(craiPath); err == nil && fi.Mode().IsRegular() {
			f.IndexPresent = true
			f.IndexSize = fi.Size()
			break
		}
	}
	return f
}

// =============================================================================
// Action planning
// =============================================================================

// PipelineAction describes what needs to happen next for a sample.
type PipelineAction struct {
	Sample string
	Stage  AlignmentStage
	Steps  []Step
}

// Step is a single executable unit of work in the pipeline.
type Step struct {
	Program string
	Input   string // primary input file or sample identifier
}

// PlanActions inspects the scan state and returns the ordered steps
// needed to bring that sample to StageComplete.
func PlanActions(s SampleAlignState) PipelineAction {
	action := PipelineAction{Sample: s.Sample, Stage: s.Stage()}

	switch s.Stage() {

	case StageComplete:
		// Nothing to do.
		return action

	case StageScratch:
		// No usable alignment files at all — restart from raw reads.
		action.Steps = append(action.Steps, Step{Program: "align", Input: s.CleanReads})
		// The rest of the pipeline follows naturally after alignment.
		action.Steps = append(action.Steps, Step{Program: "markdup", Input: ""}) // input determined at runtime
		action.Steps = append(action.Steps, Step{Program: "bam_to_cram", Input: ""})
		action.Steps = append(action.Steps, Step{Program: "index_cram", Input: ""})
		action.Steps = append(action.Steps, Step{Program: "bqsr", Input: ""})
		action.Steps = append(action.Steps, Step{Program: "bam_to_cram", Input: ""})
		action.Steps = append(action.Steps, Step{Program: "index_cram", Input: ""})

	case StageSorted:
		// Have sorted.bam — run markdup then the CRAM/BQSR chain.
		action.Steps = append(action.Steps, Step{Program: "markdup", Input: s.SortedBam.Path})
		action.Steps = append(action.Steps,
			Step{Program: "bam_to_cram", Input: ""},
			Step{Program: "index_cram", Input: ""},
			Step{Program: "bqsr", Input: ""},
			Step{Program: "bam_to_cram", Input: ""},
			Step{Program: "index_cram", Input: ""},
		)

	case StageRgmd:
		// Have rgmd bam or cram — produce rgmd.cram if missing, then run BQSR.
		rgmdInput := chooseBestRgmd(s)
		if !s.RgmdCram.IsReady() {
			action.Steps = append(action.Steps, Step{Program: "bam_to_cram", Input: rgmdInput})
		}
		if !s.RgmdCram.IsIndexed() {
			action.Steps = append(action.Steps, Step{Program: "index_cram", Input: s.RgmdCram.Path})
		}
		if !s.BqsrCram.IsComplete() {
			action.Steps = append(action.Steps, Step{Program: "bqsr", Input: rgmdInput})
			action.Steps = append(action.Steps, Step{Program: "bam_to_cram", Input: ""})
			action.Steps = append(action.Steps, Step{Program: "index_cram", Input: ""})
		}

	case StageBqsr:
		// Have bqsr bam or cram — just ensure we have a valid indexed bqsr.cram.
		if !s.BqsrCram.IsReady() {
			action.Steps = append(action.Steps, Step{Program: "bam_to_cram", Input: s.BqsrBam.Path})
		}
		if !s.BqsrCram.IsIndexed() {
			action.Steps = append(action.Steps, Step{Program: "index_cram", Input: s.BqsrCram.Path})
		}
	}

	// Cleanup: once both CRAMs are complete, intermediate BAMs can be removed.
	// We schedule this regardless of what steps precede it; the executor will
	// only run cleanup after the pipeline steps succeed.
	action.Steps = append(action.Steps, intermediateCleanupSteps(s)...)

	return action
}

// chooseBestRgmd returns the path of the best available rgmd file.
// Prefer the BAM (needed for BQSR tools that can't read CRAM), fall back to CRAM.
func chooseBestRgmd(s SampleAlignState) string {
	if s.RgmdBam.IsReady() {
		return s.RgmdBam.Path
	}
	return s.RgmdCram.Path
}

// intermediateCleanupSteps returns rm Steps for all intermediate BAM files
// that can be deleted once the final CRAMs are in place.
func intermediateCleanupSteps(s SampleAlignState) []Step {
	var steps []Step
	for _, f := range []FileInfo{s.SortedBam, s.RgmdBam, s.BqsrBam} {
		if f.Present {
			steps = append(steps, Step{Program: "rm", Input: f.Path})
		}
	}
	for _, f := range s.OtherFiles {
		steps = append(steps, Step{Program: "rm", Input: f.Path})
	}
	return steps
}

// =============================================================================
// Helpers
// =============================================================================

// hasSuffix returns true if name ends with any of the given suffixes.
func hasSuffix(name string, suffixes ...string) bool {
	for _, s := range suffixes {
		if strings.HasSuffix(name, s) {
			return true
		}
	}
	return false
}

func printSampleScanSummary(s SampleAlignState) {
	label := func(f FileInfo, name string) string {
		if !f.Present {
			return color.HiBlackString("  %-14s absent", name)
		}
		validity := color.GreenString("valid")
		if !f.Valid {
			validity = color.RedString("INVALID")
		}
		indexed := color.GreenString("indexed")
		if !f.IsIndexed() {
			indexed = color.YellowString("no-index")
		}
		return fmt.Sprintf("  %-14s %s  %s", name, validity, indexed)
	}

	stageColors := map[AlignmentStage]func(string, ...interface{}) string{
		StageComplete: color.GreenString,
		StageBqsr:     color.YellowString,
		StageRgmd:     color.YellowString,
		StageSorted:   color.YellowString,
		StageScratch:  color.RedString,
	}
	stageNames := map[AlignmentStage]string{
		StageComplete: "COMPLETE",
		StageBqsr:     "needs bqsr→cram",
		StageRgmd:     "needs bqsr",
		StageSorted:   "needs markdup",
		StageScratch:  "needs full alignment",
	}

	stage := s.Stage()
	colorFn := stageColors[stage]
	fmt.Printf("\n[%s]  stage: %s\n", color.CyanString(s.Sample), colorFn(stageNames[stage]))
	fmt.Println(label(s.SortedBam, "sorted.bam"))
	fmt.Println(label(s.RgmdBam, "rgmd.bam"))
	fmt.Println(label(s.RgmdCram, "rgmd.cram"))
	fmt.Println(label(s.BqsrBam, "bqsr.bam"))
	fmt.Println(label(s.BqsrCram, "bqsr.cram"))
	if len(s.OtherFiles) > 0 {
		color.Blue("  %d unexpected file(s) will be removed\n", len(s.OtherFiles))
	}
}
