package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
)

// AlignConfig holds all parameters needed to run the alignment pipeline.
type AlignConfig struct {
	DataDir      string
	Species      string
	RefVer       string
	RefFasta     string // explicit fasta path; mutually exclusive with GenomesDir
	GenomesDir   string // auto-discovery dir; mutually exclusive with RefFasta
	Aligner      string // "bwa-mem2" or "pbmm2"
	GatkLogLevel string
	Threads      int
	Verbose      bool
	Quick        bool // quick validation (gzip integrity only)
	SkipVer      bool
	BQSR         bool
	Bootstrap    bool
	KnownSites   []string
}

// RunAlignReadsDir is the top-level entry point.
// It validates inputs, scans each sample's bam directory, plans what work
// remains, and executes the outstanding pipeline steps concurrently.
func RunAlignReadsDir(cfg AlignConfig) {
	color.Cyan("═══════════════════════════════════════════════════\n")
	color.Cyan("   genome-whisperer  ·  alignment pipeline\n")
	color.Cyan("═══════════════════════════════════════════════════\n\n")

	// ── 1. Validate config ────────────────────────────────────────────────── //
	if err := validateConfig(cfg); err != nil {
		color.Red("Configuration error: %v\n", err)
		return
	}

	// ── 2. Resolve reference fasta ────────────────────────────────────────── //
	resolvedFasta, err := resolveRefFasta(cfg.RefFasta, cfg.GenomesDir, cfg.Species, cfg.RefVer)
	if err != nil {
		color.Red("Cannot resolve reference fasta: %v\n", err)
		return
	}

	// ── 3. Scan all sample bam directories ────────────────────────────────── //
	color.Cyan("\n─── Scanning alignment directories ───\n")
	states, err := ScanAlignments(
		cfg.DataDir, cfg.Species, cfg.RefVer,
		resolvedFasta,
		cfg.Verbose, cfg.Quick,
	)
	if err != nil {
		color.Red("Scan failed: %v\n", err)
		return
	}

	// ── 4. Plan actions for each sample ───────────────────────────────────── //
	var pending []PipelineAction
	completeCount := 0

	for _, state := range states {
		action := PlanActions(state)
		if action.Stage == StageComplete {
			completeCount++
		} else {
			pending = append(pending, action)
		}
	}

	color.Green("\n%d/%d sample(s) already complete.\n", completeCount, len(states))
	if len(pending) == 0 {
		color.Green("Nothing to do — all samples are complete. ✅\n")
		return
	}
	color.Yellow("%d sample(s) need work — starting pipeline...\n\n", len(pending))

	// ── 5. Execute pipeline concurrently, bounded by thread count ─────────── //
	threads := cfg.Threads
	if threads <= 0 {
		threads = runtime.NumCPU()
	}

	sem := make(chan struct{}, threads)
	var wg sync.WaitGroup
	results := make([]sampleResult, len(pending))

	for i, action := range pending {
		wg.Add(1)
		go func(idx int, a PipelineAction) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()

			results[idx] = executeSample(a, cfg, resolvedFasta)
		}(i, action)
	}
	wg.Wait()

	// ── 6. Print final summary ────────────────────────────────────────────── //
	printFinalSummary(results)
}

// =============================================================================
// Execution
// =============================================================================

type sampleResult struct {
	sample  string
	success bool
	err     error
}

// executeSample runs all pipeline steps for one sample in order.
// Steps are run sequentially within a sample; samples run concurrently.
func executeSample(action PipelineAction, cfg AlignConfig, resolvedFasta string) sampleResult {
	color.Cyan("[%s] Starting — stage: %s\n", action.Sample, stageName(action.Stage))

	for _, step := range action.Steps {
		color.Yellow("[%s] ▶ %s  %s\n", action.Sample, step.Program, step.Input)

		if err := executeStep(step, action.Sample, cfg, resolvedFasta); err != nil {
			color.Red("[%s] ✗ %s failed: %v\n", action.Sample, step.Program, err)
			return sampleResult{sample: action.Sample, success: false, err: err}
		}
		color.Green("[%s] ✓ %s\n", action.Sample, step.Program)
	}

	color.Green("[%s] ✅ Complete\n", action.Sample)
	return sampleResult{sample: action.Sample, success: true}
}

// executeStep dispatches a single Step to the appropriate tool function.
func executeStep(step Step, sample string, cfg AlignConfig, resolvedFasta string) error {
	switch step.Program {

	case "align":
		// Full alignment from raw reads. Finds paired FASTQ files under step.Input.
		fwd, rev, err := GetReadsPE(step.Input)
		if err != nil {
			return fmt.Errorf("finding paired reads in %s: %w", step.Input, err)
		}
		if len(fwd) == 0 || len(rev) == 0 {
			return fmt.Errorf("no paired reads found in %s", step.Input)
		}
		bamDir := filepath.Join(
			filepath.Dir(step.Input),
			"reference_genomes", cfg.RefVer, "bams",
		)
		if err := os.MkdirAll(bamDir, 0o755); err != nil {
			return fmt.Errorf("creating bam dir: %w", err)
		}
		return alignment.AlignPE(
			fwd, rev, resolvedFasta, bamDir, sample,
			cfg.Aligner, cfg.Threads, cfg.Verbose, cfg.SkipVer,
		)

	case "markdup":
		outBam := markdupOutputPath(step.Input)
		return alignment.MarkDuplicates(step.Input, outBam, cfg.Threads, cfg.Verbose)

	case "bqsr":
		outBam := bqsrOutputPath(step.Input)
		return variants.RunBQSR(
			step.Input, outBam, resolvedFasta,
			cfg.KnownSites, cfg.Bootstrap,
			cfg.GatkLogLevel, cfg.Threads, cfg.Verbose,
		)

	case "bam_to_cram":
		outCram := cramOutputPath(step.Input)
		return alignment.BamToCram(step.Input, outCram, resolvedFasta, cfg.Threads, cfg.Verbose)

	case "index_cram":
		return alignment.IndexCram(step.Input, cfg.Threads, cfg.Verbose)

	case "rm":
		color.HiBlack("[%s] Removing intermediate file: %s\n", sample, step.Input)
		return os.Remove(step.Input)

	default:
		return fmt.Errorf("unknown step program %q", step.Program)
	}
}

// =============================================================================
// Output path helpers
// =============================================================================

// markdupOutputPath derives the rgmd.bam path from a sorted.bam path.
// e.g. /data/.../bams/sample.sorted.bam → /data/.../bams/sample.rgmd.bam
func markdupOutputPath(sortedBam string) string {
	base := filepath.Base(sortedBam)
	stem := base[:len(base)-len("sorted.bam")]
	return filepath.Join(filepath.Dir(sortedBam), stem+"rgmd.bam")
}

// bqsrOutputPath derives the bqsr.bam path from an rgmd bam or cram path.
func bqsrOutputPath(rgmdPath string) string {
	dir := filepath.Dir(rgmdPath)
	base := filepath.Base(rgmdPath)
	ext := filepath.Ext(base)
	stem := base[:len(base)-len(ext)]
	// Strip trailing .rgmd / .RGMD suffix from stem.
	stem = trimSuffix(stem, ".rgmd", ".RGMD")
	return filepath.Join(dir, stem+".bqsr.bam")
}

// cramOutputPath derives a .cram path from a .bam path.
func cramOutputPath(bamPath string) string {
	return bamPath[:len(bamPath)-len(filepath.Ext(bamPath))] + ".cram"
}

// trimSuffix removes the first matching suffix from s.
func trimSuffix(s string, suffixes ...string) string {
	for _, sfx := range suffixes {
		if len(s) > len(sfx) && s[len(s)-len(sfx):] == sfx {
			return s[:len(s)-len(sfx)]
		}
	}
	return s
}

// =============================================================================
// Validation & summary
// =============================================================================

func validateConfig(cfg AlignConfig) error {
	if cfg.DataDir == "" {
		return fmt.Errorf("dataDir is required")
	}
	info, err := os.Stat(cfg.DataDir)
	if err != nil {
		return fmt.Errorf("dataDir %q: %w", cfg.DataDir, err)
	}
	if !info.IsDir() {
		return fmt.Errorf("dataDir %q is not a directory", cfg.DataDir)
	}
	if cfg.Species == "" {
		return fmt.Errorf("species is required")
	}
	if cfg.RefVer == "" {
		return fmt.Errorf("refVer is required")
	}
	if cfg.BQSR {
		if cfg.Aligner == "pbmm2" {
			return fmt.Errorf("BQSR is not supported with pbmm2; use bwa-mem2 or disable BQSR")
		}
		if len(cfg.KnownSites) == 0 && !cfg.Bootstrap {
			return fmt.Errorf("BQSR requires either --known-sites or --bootstrap")
		}
		if len(cfg.KnownSites) > 0 && cfg.Bootstrap {
			return fmt.Errorf("provide either --known-sites or --bootstrap, not both")
		}
		for _, ks := range cfg.KnownSites {
			if _, err := os.Stat(ks); err != nil {
				return fmt.Errorf("known-sites file not found: %s", ks)
			}
		}
	}
	return nil
}

func printFinalSummary(results []sampleResult) {
	ok, failed := 0, 0
	for _, r := range results {
		if r.success {
			ok++
		} else {
			failed++
		}
	}

	fmt.Println()
	color.Cyan("═══════════════════ PIPELINE SUMMARY ═══════════════════\n")
	color.Green("  ✅ Completed: %d\n", ok)
	if failed > 0 {
		color.Red("  ✗  Failed:    %d\n", failed)
		fmt.Println()
		for _, r := range results {
			if !r.success {
				color.Red("     [%s] %v\n", r.sample, r.err)
			}
		}
	}
	color.Cyan("═══════════════════════════════════════════════════════\n")
}

func stageName(s AlignmentStage) string {
	return map[AlignmentStage]string{
		StageComplete: "complete",
		StageBqsr:     "bqsr → cram",
		StageRgmd:     "bqsr needed",
		StageSorted:   "markdup needed",
		StageScratch:  "full alignment",
	}[s]
}

// =============================================================================
// FASTQ validation (kept here as it feeds into the align step)
// =============================================================================

func validateFastqGz(fastq string, verbose bool, quick bool) error {
	var cmd string
	if quick {
		cmd = fmt.Sprintf("gzip -t %s", fastq)
	} else {
		cmd = fmt.Sprintf(
			`bash -c 'gzip -cd %s | awk "NR%%4==1 && !/^@/ { print \"Bad header at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } NR%%4==3 && !/^\+/ { print \"Bad separator at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } END { if(NR%%4!=0) { print \"Truncated: \", NR, \" lines\" > \"/dev/stderr\"; exit 1 } }"'`,
			fastq,
		)
	}
	if verbose {
		return utils.RunBashCmdVerbose(cmd)
	}
	return utils.RunBashCmd(cmd)
}
