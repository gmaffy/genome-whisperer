package alignmentdir

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"

	//"runtime"
	"sort"
	"strings"
	//"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
)

func isFwd(filename string) bool {
	return strings.HasSuffix(filename, "1.fastq.gz") ||
		strings.HasSuffix(filename, "1.fq.gz") ||
		strings.HasSuffix(filename, "1.fastq") ||
		strings.HasSuffix(filename, "1.fq") ||
		strings.HasPrefix(filename, "Forward") ||
		strings.Contains(filename, "_R1_") ||
		strings.Contains(filename, "_1_") ||
		strings.Contains(filename, "1.fq")
}

func isRev(filename string) bool {
	return strings.HasSuffix(filename, "2.fastq.gz") ||
		strings.HasSuffix(filename, "2.fq.gz") ||
		strings.HasSuffix(filename, "2.fastq") ||
		strings.HasSuffix(filename, "2.fq") ||
		strings.HasPrefix(filename, "Reverse") ||
		strings.Contains(filename, "_R2_") ||
		strings.Contains(filename, "_2_") ||
		strings.Contains(filename, "_2.fq")
}

func GetReadsPE(cleanReadsDir string) ([]string, []string, error) {
	var fwdReads, revReads []string

	entries, err := os.ReadDir(cleanReadsDir)
	if err != nil {
		return nil, nil, err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}

		name := entry.Name()
		fullPath := filepath.Join(cleanReadsDir, name)

		if isFwd(name) {
			fwdReads = append(fwdReads, fullPath)
		} else if isRev(name) {
			revReads = append(revReads, fullPath)
		}
	}

	sort.Strings(fwdReads)
	sort.Strings(revReads)

	return fwdReads, revReads, nil
}

func validateFastqGz(fastq string, verbose bool, quick bool) error {
	quoted := shQuote(fastq)
	isGzip := strings.HasSuffix(strings.ToLower(fastq), ".gz")

	if quick {
		if isGzip {
			return runBash(fmt.Sprintf("gzip -t %s", quoted), verbose)
		}
		return runBash(fmt.Sprintf("head -n 4 %s > /dev/null", quoted), verbose)
	}

	reader := fmt.Sprintf("cat %s", quoted)
	if isGzip {
		reader = fmt.Sprintf("gzip -cd %s", quoted)
	}

	cmd := reader + ` | awk 'NR%4==1 && $0 !~ /^@/ { print "Bad header at record", int(NR/4)+1 > "/dev/stderr"; exit 1 } NR%4==3 && $0 !~ /^\+/ { print "Bad separator at record", int(NR/4)+1 > "/dev/stderr"; exit 1 } END { if (NR%4 != 0) { print "Truncated:", NR, "lines" > "/dev/stderr"; exit 1 } }' > /dev/null`
	return runBash(cmd, verbose)
}

type GenomeRef struct {
	RefVer    string
	FastaPath string
	DictPath  string
}

func GetValidGenomesFromDisk(genomesDir string) (map[string][]GenomeRef, error) {
	result := make(map[string][]GenomeRef)

	species, err := os.ReadDir(genomesDir)
	if err != nil {
		return nil, fmt.Errorf("reading genomes dir %s: %w", genomesDir, err)
	}

	for _, sp := range species {
		if !sp.IsDir() {
			continue
		}
		spKey := strings.ToUpper(sp.Name())
		spDir := filepath.Join(genomesDir, sp.Name())

		refs, err := os.ReadDir(spDir)
		if err != nil {
			continue
		}

		for _, ref := range refs {
			if !ref.IsDir() {
				continue
			}
			assemblyDir := filepath.Join(spDir, ref.Name(), "assembly")
			info, err := os.Stat(assemblyDir)
			if err != nil || !info.IsDir() {
				continue
			}

			entries, err := os.ReadDir(assemblyDir)
			if err != nil {
				continue
			}

			var fastaPath, dictPath string
			for _, f := range entries {
				if f.IsDir() {
					continue
				}
				name := f.Name()
				switch {
				case strings.HasSuffix(name, ".fa"),
					strings.HasSuffix(name, ".fasta"),
					strings.HasSuffix(name, ".fna"):
					fastaPath = filepath.Join(assemblyDir, name)
				case strings.HasSuffix(name, ".dict"):
					dictPath = filepath.Join(assemblyDir, name)
				}
			}

			if fastaPath != "" && dictPath != "" {
				result[spKey] = append(result[spKey], GenomeRef{
					RefVer:    strings.ToUpper(ref.Name()),
					FastaPath: fastaPath,
					DictPath:  dictPath,
				})
			}
		}
	}

	return result, nil
}

func resolveRefFasta(refFasta, genomesDir, species, refVer string) (string, error) {
	if refFasta != "" {
		if _, err := os.Stat(refFasta); err != nil {
			return "", fmt.Errorf("provided refFasta %q not accessible: %w", refFasta, err)
		}
		return refFasta, nil
	}
	if genomesDir == "" {
		return "", fmt.Errorf("either refFasta or genomesDir must be provided")
	}
	genomes, err := GetValidGenomesFromDisk(genomesDir)
	if err != nil {
		return "", fmt.Errorf("auto-discovering genomes: %w", err)
	}
	spKey := strings.ToUpper(species)
	refs, ok := genomes[spKey]
	if !ok {
		return "", fmt.Errorf("no genomes found for species %q in %s", species, genomesDir)
	}
	verKey := strings.ToUpper(refVer)
	for _, r := range refs {
		if r.RefVer == verKey {
			color.Green("Auto-resolved reference: %s\n", r.FastaPath)
			return r.FastaPath, nil
		}
	}
	return "", fmt.Errorf("species %q found but no assembly matching refVer %q in %s", species, refVer, genomesDir)
}

type sampleTask struct {
	sample        string
	bamDir        string
	cleanReadsDir string
	fwd           string
	rev           string
	state         SampleBamState
	//reason        string
}

type sampleResult struct {
	sample  string
	success bool
	err     error
}

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}

	if species == "" {
		fmt.Println("Please provide species name")
		return
	}
	if refVer == "" {
		fmt.Println("Please provide reference version name")
		return
	}

	resolvedFasta, resolveErr := resolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		color.Red("Could not resolve reference fasta: %v\n", resolveErr)
		return
	}

	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", resolvedFasta)
		return
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", resolvedFasta)
		return
	}

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	if bqsr {
		if aligner == "pbmm2" {
			color.Red("BQSR is not supported for pbmm2. Use bwa-mem2 or disable BQSR.\n")
			return
		}
		if len(knownSites) == 0 && !bootstrap {
			color.Red("BQSR requested: either provide known-sites or enable bootstrap.\n")
			return
		}
		if len(knownSites) > 0 && bootstrap {
			color.Red("Choose either known-sites OR bootstrap, not both.\n")
			return
		}
		for _, ks := range knownSites {
			if _, err := os.Stat(ks); err != nil {
				color.Red("Known-sites file not found: %s\n", ks)
				return
			}
		}
	}

	color.Green("===================== Getting directories with the clean_reads directories =====================\n")
	cleanReadsDirs, err := getPathsWithCleanReadsDir(dataDir, species)
	if err != nil {
		color.Red("Error finding samples dirs: %v\n", err)
		return
	}
	if len(cleanReadsDirs) == 0 {
		color.Yellow("No samples dirs found under %s/%s\n", dataDir, species)
		return
	}
	color.Green("Sample dirs found: %d\n", len(cleanReadsDirs))
	var (
		tasks           []sampleTask
		alreadyComplete int
		scanFailures    []string
		missingReads    []string
		longReadSamples []string
	)

	for _, cleanReadsDir := range cleanReadsDirs {
		sampleDir := filepath.Dir(cleanReadsDir)
		sample := filepath.Base(sampleDir)
		bamDir := filepath.Join(sampleDir, "reference_genomes", refVer, "bams")
		//cleanReadsDir := filepath.Join(sampleRoot, "clean_reads")

		state, scanErr := inspectSampleBamDir(sample, bamDir, resolvedFasta, verbose, quick)
		if scanErr != nil {
			scanFailures = append(scanFailures, sample)
			continue
		}

		rgmdCramOK := state.RgmdCram.Present && state.RgmdCram.Valid && state.RgmdCram.IndexPresent && state.RgmdCram.IndexSize > 0
		bqsrCramOK := state.BqsrCram.Present && state.BqsrCram.Valid && state.BqsrCram.IndexPresent && state.BqsrCram.IndexSize > 0

		if rgmdCramOK && bqsrCramOK && !state.SortedBam.Present && !state.RgmdBam.Present && !state.BqsrBam.Present && len(state.OtherFiles) == 0 {
			alreadyComplete++
			color.Green("[%s]- ✅ PASS. Skipping ...\n", sample)
			continue
		}

		//if sampleIsComplete(state, bqsr) {
		//	alreadyComplete++
		//	continue
		//}
		//
		//if canCompleteWithoutReads(state, bqsr) {
		//	tasks = append(tasks, sampleTask{
		//		sample:        sample,
		//		bamDir:        bamDir,
		//		cleanReadsDir: cleanReadsDir,
		//		state:         state,
		//		reason:        sampleWorkReason(state, bqsr),
		//	})
		//	continue
		//}

		if strings.HasSuffix(strings.ToUpper(sample), "LR") {
			color.Blue("[%s] Long reads sample ..\n", sample)
			longReadSamples = append(longReadSamples, sample)
			continue
		}

		fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
		if readsErr != nil || len(fwdReads) != 1 || len(revReads) != 1 {
			missingReads = append(missingReads, sample)
			continue
		}
		//tasks = append(tasks, sample)
		tasks = append(tasks, sampleTask{
			sample:        sample,
			bamDir:        bamDir,
			cleanReadsDir: cleanReadsDir,
			fwd:           fwdReads[0],
			rev:           revReads[0],
			state:         state,
			//reason:        "align from reads",
		})
	}

	color.Green("Samples complete: %d\n", alreadyComplete)
	color.Green("Samples queued:   %d\n", len(tasks))
	if len(scanFailures) > 0 {
		color.Red("Scan failures:    %d\n", len(scanFailures))
	}
	//if len(missingReads) > 0 {
	//	color.Yellow("Missing reads:    %d\n", len(missingReads))
	//}
	if len(longReadSamples) > 0 {
		color.Yellow("Long read samples:      %d\n", len(longReadSamples))
	}

	if len(tasks) == 0 {
		return
	}

	totalCores := runtime.NumCPU()
	if threads <= 0 {
		threads = totalCores
	}
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}

	sem := make(chan struct{}, maxParallelJobs)
	results := make([]sampleResult, len(tasks))
	var wg sync.WaitGroup

	for i, task := range tasks {
		wg.Add(1)
		go func(idx int, task sampleTask) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()

			//fmt.Printf("[%s] %s\n", task.sample, task.reason)
			fmt.Printf("[%s] ...\n", task.sample)
			results[idx] = processSample(task, resolvedFasta, gatkLogLevel, aligner, verbose, quick, skipVer, bqsr, bootstrap, knownSites, threads)
			if results[idx].success {
				color.Green("[%s] done\n", task.sample)
				return
			}
			color.Red("[%s] failed: %v\n", task.sample, results[idx].err)
		}(i, task)
	}

	wg.Wait()

	var failed []string
	successful := 0
	for _, result := range results {
		if result.success {
			successful++
			continue
		}
		failed = append(failed, result.sample)
	}

	fmt.Println()
	color.Green("Processed successfully: %d\n", successful)
	if len(failed) > 0 {
		color.Red("Failed:                %d\n", len(failed))
		color.Red("Failed samples:        %s\n", strings.Join(failed, ", "))
	}
	if len(missingReads) > 0 {
		color.Yellow("Skipped missing reads: %s\n", strings.Join(missingReads, ", "))
	}
	if len(longReadSamples) > 0 {
		color.Yellow("Skipped unsupported:   %s\n", strings.Join(longReadSamples, ", "))
	}
	if len(scanFailures) > 0 {
		color.Red("Scan failures:         %s\n", strings.Join(scanFailures, ", "))
	}
}

func getPathsWithCleanReadsDir(dataDir, species string) ([]string, error) {
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, err
	}

	var roots []string
	for _, match := range matches {
		info, err := os.Stat(match)
		if err != nil || !info.IsDir() {
			continue
		}
		roots = append(roots, match)
	}

	sort.Strings(roots)
	return roots, nil
}

func processSample(task sampleTask, refFasta, gatkLogLevel, aligner string, verbose, quick, skipVer, bqsr, bootstrap bool, knownSites []string, threads int) sampleResult {
	if err := os.MkdirAll(task.bamDir, 0o755); err != nil {
		return sampleResult{sample: task.sample, err: err}
	}

	state := task.state

	if !isUsable(state.RgmdCram) {
		switch {
		case isUsable(state.RgmdBam):
			color.Cyan("[%s] rgmd.cram is missing or invalid, converting rgmd.bam to rgmd.cram ...\n", task.sample)
			if err := alignment.BamToCram(state.RgmdBam.Path, refFasta, verbose); err != nil { //convertToCram(state.RgmdBam.Path, refFasta, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}

		case isUsable(state.SortedBam):
			//rgmdBam := filepath.Join(task.bamDir, task.sample+".rgmd.bam") //rgmdBamPath(task.sample, task.bamDir)
			color.Cyan("[%s] RGMD.bam is missing or invalid, MarkDuplicates sorted.bam to RGMD.bam ...\n", task.sample)
			rgmdBam := strings.TrimSuffix(state.SortedBam.Path, ".bam") + ".RGMD.bam"
			if err := alignment.MarkDuplicates(refFasta, state.SortedBam.Path, verbose, aligner, gatkLogLevel); err != nil { //runMarkDuplicates(refFasta, state.SortedBam.Path, rgmdBam, rgmdMetricsPath(task.sample, task.bamDir), gatkLogLevel, aligner, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}

			color.Cyan("[%s] Converting RGMD.bam to RGMD.cram ...\n", task.sample)
			if err := alignment.BamToCram(rgmdBam, refFasta, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}

		default:

			if task.fwd == "" || task.rev == "" {
				return sampleResult{sample: task.sample, err: fmt.Errorf("no usable alignment intermediates or paired reads found")}
			}
			color.Cyan("[%s] Aligning using reads: %s and %s ...\n", task.sample, task.fwd, task.rev)
			//color.Cyan("[%s] Aligning reads ...\n", task.sample)
			if !skipVer {

				color.Cyan("[%s] Forward: %s\n", task.sample, task.fwd)
				if err := utils.ValidateFastqGz(task.fwd, verbose, quick); err != nil { //validateFastqGz(task.fwd, verbose, quick); err != nil {
					return sampleResult{sample: task.sample, err: fmt.Errorf("forward read validation failed: %w", err)}
				}

				color.Cyan("[%s] Reverse: %s\n", task.sample, task.rev)
				if err := utils.ValidateFastqGz(task.rev, verbose, quick); err != nil { //validateFastqGz(task.rev, verbose, quick); err != nil {
					return sampleResult{sample: task.sample, err: fmt.Errorf("reverse read validation failed: %w", err)}
				}
			}

			color.Cyan("[%s] Aligning PE reads to %s using %s ...\n", task.sample, refFasta, aligner)
			sortedBam := filepath.Join(task.bamDir, task.sample+".sorted.bam") //sortedBamPath(task.sample, task.bamDir)
			if err := alignPairedReads(task.fwd, task.rev, refFasta, sortedBam, task.sample, aligner, threads, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}

			color.Cyan("[%s] Converting sorted.bam to RGMD.bam ...\n", task.sample)

			rgmdBam := filepath.Join(task.bamDir, task.sample+".rgmd.bam")                                        //rgmdBamPath(task.sample, task.bamDir)
			if err := alignment.MarkDuplicates(refFasta, sortedBam, verbose, aligner, gatkLogLevel); err != nil { //runMarkDuplicates(refFasta, sortedBam, rgmdBam, rgmdMetricsPath(task.sample, task.bamDir), gatkLogLevel, aligner, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}
			if err := alignment.BamToCram(rgmdBam, refFasta, verbose); err != nil { // convertToCram(rgmdBam, refFasta, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}
		}

		var err error
		state, err = inspectSampleBamDir(task.sample, task.bamDir, refFasta, verbose, quick)
		if err != nil {
			return sampleResult{sample: task.sample, err: err}
		}
	}

	if !isUsable(state.RgmdCram) {
		return sampleResult{sample: task.sample, err: fmt.Errorf("rgmd.cram is still missing or invalid")}
	}

	if !hasIndex(state.RgmdCram) {
		color.Cyan("[%s] Creating rgmd.cram index ...\n", task.sample)
		if err := alignment.BamIndex(state.RgmdCram.Path, verbose); err != nil {
			return sampleResult{sample: task.sample, err: err}
		}
		state, _ = inspectSampleBamDir(task.sample, task.bamDir, refFasta, verbose, quick)
	}

	if bqsr {
		if !isUsable(state.BqsrCram) {
			switch {
			case isUsable(state.BqsrBam):
				if err := alignment.BamToCram(state.BqsrBam.Path, refFasta, verbose); err != nil {
					return sampleResult{sample: task.sample, err: err}
				}

			default:
				sampleKnownSites := append([]string(nil), knownSites...)
				if len(sampleKnownSites) == 0 && bootstrap {
					var err error
					sampleKnownSites, err = createBootstrapKnownSites(refFasta, state.RgmdCram.Path, gatkLogLevel, verbose)
					if err != nil {
						return sampleResult{sample: task.sample, err: err}
					}
				}

				if _, err := runBQSR(refFasta, state.RgmdCram.Path, sampleKnownSites, verbose); err != nil {
					return sampleResult{sample: task.sample, err: err}
				}
			}

			var err error
			state, err = inspectSampleBamDir(task.sample, task.bamDir, refFasta, verbose, quick)
			if err != nil {
				return sampleResult{sample: task.sample, err: err}
			}
		}

		if !isUsable(state.BqsrCram) {
			return sampleResult{sample: task.sample, err: fmt.Errorf("bqsr.cram is still missing or invalid")}
		}
		if !hasIndex(state.BqsrCram) {
			if err := ensureIndex(state.BqsrCram.Path, verbose); err != nil {
				return sampleResult{sample: task.sample, err: err}
			}
			state, _ = inspectSampleBamDir(task.sample, task.bamDir, refFasta, verbose, quick)
		}
	}

	if err := cleanupSampleOutputs(task.bamDir, state, bqsr); err != nil {
		return sampleResult{sample: task.sample, err: err}
	}

	finalState, err := inspectSampleBamDir(task.sample, task.bamDir, refFasta, verbose, quick)
	if err != nil {
		return sampleResult{sample: task.sample, err: err}
	}
	if !sampleIsComplete(finalState, bqsr) {
		return sampleResult{sample: task.sample, err: fmt.Errorf("sample outputs are still incomplete after processing")}
	}

	return sampleResult{sample: task.sample, success: true}
}

func inspectSampleBamDir(sample, bamDir, refFasta string, verbose, quick bool) (SampleBamState, error) {
	state := SampleBamState{Sample: sample}

	info, err := os.Stat(bamDir)
	if err != nil {
		if os.IsNotExist(err) {
			return state, nil
		}
		return state, err
	}
	if !info.IsDir() {
		return state, fmt.Errorf("%s is not a directory", bamDir)
	}

	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return state, err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}

		name := entry.Name()
		lowerName := strings.ToLower(name)
		fullPath := filepath.Join(bamDir, name)

		switch {
		case strings.HasSuffix(lowerName, "sorted.bam"):
			state.SortedBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "rgmd.bam"):
			state.RgmdBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "rgmd.cram"):
			state.RgmdCram = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "bqsr.bam"):
			state.BqsrBam = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, "bqsr.cram"):
			state.BqsrCram = getAlignmentFileInfo(fullPath, refFasta, verbose, quick)
		case strings.HasSuffix(lowerName, ".bai"), strings.HasSuffix(lowerName, ".crai"), strings.HasSuffix(lowerName, ".pdf"):
			continue
		default:
			fileInfo, err := entry.Info()
			size := int64(0)
			if err == nil {
				size = fileInfo.Size()
			}
			state.OtherFiles = append(state.OtherFiles, FileInfo{
				Path:    fullPath,
				Size:    size,
				Present: true,
			})
		}
	}

	return state, nil
}

func getAlignmentFileInfo(path, refFasta string, verbose, quick bool) FileInfo {
	info := FileInfo{
		Path:    path,
		Present: true,
	}

	stat, err := os.Stat(path)
	if err != nil || !stat.Mode().IsRegular() {
		return info
	}

	info.Size = stat.Size()
	info.ValidateErr = utils.ValidateBam(path, refFasta, verbose, quick) //validateAlignmentFile(path, refFasta, verbose, quick)
	info.Valid = info.ValidateErr == nil

	if strings.HasSuffix(path, ".bam") {
		baiFile := strings.TrimSuffix(path, filepath.Ext(path)) + ".bai"
		baiFile2 := path + ".bai"
		if bInfo, err := os.Stat(baiFile); err == nil {
			info.IndexPresent = true
			info.IndexSize = bInfo.Size()
		} else if b2Info, err := os.Stat(baiFile2); err == nil {
			info.IndexPresent = true
			info.IndexSize = b2Info.Size()
		} else {
			info.IndexPresent = false
			info.IndexSize = 0
		}
	} else if strings.HasSuffix(path, ".cram") {
		craiFile := strings.TrimSuffix(path, filepath.Ext(path)) + ".crai"
		craiFile2 := path + ".crai"
		if cInfo, err := os.Stat(craiFile); err == nil {
			info.IndexPresent = true
			info.IndexSize = cInfo.Size()
		} else if c2Info, err := os.Stat(craiFile2); err == nil {
			info.IndexPresent = true
			info.IndexSize = c2Info.Size()
		} else {
			info.IndexPresent = false
			info.IndexSize = 0
		}

	} else {
		info.IndexPresent = false
		info.IndexSize = 0
	}

	return info
}

func indexCandidates(path string) []string {
	if strings.HasSuffix(strings.ToLower(path), ".bam") {
		return []string{strings.TrimSuffix(path, filepath.Ext(path)) + ".bai"}
	}
	if strings.HasSuffix(strings.ToLower(path), ".cram") {
		stem := strings.TrimSuffix(path, filepath.Ext(path))
		return []string{stem + ".crai", path + ".crai"}
	}
	return nil
}

func sampleHasRequiredFinals(state SampleBamState, requireBQSR bool) bool {
	rgmdOK := isUsable(state.RgmdCram) && hasIndex(state.RgmdCram)
	if !rgmdOK {
		return false
	}
	if !requireBQSR {
		return !state.BqsrCram.Present || (isUsable(state.BqsrCram) && hasIndex(state.BqsrCram))
	}
	return isUsable(state.BqsrCram) && hasIndex(state.BqsrCram)
}

func sampleIsComplete(state SampleBamState, requireBQSR bool) bool {
	if !sampleHasRequiredFinals(state, requireBQSR) {
		return false
	}
	return !state.SortedBam.Present &&
		!state.RgmdBam.Present &&
		!state.BqsrBam.Present &&
		len(state.OtherFiles) == 0
}

func canBuildRgmdWithoutReads(state SampleBamState) bool {
	return isUsable(state.SortedBam) ||
		isUsable(state.RgmdBam) ||
		isUsable(state.RgmdCram)
}

func canCompleteWithoutReads(state SampleBamState, requireBQSR bool) bool {
	if sampleHasRequiredFinals(state, requireBQSR) {
		return true
	}
	return canBuildRgmdWithoutReads(state) ||
		isUsable(state.BqsrBam) ||
		isUsable(state.BqsrCram)
}

func sampleWorkReason(state SampleBamState, requireBQSR bool) string {
	switch {
	case sampleHasRequiredFinals(state, requireBQSR):
		return "cleaning intermediates"
	case isUsable(state.BqsrBam):
		return "resuming from bqsr.bam"
	case isUsable(state.BqsrCram):
		return "finalizing bqsr.cram"
	case isUsable(state.RgmdCram):
		if requireBQSR {
			return "resuming from rgmd.cram"
		}
		return "finalizing rgmd.cram"
	case isUsable(state.RgmdBam):
		return "resuming from rgmd.bam"
	case isUsable(state.SortedBam):
		return "resuming from sorted.bam"
	default:
		return "align from reads"
	}
}

func supportsPairedReadDirAlign(aligner string) bool {
	switch aligner {
	case "bwa-mem", "bwa-mem2", "bowtie2":
		return true
	default:
		return false
	}
}

func sortedBamPath(sample, bamDir string) string {
	return filepath.Join(bamDir, sample+".sorted.bam")
}

func rgmdBamPath(sample, bamDir string) string {
	return filepath.Join(bamDir, sample+".rgmd.bam")
}

func rgmdMetricsPath(sample, bamDir string) string {
	return filepath.Join(bamDir, sample+".rgmd.metrics.txt")
}

func bqsrOutputPath(input string) string {
	switch {
	case strings.HasSuffix(strings.ToLower(input), ".bam"):
		return strings.TrimSuffix(input, filepath.Ext(input)) + "_bqsr.bam"
	case strings.HasSuffix(strings.ToLower(input), ".cram"):
		return strings.TrimSuffix(input, filepath.Ext(input)) + "_bqsr.cram"
	default:
		return input + "_bqsr"
	}
}

func convertToCram(src, refFasta string, verbose bool) (string, error) {
	dst := strings.TrimSuffix(src, filepath.Ext(src)) + ".cram"
	if err := removeIfExists(dst); err != nil {
		return "", err
	}
	if err := removeIndexFiles(dst); err != nil {
		return "", err
	}

	cmd := fmt.Sprintf("samtools view -C -T %s -o %s %s", shQuote(refFasta), shQuote(dst), shQuote(src))
	if err := runBash(cmd, verbose); err != nil {
		return "", err
	}
	if err := ensureIndex(dst, verbose); err != nil {
		return "", err
	}
	return dst, nil
}

func runMarkDuplicates(refFasta, sortedBam, rgmdBam, metricsPath, gatkLogLevel, aligner string, verbose bool) error {
	if err := removeIfExists(rgmdBam, metricsPath); err != nil {
		return err
	}
	if err := removeIndexFiles(rgmdBam); err != nil {
		return err
	}

	var cmd string
	switch aligner {
	case "pbmm2":
		cmd = fmt.Sprintf("pbmm2 markdup %s %s", shQuote(sortedBam), shQuote(rgmdBam))
	default:
		cmd = fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
			shQuote(refFasta), shQuote(sortedBam), shQuote(rgmdBam), shQuote(metricsPath), gatkLogLevel)
	}
	return runBash(cmd, verbose)
}

func alignPairedReads(fwd, rev, refFasta, sortedBam, sample, aligner string, threads int, verbose bool) error {
	if err := removeIfExists(sortedBam); err != nil {
		return err
	}
	if err := removeIndexFiles(sortedBam); err != nil {
		return err
	}

	lib := fmt.Sprintf("%s_1", sample)
	var err error
	switch aligner {
	case "bwa-mem2":
		err = alignment.BwaMem2Align(fwd, rev, refFasta, sample, lib, threads, sortedBam, verbose)
	case "bwa-mem":
		err = alignment.BwaMemAlign(fwd, rev, refFasta, sample, lib, threads, sortedBam, verbose)
	case "bowtie2":
		err = alignment.Bowtie2Align(fwd, rev, refFasta, sortedBam, sample, lib, threads, verbose)
	default:
		return fmt.Errorf("aligner %s is not supported for paired-read directory alignment", aligner)
	}
	return err
}

func ensureIndex(path string, verbose bool) error {
	if err := removeIndexFiles(path); err != nil {
		return err
	}
	cmd := fmt.Sprintf("samtools index %s", shQuote(path))
	return runBash(cmd, verbose)
}

func runBQSR(refFasta, input string, knownSites []string, verbose bool) (string, error) {
	output := bqsrOutputPath(input)
	base := strings.TrimSuffix(input, filepath.Ext(input))
	recalTable := base + "recal_table.txt"
	recalTable2 := base + "recal_table2.txt"
	plots := base + "recal_table_plots.pdf"

	if err := removeIfExists(output, recalTable, recalTable2, plots); err != nil {
		return "", err
	}
	if err := removeIndexFiles(output); err != nil {
		return "", err
	}

	var args []string
	for _, site := range knownSites {
		args = append(args, "--known-sites "+shQuote(site))
	}
	knownSitesArgs := strings.Join(args, " ")

	cmd1 := fmt.Sprintf("gatk BaseRecalibrator -R %s -I %s %s -O %s",
		shQuote(refFasta), shQuote(input), knownSitesArgs, shQuote(recalTable))
	if err := runBash(cmd1, verbose); err != nil {
		return "", err
	}

	cmd2 := fmt.Sprintf("gatk ApplyBQSR -R %s -I %s -bqsr %s -O %s",
		shQuote(refFasta), shQuote(input), shQuote(recalTable), shQuote(output))
	if err := runBash(cmd2, verbose); err != nil {
		return "", err
	}

	cmd3 := fmt.Sprintf("gatk BaseRecalibrator -R %s -I %s %s -O %s",
		shQuote(refFasta), shQuote(output), knownSitesArgs, shQuote(recalTable2))
	if err := runBash(cmd3, verbose); err != nil {
		return "", err
	}

	cmd4 := fmt.Sprintf("gatk AnalyzeCovariates -before %s -after %s -plots %s",
		shQuote(recalTable), shQuote(recalTable2), shQuote(plots))
	if err := runBash(cmd4, verbose); err != nil {
		return "", err
	}

	if err := ensureIndex(output, verbose); err != nil {
		return "", err
	}

	return output, nil
}

func createBootstrapKnownSites(refFasta, input, gatkLogLevel string, verbose bool) ([]string, error) {
	base := strings.TrimSuffix(input, filepath.Ext(input))
	rawVCF := base + ".raw.vcf.gz"
	snpVCF := base + ".raw.SNP.vcf.gz"
	indelVCF := base + ".raw.INDEL.vcf.gz"
	snpColumns := base + ".raw.SNP.columns.vcf.gz"
	indelColumns := base + ".raw.INDEL.columns.vcf.gz"
	filteredSNP := base + ".raw.SNP.hard_filtered.vcf.gz"
	filteredINDEL := base + ".raw.INDEL.hard_filtered.vcf.gz"

	if err := removeIfExists(rawVCF, snpVCF, indelVCF, snpColumns, indelColumns, filteredSNP, filteredINDEL); err != nil {
		return nil, err
	}

	cmds := []string{
		fmt.Sprintf("gatk HaplotypeCaller -R %s -I %s -O %s", shQuote(refFasta), shQuote(input), shQuote(rawVCF)),
		fmt.Sprintf("gatk SelectVariants -V %s --select-type-to-include SNP -O %s --verbosity %s", shQuote(rawVCF), shQuote(snpVCF), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants -V %s --select-type-to-include INDEL -O %s --verbosity %s", shQuote(rawVCF), shQuote(indelVCF), gatkLogLevel),
		fmt.Sprintf(`gatk VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O %s --verbosity %s`, shQuote(snpVCF), shQuote(snpColumns), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants --exclude-filtered -V %s -O %s --verbosity %s", shQuote(snpColumns), shQuote(filteredSNP), gatkLogLevel),
		fmt.Sprintf(`gatk VariantFiltration -V %s -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O %s --verbosity %s`, shQuote(indelVCF), shQuote(indelColumns), gatkLogLevel),
		fmt.Sprintf("gatk SelectVariants --exclude-filtered -V %s -O %s --verbosity %s", shQuote(indelColumns), shQuote(filteredINDEL), gatkLogLevel),
	}

	for _, cmd := range cmds {
		if err := runBash(cmd, verbose); err != nil {
			return nil, err
		}
	}

	return []string{filteredSNP, filteredINDEL}, nil
}

func cleanupSampleOutputs(bamDir string, state SampleBamState, requireBQSR bool) error {
	keep := make(map[string]struct{})

	if isUsable(state.RgmdCram) && hasIndex(state.RgmdCram) {
		keep[state.RgmdCram.Path] = struct{}{}
		for _, idx := range existingIndexFiles(state.RgmdCram.Path) {
			keep[idx] = struct{}{}
		}
	}

	if isUsable(state.BqsrCram) && hasIndex(state.BqsrCram) {
		keep[state.BqsrCram.Path] = struct{}{}
		for _, idx := range existingIndexFiles(state.BqsrCram.Path) {
			keep[idx] = struct{}{}
		}
	}

	if requireBQSR && len(keep) < 4 {
		return fmt.Errorf("cannot clean up before both final CRAMs are ready")
	}
	if !requireBQSR && len(keep) == 0 {
		return fmt.Errorf("cannot clean up before rgmd.cram is ready")
	}

	entries, err := os.ReadDir(bamDir)
	if err != nil {
		return err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		fullPath := filepath.Join(bamDir, entry.Name())
		if _, ok := keep[fullPath]; ok {
			continue
		}
		if err := os.Remove(fullPath); err != nil && !os.IsNotExist(err) {
			return err
		}
	}

	return nil
}

func existingIndexFiles(path string) []string {
	var files []string
	for _, candidate := range indexCandidates(path) {
		if _, err := os.Stat(candidate); err == nil {
			files = append(files, candidate)
		}
	}
	return files
}

func removeIfExists(paths ...string) error {
	for _, path := range paths {
		if path == "" {
			continue
		}
		if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
			return err
		}
	}
	return nil
}

func removeIndexFiles(path string) error {
	return removeIfExists(indexCandidates(path)...)
}

func runBash(cmdStr string, verbose bool) error {
	cmd := exec.Command("bash", "-c", cmdStr)
	if verbose {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
	}
	return cmd.Run()
}

func shQuote(value string) string {
	if value == "" {
		return "''"
	}
	return "'" + strings.ReplaceAll(value, `'`, `'\''`) + "'"
}
