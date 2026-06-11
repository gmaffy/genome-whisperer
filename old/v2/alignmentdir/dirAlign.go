package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
)

type sampleWork struct {
	sample  string
	cram    string
	cramDir string // filepath.Dir of the cram, used to derive gvcf output path
}

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

	return fwdReads, revReads, nil
}

func validateFastqGz(fastq string, verbose bool, quick bool) error {
	if quick {
		// Just check the gzip magic bytes and EOF integrity
		valStr := fmt.Sprintf("gzip -t %s", fastq)
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
		if verbose {
			return utils.RunBashCmdVerbose(valStr)
		}
		return utils.RunBashCmd(valStr)
	}

	// Full validation: gzip integrity + FASTQ format check (4-line records, quality scores)
	valStr := fmt.Sprintf(
		`bash -c 'gzip -cd %s | awk "NR%%4==1 && !/^@/ { print \"Bad header at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } NR%%4==3 && !/^\+/ { print \"Bad separator at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } END { if(NR%%4!=0) { print \"Truncated: \", NR, \"lines\" > \"/dev/stderr\"; exit 1 } }" '`,
		fastq,
	)
	fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
	if verbose {
		return utils.RunBashCmdVerbose(valStr)
	}
	return utils.RunBashCmd(valStr)
}

type GenomeRef struct {
	RefVer    string // e.g. "GRCh38"
	FastaPath string // path to .fa / .fasta / .fna
	DictPath  string // path to .dict
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
				case strings.HasSuffix(name, ".fa") ||
					strings.HasSuffix(name, ".fasta") ||
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

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
	// ============================================= Validate paths ================================================ //
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		fmt.Printf("Error accessing data directory: %s\n", dataDir)
		return
	}
	if !dInfo.IsDir() {
		fmt.Printf("Data directory %s is not a directory\n", dataDir)
		return
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		fmt.Printf("Error getting absolute path for data directory: %s\n", dataDir)
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

	// ---- Resolve reference fasta (explicit path or auto-discover) ---- //
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

	// ---- Validate BQSR parameters up-front ---- //
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

	color.Green("All file paths valid\n....................................................\n\n")
	fmt.Println("data dir abs:", dataDirAbs)

	// ========================================= Discover samples =================================================== //
	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "clean_reads")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		panic(err)
	}

	color.Green("SAMPLES FOUND:\n---------------------------------------------------------------------\n\n ")
	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, match)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n", len(samples), color.GreenString(species))

	// ============================================= Scan alignment files =========================================== //

	bamsStatResults, _ := ScanAlignments(dataDir, species, refVer, genomesDir, resolvedFasta, verbose, quick)

	var failedSamples []string
	for _, s := range bamsStatResults {
		sampleAction := GetSampleActions(s)
		color.Cyan("[%s] Things to do\n----------------------------------------\n", s.Sample)
		if sampleAction.Perfect {
			color.Green("Sample %s is already aligned and indexed. Skipping.\n", s.Sample)
		} else {

			for _, step := range sampleAction.Steps {
				color.Blue("[%s] %s on %s\n", s.Sample, step.Program, strings.Join(step.Inputs, ", "))
				prog := step.Program
				input := step.Inputs[0]

				switch prog {
				case "bam_to_cram":
					cramPath := strings.Replace(input, ".bam", ".cram", 1)
					err = alignment.BamToCram(input, cramPath, resolvedFasta, verbose)
					if err != nil {
						failedSamples = append(failedSamples, s.Sample)
						color.Red("Error converting %s to cram: %v\n", input, err)
						continue
					}
				case "markdup":
					//rgmdBam := strings.Replace(input, ".bam", ".rgmd.bam", 1)
					err = alignment.MarkDuplicates(resolvedFasta, input, verbose, aligner, gatkLogLevel)
					if err != nil {
						failedSamples = append(failedSamples, s.Sample)
						color.Red("[%s] Error marking duplicates: %v\n", s.Sample, err)
						continue

					}
				case "indexBam":
					err = alignment.BamIndex(input, verbose)
					if err != nil {
						failedSamples = append(failedSamples, s.Sample)
						color.Red("[%s] Error indexing %s: %v\n", s.Sample, input, err)
						continue
					}

				case "scratch":
					color.Yellow("Aligning %s\n", input)

				case "rm":
					err = os.Remove(input)
					if err != nil {
						failedSamples = append(failedSamples, s.Sample)
					}
				}

			}
		}

	}
}
