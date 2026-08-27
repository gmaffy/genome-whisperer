package utils

import (
	"bufio"
	"encoding/json"
	"fmt"
	"hash/fnv"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
)

type LogEntry struct {
	Timestamp  string
	Tool       string
	Program    string
	Sample     string
	Chromosome string
	Status     string
	Cmd        string
}

type VariantCallingConfig struct {
	dataDir    string
	Reference  string
	Species    string
	RefVersion string
	Bam        string
	OutputDir  string
}

type Config struct {
	dataDir     string
	Reference   string
	GFF         string
	Proteins    string
	CDS         string
	Species     string
	OutputDir   string
	BaseName    string
	Bams        []string
	ReadPairs   [][]string
	SeReads     [][]string
	BSAseqBams  [][]string
	VCF         string
	Version     string
	VCFs        []string
	SelectChrom string
	SelectStart string
	SelectStop  string
	SelectVCF   string
	KnownSites  []string
	Threads     string
	InputDir    string
	DataType    string
	GVCFsDir    string
	CallerName  string
	GVCFs       []string
}

type HardFilterConfig struct {
	LightFilter              bool
	SNP_QD_Min               float64
	SNP_QUAL_Min             float64
	SNP_SOR_Max              float64
	SNP_FS_Max               float64
	SNP_MQ_Min               float64
	SNP_MQRankSum_Min        float64
	SNP_ReadPosRankSum_Min   float64
	INDEL_QD_Min             float64
	INDEL_QUAL_Min           float64
	INDEL_FS_Max             float64
	INDEL_ReadPosRankSum_Min float64
	INDEL_SOR_Max            float64
}

func ReadConfig(configPath string) (Config, error) {
	configFile, err := os.Open(configPath)
	if err != nil {
		return Config{}, err
	}
	defer configFile.Close()
	var cfg Config

	scanner := bufio.NewScanner(configFile)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if line == "" {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])

		switch key {
		case "Reference":
			cfg.Reference = value
		case "gff":
			cfg.GFF = value
		case "proteins":
			cfg.Proteins = value
		case "cds":
			cfg.CDS = value
		case "Species":
			cfg.Species = value
		case "OutputDir":
			cfg.OutputDir = value
		case "bam":
			cfg.Bams = append(cfg.Bams, value)
		case "BSAseqBam":
			bsaSeqBam := strings.Fields(value)
			cfg.BSAseqBams = append(cfg.BSAseqBams, bsaSeqBam)
		case "SingleEndReads":
			seRead := strings.Fields(value)
			cfg.SeReads = append(cfg.SeReads, seRead)
		case "ReadPair":
			pairs := strings.Fields(value)
			cfg.ReadPairs = append(cfg.ReadPairs, pairs)
		case "VCF":
			cfg.VCF = value
			cfg.VCFs = append(cfg.VCFs, value)
		case "gvcf":
			cfg.GVCFs = append(cfg.GVCFs, value)
		case "select_chrom":
			cfg.SelectChrom = value
		case "select_start":
			cfg.SelectStart = value
		case "select_stop":
			cfg.SelectStop = value
		case "select_vcf":
			cfg.SelectVCF = value
		case "known-sites":
			cfg.KnownSites = append(cfg.KnownSites, value)
		case "InputDir":
			cfg.InputDir = value
		case "DATA_TYPE":
			cfg.DataType = value
		case "GVCFS_Dir":
			cfg.GVCFsDir = value
		case "Caller_name":
			cfg.CallerName = value
		}
	}

	if err := scanner.Err(); err != nil {
		return cfg, err
	}

	return cfg, nil

}

func CheckDeps(deps []string) error {
	paths := make(map[string]string)
	missing := make([]string, 0)

	for _, dep := range deps {
		path, err := exec.LookPath(dep)
		if err != nil {
			missing = append(missing, dep)
			fmt.Printf("%s not found!\n", dep)
		} else {
			paths[dep] = path

		}
	}

	if len(missing) > 0 {
		return fmt.Errorf("\n\nmissing dependencies: %s", strings.Join(missing, ", "))
	}

	fmt.Println("\nDependency locations:")
	for dep, path := range paths {
		fmt.Printf("Using %s at %s\n", dep, path)
	}
	fmt.Println("")
	return nil
}

func RunBashCmdVerbose(cmdStr string) error {
	cmd := exec.Command("bash", "-c", cmdStr)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	if err != nil {
		fmt.Println("CMD error:", err)
		return err
	}
	return nil
}

func RunBashCmd(cmdStr string) error {
	cmd := exec.Command("bash", "-c", cmdStr)
	err := cmd.Run()
	if err != nil {
		fmt.Println("CMD error:", err)
		return err
	}
	return nil
}

func ParseLogFile(logFilePath string) []LogEntry {
	var data []LogEntry
	file, err := os.Open(logFilePath)
	if err != nil {
		fmt.Printf("Log file '%s' not found, starting fresh or assuming no previous runs.\n", logFilePath)
		return data
	}
	defer file.Close()
	//fmt.Println("Parsing log file ...")

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		var entry map[string]interface{}
		lineBytes := scanner.Bytes()
		if len(lineBytes) == 0 { // Skip empty lines
			continue
		}
		if err := json.Unmarshal(lineBytes, &entry); err != nil {
			fmt.Printf("Warning: Skipping malformed log line: %v - Line: %s\n", err, string(lineBytes))
			continue // skip malformed line
		}

		timestamp, timeOk := entry["time"]
		tool, toolOk := entry["msg"]
		chromVal, chromOk := entry["CHROMOSOME"]
		statusVal, statusOk := entry["STATUS"]
		sampleVal, sampleOk := entry["SAMPLE"]
		programVal, programOk := entry["PROGRAM"]
		if chromOk && statusOk && sampleOk && programOk && timeOk && toolOk {
			r := LogEntry{
				Timestamp:  timestamp.(string),
				Tool:       tool.(string),
				Program:    programVal.(string),
				Sample:     sampleVal.(string),
				Chromosome: chromVal.(string),
				Status:     statusVal.(string),
			}
			data = append(data, r)
		}

	}

	return data
}

func StageHasCompleted(logEntries []LogEntry, prog string, sample string, chrom string) bool {
	for _, entry := range logEntries {
		if entry.Program == prog && entry.Sample == sample && entry.Chromosome == chrom && entry.Status == "COMPLETED" {
			return true
		}
	}
	return false
}

func CopyFile(src, dst string) error {
	sourceFile, sErr := os.Open(src)
	if sErr != nil {
		return fmt.Errorf("couldn't open source file %s: %w", src, sErr)
	}
	defer sourceFile.Close()

	dstFile, dErr := os.Create(dst)
	if dErr != nil {
		return fmt.Errorf("couldn't create destination file %s: %w", dst, dErr)
	}
	defer dstFile.Close()

	_, err := io.Copy(dstFile, sourceFile)
	if err != nil {
		return fmt.Errorf("failed to copy file contents: %w", err)
	}

	return nil
}

func PrepareFasta(ref, aligner string, verbose bool) error {

	fmt.Printf("Starting prepare fasta with %s................\n\n", aligner)
	indexStr := ``
	indexFile := ""
	if aligner == "bwa-mem" {
		indexStr = fmt.Sprintf(`bwa index %s`, ref)
		indexFile = fmt.Sprintf("%s.bwt", ref)
	} else if aligner == "bwa-mem2" {
		indexStr = fmt.Sprintf(`bwa-mem2 index %s`, ref)
		indexFile = fmt.Sprintf("%s.bwt.2bit.64", ref)
	} else if aligner == "bowtie2" {
		indexStr = fmt.Sprintf(`bowtie2-build %s %s`, ref, ref)
		indexFile = fmt.Sprintf("%s.1.bt2", ref)
	} else if aligner == "pbmm2" {
		indexStr = fmt.Sprintf(`pbmm2 index %s`, ref)
		indexFile = fmt.Sprintf("%s.mmi", ref)
	} else {
		return fmt.Errorf("unsupported aligner: %s. Supported aligners are 'bwa-mem', bwa-mem2, 'bowtie2', 'pbmm2'", aligner)
	}

	_, bwErr := os.Stat(indexFile)
	if bwErr != nil {
		fmt.Printf("\nreference file %s is not indexed using  %s....\n\n", ref, aligner)
		fmt.Printf("Running %s ...\n\n", indexStr)

		var indexErr error
		if verbose {
			indexErr = RunBashCmdVerbose(indexStr)
		} else {
			indexErr = RunBashCmd(indexStr)
		}

		if indexErr != nil {
			return fmt.Errorf("%s indexing failed: %v", aligner, indexErr)

		}
	}

	// samtools writes both the .fai and, for a bgzipped reference, the .gzi
	// offset index in one pass, so a missing .gzi has to re-run faidx too:
	// without it nothing can seek into the compressed reference.
	_, faiErr := os.Stat(FaiPath(ref))
	gziErr := error(nil)
	if gzi := GziPath(ref); gzi != "" {
		_, gziErr = os.Stat(gzi)
	}
	if faiErr != nil || gziErr != nil {
		samIndexStr := fmt.Sprintf(`samtools faidx %s`, ref)
		fmt.Printf("Running %s ...\n\n", samIndexStr)

		fmt.Printf("\nRunning: %s ...\n\n", samIndexStr)

		var samIndexErr error
		if verbose {
			samIndexErr = RunBashCmdVerbose(samIndexStr)
		} else {
			samIndexErr = RunBashCmd(samIndexStr)
		}
		//samIndexErr := RunBashCmd(samIndexStr)
		if samIndexErr != nil {
			return fmt.Errorf("samtools indexing failed: %v", samIndexErr)
		}
	}
	// -O is passed explicitly rather than left to GATK, which derives the name
	// by stripping extensions and would write genome.dict for a genome.fa.gz —
	// not the genome.fa.gz.dict every reader here looks for.
	dictPath := DictPath(ref)
	_, dicfErr := os.Stat(dictPath)
	if dicfErr != nil {
		fmt.Printf("Running  gatk CreateSequenceDictionary -R %s -O %s ...\n\n", ref, dictPath)
		dicStr := fmt.Sprintf(`gatk CreateSequenceDictionary -R %s -O %s`, ref, dictPath)

		var dicErr error
		if verbose {
			dicErr = RunBashCmdVerbose(dicStr)
		} else {
			dicErr = RunBashCmd(dicStr)
		}

		if dicErr != nil {
			return fmt.Errorf("gatk creating sequence dictionary failed: %s", dicErr)
		}

	}

	// makeblastdb reads plain text only, so a bgzipped reference gets no blast
	// database. Nothing on the alignment path needs one; the step is skipped
	// out loud rather than failing the whole preparation over it.
	_, nsqErr := os.Stat(ref + ".nsq")
	if IsBgzippedFasta(ref) {
		if nsqErr != nil {
			fmt.Printf("Skipping makeblastdb: %s is block-compressed and makeblastdb cannot read it\n\n", ref)
		}
	} else if nsqErr != nil {
		mbdbStr := fmt.Sprintf(`makeblastdb -in %s -dbtype nucl -out %s`, ref, ref)
		fmt.Printf("Running %s ...\n\n", mbdbStr)

		fmt.Printf("\nRunning: %s ...\n\n", mbdbStr)

		var mbdbErr error
		if verbose {
			mbdbErr = RunBashCmdVerbose(mbdbStr)
		} else {
			mbdbErr = RunBashCmd(mbdbStr)
		}

		if mbdbErr != nil {
			return fmt.Errorf("makeblastdb creation failed: %v", mbdbErr)
		}
	}
	fmt.Printf("Fasta file prep done ...\n\n")

	return nil
}

func ValidateGvcf(vcf string, verbose bool, quick bool) error {
	// Which index belongs beside a VCF is decided by the VCF: a bgzipped one
	// carries the .tbi bcftools writes, an uncompressed one the Tribble .idx
	// GATK writes. On a reference with contigs past the BAI ceiling only the
	// second can be built at all.
	compressed := strings.HasSuffix(strings.ToLower(vcf), ".gz")
	index := vcf + ".idx"
	if compressed {
		index = vcf + ".tbi"
	}

	if quick {
		if _, err := os.Stat(index); err != nil {
			return fmt.Errorf("%s index missing for %s", filepath.Ext(index), vcf)
		}
		valStr := fmt.Sprintf("bcftools view -h %s > /dev/null", vcf)
		if verbose {
			fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
			return RunBashCmdVerbose(valStr)
		}
		return RunBashCmd(valStr)
	}

	// The thorough check rebuilds the index, which bcftools can only do for a
	// bgzipped file; for the rest, reading every record end to end is the
	// equivalent test of whether the file is whole.
	valStr := fmt.Sprintf("bcftools index --tbi --force %s", vcf)
	if !compressed {
		valStr = fmt.Sprintf("bcftools view %s > /dev/null", vcf)
	}
	if verbose {
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
		return RunBashCmdVerbose(valStr)
	}
	return RunBashCmd(valStr)
}

func ValidateBam(bam string, ref string, verbose bool, quick bool) error {
	var valStr string
	if quick {
		valStr = fmt.Sprintf("samtools quickcheck %s > /dev/null", bam)
	} else {
		valStr = fmt.Sprintf(`samtools view -T %s -h %s > /dev/null`, ref, bam)
	}
	if verbose {
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
		return RunBashCmdVerbose(valStr)
	}
	return RunBashCmd(valStr)
}

func ValidateFastqGz(fastq string, verbose bool, quick bool) error {
	if quick {
		valStr := fmt.Sprintf("gzip -t %s", fastq)
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
		if verbose {
			return RunBashCmdVerbose(valStr)
		}
		return RunBashCmd(valStr)
	}

	valStr := fmt.Sprintf(
		`bash -c 'gzip -cd %s | awk "NR%%4==1 && !/^@/ { print \"Bad header at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } NR%%4==3 && !/^\+/ { print \"Bad separator at record\", int(NR/4)+1 > \"/dev/stderr\"; exit 1 } END { if(NR%%4!=0) { print \"Truncated: \", NR, \"lines\" > \"/dev/stderr\"; exit 1 } }" '`,
		fastq,
	)
	fmt.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", valStr)
	if verbose {
		return RunBashCmdVerbose(valStr)
	}
	return RunBashCmd(valStr)
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

			var fastaPath, gzFastaPath, dictPath string
			for _, f := range entries {
				if f.IsDir() {
					continue
				}
				name := f.Name()
				switch {
				case IsBgzippedFasta(name):
					gzFastaPath = filepath.Join(assemblyDir, name)
				case IsFasta(name):
					fastaPath = filepath.Join(assemblyDir, name)
				case strings.HasSuffix(name, ".dict"):
					dictPath = filepath.Join(assemblyDir, name)
				}
			}

			// An assembly directory holding both forms is used uncompressed:
			// GATK's recalibration tools refuse a block-compressed reference,
			// and ReadDir order alone would pick whichever name sorts last.
			if fastaPath == "" {
				fastaPath = gzFastaPath
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

	fmt.Printf("Discovered genomes: %d\n", len(result))

	return result, nil
}

func ResolveRefFasta(refFasta, genomesDir, species, refVer string) (string, error) {
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

func GetFloat(v *vcfgo.Variant, key string) (float64, bool) {
	raw, err := v.Info().Get(key)
	if err != nil || raw == nil {
		return 0, false
	}
	switch val := raw.(type) {
	case float32:
		return float64(val), true
	case float64:
		return val, true
	case int:
		return float64(val), true
	default:
		return 0, false
	}
}

func RunCmd(cmdStr string, verbose bool) error {
	if verbose {
		return RunBashCmdVerbose(cmdStr)
	}
	return RunBashCmd(cmdStr)
}

func CalculateGATKMemory(totalSystemRamGB int, maxWorkers int) string {
	if maxWorkers <= 0 {
		maxWorkers = 1
	}

	safeRam := totalSystemRamGB - 4
	if safeRam < 4 {
		safeRam = 4 // Absolute minimum floor
	}

	memPerWorker := safeRam / maxWorkers

	if memPerWorker < 4 {
		memPerWorker = 4
	}

	return fmt.Sprintf("-Xmx%dg -Xms%dg", memPerWorker, memPerWorker)
}

// envTmpDir names the environment variable that moves every scratch directory
// this package hands out — the small ones and the large spills alike — off the
// default base in one go.
const envTmpDir = "GENOME_WHISPERER_TMPDIR"

// TmpBase returns the directory scratch goes under: $GENOME_WHISPERER_TMPDIR if
// it is set, otherwise os.TempDir(), which honours $TMPDIR and is /tmp on a
// default Linux host.
func TmpBase() string {
	if explicit := os.Getenv(envTmpDir); explicit != "" {
		return explicit
	}
	return os.TempDir()
}

// WorkTmpDir returns, creating it if needed, the scratch directory the external
// tools are told to spill into: <TmpBase()>/genome-whisperer/<key>.
//
// Everything the pipeline shells out to — samtools sort, the Picard-style GATK
// tools, the GATK walkers — is pointed here rather than at a directory beside
// its own output. Outputs live on the mounted data drives, where the small-file
// random I/O a spill is made of is the slowest thing a step does, while
// os.TempDir() is local. The base is TmpBase() rather than a hard-coded /tmp, so
// a run whose spill will not fit there can be sent elsewhere with
// $GENOME_WHISPERER_TMPDIR (or $TMPDIR) without touching the code — worth
// knowing on a host where /tmp is a RAM-backed tmpfs.
//
// The key is derived from the output's own directory, so two samples, or two
// concurrent processes, never share a scratch directory.
func WorkTmpDir(outPath string) string {
	dir := WorkTmpDirFor(filepath.Dir(outPath))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		// Nowhere to spill but beside the output, which is where these tools
		// would have gone on their own.
		return filepath.Dir(outPath)
	}
	return dir
}

// WorkTmpDirFor returns the scratch path WorkTmpDir hands out for outputs
// written to dir, without creating it. Cleanup needs the path of a directory it
// may never have created.
func WorkTmpDirFor(dir string) string {
	return filepath.Join(TmpBase(), "genome-whisperer", scratchKey(dir))
}

// SpillTmpDir is WorkTmpDir for the two steps that spill on the order of the
// read set rather than a few MB: the samtools sort at the end of an aligner pipe
// and MarkDuplicates' sorting collection. A whole-genome sort of a 40x sample
// writes tens of GB of chunks, so those two cannot go wherever the small
// scratch goes when os.TempDir() is a RAM-backed tmpfs — that is a full tmpfs
// and an OOM, not a slow step.
//
// The base is the first of these that is a directory:
//
//	$GENOME_WHISPERER_TMPDIR   taken as given, RAM-backed or not
//	os.TempDir()               unless it is a tmpfs/ramfs (honours $TMPDIR)
//	/var/tmp                   unless it is a tmpfs/ramfs
//
// and if none of them can be created, the spill goes back beside the output as
// <output dir>/tmp — slow on a mounted data drive, but never out of space.
//
// The joint-genotyping scratch in the variants package is deliberately not
// routed through here: GenomicsDBImport's temp files stay on os.TempDir().
func SpillTmpDir(outPath string) string {
	dir := SpillTmpDirFor(filepath.Dir(outPath))
	if err := os.MkdirAll(dir, 0o755); err == nil {
		return dir
	}
	beside := filepath.Join(filepath.Dir(outPath), "tmp")
	if err := os.MkdirAll(beside, 0o755); err != nil {
		return filepath.Dir(outPath)
	}
	return beside
}

// SpillTmpDirFor returns the large-spill scratch path for outputs written to
// dir, without creating it, so cleanup can find what SpillTmpDir created.
func SpillTmpDirFor(dir string) string {
	for _, base := range spillBases() {
		if info, err := os.Stat(base); err == nil && info.IsDir() {
			return filepath.Join(base, "genome-whisperer", scratchKey(dir))
		}
	}
	return filepath.Join(dir, "tmp")
}

func spillBases() []string {
	var bases []string
	if explicit := os.Getenv(envTmpDir); explicit != "" {
		bases = append(bases, explicit)
	}
	if tmp := os.TempDir(); !isRAMBacked(tmp) {
		bases = append(bases, tmp)
	}
	if !isRAMBacked("/var/tmp") {
		bases = append(bases, "/var/tmp")
	}
	return bases
}

// isRAMBacked reports whether path sits on a filesystem that is really memory:
// filling it costs the run its RAM, not its disk. Unknown filesystems — and any
// host with no /proc/mounts to read — count as real disks, so the answer only
// ever moves a spill away from a tmpfs it can prove is there.
func isRAMBacked(path string) bool {
	mounts, err := os.ReadFile("/proc/mounts")
	if err != nil {
		return false
	}
	abs, err := filepath.Abs(path)
	if err != nil {
		abs = path
	}
	abs = filepath.Clean(abs)

	var (
		bestPoint string
		bestType  string
	)
	for _, line := range strings.Split(string(mounts), "\n") {
		fields := strings.Fields(line)
		if len(fields) < 3 {
			continue
		}
		point := filepath.Clean(fields[1])
		if point != "/" && !strings.HasPrefix(point, "/") {
			continue
		}
		if abs != point && point != "/" && !strings.HasPrefix(abs, point+string(os.PathSeparator)) {
			continue
		}
		if len(point) >= len(bestPoint) {
			bestPoint, bestType = point, fields[2]
		}
	}

	switch bestType {
	case "tmpfs", "ramfs", "devtmpfs":
		return true
	}
	return false
}

// scratchKey names one output directory's scratch: a readable tail of the path
// so /tmp stays legible mid-run, plus a hash of the whole path so flattening the
// separators out of that tail cannot make two directories collide.
func scratchKey(dir string) string {
	abs, err := filepath.Abs(dir)
	if err != nil {
		abs = dir
	}
	h := fnv.New64a()
	_, _ = h.Write([]byte(abs))

	label := strings.Trim(strings.Map(func(r rune) rune {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z', r >= '0' && r <= '9', r == '.', r == '-':
			return r
		}
		return '_'
	}, abs), "_")
	if len(label) > 60 {
		label = label[len(label)-60:]
	}

	return fmt.Sprintf("%s-%x", label, h.Sum64())
}
