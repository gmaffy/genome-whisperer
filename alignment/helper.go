package alignment

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// samThreads bounds an -@ value. samtools' BGZF and CRAM codecs stop scaling
// much past 8 extra threads, and threads here is a per-sample budget that
// several samples spend concurrently, so it is a ceiling rather than a target.
func samThreads(threads int) int {
	if threads < 1 {
		return 1
	}
	if threads > 8 {
		return 8
	}
	return threads
}

// MaxCycleValue is the --maximum-cycle-value passed to BaseRecalibrator.
//
// The Cycle covariate keys on the base's position within its read, and GATK
// rejects the whole run when a read is longer than this bound — its own default
// of 500 dies on any library whose read pairs were stitched or merged before
// alignment (a CLC-merged 2x250 run reaches 512). The bound only sizes the
// covariate's key space, which GATK allocates lazily per observed quality, so
// raising it to cover any short-read library costs a few MB of heap rather than
// the "more memory" the error message warns about. BaseRecalibrator writes the
// value into the recalibration table's argument section, so ApplyBQSR,
// GatherBQSRReports and AnalyzeCovariates all pick it up from there.
const MaxCycleValue = 10000

// WorkTmpDir is utils.WorkTmpDir: the scratch directory under os.TempDir() that
// every tool the pipeline shells out to is told to spill into. It stays exported
// here because most of its callers are in this package and in alignmentdir.
func WorkTmpDir(outPath string) string {
	return utils.WorkTmpDir(outPath)
}

// IndexPath returns the index file samtools writes for an alignment file.
//
// BAM is indexed as CSI, not BAI. A BAI cannot address a position beyond
// 2^29-1 (536,870,911) on a single reference, and large genomes exceed that —
// onion's chromosomes arrive as ~1.06 Gb pieces — so both samtools and htsjdk
// refuse to build a BAI for them.
//
// CRAM's .crai has no such limit, but that only helps samtools. htsjdk has no
// CRAM index reader of its own: it converts a .crai into an in-memory BAI when
// it opens the file, so GATK inherits the 2^29 ceiling from the index it builds
// rather than from the one on disk. That is why the GATK steps are handed a
// CSI-indexed BAM on such references — see gatkReadsNeedBam in
// alignmentdir/scatter.go.
func IndexPath(path string) string {
	if strings.HasSuffix(strings.ToLower(path), ".cram") {
		return path + ".crai"
	}
	return path + ".csi"
}

// BamToCram encodes a BAM as the CRAM beside it.
//
// Written under a temporary name and renamed, like CramToBam and BamIndex. The
// cram is the sample's durable artefact, so an interrupted encode leaving a
// truncated file under that name is the worst shape this can fail in: the next
// run finds something that looks finished, and where the input was an adopted
// external alignment, recovering means normalising and re-marking the whole
// thing again.
func BamToCram(bamPath, refFasta string, threads int, verbose bool) error {
	cramPath := strings.TrimSuffix(bamPath, filepath.Ext(bamPath)) + ".cram"
	tmpPath := cramPath + ".partial"

	bamToCramStr := fmt.Sprintf(`samtools view -@ %d -T %s -C --output-fmt cram,version=3.0 -o %s %s`, samThreads(threads), refFasta, tmpPath, bamPath)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", bamToCramStr)
	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(bamToCramStr)
	} else {
		err = utils.RunBashCmd(bamToCramStr)
	}
	if err != nil {
		_ = os.Remove(tmpPath)
		return err
	}
	return os.Rename(tmpPath, cramPath)
}

// CramToBam decodes a CRAM back to a BAM beside it and returns the new path.
//
// This exists for one reason: on a reference with contigs past the BAI limit
// GATK cannot read an indexed CRAM at all, so a sample whose only surviving
// alignment is the cram has to be decoded before BQSR can run on it. The BAM is
// an intermediate — it is written under a temporary name so an interrupted
// decode cannot be mistaken for a finished one, and the sample cleanup removes
// it once the crams are in place.
func CramToBam(cramPath, refFasta string, threads int, verbose bool) (string, error) {
	bamPath := strings.TrimSuffix(cramPath, filepath.Ext(cramPath)) + ".bam"
	tmpPath := bamPath + ".partial"

	cmdStr := fmt.Sprintf(`samtools view -@ %d -T %s -b -o %s %s`, samThreads(threads), refFasta, tmpPath, cramPath)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	if err != nil {
		_ = os.Remove(tmpPath)
		return "", err
	}
	if err := os.Rename(tmpPath, bamPath); err != nil {
		return "", err
	}
	return bamPath, nil
}

func BwaMem2Align(forwardPath string, reversePath string, referencePath string, sampleName string, libName string, threads int, sortedBam string, verbose bool) error {

	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	cmdStr := fmt.Sprintf(`bwa-mem2 mem -t %v -M -Y -R '%s' %s %s %s | %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortCmd(sortedBam, threads))
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	return err
}

func BwaMemAlign(forwardPath string, reversePath string, referencePath string, sampleName string, libName string, threads int, sortedBam string, verbose bool) error {
	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	cmdStr := fmt.Sprintf(`bwa mem -t %v -M -Y -R '%s' %s %s %s | %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortCmd(sortedBam, threads))
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	return err
}

func Bowtie2Align(forwardPath string, reversePath string, referencePath string, sortedBam string, sampleName string, libName string, threads int, verbose bool) error {
	cmdStr := fmt.Sprintf(`bowtie2 -I 0 -X 1000 -x %s -1 %s -2 %s --end-to-end --sensitive --threads %v  --rg-id %s.1 --rg PL:BGISEQ --rg SM:%s --rg LB:%s | %s`, referencePath, forwardPath, reversePath, threads, sampleName, sampleName, libName, sortCmd(sortedBam, threads))
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", cmdStr)
	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	return err
}

// sortPrefix is the -T a samtools sort is given: the large-spill scratch
// directory, plus the output's own name so that two sorts sharing that
// directory cannot pick up each other's chunks.
//
// This is utils.SpillTmpDir rather than WorkTmpDir because a whole-genome sort
// writes chunks totalling roughly the size of the BAM, which is more than a
// tmpfs /tmp can hold.
func sortPrefix(sortedBam string) string {
	return filepath.Join(utils.SpillTmpDir(sortedBam), "sort_"+filepath.Base(sortedBam))
}

// sortCmd builds the samtools sort half of an aligner pipe.
//
// The default -m of 768 MiB on a single thread is the reason a whole-genome
// sort takes days here: a 15 Gb reference at 15x produces a BAM large enough to
// spill into hundreds of temporary chunks, merged by one thread. -@/-m cut both
// the chunk count and the merge time, and -T puts the chunks on a local disk
// instead of wherever samtools would default to.
func sortCmd(sortedBam string, threads int) string {
	return fmt.Sprintf(`samtools sort -@ %d -m 1G -T %s -o %s`,
		samThreads(threads), sortPrefix(sortedBam), sortedBam)
}

// EnsurePbmm2Index returns the reference to hand pbmm2, building the minimizer
// index once if it is missing.
//
// Callers build this before starting any concurrent work. Pbmm2Align used to
// build it lazily itself, which under the directory pipeline meant N long-read
// samples racing to write one file and each reading whatever the others had
// half-written.
//
// The name carries the preset because pbmm2 index bakes -k and -w into the
// index and pbmm2 align does not re-derive them from --preset: aligning HIFI
// against a SUBREAD-built index does not fail, it warns and seeds wrongly. A
// silent change in alignment quality is worth two files.
//
// A reference this cannot index is not a failed run: pbmm2 accepts a FASTA, so
// a block-compressed reference or a read-only assembly store falls back to
// handing it the FASTA and paying an in-memory index build per sample.
func EnsurePbmm2Index(refFasta, preset string, verbose bool) (string, error) {
	if utils.IsBgzippedFasta(refFasta) {
		fmt.Printf("Reference %s is block-compressed: handing pbmm2 the FASTA instead of a .mmi\n", refFasta)
		return refFasta, nil
	}

	mmi := fmt.Sprintf("%s.%s.mmi", refFasta, strings.ToUpper(preset))
	if info, err := os.Stat(mmi); err == nil && info.Size() > 0 {
		fmt.Printf("Using pbmm2 index: %s\n", mmi)
		return mmi, nil
	}

	// An index built before the name carried a preset, or by hand. Reuse it
	// rather than spending an hour rebuilding what is on disk — indexing a
	// plant genome is not cheap, and for SUBREAD, CCS and HIFI the minimizer
	// parameters are the same anyway (-k 19 -w 10). Say that its preset is
	// unknown, because for ISOSEQ and UNROLLED they are not.
	legacy := refFasta + ".mmi"
	if info, err := os.Stat(legacy); err == nil && info.Size() > 0 {
		fmt.Printf("Using the existing pbmm2 index: %s\n", legacy)
		fmt.Printf("(built before indexes were named per preset — it may not carry %s's -k/-w;\n", strings.ToUpper(preset))
		fmt.Printf(" delete it to have one built for this preset)\n")
		return legacy, nil
	}

	// Built under a process-unique name and renamed, so two runs sharing a
	// genomes directory produce identical bytes and the last rename wins
	// atomically, rather than each writing into the name the other is reading.
	partial := fmt.Sprintf("%s.partial.%d", mmi, os.Getpid())
	cmdStr := fmt.Sprintf(`pbmm2 index --preset %s %s %s`, strings.ToUpper(preset), refFasta, partial)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	if err != nil {
		_ = os.Remove(partial)
		fmt.Printf("Could not build %s: %v\nHanding pbmm2 the FASTA instead; each long-read sample will index in memory\n", mmi, err)
		return refFasta, nil
	}
	if renameErr := os.Rename(partial, mmi); renameErr != nil {
		_ = os.Remove(partial)
		return "", fmt.Errorf("installing %s: %w", mmi, renameErr)
	}
	return mmi, nil
}

// pbmm2AlignCmd builds the pbmm2 half of the long-read pipe, plus the sort.
//
// Split out so the command string is testable: pbmm2 is not installed
// everywhere this builds, and the read-group and reference arguments are
// exactly the things that are wrong in silence rather than loudly.
func pbmm2AlignCmd(sePath, refIndex, sortedBam, sampleName, libName, preset string, threads int) string {
	// Stamped at align time, the way the short-read aligners do it, rather than
	// by a second gatk AddOrReplaceReadGroups pass over a whole-genome BAM.
	// pbmm2 documents --rg for FASTA/Q input, which is the only input this path
	// takes. It is not cosmetic: HaplotypeCaller names the gVCF's sample column
	// from SM.
	readGroup := fmt.Sprintf(`@RG\tID:%s.1\tSM:%s\tLB:%s\tPL:PACBIO`, sampleName, sampleName, libName)

	// Piped into the repo's own sort rather than pbmm2 --sort. pbmm2 sorts to
	// its own temporary files in a location this pipeline does not control,
	// while every GATK step here spills under utils.TmpBase — so concurrent
	// long-read sorts would fill the directory the short-read MarkDuplicates
	// and BaseRecalibrator jobs depend on, and scratchDirs could not clean up
	// what was left behind. Piping also makes -j the whole thread budget:
	// pbmm2 --sort spends -J sort threads on top of -j.
	return fmt.Sprintf(`pbmm2 align -j %d --preset %s --rg '%s' %s %s | %s`,
		threads, strings.ToUpper(preset), readGroup, refIndex, sePath, sortCmd(sortedBam, threads))
}

// Pbmm2Align aligns one long-read FASTQ to sortedBam.
//
// refIndex is the reference to align against — a preset-keyed .mmi, or the
// FASTA where one could not be built. It is passed in rather than derived:
// this used to build <ref>.mmi lazily if it was missing, which under the
// directory pipeline meant every concurrent long-read sample racing to write
// one file. See alignmentdir.ensurePbmm2Index, which builds it once per run.
func Pbmm2Align(sePath, refIndex, sortedBam, sampleName, libName, preset string, threads int, verbose bool) error {
	pbmm2CmdStr := pbmm2AlignCmd(sePath, refIndex, sortedBam, sampleName, libName, preset, threads)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", pbmm2CmdStr)
	var pbmm2Err error
	if verbose {
		pbmm2Err = utils.RunBashCmdVerbose(pbmm2CmdStr)
	} else {
		pbmm2Err = utils.RunBashCmd(pbmm2CmdStr)
	}

	return pbmm2Err
}

func ReadGroups(sortedBam string, rgBam string, sampleName string, libName string, verbose bool) error {
	rgCmdStr := fmt.Sprintf(`gatk AddOrReplaceReadGroups -I %s -O %s -ID %s.1 -LB %s -PL PACBIO -PU BKD -SM %s --TMP_DIR %s`, sortedBam, rgBam, sampleName, libName, sampleName, WorkTmpDir(rgBam))
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", rgCmdStr)
	var rgErr error
	if verbose {
		rgErr = utils.RunBashCmdVerbose(rgCmdStr)
	} else {
		rgErr = utils.RunBashCmd(rgCmdStr)
	}
	return rgErr
}

// BamIndex indexes a BAM or CRAM.
//
// The index is built under a temporary name and renamed on success. An index
// interrupted part-way through is otherwise left behind, and a partial index is
// still present and non-empty as far as a directory scan can tell, so the next
// run reuses it and hands a truncated index to GATK.
func BamIndex(bam string, threads int, verbose bool) error {
	idxPath := IndexPath(bam)
	tmpPath := idxPath + ".partial"

	csiFlag := ""
	if !strings.HasSuffix(strings.ToLower(bam), ".cram") {
		csiFlag = "-c "
	}

	indexCmdStr := fmt.Sprintf(`samtools index %s-@ %d -o %s %s`, csiFlag, samThreads(threads), tmpPath, bam)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", indexCmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(indexCmdStr)
	} else {
		err = utils.RunBashCmd(indexCmdStr)
	}
	if err != nil {
		_ = os.Remove(tmpPath)
		return err
	}
	return os.Rename(tmpPath, idxPath)
}

// RgmdBamPath returns the file MarkDuplicates writes for a given sorted BAM.
//
// The ".sorted" tag is dropped on the way through so the durable artefacts stay
// <sample>.RGMD.*, which is the name the data directories, the resume scan and
// the publishing scripts all expect. Both MarkDuplicates and its caller derive
// the name from here; deriving it twice is how they came to disagree.
func RgmdBamPath(sortedBam string) string {
	base := strings.TrimSuffix(sortedBam, filepath.Ext(sortedBam))
	base = strings.TrimSuffix(base, ".sorted")
	return base + ".RGMD.bam"
}

// Duplicate markers, as chosen by the caller rather than inferred from which
// aligner produced the file.
const (
	// DupMarkerGatk is GATK MarkDuplicates: works on any coordinate-sorted BAM,
	// flags duplicates in place and writes a metrics file.
	DupMarkerGatk = "gatk"
	// DupMarkerPbmarkdup is pbmarkdup, for NATIVE PacBio BAM only.
	//
	// It needs the per-read tags a PacBio instrument writes — zm, np, rq — and
	// a BAM aligned from a FASTQ has none of them, because a FASTQ cannot carry
	// them. Given one it still reports the right duplicate counts and then dies
	// writing the output ("ERROR: stoul"), leaving a header-only BAM. So this is
	// only correct where the reads arrived as a PacBio uBAM; the directory
	// pipeline, whose long-read input is a FASTQ, uses DupMarkerGatk.
	DupMarkerPbmarkdup = "pbmarkdup"
)

// markDupCmd builds the duplicate-marking command and returns it with the file
// it writes.
//
// Split out from MarkDuplicates so the command string can be tested without a
// subprocess. That is worth a function on its own here: the tool name lives
// inside a format string, where a wrong one is invisible to every other kind of
// test and shows up only as a failed run hours in. This branch shipped for a
// while as "pbmm2 markdup", which is not a pbmm2 subcommand at all — pbmarkdup
// is a separate binary, which is why cmd/AlignReads.go adds it to the
// dependency list for that aligner.
func markDupCmd(referencePath, sortedBam, dupMarker, gatkLogLevel, javaOpts string, threads int) (cmd, rgmdBam string) {
	rgmdBam = RgmdBamPath(sortedBam)

	if dupMarker == DupMarkerPbmarkdup {
		// No metrics file: pbmarkdup does not write one. Anything reading
		// <sample>.RGMD.metrics.txt must tolerate its absence.
		return fmt.Sprintf(`pbmarkdup -j %d %s %s`, samThreads(threads), sortedBam, rgmdBam), rgmdBam
	}

	rgmdMetrics := strings.TrimSuffix(rgmdBam, ".bam") + ".metrics.txt"

	// --TMP_DIR matters as much as the heap: duplicate marking a
	// whole-genome BAM spills a sorting collection roughly the size of the
	// read set. SpillTmpDir keeps that off both the mounted data drive the
	// output lives on and a RAM-backed /tmp — see utils.SpillTmpDir for the
	// order it picks a disk in.
	return fmt.Sprintf(`gatk --java-options "%s" MarkDuplicates -R %s -I %s -O %s -M %s --TMP_DIR %s --VERBOSITY %s`,
		javaOpts, referencePath, sortedBam, rgmdBam, rgmdMetrics, utils.SpillTmpDir(rgmdBam), gatkLogLevel), rgmdBam
}

// MarkDuplicates marks duplicates in sortedBam, writing RgmdBamPath(sortedBam).
//
// dupMarker is one of the DupMarker* constants. It is passed in rather than
// derived from the aligner because the two are not the same question: what
// aligned a file says nothing about whether the file carries the native PacBio
// tags pbmarkdup needs.
func MarkDuplicates(referencePath string, sortedBam string, verbose bool, dupMarker string, gatkLogLevel string, javaOpts string, threads int) error {
	cmd, rgmdBam := markDupCmd(referencePath, sortedBam, dupMarker, gatkLogLevel, javaOpts, threads)

	// pbmarkdup refuses to overwrite an existing output, so a resumed sample
	// would fail on the file its own previous attempt left behind. GATK
	// MarkDuplicates overwrites happily, but clearing the target first is the
	// right thing on both paths: a half-written file from an interrupted run
	// must never be mistaken for this run's output.
	if err := os.Remove(rgmdBam); err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("clearing %s before marking duplicates: %w", rgmdBam, err)
	}

	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", cmd)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmd)
	} else {
		err = utils.RunBashCmd(cmd)
	}
	return err

}

func BQSR(bam string, verbose bool, bootrap bool, knownSites []string) error {
	return nil
}
