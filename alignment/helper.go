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

func BamToCram(bamPath, refFasta string, threads int, verbose bool) error {
	cramPath := strings.TrimSuffix(bamPath, filepath.Ext(bamPath)) + ".cram"
	bamToCramStr := fmt.Sprintf(`samtools view -@ %d -T %s -C --output-fmt cram,version=3.0 -o %s %s`, samThreads(threads), refFasta, cramPath, bamPath)
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s ...\n------------------------------------------------------------------\n\n", bamToCramStr)
	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(bamToCramStr)
	} else {
		err = utils.RunBashCmd(bamToCramStr)
	}
	return err
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

func Pbmm2Align(sePath string, referencePath string, sortedBam string, preset string, threads int, verbose bool) error {
	_, refIndexErr := os.Stat(referencePath + ".mmi")

	if refIndexErr != nil {
		fmt.Println("Reference index not found")
		fmt.Println("Indexing reference ...")
		indexCmdStr := fmt.Sprintf(`pbmm2 index %s %s`, referencePath, referencePath+".mmi")
		fmt.Println(indexCmdStr)
		var indexErr error
		if verbose {
			indexErr = utils.RunBashCmdVerbose(indexCmdStr)
		} else {
			indexErr = utils.RunBashCmd(indexCmdStr)
		}
		if indexErr != nil {
			fmt.Println("Indexing reference failed")
			return indexErr
		}
	}

	pbmm2CmdStr := fmt.Sprintf(`pbmm2 align --sort -j %v --preset %s %s %s %s`, threads, preset, referencePath+".mmi", sePath, sortedBam)
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
	rgCmdStr := fmt.Sprintf(`gatk AddOrReplaceReadGroups -I %s -O %s -ID %s.1 -LB %s -PL PACBIO -PU BKD -SM %s --tmp-dir %s`, sortedBam, rgBam, sampleName, libName, sampleName, WorkTmpDir(rgBam))
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

func MarkDuplicates(referencePath string, sortedBam string, verbose bool, aligner string, gatkLogLevel string, javaOpts string) error {

	rgmdBam := RgmdBamPath(sortedBam)
	rgmdMetrics := strings.TrimSuffix(rgmdBam, ".bam") + ".metrics.txt"

	var cmd string
	if aligner == "pbmm2" {
		cmd = fmt.Sprintf(`pbmm2 markdup %s %s`, sortedBam, rgmdBam)
	} else {
		// --TMP_DIR matters as much as the heap: duplicate marking a
		// whole-genome BAM spills a sorting collection roughly the size of the
		// read set. SpillTmpDir keeps that off both the mounted data drive the
		// output lives on and a RAM-backed /tmp — see utils.SpillTmpDir for the
		// order it picks a disk in.
		cmd = fmt.Sprintf(`gatk --java-options "%s" MarkDuplicates -R %s -I %s -O %s -M %s --TMP_DIR %s --VERBOSITY %s`,
			javaOpts, referencePath, sortedBam, rgmdBam, rgmdMetrics, utils.SpillTmpDir(rgmdBam), gatkLogLevel)

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
