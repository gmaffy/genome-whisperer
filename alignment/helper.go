package alignment

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

func BamToCram(bamPath, refFasta string, verbose bool) error {
	cramPath := strings.TrimSuffix(bamPath, filepath.Ext(bamPath)) + ".cram"
	bamToCramStr := fmt.Sprintf(`samtools view -T %s -C -o %s %s`, refFasta, cramPath, bamPath)
	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(bamToCramStr)
	} else {
		err = utils.RunBashCmd(bamToCramStr)
	}
	return err
}

func BwaMem2Align(forwardPath string, reversePath string, referencePath string, sampleName string, libName string, threads int, sortedBam string, verbose bool) error {

	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	cmdStr := fmt.Sprintf(`bwa-mem2 mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
	fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

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
	cmdStr := fmt.Sprintf(`bwa mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
	fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	return err
}

func Bowtie2Align(forwardPath string, reversePath string, referencePath string, sortedBam string, sampleName string, libName string, threads int, verbose bool) error {
	cmdStr := fmt.Sprintf(`bowtie2 -I 0 -X 1000 -x %s -1 %s -2 %s --end-to-end --sensitive --threads %v  --rg-id %s.1 --rg PL:BGISEQ --rg SM:%s --rg LB:%s | samtools sort -o %s`, referencePath, forwardPath, reversePath, threads, sampleName, sampleName, libName, sortedBam)
	fmt.Println(cmdStr)
	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(cmdStr)
	} else {
		err = utils.RunBashCmd(cmdStr)
	}
	return err
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
	fmt.Println(pbmm2CmdStr)
	var pbmm2Err error
	if verbose {
		pbmm2Err = utils.RunBashCmdVerbose(pbmm2CmdStr)
	} else {
		pbmm2Err = utils.RunBashCmd(pbmm2CmdStr)
	}

	return pbmm2Err
}

func ReadGroups(sortedBam string, rgBam string, sampleName string, libName string, verbose bool) error {
	rgCmdStr := fmt.Sprintf(`gatk AddOrReplaceReadGroups -I %s -O %s -ID %s.1 -LB %s -PL PACBIO -PU BKD -SM %s`, sortedBam, rgBam, sampleName, libName, sampleName)
	fmt.Printf("%s\n-----------------------------------------------\n\n", rgCmdStr)
	var rgErr error
	if verbose {
		rgErr = utils.RunBashCmdVerbose(rgCmdStr)
	} else {
		rgErr = utils.RunBashCmd(rgCmdStr)
	}
	return rgErr
}

func BamIndex(bam string, verbose bool) error {
	indexCmdStr := fmt.Sprintf(`samtools index %s`, bam)
	fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

	var err error
	if verbose {
		err = utils.RunBashCmdVerbose(indexCmdStr)
	} else {
		err = utils.RunBashCmd(indexCmdStr)
	}
	return err
}

func MarkDuplicates(referencePath string, sortedBam string, verbose bool, aligner string, gatkLogLevel string) error {

	baseName := ""
	if strings.HasSuffix(sortedBam, ".bam") {
		baseName = strings.TrimSuffix(sortedBam, ".bam")
	} else {
		baseName = strings.TrimSuffix(sortedBam, ".cram")
	}
	rgmdBam := baseName + ".RGMD.bam"
	rgmdMetrics := baseName + ".RGMD.metrics.txt"

	var cmd string
	switch aligner {
	case "bwa-mem":
		cmd = fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`, referencePath, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)

	case "bwa-mem2":
		cmd = fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`, referencePath, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)

	case "bowtie2":
		cmd = fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`, referencePath, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)

	case "pbmm2":
		cmd = fmt.Sprintf(`pbmm2 markdup %s %s`, sortedBam, rgmdBam)
	}

	fmt.Printf("%s\n-----------------------------------------------\n\n", cmd)

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
