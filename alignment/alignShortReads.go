package alignment

import (
	"fmt"
	"log"
	"os"
	"os/exec"
)

func AlignShortReadsMem(referencePath string, forwardPath string, reversePath string, sampleName string, libName string, outputDir string) {
	fmt.Println("Reading ...")
	fmt.Printf("ReferencePath: %s\nFwd: %s\nRev: %s\nSample Name: %s\nlib Name: %s\nOuDir: %s\n", referencePath, forwardPath, reversePath, sampleName, libName, outputDir)

	lineDir := fmt.Sprintf("%s/%s", outputDir, sampleName)
	bErr := os.MkdirAll(lineDir, 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
	}
	readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
	sortedBam := fmt.Sprintf("%s/%s.sorted.bam", lineDir, sampleName)
	rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sampleName)
	rgmdMetrics := fmt.Sprintf("%s/%s.RGMD.metrics.txt", lineDir, sampleName)
	rgmdIndex := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sampleName)

	// Construct the one-liner shell command
	cmdStr := fmt.Sprintf(`bwa mem -t 8 -M -Y -R '%s' %s %s %s | samtools sort -o %s`, readGroup, referencePath, forwardPath, reversePath, sortedBam)
	fmt.Println(cmdStr)
	cmd := exec.Command("bash", "-c", cmdStr)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	if err != nil {
		return
	}
	fmt.Printf("Marking duplicates ....")
	mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -I %s -O %s -M %s`, sortedBam, rgmdBam, rgmdMetrics)
	mDupCmd := exec.Command("bash", "-c", mDupCmdStr)
	mDupCmd.Stdout = os.Stdout
	mDupCmd.Stderr = os.Stderr

	mErr := mDupCmd.Run()
	if mErr != nil {
		return
	}

	fmt.Printf("Index Bam ....")
	indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
	indexCmd := exec.Command("bash", "-c", indexCmdStr)
	indexCmd.Stdout = os.Stdout
	indexCmd.Stderr = os.Stderr

	iErr := indexCmd.Run()
	if iErr != nil {
		return
	}

}
