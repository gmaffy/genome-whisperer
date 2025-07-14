package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"strings"
)

func AlignmentStats(ref, bam string) error {
	cmdStr := fmt.Sprintf(`gatk CollectAlignmentSummaryMetrics -R %s -I %s -O %s`, ref, bam, strings.TrimSuffix(bam, ".bam")+"_alignment_metrics.txt")

	cmdStr2 := fmt.Sprintf(`gatk CollectInsertSizeMetrics -R %s -I %s -O %s -H %s`, ref, bam, strings.TrimSuffix(bam, ".bam")+"_insert_metrics.txt", strings.TrimSuffix(bam, ".bam")+"_insert_size_histogram.pdf")

	cmdStr3 := fmt.Sprintf(`samtools flagstat %s > %s`, bam, strings.TrimSuffix(bam, ".bam")+"_flagstats.txt")

	fmt.Println(cmdStr)
	err := utils.RunBashCmdVerbose(cmdStr)
	if err != nil {
		fmt.Printf("Error running CollectAlignmentSummaryMetrics: %v\n", err)
		return err
	}

	fmt.Println(cmdStr2)
	if err2 := utils.RunBashCmdVerbose(cmdStr2); err2 != nil {
		log.Fatalf("Error running CollectInsertSizeMetrics: %v\n", err2)
		return err2
	}

	fmt.Println(cmdStr3)
	if err3 := utils.RunBashCmdVerbose(cmdStr3); err3 != nil {
		log.Fatalf("Error running samtools flagstat: %v\n", err3)
		return err3
	}
	return nil

}
