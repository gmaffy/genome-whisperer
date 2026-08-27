package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"strings"
)

func Stats(ref, bam string, verbose bool) error {
	cmdStr := fmt.Sprintf(`gatk CollectAlignmentSummaryMetrics -R %s -I %s -O %s --tmp-dir %s`, ref, bam, strings.TrimSuffix(bam, ".bam")+"_alignment_metrics.txt", WorkTmpDir(bam))

	cmdStr2 := fmt.Sprintf(`gatk CollectInsertSizeMetrics -R %s -I %s -O %s -H %s --tmp-dir %s`, ref, bam, strings.TrimSuffix(bam, ".bam")+"_insert_metrics.txt", strings.TrimSuffix(bam, ".bam")+"_insert_size_histogram.pdf", WorkTmpDir(bam))

	cmdStr3 := fmt.Sprintf(`samtools flagstat %s > %s`, bam, strings.TrimSuffix(bam, ".bam")+"_flagstats.txt")

	if verbose {
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
	} else {
		err := utils.RunBashCmd(cmdStr)
		if err != nil {
			fmt.Printf("Error running CollectAlignmentSummaryMetrics: %v\n", err)
			return err
		}

		if err2 := utils.RunBashCmd(cmdStr2); err2 != nil {
			log.Fatalf("Error running CollectInsertSizeMetrics: %v\n", err2)
			return err2
		}

		if err3 := utils.RunBashCmd(cmdStr3); err3 != nil {
			log.Fatalf("Error running samtools flagstat: %v\n", err3)
			return err3
		}
	}
	return nil

}
