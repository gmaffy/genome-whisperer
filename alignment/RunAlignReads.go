package alignment

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
)

func RunAlignReads(referencePath string, forwardPath string, reversePath string, sePath string, sampleName string, libName string, outputDir string, threads int, aligner string, knownSites []string, bqsr bool, bootstrap bool, logFilePath string, preset string, gatkLogLevel string, verbose bool, outputFmt string) (string, error) {

	// ----------------------------------------- Log file ----------------------------------------------------------- //
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)

	jlog := slog.New(jsonHandler)
	logged := utils.ParseLogFile(logFilePath)

	// ----------------------------------------- Output Paths ------------------------------------------------------- //
	lineDir := fmt.Sprintf("%s/%s", outputDir, sampleName)
	rgBam := ""
	rgmdBam := ""
	rgmdIndex := ""
	sortedBam := ""
	sortedBai := ""
	if outputFmt == "bam" {

		rgmdBam = fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sampleName)
		rgmdIndex = fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sampleName)

	} else if outputFmt == "cram" {

		rgmdBam = fmt.Sprintf("%s/%s.RGMD.cram", lineDir, sampleName)
		rgmdIndex = fmt.Sprintf("%s/%s.RGMD.cram.crai", lineDir, sampleName)

	} else {
		return "", fmt.Errorf("invalid output format: %s", outputFmt)
	}
	rgBam = fmt.Sprintf("%s/%s.RG.bam", lineDir, sampleName)
	sortedBam = fmt.Sprintf("%s/%s.sorted.bam", lineDir, sampleName)
	sortedBai = fmt.Sprintf("%s/%s.sorted.bai", lineDir, sampleName)

	rgmdMetrics := fmt.Sprintf("%s/%s.RGMD.metrics.txt", lineDir, sampleName)

	finalBam := ""

	bErr := os.MkdirAll(lineDir, 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
		return "", bErr
	}

	// ========================== CREATING RGMD (depending on aligner ) STARTS ====================================== //

	if aligner == "bwa-mem2" {
		fmt.Printf("\n\n ----------------------------------------- Running aligner: bwa mem -------------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BWA_MEM", sampleName, "ALL") {
			msg := fmt.Sprintf("BWA_MEM already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")

			readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
			cmdStr := fmt.Sprintf(`bwa-mem2 mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
			fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

			var memErr error
			if verbose {
				memErr = utils.RunBashCmdVerbose(cmdStr)
			} else {
				memErr = utils.RunBashCmd(cmdStr)
			}
			if memErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr))
				return "", memErr
			}

			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- Mark Duplicates -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "MARK_DUPLICATES", sampleName, "ALL") {
			msg := fmt.Sprintf("MARK_DUPLICATES already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {

			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`, referencePath, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
			fmt.Printf("%s\n-----------------------------------------------\n\n", mDupCmdStr)

			var mdupErr error
			if verbose {
				mdupErr = utils.RunBashCmdVerbose(mDupCmdStr)
			} else {
				mdupErr = utils.RunBashCmd(mDupCmdStr)
			}
			if mdupErr != nil {
				jlog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr))
				slog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr))
				return "", mdupErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- BAM/CRAM INDEXING -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BAM_INDEX", sampleName, "ALL") {
			msg := fmt.Sprintf("BAM INDEXING already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			indexCmdStr := fmt.Sprintf(`samtools index %s`, rgmdBam)
			fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

			var indErr error
			if verbose {
				indErr = utils.RunBashCmdVerbose(indexCmdStr)
			} else {
				indErr = utils.RunBashCmd(indexCmdStr)
			}
			if indErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				return "", indErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		}

	} else if aligner == "bwa-mem" {
		fmt.Printf("\n\n ----------------------------------------- Running aligner: bwa mem -------------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BWA_MEM", sampleName, "ALL") {
			msg := fmt.Sprintf("BWA_MEM already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")

			readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", sampleName, sampleName, libName)
			cmdStr := fmt.Sprintf(`bwa mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, referencePath, forwardPath, reversePath, sortedBam)
			fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

			var memErr error
			if verbose {
				memErr = utils.RunBashCmdVerbose(cmdStr)
			} else {
				memErr = utils.RunBashCmd(cmdStr)
			}
			if memErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr))
				return "", memErr
			}

			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- Mark Duplicates -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "MARK_DUPLICATES", sampleName, "ALL") {
			msg := fmt.Sprintf("MARK_DUPLICATES already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {

			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -I %s -O %s -M %s --VERBOSITY %s`, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
			fmt.Printf("%s\n-----------------------------------------------\n\n", mDupCmdStr)

			var mdupErr error
			if verbose {
				mdupErr = utils.RunBashCmdVerbose(mDupCmdStr)
			} else {
				mdupErr = utils.RunBashCmd(mDupCmdStr)
			}
			if mdupErr != nil {
				jlog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr))
				slog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr))
				return "", mdupErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- BAM/CRAM INDEXING -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BAM_INDEX", sampleName, "ALL") {
			msg := fmt.Sprintf("BAM INDEXING already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
			fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

			var indErr error
			if verbose {
				indErr = utils.RunBashCmdVerbose(indexCmdStr)
			} else {
				indErr = utils.RunBashCmd(indexCmdStr)
			}
			if indErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				return "", indErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		}

	} else if aligner == "bowtie2" {
		fmt.Printf("\n\n ----------------------------------------- Running aligner: Bowtie2 -------------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BOWTIE2", sampleName, "ALL") {
			msg := fmt.Sprintf("BOWTIE2 already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			cmdStr := fmt.Sprintf(`bowtie2 -I 0 -X 1000 -x %s -1 %s -2 %s --end-to-end --sensitive --threads %v  --rg-id %s.1 --rg PL:BGISEQ --rg SM:%s --rg LB:%s | samtools sort -o %s`, referencePath, forwardPath, reversePath, threads, sampleName, sampleName, libName, sortedBam)
			fmt.Println(cmdStr)
			var bowErr error
			if verbose {
				bowErr = utils.RunBashCmdVerbose(cmdStr)
			} else {
				bowErr = utils.RunBashCmd(cmdStr)
			}
			if bowErr != nil {
				jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bowErr))
				slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bowErr))
				return "", bowErr
			}

			fmt.Printf("Indexing Bam ....")
			mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, sortedBam, sortedBai)
			var mErr error
			if verbose {
				mErr = utils.RunBashCmdVerbose(mDupCmdStr)
			} else {
				mErr = utils.RunBashCmd(mDupCmdStr)
			}
			if mErr != nil {
				jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mErr))
				slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mErr))
				return "", mErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

	} else if aligner == "pbmm2" {
		fmt.Printf("\n\n ----------------------------------------- Running aligner: PBMM2 -------------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "PBMM2", sampleName, "ALL") {
			msg := fmt.Sprintf("PBMM2 already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
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
					return "", indexErr
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
			if pbmm2Err != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", pbmm2Err))
				slog.Error("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", pbmm2Err))
				//fmt.Println("pbmm2 alignment failed")
				return "", pbmm2Err
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- Adding read groups -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "RG", sampleName, "ALL") {
			msg := fmt.Sprintf("RG already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			rgCmdStr := fmt.Sprintf(`gatk AddOrReplaceReadGroups -I %s -O %s -ID %s.1 -LB %s -PL PACBIO -PU BKD -SM %s`, sortedBam, rgBam, sampleName, libName, sampleName)
			fmt.Printf("%s\n-----------------------------------------------\n\n", rgCmdStr)
			var rgErr error
			if verbose {
				rgErr = utils.RunBashCmdVerbose(rgCmdStr)
			} else {
				rgErr = utils.RunBashCmd(rgCmdStr)
			}
			if rgErr != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", rgErr))
				slog.Error("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", rgErr))
				return "", rgErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		}

		fmt.Printf("\n\n------------------------------------------- Running pbmarkdup -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "PBMARKDUP", sampleName, "ALL") {
			msg := fmt.Sprintf("PBMARKDUP already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			mkdpCmdStr := fmt.Sprintf(`pbmm2 markdup %s %s`, rgBam, rgmdBam)
			fmt.Printf("%s\n-----------------------------------------------\n\n", mkdpCmdStr)
			var mkdpErr error
			if verbose {
				mkdpErr = utils.RunBashCmdVerbose(mkdpCmdStr)
			} else {
				mkdpErr = utils.RunBashCmd(mkdpCmdStr)
			}
			if mkdpErr != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mkdpErr))
				slog.Error("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mkdpErr))
				return "", mkdpErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- BAM INDEX -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BAM_INDEX", sampleName, "ALL") {
			msg := fmt.Sprintf("BAM INDEXING already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
			fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

			var indErr error
			if verbose {
				indErr = utils.RunBashCmdVerbose(indexCmdStr)
			} else {
				indErr = utils.RunBashCmd(indexCmdStr)
			}
			if indErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				slog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				return "", indErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

	} else {
		fmt.Println("Aligner not supported")
		return "", fmt.Errorf("aligner not supported\n\n Only bwa-mem, bowtie2 and pbmm2 are supported\n\n")
	}

	// ============================================== RGMD CREATION ENDS ============================================ //

	// ==================================== BSQR CREATION STARTS (IT TRUE) ========================================== //
	if bqsr {
		fmt.Println("Running BQSR")
		//fmt.Println("Check if mark duplicates was successful")
		//_, mdErr := os.Stat(rgmdBam)
		//_, mdiErr := os.Stat(rgmdIndex)
		//
		//if mdErr != nil || mdiErr != nil {
		//	fmt.Println("Mark duplicates failed")
		//	return "", fmt.Errorf("mark duplicates failed\n\n")
		//}
		//fmt.Println("Mark duplicates successful")

		if len(knownSites) == 0 && bootstrap == true {
			fmt.Println("Running with bootstrap method")
			fmt.Println("Create Known variants ..............")
			if utils.StageHasCompleted(logged, "BOOTSTRAP_CKV", sampleName, "ALL") {
				msg := fmt.Sprintf("BOOTSTRAP already completed for bam file: %s. Reusing known-sites.\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
				// Derive expected known-sites paths from the BAM name
				baseName := ""
				if strings.HasSuffix(rgmdBam, ".bam") {
					baseName = strings.TrimSuffix(rgmdBam, ".bam")
				} else {
					baseName = strings.TrimSuffix(rgmdBam, ".cram")
				}
				snpVCF := baseName + ".raw.SNP.hard_filtered.vcf.gz"
				indelVCF := baseName + ".raw.INDEL.hard_filtered.vcf.gz"
				// Use existing known-sites if present
				if _, err := os.Stat(snpVCF); err == nil {
					knownSites = append(knownSites, snpVCF)
				}
				if _, err := os.Stat(indelVCF); err == nil {
					knownSites = append(knownSites, indelVCF)
				}
				if len(knownSites) == 0 {
					// Fallback to recomputing if files are missing
					jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "RETRYING_CREATE_KNOWN_VARIANTS")
					newKnownSites, ksErr := CreateKnownVariants(referencePath, rgmdBam, logFilePath, gatkLogLevel, verbose)
					if ksErr != nil {
						jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
						slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
						return "", fmt.Errorf("create known variants failed\n\n")
					}
					knownSites = newKnownSites
				}
			} else {
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				newKnownSites, ksErr := CreateKnownVariants(referencePath, rgmdBam, logFilePath, gatkLogLevel, verbose)
				if ksErr != nil {
					jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
					slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
					return "", fmt.Errorf("create known variants failed\n\n")
				}
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				knownSites = newKnownSites
				fmt.Printf("Known variants created at %v\n...................................................\n", knownSites)

			}

			if utils.StageHasCompleted(logged, "BOOTSTRAP_BQSR", sampleName, "ALL") {
				msg := fmt.Sprintf("BOOTSTRAP_BQSR already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
			} else {
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				bqsrBam, recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath, verbose)
				if recErr != nil {
					jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					return "", fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
				}
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("COMPLETED - BQSR BAM: %s", bqsrBam))
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			}

		} else if len(knownSites) > 0 && bootstrap == false {
			if utils.StageHasCompleted(logged, "DBSNP_BQSR", sampleName, "ALL") {
				msg := fmt.Sprintf("DBSNP_BQSR already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
			} else {
				jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				bqsrBam, recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath, verbose)
				if recErr != nil {
					jlog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					slog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					return "", fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
				}
				jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("COMPLETED - BQSR BAM: %s", bqsrBam))
				jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				finalBam = bqsrBam

			}

		} else {
			finalBam = rgmdBam
		}

	}
	return finalBam, nil
}

func RunAlignReadsConfig(configPath string, threadsPerSample int, bqsr bool, bootstrap bool, aligner string, preset string, gatkLogLevel string, verbose bool, outputFmt string) ([]string, error) {

	fmt.Printf(" ---------------------------------------- Checking File Paths --------------------------------------------------------- \n\n")
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return nil, fmt.Errorf("error reading config: %v", err)
	}

	ref := cfg.Reference
	_, refErr := os.Stat(ref)
	if refErr != nil {
		fmt.Printf("Reference genome path: %s, is not valid\n", ref)
		return nil, fmt.Errorf("reference genome path: %s, is not valid", ref)
	}

	out := cfg.OutputDir
	outInfo, outErr := os.Stat(out)

	if outErr != nil {

		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory: %s does not exist. Attempting to create it.", out)
			if createErr := os.MkdirAll(out, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v", out, createErr)
				return nil, fmt.Errorf("failed to create output directory %s: %v", out, createErr)
			}
			fmt.Printf("Output directory %s created successfully.\n", out)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", out, outErr)
			return nil, fmt.Errorf("error accessing output directory %s: %v", out, outErr)
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output Directory %s file path is not a directory\n", out)
		return nil, fmt.Errorf("output directory %s file path is not a directory", out)
	}

	fmt.Println("Checking reference file index ...")

	pErr := utils.PrepareFasta(ref, aligner, verbose)

	if pErr != nil {
		log.Fatalf("Error preparing reference file: %v", pErr)
	}

	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(out, "alignment.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	fmt.Println("Log file created.")

	jsonHandler := slog.NewJSONHandler(logFile, nil)

	jlog := slog.New(jsonHandler)

	logged := utils.ParseLogFile(logFilePath)

	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Output directory:", cfg.OutputDir)

	// ----------------------------------------------- Check Paths if bqsr ------------------------------------------ //
	knownSites := cfg.KnownSites
	if bqsr {
		fmt.Println("----------------------------------------------- Check Paths if bqsr ------------------------------------------")
		if aligner == "pbmm2" {
			fmt.Println("We do not support BQSR for pbmm2 aligner. Please use bwa-mem or bowtie2 aligner or disable BQSR")
			return nil, fmt.Errorf("we do not support BQSR for pbmm2 aligner. Please use bwa-mem or bowtie2 aligner or disable BQSR")
		}
		if len(knownSites) == 0 && !bootstrap {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return nil, fmt.Errorf("either pass a known-sites file or enable bootstrap method")
		} else if len(knownSites) > 0 {
			fmt.Println("Running with known-sites flag")
			// ---------------------------- Checking Known sites file paths ----------------------------------------- //
			for j := range knownSites {
				_, err := os.Stat(knownSites[j])
				if err != nil {
					fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
					return nil, fmt.Errorf("known-sites file: %s is not a valid file path", knownSites[j])
				}
			}
			if bootstrap {
				fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
				return nil, fmt.Errorf("choose either pass a known-sites file or enable bootstrap method, but not both")
			}
		}
	}

	//----------------------------------------- Check reads based on SE OR PE ----------------------------------------//

	i := 0
	if aligner == "bwa-mem" || aligner == "bowtie2" || aligner == "bwa-mem2" {
		fmt.Printf("\n---------------------------------------------Checking Read pairs---------------------------------------------------\n\n")
		for _, pair := range cfg.ReadPairs {
			if len(pair) < 4 {
				fmt.Printf("This read pair is wrongly formated %s\n", pair)
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				continue
			}

			fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]

			_, fwdErr := os.Stat(fwd)
			_, revErr := os.Stat(rev)

			if fwdErr != nil {
				fmt.Printf("Forward reads path %s, is not valid\n", fwd)
				return nil, fmt.Errorf("forward reads path %s, is not valid\n", fwd)
			}

			if revErr != nil {
				fmt.Printf("Reverse reads path %s, is not valid\n", rev)
				return nil, fmt.Errorf("reverse reads path %s, is not valid\n", rev)
			}

			if sn == "" {
				fmt.Println("Please provide sample name ")
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				return nil, fmt.Errorf("please provide sample name. Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			}
			if lb == "" {
				fmt.Println("Please provide library name  ")
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				return nil, fmt.Errorf("please provide library name. Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
			}
			i++
		}
		fmt.Printf("There are  %v read pairs", i)
	}

	if aligner == "pbmm2" {
		fmt.Printf("\n---------------------------------------------Checking SingleRead pairs---------------------------------------------------\n\n")
		for _, se := range cfg.SeReads {
			if len(se) < 3 {
				fmt.Printf("This read pair is wrongly formated %s\n", se)
				fmt.Println("Supply reads in this format: SingleRead: <reads> <sample name> <library name> ")
				continue
			}
			seRead, sn, lb := se[0], se[1], se[2]
			_, seReadErr := os.Stat(seRead)
			if seReadErr != nil {
				fmt.Printf("Single reads path %s, is not valid\n", seRead)
				return nil, fmt.Errorf("single reads path %s, is not valid\n", seRead)
			}
			if sn == "" {
				fmt.Println("Please provide sample name ")
				return nil, fmt.Errorf("please provide sample name. Supply reads in this format: SingleRead: <reads> <sample name> <library name> ")
			}
			if lb == "" {
				fmt.Println("Please provide library name  ")
				return nil, fmt.Errorf("please provide library name. Supply reads in this format: SingleRead: <reads> <sample name> <library name> ")
			}
			i++
		}
		fmt.Printf("There are  %v reads", i)
	}

	// ============================================== Run Alignments ================================================ //

	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threadsPerSample = totalCores
	}

	fmt.Printf("Running up to %d jobs in parallel with %d threads each\n", maxParallelJobs, threadsPerSample)

	var wg sync.WaitGroup
	sem := make(chan struct{}, maxParallelJobs)

	//------------------------------------------ Run alignment ------------------------------------------------------ //

	// For each sample, this outputs bam. Either bqsr bam or rgmd bams
	var bams []string
	if aligner == "bwa-mem" || aligner == "bowtie2" || aligner == "bwa-mem2" {

		// ----------------------------------------- RUNNING STATS HERE --------------------------------------------- //
		for _, pair := range cfg.ReadPairs {
			if len(pair) < 4 {
				fmt.Printf("This read pair is wrongly formated %s\n", pair)
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name>")
				continue
			}

			wg.Add(1)
			sem <- struct{}{}
			go func(pair []string) {
				defer wg.Done()
				defer func() { <-sem }()

				fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]

				jlog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", "STARTED")

				isDone := utils.StageHasCompleted(logged, "PE_ALIGNMENT", sn, "ALL")
				if isDone {
					msg := fmt.Sprintf("%s and MarkDuplicates already completed for %s. Skipping.\n\n-------------------------------------------------------\n\n", aligner, sn)
					slog.Info(msg)

				} else {
					bam, alErr := RunAlignReads(cfg.Reference, fwd, rev, "", sn, lb, cfg.OutputDir, threadsPerSample, aligner, knownSites, bqsr, bootstrap, logFilePath, preset, gatkLogLevel, verbose, outputFmt)
					if alErr != nil {
						jlog.Error("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", alErr))
						slog.Error("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED - %v", alErr))
						return
					}
					jlog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
					slog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", "COMPLETED")
					bams = append(bams, bam)
				}
			}(pair)

		}
		wg.Wait()

	} else if aligner == "pbmm2" {

		for _, se := range cfg.SeReads {
			if len(se) < 3 {
				fmt.Printf("This read pair is wrongly formated %s\n", se)
				fmt.Println("Supply reads in this format: SingleRead: <reads> <sample name> <library name>")
				continue
			}
			wg.Add(1)
			sem <- struct{}{}
			go func(pair []string) {
				defer wg.Done()
				defer func() { <-sem }()
				seRead, sn, lb := se[0], se[1], se[2]

				jlog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "STATUS", "STARTED")

				isDone := utils.StageHasCompleted(logged, "PB_ALIGNMENT", sn, "ALL")
				if isDone {
					msg := fmt.Sprintf("%s and MarkDuplicates already completed for %s. Skipping.\n\n-------------------------------------------------------\n\n", aligner, sn)
					slog.Info(msg)

				} else {
					bam, alErr := RunAlignReads(cfg.Reference, "", "", seRead, sn, lb, cfg.OutputDir, threadsPerSample, aligner, knownSites, bqsr, bootstrap, logFilePath, preset, gatkLogLevel, verbose, outputFmt)
					if alErr != nil {
						jlog.Error("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", alErr))
						slog.Error("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED - %v", alErr))

						return
					}
					jlog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
					slog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "STATUS", "COMPLETED")
					bams = append(bams, bam)
				}

			}(se)
		}
		wg.Wait()

	} else {
		Msg := fmt.Sprintf("Aligner %s is not supported", aligner)
		slog.Error(Msg)
		return nil, fmt.Errorf(Msg)
	}
	return bams, nil

}

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

func RunAlignReadsDir(dataDir string, species string, refVer string, refFasta string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, threads int) {
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
	if refFasta == "" {
		fmt.Println("Please provide reference name")
		return
	}

	fastaInfo, err := os.Stat(refFasta)
	if err != nil {
		fmt.Printf("Error accessing reference fasta file: %s\n", refFasta)
		return
	}
	if !fastaInfo.Mode().IsRegular() {
		fmt.Printf("Reference fasta file: %s is not a regular file\n", refFasta)
		return
	}

	dictFilePath := refFasta[:len(refFasta)-len(filepath.Ext(refFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
		return
	}

	color.Green("All file paths valid\n....................................................\n\n")
	fmt.Println("data dir abs:", dataDirAbs)
	// ================================== Discover samples ===================================== //
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

	// ======================================= Resolve valid samples up-front ======================================= //
	fmt.Println("Looking for bam files ....")

	var (
		validSamples []sampleWork
		missingBams  []string
		multipleBams []string
	)

	for _, sample := range samples {
		_, cramFiles := variants.FindBamOrVcfs(dataDirAbs, species, sample, refVer, "bams", "")
		switch len(cramFiles) {
		case 0:
			missingBams = append(missingBams, sample)
			color.Red("[%s] bqsr.cram file MISSING! 😒\n", sample)
		case 1:
			cram := cramFiles[0]
			//fmt.Printf("CRAM file found for %s: %s\n", sample, cram)
			color.Green("[%s] bqsr.cram file FOUND! 😊: %s\n", sample, cram)
			err := variants.ValidateBam(cram, refFasta, verbose, quick)
			if err != nil {
				color.Red("[%s] bqsr.cram file is not valid: %v\n", sample, err)
				missingBams = append(missingBams, sample)
			} else {

				validSamples = append(validSamples, sampleWork{sample: sample, cram: cram, cramDir: filepath.Dir(filepath.Dir(filepath.Dir(cram)))})
			}

		default:
			color.Red("[%s] has multiple cram files.⚠️.\n%s\n\nI don't know which one to use.  Remove bad ones: %v\n\n", sample, cramFiles)
			multipleBams = append(multipleBams, sample)
		}
	}

	color.Green("\nValid: %d\n", len(validSamples))
	color.Red("Missing: %d\n", len(missingBams))
	color.Yellow("Multiple: %d\n", len(multipleBams))

	fmt.Printf("\n==============================================================\n\n")

	type samplePair struct {
		sample string
		fwd    string
		rev    string
		bamDir string
	}

	var missingCleanReads []string
	var multipleCleanReads []string
	var validCleanReads []samplePair

	if len(missingBams) > 0 {
		color.Green("Aligning reads for %d samples:\n---------------------------------------------------------------------\n\n ", len(missingBams))
		for _, cleanReadsDir := range missingBams {
			sample := filepath.Base(filepath.Dir(cleanReadsDir))
			fmt.Println(sample)
			if strings.HasSuffix(sample, "LR") {
				fmt.Println("This is LR sample. Using pbmm2")
			} else {
				fmt.Println("This is not LR sample. Using bwa-mem2")
				fmt.Println("Check if reference_genomes dir is present in the data directory")
				refDir := filepath.Join(filepath.Dir(cleanReadsDir), "reference_genomes")
				fmt.Println("refDir:", refDir)
				_, err := os.Stat(refDir)
				if err != nil {
					err := os.MkdirAll(refDir, 0755)
					if err != nil {
						color.Red("Error creating reference_genomes dir in the data directory: %v\n", err)
						missingCleanReads = append(missingCleanReads, sample)
					}
				} else {
					fmt.Println("Check if reference_genomes species bams dir is present in the data directory")
					verDir := filepath.Join(refDir, refVer, "bams")
					_, err := os.Stat(verDir)
					if err != nil {
						err := os.MkdirAll(verDir, 0755)
						if err != nil {
							color.Red("Error creating reference_genomes bams dir in the data directory: %v\n", err)
							missingCleanReads = append(missingCleanReads, sample)
						}
					} else {
						fmt.Println("Check for PE reads ..")

						fwdReads, revReads, err := GetReadsPE(cleanReadsDir)
						fmt.Println("fwd reads:", fwdReads)
						fmt.Println("rev reads:", revReads)

						if err != nil {
							Msg := fmt.Sprintf("Error getting reads from %s: %v", cleanReadsDir, err)
							color.Red(Msg)
							missingCleanReads = append(missingCleanReads, sample)

						} else {
							if len(fwdReads) == 0 && len(revReads) == 0 {
								Msg := fmt.Sprintf("No reads found in %s", cleanReadsDir)
								color.Red(Msg)
								missingCleanReads = append(missingCleanReads, sample)
							} else if len(fwdReads) == 1 && len(revReads) == 1 {

								if skipVer {
									color.Yellow("Skip fastq verification for %s\n", sample)
									validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir})
									color.Green("[%s] PE reads found", sample)
								} else {
									color.Blue("Validating fwd PE read: %s\n", fwdReads[0])
									fwdErr := validateFastqGz(fwdReads[0], verbose, quick)
									revErr := validateFastqGz(revReads[0], verbose, quick)

									if fwdErr != nil || revErr != nil {
										color.Red("Error validating fwd PE read: %v\n", err)
										missingCleanReads = append(missingCleanReads, sample)
									} else {
										color.Green("Valid PE reads found")
										validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir})
									}

								}
							} else {
								multipleCleanReads = append(multipleCleanReads, cleanReadsDir)
								Msg := fmt.Sprintf("Multiple PE reads found in %s. Can not align", cleanReadsDir)
								color.Yellow(Msg)
								missingCleanReads = append(missingCleanReads, sample)

							}
						}

					}

				}
			}
		}

	}

	var failedAlignments []string
	color.Green("\nValid: %d\n", len(validCleanReads))
	if len(validCleanReads) > 0 {
		for _, readPair := range validCleanReads {
			switch aligner {
			case "bwa-mem":
				color.Green("Aligning reads for %s to reference genome %s with bwa-mem\n", readPair.sample, refFasta)

			default:

				color.Green("Aligning reads for %s to reference genome %s. with bwa-mem2\n", readPair.sample, refFasta)
				sortedBam := filepath.Join(readPair.bamDir, readPair.sample+".sorted.bam")
				readGroup := fmt.Sprintf("@RG\\tID:%s.1\\tSM:%s\\tLB:%s\\tPL:BGISEQ", readPair.sample, readPair.sample, fmt.Sprintf("%s_1", readPair.sample))
				cmdStr := fmt.Sprintf(`bwa-mem2 mem -t %v -M -Y -R '%s' %s %s %s | samtools sort -o %s`, threads, readGroup, refFasta, readPair.fwd, readPair.rev, sortedBam)
				fmt.Printf("%s\n--------------------------------------------\n\n", cmdStr)

				var memErr error
				if verbose {
					memErr = utils.RunBashCmdVerbose(cmdStr)
				} else {
					memErr = utils.RunBashCmd(cmdStr)
				}
				if memErr != nil {
					failedAlignments = append(failedAlignments, readPair.sample)
				} else {
					color.Green("Aligned reads for %s to reference genome %s. with bwa-mem2\n", readPair.sample, refFasta)
					rgmdBam := filepath.Join(readPair.bamDir, readPair.sample+".rgmd.bam")
					rgmdMetrics := filepath.Join(readPair.bamDir, readPair.sample+".rgmd.metrics.txt")
					mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`, refFasta, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
					fmt.Printf("%s\n-----------------------------------------------\n\n", mDupCmdStr)

					var mdupErr error
					if verbose {
						mdupErr = utils.RunBashCmdVerbose(mDupCmdStr)
					} else {
						mdupErr = utils.RunBashCmd(mDupCmdStr)
					}
					if mdupErr != nil {
						failedAlignments = append(failedAlignments, readPair.sample)

					} else {
						indexCmdStr := fmt.Sprintf(`samtools index %s`, rgmdBam)
						fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

						var indErr error
						if verbose {
							indErr = utils.RunBashCmdVerbose(indexCmdStr)
						} else {
							indErr = utils.RunBashCmd(indexCmdStr)
						}
						if indErr != nil {
							failedAlignments = append(failedAlignments, readPair.sample)
						} else {
							if bqsr {
								fmt.Println("Running BQSR")
								//fmt.Println("Check if mark duplicates was successful")
								//_, mdErr := os.Stat(rgmdBam)
								//_, mdiErr := os.Stat(rgmdIndex)
								//
								//if mdErr != nil || mdiErr != nil {
								//	fmt.Println("Mark duplicates failed")
								//	return "", fmt.Errorf("mark duplicates failed\n\n")
								//}
								//fmt.Println("Mark duplicates successful")

								if len(knownSites) == 0 && bootstrap == true {
									fmt.Println("Running with bootstrap method")
									fmt.Println("Create Known variants ..............")
									if utils.StageHasCompleted(logged, "BOOTSTRAP_CKV", sampleName, "ALL") {
										msg := fmt.Sprintf("BOOTSTRAP already completed for bam file: %s. Reusing known-sites.\n\n------------------------------------------------------------------------\n\n", sampleName)
										slog.Info(msg)
										// Derive expected known-sites paths from the BAM name
										baseName := ""
										if strings.HasSuffix(rgmdBam, ".bam") {
											baseName = strings.TrimSuffix(rgmdBam, ".bam")
										} else {
											baseName = strings.TrimSuffix(rgmdBam, ".cram")
										}
										snpVCF := baseName + ".raw.SNP.hard_filtered.vcf.gz"
										indelVCF := baseName + ".raw.INDEL.hard_filtered.vcf.gz"
										// Use existing known-sites if present
										if _, err := os.Stat(snpVCF); err == nil {
											knownSites = append(knownSites, snpVCF)
										}
										if _, err := os.Stat(indelVCF); err == nil {
											knownSites = append(knownSites, indelVCF)
										}
										if len(knownSites) == 0 {
											// Fallback to recomputing if files are missing
											jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "RETRYING_CREATE_KNOWN_VARIANTS")
											newKnownSites, ksErr := CreateKnownVariants(referencePath, rgmdBam, logFilePath, gatkLogLevel, verbose)
											if ksErr != nil {
												jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
												slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
												return "", fmt.Errorf("create known variants failed\n\n")
											}
											knownSites = newKnownSites
										}
									} else {
										jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										newKnownSites, ksErr := CreateKnownVariants(referencePath, rgmdBam, logFilePath, gatkLogLevel, verbose)
										if ksErr != nil {
											jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
											slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
											return "", fmt.Errorf("create known variants failed\n\n")
										}
										jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
										slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CKV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
										knownSites = newKnownSites
										fmt.Printf("Known variants created at %v\n...................................................\n", knownSites)

									}

									if utils.StageHasCompleted(logged, "BOOTSTRAP_BQSR", sampleName, "ALL") {
										msg := fmt.Sprintf("BOOTSTRAP_BQSR already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
										slog.Info(msg)
									} else {
										jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										bqsrBam, recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath, verbose)
										if recErr != nil {
											jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
											slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
											return "", fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
										}
										jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("COMPLETED - BQSR BAM: %s", bqsrBam))
										jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
										slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
									}

								} else if len(knownSites) > 0 && bootstrap == false {
									if utils.StageHasCompleted(logged, "DBSNP_BQSR", sampleName, "ALL") {
										msg := fmt.Sprintf("DBSNP_BQSR already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
										slog.Info(msg)
									} else {
										jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										slog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
										bqsrBam, recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath, verbose)
										if recErr != nil {
											jlog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
											slog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
											return "", fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
										}
										jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("COMPLETED - BQSR BAM: %s", bqsrBam))
										jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
										slog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
										finalBam = bqsrBam

									}

								} else {
									finalBam = rgmdBam
								}

							}
							return finalBam, nil
						}
					}
				}

			}

		}
	}
}
