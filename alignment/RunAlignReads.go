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

// resumeStage records where in the pipeline a sample should resume.
type resumeStage int

const (
	resumeFromAlign            resumeStage = iota // no intermediate files found – start from scratch
	resumeFromMarkDupAndBQSR                      // sorted.bam exists – run MarkDuplicates → BQSR → convert
	resumeFromBQSR                                // rgmd.cram exists – run BQSR directly
	resumeFromCramRgmdThenBQSR                    // rgmd.bam exists – convert to cram, then BQSR
	resumeFromBamToCram                           // bqsr.bam exists – convert BAM→CRAM only
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

// GenomeRef holds a discovered reference genome on disk.
type GenomeRef struct {
	RefVer    string // e.g. "GRCh38"
	FastaPath string // path to .fa / .fasta / .fna
	DictPath  string // path to .dict
}

// GetValidGenomesFromDisk scans genomesDir for species → reference assemblies.
// Expected layout: <genomesDir>/<species>/<refVer>/assembly/  containing a fasta + dict file.
// Returns map[SPECIES_UPPER] → []GenomeRef.
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

// resolveRefFasta returns refFasta as-is when provided, otherwise auto-discovers
// it from genomesDir using species and refVer. Returns the fasta path or an error.
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

func RunAlignReadsDir1(dataDir string, species string, refVer string, refFasta string, genomesDir string, verbose bool, gatkLogLevel string, aligner string, quick bool, skipVer bool, bqsr bool, bootstrap bool, knownSites []string, threads int) {
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
	fmt.Println("Looking for existing bam/cram files ....")

	var (
		validSamples []sampleWork
		missingBams  []string
		multipleBams []string
	)

	for _, sample := range samples {
		sampleName := filepath.Base(filepath.Dir(sample))
		sampleBaseDir := filepath.Dir(sample)
		verDir := filepath.Join(sampleBaseDir, "reference_genomes", refVer, "bams")

		// Helper: check if a file exists in verDir.
		fileExists := func(name string) (string, bool) {
			p := filepath.Join(verDir, name)
			if _, err := os.Stat(p); err == nil {
				return p, true
			}
			return "", false
		}

		bqsrCram := sampleName + ".bqsr.cram"
		bqsrBam := sampleName + ".bqsr.bam"
		rgmdCram := sampleName + ".rgmd.cram"
		rgmdBamFile := sampleName + ".rgmd.bam"
		sortedBamFile := sampleName + ".sorted.bam"

		// ---- Stage 1: bqsr.cram already present ---- //
		if cramPath, ok := fileExists(bqsrCram); ok {
			color.Green("[%s] bqsr.cram FOUND 😊: %s\n", sampleName, cramPath)
			if err := variants.ValidateBam(cramPath, resolvedFasta, verbose, quick); err != nil {
				color.Red("[%s] bqsr.cram not valid: %v — will re-process\n", sampleName, err)
				missingBams = append(missingBams, sample)
			} else {
				validSamples = append(validSamples, sampleWork{
					sample:  sample,
					cram:    cramPath,
					cramDir: filepath.Dir(filepath.Dir(filepath.Dir(cramPath))),
				})
			}
			continue
		}

		// ---- Stage 2: bqsr.bam present → convert BAM→CRAM only ---- //
		if bamPath, ok := fileExists(bqsrBam); ok {
			color.Yellow("[%s] bqsr.bam found – will convert BAM→CRAM only.\n", sampleName)
			missingBams = append(missingBams, sample)
			// Store intermediate path; the goroutine will pick it up via samplePair.
			// We encode this by appending a special marker; handled below in the work list.
			_ = bamPath // used when building samplePair further down
			continue
		}

		// ---- Stage 3: rgmd.cram present → run BQSR on it ---- //
		if cramPath, ok := fileExists(rgmdCram); ok {
			color.Yellow("[%s] rgmd.cram found – will run BQSR.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = cramPath
			continue
		}

		// ---- Stage 4: rgmd.bam present → convert to CRAM, then BQSR ---- //
		if bamPath, ok := fileExists(rgmdBamFile); ok {
			color.Yellow("[%s] rgmd.bam found – will convert to CRAM then run BQSR.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = bamPath
			continue
		}

		// ---- Stage 5: sorted.bam present → run MarkDuplicates + BQSR + convert ---- //
		if bamPath, ok := fileExists(sortedBamFile); ok {
			color.Yellow("[%s] sorted.bam found – will run MarkDuplicates → BQSR → convert.\n", sampleName)
			missingBams = append(missingBams, sample)
			_ = bamPath
			continue
		}

		// ---- No intermediate files found – needs full alignment ---- //
		color.Red("[%s] No bqsr.cram or intermediate files found 😒 – needs full alignment.\n", sampleName)
		missingBams = append(missingBams, sample)
	}

	color.Green("\nAlready aligned (valid): %d\n", len(validSamples))
	color.Red("Need alignment (missing): %d\n", len(missingBams))
	color.Yellow("Skipped (multiple crам files): %d\n", len(multipleBams))
	fmt.Printf("\n==============================================================\n\n")

	// ======================================= Build work list for missing samples ================================== //
	type samplePair struct {
		sample          string
		fwd             string
		rev             string
		bamDir          string
		resume          resumeStage // where to start in the pipeline
		intermediateBam string      // pre-existing BAM/CRAM to resume from (if any)
	}

	var missingCleanReads []string
	var multipleCleanReads []string
	var validCleanReads []samplePair

	if len(missingBams) > 0 {
		color.Green("Preparing to align %d sample(s):\n---------------------------------------------------------------------\n\n", len(missingBams))
		for _, cleanReadsDir := range missingBams {
			sample := filepath.Base(filepath.Dir(cleanReadsDir))
			fmt.Println(sample)

			if strings.HasSuffix(sample, "LR") {
				// Long-read samples: pbmm2 path — not handled in this function yet
				color.Yellow("[%s] LR sample detected (pbmm2). Skipping in dir-mode for now.\n", sample)
				missingCleanReads = append(missingCleanReads, sample)
				continue
			}

			// Ensure output bams dir exists: <sampleDir>/reference_genomes/<refVer>/bams/
			sampleBaseDir := filepath.Dir(cleanReadsDir)
			verDir := filepath.Join(sampleBaseDir, "reference_genomes", refVer, "bams")
			if mkErr := os.MkdirAll(verDir, 0755); mkErr != nil {
				color.Red("[%s] Error creating bams dir %s: %v\n", sample, verDir, mkErr)
				missingCleanReads = append(missingCleanReads, sample)
				continue
			}

			// ---- Determine resume stage by re-checking intermediate files ---- //
			fileInVerDir := func(name string) (string, bool) {
				p := filepath.Join(verDir, name)
				if _, err := os.Stat(p); err == nil {
					return p, true
				}
				return "", false
			}

			bqsrBamPath, hasBqsrBam := fileInVerDir(sample + ".bqsr.bam")
			rgmdCramPath, hasRgmdCram := fileInVerDir(sample + ".rgmd.cram")
			rgmdBamPath, hasRgmdBam := fileInVerDir(sample + ".rgmd.bam")
			sortedBamPath, hasSortedBam := fileInVerDir(sample + ".sorted.bam")

			switch {
			case hasBqsrBam:
				color.Yellow("[%s] Resuming: converting bqsr.bam → bqsr.cram\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromBamToCram,
					intermediateBam: bqsrBamPath,
				})

			case hasRgmdCram:
				color.Yellow("[%s] Resuming: running BQSR on rgmd.cram\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromBQSR,
					intermediateBam: rgmdCramPath,
				})

			case hasRgmdBam:
				color.Yellow("[%s] Resuming: converting rgmd.bam → cram, then BQSR\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromCramRgmdThenBQSR,
					intermediateBam: rgmdBamPath,
				})

			case hasSortedBam:
				color.Yellow("[%s] Resuming: sorted.bam found – MarkDuplicates → BQSR → convert\n", sample)
				validCleanReads = append(validCleanReads, samplePair{
					sample:          sample,
					bamDir:          verDir,
					resume:          resumeFromMarkDupAndBQSR,
					intermediateBam: sortedBamPath,
				})

			default:
				// No intermediate files – need reads for full alignment
				fwdReads, revReads, readsErr := GetReadsPE(cleanReadsDir)
				fmt.Printf("[%s] fwd: %v  rev: %v\n", sample, fwdReads, revReads)
				if readsErr != nil {
					color.Red("[%s] Error reading clean_reads dir: %v\n", sample, readsErr)
					missingCleanReads = append(missingCleanReads, sample)
					continue
				}

				switch {
				case len(fwdReads) == 0 && len(revReads) == 0:
					color.Red("[%s] No PE reads found in %s\n", sample, cleanReadsDir)
					missingCleanReads = append(missingCleanReads, sample)
				case len(fwdReads) == 1 && len(revReads) == 1:
					if skipVer {
						color.Yellow("[%s] Skipping fastq verification.\n", sample)
						validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir, resume: resumeFromAlign})
					} else {
						color.Blue("[%s] Validating fastq files...\n", sample)
						fwdErr := validateFastqGz(fwdReads[0], verbose, quick)
						revErr := validateFastqGz(revReads[0], verbose, quick)
						if fwdErr != nil || revErr != nil {
							color.Red("[%s] Fastq validation failed (fwd=%v rev=%v)\n", sample, fwdErr, revErr)
							missingCleanReads = append(missingCleanReads, sample)
						} else {
							color.Green("[%s] PE reads valid.\n", sample)
							validCleanReads = append(validCleanReads, samplePair{sample: sample, fwd: fwdReads[0], rev: revReads[0], bamDir: verDir, resume: resumeFromAlign})
						}
					}
				default:
					color.Yellow("[%s] Multiple PE reads found — cannot determine which to use.\n", sample)
					multipleCleanReads = append(multipleCleanReads, cleanReadsDir)
					missingCleanReads = append(missingCleanReads, sample)
				}
			}
		}
	}

	color.Green("\nReady to align: %d\n", len(validCleanReads))
	color.Red("Missing/invalid reads: %d\n", len(missingCleanReads))
	color.Yellow("Skipped (multiple read files): %d\n", len(multipleCleanReads))
	fmt.Printf("\n==============================================================\n\n")

	if len(validCleanReads) == 0 {
		color.Yellow("No samples to align. Exiting.\n")
		return
	}

	// ======================================= Parallel alignment ================================================== //
	totalCores := runtime.NumCPU()
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}
	color.Green("Aligning %d sample(s) — up to %d jobs in parallel, %d threads each.\n\n",
		len(validCleanReads), maxParallelJobs, threads)

	var (
		mu               sync.Mutex
		failedAlignments []string
		wg               sync.WaitGroup
		sem              = make(chan struct{}, maxParallelJobs)
	)

	for _, rp := range validCleanReads {
		wg.Add(1)
		sem <- struct{}{}

		go func(readPair samplePair) {
			defer wg.Done()
			defer func() { <-sem }()

			sampleName := readPair.sample

			addFailure := func(reason string) {
				color.Red("[%s] FAILED: %s\n", sampleName, reason)
				mu.Lock()
				failedAlignments = append(failedAlignments, sampleName)
				mu.Unlock()
			}

			// ---- Paths ---- //
			sortedBam := filepath.Join(readPair.bamDir, sampleName+".sorted.bam")
			rgmdBam := filepath.Join(readPair.bamDir, sampleName+".rgmd.bam")
			rgmdCram := filepath.Join(readPair.bamDir, sampleName+".rgmd.cram")
			rgmdMetrics := filepath.Join(readPair.bamDir, sampleName+".rgmd.metrics.txt")
			bqsrCram := filepath.Join(readPair.bamDir, sampleName+".bqsr.cram")

			// bamToCram converts src BAM/CRAM to dst CRAM using samtools.
			bamToCram := func(src, dst string) error {
				cmd := fmt.Sprintf(`samtools view -C -T %s -o %s %s && samtools index %s`, resolvedFasta, dst, src, dst)
				color.Green("[%s] Converting to CRAM: %s\n", sampleName, cmd)
				if verbose {
					return utils.RunBashCmdVerbose(cmd)
				}
				return utils.RunBashCmd(cmd)
			}

			// runBQSR runs bootstrap/known-sites recalibration on inputBam and returns the output path.
			runBQSR := func(inputBam string) (string, error) {
				sampleKnownSites := make([]string, len(knownSites))
				copy(sampleKnownSites, knownSites)

				if len(sampleKnownSites) == 0 && bootstrap {
					color.Green("[%s] Creating known variants (bootstrap)...\n", sampleName)
					newKS, ksErr := CreateKnownVariants(resolvedFasta, inputBam, "", gatkLogLevel, verbose)
					if ksErr != nil {
						return "", fmt.Errorf("CreateKnownVariants: %v", ksErr)
					}
					sampleKnownSites = newKS
					color.Green("[%s] Known variants created: %v\n", sampleName, sampleKnownSites)
				}

				color.Green("[%s] Running Recalibrate...\n", sampleName)
				out, recErr := Recalibrate(resolvedFasta, inputBam, sampleKnownSites, "", verbose)
				if recErr != nil {
					return "", fmt.Errorf("Recalibrate: %v", recErr)
				}
				return out, nil
			}

			switch readPair.resume {

			// ------------------------------------------------------------------ //
			// resumeFromBamToCram: bqsr.bam exists → convert to bqsr.cram only   //
			// ------------------------------------------------------------------ //
			case resumeFromBamToCram:
				color.Green("[%s] Converting bqsr.bam → bqsr.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, bqsrCram); err != nil {
					addFailure(fmt.Sprintf("bqsr BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] Done. Final CRAM: %s\n", sampleName, bqsrCram)
				return

			// ------------------------------------------------------------------ //
			// resumeFromBQSR: rgmd.cram exists → run BQSR directly              //
			// ------------------------------------------------------------------ //
			case resumeFromBQSR:
				if !bqsr {
					color.Green("[%s] BQSR disabled — nothing to do (rgmd.cram already present).\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(readPair.intermediateBam)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromCramRgmdThenBQSR: rgmd.bam → cram it, then BQSR         //
			// ------------------------------------------------------------------ //
			case resumeFromCramRgmdThenBQSR:
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] rgmd.cram ready: %s\n", sampleName, rgmdCram)

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete.\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromMarkDupAndBQSR: sorted.bam → MarkDuplicates → BQSR → CRAM //
			// ------------------------------------------------------------------ //
			case resumeFromMarkDupAndBQSR:
				color.Green("[%s] Marking duplicates on sorted.bam...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, readPair.intermediateBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates complete.\n", sampleName)

				color.Green("[%s] Indexing rgmd.bam...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}

				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete. Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromAlign (default): full alignment from clean reads          //
			// ------------------------------------------------------------------ //
			default:
				// ----------------------------------------- BWA-MEM2 / BWA-MEM ----------------------------------------- //
				alignCmd := fmt.Sprintf(
					`bwa-mem2 mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
					threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				if aligner == "bwa-mem" {
					alignCmd = fmt.Sprintf(
						`bwa mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
						threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				}

				color.Green("[%s] Aligning with %s...\n", sampleName, aligner)
				fmt.Printf("[%s] %s\n", sampleName, alignCmd)
				var memErr error
				if verbose {
					memErr = utils.RunBashCmdVerbose(alignCmd)
				} else {
					memErr = utils.RunBashCmd(alignCmd)
				}
				if memErr != nil {
					addFailure(fmt.Sprintf("alignment: %v", memErr))
					return
				}
				color.Green("[%s] Alignment completed.\n", sampleName)

				// ----------------------------------------- Mark Duplicates -------------------------------------------- //
				color.Green("[%s] Marking duplicates...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates completed.\n", sampleName)

				// ----------------------------------------- BAM Index -------------------------------------------------- //
				color.Green("[%s] Indexing BAM...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				fmt.Printf("[%s] %s\n", sampleName, indexCmd)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}
				color.Green("[%s] BAM index completed.\n", sampleName)

				// ----------------------------------------- Convert to CRAM -------------------------------------------- //
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				// ============================================= BQSR ==================================================== //
				if !bqsr {
					color.Green("[%s] Pipeline complete (no BQSR). Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}

				color.Green("[%s] Starting BQSR...\n", sampleName)
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Final: %s\n", sampleName, bqsrOut)
			}

		}(rp)
	}

	wg.Wait()

	// ======================================= Final summary ======================================================= //
	fmt.Printf("\n==============================================================\n")
	color.Green("Already aligned (skipped): %d\n", len(validSamples))
	color.Green("Successfully processed:    %d\n", len(validCleanReads)-len(failedAlignments))
	if len(failedAlignments) > 0 {
		color.Red("Failed:                    %d\n", len(failedAlignments))
		for _, s := range failedAlignments {
			color.Red("  - %s\n", s)
		}
	}
	if len(missingCleanReads) > 0 {
		color.Yellow("Skipped (no valid reads):  %d\n", len(missingCleanReads))
		for _, s := range missingCleanReads {
			color.Yellow("  - %s\n", s)
		}
	}
	fmt.Printf("==============================================================\n\n")
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

	// ============================================= Scan alignment files =========================================== //

	bamsStatResults, _ := ScanAlignments(dataDir, species, refVer, genomesDir, resolvedFasta, verbose, quick)

	for _, s := range bamsStatResults {
		bqsrCramOK := isUsable(s.BqsrCram) && hasIndex(s.BqsrCram)
		rgmdCramOK := isUsable(s.RgmdCram) && hasIndex(s.RgmdCram)

		// --- If already perfect ---
		if bqsrCramOK && rgmdCramOK &&
			!s.SortedBam.Present &&
			!s.RgmdBam.Present &&
			!s.BqsrBam.Present &&
			len(s.OtherFiles) == 0 {

			color.Green("[%s] PASS ✅ - Skipping'\n", s.Sample)
		} else if !rgmdCramOK {
			if isUsable(s.RgmdBam) {
				color.Green("[%s] Converting rgmd.bam → rgmd.cram\n", s.Sample)
				bamToCramStr := fmt.Sprintf(`echo "samtools view -T %s -C -o %s %s"
    samtools view -T %s -C -o %s %s
    echo "samtools index %s"
    samtools index %s
`, resolvedFasta, strings.Replace(s.RgmdBam.Path, ".bam", ".cram", 1), s.RgmdBam.Path, resolvedFasta, strings.Replace(s.RgmdBam.Path, ".bam", ".cram", 1), s.RgmdBam.Path, strings.Replace(s.RgmdBam.Path, ".bam", ".cram", 1), strings.Replace(s.RgmdBam.Path, ".bam", ".cram", 1))
				if verbose {
					utils.RunBashCmdVerbose(bamToCramStr)
				} else {
					utils.RunBashCmd(bamToCramStr)
				}
			} else if isUsable(s.SortedBam) {
				steps = append(steps, "Run markdup: sorted.bam → rgmd.bam")
				steps = append(steps, "Convert rgmd.bam → rgmd.cram")
			} else {
				steps = append(steps, "Missing sorted.bam → must rerun alignment (bwa-mem)")
			}
		}

		// --- Ensure rgmd.cram index ---
		if isUsable(s.RgmdCram) && !hasIndex(s.RgmdCram) {
			steps = append(steps, "Index rgmd.cram")
		}

		// --- Ensure bqsr.cram ---
		if !bqsrCramOK {
			if isUsable(s.RgmdCram) || isUsable(s.RgmdBam) {
				steps = append(steps, "Run BQSR → bqsr.bam")
				steps = append(steps, "Convert bqsr.bam → bqsr.cram")
			} else {
				steps = append(steps, "Cannot run BQSR: missing rgmd input")
			}
		}

		// --- Ensure bqsr.cram index ---
		if isUsable(s.BqsrCram) && !hasIndex(s.BqsrCram) {
			steps = append(steps, "Index bqsr.cram")
		}

		// --- Cleanup phase ---
		if bqsrCramOK && rgmdCramOK {
			if s.SortedBam.Present {
				steps = append(steps, "Delete sorted.bam")
			}
			if s.RgmdBam.Present {
				steps = append(steps, "Delete rgmd.bam")
			}
			if s.BqsrBam.Present {
				steps = append(steps, "Delete bqsr.bam")
			}
			if len(s.OtherFiles) > 0 {
				steps = append(steps, "Remove unexpected extra files")
			}
		}

		return Action{
			Sample:  s.Sample,
			Steps:   steps,
			Perfect: false,
		}

	}

	color.Green("\nReady to align: %d\n", len(validCleanReads))
	color.Red("Missing/invalid reads: %d\n", len(missingCleanReads))
	color.Yellow("Skipped (multiple read files): %d\n", len(multipleCleanReads))
	fmt.Printf("\n==============================================================\n\n")

	if len(validCleanReads) == 0 {
		color.Yellow("No samples to align. Exiting.\n")
		return
	}

	// ======================================= Parallel alignment ================================================== //
	totalCores := runtime.NumCPU()
	maxParallelJobs := totalCores / threads
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
		threads = totalCores
	}
	color.Green("Aligning %d sample(s) — up to %d jobs in parallel, %d threads each.\n\n",
		len(validCleanReads), maxParallelJobs, threads)

	var (
		mu               sync.Mutex
		failedAlignments []string
		wg               sync.WaitGroup
		sem              = make(chan struct{}, maxParallelJobs)
	)

	for _, rp := range validCleanReads {
		wg.Add(1)
		sem <- struct{}{}

		go func(readPair samplePair) {
			defer wg.Done()
			defer func() { <-sem }()

			sampleName := readPair.sample

			addFailure := func(reason string) {
				color.Red("[%s] FAILED: %s\n", sampleName, reason)
				mu.Lock()
				failedAlignments = append(failedAlignments, sampleName)
				mu.Unlock()
			}

			// ---- Paths ---- //
			sortedBam := filepath.Join(readPair.bamDir, sampleName+".sorted.bam")
			rgmdBam := filepath.Join(readPair.bamDir, sampleName+".rgmd.bam")
			rgmdCram := filepath.Join(readPair.bamDir, sampleName+".rgmd.cram")
			rgmdMetrics := filepath.Join(readPair.bamDir, sampleName+".rgmd.metrics.txt")
			bqsrCram := filepath.Join(readPair.bamDir, sampleName+".bqsr.cram")

			// bamToCram converts src BAM/CRAM to dst CRAM using samtools.
			bamToCram := func(src, dst string) error {
				cmd := fmt.Sprintf(`samtools view -C -T %s -o %s %s && samtools index %s`, resolvedFasta, dst, src, dst)
				color.Green("[%s] Converting to CRAM: %s\n", sampleName, cmd)
				if verbose {
					return utils.RunBashCmdVerbose(cmd)
				}
				return utils.RunBashCmd(cmd)
			}

			// runBQSR runs bootstrap/known-sites recalibration on inputBam and returns the output path.
			runBQSR := func(inputBam string) (string, error) {
				sampleKnownSites := make([]string, len(knownSites))
				copy(sampleKnownSites, knownSites)

				if len(sampleKnownSites) == 0 && bootstrap {
					color.Green("[%s] Creating known variants (bootstrap)...\n", sampleName)
					newKS, ksErr := CreateKnownVariants(resolvedFasta, inputBam, "", gatkLogLevel, verbose)
					if ksErr != nil {
						return "", fmt.Errorf("CreateKnownVariants: %v", ksErr)
					}
					sampleKnownSites = newKS
					color.Green("[%s] Known variants created: %v\n", sampleName, sampleKnownSites)
				}

				color.Green("[%s] Running Recalibrate...\n", sampleName)
				out, recErr := Recalibrate(resolvedFasta, inputBam, sampleKnownSites, "", verbose)
				if recErr != nil {
					return "", fmt.Errorf("Recalibrate: %v", recErr)
				}
				return out, nil
			}

			switch readPair.resume {

			// ------------------------------------------------------------------ //
			// resumeFromBamToCram: bqsr.bam exists → convert to bqsr.cram only   //
			// ------------------------------------------------------------------ //
			case resumeFromBamToCram:
				color.Green("[%s] Converting bqsr.bam → bqsr.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, bqsrCram); err != nil {
					addFailure(fmt.Sprintf("bqsr BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] Done. Final CRAM: %s\n", sampleName, bqsrCram)
				return

			// ------------------------------------------------------------------ //
			// resumeFromBQSR: rgmd.cram exists → run BQSR directly              //
			// ------------------------------------------------------------------ //
			case resumeFromBQSR:
				if !bqsr {
					color.Green("[%s] BQSR disabled — nothing to do (rgmd.cram already present).\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(readPair.intermediateBam)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromCramRgmdThenBQSR: rgmd.bam → cram it, then BQSR         //
			// ------------------------------------------------------------------ //
			case resumeFromCramRgmdThenBQSR:
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(readPair.intermediateBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM conversion: %v", err))
					return
				}
				color.Green("[%s] rgmd.cram ready: %s\n", sampleName, rgmdCram)

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete.\n", sampleName)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromMarkDupAndBQSR: sorted.bam → MarkDuplicates → BQSR → CRAM //
			// ------------------------------------------------------------------ //
			case resumeFromMarkDupAndBQSR:
				color.Green("[%s] Marking duplicates on sorted.bam...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, readPair.intermediateBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates complete.\n", sampleName)

				color.Green("[%s] Indexing rgmd.bam...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}

				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				if !bqsr {
					color.Green("[%s] BQSR disabled — pipeline complete. Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Output: %s\n", sampleName, bqsrOut)
				return

			// ------------------------------------------------------------------ //
			// resumeFromAlign (default): full alignment from clean reads          //
			// ------------------------------------------------------------------ //
			default:
				// ----------------------------------------- BWA-MEM2 / BWA-MEM ----------------------------------------- //
				alignCmd := fmt.Sprintf(
					`bwa-mem2 mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
					threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				if aligner == "bwa-mem" {
					alignCmd = fmt.Sprintf(
						`bwa mem -t %v -M -Y -R '@RG\tID:%s.1\tSM:%s\tLB:%s_1\tPL:BGISEQ' %s %s %s | samtools sort -o %s`,
						threads, sampleName, sampleName, sampleName, resolvedFasta, readPair.fwd, readPair.rev, sortedBam)
				}

				color.Green("[%s] Aligning with %s...\n", sampleName, aligner)
				fmt.Printf("[%s] %s\n", sampleName, alignCmd)
				var memErr error
				if verbose {
					memErr = utils.RunBashCmdVerbose(alignCmd)
				} else {
					memErr = utils.RunBashCmd(alignCmd)
				}
				if memErr != nil {
					addFailure(fmt.Sprintf("alignment: %v", memErr))
					return
				}
				color.Green("[%s] Alignment completed.\n", sampleName)

				// ----------------------------------------- Mark Duplicates -------------------------------------------- //
				color.Green("[%s] Marking duplicates...\n", sampleName)
				mDupCmd := fmt.Sprintf(
					`gatk --java-options "-Xmx8G" MarkDuplicates -R %s -I %s -O %s -M %s --VERBOSITY %s`,
					resolvedFasta, sortedBam, rgmdBam, rgmdMetrics, gatkLogLevel)
				fmt.Printf("[%s] %s\n", sampleName, mDupCmd)
				var mdupErr error
				if verbose {
					mdupErr = utils.RunBashCmdVerbose(mDupCmd)
				} else {
					mdupErr = utils.RunBashCmd(mDupCmd)
				}
				if mdupErr != nil {
					addFailure(fmt.Sprintf("MarkDuplicates: %v", mdupErr))
					return
				}
				color.Green("[%s] MarkDuplicates completed.\n", sampleName)

				// ----------------------------------------- BAM Index -------------------------------------------------- //
				color.Green("[%s] Indexing BAM...\n", sampleName)
				indexCmd := fmt.Sprintf(`samtools index %s`, rgmdBam)
				fmt.Printf("[%s] %s\n", sampleName, indexCmd)
				var indErr error
				if verbose {
					indErr = utils.RunBashCmdVerbose(indexCmd)
				} else {
					indErr = utils.RunBashCmd(indexCmd)
				}
				if indErr != nil {
					addFailure(fmt.Sprintf("BAM index: %v", indErr))
					return
				}
				color.Green("[%s] BAM index completed.\n", sampleName)

				// ----------------------------------------- Convert to CRAM -------------------------------------------- //
				color.Green("[%s] Converting rgmd.bam → rgmd.cram...\n", sampleName)
				if err := bamToCram(rgmdBam, rgmdCram); err != nil {
					addFailure(fmt.Sprintf("rgmd BAM→CRAM: %v", err))
					return
				}

				// ============================================= BQSR ==================================================== //
				if !bqsr {
					color.Green("[%s] Pipeline complete (no BQSR). Final CRAM: %s\n", sampleName, rgmdCram)
					return
				}

				color.Green("[%s] Starting BQSR...\n", sampleName)
				bqsrOut, err := runBQSR(rgmdCram)
				if err != nil {
					addFailure(err.Error())
					return
				}
				color.Green("[%s] BQSR complete. Final: %s\n", sampleName, bqsrOut)
			}

		}(rp)
	}

	wg.Wait()

	// ======================================= Final summary ======================================================= //
	fmt.Printf("\n==============================================================\n")
	color.Green("Already aligned (skipped): %d\n", len(validSamples))
	color.Green("Successfully processed:    %d\n", len(validCleanReads)-len(failedAlignments))
	if len(failedAlignments) > 0 {
		color.Red("Failed:                    %d\n", len(failedAlignments))
		for _, s := range failedAlignments {
			color.Red("  - %s\n", s)
		}
	}
	if len(missingCleanReads) > 0 {
		color.Yellow("Skipped (no valid reads):  %d\n", len(missingCleanReads))
		for _, s := range missingCleanReads {
			color.Yellow("  - %s\n", s)
		}
	}
	fmt.Printf("==============================================================\n\n")
}
