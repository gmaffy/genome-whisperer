package alignment

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sync"

	"github.com/gmaffy/genome-whisperer/utils"
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
				msg := fmt.Sprintf("BOOTSTRAP already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
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
				fmt.Printf("Known variants created at %s\n...................................................\n", knownSites)

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
