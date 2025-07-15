package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sync"
)

func RunAlignReads(referencePath string, forwardPath string, reversePath string, sePath string, sampleName string, libName string, outputDir string, threads int, aligner string, knownSites []string, bqsr bool, bootstrap bool, logFilePath string, preset string) error {

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
	rgBam := fmt.Sprintf("%s/%s.RG.bam", lineDir, sampleName)
	rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sampleName)
	rgmdMetrics := fmt.Sprintf("%s/%s.RGMD.metrics.txt", lineDir, sampleName)
	rgmdIndex := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sampleName)
	sortedBam := fmt.Sprintf("%s/%s.sorted.bam", lineDir, sampleName)
	sortedBai := fmt.Sprintf("%s/%s.sorted.bai", lineDir, sampleName)

	bErr := os.MkdirAll(lineDir, 0755)
	if bErr != nil {
		//jlog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bErr))    //, "CMD", "ALL")
		//slog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bErr))
		log.Fatalf("Error creating results directory: %s\n", bErr)
		return bErr
	}

	if aligner == "bwa-mem" {
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

			memErr := utils.RunBashCmdVerbose(cmdStr)
			if memErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", memErr))
				return memErr
			}

			jlog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BWA_MEM", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- Mark Duplicates -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "MARK_DUPLICATES", sampleName, "ALL") {
			msg := fmt.Sprintf("MARK_DUPLICATES already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {

			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" MarkDuplicates -I %s -O %s -M %s`, sortedBam, rgmdBam, rgmdMetrics)
			fmt.Printf("%s\n-----------------------------------------------\n\n", mDupCmdStr)

			mdupErr := utils.RunBashCmdVerbose(mDupCmdStr)
			if mdupErr != nil {
				jlog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mdupErr))
				return mdupErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "MARK_DUPLICATES", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

		fmt.Printf("\n\n------------------------------------------- BAM INDEXING -------------------------------------------------- \n\n")

		if utils.StageHasCompleted(logged, "BAM_INDEX", sampleName, "ALL") {
			msg := fmt.Sprintf("BAM INDEXING already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
			slog.Info(msg)
		} else {
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")

			indexCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, rgmdBam, rgmdIndex)
			fmt.Printf("%s\n-----------------------------------------------\n\n", indexCmdStr)

			indErr := utils.RunBashCmdVerbose(indexCmdStr)
			if indErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr)) //, "CMD", "ALL")
				slog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				return indErr
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
			bowErr := utils.RunBashCmdVerbose(cmdStr)
			if bowErr != nil {
				jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bowErr))
				slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", bowErr))
				return bowErr
			}

			fmt.Printf("Indexing Bam ....")
			mDupCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8G" BuildBamIndex -I %s -O %s`, sortedBam, sortedBai)
			mErr := utils.RunBashCmdVerbose(mDupCmdStr)
			if mErr != nil {
				jlog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mErr))
				slog.Info("ALIGNMENT", "PROGRAM", "BOWTIE2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mErr))
				return mErr
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
				indexErr := utils.RunBashCmdVerbose(indexCmdStr)
				if indexErr != nil {
					fmt.Println("Indexing reference failed")
					return indexErr
				}
			}
			pbmm2CmdStr := fmt.Sprintf(`pbmm2 align --sort -j %v --preset %s %s %s %s`, threads, preset, referencePath+".mmi", sePath, sortedBam)
			fmt.Println(pbmm2CmdStr)
			pbmm2Err := utils.RunBashCmdVerbose(pbmm2CmdStr)
			if pbmm2Err != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", pbmm2Err))
				slog.Error("ALIGNMENT", "PROGRAM", "PBMM2", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", pbmm2Err))
				//fmt.Println("pbmm2 alignment failed")
				return pbmm2Err
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
			rgErr := utils.RunBashCmdVerbose(rgCmdStr)
			if rgErr != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", rgErr))
				slog.Error("ALIGNMENT", "PROGRAM", "RG", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", rgErr))
				return rgErr
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
			mkdpErr := utils.RunBashCmdVerbose(mkdpCmdStr)
			if mkdpErr != nil {
				jlog.Error("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mkdpErr))
				slog.Error("ALIGNMENT", "PROGRAM", "PBMARKDUP", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", mkdpErr))
				return mkdpErr
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

			indErr := utils.RunBashCmdVerbose(indexCmdStr)
			if indErr != nil {
				jlog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				slog.Error("BQSR", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", indErr))
				return indErr
			}
			jlog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("ALIGNMENT", "PROGRAM", "BAM_INDEX", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

	} else {
		fmt.Println("Aligner not supported")
		return fmt.Errorf("aligner not supported\n\n Only bwa-mem, bowtie2 and pbmm2 are supported\n\n")
	}

	if bqsr {

		fmt.Println("Running BQSR")
		fmt.Println("Check if mark duplicates was successful", knownSites)
		_, mdErr := os.Stat(rgmdBam)
		_, mdiErr := os.Stat(rgmdIndex)

		if mdErr != nil || mdiErr != nil {
			fmt.Println("Mark duplicates failed")
			return fmt.Errorf("mark duplicates failed\n\n")
		}

		if len(knownSites) == 0 && bootstrap == true {
			fmt.Println("Running with bootstrap method")
			fmt.Println("Create Known variants", knownSites)
			if utils.StageHasCompleted(logged, "BOOTSTRAP_CKV", sampleName, "ALL") {
				msg := fmt.Sprintf("BOOTSTRAP already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
			} else {
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				newKnownSites, ksErr := CreateKnownVariants(referencePath, rgmdBam, logFilePath)
				if ksErr != nil {
					jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
					slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", ksErr))
					return fmt.Errorf("create known variants failed\n\n")
				}
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_CSV", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				knownSites = newKnownSites
				fmt.Printf("Known variants created at %s\n...................................................\n", knownSites)

			}

			if utils.StageHasCompleted(logged, "BOOTSTRAP_BQSR", sampleName, "ALL") {
				msg := fmt.Sprintf("BOOTSTRAP_BQSR already completed for bam file: %s. Skipping\n\n------------------------------------------------------------------------\n\n", sampleName)
				slog.Info(msg)
			} else {
				jlog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath)
				if recErr != nil {
					jlog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					slog.Error("BQSR", "PROGRAM", "BOOTSTRAP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					return fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
				}
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
				recErr := Recalibrate(referencePath, rgmdBam, knownSites, logFilePath)
				if recErr != nil {
					jlog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					slog.Error("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %s", recErr))
					return fmt.Errorf("BQSR failed - %s ................................\n\n", recErr)
				}
				jlog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
				slog.Info("BQSR", "PROGRAM", "DBSNP_BQSR", "SAMPLE", sampleName, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

			}

		}
	}
	return nil
}

func RunAlignReadsConfig(configPath string, threadsPerSample int, bqsr bool, bootstrap bool, aligner string, preset string) {

	fmt.Printf(" ---------------------------------------- Checking File Paths --------------------------------------------------------- \n\n")
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}

	ref := cfg.Reference
	_, refErr := os.Stat(ref)
	if refErr != nil {
		fmt.Printf("Reference genome path: %s, is not valid\n", ref)
		return
	}

	out := cfg.OutputDir
	outInfo, outErr := os.Stat(out)

	if outErr != nil {

		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", out)
			if createErr := os.MkdirAll(out, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", out, createErr)
				return
			}
			fmt.Printf("Output directory %s created successfully.\n", out)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", out, outErr)
			return
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output Directory %s file path is not a directory\n", out)
		return
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
			return
		}
		if len(knownSites) == 0 && bootstrap == false {
			fmt.Println("Either pass a known-sites file or enable bootstrap method")
			return
		} else if len(knownSites) > 0 {
			fmt.Println("Running with known-sites flag")
			// ---------------------------- Checking Known sites file paths ----------------------------------------- //
			for j, _ := range knownSites {
				_, err := os.Stat(knownSites[j])
				if err != nil {
					fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
					log.Fatal(err)
				}
			}
			if bootstrap == true {
				fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
				return
			}
		}
	}

	i := 0
	if aligner == "bwa-mem" || aligner == "bowtie2" {
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
				return
			}

			if revErr != nil {
				fmt.Printf("Reverse reads path %s, is not valid\n", rev)
				return
			}

			if sn == "" {
				fmt.Println("Please provide sample name ")
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				return
			}
			if lb == "" {
				fmt.Println("Please provide library name  ")
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				return
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
				return
			}
			if sn == "" {
				fmt.Println("Please provide sample name ")
			}
			if lb == "" {
				fmt.Println("Please provide library name  ")
			}

		}

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
	var bams []string
	if aligner == "bwa-mem" || aligner == "bowtie2" {
		if aligner == "bwa-mem" {
			depErr := utils.CheckDeps([]string{"bwa", "gatk", "samtools"})
			if depErr != nil {
				fmt.Printf("Dependency check failed!\n")
				return
			}
		} else {
			depErr := utils.CheckDeps([]string{"bowtie2", "samtools", "gatk"})
			if depErr != nil {
				fmt.Printf("Dependency check failed!\n")
			}
			if depErr != nil {
				fmt.Printf("Dependency check failed!\n")
			}
		}

		for _, pair := range cfg.ReadPairs {
			if len(pair) < 4 {
				fmt.Printf("This read pair is wrongly formated %s\n", pair)
				fmt.Println("Supply reads in this format: ReadPair: <fwd reads> <rev reads> <sample name> <library name> ")
				continue
			}

			wg.Add(1)
			sem <- struct{}{}
			go func(pair []string) {
				defer wg.Done()
				defer func() { <-sem }()

				fwd, rev, sn, lb := pair[0], pair[1], pair[2], pair[3]
				lineDir := fmt.Sprintf("%s/%s", out, sn)
				rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sn)
				rgmdBai := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sn)

				jlog.Info("ALIGNMENT", "PROGRAM", "SR_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("ALIGNMENT", "PROGRAM", "SR_ALIGNMENT", "SAMPLE", sn, "STATUS", "STARTED")

				_, rgmdBamErr := os.Stat(rgmdBam)
				_, rgmdBaiErr := os.Stat(rgmdBai)

				isDone := utils.StageHasCompleted(logged, "PE_ALIGNMENT", sn, "ALL")
				if isDone && rgmdBamErr == nil && rgmdBaiErr == nil {
					msg := fmt.Sprintf("%s and MarkDuplicates already completed for %s. Skipping.\n\n-------------------------------------------------------\n\n", aligner, sn)
					slog.Info(msg)

				} else {
					alErr := RunAlignReads(cfg.Reference, fwd, rev, "", sn, lb, cfg.OutputDir, threadsPerSample, aligner, knownSites, bqsr, bootstrap, logFilePath, preset)
					if alErr != nil {
						jlog.Error("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", alErr))
						slog.Error("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED - %v", alErr))

						return
					}
					jlog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
					slog.Info("ALIGNMENT", "PROGRAM", "PE_ALIGNMENT", "SAMPLE", sn, "STATUS", "COMPLETED")
				}
				bams = append(bams, rgmdBam)

			}(pair)

		}
		wg.Wait()

	} else if aligner == "pbmm2" {
		depErr := utils.CheckDeps([]string{"pbmm2", "pbmarkdup", "samtools", "gatk"})
		if depErr != nil {
			fmt.Printf("Dependency check failed!\n")
			return
		}
		for _, se := range cfg.SeReads {
			if len(se) < 3 {
				fmt.Printf("This read pair is wrongly formated %s\n", se)
				fmt.Println("Supply reads in this format: SingleRead: <reads> <sample name> <library name> ")
				continue
			}
			wg.Add(1)
			sem <- struct{}{}
			go func(pair []string) {
				defer wg.Done()
				defer func() { <-sem }()
				seRead, sn, lb := se[0], se[1], se[2]
				lineDir := fmt.Sprintf("%s/%s", out, sn)
				rgmdBam := fmt.Sprintf("%s/%s.RGMD.bam", lineDir, sn)
				rgmdBai := fmt.Sprintf("%s/%s.RGMD.bai", lineDir, sn)

				jlog.Info("ALIGNMENT", "PROGRAM", "SR_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "STARTED")
				slog.Info("ALIGNMENT", "PROGRAM", "SR_ALIGNMENT", "SAMPLE", sn, "STATUS", "STARTED")

				_, rgmdBamErr := os.Stat(rgmdBam)
				_, rgmdBaiErr := os.Stat(rgmdBai)

				isDone := utils.StageHasCompleted(logged, "PB_ALIGNMENT", sn, "ALL")
				if isDone && rgmdBamErr == nil && rgmdBaiErr == nil {
					msg := fmt.Sprintf("%s and MarkDuplicates already completed for %s. Skipping.\n\n-------------------------------------------------------\n\n", aligner, sn)
					slog.Info(msg)

				} else {
					alErr := RunAlignReads(cfg.Reference, "", "", seRead, sn, lb, cfg.OutputDir, threadsPerSample, aligner, knownSites, bqsr, bootstrap, logFilePath, preset)
					if alErr != nil {
						jlog.Error("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", alErr))
						slog.Error("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "STATUS", fmt.Sprintf("FAILED - %v", alErr))

						return
					}
					jlog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
					slog.Info("ALIGNMENT", "PROGRAM", "PB_ALIGNMENT", "SAMPLE", sn, "STATUS", "COMPLETED")
				}
				bams = append(bams, rgmdBam)
			}(se)
		}
		wg.Wait()

	} else {
		Msg := fmt.Sprintf("Aligner %s is not supported", aligner)
		slog.Error(Msg)
		return
	}

}
