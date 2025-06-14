package pangenome

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func LatestRef(outputDir, ref string) (int, string, error) {
	var pansList []string

	// Check if outputDir exists and is a directory
	info, err := os.Stat(outputDir)
	if err != nil || !info.IsDir() {
		return 0, ref, nil // If not a dir, return default
	}

	// Walk through outputDir
	entries, err := os.ReadDir(outputDir)
	if err != nil {
		return 0, ref, err
	}

	for _, entry := range entries {
		if entry.IsDir() {
			subDir := filepath.Join(outputDir, entry.Name())
			files, err := os.ReadDir(subDir)
			if err != nil {
				continue
			}
			for _, f := range files {
				if strings.HasSuffix(f.Name(), "_updated_ref.fa") {
					fullPath := filepath.Join(subDir, f.Name())
					pansList = append(pansList, fullPath)
				}
			}
		}
	}

	sort.Strings(pansList) // Ensure consistent ordering
	fmt.Println("PAN LIST:", pansList)

	if len(pansList) > 0 {
		latestRef := pansList[len(pansList)-1]
		latestBase := filepath.Base(latestRef)
		latestIndexStr := strings.SplitN(latestBase, "_", 2)[0]
		latestIndex, err := strconv.Atoi(latestIndexStr)
		if err != nil {
			return 0, ref, fmt.Errorf("invalid prefix in filename: %s", latestBase)
		}

		if latestIndex == 0 {
			fmt.Printf("Re-doing the first read pair %s — use initial reference!\n", latestRef)
			return 0, ref, nil
		} else {
			fmt.Printf("Sample %s is done! Starting with the next\n", latestRef)
			return latestIndex + 1, latestRef, nil
		}
	}

	return 0, ref, nil
}

func ExtractUnmappedReads(bamFile, outDir string) error {
	//samtools flagstat $bam_file > ${bam_file%.bam}_flagstats.txt

	cmdStr := fmt.Sprintf(`samtools view -u  -f 4 -F264 %s  > %s/tmps1.bam`, bamFile, outDir)
	cmdStr1 := fmt.Sprintf(`samtools view -u -f 8 -F 260 %s  > %s/tmps2.bam`, bamFile, outDir)
	cmdStr2 := fmt.Sprintf(`samtools view -u -f 12 -F 256 %s > %s/tmps3.bam`, bamFile, outDir)
	cmdStr3 := fmt.Sprintf(`samtools merge -u - %s/tmps[123].bam > %s/temps_unsorted.bam`, outDir, outDir)
	cmdStr4 := fmt.Sprintf(`samtools sort -n %s/temps_unsorted.bam -o %s`, outDir, strings.TrimSuffix(bamFile, ".bam")+"_unmapped.bam")

	fmt.Println(cmdStr)
	err := utils.RunBashCmdVerbose(cmdStr)
	if err != nil {
		fmt.Printf("Error running samtools view: %v\n", err)
		return err
	}
	fmt.Println(cmdStr1)
	err = utils.RunBashCmdVerbose(cmdStr1)
	if err != nil {
		fmt.Printf("Error running samtools view: %v\n", err)
		return err
	}
	fmt.Println(cmdStr2)
	err = utils.RunBashCmdVerbose(cmdStr2)
	if err != nil {
		fmt.Printf("Error running samtools view: %v\n", err)
		return err
	}
	fmt.Println(cmdStr3)
	err = utils.RunBashCmdVerbose(cmdStr3)
	if err != nil {
		fmt.Printf("Error running samtools merge: %v\n", err)
	}

	fmt.Println(cmdStr4)
	err = utils.RunBashCmdVerbose(cmdStr4)
	if err != nil {
		fmt.Printf("Error running samtools sort: %v\n", err)
		return err
	}

	fmt.Println("removing temporary files")
	err = os.RemoveAll(outDir + "/tmps1.bam")
	if err != nil {
		return err
	}
	err = os.RemoveAll(outDir + "/tmps2.bam")
	if err != nil {
		return err
	}
	err = os.RemoveAll(outDir + "/tmps3.bam")
	if err != nil {
		return err
	}
	err = os.RemoveAll(outDir + "/temps_unsorted.bam")
	if err != nil {
		return err
	}
	return nil

}

func GoPan(config string, assembler string, threads int) {

	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(config)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Species:", cfg.Species)
	fmt.Println("Read Pairs:", cfg.ReadPairs)
	readPairs := cfg.ReadPairs
	ref := cfg.Reference
	out := cfg.OutputDir

	//--------------------------------------------- Check Paths ----------------------------------------------------- //
	_, refErr := os.Stat(ref)
	if refErr != nil {
		fmt.Printf("Reference genome path: %s, is not valid\n", ref)
		return
	}
	info, err := os.Stat(out)
	if err != nil || !info.IsDir() {
		fmt.Printf("Output directory: %s is not a valid path\n", out)
		return
	}

	for _, pair := range readPairs {
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
			fmt.Println("Please provide sample name")
			return
		}
		if lb == "" {
			fmt.Println("Please provide library name")
			return
		}
	}

	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(out, "pangenome.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	mw := io.MultiWriter(logFile, os.Stdout)
	log.SetOutput(mw)
	fmt.Println("Log file created.")
	//--------------------------------------- Start goPan ----------------------------------------------------------- //
	i, latestRef, err := LatestRef(out, ref)
	fmt.Println("Check log for latest reference genome ...")
	for i < len(readPairs) {
		fmt.Println(i, readPairs[i])
		// --------------------------- Align reads to latest reference ---------------------------------------------- //
		fmt.Println("Aligning reads to latest reference")
		fwd, rev, sn, lb := readPairs[i][0], readPairs[i][1], readPairs[i][2], readPairs[i][3]

		sampleDir, _ := filepath.Abs(filepath.Join(out, sn))

		bamFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted.bam"))
		baiFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted.bai"))

		_, bamErr := os.Stat(bamFile)
		_, baiErr := os.Stat(baiFile)

		fmt.Printf("Check log to see if alignment for %s is finished. If not, run alignment again. Otherwise skip\n", sn)

		if bamErr == nil && baiErr == nil {
			fmt.Println("Bam file and bai file exist for ", sn, "skip alignment")
			continue
		} else {
			fmt.Printf("Bam file and bai file do not exist for %s \n", sn)
			fmt.Printf("Align %s and %s to  %s using bowtie2 \n", fwd, rev, sn)
			log.Printf("%s\tBOWTIE2\t%v\t%sSTARTED", sn, i, latestRef)

			aErr := alignment.AlignShortReadsBt(latestRef, fwd, rev, sn, lb, sampleDir, 8)
			if aErr != nil {
				fmt.Printf("Error aligning reads: %v\n", aErr)
				log.Printf("%s\tBOWTIE2\t%v\t%sFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tBOWTIE2\t%v\t%sFINISHED", sn, i, latestRef)
		}

		// ---------------------------------- Alignment Statistics -------------------------------------------------- //

		pdfFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_insert_size_histogram.pdf"))
		alignMetrics, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_align_metrics.txt"))
		insertMetrics, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_insert_metrics.txt"))
		flagStats, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_flagstats.txt"))

		_, pdfErr := os.Stat(pdfFile)
		_, alignErr := os.Stat(alignMetrics)
		_, insertErr := os.Stat(insertMetrics)
		_, flagErr := os.Stat(flagStats)
		if pdfErr == nil && alignErr == nil && insertErr == nil && flagErr == nil {
			fmt.Println("Alignment stats for ", sn, "exist. skip")
			continue
		} else {
			fmt.Printf("Alignment stats files do not exist for %s \n", sn)
			fmt.Printf("Getting alignment stats for %s ... \n", sn)
			log.Printf("%s\tALIGNMENT_STATS\t%v\t%sSTARTED", sn, i, latestRef)
			sErr := alignment.AlignmentStats(ref, bamFile)
			if sErr != nil {
				fmt.Printf("Error getting alignment stats: %v\n", sErr)
				log.Printf("%s\tALIGNMENT_STATS\t%v\t%sFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tALIGNMENT_STATS\t%v\t%sFINISHED", sn, i, latestRef)

		}

		// -------------------------------------- Extract reads ----------------------------------------------------- //
		unmappedBam := strings.TrimSuffix(bamFile, ".bam") + "_unmapped.bam"
		_, unmappedBamErr := os.Stat(unmappedBam)

		if unmappedBamErr == nil {
			fmt.Println("Unmapped bam and bai files exist for ", sn, "skip")
			continue
		} else {
			fmt.Printf("Unmapped bam and bai files do not exist for %s \n", sn)
			fmt.Printf("Extract unmapped reads for %s ... \n", sn)
			log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%sSTARTED", sn, i, latestRef)
			eErr := ExtractUnmappedReads(unmappedBam, sampleDir)
			if eErr != nil {
				fmt.Printf("Error extracting unmapped reads: %v\n", eErr)
				log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%sFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%sFINISHED", sn, i, latestRef)
		}

		// ---------------------------------- Convert unmapped bams to fasta ---------------------------------------- //
		unmappedFwdReads, _ := filepath.Abs(strings.Replace(unmappedBam, ".bam", "_1.fastq", 1))
		unmappedRevReads, _ := filepath.Abs(strings.Replace(unmappedBam, ".bam", "_2.fastq", 1))

		_, unmappedFwdReadsErr := os.Stat(unmappedFwdReads)
		_, unmappedRevReadsErr := os.Stat(unmappedRevReads)

		if unmappedFwdReadsErr == nil && unmappedRevReadsErr == nil {
			fmt.Println("Unmapped fastq files exist for ", sn, "skip")
			continue
		} else {
			fmt.Printf("Unmapped bam and bai files do not exist for %s \n", sn)
			fmt.Printf("Extract unmapped reads for %s ... \n", sn)
			log.Printf("%s\tBAM_TO_FASTQ\t%v\t%sSTARTED", sn, i, latestRef)
			cmdStr := fmt.Sprintf(`bedtools bamtofastq -i %s -fq %s -fq2 %s`, unmappedBam, unmappedFwdReads, unmappedRevReads)
			fmt.Println(cmdStr)
			bErr := utils.RunBashCmdVerbose(cmdStr)
			if bErr != nil {
				fmt.Printf("Error running bedtools bamtofastq: %v\n", bErr)
				log.Printf("%s\tBAM_TO_FASTQ\t%v\t%sFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tBAM_TO_FASTQ\t%v\t%sFINISHED", sn, i, latestRef)
		}

		//---------------------------------- Assemble unmapped reads ------------------------------------------------ //
		if assembler == "masurca" {

		} else if assembler == "megahit" {
			megahitDir := filepath.Join(sampleDir, "MegaHit")
			megahitAssembly := filepath.Join(megahitDir, "final.contigs.fa")

			_, megahitErr := os.Stat(megahitAssembly)
			if megahitErr == nil {
				fmt.Println("Megahit assembly exists for ", sn, "skip")
				continue
			} else {
				err := os.RemoveAll(megahitDir)
				if err != nil {
					fmt.Printf("Error removing megahit directory: %v\n", err)
					return
				}
				fmt.Printf("Megahit assembly does not exist for %s \n", sn)
				fmt.Printf("Assemble unmapped reads for %s ... \n", sn)
				log.Printf("%s\tMEGAHIT\t%v\t%sSTARTED", sn, i, latestRef)
				cmdStr := fmt.Sprintf(`megahit  -1 %s -2 %s -o %s`, unmappedFwdReads, unmappedRevReads, megahitDir)
				fmt.Println(cmdStr)
				mErr := utils.RunBashCmdVerbose(cmdStr)
				if mErr != nil {
					fmt.Printf("Error running megahit: %v\n", mErr)
					log.Printf("%s\tMEGAHIT\t%v\t%sFAILED", sn, i, latestRef)
					return
				}
				log.Printf("%s\tMEGAHIT\t%v\t%sFINISHED", sn, i, latestRef)
			}

		} else if assembler == "both" {
			fmt.Println("Invalid assembler")

		} else {
			fmt.Println("Invalid assembler")
		}

	}

}
