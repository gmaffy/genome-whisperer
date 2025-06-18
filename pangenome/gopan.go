package pangenome

import (
	"bufio"
	"fmt"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func LatestRef(outputDir, ref string) (int, string, error) {
	var pansList []string

	info, err := os.Stat(outputDir)
	if err != nil || !info.IsDir() {
		return 0, ref, nil
	}

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

	sort.Strings(pansList)
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

func RenameScaffs(trimedFasta string, renamedFasta string, sn string) error {
	fmt.Println("Renaming scaffolds in fasta file")
	f, err := os.Open(trimedFasta)
	if err != nil {
		return err
	}
	defer f.Close()

	r, rErr := os.Create(renamedFasta)
	if rErr != nil {
		return rErr
	}
	defer r.Close()

	lineScanner := bufio.NewScanner(f)
	i := 0
	totalBytes := 0
	for lineScanner.Scan() {
		line := lineScanner.Text()
		if strings.HasPrefix(line, ">") {
			line = fmt.Sprintf(">%s_scaff_%d", sn, i+1)
			bytes, _ := fmt.Fprintln(r, line)
			totalBytes += bytes
		} else {
			bytes, _ := fmt.Fprintln(r, line)
			totalBytes += bytes
		}
		i++

	}
	if err := lineScanner.Err(); err != nil {
		return err
	}
	fmt.Printf("Total bytes written: %d\n", totalBytes)
	return nil
}

func TrimFasta(fasta string, trimmedFasta string) error {
	fmt.Println("Trimming fasta file")
	f, err := os.Open(fasta)
	if err != nil {
		return err
	}
	defer f.Close()

	return nil
}

func ConcatFasta(fastas []string, outFasta string) error {

	outFile, err := os.Create(outFasta)
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	defer outFile.Close()

	fmt.Println("Concatenating fasta files")
	for _, f := range fastas {
		inFile, err := os.Open(f)
		if err != nil {
			log.Fatalf("Failed to open input file %s: %v", f, err)
		}

		_, err = io.Copy(outFile, inFile)
		if err != nil {
			log.Fatalf("Failed to copy contents from %s: %v", f, err)
		}
		inFile.Close()
	}
	return nil
}

func GoPan(config string, assembler string) {

	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(config)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}

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
	fmt.Println("Starting goPan")
	i, latestRef, err := LatestRef(out, ref)
	fmt.Printf("Latest reference is: %s\n\n------------------------------------------\n\n", latestRef)
	for i < len(readPairs) {

		// --------------------------- Align reads to latest reference ---------------------------------------------- //
		fmt.Println("Aligning reads to latest reference")
		fwd, rev, sn, lb := readPairs[i][0], readPairs[i][1], readPairs[i][2], readPairs[i][3]

		sampleDir, _ := filepath.Abs(filepath.Join(out, sn))

		bamFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted.bam"))
		baiFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted.bai"))

		_, bamErr := os.Stat(bamFile)
		_, baiErr := os.Stat(baiFile)

		fmt.Printf("Checking log to see if alignment for %s is finished... \n\n", sn)

		if bamErr == nil && baiErr == nil {
			fmt.Printf("Bam file and bai file exist for: %s. Skip alignment....\n\n-------------------------------------------\n\n", sn)
			//continue
		} else {
			fmt.Printf("Bam file and bai file do not exist for %s \n", sn)
			fmt.Printf("Align %s and %s to  %s using bowtie2 \n", fwd, rev, sn)
			log.Printf("%s\tBOWTIE2\t%v\t%s\tSTARTED", sn, i, latestRef)

			aErr := alignment.AlignShortReadsBt(latestRef, fwd, rev, sn, lb, sampleDir, 8)
			if aErr != nil {
				fmt.Printf("Error aligning reads: %v\n", aErr)
				log.Printf("%s\tBOWTIE2\t%v\t%s\tFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tBOWTIE2\t%v\t%s\tFINISHED", sn, i, latestRef)
		}

		// ---------------------------------- Alignment Statistics -------------------------------------------------- //

		pdfFile, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_insert_size_histogram.pdf"))
		alignMetrics, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_alignment_metrics.txt"))
		insertMetrics, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_insert_metrics.txt"))
		flagStats, _ := filepath.Abs(filepath.Join(sampleDir, sn+".sorted_flagstats.txt"))
		//fmt.Printf("pdfFile: %s\n, alignMetrics:%s\n, insertMetrics:%s\n, flagStats:%s\n ", pdfFile, alignMetrics, insertMetrics, flagStats)

		_, pdfErr := os.Stat(pdfFile)
		_, alignErr := os.Stat(alignMetrics)
		_, insertErr := os.Stat(insertMetrics)
		_, flagErr := os.Stat(flagStats)
		//fmt.Printf("pdfErr: %s\n, alignErr:%s\n, insertErr:%s\n, flagErr:%s\n\n\n----------------\n\n ", pdfErr, alignErr, insertErr, flagErr)
		if pdfErr == nil && alignErr == nil && insertErr == nil && flagErr == nil {
			fmt.Println("Alignment stats for ", sn, "exist. skip")
			//continue
		} else {
			fmt.Printf("Alignment stats files do not exist for %s \n\n", sn)
			fmt.Printf("Getting alignment stats for %s ... \n", sn)
			log.Printf("%s\tALIGNMENT_STATS\t%v\t%s\tSTARTED", sn, i, latestRef)
			sErr := alignment.AlignmentStats(ref, bamFile)
			if sErr != nil {
				fmt.Printf("Error getting alignment stats: %v\n", sErr)
				log.Printf("%s\tALIGNMENT_STATS\t%v\t%s\tFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tALIGNMENT_STATS\t%v\t%s\tFINISHED", sn, i, latestRef)

		}

		// -------------------------------------- Extract reads ----------------------------------------------------- //
		unmappedBam := strings.TrimSuffix(bamFile, ".bam") + "_unmapped.bam"
		_, unmappedBamErr := os.Stat(unmappedBam)

		if unmappedBamErr == nil {
			fmt.Println("Unmapped bam and bai files exist for ", sn, "skip")
			//continue
		} else {
			fmt.Printf("Unmapped bam and bai files do not exist for %s \n", sn)
			fmt.Printf("Extract unmapped reads for %s ... \n", sn)
			log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%s\tSTARTED", sn, i, latestRef)
			eErr := ExtractUnmappedReads(bamFile, sampleDir)
			if eErr != nil {
				fmt.Printf("Error extracting unmapped reads: %v\n", eErr)
				log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%s\tFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tEXTRACT_UNMAPPED_READS\t%v\t%s\tFINISHED", sn, i, latestRef)
		}

		// ---------------------------------- Convert unmapped bams to fasta ---------------------------------------- //
		unmappedFwdReads, _ := filepath.Abs(strings.Replace(unmappedBam, ".bam", "_1.fastq", 1))
		unmappedRevReads, _ := filepath.Abs(strings.Replace(unmappedBam, ".bam", "_2.fastq", 1))

		_, unmappedFwdReadsErr := os.Stat(unmappedFwdReads)
		_, unmappedRevReadsErr := os.Stat(unmappedRevReads)

		if unmappedFwdReadsErr == nil && unmappedRevReadsErr == nil {
			fmt.Println("Unmapped fastq files exist for ", sn, "skip")
			//continue
		} else {
			fmt.Printf("Unmapped fastq files do not exist for %s \n", sn)
			fmt.Printf("Extract unmapped reads for %s ... \n", sn)
			log.Printf("%s\tBAM_TO_FASTQ\t%v\t%s\tSTARTED", sn, i, latestRef)
			cmdStr := fmt.Sprintf(`bedtools bamtofastq -i %s -fq %s -fq2 %s`, unmappedBam, unmappedFwdReads, unmappedRevReads)
			fmt.Println(cmdStr)
			bErr := utils.RunBashCmdVerbose(cmdStr)
			if bErr != nil {
				fmt.Printf("Error running bedtools bamtofastq: %v\n", bErr)
				log.Printf("%s\tBAM_TO_FASTQ\t%v\t%s\tFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tBAM_TO_FASTQ\t%v\t%s\tFINISHED", sn, i, latestRef)
		}

		//---------------------------------- Assemble unmapped reads ------------------------------------------------ //
		var assembledContigs string
		if assembler == "masurca" {
			fmt.Printf("Assemble unmapped reads for %s with MASURCA only... \n", sn)
			masurcaDir := filepath.Join(sampleDir, "MASURCA")
			masurcaAssembly := filepath.Join(masurcaDir, "CA", "primary.genome.scf.fasta")
			_, masurcaErr := os.Stat(masurcaAssembly)
			if masurcaErr == nil {
				fmt.Println("MASURCA assembly exists for ", sn, "skip")
				//continue
			} else {
				log.Printf("%s\tMASURCA\t%v\t%s\tSTARTED", sn, i, latestRef)
				err := os.RemoveAll(masurcaDir)
				if err != nil {
					fmt.Printf("Error removing masurca directory: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				err = os.MkdirAll(masurcaDir, os.ModePerm)
				if err != nil {
					fmt.Printf("Error creating masurca directory: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				cmd := exec.Command("masurca", "-t", "32", "-i", fmt.Sprintf("%s,%s", unmappedFwdReads, unmappedRevReads))
				cmd.Dir = masurcaDir
				fmt.Println(cmd.String())
				err = cmd.Run()
				if err != nil {
					fmt.Printf("Error running masurca: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
				}
				log.Printf("%s\tMASURCA\t%v\t%s\tFINISHED", sn, i, latestRef)
				assembledContigs = masurcaAssembly
			}

		} else if assembler == "megahit" {

			fmt.Printf("Assemble unmapped reads for %s with MEGAHIT only... \n", sn)
			megahitDir := filepath.Join(sampleDir, "MegaHit")
			megahitAssembly := filepath.Join(megahitDir, "final.contigs.fa")

			_, megahitErr := os.Stat(megahitAssembly)
			if megahitErr == nil {
				fmt.Println("Megahit assembly exists for ", sn, "skip")
				//continue
			} else {
				log.Printf("%s\tMEGAHIT\t%v\t%s\tSTARTED", sn, i, latestRef)
				err := os.RemoveAll(megahitDir)
				if err != nil {
					fmt.Printf("Error removing megahit directory: %v\n", err)
					log.Fatalf("%s\tMEGAHIT\t%v\t%s\tFAILED", sn, i, latestRef)
					return
				}
				fmt.Printf("Megahit assembly does not exist for %s \n", sn)
				fmt.Printf("Assemble unmapped reads for %s ... \n", sn)
				log.Printf("%s\tMEGAHIT\t%v\t%s\tSTARTED", sn, i, latestRef)
				cmdStr := fmt.Sprintf(`megahit  -1 %s -2 %s -o %s`, unmappedFwdReads, unmappedRevReads, megahitDir)
				fmt.Println(cmdStr)
				mErr := utils.RunBashCmdVerbose(cmdStr)
				if mErr != nil {
					fmt.Printf("Error running megahit: %v\n", mErr)
					log.Fatalf("%s\tMEGAHIT\t%v\t%s\tFAILED", sn, i, latestRef)
					return
				}
				log.Printf("%s\tMEGAHIT\t%v\t%s\tFINISHED", sn, i, latestRef)
				assembledContigs = megahitAssembly
			}

		} else {
			// --------------------------------------- BOTH ASSEMBLERS ---------------------------------------------- //

			// --------------------------------------- Run Masurca first -------------------------------------------- //
			fmt.Printf("Assemble unmapped reads for %s with MASURCA and MEGAHIT ... \n", sn)
			masurcaDir := filepath.Join(sampleDir, "MASURCA")
			masurcaAssembly := filepath.Join(masurcaDir, "CA", "primary.genome.scf.fasta")

			_, masurcaErr := os.Stat(masurcaAssembly)
			if masurcaErr == nil {
				fmt.Println("MASURCA assembly exists for ", sn, "skip")
				//continue
			} else {
				log.Printf("%s\tMASURCA\t%v\t%s\tSTARTED", sn, i, latestRef)
				err := os.RemoveAll(masurcaDir)
				if err != nil {
					fmt.Printf("Error removing masurca directory: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				err = os.MkdirAll(masurcaDir, os.ModePerm)
				if err != nil {
					fmt.Printf("Error creating masurca directory: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				cmd := exec.Command("masurca", "-t", "32", "-i", fmt.Sprintf("%s,%s", unmappedFwdReads, unmappedRevReads))
				cmd.Dir = masurcaDir
				fmt.Println(cmd.String())
				err = cmd.Run()
				if err != nil {
					fmt.Printf("Error running masurca: %v\n", err)
					log.Fatalf("%s\tMASURCA\t%v\t%s\tFAILED", sn, i, latestRef)
				}
			}
			log.Printf("%s\tMASURCA\t%v\t%s\tFINISHED", sn, i, latestRef)

			//--------------------------------------- Run MegaHit 2nd --------------------------------------------------- //

			fmt.Printf("Assemble unmapped reads for %s with MEGAHIT only... \n", sn)
			megahitDir := filepath.Join(sampleDir, "MegaHit")
			megahitAssembly := filepath.Join(megahitDir, "final.contigs.fa")

			_, megahitErr := os.Stat(megahitAssembly)
			if megahitErr == nil {
				fmt.Println("Megahit assembly exists for ", sn, "skip")
				//continue
			} else {
				log.Printf("%s\tMEGAHIT\t%v\t%s\tSTARTED", sn, i, latestRef)
				err := os.RemoveAll(megahitDir)
				if err != nil {
					fmt.Printf("Error removing megahit directory: %v\n", err)
					log.Fatalf("%s\tMEGAHIT\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				fmt.Printf("Megahit assembly does not exist for %s \n", sn)
				fmt.Printf("Assemble unmapped reads for %s ... \n", sn)
				//log.Printf("%s\tMEGAHIT\t%v\t%sSTARTED", sn, i, latestRef)
				cmdStr := fmt.Sprintf(`megahit  -1 %s -2 %s -o %s`, unmappedFwdReads, unmappedRevReads, megahitDir)
				fmt.Println(cmdStr)
				mErr := utils.RunBashCmdVerbose(cmdStr)
				if mErr != nil {
					fmt.Printf("Error running megahit: %v\n", mErr)
					log.Fatalf("%s\tMEGAHIT\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				log.Printf("%s\tMEGAHIT\t%v\t%s\tFINISHED", sn, i, latestRef)
			}

			// -------------------------------------------- MAC --------------------------------------------------------- //
			macDir, err := filepath.Abs(filepath.Join(sampleDir, "MAC"))
			if err != nil {
				fmt.Printf("Error determining absolute path for MAC directory: %v\n", err)
				log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
				//return
			}

			macInputDir := filepath.Join(macDir, "input")
			macOutputDir := filepath.Join(macDir, "output")
			macTempDir := filepath.Join(macDir, "temp")
			macAssembly := filepath.Join(macDir, "output", "scaffold.fasta")

			_, macErr := os.Stat(macAssembly)
			if macErr == nil {
				fmt.Printf("MAC assembly exists for %s. Skipping...\n", sn)
				//continue
			} else {
				log.Printf("%s\tMAC\t%v\t%s\tSTARTED", sn, i, latestRef)
				err = os.RemoveAll(macDir)
				if err != nil {
					fmt.Printf("Error removing MAC directory: %v\n", err)
					log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}

				macDirs := []string{macDir, macInputDir, macOutputDir, macTempDir}
				for _, d := range macDirs {
					err = os.MkdirAll(d, os.ModePerm)
					if err != nil {
						fmt.Printf("Error creating directory %s: %v\n", d, err)
						log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
						//return
					}
				}

				err = utils.CopyFile(megahitAssembly, filepath.Join(macInputDir, "megahit_assembly.fa"))
				if err != nil {
					fmt.Printf("Error copying megahit assembly to MAC directory: %v\n", err)
					log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				err = utils.CopyFile(masurcaAssembly, filepath.Join(macInputDir, "masurca_assembly.fa"))
				if err != nil {
					fmt.Printf("Error copying masurca assembly to MAC directory: %v\n", err)
					log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}

				cmdMac := exec.Command("MAC2.0", "megahit_assembly.fa", "masurca_assembly.fa")
				cmdMac.Dir = macDir
				fmt.Println(cmdMac.String())
				err = cmdMac.Run()
				if err != nil {
					fmt.Printf("Error running MAC: %v\n", err)
					log.Fatalf("%s\tMAC\t%v\t%s\tFAILED", sn, i, latestRef)
					//return
				}
				log.Printf("%s\tMAC\t%v\t%s\tFINISHED", sn, i, latestRef)
				assembledContigs = macAssembly
			}

		}

		// ----------------------------------------- TRIM FASTA ----------------------------------------------------- //
		fmt.Printf("Removing contigs less than 200bp from fasta: %s ... \n", assembledContigs)

		trimmedContigs := filepath.Join(sampleDir, sn+".trimmed.fasta")

		_, trimmedContigsErr := os.Stat(trimmedContigs)
		if trimmedContigsErr == nil {
			fmt.Printf("Trimmed contigs exist for %s. Skipping...\n", sn)
		} else {
			log.Printf("%s\tTRIM\t%v\t%s\tSTARTED", sn, i, latestRef)
			trimCmdStr := fmt.Sprintf("seqtk seq -L 200 %s > %s", assembledContigs, trimmedContigs)
			fmt.Println(trimCmdStr)
			trimErr := utils.RunBashCmdVerbose(trimCmdStr)
			if trimErr != nil {
				fmt.Printf("Error trimming contigs: %v\n", trimErr)
				log.Fatalf("%s\tTRIM\t%v\t%sFAILED", sn, i, latestRef)
				return
			}
			log.Printf("%s\tTRIM\t%v\t%sFINISHED", sn, i, latestRef)

		}

		// ---------------------------------------- RENAME CONTIGS -------------------------------------------------- //
		fmt.Printf("Renaming contigs in fasta: %s ... \n", trimmedContigs)
		log.Printf("%s\tRENAME\t%v\t%s\tSTARTED", sn, i, latestRef)
		renamedContigs := filepath.Join(sampleDir, sn+".renamed.fasta")
		rErr := RenameScaffs(trimmedContigs, renamedContigs, sn)
		if rErr != nil {
			fmt.Printf("Error renaming contigs: %v\n", rErr)
			log.Fatalf("%s\tRENAME\t%v\t%s\tFAILED", sn, i, latestRef)
			return
		}
		log.Printf("%s\tRENAME\t%v\t%s\tFINISHED", sn, i, latestRef)

		// ---------------------------------------- Update Reference ------------------------------------------------ //
		log.Printf("%s\tUPDATE_REF\t%v\t%s\tSTARTED", sn, i, latestRef)
		fastas := []string{latestRef, renamedContigs}
		newRef := filepath.Join(sampleDir, fmt.Sprintf("%v_%s_updated_ref.fa", i, sn))
		fErr := ConcatFasta(fastas, newRef)
		if fErr != nil {
			fmt.Printf("Error concatenating fasta files: %v\n", fErr)
			log.Fatalf("%s\tUPDATE_REF\t%v\t%s\tFAILED", sn, i, latestRef)

			return
		}
		latestRef = newRef
		log.Printf("%s\tUPDATE_REF\t%v\t%s\tFINISHED", sn, i, latestRef)
		fmt.Printf("Updated reference: %s\n", latestRef)
		fmt.Println("-------------------------------------------------------")
		i++
	}

}
