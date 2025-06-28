package variants

import (
	"compress/gzip"
	"fmt"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"strings"
	"sync"
)

func VariantCalling(refFile string, bams []string, out string, species string, maxParallelJobs int, verbosity string) {

	// --------------------------------------- Opening fasta file --------------------------------------------------- //
	fmt.Println("Working on FASTA file ...")
	fna, err := os.Open(refFile)
	if err != nil {
		log.Fatalf("Failed to open FASTA file: %v", err)
	}
	defer func(fna *os.File) {
		err := fna.Close()
		if err != nil {
			panic(err)
		}
	}(fna)

	var reader io.Reader = fna
	if strings.HasSuffix(refFile, ".gz") {
		gzReader, err := gzip.NewReader(fna)
		if err != nil {
			log.Fatalf("Failed to create gzip reader: %v", err)
		}
		defer gzReader.Close()
		reader = gzReader
	}

	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")
	logFilePath := filepath.Join(out, "variant_calling.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	jlog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")
	slog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", "ALL")

	//-------------------------- If resuming (read logfile and check for completed stages) -------------------------- //
	logged := utils.ParseLogFile(logFilePath)

	if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", "ALL") {
		fmt.Println("All stages  completed. Skipping.")
		return
	}

	// --------------------------------------- Reading fasta -------------------------------------------------------- //

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	var wg sync.WaitGroup
	sem := make(chan struct{}, maxParallelJobs)
	var jointvSlice []string

	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		chromDir := strings.ReplaceAll(seq.ID, ".", "_")
		chromDirPath := filepath.Join(out, chromDir)
		gvcfPath := filepath.Join(chromDirPath, "gvcfs")
		tmpPath := filepath.Join(chromDirPath, "tmp")
		tmp2Path := filepath.Join(chromDirPath, "tmp2")
		vcfPath := filepath.Join(chromDirPath, "VCFs")

		jointVCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.vcf.gz")
		snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
		indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
		hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		theDB := filepath.Join(chromDirPath, chromDir+"DB")

		//if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
		//	slog.Info(fmt.Sprintf("Chromosome %s already processed all steps. Skipping\n", seq.ID))
		//	jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)
		//	continue
		//}

		sem <- struct{}{}
		wg.Add(1)

		go func(seq *linear.Seq) {
			defer func() {
				wg.Done()
				<-sem
			}()

			fmt.Println(seq.ID)

			slog.Info("Creating directories ...")
			dirsToCreate := []string{chromDirPath, gvcfPath, tmpPath, tmp2Path, vcfPath}
			for _, dir := range dirsToCreate {
				if _, err := os.Stat(dir); os.IsNotExist(err) {
					cErr := os.MkdirAll(dir, 0755)
					if cErr != nil {

						slog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", cErr))
						jlog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", cErr))

						return
					}
				}
			}

			// ------------------------------------ HAPLOTYPE CALLER (Skip completed) ------------------------------- //
			var vSlice []string

			for _, bam := range bams {
				bamName := filepath.Base(bam)
				theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))

				if utils.StageHasCompleted(logged, "HaplotypeCaller", bamName, seq.ID) {
					msg := fmt.Sprintf("HaplotypeCaller already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
					slog.Info(msg)
				} else {
					hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, seq.ID, theGVCF, verbosity)
					vSlice = append(vSlice, "-V "+theGVCF)

					jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
					slog.Info(hapCmdStr)
					hapErr := utils.RunBashCmdVerbose(hapCmdStr)

					if hapErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
						log.Fatalf("FAILED: %v", hapErr)
						return
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED") //, "CMD", hapCmdStr)
					slog.Info("VARIANT CALLING", "CMD", hapCmdStr, "STATUS", "COMPLETED")

				}

			}

			// ---------------------------------- GENOMICS DB IMPORT (Skip completed) ------------------------------- //

			if utils.StageHasCompleted(logged, "GenomicsDBImport", "ALL", seq.ID) {
				msg := fmt.Sprintf("GenomicsDBImport already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)
			} else {
				vArgs := strings.Join(vSlice, " ")

				//----------------------------- Delete DB if present and delete ------------------------------------- //

				dErr := os.RemoveAll(theDB)
				if dErr != nil {
					fmt.Println("Error removing directory:", dErr)
					slog.Error("VARIANT CALLING", "PROGRAM", "rm ", "SAMPLE", theDB, "CHROMOSOME", seq.ID, "STATUS", "ERROR", "CMD", fmt.Sprintf("%v", dErr))
					log.Fatalf("Error removing directory: %v", dErr)
					return
				} else {
					fmt.Println("Directory removed successfully (if it existed).")
				}

				// -------------------------------------------------------------------------------------------------- //

				gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader`, vArgs, theDB, tmpPath, seq.ID)

				jlog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
				slog.Info(gDBImpCmdStr)

				gErr := utils.RunBashCmdVerbose(gDBImpCmdStr)
				if gErr != nil {

					jlog.Error("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", gErr))
					log.Fatalf("FAILED: %v", gErr)
					return
				}

				jlog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", gDBImpCmdStr, "STATUS", "COMPLETED")

			}

			// --------------------------------------- GENOTYPE GVCFS (Skip completed) ------------------------------ //

			if utils.StageHasCompleted(logged, "GenotypeGVCFs", "ALL", seq.ID) {
				msg := fmt.Sprintf("GenotypeGVCFs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)

			} else {
				genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s`, refFile, theDB, jointVCF, tmpPath)

				jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
				slog.Info(genoCmdStr)

				gtErr := utils.RunBashCmdVerbose(genoCmdStr)
				if gtErr != nil {

					jlog.Error("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
					log.Fatalf("FAILED: %v", gtErr)
					return
				}
				jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", genoCmdStr, "STATUS", "COMPLETED")

			}

			// -------------------------------------- SELECT SNPs (Skip completed) ---------------------------------- //

			fmt.Println("Hard filtered  joint VCF ...")

			if utils.StageHasCompleted(logged, "SelectSNPs", "ALL", seq.ID) {
				msg := fmt.Sprintf("SELECT_SNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)

			} else {

				jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
				slog.Info("gatk SelectVariants --select-type-to-include SNP")

				sErr := GetVariantType(jointVCF, "SNP")
				if sErr != nil {

					jlog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", sErr))
					log.Fatalf("FAILED: %v", sErr)
					return
				}
				jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", "gatk SelectVariants --select-type-to-include SNP", "STATUS", "COMPLETED")

			}

			// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //
			if utils.StageHasCompleted(logged, "SelectINDELs", "ALL", seq.ID) {
				msg := fmt.Sprintf("SELECT_INDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)

			} else {

				jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
				slog.Info("gatk SelectVariants --select-type-to-include INDEL")

				iErr := GetVariantType(jointVCF, "INDEL")
				if iErr != nil {
					jlog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", iErr))
					log.Fatalf("FAILED: %v", iErr)
					return
				}
				jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", "gatk SelectVariants --select-type-to-include INDEL", "STATUS", "COMPLETED")

			}

			// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

			//if _, exists := completed["HardFilteringSNPS"][seq.ID]; exists {
			//	fmt.Printf("HardFilteringSNPS already completed for %s. Skipping.\n", seq.ID)
			if utils.StageHasCompleted(logged, "HardFilteringSNPS", "ALL", seq.ID) {
				msg := fmt.Sprintf("HardFilteringSNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)

			} else {
				//log.Printf("%s\tHardFilteringSNPS\tALL\tSTARTED\tSNPs\n", seq.ID)
				jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
				slog.Info("gatk HardFilteringSNPS ")

				hsErr := HardFilterSNPs(snpVCF)
				if hsErr != nil {
					//log.Printf("%s\tHardFilteringSNPS\tALL\tFAILED\tSNPs - Error: %v\n", seq.ID, hsErr)
					jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
					log.Fatalf("FAILED: %v", hsErr)
					return
				}
				//log.Printf("%s\tHardFilteringSNPS\tALL\tFINISHED\tSNPs\n", seq.ID)
				jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", "gatk HardFilteringSNPS ", "STATUS", "COMPLETED")

			}

			// -------------------------------- HARD FILTERING INDELS (Skip completed) ------------------------------ //

			//if _, exists := completed["HardFilteringINDELS"][seq.ID]; exists {
			//	fmt.Printf("HardFilteringINDELS already completed for %s. Skipping.\n", seq.ID)
			if utils.StageHasCompleted(logged, "HardFilteringINDELs", "ALL", seq.ID) {
				msg := fmt.Sprintf("HardFilteringINDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)

			} else {
				//log.Printf("%s\tHardFilteringINDELS\tALL\tSTARTED\tINDELs\n", seq.ID)
				jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
				slog.Info("gatk HardFilteringINDELS ")

				hiErr := HardFilterINDELs(indelVCF)
				if hiErr != nil {
					//log.Printf("%s\tHardFilteringINDELS\tALL\tFAILED\tINDELs - Error: %v\n", seq.ID, hiErr)
					jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
					log.Fatalf("FAILED: %v", hiErr)
					return
				}
				//log.Printf("%s\tHardFilteringINDELS\tALL\tFINISHED\tINDELs\n", seq.ID)

				jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("CMD", "gatk HardFilteringINDELS ", "STATUS", "COMPLETED")

			}

			// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

			//if _, exists := completed["MergeVcfs"][seq.ID]; exists {
			//	fmt.Printf("MergeVcfs already completed for %s. Skipping.\n", seq.ID)
			if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
				msg := fmt.Sprintf("MergeVcfs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
				slog.Info(msg)
			} else {
				jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

				mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
				slog.Info(mergeCmdStr)
				//log.Printf("%s\tMergeVcfs\tALL\tSTART\t%s\n", seq.ID, mergeCmdStr)
				//fmt.Println(mergeCmdStr)
				mErr := utils.RunBashCmdVerbose(mergeCmdStr)
				if mErr != nil {
					//log.Printf("%s\tMergeVcfs\tALL\tFAILED\t%s - Error: %v\n", seq.ID, mergeCmdStr, mErr)
					jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
					slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
					log.Fatalf("FAILED: %v", mErr)
					return

				}
				//log.Printf("%s\tMergeVcfs\tALL\tFINISHED\t%s\n", seq.ID, mergeCmdStr)
				jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				slog.Info("VARIANT CALLING", "STATUS", "COMPLETED")
				//jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)

			}
			jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)

		}(seq)

	}
	wg.Wait()
	//fmt.Println("All jobs completed.")
	//slog.Info("VARIANT CALLING", "STATUS", "COMPLETED")

	fmt.Println("Merging ALL Hard filtered VCFs ...")

	finalVcf := filepath.Join(out, species+".joint_hard_filtered.vcf.gz")
	mergeCmdStr2 := fmt.Sprintf(`gatk MergeVcfs %s -O %s`, strings.Join(jointvSlice, " "), finalVcf)
	slog.Info(mergeCmdStr2)
	jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

	mErr := utils.RunBashCmdVerbose(mergeCmdStr2)
	if mErr != nil {
		jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		return
	}

	jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("VARIANT CALLING", "STATUS", "COMPLETED")

}

func VariantCallingConfig(configFile string, species string, maxParallelJobs int, verbosity string) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Bams", cfg.Bams)
	fmt.Println("Ouput", cfg.OutputDir)

	refFile := cfg.Reference

	_, rErr := os.Stat(refFile)
	if rErr != nil {
		fmt.Printf("Reference file: %s does not exist\n", refFile)
		return
	}
	bams := cfg.Bams

	fmt.Printf("bams: %v\n", bams)
	if len(bams) == 0 {
		fmt.Println("You must provide at least one bam file")
		return
	} else {
		for i := range bams {
			_, err := os.Stat(bams[i])
			if err != nil {
				fmt.Printf("Bam file: %s is not a valid file path: %v\n", bams[i], err)
				return
			}
		}
	}

	outputDir := cfg.OutputDir

	outInfo, outErr := os.Stat(outputDir)

	if outErr != nil {

		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory: %s does not exist. Attempting to create it.\n", outputDir)
			if createErr := os.MkdirAll(outputDir, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", outputDir, createErr)
				return
			}
			fmt.Printf("Output directory %s created successfully.\n", outputDir)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", outputDir, outErr)
			return
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output Directory %s file path is not a directory\n", outputDir)
		return
	}

	VariantCalling(refFile, bams, outputDir, species, maxParallelJobs, verbosity)
}
