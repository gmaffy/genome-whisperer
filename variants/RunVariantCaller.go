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

func VariantCalling(refFile string, bams []string, out string, species string, maxParallelJobs int, verbosity string, caller string, merger string, logFilePath string, dvVer string, modelType string) {

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
		jointBCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.bcf")
		snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
		indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
		hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		theDB := filepath.Join(chromDirPath, chromDir+"DB")

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

		sem <- struct{}{}
		wg.Add(1)

		if caller == "gatk" {

			// ------------------------------------------------ Run ------------------------------------------------- //
			go func(seq *linear.Seq) {
				defer func() {
					wg.Done()
					<-sem
				}()

				slog.Info(fmt.Sprintf("Working on chromosome %s................\n\n", seq.ID))

				// -------------------------------- HAPLOTYPE CALLER (Skip completed) ------------------------------- //
				var vSlice []string

				for _, bam := range bams {
					bamName := filepath.Base(bam)
					theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
					vSlice = append(vSlice, "-V "+theGVCF)

					if utils.StageHasCompleted(logged, "HaplotypeCaller", bamName, seq.ID) {
						msg := fmt.Sprintf("HaplotypeCaller already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
						slog.Info(msg)

					} else {
						hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, seq.ID, theGVCF, verbosity)

						jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
						//slog.Info(hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
						hapErr := utils.RunBashCmdVerbose(hapCmdStr)

						if hapErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							log.Fatalf("FAILED: %v", hapErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED") //, "CMD", hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					}

				}

				// ---------------------------------- MERGING (Skip completed) ------------------------------- //

				if merger == "gatk" {

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
						} else {
							fmt.Println("Directory removed successfully (if it existed).")
						}

						// -------------------------------------------------------------------------------------------------- //

						gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader --verbosity %s`, vArgs, theDB, tmpPath, seq.ID, verbosity)

						jlog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

						gErr := utils.RunBashCmdVerbose(gDBImpCmdStr)
						if gErr != nil {

							jlog.Error("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", gErr))
							log.Fatalf("FAILED: %v", gErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					}

					// --------------------------------------- GENOTYPE GVCFS (Skip completed) ------------------------------ //

					if utils.StageHasCompleted(logged, "GenotypeGVCFs", "ALL", seq.ID) {
						msg := fmt.Sprintf("GenotypeGVCFs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`, refFile, theDB, jointVCF, tmpPath, verbosity)

						jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

						gtErr := utils.RunBashCmdVerbose(genoCmdStr)
						if gtErr != nil {

							jlog.Error("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
							log.Fatalf("FAILED: %v", gtErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					}

				} else if merger == "glnexus" {

					if utils.StageHasCompleted(logged, "GLNEXUS", "ALL", seq.ID) {
						msg := fmt.Sprintf("glnexus_cli already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						glnexusCmdStr := fmt.Sprintf(`glnexus_cli --config gatk  %s/*.vcf.gz > %s`, gvcfPath, jointBCF)
						slog.Info(glnexusCmdStr)
						glErr := utils.RunBashCmdVerbose(glnexusCmdStr)
						if glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							log.Fatalf("FAILED: %v", glErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("CMD", glnexusCmdStr, "STATUS", "COMPLETED")
					}

					if utils.StageHasCompleted(logged, "GLNEXUS_BCFTOOLS", "ALL", seq.ID) {
						msg := fmt.Sprintf("glnexus_cli bcftools already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						bcftoolsCmdStr := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
						slog.Info(bcftoolsCmdStr)
						glErr := utils.RunBashCmdVerbose(bcftoolsCmdStr)
						if glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							log.Fatalf("FAILED: %v", glErr)
						}

						rmErr := os.RemoveAll("GLnexus.DB")
						if rmErr != nil {
							fmt.Println("Error removing directory:", rmErr)
							//slog.Error("VARIANT CALLING", "PROGRAM", "rm ", "SAMPLE", "GLnexus.DB", "CHROMOSOME", seq.ID, "STATUS", "ERROR", "CMD", fmt.Sprintf("%v", rmErr))
							log.Fatalf("Error removing directory: %v", rmErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("CMD", bcftoolsCmdStr, "STATUS", "COMPLETED")
					}

				} else {
					fmt.Println("Merger should either be gatk or glnexus")
					os.Exit(1)
				}

				// -------------------------------------- SELECT SNPs (Skip completed) ---------------------------------- //

				fmt.Println("Hard filtered  joint VCF ...")

				if utils.StageHasCompleted(logged, "SelectSNPs", "ALL", seq.ID) {
					msg := fmt.Sprintf("SELECT_SNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					sErr := GetVariantType(jointVCF, "SNP")
					if sErr != nil {

						jlog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
						slog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
						log.Fatalf("FAILED: %v", sErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //
				if utils.StageHasCompleted(logged, "SelectINDELs", "ALL", seq.ID) {
					msg := fmt.Sprintf("SELECT_INDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					iErr := GetVariantType(jointVCF, "INDEL")
					if iErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
						slog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
						log.Fatalf("FAILED: %v", iErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

				if utils.StageHasCompleted(logged, "HardFilteringSNPS", "ALL", seq.ID) {
					msg := fmt.Sprintf("HardFilteringSNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {
					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					hsErr := HardFilterSNPs(snpVCF, verbosity)
					if hsErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
						slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
						log.Fatalf("FAILED: %v", hsErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------- HARD FILTERING INDELS (Skip completed) ------------------------------ //

				if utils.StageHasCompleted(logged, "HardFilteringINDELS", "ALL", seq.ID) {
					msg := fmt.Sprintf("HardFilteringINDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					hiErr := HardFilterINDELs(indelVCF, verbosity)
					if hiErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
						slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
						log.Fatalf("FAILED: %v", hiErr)
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

				if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
					msg := fmt.Sprintf("MergeVcfs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)
				} else {
					jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
					mErr := utils.RunBashCmdVerbose(mergeCmdStr)
					if mErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
						log.Fatalf("FAILED: %v", mErr)
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}
				jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)

			}(seq)
		} else if caller == "deepvariant" {

			// --------------------------------------- Check Deps -------------------------------------------------- //

			go func(seq *linear.Seq) {
				defer func() {
					wg.Done()
					<-sem
				}()

				fmt.Println(seq.ID)

				// ------------------------------------ DEEP VARIANT CALLING (Skip completed) ------------------------------- //

				for _, bam := range bams {

					refDir := filepath.Dir(refFile)
					refName := filepath.Base(refFile)
					bamDir := filepath.Dir(bam)
					bamName := filepath.Base(bam)
					gvcfName := strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1)
					vcfName := strings.Replace(bamName, "bam", fmt.Sprintf("%s.vcf.gz", chromDir), 1)

					if utils.StageHasCompleted(logged, "DEEPVARIANT", bamName, seq.ID) {
						msg := fmt.Sprintf("DEEPVARIANT already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						dvCmdStr := fmt.Sprintf(`docker run -v "%s":"/bam" -v "%s":"/ref" -v "%s":"/output" google/deepvariant:"%s" /opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s --regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s --intermediate_results_dir /output/tmp`, bamDir, refDir, gvcfPath, dvVer, modelType, refName, bamName, seq.ID, vcfName, gvcfName)
						slog.Info(dvCmdStr)
						dvErr := utils.RunBashCmdVerbose(dvCmdStr)
						if dvErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							slog.Error(dvCmdStr, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							log.Fatalf("FAILED: %v", dvErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "STATUS", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "COMPLETED")
					}
				}

				if merger == "glnexus" {

					if utils.StageHasCompleted(logged, "GLNEXUS", "ALL", seq.ID) {
						msg := fmt.Sprintf("glnexus_cli already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						glnexusCmdStr := fmt.Sprintf(`glnexus_cli --config gatk  %s/*.vcf.gz > %s`, gvcfPath, jointBCF)
						slog.Info(glnexusCmdStr)
						glErr := utils.RunBashCmdVerbose(glnexusCmdStr)
						if glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							log.Fatalf("FAILED: %v", glErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("CMD", glnexusCmdStr, "STATUS", "COMPLETED")
					}

					if utils.StageHasCompleted(logged, "GLNEXUS_BCFTOOLS", "ALL", seq.ID) {
						msg := fmt.Sprintf("glnexus_cli bcftools already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						bcftoolsCmdStr := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
						slog.Info(bcftoolsCmdStr)
						glErr := utils.RunBashCmdVerbose(bcftoolsCmdStr)
						if glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							log.Fatalf("FAILED: %v", glErr)
						}

						rmErr := os.RemoveAll("GLnexus.DB")
						if rmErr != nil {
							fmt.Println("Error removing directory:", rmErr)
							//slog.Error("VARIANT CALLING", "PROGRAM", "rm ", "SAMPLE", "GLnexus.DB", "CHROMOSOME", seq.ID, "STATUS", "ERROR", "CMD", fmt.Sprintf("%v", rmErr))
							log.Fatalf("Error removing directory: %v", rmErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("CMD", bcftoolsCmdStr, "STATUS", "COMPLETED")
					}

				} else {
					fmt.Println("Merger can only be glnexus if deep variant is the variant caller")
					os.Exit(1)
				}
				// -------------------------------------- SELECT SNPs (Skip completed) ---------------------------------- //

				fmt.Println("Hard filtered  joint VCF ...")

				if utils.StageHasCompleted(logged, "SelectSNPs", "ALL", seq.ID) {
					msg := fmt.Sprintf("SELECT_SNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED") //, "CMD", hapCmdStr)
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					sErr := GetVariantType(jointVCF, "SNP")
					if sErr != nil {

						jlog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", sErr))
						log.Fatalf("FAILED: %v", sErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //
				if utils.StageHasCompleted(logged, "SelectINDELs", "ALL", seq.ID) {
					msg := fmt.Sprintf("SELECT_INDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					iErr := GetVariantType(jointVCF, "INDEL")
					if iErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", iErr))
						log.Fatalf("FAILED: %v", iErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

				if utils.StageHasCompleted(logged, "HardFilteringSNPS", "ALL", seq.ID) {
					msg := fmt.Sprintf("HardFilteringSNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {
					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					hsErr := HardFilterSNPs(snpVCF, verbosity)
					if hsErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
						log.Fatalf("FAILED: %v", hsErr)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------- HARD FILTERING INDELS (Skip completed) ------------------------------ //

				if utils.StageHasCompleted(logged, "HardFilteringINDELS", "ALL", seq.ID) {
					msg := fmt.Sprintf("HardFilteringINDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)

				} else {

					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					//fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

					hiErr := HardFilterINDELs(indelVCF, verbosity)
					if hiErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
						log.Fatalf("FAILED: %v", hiErr)
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

				}

				// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

				if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
					msg := fmt.Sprintf("MergeVcfs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
					slog.Info(msg)
				} else {
					jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

					mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s --verbosity %s`, snpVCF, indelVCF, hardFilteredVCF, verbosity)
					//slog.Info(mergeCmdStr)
					mErr := utils.RunBashCmdVerbose(mergeCmdStr)
					if mErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
						slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
						log.Fatalf("FAILED: %v", mErr)
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					//slog.Info("VARIANT CALLING", "STATUS", "COMPLETED")

				}
				jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)

			}(seq)

		}
	}
	wg.Wait()

	fmt.Println("Merging ALL Hard filtered VCFs ...")

	finalVcf := filepath.Join(out, species+".joint_hard_filtered.vcf.gz")
	mergeCmdStr2 := fmt.Sprintf(`gatk MergeVcfs %s -O %s`, strings.Join(jointvSlice, " "), finalVcf)
	//slog.Info(mergeCmdStr2)
	jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

	mErr := utils.RunBashCmdVerbose(mergeCmdStr2)
	if mErr != nil {
		jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		log.Fatalf("FAILED: %v", mErr)
	}

	jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	print("-------------------------------------------------------------------------------------------\n\n")
	//slog.Info("ALL VARIANT CALLING STEPS", "STATUS", "COMPLETED")

}

func VariantCallingConfig(configFile string, species string, maxParallelJobs int, verbosity string, caller string, merger string, dvVer string, modelType string) {
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
	logFilePath := filepath.Join(outputDir, "variant_calling.log")
	VariantCalling(refFile, bams, outputDir, species, maxParallelJobs, verbosity, caller, merger, logFilePath, dvVer, modelType)
}
