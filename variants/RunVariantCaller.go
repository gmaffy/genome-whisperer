package variants

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

type SeqInfo struct {
	ID  string
	Len int
}

func getChromsAndContigs(dictFilePath string) ([]SeqInfo, []SeqInfo, error) {
	dictFile, err := os.Open(dictFilePath)
	if err != nil {
		fmt.Printf("Error opening reference dict file: %s: %v\n", dictFilePath, err)
		return nil, nil, err
	}
	defer dictFile.Close()

	scanner := bufio.NewScanner(dictFile)
	var seqs []SeqInfo

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "@SQ") {
			parts := strings.Split(line, "\t")
			seqID := strings.TrimPrefix(parts[1], "SN:")
			seqLenStr := strings.TrimPrefix(parts[2], "LN:")
			seqLen, _ := strconv.Atoi(seqLenStr)

			seqs = append(seqs, SeqInfo{ID: seqID, Len: seqLen})
		}
	}

	// Sort by length descending
	sort.Slice(seqs, func(i, j int) bool {
		return seqs[i].Len > seqs[j].Len
	})

	// Split into chroms and contigs
	var chroms, contigs []SeqInfo
	if len(seqs) > 21 {
		chroms = append(chroms, seqs[:21]...)
		contigs = append(contigs, seqs[21:]...)
	} else {
		chroms = append(chroms, seqs...)
	}

	// Ensure MT and Pltd are included in chroms
	for i := 0; i < len(contigs); {
		if contigs[i].ID == "MT" || contigs[i].ID == "Pltd" {
			chroms = append(chroms, contigs[i])
			contigs = append(contigs[:i], contigs[i+1:]...)
		} else {
			i++
		}
	}

	return chroms, contigs, nil
}

func FindBamOrVcfs(dataDirAbs string, species string, sample string, refVer string, bamOrVcf string, chrom string) (error, []string) {
	var pattern string
	switch bamOrVcf {
	case "bams":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*bqsr.cram", dataDirAbs, strings.ToLower(species), sample, refVer)
	case "gvcfs":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/gvcfs/*%s.g.vcf.gz", dataDirAbs, strings.ToLower(species), sample, refVer, chrom)
	}

	matches, err := filepath.Glob(pattern)
	if err != nil {
		fmt.Println("glob error:", err)
		return err, []string{}
	}

	if len(matches) == 0 {
		//color.Red("No %s files found in %s. Skipping ....\n\n", bamOrVcf, sample)
		return fmt.Errorf("no %s files found in %s ", bamOrVcf, sample), matches
	}
	if len(matches) > 1 {
		//color.Red("Multiple %s files found in %s: Skipping ....\n\n", bamOrVcf, sample)
		return fmt.Errorf("multiple %s files found in %s", bamOrVcf, sample), matches
	}

	//color.Green("Found %s file for sample %s: %s\n\n", bamOrVcf, color.BlueString(sample), matches[0])
	return nil, matches
}

func CreateGvcf(bam string, refFile string, chroms []SeqInfo, theGVCF string, gatkLogLevel string, caller string, dvVer string, modelType string, verbose bool) (string, error) {
	var sID string
	var err error
	switch caller {
	case "gatk":
		if len(chroms) == 1 {
			sID = chroms[0].ID
		} else {
			sID = filepath.Join(os.TempDir(), "gatk_intervals_contigs.list")
			f, err1 := os.Create(sID)
			if err1 != nil {
				return "", err1
			}
			defer f.Close()

			for _, chrom := range chroms {
				fmt.Fprintln(f, chrom.ID)
			}
		}
		hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, sID, theGVCF, gatkLogLevel)

		color.Green("\n---------------------------------------------------------------------------------\n%s\n\n----------------------------------------------------------------------------------------\n\n", hapCmdStr)
		if verbose {
			err = utils.RunBashCmdVerbose(hapCmdStr)
		} else {
			err = utils.RunBashCmd(hapCmdStr)
		}

	case "DeepVariant":
		if len(chroms) == 1 {
			sID = chroms[0].ID
		} else {
			sID = filepath.Join(os.TempDir(), "deepvariant_intervals_contigs.list")
			f, err1 := os.Create(sID)
			if err1 != nil {
				return "", err1
			}
			defer f.Close()

			for _, chrom := range chroms {
				fmt.Fprintf(f, "%s\t0\t%d\n", chrom.ID, chrom.Len)
			}
		}
		bamDir := filepath.Dir(bam)
		bamName := filepath.Base(bam)
		refName := filepath.Base(refFile)
		refDir := filepath.Dir(refFile)
		gvcfPath := filepath.Dir(theGVCF)
		gvcfName := filepath.Base(theGVCF)
		vcfName := strings.Replace(gvcfName, ".g.vcf.gz", ".vcf.gz", 1)
		dvCmdStr := fmt.Sprintf(`docker run -v "%s":"/bam" -v "%s":"/ref" -v "%s":"/output" google/deepvariant:"%s" /opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s --regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s --intermediate_results_dir /output/tmp`, bamDir, refDir, gvcfPath, dvVer, modelType, refName, bamName, sID, vcfName, gvcfName)

		if verbose {
			err = utils.RunBashCmdVerbose(dvCmdStr)
		} else {
			err = utils.RunBashCmd(dvCmdStr)
		}
	}

	return theGVCF, err
}

func VariantCalling(refFile string, bams []string, out string, species string, threadsPerSample int, gatkLogLevel string, caller string, merger string, logFilePath string, dvVer string, modelType string, verbose bool, noMerging bool) (string, error) {
	// --------------------------------------- Opening fasta file --------------------------------------------------- //
	fmt.Println("Working on FASTA file ...")
	fna, err := os.Open(refFile)
	if err != nil {
		return "", fmt.Errorf("failed to open FASTA file: %v", err)
	}
	defer func(fna *os.File) {
		err := fna.Close()
		if err != nil {
			slog.Error("VARIANT CALLING", "PROGRAM", "FileClose", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
			panic(err)
		}
	}(fna)

	var reader io.Reader = fna
	if strings.HasSuffix(refFile, ".gz") {
		gzReader, err := gzip.NewReader(fna)
		if err != nil {
			return "", fmt.Errorf("failed to create gzip reader: %v", err)
		}
		defer gzReader.Close()
		reader = gzReader
	}

	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")

	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		return "", fmt.Errorf("failed to open log file %s - %s", logFilePath, err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	finalVcf := filepath.Join(out, species+".joint_hard_filtered.vcf.gz")

	jlog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

	//-------------------------- If resuming (read logfile and check for completed stages) -------------------------- //
	logged := utils.ParseLogFile(logFilePath)

	if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", "ALL") {
		fmt.Println("All stages  completed. Skipping.")
		return finalVcf, nil
	}

	// --------------------------------------- Reading fasta -------------------------------------------------------- //

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
	}

	var wg sync.WaitGroup
	var mu sync.Mutex
	// Mutex to ensure GLnexus (and its follow-ups) run strictly one-at-a-time across chromosomes
	var glnxMu sync.Mutex
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
		hardFilteredSNPs := strings.TrimSuffix(snpVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		hardFilteredINDELs := strings.TrimSuffix(indelVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		theDB := filepath.Join(chromDirPath, chromDir+"DB")

		dirsToCreate := []string{chromDirPath, gvcfPath, tmpPath, tmp2Path, vcfPath}
		for _, dir := range dirsToCreate {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				cErr := os.MkdirAll(dir, 0755)
				if cErr != nil {
					slog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", cErr))
					jlog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", cErr))
					return "", fmt.Errorf("failed to created directory %s - %s", dir, cErr)
				}
			}
		}

		sem <- struct{}{}
		wg.Add(1)

		switch caller {
		case "gatk":

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
					theGVCF := ""
					if strings.HasSuffix(bamName, ".bam") {
						theGVCF = filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
					} else if strings.HasSuffix(bamName, ".cram") {
						theGVCF = filepath.Join(gvcfPath, strings.Replace(bamName, "cram", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
					} else {
						fmt.Println("BAM file should have either .bam or .cram extension")
						os.Exit(1)
					}

					vSlice = append(vSlice, "-V "+theGVCF)

					if utils.StageHasCompleted(logged, "HaplotypeCaller", bamName, seq.ID) {
						msg := fmt.Sprintf("HaplotypeCaller already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
						slog.Info(msg)

					} else {
						hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, seq.ID, theGVCF, gatkLogLevel)

						jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						var hapErr error
						if verbose {
							hapErr = utils.RunBashCmdVerbose(hapCmdStr)
						} else {
							hapErr = utils.RunBashCmd(hapCmdStr)
						}

						if hapErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							log.Fatalf("FAILED: %v", hapErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HaplotypeCaller", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}
				}

				// ---------------------------------- MERGING (Skip completed) ------------------------------- //

				if !noMerging {
					switch merger {
					case "gatk":

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

							// ------------------------------------------------------------------------------------------ //
							gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader --verbosity %s`, vArgs, theDB, tmpPath, seq.ID, gatkLogLevel)

							jlog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GenomicsDBImport", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var gErr error
							if verbose {
								gErr = utils.RunBashCmdVerbose(gDBImpCmdStr)
							} else {
								gErr = utils.RunBashCmd(gDBImpCmdStr)
							}

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
							genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`, refFile, theDB, jointVCF, tmpPath, gatkLogLevel)

							jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var gtErr error
							if verbose {
								gtErr = utils.RunBashCmdVerbose(genoCmdStr)
							} else {
								gtErr = utils.RunBashCmd(genoCmdStr)
							}

							if gtErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
								log.Fatalf("FAILED: %v", gtErr)
							}
							jlog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GenotypeGVCFs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
						}

					case "glnexus":

						if utils.StageHasCompleted(logged, "GLNEXUS", "ALL", seq.ID) {
							msg := fmt.Sprintf("glnexus_cli already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
							slog.Info(msg)
						} else {
							glnxMu.Lock()
							glnexusCmdStr := fmt.Sprintf(`glnexus_cli --config gatk  %s/*.vcf.gz > %s`, gvcfPath, jointBCF)
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", glnexusCmdStr)
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
							rmErr := os.RemoveAll("GLnexus.DB")
							if rmErr != nil {
								fmt.Println("Error removing the GLnexus.DB directory:", rmErr)
							}
							var glErr error
							if verbose {
								glErr = utils.RunBashCmdVerbose(glnexusCmdStr)
							} else {
								glErr = utils.RunBashCmd(glnexusCmdStr)
							}

							rErr := os.RemoveAll("GLnexus.DB")
							if rErr != nil {
								fmt.Println("Error removing the GLnexus.DB directory:", rErr)
								log.Fatalf("Error removing directory: %v", rErr)
							}
							glnxMu.Unlock()

							if glErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								log.Fatalf("FAILED: %v", glErr)
							}
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						}

						if utils.StageHasCompleted(logged, "GLNEXUS_BCFTOOLS", "ALL", seq.ID) {
							msg := fmt.Sprintf("glnexus_cli bcftools already completed for %s. Skipping.\n\n---------------------------------------------------------\n\n", seq.ID)
							slog.Info(msg)
						} else {
							bcftoolsCmdStr := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", bcftoolsCmdStr)
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var glErr error
							if verbose {
								glErr = utils.RunBashCmdVerbose(bcftoolsCmdStr)
							} else {
								glErr = utils.RunBashCmd(bcftoolsCmdStr)
							}

							if glErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								log.Fatalf("FAILED: %v", glErr)
							}

							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						}

						if utils.StageHasCompleted(logged, "GLNEXUS_INDEX", "ALL", seq.ID) {
							msg := fmt.Sprintf("VCF INDEXING already completed for %s. Skipping.\n\n---------------------------------------------------------\n\n", seq.ID)
							slog.Info(msg)
						} else {
							indexCmdStr := fmt.Sprintf(`gatk IndexFeatureFile -I %s`, jointVCF)
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", indexCmdStr)
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var indErr error
							if verbose {
								indErr = utils.RunBashCmdVerbose(indexCmdStr)
							} else {
								indErr = utils.RunBashCmd(indexCmdStr)
							}

							if indErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", indErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", indErr))
							}
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_INDEX", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						}

					default:
						fmt.Println("Merger should either be gatk or glnexus")
						os.Exit(1)
					}

					// -------------------------------------- SELECT SNPs (Skip completed) ---------------------------------- //

					fmt.Println("Hard filtered  joint VCF ...")

					if utils.StageHasCompleted(logged, "SelectSNPs", "ALL", seq.ID) {
						msg := fmt.Sprintf("SELECT_SNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {

						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newSnpVCF, sErr := GetVariantType(jointVCF, "SNP", gatkLogLevel, verbose)
						if sErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
							log.Fatalf("FAILED: %v", sErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						snpVCF = newSnpVCF
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //
					if utils.StageHasCompleted(logged, "SelectINDELs", "ALL", seq.ID) {
						msg := fmt.Sprintf("SELECT_INDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newIndelVCF, iErr := GetVariantType(jointVCF, "INDEL", gatkLogLevel, verbose)
						if iErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
							log.Fatalf("FAILED: %v", iErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						indelVCF = newIndelVCF
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

					if utils.StageHasCompleted(logged, "HardFilteringSNPS", "ALL", seq.ID) {
						msg := fmt.Sprintf("HardFilteringSNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newHardFilteredSNPs, hsErr := HardFilterSNPs(snpVCF, gatkLogLevel, verbose)
						if hsErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
							log.Fatalf("FAILED: %v", hsErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						hardFilteredSNPs = newHardFilteredSNPs
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// -------------------------------- HARD FILTERING INDELS (Skip completed) ------------------------------ //

					if utils.StageHasCompleted(logged, "HardFilteringINDELS", "ALL", seq.ID) {
						msg := fmt.Sprintf("HardFilteringINDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newHardFilteredINDELs, hiErr := HardFilterINDELs(indelVCF, gatkLogLevel, verbose)
						if hiErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
							log.Fatalf("FAILED: %v", hiErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						hardFilteredINDELs = newHardFilteredINDELs
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

					if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
						msg := fmt.Sprintf("MergeVcfs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s --VERBOSITY %s`, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF, gatkLogLevel)
						var mErr error
						if verbose {
							mErr = utils.RunBashCmdVerbose(mergeCmdStr)
						} else {
							mErr = utils.RunBashCmd(mergeCmdStr)
						}

						if mErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
							log.Fatalf("FAILED: %v", mErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					}
					mu.Lock()
					jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)
					mu.Unlock()
				}
			}(seq)

		case "DeepVariant":
			// ---------------------------------------- Check Deps -------------------------------------------------- //
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

					gvcfName := ""
					vcfName := ""
					if strings.HasSuffix(bamName, ".bam") {
						gvcfName = strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1)
						vcfName = strings.Replace(bamName, "bam", fmt.Sprintf("%s.vcf.gz", chromDir), 1)
					} else if strings.HasSuffix(bamName, ".cram") {
						gvcfName = strings.Replace(bamName, "cram", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1)
						vcfName = strings.Replace(bamName, "cram", fmt.Sprintf("%s.vcf.gz", chromDir), 1)
					} else {
						fmt.Println("BAM file should have either .bam or .cram extension")
						os.Exit(1)
					}

					if utils.StageHasCompleted(logged, "DEEPVARIANT", bamName, seq.ID) {
						msg := fmt.Sprintf("DEEPVARIANT already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						dvCmdStr := fmt.Sprintf(`docker run -v "%s":"/bam" -v "%s":"/ref" -v "%s":"/output" google/deepvariant:"%s" /opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s --regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s --intermediate_results_dir /output/tmp`, bamDir, refDir, gvcfPath, dvVer, modelType, refName, bamName, seq.ID, vcfName, gvcfName)

						var dvErr error
						if verbose {
							dvErr = utils.RunBashCmdVerbose(dvCmdStr)
						} else {
							dvErr = utils.RunBashCmd(dvCmdStr)
						}

						if dvErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							log.Fatalf("FAILED: %v", dvErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "DEEPVARIANT", "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					}
				}

				if !noMerging {
					if merger == "glnexus" {

						if utils.StageHasCompleted(logged, "GLNEXUS", "ALL", seq.ID) {
							msg := fmt.Sprintf("glnexus_cli already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
							slog.Info(msg)

						} else {
							glnxMu.Lock()
							glnexusCmdStr := fmt.Sprintf(`glnexus_cli --config %s  %s/*.vcf.gz > %s`, caller, gvcfPath, jointBCF)
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", glnexusCmdStr)
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var glErr error
							if verbose {
								glErr = utils.RunBashCmdVerbose(glnexusCmdStr)
							} else {
								glErr = utils.RunBashCmd(glnexusCmdStr)
							}
							glnxMu.Unlock()

							if glErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								log.Fatalf("FAILED: %v", glErr)
							}
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						}

						if utils.StageHasCompleted(logged, "GLNEXUS_BCFTOOLS", "ALL", seq.ID) {
							msg := fmt.Sprintf("glnexus_cli bcftools already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
							slog.Info(msg)

						} else {
							bcftoolsCmdStr := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", bcftoolsCmdStr)
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

							var glErr error
							if verbose {
								glErr = utils.RunBashCmdVerbose(bcftoolsCmdStr)
							} else {
								glErr = utils.RunBashCmd(bcftoolsCmdStr)
							}

							if glErr != nil {
								jlog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								slog.Error("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
								log.Fatalf("FAILED: %v", glErr)
							}

							rmErr := os.RemoveAll("GLnexus.DB")
							if rmErr != nil {
								fmt.Println("Error removing GLnexus.DB directory:", rmErr)
								log.Fatalf("Error removing directory: %v", rmErr)
							}
							jlog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
							slog.Info("VARIANT CALLING", "PROGRAM", "GLNEXUS_BCFTOOLS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
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
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newSnpVCF, sErr := GetVariantType(jointVCF, "SNP", gatkLogLevel, verbose)
						if sErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
							log.Fatalf("FAILED: %v", sErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectSNPs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						snpVCF = newSnpVCF
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //
					if utils.StageHasCompleted(logged, "SelectINDELs", "ALL", seq.ID) {
						msg := fmt.Sprintf("SELECT_INDELS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newIndelVCF, iErr := GetVariantType(jointVCF, "INDEL", gatkLogLevel, verbose)
						if iErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
							log.Fatalf("FAILED: %v", iErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "SelectINDELs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						indelVCF = newIndelVCF
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}

					// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

					if utils.StageHasCompleted(logged, "HardFilteringSNPS", "ALL", seq.ID) {
						msg := fmt.Sprintf("HardFilteringSNPS already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)

					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						newHardFilteredSNPs, hsErr := HardFilterSNPs(snpVCF, gatkLogLevel, verbose)
						if hsErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
							log.Fatalf("FAILED: %v", hsErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringSNPS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						hardFilteredSNPs = newHardFilteredSNPs
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

						newHardFilteredINDELs, hiErr := HardFilterINDELs(indelVCF, gatkLogLevel, verbose)
						if hiErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
							slog.Error("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
							log.Fatalf("FAILED: %v", hiErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "HardFilteringINDELS", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						hardFilteredINDELs = newHardFilteredINDELs
					}

					// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

					if utils.StageHasCompleted(logged, "MergeVcfs", "ALL", seq.ID) {
						msg := fmt.Sprintf("MergeVcfs already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)
						slog.Info(msg)
					} else {
						jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s --VERBOSITY %s`, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF, gatkLogLevel)
						var mErr error
						if verbose {
							mErr = utils.RunBashCmdVerbose(mergeCmdStr)
						} else {
							mErr = utils.RunBashCmd(mergeCmdStr)
						}

						if mErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
							slog.Error("VARIANT CALLING", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
							log.Fatalf("FAILED: %v", mErr)
						}

						jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
					}
					mu.Lock()
					jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)
					mu.Unlock()
				}
			}(seq)
		}
	}
	wg.Wait()

	if !noMerging {
		fmt.Println("Merging ALL Hard filtered VCFs ...")

		mergeCmdStr2 := fmt.Sprintf(`gatk MergeVcfs %s -O %s --VERBOSITY %s`, strings.Join(jointvSlice, " "), finalVcf, gatkLogLevel)

		jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", mergeCmdStr2)
		slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")

		var mErr error

		if verbose {
			mErr = utils.RunBashCmdVerbose(mergeCmdStr2)
		} else {
			mErr = utils.RunBashCmd(mergeCmdStr2)
		}

		if mErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
			slog.Error("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
			log.Fatalf("FAILED: %v", mErr)
		}

		jlog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", "MergeVcfs", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		print("-------------------------------------------------------------------------------------------\n\n")

		return finalVcf, nil
	}

	return finalVcf, nil
}

func VariantCallingConfig(configFile string, species string, maxParallelJobs int, gatkLogLevel string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool) {
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
	finalVcf, fErr := VariantCalling(refFile, bams, outputDir, species, maxParallelJobs, gatkLogLevel, caller, merger, logFilePath, dvVer, modelType, verbose, noMerging)
	if fErr != nil {
		fmt.Printf("Error calling variants: %v\nNo multisample vcf: %s", fErr, finalVcf)
		return
	}
}

// SampleWork bundles per-sample data needed by the worker goroutine.
type SampleWork struct {
	sample  string
	cram    string
	cramDir string // filepath.Dir of the cram, used to derive gvcf output path
}

func validateGvcf(vcf string, verbose bool, quick bool) error {
	if quick {
		tbi := vcf + ".tbi"
		if _, err := os.Stat(tbi); err != nil {
			return fmt.Errorf("TBI index missing for %s", vcf)
		}
		valStr := fmt.Sprintf("bcftools view -h %s > /dev/null", vcf)
		fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
		if verbose {
			return utils.RunBashCmdVerbose(valStr)
		}
		return utils.RunBashCmd(valStr)
	}

	valStr := fmt.Sprintf("bcftools index --tbi --force %s", vcf)
	fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
	if verbose {
		return utils.RunBashCmdVerbose(valStr)
	}
	return utils.RunBashCmd(valStr)
}

func ValidateBam(bam string, ref string, verbose bool, quick bool) error {
	var valStr string
	if quick {
		valStr = fmt.Sprintf("samtools quickcheck %s > /dev/null", bam)
	} else {
		valStr = fmt.Sprintf(`samtools view -T %s -h %s > /dev/null`, ref, bam)
	}
	fmt.Printf("\n-------------------------------------------------------------------\n%s\n------------------------------------------------------------------\n\n", valStr)
	if verbose {
		return utils.RunBashCmdVerbose(valStr)
	}
	return utils.RunBashCmd(valStr)
}
func VariantCallingDir(dataDir string, species string, refVer string, refFasta string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool, gatkLogLevel string, quick bool) {

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

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		fmt.Printf("Error getting chromosomes and IDs: %v\n", err)
		return
	}

	color.Green("All file paths valid\n....................................................\n\n")

	// ================================== Discover samples ===================================== //
	color.Green("Checking Samples in dir structure ...\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "reference_genomes")
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
			samples = append(samples, s)
			fmt.Println(s)
		}
	}
	color.Green("\nFound %d sample(s) in the data directory for %s\n==================================\n\n", len(samples), color.GreenString(species))

	// ======================================= Resolve valid samples up-front ======================================= //
	fmt.Println("Looking for bam files ....")

	var (
		validSamples []SampleWork
		missingBams  []string
		multipleBams []string
	)

	for _, sample := range samples {
		_, cramFiles := FindBamOrVcfs(dataDirAbs, species, sample, refVer, "bams", "")
		switch len(cramFiles) {
		case 0:
			missingBams = append(missingBams, sample)
			color.Red("[%s] bqsr.cram file MISSING! 😒\n", sample)
		case 1:
			cram := cramFiles[0]
			//fmt.Printf("CRAM file found for %s: %s\n", sample, cram)
			color.Green("[%s] bqsr.cram file FOUND! 😊: %s\n", sample, cram)
			err := ValidateBam(cram, refFasta, verbose, quick)
			if err != nil {
				color.Red("[%s] bqsr.cram file is not valid: %v\n", sample, err)
				missingBams = append(missingBams, sample)
			} else {
				validSamples = append(validSamples, SampleWork{sample: sample, cram: cram, cramDir: filepath.Dir(filepath.Dir(filepath.Dir(cram)))})
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

	color.Green("Creating gVCFs for %d valid samples ...................................\n\n", len(validSamples))

	if len(validSamples) == 0 {
		color.Red("No valid samples found. Exiting.")
		return
	}

	// ==================================== Concurrent per-chrom GVCF creation ====================================== //
	var failedSamples []string
	maxWorkers := runtime.NumCPU()
	if maxWorkers < 1 {
		maxWorkers = 1
	}

	type workItem struct {
		sw    SampleWork
		chrom SeqInfo
	}

	workCh := make(chan workItem)
	var wg sync.WaitGroup

	// Fan out workers.
	for i := 0; i < maxWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for item := range workCh {
				sw := item.sw
				chrom := item.chrom

				cramName := filepath.Base(sw.cram)
				gvcfName := strings.Replace(cramName, ".cram", fmt.Sprintf(".%s.g.vcf.gz", chrom.ID), 1)
				gvcfPath := filepath.Join(sw.cramDir, refVer, "gvcfs", gvcfName)

				_, vcfFiles := FindBamOrVcfs(dataDirAbs, species, sw.sample, refVer, "gvcfs", chrom.ID)

				switch len(vcfFiles) {
				case 0:
					color.Red("[%s] gVCF not found for chrom %s — creating ........\n\n", sw.sample, color.GreenString(chrom.ID))
					if _, err := CreateGvcf(sw.cram, refFasta, []SeqInfo{chrom}, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
						color.Red("[%s] Error creating gVCF for chrom %s: %v\n\n", sw.sample, chrom.ID, err)
						failedSamples = append(failedSamples, []string{sw.sample, chrom.ID}...)
					}

				case 1:
					vcf := vcfFiles[0]
					color.Green("[%s] gVCF file for chrom %s exists: %s\n\n", sw.sample, color.BlueString(chrom.ID), vcf)
					fmt.Printf("[%s] checking integrity of gVCF file: %s ..........\n", sw.sample, color.BlueString(vcf))
					if vErr := validateGvcf(vcf, verbose, quick); vErr != nil {
						color.Red("[%s] gVCF %s corrupted. Error: (%v) — re-creating\n", sw.sample, color.BlueString(vcf), vErr)
						if _, err := CreateGvcf(sw.cram, refFasta, []SeqInfo{chrom}, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
							color.Red("[%s] Error re-creating gVCF %s: %v\n", sw.sample, color.BlueString(vcf), err)
							failedSamples = append(failedSamples, []string{sw.sample, chrom.ID}...)
						} else {
							color.Green("[%s] gVCF %s re-created successfully\n", sw.sample, color.BlueString(vcf))
						}
					} else {
						color.Green("[%s] gVCF %s is valid!!\n\n", sw.sample, vcf)
					}

				default:
					color.Red("[%s] Multiple gVCF files found for chrom %s. Please remove extra gVCF in here.\n\n", sw.sample, chrom.ID)
				}
			}
		}()
	}

	// Enqueue (sample × chrom) work items.
	for _, sw := range validSamples {
		for _, chrom := range chroms {
			workCh <- workItem{sw: sw, chrom: chrom}
		}
	}
	close(workCh)
	wg.Wait()

	// ============================== Contig gVCFs (skip if already present) ============================== //
	if len(contigs) == 0 {
		color.Yellow("No contigs found in reference dict. Skipping contig gVCF creation.\n\n")
	} else {
		color.Green("Creating contig gVCFs for %d valid samples ...................................\n\n", len(validSamples))
		for _, sw := range validSamples {
			cramName := filepath.Base(sw.cram)
			contigGvcfName := strings.Replace(cramName, ".cram", ".contigs.g.vcf.gz", 1)
			contigGvcfPath := filepath.Join(sw.cramDir, refVer, "gvcfs", contigGvcfName)

			// Skip if the contig gVCF already exists and is readable.
			if _, statErr := os.Stat(contigGvcfPath); statErr == nil {

				fmt.Printf("[%s] checking integrity of gVCF file: %s .......................\n\n", sw.sample, color.BlueString(contigGvcfPath))
				if vErr := validateGvcf(contigGvcfPath, verbose, quick); vErr != nil {
					color.Red("[%s] gVCF %s corrupted. Error: (%v) — re-creating\n", sw.sample, color.BlueString(contigGvcfPath), vErr)
					if _, err := CreateGvcf(sw.cram, refFasta, contigs, contigGvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
						color.Red("[%s] Error re-creating gVCF %s: %v\n", sw.sample, color.BlueString(contigGvcfPath), err)
					} else {
						color.Green("[%s] gVCF %s re-created successfully\n", sw.sample, color.BlueString(contigGvcfPath))
					}
				} else {
					color.Green("[%s] gVCF %s is valid!!\n\n", sw.sample, color.BlueString(contigGvcfPath))
				}
			} else {
				color.Green("[%s] Creating contig gVCF: %s\n", sw.sample, color.BlueString(contigGvcfPath))
				if _, err := CreateGvcf(sw.cram, refFasta, contigs, contigGvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
					fmt.Printf("[%s] Error creating contig gVCF: %v\n", sw.sample, err)
				}
			}
		}
	}

	fmt.Printf("Missing bams (%d): %v\n==================================\n", len(missingBams), missingBams)

	if !noMerging {
		if len(missingBams) == 0 && len(failedSamples) == 0 {
			fmt.Println("no missing bam files and no failed samples. Calling MergeVcfs")

		}
	}
}
