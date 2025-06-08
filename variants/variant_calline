package variants

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func VariantCalling(refFile string, bams []string, out string, species string) {
	// --------------------------------------- Opening fasta file --------------------------------------------------- //
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
	logFilePath := filepath.Join(out, "variant_calling.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	// Initialize slog with JSON handler
	logger := slog.New(slog.NewJSONHandler(logFile, &slog.HandlerOptions{
		AddSource: true,
		Level:     slog.LevelInfo,
	}))

	// Read existing log entries to check completed processes
	completedProcesses := readCompletedProcesses(logFilePath)

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	var wg sync.WaitGroup
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		wg.Add(1)
		go func(seq *linear.Seq) {
			defer wg.Done()
			logger.Info("Processing sequence", "seqID", seq.ID)

			chromDir := strings.ReplaceAll(seq.ID, ".", "_")
			chromDirPath := filepath.Join(out, chromDir)
			gvcfPath := filepath.Join(chromDirPath, "gvcfs")
			tmpPath := filepath.Join(chromDirPath, "tmp")
			tmp2Path := filepath.Join(chromDirPath, "tmp2")
			vcfPath := filepath.Join(chromDirPath, "VCFs")

			if err := os.MkdirAll(chromDirPath, 0755); err != nil {
				logger.Error("Error creating directory", "path", chromDirPath, "error", err)
				log.Fatalf("Error creating directory: %v", err)
			}
			if err := os.MkdirAll(gvcfPath, 0755); err != nil {
				logger.Error("Error creating directory", "path", gvcfPath, "error", err)
				log.Fatalf("Error creating directory: %v", err)
			}
			if err := os.MkdirAll(tmpPath, 0755); err != nil {
				logger.Error("Error creating directory", "path", tmpPath, "error", err)
				log.Fatalf("Error creating directory: %v", err)
			}
			if err := os.MkdirAll(tmp2Path, 0755); err != nil {
				logger.Error("Error creating directory", "path", tmp2Path, "error", err)
				log.Fatalf("Error creating directory: %v", err)
			}
			if err := os.MkdirAll(vcfPath, 0755); err != nil {
				logger.Error("Error creating directory", "path", vcfPath, "error", err)
				log.Fatalf("Error creating directory: %v", err)
			}

			var vSlice []string
			for _, bam := range bams {
				bamName := filepath.Base(bam)
				theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
				processKey := fmt.Sprintf("HaplotypeCaller:%s:%s", seq.ID, bamName)

				if _, completed := completedProcesses[processKey]; !completed {
					logger.Info("Starting HaplotypeCaller", "seqID", seq.ID, "bam", bamName)
					hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF`, refFile, bam, seq.ID, theGVCF)
					fmt.Println(hapCmdStr)
					if err := utils.RunBashCmdVerbose(hapCmdStr); err != nil {
						log.Fatalf("HaplotypeCaller failed for %s: %v", seq.ID, err)
					}
					logger.Info("Completed HaplotypeCaller", "seqID", seq.ID, "bam", bamName, "processKey", processKey)
				} else {
					logger.Info("Skipping HaplotypeCaller (already completed)", "seqID", seq.ID, "bam", bamName)
				}
				vSlice = append(vSlice, "-V "+theGVCF)
			}

			bName := strings.Replace(filepath.Base(bams[0]), ".bam", "", 1)
			vArgs := strings.Join(vSlice, " ")
			theDB := filepath.Join(chromDirPath, chromDir+"DB")
			jointVCF := filepath.Join(vcfPath, bName+"_"+chromDir+".joint.vcf.gz")
			gdbProcessKey := fmt.Sprintf("GenomicsDBImport:%s", seq.ID)

			if _, completed := completedProcesses[gdbProcessKey]; !completed {
				logger.Info("Starting GenomicsDBImport", "seqID", seq.ID)
				gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50 --bypass-feature-reader`, vArgs, theDB, tmpPath, seq.ID)
				fmt.Println(gDBImpCmdStr)
				if err := utils.RunBashCmdVerbose(gDBImpCmdStr); err != nil {
					log.Fatalf("GenomicsDBImport failed for %s: %v", seq.ID, err)
				}
				logger.Info("Completed GenomicsDBImport", "seqID", seq.ID, "processKey", gdbProcessKey)
			} else {
				logger.Info("Skipping GenomicsDBImport (already completed)", "seqID", seq.ID)
			}

			genoProcessKey := fmt.Sprintf("GenotypeGVCFs:%s", seq.ID)
			if _, completed := completedProcesses[genoProcessKey]; !completed {
				logger.Info("Starting GenotypeGVCFs", "seqID", seq.ID)
				genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s`, refFile, theDB, jointVCF, tmpPath)
				fmt.Println(genoCmdStr)
				if err := utils.RunBashCmdVerbose(genoCmdStr); err != nil {
					log.Fatalf("GenotypeGVCFs failed for %s: %v", seq.ID, err)
				}
				logger.Info("Completed GenotypeGVCFs", "seqID", seq.ID, "processKey", genoProcessKey)
			} else {
				logger.Info("Skipping GenotypeGVCFs (already completed)", "seqID", seq.ID)
			}

			snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
			indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
			hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
			hardFilterProcessKey := fmt.Sprintf("HardFilter:%s", seq.ID)

			if _, completed := completedProcesses[hardFilterProcessKey]; !completed {
				logger.Info("Starting Hard Filtering", "seqID", seq.ID)
				fmt.Println("Hard filtered joint VCF ...")
				if err := GetVariantType(jointVCF, "INDEL"); err != nil {
					log.Fatalf("GetVariantType for INDEL failed for %s: %v", seq.ID, err)
				}

				if err := GetVariantType(jointVCF, "SNP"); err != nil { // Assuming GetVariantType returns an error
					log.Fatalf("GetVariantType for SNP failed for %s: %v", seq.ID, err)
				}

				if err := HardFilterSNPs(snpVCF); err != nil {
					log.Fatalf("Hard Filtering for SNPS failed for %s: %v", seq.ID, err)
				}
				if err := HardFilterINDELs(indelVCF); err != nil {
					log.Fatalf("Hard Filtering for INDELS failed for %s: %v", seq.ID, err)
				}
				mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
				fmt.Println(mergeCmdStr)
				if err := utils.RunBashCmdVerbose(mergeCmdStr); err != nil {
					log.Fatalf("MergeVcfs failed for %s: %v", seq.ID, err)
				}
				logger.Info("Completed Hard Filtering", "seqID", seq.ID, "processKey", hardFilterProcessKey)
			} else {
				logger.Info("Skipping Hard Filtering (already completed)", "seqID", seq.ID)
			}
		}(seq)
	}
	wg.Wait()
}

// readCompletedProcesses reads the log file and returns a map of completed process keys
func readCompletedProcesses(logFilePath string) map[string]struct{} {
	completed := make(map[string]struct{})
	file, err := os.Open(logFilePath)
	if err != nil {
		if os.IsNotExist(err) {
			return completed
		}
		log.Fatalf("Failed to read log file: %v", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		var logEntry struct {
			Level      string `json:"level"`
			Msg        string `json:"msg"`
			ProcessKey string `json:"processKey"`
		}
		if err := json.Unmarshal([]byte(line), &logEntry); err == nil {
			if logEntry.Level == "INFO" && strings.HasPrefix(logEntry.Msg, "Completed") && logEntry.ProcessKey != "" {
				completed[logEntry.ProcessKey] = struct{}{}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		log.Fatalf("Error scanning log file: %v", err)
	}

	return completed
}
