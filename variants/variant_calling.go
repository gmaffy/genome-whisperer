package variants

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/gmaffy/genome-whisperer/utils"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
	"sync"
)

type LogEntry struct {
	Timestamp  string
	Chromosome string
	Program    string
	Bam        string
	Status     string
	Cmd        string
}

func VariantCalling(refFile string, bams []string, out string, species string) {

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
	logFilePath := filepath.Join(out, "variant_calling.log")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	log.SetOutput(logFile)

	log.Printf("CHROM\tPROGRAM\tBAM\tSTATUS\tCMD\n")

	// --------------------------------------- Reading fasta -------------------------------------------------------- //

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	var wg sync.WaitGroup
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		wg.Add(1)
		go func(seq *linear.Seq) {
			defer wg.Done()
			fmt.Println(seq.ID)

			chromDir := strings.ReplaceAll(seq.ID, ".", "_")
			chromDirPath := filepath.Join(out, chromDir)
			gvcfPath := filepath.Join(chromDirPath, "gvcfs")
			tmpPath := filepath.Join(chromDirPath, "tmp")
			tmp2Path := filepath.Join(chromDirPath, "tmp2")
			vcfPath := filepath.Join(chromDirPath, "VCFs")

			dirsToCreate := []string{chromDirPath, gvcfPath, tmpPath, tmp2Path, vcfPath}
			for _, dir := range dirsToCreate {
				if _, err := os.Stat(dir); os.IsNotExist(err) {
					cErr := os.MkdirAll(dir, 0755)
					if cErr != nil {
						log.Fatalf("Dir:%s\tDir creation\tNA\tFAILED\tNA")
						return
					}
				}
			}

			var vSlice []string
			for _, bam := range bams {
				bamName := filepath.Base(bam)
				theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
				hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF`, refFile, bam, seq.ID, theGVCF)
				vSlice = append(vSlice, "-V "+theGVCF)
				fmt.Println(hapCmdStr)
				log.Printf("%s\tHaplotypeCaller\t%s\tSTARTED\t%s\n", seq.ID, bamName, hapCmdStr)
				hapErr := utils.RunBashCmdVerbose(hapCmdStr)
				fmt.Println(hapErr)
				if hapErr != nil {
					log.Fatalf("%s\tHaplotypeCaller\t%s\tFAILED\t%s\n", seq.ID, bamName, hapErr)
					return
				}
				log.Printf("%s\tHaplotypeCaller\t%s\tFINISHED\t%s\n", seq.ID, bamName, hapCmdStr)

			}

			//bName := strings.Replace(filepath.Base(bams[0]), ".bam", "", 1)

			vArgs := strings.Join(vSlice, " ")
			theDB := filepath.Join(chromDirPath, chromDir+"DB")
			jointVCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.vcf.gz")

			gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader`, vArgs, theDB, tmpPath, seq.ID)
			fmt.Println(gDBImpCmdStr)
			log.Printf("%s\tGenomicsDBImport\tALL\tSTARTED\t%s\n", seq.ID, gDBImpCmdStr)

			gErr := utils.RunBashCmdVerbose(gDBImpCmdStr)
			if gErr != nil {
				log.Fatalf("%s\tGenomicsDBImport\tALL\tFAILED\t%s\n", seq.ID, gDBImpCmdStr)
				return
			}
			log.Printf("%s\tGenomicsDBImport\tALL\tFINISHED\t%s\n", seq.ID, gDBImpCmdStr)

			genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s`, refFile, theDB, jointVCF, tmpPath)
			fmt.Println(genoCmdStr)
			log.Printf("%s\tGenotypeGVCFs\tALL\tSTARTED\t%s\n", seq.ID, genoCmdStr)
			gtErr := utils.RunBashCmdVerbose(genoCmdStr)
			if gtErr != nil {
				fmt.Println(gtErr)
				log.Fatalf("%s\tGenotypeGVCFs\tALL\tFAILED\t%s\n", seq.ID, genoCmdStr)
				return
			}
			log.Printf("%s\tGenotypeGVCFs\tALL\tFINISED\t%s\n", seq.ID, genoCmdStr)

			snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
			indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
			hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"

			fmt.Println("Hard filtered joint VCF ...")
			log.Printf("%s\tSELECT_SNPS\tALL\tSTARTED\tSNPS\n", seq.ID)
			sErr := GetVariantType(jointVCF, "SNP")
			if sErr != nil {
				log.Fatalf("%s\tSELECT_SNPS\tALL\tFAILED\tSNPS\n", seq.ID)
				return
			}
			log.Printf("%s\tSelectVariants\tALL\tFINISHED\tSNPS\n", seq.ID)

			log.Printf("%s\tSelectVariants\tALL\tSTARTED\tINDELs\n", seq.ID)
			iErr := GetVariantType(jointVCF, "INDEL")
			if iErr != nil {
				log.Fatalf("%s\tSelectVariants\tALL\tFAILED\tINDELs\n", seq.ID)
				return
			}
			log.Printf("%s\tSelectVariants\tALL\tFINISHED\tINDELs\n", seq.ID)

			log.Printf("%s\tHardFiltering\tALL\tSTARTED\tSNPs\n", seq.ID)
			hsErr := HardFilterSNPs(snpVCF)
			if hsErr != nil {
				log.Fatalf("%s\tHardFiltering\tALL\tFAILED\tSNPs\n", seq.ID)
				return
			}
			log.Printf("%s\tHardFiltering\tALL\tFINISHED\tSNPs\n", seq.ID)

			log.Printf("%s\tHardFiltering\tALL\tSTARTED\tINDELs\n", seq.ID)
			hiErr := HardFilterINDELs(indelVCF)
			if hiErr != nil {
				log.Fatalf("%s\tHardFiltering\tALL\tFAILED\tINDELs\n", seq.ID)
				return
			}
			log.Printf("%s\tHardFiltering\tALL\tFINISHED\tINDELs\n", seq.ID)

			mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
			log.Printf("%s\tMergeVcfs\tALL\tSTART\t%s\n", seq.ID, mergeCmdStr)
			fmt.Println(mergeCmdStr)
			mErr := utils.RunBashCmdVerbose(mergeCmdStr)
			if mErr != nil {
				log.Fatalf("%s\tMergeVcfs\tALL\tFAILED\t%s\n", seq.ID, mergeCmdStr)
				return
			}
			log.Printf("%s\tMergeVcfs\tALL\tFINISHED\t%s\n", seq.ID, mergeCmdStr)

		}(seq)

	}
	wg.Wait()

}

func parseLogFile(logFilePath string) []LogEntry {

	file, err := os.Open(logFilePath)
	if err != nil {
		log.Fatalf("Failed to open log file for parsing: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.TrimLeadingSpace = true

	records, err := reader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}
	if len(records) < 2 {
		log.Fatalf("Expected at least 2 rows, got %d", len(records))
	}

	var data []LogEntry
	for _, row := range records[1:] {
		if len(row) != 6 {
			log.Fatalf("Expected 6 columns, got %d", len(row))
		}
		timeStamp := row[0]
		chromID := row[1]
		stage := row[2]
		bam := row[3]
		status := row[4]
		details := row[5]
		r := LogEntry{
			Timestamp:  timeStamp,
			Chromosome: chromID,
			Program:    stage,
			Bam:        bam,
			Status:     status,
			Cmd:        details,
		}
		data = append(data, r)
	}

	//scanner := bufio.NewScanner(file)
	//for scanner.Scan() {
	//	line := scanner.Text()
	//	fields := strings.Split(line, "\t")
	//	fmt.Println(fields)
	//}
	return data
}

// Define log messages for different stages
//const (
//	logStageHaplotypeCaller  = "HaplotypeCaller completed for %s"
//	logStageGenomicsDBImport = "GenomicsDBImport completed for %s"
//	logStageGenotypeGVCFs    = "GenotypeGVCFs completed for %s"
//	logStageSNPsFiltered     = "SNP hard filtering completed for %s"
//	logStageINDELsFiltered   = "INDEL hard filtering completed for %s"
//	logStageMergeVcfs        = "MergeVcfs completed for %s"
//	logStageChromosomeDone   = "ALL STAGES COMPLETED FOR CHROMOSOME %s"
//)
//
//// completedStages tracks which stages are done for each chromosome
//type completedStages map[string]map[string]bool
//
//func VariantCalling(refFile string, bams []string, out string, species string) {
//
//	// --------------------------------------- Opening fasta file --------------------------------------------------- //
//	fmt.Println("Working on FASTA file ...")
//	fna, err := os.Open(refFile)
//	if err != nil {
//		log.Fatalf("Failed to open FASTA file: %v", err)
//	}
//	defer func(fna *os.File) {
//		err := fna.Close()
//		if err != nil {
//			panic(err)
//		}
//	}(fna)
//
//	var reader io.Reader = fna
//	if strings.HasSuffix(refFile, ".gz") {
//		gzReader, err := gzip.NewReader(fna)
//		if err != nil {
//			log.Fatalf("Failed to create gzip reader: %v", err)
//		}
//		defer gzReader.Close()
//		reader = gzReader
//	}
//
//	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
//	logFilePath := filepath.Join(out, "variant_calling.log")
//	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
//	if err != nil {
//		log.Fatalf("Failed to open log file: %v", err)
//	}
//	defer logFile.Close()
//
//	// Direct standard log output to the file
//	log.SetOutput(logFile)
//
//	// Read existing log to determine completed chromosomes/stages
//	completed := parseLogFile(logFilePath)
//
//	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
//	sc := seqio.NewScanner(r)
//
//	var wg sync.WaitGroup
//	for sc.Next() {
//		seq := sc.Seq().(*linear.Seq)
//		chromID := seq.ID // Use original seq.ID for logging consistency
//
//		// Check if this chromosome is already fully processed
//		if completed[chromID] != nil && completed[chromID][fmt.Sprintf(logStageChromosomeDone, chromID)] {
//			log.Printf("Chromosome %s already completed. Skipping.\n", chromID)
//			continue
//		}
//
//		wg.Add(1)
//		go func(seq *linear.Seq) {
//			defer wg.Done()
//
//			fmt.Printf("Processing chromosome: %s\n", seq.ID)
//			log.Printf("Starting processing for chromosome: %s\n", seq.ID)
//
//			chromDir := strings.ReplaceAll(seq.ID, ".", "_")
//			chromDirPath := filepath.Join(out, chromDir)
//			gvcfPath := filepath.Join(chromDirPath, "gvcfs")
//			tmpPath := filepath.Join(chromDirPath, "tmp")
//			tmp2Path := filepath.Join(chromDirPath, "tmp2")
//			vcfPath := filepath.Join(chromDirPath, "VCFs")
//
//			// Create directories if they don't exist
//			dirsToCreate := []string{chromDirPath, gvcfPath, tmpPath, tmp2Path, vcfPath}
//			for _, dir := range dirsToCreate {
//				if _, err := os.Stat(dir); os.IsNotExist(err) {
//					cErr := os.MkdirAll(dir, 0755)
//					if cErr != nil {
//						log.Fatalf("Error creating directory %s: %v", dir, cErr)
//						return
//					}
//				}
//			}
//
//			// --- HaplotypeCaller Stage ---
//			haplotypeCallerKey := fmt.Sprintf(logStageHaplotypeCaller, chromID)
//			if completed[chromID] == nil || !completed[chromID][haplotypeCallerKey] {
//				var vSlice []string
//				for _, bam := range bams {
//					bamName := filepath.Base(bam)
//					theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, ".bam", fmt.Sprintf("_%s.g.vcf.gz", chromDir), 1))
//					hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF`, refFile, bam, seq.ID, theGVCF)
//					vSlice = append(vSlice, "-V "+theGVCF)
//					fmt.Printf("Running HaplotypeCaller for %s: %s\n", bamName, hapCmdStr)
//					log.Printf("Running HaplotypeCaller for %s: %s\n", bamName, hapCmdStr)
//					if err := utils.RunBashCmdVerbose(hapCmdStr); err != nil {
//						log.Fatalf("HaplotypeCaller failed for %s: %v", bamName, err)
//					}
//				}
//				log.Printf(haplotypeCallerKey + "\n") // Log completion of HaplotypeCaller for this chromosome
//			} else {
//				log.Printf("HaplotypeCaller already completed for %s. Skipping.\n", chromID)
//				// Reconstruct vSlice if skipping HaplotypeCaller, necessary for subsequent steps
//
//			}
//
//			// --- GenomicsDBImport Stage ---
//			var vSlice []string
//			for _, bam := range bams {
//				bamName := filepath.Base(bam)
//				theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, ".bam", fmt.Sprintf("_%s.g.vcf.gz", chromDir), 1))
//				vSlice = append(vSlice, "-V "+theGVCF)
//			}
//
//			genomicsDBImportKey := fmt.Sprintf(logStageGenomicsDBImport, chromID)
//			if completed[chromID] == nil || !completed[chromID][genomicsDBImportKey] {
//				vArgs := strings.Join(vSlice, " ")
//				theDB := filepath.Join(chromDirPath, chromDir+"DB")
//				gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader`, vArgs, theDB, tmpPath, seq.ID)
//				fmt.Printf("Running GenomicsDBImport for %s: %s\n", chromID, gDBImpCmdStr)
//				log.Printf("Running GenomicsDBImport for %s: %s\n", chromID, gDBImpCmdStr)
//				if err := utils.RunBashCmdVerbose(gDBImpCmdStr); err != nil {
//					log.Fatalf("GenomicsDBImport failed for %s: %v", chromID, err)
//				}
//				log.Printf(genomicsDBImportKey + "\n") // Log completion
//			} else {
//				log.Printf("GenomicsDBImport already completed for %s. Skipping.\n", chromID)
//			}
//
//			// --- GenotypeGVCFs Stage ---
//			genotypeGVCFsKey := fmt.Sprintf(logStageGenotypeGVCFs, chromID)
//			if completed[chromID] == nil || !completed[chromID][genotypeGVCFsKey] {
//				bName := strings.Replace(filepath.Base(bams[0]), ".bam", "", 1)
//				theDB := filepath.Join(chromDirPath, chromDir+"DB")
//				jointVCF := filepath.Join(vcfPath, bName+"_"+chromDir+".joint.vcf.gz")
//				genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s`, refFile, theDB, jointVCF, tmpPath)
//				fmt.Printf("Running GenotypeGVCFs for %s: %s\n", chromID, genoCmdStr)
//				log.Printf("Running GenotypeGVCFs for %s: %s\n", chromID, genoCmdStr)
//				if err := utils.RunBashCmdVerbose(genoCmdStr); err != nil {
//					log.Fatalf("GenotypeGVCFs failed for %s: %v", chromID, err)
//				}
//				log.Printf(genotypeGVCFsKey + "\n") // Log completion
//			} else {
//				log.Printf("GenotypeGVCFs already completed for %s. Skipping.\n", chromID)
//			}
//
//			bName := strings.Replace(filepath.Base(bams[0]), ".bam", "", 1)
//			jointVCF := filepath.Join(vcfPath, bName+"_"+chromDir+".joint.vcf.gz")
//			snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
//			indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
//			hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
//
//			// --- HardFilterSNPs Stage ---
//			snpsFilteredKey := fmt.Sprintf(logStageSNPsFiltered, chromID)
//			if completed[chromID] == nil || !completed[chromID][snpsFilteredKey] {
//				fmt.Printf("Hard filtering SNPs for %s...\n", chromID)
//				log.Printf("Hard filtering SNPs for %s...\n", chromID)
//				// Assuming GetVariantType and HardFilterSNPs write to snpVCF
//				if err := GetVariantType(jointVCF, "SNP"); err != nil { // Assuming GetVariantType returns an error
//					log.Fatalf("GetVariantType for SNP failed for %s: %v", chromID, err)
//				}
//				if err := HardFilterSNPs(snpVCF); err != nil { // Assuming HardFilterSNPs returns an error
//					log.Fatalf("HardFilterSNPs failed for %s: %v", chromID, err)
//				}
//				log.Printf(snpsFilteredKey + "\n") // Log completion
//			} else {
//				log.Printf("SNP hard filtering already completed for %s. Skipping.\n", chromID)
//			}
//
//			// --- HardFilterINDELs Stage ---
//			indelsFilteredKey := fmt.Sprintf(logStageINDELsFiltered, chromID)
//			if completed[chromID] == nil || !completed[chromID][indelsFilteredKey] {
//				fmt.Printf("Hard filtering INDELs for %s...\n", chromID)
//				log.Printf("Hard filtering INDELs for %s...\n", chromID)
//				// Assuming GetVariantType and HardFilterINDELs write to indelVCF
//				if err := GetVariantType(jointVCF, "INDEL"); err != nil { // Assuming GetVariantType returns an error
//					log.Fatalf("GetVariantType for INDEL failed for %s: %v", chromID, err)
//				}
//				if err := HardFilterINDELs(indelVCF); err != nil { // Assuming HardFilterINDELs returns an error
//					log.Fatalf("HardFilterINDELs failed for %s: %v", chromID, err)
//				}
//				log.Printf(indelsFilteredKey + "\n") // Log completion
//			} else {
//				log.Printf("INDEL hard filtering already completed for %s. Skipping.\n", chromID)
//			}
//
//			// --- MergeVcfs Stage ---
//			mergeVcfsKey := fmt.Sprintf(logStageMergeVcfs, chromID)
//			if completed[chromID] == nil || !completed[chromID][mergeVcfsKey] {
//				mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
//				fmt.Printf("Running MergeVcfs for %s: %s\n", chromID, mergeCmdStr)
//				log.Printf("Running MergeVcfs for %s: %s\n", chromID, mergeCmdStr)
//				if err := utils.RunBashCmdVerbose(mergeCmdStr); err != nil {
//					log.Fatalf("MergeVcfs failed for %s: %v", chromID, err)
//				}
//				log.Printf(mergeVcfsKey + "\n") // Log completion
//			} else {
//				log.Printf("MergeVcfs already completed for %s. Skipping.\n", chromID)
//			}
//
//			log.Printf(fmt.Sprintf(logStageChromosomeDone, chromID) + "\n") // Mark chromosome as fully done
//			fmt.Printf("Finished processing chromosome: %s\n", seq.ID)
//
//		}(seq)
//	}
//	wg.Wait()
//	log.Println("All variant calling processes completed.")
//	fmt.Println("All variant calling processes completed.")
//}
//
//// parseLogFile reads the log file and returns a map of completed stages for each chromosome.
//func parseLogFile(logFilePath string) completedStages {
//	completed := make(completedStages)
//
//	file, err := os.Open(logFilePath)
//	if err != nil {
//		if os.IsNotExist(err) {
//			return completed // Log file doesn't exist yet, nothing completed
//		}
//		log.Fatalf("Failed to open log file for parsing: %v", err)
//	}
//	defer file.Close()
//
//	scanner := bufio.NewScanner(file)
//	for scanner.Scan() {
//		line := scanner.Text()
//		// Parse based on the log messages defined
//		// Example: "HaplotypeCaller completed for chr1"
//		// Example: "ALL STAGES COMPLETED FOR CHROMOSOME chr1"
//
//		// Extract chromosome ID from the log message
//		var chromID string
//		if strings.Contains(line, "for chromosome ") {
//			parts := strings.Split(line, "for chromosome ")
//			if len(parts) > 1 {
//				chromID = strings.TrimSpace(strings.TrimSuffix(parts[1], "."))
//			}
//		} else if strings.Contains(line, "for ") {
//			parts := strings.Split(line, "for ")
//			if len(parts) > 1 {
//				chromID = strings.TrimSpace(strings.Split(parts[1], ":")[0]) // For messages like "Running HaplotypeCaller for chr1: ..."
//				if strings.Contains(chromID, " ") {                          // Clean up if it captures more than just the ID
//					chromID = strings.Split(chromID, " ")[0]
//				}
//			}
//		}
//
//		if chromID == "" {
//			continue // Couldn't parse chromosome ID
//		}
//
//		if completed[chromID] == nil {
//			completed[chromID] = make(map[string]bool)
//		}
//
//		// Check for stage completion markers
//		if strings.Contains(line, fmt.Sprintf(logStageHaplotypeCaller, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageHaplotypeCaller, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageGenomicsDBImport, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageGenomicsDBImport, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageGenotypeGVCFs, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageGenotypeGVCFs, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageSNPsFiltered, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageSNPsFiltered, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageINDELsFiltered, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageINDELsFiltered, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageMergeVcfs, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageMergeVcfs, chromID)] = true
//		} else if strings.Contains(line, fmt.Sprintf(logStageChromosomeDone, chromID)) {
//			completed[chromID][fmt.Sprintf(logStageChromosomeDone, chromID)] = true
//		}
//	}
//
//	if err := scanner.Err(); err != nil {
//		log.Fatalf("Error reading log file: %v", err)
//	}
//
//	return completed
//}
