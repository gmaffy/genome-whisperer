package variants

import (
	"bufio"
	"compress/gzip"
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

type ChromosomeSamplePair struct {
	Chromosome string
	Sample     string
}

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

	mw := io.MultiWriter(logFile, os.Stdout)
	log.SetOutput(mw)
	fmt.Println("Log file created.")

	//-------------------------- If resuming (read logfile and check for completed stages) -------------------------- //

	logged := parseLogFile(logFilePath)
	completed := getCompletedStages(logged)
	hapPairCompleted := getHapFinished(logged)

	// --------------------------------------- Reading fasta -------------------------------------------------------- //

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	var wg sync.WaitGroup
	sem := make(chan struct{}, maxParallelJobs)
	var jointvSlice []string

	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)

		// ----------------------- Check if MergeVcfs already completed for this chromosome ------------------------- //

		if _, exists := completed["MergeVcfs"][seq.ID]; exists {
			fmt.Printf("Chromosome %s already processed all steps. Skipping\n", seq.ID)
			continue
		}

		sem <- struct{}{}
		wg.Add(1)

		go func(seq *linear.Seq) {
			defer func() {
				wg.Done()
				<-sem
			}()

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
						log.Printf("Dir:%s\tDir creation\tNA\tFAILED\tNA - Error: %v\n", seq.ID, cErr)
						return
					}
				}
			}

			// ------------------------------------ HAPLOTYPE CALLER (Skip completed) ------------------------------- //
			var vSlice []string
			hapCompletedMap := make(map[string]map[string]bool)
			for _, hap := range hapPairCompleted {
				if hapCompletedMap[hap.Chromosome] == nil {
					hapCompletedMap[hap.Chromosome] = make(map[string]bool)
				}
				hapCompletedMap[hap.Chromosome][hap.Sample] = true
			}

			for _, bam := range bams {
				bamName := filepath.Base(bam)
				fmt.Println(bamName)
				theGVCF := filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))

				if hapCompletedMap[seq.ID] != nil && hapCompletedMap[seq.ID][bamName] {
					fmt.Printf("HaplotypeCaller already completed for BAM FILE %s, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", bamName, seq.ID)
					vSlice = append(vSlice, "-V "+theGVCF)
					continue
				}

				hapCmdStr := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, seq.ID, theGVCF, verbosity)
				vSlice = append(vSlice, "-V "+theGVCF)
				fmt.Println(hapCmdStr)
				log.Printf("%s\tHaplotypeCaller\t%s\tSTARTED\t%s\n", seq.ID, bamName, hapCmdStr)
				hapErr := utils.RunBashCmdVerbose(hapCmdStr)
				fmt.Println(hapErr)
				if hapErr != nil {
					log.Printf("%s\tHaplotypeCaller\t%s\tFAILED\t%s - Error: %v\n", seq.ID, bamName, hapErr) // Use Printf
					return
				}
				log.Printf("%s\tHaplotypeCaller\t%s\tFINISHED\t%s\n", seq.ID, bamName, hapCmdStr)

			}

			jointVCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.vcf.gz")
			snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
			indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
			hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
			theDB := filepath.Join(chromDirPath, chromDir+"DB")

			// ------------------------------------ GENOMICSDBIMPORT (Skip completed) ------------------------------- //

			if _, exists := completed["GenomicsDBImport"][seq.ID]; exists {
				fmt.Printf("GenomicsDBImport already completed for %s. Skipping.\n\n----------------------------\n\n", seq.ID)

			} else {
				vArgs := strings.Join(vSlice, " ")

				//----------------------------- Delete DB if present and delete ------------------------------------- //

				dErr := os.RemoveAll(theDB)
				if dErr != nil {
					fmt.Println("Error removing directory:", dErr)
				} else {
					fmt.Println("Directory removed successfully (if it existed).")
				}

				// -------------------------------------------------------------------------------------------------- //

				gDBImpCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s --genomicsdb-workspace-path %s --tmp-dir %s -L %s --genomicsdb-shared-posixfs-optimizations true --batch-size 50  --bypass-feature-reader`, vArgs, theDB, tmpPath, seq.ID)
				fmt.Println(gDBImpCmdStr)
				log.Printf("%s\tGenomicsDBImport\tALL\tSTARTED\t%s\n", seq.ID, gDBImpCmdStr)

				gErr := utils.RunBashCmdVerbose(gDBImpCmdStr)
				if gErr != nil {
					log.Printf("%s\tGenomicsDBImport\tALL\tFAILED\t%s - Error: %v\n", seq.ID, gDBImpCmdStr, gErr)
					return
				}
				log.Printf("%s\tGenomicsDBImport\tALL\tFINISHED\t%s\n", seq.ID, gDBImpCmdStr)

			}

			// --------------------------------------- GENOTYPE GVCFS (Skip completed) ------------------------------ //

			if _, exists := completed["GenotypeGVCFs"][seq.ID]; exists {
				fmt.Printf("GenomicsDBImport already completed for %s. Skipping.\n\n------------------------------\n\n", seq.ID)

			} else {
				genoCmdStr := fmt.Sprintf(`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s`, refFile, theDB, jointVCF, tmpPath)
				fmt.Println(genoCmdStr)
				log.Printf("%s\tGenotypeGVCFs\tALL\tSTARTED\t%s\n", seq.ID, genoCmdStr)
				gtErr := utils.RunBashCmdVerbose(genoCmdStr)
				if gtErr != nil {
					fmt.Println(gtErr)
					log.Printf("%s\tGenotypeGVCFs\tALL\tFAILED\t%s - Error: %v\n", seq.ID, genoCmdStr, gtErr)
					return
				}
				log.Printf("%s\tGenotypeGVCFs\tALL\tFINISED\t%s\n", seq.ID, genoCmdStr)

			}

			// -------------------------------------- SELECT SNPs (Skip completed) ---------------------------------- //

			fmt.Println("Hard filtered  joint VCF ...")

			if _, exists := completed["SELECT_SNPS"][seq.ID]; exists {
				fmt.Printf("SELECT_SNPS already completed for %s. Skipping.\n\n", seq.ID)

			} else {

				log.Printf("%s\tSELECT_SNPS\tALL\tSTARTED\tSNPS\n", seq.ID)
				sErr := GetVariantType(jointVCF, "SNP")
				if sErr != nil {
					log.Printf("%s\tSELECT_SNPS\tALL\tFAILED\tSNPS - Error: %v\n", seq.ID, sErr)
					return
				}
				log.Printf("%s\tSELECT_SNPS\tALL\tFINISHED\tSNPS\n", seq.ID)

			}

			// -------------------------------------- SELECT INDELS (Skip completed) -------------------------------- //

			if _, exists := completed["SELECT_INDELS"][seq.ID]; exists {
				fmt.Printf("SELECT_INDELS already completed for %s. Skipping.\n", seq.ID)

			} else {
				log.Printf("%s\tSELECT_INDELS\tALL\tSTARTED\tINDELs\n", seq.ID)
				iErr := GetVariantType(jointVCF, "INDEL")
				if iErr != nil {
					log.Printf("%s\tSELECT_INDELS\tALL\tFAILED\tINDELs - Error: %v\n", seq.ID, iErr)
					return
				}
				log.Printf("%s\tSELECT_INDELS\tALL\tFINISHED\tINDELs\n", seq.ID)

			}

			// --------------------------------- HARD FILTERING SNPs (Skip completed) ------------------------------- //

			if _, exists := completed["HardFilteringSNPS"][seq.ID]; exists {
				fmt.Printf("HardFilteringSNPS already completed for %s. Skipping.\n", seq.ID)

			} else {
				log.Printf("%s\tHardFilteringSNPS\tALL\tSTARTED\tSNPs\n", seq.ID)
				hsErr := HardFilterSNPs(snpVCF)
				if hsErr != nil {
					log.Printf("%s\tHardFilteringSNPS\tALL\tFAILED\tSNPs - Error: %v\n", seq.ID, hsErr)
					return
				}
				log.Printf("%s\tHardFilteringSNPS\tALL\tFINISHED\tSNPs\n", seq.ID)

			}

			// -------------------------------- HARD FILTERING INDELS (Skip completed) ------------------------------ //

			if _, exists := completed["HardFilteringINDELS"][seq.ID]; exists {
				fmt.Printf("HardFilteringINDELS already completed for %s. Skipping.\n", seq.ID)

			} else {
				log.Printf("%s\tHardFilteringINDELS\tALL\tSTARTED\tINDELs\n", seq.ID)
				hiErr := HardFilterINDELs(indelVCF)
				if hiErr != nil {
					log.Printf("%s\tHardFilteringINDELS\tALL\tFAILED\tINDELs - Error: %v\n", seq.ID, hiErr)
					return
				}
				log.Printf("%s\tHardFilteringINDELS\tALL\tFINISHED\tINDELs\n", seq.ID)

			}

			// -------------------------------------- MERGE VCFs (Skip completed) ----------------------------------- //

			if _, exists := completed["MergeVcfs"][seq.ID]; exists {
				fmt.Printf("MergeVcfs already completed for %s. Skipping.\n", seq.ID)

			} else {
				mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s`, snpVCF, indelVCF, hardFilteredVCF)
				log.Printf("%s\tMergeVcfs\tALL\tSTART\t%s\n", seq.ID, mergeCmdStr)
				fmt.Println(mergeCmdStr)
				mErr := utils.RunBashCmdVerbose(mergeCmdStr)
				if mErr != nil {
					log.Printf("%s\tMergeVcfs\tALL\tFAILED\t%s - Error: %v\n", seq.ID, mergeCmdStr, mErr)
					return
				}
				log.Printf("%s\tMergeVcfs\tALL\tFINISHED\t%s\n", seq.ID, mergeCmdStr)
				jointvSlice = append(jointvSlice, "-V "+hardFilteredVCF)

			}

		}(seq)

	}
	wg.Wait()
	fmt.Println("All jobs completed.")
	fmt.Println("Merging VCFs ...")
	finalVcf := filepath.Join(out, species+"_"+".joint_hard_filtered.vcf.gz")
	mergeCmdStr := fmt.Sprintf(`gatk MergeVcfs -I %s -O %s`, strings.Join(jointvSlice, " "), finalVcf)
	fmt.Println(mergeCmdStr)
	mErr := utils.RunBashCmdVerbose(mergeCmdStr)
	if mErr != nil {
		fmt.Println(mErr)
		return
	}
	fmt.Println("VCFs merged successfully.")

}

func parseLogFile(logFilePath string) []LogEntry {
	var data []LogEntry
	file, err := os.Open(logFilePath)
	if err != nil {
		fmt.Printf("Log file '%s' not found, starting fresh or assuming no previous runs.\n", logFilePath)
		return data
	}
	defer file.Close()
	fmt.Println("Parsing log file ...")

	scanner := bufio.NewScanner(file)
	if !scanner.Scan() {
		return data
	}
	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")
		if len(line) != 5 {
			fmt.Printf("Skipping malformed log line: '%s', Expected 5 columns, got %d\n", scanner.Text(), len(line))
			continue
		}

		timeStampChrom := strings.Split(line[0], " ")
		timeStamp := strings.Join(timeStampChrom[:len(timeStampChrom)-1], " ")

		chromID := timeStampChrom[len(timeStampChrom)-1]
		stage := line[1]
		bam := line[2]
		status := line[3]
		details := line[4]

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
	return data
}

func getHapFinished(logEntries []LogEntry) []ChromosomeSamplePair {
	var hapFinished []ChromosomeSamplePair
	for _, entry := range logEntries {
		if entry.Program == "HaplotypeCaller" && entry.Status == "FINISHED" {
			pair := ChromosomeSamplePair{Chromosome: entry.Chromosome, Sample: entry.Bam}
			hapFinished = append(hapFinished, pair)
		}
	}
	return hapFinished
}

func getCompletedStages(logEntries []LogEntry) map[string]map[string]bool {
	cs := make(map[string]map[string]bool)
	for _, entry := range logEntries {

		if entry.Program == "GenomicsDBImport" && entry.Status == "FINISHED" {
			if _, ok := cs["GenomicsDBImport"]; !ok {
				cs["GenomicsDBImport"] = make(map[string]bool)
			}
			cs["GenomicsDBImport"][entry.Chromosome] = true
		} else if entry.Program == "GenotypeGVCFs" && entry.Status == "FINISHED" {
			if _, ok := cs["GenotypeGVCFs"]; !ok {
				cs["GenotypeGVCFs"] = make(map[string]bool)
			}
			cs["GenotypeGVCFs"][entry.Chromosome] = true
		} else if entry.Program == "SELECT_SNPS" && entry.Status == "FINISHED" {
			if _, ok := cs["SELECT_SNPS"]; !ok {
				cs["SELECT_SNPS"] = make(map[string]bool)
			}
			cs["SELECT_SNPS"][entry.Chromosome] = true

		} else if entry.Program == "SELECT_INDELS" && entry.Status == "FINISHED" {
			if _, ok := cs["SELECT_INDELS"]; !ok {
				cs["SELECT_INDELS"] = make(map[string]bool)
			}
			cs["SELECT_INDELS"][entry.Chromosome] = true
		} else if entry.Program == "HardFilteringSNPS" && entry.Status == "FINISHED" {
			if _, ok := cs["HardFilteringSNPS"]; !ok {
				cs["HardFilteringSNPS"] = make(map[string]bool)
			}
			cs["HardFilteringSNPS"][entry.Chromosome] = true
		} else if entry.Program == "HardFilteringINDELS" && entry.Status == "FINISHED" {
			if _, ok := cs["HardFilteringINDELS"]; !ok {
				cs["HardFilteringINDELS"] = make(map[string]bool)
			}
			cs["HardFilteringINDELS"][entry.Chromosome] = true
		} else if entry.Program == "MergeVcfs" && entry.Status == "FINISHED" {
			if _, ok := cs["MergeVcfs"]; !ok {
				cs["MergeVcfs"] = make(map[string]bool)
			}
			cs["MergeVcfs"][entry.Chromosome] = true
		}
	}
	return cs
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
