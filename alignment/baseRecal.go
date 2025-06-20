package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"log/slog"
	"os"
	"strings"
	"sync"
)

func Recalibrate(ref string, bam string, knownSites []string, logFilePath string) error {

	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	logged := utils.ParseLogFile(logFilePath)

	var kslice []string
	for _, site := range knownSites {
		kslice = append(kslice, "--known-sites "+site)
	}
	ks := strings.Join(kslice, " ")

	recalTable := strings.TrimSuffix(bam, ".bam") + "recal_table.txt"
	recalTable2 := strings.TrimSuffix(bam, ".bam") + "recal_table2.txt"
	plots := strings.TrimSuffix(bam, ".bam") + "recal_table_plots.pdf"
	bqsrBam := strings.TrimSuffix(bam, ".bam") + "_bqsr.bam"

	// ----------------------------------- First Recalibration table ------------------------------------------------ //
	if utils.StageHasCompleted(logged, "BaseRecalibrator", bam, "ALL") {
		msg := fmt.Sprintf("BaseRecalibrator already completed for bam file: %s. Skipping ....\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		cmdStr := fmt.Sprintf(`gatk BaseRecalibrator -R %s -I %s %s -O %s`, ref, bam, ks, recalTable)
		slog.Info(fmt.Sprintf("%s\n-------------------------------------------------\n\n", cmdStr))

		err = utils.RunBashCmdVerbose(cmdStr)
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED- %v", err))
			slog.Error("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED- %v", err))
			return err
		}
		jlog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	// ------------------------------------------------ Apply BQSR -------------------------------------------------- //
	if utils.StageHasCompleted(logged, "ApplyBQSR", bam, "ALL") {
		msg := fmt.Sprintf("ApplyBQSR already completed for bam file: %s. Skipping ....\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "ApplyBQSR", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		aCmdStr := fmt.Sprintf(`gatk ApplyBQSR -R %s -I %s -bqsr %s -O %s`, ref, bam, recalTable, bqsrBam)
		slog.Info(fmt.Sprintf("%s\n-------------------------------------------------\n\n", aCmdStr))
		aErr := utils.RunBashCmdVerbose(aCmdStr)
		if aErr != nil {
			jlog.Error("BQSR", "PROGRAM", "ApplyBQSR", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED- %v", aErr))
			slog.Error("BQSR", "PROGRAM", "ApplyBQSR", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED- %v", aErr))
			return err
		}
		jlog.Info("BQSR", "PROGRAM", "ApplyBQSR", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "ApplyBQSR", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	// -------------------------------------- 2nd Recalibration table ----------------------------------------------- //
	if utils.StageHasCompleted(logged, "BaseRecalibrator", bqsrBam, "ALL") {
		msg := fmt.Sprintf("BaseRecalibrator already completed for bam file: %s. Skipping ....\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		cmdStr2 := fmt.Sprintf(`gatk BaseRecalibrator -R %s -I %s %s -O %s`, ref, bqsrBam, ks, recalTable2)
		slog.Info(fmt.Sprintf("%s\n-------------------------------------------------\n\n", cmdStr2))

		bErr := utils.RunBashCmdVerbose(cmdStr2)
		if bErr != nil {
			jlog.Error("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED- %v", err))
			slog.Error("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bqsrBam, "STATUS", fmt.Sprintf("FAILED- %v", err))
			return err
		}
		jlog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "BaseRecalibrator", "SAMPLE", bqsrBam, "STATUS", "COMPLETED")
	}

	// ------------------------------------------ Analyse Covariates ------------------------------------------------ //
	if utils.StageHasCompleted(logged, "AnalyzeCovariates", bqsrBam, "ALL") {
		msg := fmt.Sprintf("AnalyzeCovariates already completed for bam file: %s. Skipping ....\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "AnalyzeCovariates", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		cmdStrA := fmt.Sprintf(`gatk AnalyzeCovariates -before %s -after %s -plots %s `, recalTable, recalTable2, plots)
		slog.Info(fmt.Sprintf("%s\n-------------------------------------------------\n\n", cmdStrA))

		A2err := utils.RunBashCmdVerbose(cmdStrA)
		if A2err != nil {
			jlog.Error("BQSR", "PROGRAM", "AnalyzeCovariates", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED- %v", err))
			slog.Error("BQSR", "PROGRAM", "AnalyzeCovariates", "SAMPLE", bqsrBam, "STATUS", fmt.Sprintf("FAILED- %v", err))
			return err
		}
		jlog.Info("BQSR", "PROGRAM", "AnalyzeCovariates", "SAMPLE", bqsrBam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "AnalyzeCovariates", "SAMPLE", bqsrBam, "STATUS", "COMPLETED")

	}

	return nil

}

func DbSnpBqsr(ref string, bams []string, knownSites []string, numJobs int, logFilePath string) error {
	fmt.Println("dbSnpBqsr")
	// ------------------------------------------- Open log file ---------------------------------------------------- //
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	logged := utils.ParseLogFile(logFilePath)

	var wg sync.WaitGroup
	sem := make(chan struct{}, numJobs)

	for _, bam := range bams {
		wg.Add(1)
		sem <- struct{}{}
		go func(bam string) {
			defer wg.Done()
			defer func() { <-sem }()
			jlog.Info("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			slog.Info("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
			if utils.StageHasCompleted(logged, "RECALIBRATE", bam, "ALL") {
				msg := fmt.Sprintf("RECALIBRATE already completed for bam file: %s. Skipping ....\n", bam)
				slog.Info(msg)
			} else {
				err := Recalibrate(ref, bam, knownSites, logFilePath)
				if err != nil {
					jlog.Error("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
					slog.Error("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
					return
				}

			}

			jlog.Info("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("BQSR", "PROGRAM", "RECALIBRATE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

		}(bam)
	}

	wg.Wait()
	return nil
}

func CreateKnownVariants(ref string, bam string, logFilePath string) ([]string, error) {

	// ------------------------------------------- Open log file ---------------------------------------------------- //
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)
	jlog := slog.New(jsonHandler)

	jlog.Info("BQSR", "PROGRAM", "INITIALISE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("BQSR", "PROGRAM", "INITIALISE", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")

	logged := utils.ParseLogFile(logFilePath)

	// --------------------------------------------- Output Paths --------------------------------------------------- //
	rawVCF := strings.TrimSuffix(bam, ".bam") + ".raw.vcf.gz"
	snpVCF := strings.TrimSuffix(bam, ".bam") + ".SNP.vcf.gz"
	indelVCF := strings.TrimSuffix(bam, ".bam") + ".INDEL.vcf.gz"
	hardFilteredSnpVCF := strings.TrimSuffix(bam, ".bam") + ".hard_filtered.SNP.vcf.gz"
	hardFilteredIndelVCF := strings.TrimSuffix(bam, ".bam") + ".hard_filtered.INDEL.vcf.gz"

	// --------------------------------------------- Haplotype Caller ----------------------------------------------- //
	if utils.StageHasCompleted(logged, "HaplotypeCaller", bam, "ALL") {
		msg := fmt.Sprintf("HaplotypeCaller already completed for bam file: %s. Skipping\n", bam)
		slog.Info(msg)

	} else {
		cmdStrHap := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -O %s`, ref, bam, rawVCF)
		fmt.Println(cmdStrHap)
		jlog.Info("BQSR", "PROGRAM", "HaplotypeCaller", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info(fmt.Sprintf("%s\n-------------------------------------------------\n\n", cmdStrHap))

		err = utils.RunBashCmdVerbose(cmdStrHap)
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "HaplotypeCaller", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BQSR", "PROGRAM", "HaplotypeCaller", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
			return []string{}, err
		}
		jlog.Info("BQSR", "PROGRAM", "HaplotypeCaller", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "HaplotypeCaller", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	// --------------------------------------------- Select SNPs ---------------------------------------------------- //
	fmt.Println("Get SNPs and INDELs vcf files from raw VCF file ...")
	if utils.StageHasCompleted(logged, "SelectSNPs", bam, "ALL") {
		msg := fmt.Sprintf("SelectSNPs already completed for bam file: %s. Skipping\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "STATUS", "STARTED")
		err = variants.GetVariantType(rawVCF, "SNP")
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
			return []string{}, err
		}
		jlog.Info("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "SelectSNPs", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	// --------------------------------------------- Select INDELs -------------------------------------------------- //
	if utils.StageHasCompleted(logged, "SelectINDELs", bam, "ALL") {
		msg := fmt.Sprintf("SelectINDELs already completed for bam file: %s. Skipping\n", bam)
		slog.Info(msg)

	} else {
		jlog.Info("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "STATUS", "STARTED")

		err = variants.GetVariantType(rawVCF, "INDEL")
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
			return []string{}, err
		}
		jlog.Info("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "SelectINDELs", "SAMPLE", bam, "STATUS", "COMPLETED")

	}

	// ---------------------------------------- HARD FILTER SNPs ---------------------------------------------------- //
	if utils.StageHasCompleted(logged, "HardFilterSNPs", bam, "ALL") {
		msg := fmt.Sprintf("HardFilterSNPs already completed for bam file: %s. Skipping\n", bam)
		slog.Info(msg)
	} else {
		jlog.Error("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Error("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "STATUS", "STARTED")

		err = variants.HardFilterSNPs(snpVCF)
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
			return []string{}, err
		}
		jlog.Info("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Error("BQSR", "PROGRAM", "HardFilterSNPs", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	// -------------------------------------- HARD FILTER INDELS ---------------------------------------------------- //
	if utils.StageHasCompleted(logged, "HardFilterINDELs", bam, "ALL") {
		msg := fmt.Sprintf("HardFilterINDELs already completed for bam file: %s. Skipping\n", bam)
		slog.Info(msg)
	} else {
		jlog.Info("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Error("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "STATUS", "STARTED")
		err = variants.HardFilterINDELs(indelVCF)
		if err != nil {
			jlog.Error("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Error("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
			return []string{}, err
		}
		jlog.Info("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Error("BQSR", "PROGRAM", "HardFilterINDELs", "SAMPLE", bam, "STATUS", "COMPLETED")
	}

	_, sErr := os.Stat(hardFilteredSnpVCF)
	_, iErr := os.Stat(hardFilteredIndelVCF)
	if sErr != nil || iErr != nil {
		msg := fmt.Sprintf("Hard Filtered VCF files not found for bam file: %s. Skipping\n", bam)
		slog.Info(msg)
		return []string{}, fmt.Errorf(msg)
	}

	knownSites := []string{hardFilteredSnpVCF, hardFilteredIndelVCF}
	return knownSites, nil

}

func BootstrapBqsr(ref string, bams []string, numJobs int, logFilePath string) error {
	fmt.Println("bootstrapBqsr")
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)

	jlog := slog.New(jsonHandler)

	jlog.Info("BQSR", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
	slog.Info("BQSR", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED") //, "CMD", "ALL")
	var wg sync.WaitGroup
	sem := make(chan struct{}, numJobs)

	var knownSites []string
	// -------------------------------- Get Known variants from each bam file --------------------------------------- //
	for _, bam := range bams {
		wg.Add(1)
		sem <- struct{}{}
		go func(bam string) {
			defer wg.Done()
			defer func() { <-sem }()
			jlog.Info("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "STATUS", "COMPLETED")
			knownSites, err = CreateKnownVariants(ref, bam, logFilePath)
			if err != nil {
				jlog.Error("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
				slog.Error("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "STATUS", fmt.Sprintf("FAILED - %v", err))
				return
			}
			jlog.Info("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("BQSR", "PROGRAM", "CREATE_KNOWN_VARIANTS", "SAMPLE", bam, "STATUS", "COMPLETED")

		}(bam)
	}
	wg.Wait()

	// --------------------------------------------- Run BQSR ------------------------------------------------------- //
	jlog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", "STARTED")
	err = DbSnpBqsr(ref, bams, knownSites, numJobs, logFilePath)
	if err != nil {
		jlog.Error("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Error("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", "STARTED")
		return err
	}
	jlog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", "COMPLETED")
	return nil
}

func BQSRconfig(configPath string, bootstrap bool, jobs int, logFilePath string) {

	// ----------------------------------------- Log file ----------------------------------------------------------- //
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		log.Fatalf("Failed to open log file: %v", err)
	}
	defer logFile.Close()

	jsonHandler := slog.NewJSONHandler(logFile, nil)

	jlog := slog.New(jsonHandler)

	// ----------------------------------------- Read Config File --------------------------------------------------- //
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}

	// ------------------------------------------ Checking  Bam Paths ----------------------------------------------- //

	refFile := cfg.Reference
	bams := cfg.Bams
	knownSites := cfg.KnownSites

	fmt.Printf("bams: %v\n", bams)
	if len(bams) == 0 {
		fmt.Println("You must provide at least one bam file")
		return
	} else {
		for i, _ := range bams {
			_, err := os.Stat(bams[i])
			if err != nil {
				fmt.Printf("Bam file: %s is not a valid file path", bams[i])
				log.Fatal(err)
			}
		}
	}

	// ------------------------------------------- Create log file -------------------------------------------------- //
	// ----------------------------------- Create/Open log file ----------------------------------------------------- //
	fmt.Println("Reading log file ...")

	if len(knownSites) == 0 && bootstrap == false {
		fmt.Println("Either pass a known-sites file or enable bootstrap method")
		return
	} else if len(knownSites) == 0 && bootstrap == true {
		fmt.Println("Running with bootstrap method")
		BootstrapBqsr(refFile, bams, jobs, logFilePath)
	} else if len(knownSites) > 0 {
		fmt.Println("Running with known-sites flag")
		// -------------------------------- Checking Known sites file paths ----------------------------------------- //
		for j, _ := range knownSites {
			_, err := os.Stat(knownSites[j])
			if err != nil {
				fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
				log.Fatal(err)
			}
		}

		// ----------------------------------- Running dbSnpBQSR ---------------------------------------------------- //
		jlog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		slog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", "STARTED")
		err := DbSnpBqsr(refFile, bams, knownSites, jobs, logFilePath)
		if err != nil {
			jlog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			slog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", fmt.Sprintf("FAILED - %v", err))
			return
		}
		jlog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("BQSR", "PROGRAM", "BQSR", "SAMPLE", "ALL", "STATUS", "COMPLETED")

	} else {
		fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
		return
	}
}
