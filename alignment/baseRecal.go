package alignment

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"log"
	"os"
	"strings"
	"sync"
)

func Recalibrate(ref string, bam string, knownSites []string) {

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

	cmdStr := fmt.Sprintf(`gatk BaseRecalibrator -R %s -I %s %s -O %s`, ref, bam, ks, recalTable)
	fmt.Println(cmdStr)
	err := utils.RunBashCmdVerbose(cmdStr)
	if err != nil {
		return
	}

	// ------------------------------------------------ Apply BQSR -------------------------------------------------- //
	aCmdStr := fmt.Sprintf(`gatk ApplyBQSR -R %s -I %s -bqsr %s -O %s`, ref, bam, recalTable, bqsrBam)
	fmt.Println(aCmdStr)
	aErr := utils.RunBashCmdVerbose(aCmdStr)
	if aErr != nil {
		return
	}

	// -------------------------------------- 2nd Recalibration table ----------------------------------------------- //
	cmdStr2 := fmt.Sprintf(`gatk BaseRecalibrator -R %s -I %s %s -O %s`, ref, bqsrBam, ks, recalTable2)
	fmt.Println(cmdStr2)
	bErr := utils.RunBashCmdVerbose(cmdStr2)
	if bErr != nil {
		return
	}

	// ------------------------------------------ Analyse Covariates ------------------------------------------------ //
	cmdStrA := fmt.Sprintf(`gatk AnalyzeCovariates -before %s -after %s -plots %s `, recalTable, recalTable2, plots)
	fmt.Println(cmdStrA)
	A2err := utils.RunBashCmdVerbose(cmdStrA)
	if A2err != nil {
		return
	}

}

func DbSnpBqsr(ref string, bams []string, knownSites []string, numJobs int) {
	fmt.Println("dbSnpBqsr")

	var wg sync.WaitGroup
	sem := make(chan struct{}, numJobs)

	for _, bam := range bams {
		wg.Add(1)
		sem <- struct{}{}
		go func(bam string) {
			defer wg.Done()
			defer func() { <-sem }()
			Recalibrate(ref, bam, knownSites)

		}(bam)
	}

	wg.Wait()
}

func CreateKnownVariants(ref string, bam string) []string {
	rawVCF := strings.TrimSuffix(bam, ".bam") + ".raw.vcf.gz"
	snpVCF := strings.TrimSuffix(bam, ".bam") + ".SNP.vcf.gz"
	indelVCF := strings.TrimSuffix(bam, ".bam") + ".INDEL.vcf.gz"
	hardFilteredSnpVCF := strings.TrimSuffix(bam, ".bam") + ".hard_filtered.SNP.vcf.gz"
	hardFilteredIndelVCF := strings.TrimSuffix(bam, ".bam") + ".hard_filtered.INDEL.vcf.gz"
	cmdStrHap := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -O %s`, ref, bam, rawVCF)
	fmt.Println(cmdStrHap)
	utils.RunBashCmdVerbose(cmdStrHap)

	fmt.Println("Get SNPs and INDELs vcf files from raw VCF file ...")
	variants.GetVariantType(rawVCF, "SNP")
	variants.GetVariantType(rawVCF, "INDEL")

	variants.HardFilterSNPs(snpVCF)
	variants.HardFilterINDELs(indelVCF)

	knownSites := []string{hardFilteredSnpVCF, hardFilteredIndelVCF}
	return knownSites

}

func BootstrapBqsr(ref string, bams []string, numJobs int) {
	fmt.Println("bootstrapBqsr")
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
			knownSites = CreateKnownVariants(ref, bam)

		}(bam)
	}
	wg.Wait()

	// --------------------------------------------- Run BQSR ------------------------------------------------------- //
	DbSnpBqsr(ref, bams, knownSites, numJobs)
}

func BQSRconfig(configPath string, bootstrap bool, jobs int) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configPath)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Species:", cfg.Species)
	fmt.Println("Bams", cfg.Bams)
	fmt.Println("Known-sites", cfg.KnownSites)

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

	if len(knownSites) == 0 && bootstrap == false {
		fmt.Println("Either pass a known-sites file or enable bootstrap method")
		return
	} else if len(knownSites) == 0 && bootstrap == true {
		fmt.Println("Running with bootstrap method")
		BootstrapBqsr(refFile, bams, jobs)
	} else if len(knownSites) > 0 {
		fmt.Println("Running with known-sites flag")
		// ------------------------ Checking Known sites file paths ----------------------------------------- //
		for j, _ := range knownSites {
			_, err := os.Stat(knownSites[j])
			if err != nil {
				fmt.Printf("Known-sites file: %s is not a valid file path", knownSites[j])
				log.Fatal(err)
			}
		}

		// --------------------------- Running dbSnpBQSR ---------------------------------------------------- //
		DbSnpBqsr(refFile, bams, knownSites, jobs)

	} else {
		fmt.Println("Choose either pass a known-sites file or enable bootstrap method, but not both")
		return
	}
}
