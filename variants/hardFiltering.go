package variants

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"strings"
)

func GetVariantType(vcf string, varType string) {
	var varVCF string
	if strings.HasSuffix(vcf, ".vcf.gz") {
		suf := fmt.Sprintf(".%s.vcf.gz", varType)
		varVCF = strings.TrimSuffix(vcf, ".vcf.gz") + suf
	} else if strings.HasSuffix(vcf, ".vcf") {
		suf := fmt.Sprintf(".%s.vcf.gz", varType)
		varVCF = strings.TrimSuffix(vcf, ".vcf") + suf
	} else {
		fmt.Println("vcf file must be in vcf or vcf.gz format")
		return
	}

	cmdStrHap := fmt.Sprintf(`gatk SelectVariants -V %s --select-type-to-include %s -O %s`, vcf, varType, varVCF)
	fmt.Println(cmdStrHap)
	utils.RunBashCmdVerbose(cmdStrHap)
	fmt.Printf("%s created\n", varVCF)
}

func HardFilterINDELs(vcf string) {
	var vcfCol string
	var vcfFiltered string
	if strings.HasSuffix(vcf, ".vcf.gz") {
		vcfCol = strings.TrimSuffix(vcf, ".vcf.gz") + ".columns.vcf.gz"
		vcfFiltered = strings.TrimSuffix(vcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
	} else if strings.HasSuffix(vcf, ".vcf") {
		vcfCol = strings.TrimSuffix(vcf, ".vcf") + ".columns.vcf.gz"
		vcfFiltered = strings.TrimSuffix(vcf, ".vcf") + ".hard_filtered.vcf.gz"
	} else {
		fmt.Println("vcf file must be in vcf or vcf.gz format")
		return
	}

	cmdStr := fmt.Sprintf(`gatk VariantFiltration \ 
    -V %s \ 
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
    -O %s`, vcf, vcfCol)

	utils.RunBashCmdVerbose(cmdStr)

	sCmdStr := fmt.Sprintf(`gatk SelectVariants --exclude-filtered -V %s -O %s`, vcfCol, vcfFiltered)
	utils.RunBashCmdVerbose(sCmdStr)
	fmt.Printf("%s created\n", vcfFiltered)
}

func HardFilterSNPs(vcf string) {
	var vcfCol string
	var vcfFiltered string
	if strings.HasSuffix(vcf, ".vcf.gz") {
		vcfCol = strings.TrimSuffix(vcf, ".vcf.gz") + ".columns.vcf.gz"
		vcfFiltered = strings.TrimSuffix(vcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
	} else if strings.HasSuffix(vcf, ".vcf") {
		vcfCol = strings.TrimSuffix(vcf, ".vcf") + ".columns.vcf.gz"
		vcfFiltered = strings.TrimSuffix(vcf, ".vcf") + ".hard_filtered.vcf.gz"
	} else {
		fmt.Println("vcf file must be in vcf or vcf.gz format")
		return
	}

	cmdStr := fmt.Sprintf(`gatk VariantFiltration \
    -V %s \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O %s`, vcf, vcfCol)

	utils.RunBashCmdVerbose(cmdStr)

	sCmdStr := fmt.Sprintf(`gatk SelectVariants --exclude-filtered -V %s -O %s`, vcfCol, vcfFiltered)
	utils.RunBashCmdVerbose(sCmdStr)
}
