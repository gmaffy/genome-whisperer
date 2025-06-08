package variants

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
	"strings"
)

func GetVariantType(vcf string, varType string) error {
	var varVCF string
	if strings.HasSuffix(vcf, ".vcf.gz") {
		suf := fmt.Sprintf(".%s.vcf.gz", varType)
		varVCF = strings.TrimSuffix(vcf, ".vcf.gz") + suf
	} else if strings.HasSuffix(vcf, ".vcf") {
		suf := fmt.Sprintf(".%s.vcf.gz", varType)
		varVCF = strings.TrimSuffix(vcf, ".vcf") + suf
	} else {
		fmt.Println("vcf file must be in vcf or vcf.gz format")
		return fmt.Errorf("vcf file must be in vcf or vcf.gz format")
	}

	cmdStrHap := fmt.Sprintf(`gatk SelectVariants -V %s --select-type-to-include %s -O %s`, vcf, varType, varVCF)
	fmt.Println(cmdStrHap)
	if err := utils.RunBashCmdVerbose(cmdStrHap); err != nil {
		return err
	}
	fmt.Printf("%s created\n", varVCF)
	return nil
}

func HardFilterINDELs(vcf string) error {
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
		return fmt.Errorf("vcf file must be in vcf or vcf.gz format")
	}

	cmdStr := fmt.Sprintf(`gatk VariantFiltration \ 
    -V %s \ 
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
    -O %s`, vcf, vcfCol)

	if err := utils.RunBashCmdVerbose(cmdStr); err != nil {
		return err
	}

	sCmdStr := fmt.Sprintf(`gatk SelectVariants --exclude-filtered -V %s -O %s`, vcfCol, vcfFiltered)
	if err := utils.RunBashCmdVerbose(sCmdStr); err != nil {
		return err
	}
	fmt.Printf("%s created\n", vcfFiltered)
	return nil
}

func HardFilterSNPs(vcf string) error {
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
		return fmt.Errorf("vcf file must be in vcf or vcf.gz format")
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

	if err := utils.RunBashCmdVerbose(cmdStr); err != nil {
		return err
	}

	sCmdStr := fmt.Sprintf(`gatk SelectVariants --exclude-filtered -V %s -O %s`, vcfCol, vcfFiltered)
	if err := utils.RunBashCmdVerbose(sCmdStr); err != nil {
		return err
	}
	fmt.Printf("%s created\n", vcfFiltered)
	return nil
}
