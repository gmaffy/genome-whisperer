package alignmentdir

//
//func DirBQSR(ref string, bam string, bqsrbam string, threads int, gatkLogLevel string, verbose bool, knownSites []string, bootstrap bool, quick bool) error {
//	if len(knownSites) == 0 && bootstrap == false {
//		fmt.Println("Either pass a known-sites file or enable bootstrap method")
//		return fmt.Errorf("Either pass a known-sites file or enable bootstrap method")
//	} else if len(knownSites) == 0 && bootstrap == true {
//		fmt.Println("Running with bootstrap method")
//		// =============================== OutPut names =============================================================//
//		var rawVcf string
//		var recalTable string
//		var recalTable2 string
//		var plots string
//		var bqsrBam string
//		var snpVCF string
//		var indelVCF string
//		var hardFilteredSnpVCF string
//		var hardFilteredIndelVCF string
//
//		if strings.HasSuffix(bam, ".bam") {
//			rawVcf = strings.TrimSuffix(bam, ".bam") + ".raw.vcf.gz"
//			recalTable = strings.TrimSuffix(bam, ".bam") + ".recal_table.txt"
//			recalTable2 = strings.TrimSuffix(bam, ".bam") + "recal_table2.txt"
//			plots = strings.TrimSuffix(bam, ".bam") + "recal_table_plots.pdf"
//			bqsrBam = strings.TrimSuffix(bam, ".bam") + "_bqsr.bam"
//
//			snpVCF = strings.TrimSuffix(bam, ".bam") + ".raw.SNP.vcf.gz" //baseName + ".raw.SNP.vcf.gz"
//			indelVCF = strings.TrimSuffix(bam, ".bam") + ".raw.INDEL.vcf.gz"
//			hardFilteredSnpVCF = strings.TrimSuffix(bam, ".bam") + ".raw.SNP.hard_filtered.vcf.gz"
//			hardFilteredIndelVCF = strings.TrimSuffix(bam, ".bam") + ".raw.INDEL.hard_filtered.vcf.gz"
//		} else {
//			recalTable = strings.TrimSuffix(bam, ".cram") + "recal_table.txt"
//			recalTable2 = strings.TrimSuffix(bam, ".cram") + "recal_table2.txt"
//			plots = strings.TrimSuffix(bam, ".cram") + "recal_table_plots.pdf"
//			bqsrBam = strings.TrimSuffix(bam, ".cram") + "_bqsr.cram"
//			rawVcf = strings.TrimSuffix(bam, ".cram") + ".raw.vcf.gz"
//			snpVCF = strings.TrimSuffix(bam, ".cram") + ".raw.SNP.vcf.gz" //baseName + ".raw.SNP.vcf.gz"
//			indelVCF = strings.TrimSuffix(bam, ".cram") + ".raw.INDEL.vcf.gz"
//			hardFilteredSnpVCF = strings.TrimSuffix(bam, ".cram") + ".raw.SNP.hard_filtered.vcf.gz"
//			hardFilteredIndelVCF = strings.TrimSuffix(bam, ".cram") + ".raw.INDEL.hard_filtered.vcf.gz"
//
//		}
//
//		valErr := variants.ValidateGvcf(rawVcf, verbose, quick)
//		if valErr != nil {
//			err := variants.HapCaller(ref, bam, verbose, "raw", rawVcf)
//			if err != nil {
//				return err
//			}
//
//		} else {
//			valErr = variants.ValidateGvcf(snpVCF, verbose, quick)
//			if valErr != nil {
//				snpVCF, err := variants.GetVariantType(rawVcf, "SNP", gatkLogLevel, verbose)
//				if err != nil {
//					return err
//				}
//			}
//		}
//
//		hardFilteredSnpVCF, err = variants.HardFilterSNPs(snpVCF, gatkLogLevel, verbose)
//		if err != nil {
//			return err
//		}
//
//		indelVCF, err = variants.GetVariantType(rawVcf, "INDEL", gatkLogLevel, verbose)
//		if err != nil {
//			return err
//		}
//
//		hardFilteredIndelVCF, err = variants.HardFilterSNPs(indelVCF, gatkLogLevel, verbose)
//		if err != nil {
//			return err
//		}
//
//		knownSites = []string{hardFilteredSnpVCF, hardFilteredIndelVCF}
//
//	}
//}
