package variants

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

func MergeGvcfs(config string, gvcfs []string, dataDir string, species string, refVer string, refFasta string, outDir string, merger string) {
	fmt.Println("Merging GVCFs ...")
	if config != "" {
		fmt.Printf("Using config file: %s\n", config)
	} else if len(gvcfs) > 0 {
		fmt.Printf("Using gvcfs: %v\n", gvcfs)
	} else {
		fmt.Println("Using directory structure:")

		// ============================================= Check paths ================================================ //

		// --------------------------------------------- Data dir --------------------------------------------------- //
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
		fmt.Printf("Data directory: %s\n", dataDirAbs)

		// ---------------------------------------- Species & Version ----------------------------------------------- //

		if species == "" {
			fmt.Println("Please provide species name")
			return
		}

		if refVer == "" {
			fmt.Println("Please provide reference version name")
			return
		}

		fmt.Printf("Species: %s, ver: %s\n", species, refVer)

		// -------------------------------------- Reference fasta & Dict -------------------------------------------- //
		if refFasta == "" {
			fmt.Println("Please provide reference name")
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
		_, dicfErr := os.Stat(dictFilePath)
		if dicfErr != nil {
			fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
			return
		}

		fmt.Printf("Reference fasta: %s\n", refFasta)
		fmt.Printf("Reference dict: %s\n", dictFilePath)

		// --------------------------------- Output directory ------------------------------------------------------- //
		if outDir == "" {
			fmt.Println("No output directory provided. ")
			outDir = filepath.Join(dataDirAbs, species, "VCFs", refVer)
			fmt.Printf("Creating output directory at: %s...\n", outDir)
			mkdirErr := os.MkdirAll(outDir, 0755)
			if mkdirErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", outDir, mkdirErr)
				return

			}
			//return
		} else {
			fmt.Printf("Output directory: %s\n", outDir)
			oInfo, err := os.Stat(outDir)
			if err != nil {
				fmt.Printf("Error accessing output directory: %s\n", outDir)
				return
			}
			if !oInfo.IsDir() {
				fmt.Printf("Output directory %s is not a directory\n", outDir)
				fmt.Println("Creating output directory ...")
				mkdirErr := os.MkdirAll(outDir, 0755)
				if mkdirErr != nil {
					fmt.Printf("Failed to create output directory %s: %v\n", outDir, mkdirErr)
					return

				}
			}

		}

		// ================================== Checking Samples in dir Structure ===================================== //
		fmt.Println("Checking Samples in dir structure ...")
		pattern := filepath.Join(dataDir, species, "*", "*", "reference_genomes")
		matches, err := filepath.Glob(pattern)
		if err != nil {
			panic(err)
		}

		samples := []string{}
		for _, match := range matches {
			// match is the reference_genomes dir, parent is sampleDir
			sample := filepath.Base(filepath.Dir(match))
			//fmt.Println(sample)
			samples = append(samples, sample)
		}
		fmt.Printf("there are %v Samples in the data directory for %s\n==================================\n", len(samples), species)

		// ================================= Getting gvcfs for each sample ========================================= //

		fmt.Println("Looking for gvcfs ....")

		fmt.Printf("Reference dict file: %s\n", dictFilePath)
		dictFile, err := os.Open(dictFilePath)
		if err != nil {
			fmt.Printf("Error opening reference dict file: %s: %v\n", dictFilePath, err)
			return
		}
		defer dictFile.Close()

		scanner := bufio.NewScanner(dictFile)

		missingContigs := []string{}
		contigsMissingSamples := []string{}
		completeContigs := []string{}
		//
		//for scanner.Scan() {
		//	//fmt.Printf("%s\n", scanner.Text())
		//	if strings.HasPrefix(scanner.Text(), "@SQ") {
		//		seqID := strings.Split(scanner.Text(), "\t")[1][3:]
		//		pattern = fmt.Sprintf("%s/%s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz", dataDirAbs, strings.ToLower(species), refVer, seqID)
		//
		//		matches, err = filepath.Glob(pattern)
		//		if err != nil {
		//			fmt.Println("Error with glob pattern:", err)
		//			return
		//		}
		//
		//		if len(matches) == 0 {
		//			missingContigs = append(missingContigs, seqID)
		//		} else if len(matches) != len(samples) {
		//			contigsMissingSamples = append(contigsMissingSamples, seqID)
		//		} else {
		//			completeContigs = append(completeContigs, seqID)
		//		}
		//
		//		//seqLen := strings.Split(scanner.Text(), "\t")[2][3:]
		//		//fmt.Printf("Found seqID: %s, Len:%v\n", seqID, seqLen)
		//		//gvcfs = append(gvcfs, filepath.Join(dataDir, seqID, species+"_"+refVer+"_"+seqID+".g.vcf.gz"))
		//		//glnexusCmdStr := fmt.Sprintf(`glnexus_cli --config gatk  %s/%s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz > %s/%s.%s.joint.bcf`, dataDirAbs, strings.ToLower(species), refVer, seqID, outDir, strings.ToUpper(species), seqID)
		//		//bcftoolsCmdStr := fmt.Sprintf(`bcftools view %s/%s.%s.joint.bcf | bgzip -@ 4 -c > %s/%s.%s.joint.vcf.gz`, outDir, strings.ToUpper(species), seqID, outDir, strings.ToUpper(species), seqID)
		//		//
		//		//fmt.Println(glnexusCmdStr)
		//		//fmt.Println(bcftoolsCmdStr)
		//	}
		//}
		//
		//fmt.Printf("Found %v contigs with no gvcfs: %v\n", len(missingContigs), missingContigs)
		//fmt.Printf("Found %v contigs with some missing samples: %v\n", len(contigsMissingSamples), contigsMissingSamples)
		//fmt.Printf("Found %v contigs with all samples: %v\n", len(completeContigs), completeContigs)

		// One glob to get all gvcf files upfront
		allPattern := fmt.Sprintf("%s/%s/*/*/reference_genomes/%s/gvcfs/*.g.vcf.gz", dataDirAbs, strings.ToLower(species), refVer)
		allMatches, err := filepath.Glob(allPattern)
		if err != nil {
			fmt.Println("Error with glob pattern:", err)
			return
		}

		// Build map: seqID -> count of matching files
		contigCounts := make(map[string]int)
		for _, match := range allMatches {
			base := filepath.Base(match) // e.g. "SAMPLE.CONTIG.g.vcf.gz"
			// extract seqID from filename - adjust based on your naming convention
			// if filename is *SEQID.g.vcf.gz, extract the part before .g.vcf.gz
			seqID := strings.TrimSuffix(base, ".g.vcf.gz")
			seqID = seqID[strings.LastIndex(seqID, ".")+1:] // get last segment
			contigCounts[seqID]++
		}

		// Now scan is fast - just map lookups
		for scanner.Scan() {
			if strings.HasPrefix(scanner.Text(), "@SQ") {
				seqID := strings.Split(scanner.Text(), "\t")[1][3:]
				count := contigCounts[seqID]

				if count == 0 {
					missingContigs = append(missingContigs, seqID)
				} else if count != len(samples) {
					contigsMissingSamples = append(contigsMissingSamples, seqID)
				} else {
					completeContigs = append(completeContigs, seqID)
				}
			}
		}

		fmt.Printf("Found %v contigs with no gvcfs:\n", len(missingContigs))
		fmt.Printf("Found %v contigs with some missing samples: %v\n", len(contigsMissingSamples), contigsMissingSamples)
		fmt.Printf("Found %v contigs with all samples: %v\n", len(completeContigs), completeContigs)

		switch merger {
		case "gatk":
			fmt.Println("Using GATK MergeVcfs")
			break
		case "glnexus":
			fmt.Println("Using GLNEXUS merge")
			for _, sID := range completeContigs {
				speciesDir := filepath.Join(dataDirAbs, strings.ToLower(species))
				glnexusCmdStr := fmt.Sprintf(`cd ~/

glnexus_cli --config gatk  %s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz > %s/%s.%s.joint.bcf

bcftools view %s/%s.%s.joint.bcf | bgzip -@ 4 -c > %v/%v.%v.joint.vcf.gz`, speciesDir, refVer, sID, outDir, strings.ToUpper(species), sID, sID, outDir, strings.ToUpper(species), sID, outDir, strings.ToUpper(species), sID)
				fmt.Println(glnexusCmdStr)
			}
		default:
			fmt.Println("Please provide a valid merger: Either gatk or glnexus")

		}

	}
	//if merger == "gatk" {
	//	fmt.Println("Using GATK MergeVcfs")
	//} else if merger == "glnexus" {
	//	fmt.Println("Using GLNEXUS merge")
	//	for sID := range completeContigs {
	//
	//	}
	//} else {
	//	fmt.Println("Please provide a valid merger: Either gatk or glnexus")
	//
	//}

}
