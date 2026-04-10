package variants

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
	"golang.org/x/exp/slog"
)

func glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper string) string {
	lines := []string{
		"cd ~/", "rm -r GLnexus.DB",
		fmt.Sprintf("glnexus_cli --config gatk %s/*/*/reference_genomes/%s/gvcfs/*%s.g.vcf.gz > %s/%s.%s.joint.bcf", speciesDir, refVer, sID, outDir, speciesUpper, sID),
		fmt.Sprintf("bcftools view %s/%s.%s.joint.bcf | bgzip -@ 4 -c > %s/%s.%s.joint.vcf.gz", outDir, speciesUpper, sID, outDir, speciesUpper, sID),
	}
	return strings.Join(lines, "\n")
}

func MergeGvcfs(config string, gvcfs []string, dataDir string, species string, refVer string, refFasta string, outDir string, merger string, logFilePath string) {
	fmt.Println("Merging GVCFs ...")

	if logFilePath == "" {
		fmt.Println("No log file path specified. Set log file path")
		return
	}

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
		_, dicfErr := os.Stat(dictFilePath)
		if dicfErr != nil {
			fmt.Printf("Reference dict file: %s does not exist\n", dictFilePath)
			return
		}

		fmt.Printf("Reference fasta: %s\n", refFasta)
		fmt.Printf("Reference dict: %s\n", dictFilePath)

		// --------------------------------- Output directory ------------------------------------------------------- //
		if outDir == "" {
			fmt.Println("No output directory provided ... ")
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

		allPattern := fmt.Sprintf(
			"%s/%s/*/*/reference_genomes/%s/gvcfs/*.g.vcf.gz",
			dataDirAbs, strings.ToLower(species), refVer)

		allMatches, err := filepath.Glob(allPattern)
		if err != nil {
			fmt.Println("Error with glob pattern:", err)
			return
		}

		contigCounts := make(map[string]int)
		for _, match := range allMatches {
			base := filepath.Base(match)
			seqID := strings.TrimSuffix(base, ".g.vcf.gz")
			seqID = seqID[strings.LastIndex(seqID, ".")+1:]
			contigCounts[seqID]++
		}

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

		// ----------------------------------- Create/Open log file ----------------------------------------------------- //
		fmt.Println("Reading log file ...")

		logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
		if err != nil {
			fmt.Printf("failed to open log file %s - %s", logFilePath, err)
			return
		}
		defer logFile.Close()

		jsonHandler := slog.NewJSONHandler(logFile, nil)
		jlog := slog.New(jsonHandler)

		logged := utils.ParseLogFile(logFilePath)

		switch merger {
		case "gatk":
			fmt.Println("Using GATK MergeVcfs")

		case "glnexus":
			fmt.Println("Using GLNEXUS merge")
			for _, sID := range completeContigs {
				if utils.StageHasCompleted(logged, "GLNEXUS", "ALL", sID) {
					msg := fmt.Sprintf("GLNEXUS already completed for ALL gvcfs, CHROMOSOME %s. Skipping.\n\n------------------------------\n\n", sID)
					slog.Info(msg)

				} else {

					speciesDir := filepath.Join(dataDirAbs, strings.ToLower(species))
					speciesUpper := strings.ToUpper(species)
					glnexusCmdStr := glnexusSript(speciesDir, refVer, sID, outDir, speciesUpper)
					fmt.Println(glnexusCmdStr)
					jlog.Info("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", "STARTED", "CMD", glnexusCmdStr)

					slog.Info("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", "STARTED")

					hapErr := utils.RunBashCmdVerbose(glnexusCmdStr)

					if hapErr != nil {
						jlog.Error("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
						slog.Error("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
						log.Fatalf("FAILED: %v", hapErr)
					}

					jlog.Info("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", "COMPLETED")
					slog.Info("GVCF MERGING", "PROGRAM", "GLNEXUS", "SAMPLE", "ALL", "CHROMOSOME", sID, "STATUS", "COMPLETED")
					fmt.Printf("\n\n-------------------------------------------------------------------------------------------\n\n")
				}
			}
		default:
			fmt.Println("Please provide a valid merger: Either gatk or glnexus")

		}
	}
}
