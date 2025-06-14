package pangenome

import (
	"fmt"
	"github.com/gmaffy/genome-whisperer/utils"
)

func GoPan(config string, assembler string) {

	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(config)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Species:", cfg.Species)
	fmt.Println("Read Pairs:", cfg.ReadPairs)

	readPairs := cfg.ReadPairs

	for i, readPair := range readPairs {
		fmt.Println(i, readPair)
		// --------------------------- Align reads to latest reference ---------------------------------------------- //
		fmt.Println("Aligning reads to latest reference")
		fmt.Println("Check log for latest reference genome ...")
		fmt.Printf("Check log to see if alignment for %s is finished. If not, run alignment again. Otherwise skip\n", readPair[3])
		fmt.Println("start aligning reads for ", readPair[3], "")

		fmt.Println("Log success ", readPair[3], "if alignment is finished, return if failed")

		// ---------------------------------- Alignment Statistics -------------------------------------------------- //

		fmt.Printf("Check log to see if alignment stats for %s is finished. If not, run alignment again. Otherwise skip\n", readPair[3])
		fmt.Println("Get alignment stats for ", readPair[3], "")

		fmt.Println("Log success if stats for", readPair[3], "have been generated. return if failed")

		// -------------------------------------- Extract reads ----------------------------------------------------- //
		fmt.Println("Check log to check if reads have been extracted for ", readPair[3], "")
		fmt.Println("Extract reads for ", readPair[3], "from latest ref")
		fmt.Println("Log success if reads have been extracted for", readPair[3], "return if failed")

		// ---------------------------------- Convert unmapped bams to fasta ---------------------------------------- //
		fmt.Println("Check log to check if unmapped bams have been converted to fasta for ", readPair[3], "")
		fmt.Println("Convert unmapped bams to fasta for ", readPair[3], "from latest ref")
		fmt.Println("Log success if unmapped bams have been converted to fasta for", readPair[3], "return if failed")

		//---------------------------------- Assemble unmapped reads ------------------------------------------------ //
		fmt.Println("Check log to check if unmapped reads have been assembled for ", readPair[3], "")
		fmt.Println("Assemble unmapped reads for ", readPair[3], "from latest ref")
		fmt.Println("Log success if unmapped reads have been assembled for", readPair[3], "return if failed")

	}

}
