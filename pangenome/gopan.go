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

	}

}
