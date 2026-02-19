package variants

import (
	"fmt"
	"os"
)

func MergeGvcfs(config string, gvcfs []string, dataDir string, species string, refName string, outDir string, merger string) {
	fmt.Println("Merging GVCFs ...")
	if config != "" {
		fmt.Printf("Using config file: %s\n", config)
	} else if len(gvcfs) > 0 {
		fmt.Printf("Using gvcfs: %v\n", gvcfs)
	} else {
		fmt.Println("Using directory structure:")
		dInfo, err := os.Stat(dataDir)
		if err != nil {
			fmt.Printf("Error accessing data directory: %s\n", dataDir)
			return
		}
		if !dInfo.IsDir() {
			fmt.Printf("Data directory %s is not a directory\n", dataDir)
			return
		}
		fmt.Printf("Data directory: %s\n", dataDir)

		if species == "" {
			fmt.Println("Please provide species name")
			return
		}
		if refName == "" {
			fmt.Println("Please provide reference name")
		}

		fmt.Printf("Species: %s\n", species)
		fmt.Printf("Reference: %s\n", refName)

		fmt.Println("Scanning directory structure for samples and libraries")

		fmt.Printf("%s\n", dInfo.Name())
	}
	if merger == "gatk" {
		fmt.Println("Using GATK MergeVcfs")
	} else if merger == "glnexus" {
		fmt.Println("Using GLNEXUS merge")
		fmt.Println("")
	} else {
		fmt.Println("Please provide a valid merger: Either gatk or glnexus")

	}

}
