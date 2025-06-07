package variants

import (
	"compress/gzip"
	"fmt"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
)

func VariantCalling(refFile string, bams []string, out string) {

	fmt.Println("Working on FASTA file ...")
	fna, err := os.Open(refFile)
	if err != nil {
		log.Fatalf("Failed to open FASTA file: %v", err)
	}
	defer func(fna *os.File) {
		err := fna.Close()
		if err != nil {
			panic(err)
		}
	}(fna)

	var reader io.Reader = fna
	if strings.HasSuffix(refFile, ".gz") {
		gzReader, err := gzip.NewReader(fna)
		if err != nil {
			log.Fatalf("Failed to create gzip reader: %v", err)
		}
		defer gzReader.Close()
		reader = gzReader
	}

	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	// Iterate over sequences in the FASTA file
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		fmt.Println(seq.ID)
		chromDir := filepath.Join(out, strings.ReplaceAll(seq.ID, ".", "_"))
		gvcfPath := filepath.Join(chromDir, "gvcfs")
		tmpPath := filepath.Join(chromDir, "tmp")
		tmp2Path := filepath.Join(chromDir, "tmp2")
		vcfPath := filepath.Join(chromDir, "VCFs")

		cErr := os.MkdirAll(chromDir, 0755)
		if cErr != nil {
			log.Fatalf("Error creating directory: %v", cErr)
			return
		}

		gErr := os.MkdirAll(gvcfPath, 0755)
		if gErr != nil {
			log.Fatalf("Error creating directory: %v", gErr)
			return
		}

		tErr := os.MkdirAll(tmpPath, 0755)
		if tErr != nil {
			log.Fatalf("Error creating directory: %v", cErr)
			return
		}

		t2Err := os.MkdirAll(tmp2Path, 0755)
		if t2Err != nil {
			log.Fatalf("Error creating directory: %v", t2Err)
			return
		}

		vErr := os.MkdirAll(vcfPath, 0755)
		if vErr != nil {
			log.Fatalf("Error creating directory: %v", vErr)
			return
		}

		for _, bam := range bams {
			
		}

	}

}
