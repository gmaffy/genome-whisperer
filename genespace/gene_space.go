package gene_space

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

type GeneDetails struct {
	Gene  string
	Start int
	Stop  int
}

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}

	cleanup := func() { f.Close() }

	// Check suffix
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}

		cleanup = func() {
			gz.Close()
			f.Close()
		}

		return gz, cleanup, nil
	}

	// Plain text VCF
	return f, cleanup, nil
}

func geneStartStop(gff, chrom string, start, stop int) (map[string][2]string, error) {
	dic := make(map[string][2]string)

	f, err := os.Open(gff)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) < 9 {
			continue
		}

		pos, err := strconv.Atoi(fields[3])
		if err != nil {
			continue
		}

		if fields[0] == chrom && fields[2] == "mRNA" && pos >= start && pos < stop {
			gene := strings.SplitN(strings.SplitN(fields[8], ";", 2)[0], "=", 2)[1]
			dic[gene] = [2]string{fields[3], fields[4]}
		}
	}

	return dic, scanner.Err()
}

func GeneSpace(gffPath, vcfTable, chrom string, start, stop int, resLines, susLines []string, descFile, prgFile string) {
	f, cleanup, err := openVCF(vcfTable)
	if err != nil {
		panic(err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(f, false)
	if err != nil {
		panic(err)
	}

	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		if v.Chromosome == chrom && int(v.Pos) >= start && int(v.Pos) < stop {
			continue
		}

	}

}
