package utils

import (
	"os"
	"path/filepath"
	"testing"
)

func TestFastaExt(t *testing.T) {
	cases := map[string]string{
		"/g/genome.fa":        ".fa",
		"/g/genome.fasta":     ".fasta",
		"/g/genome.fna":       ".fna",
		"/g/genome.fa.gz":     ".fa.gz",
		"/g/genome.fasta.gz":  ".fasta.gz",
		"/g/genome.fna.gz":    ".fna.gz",
		"/g/Genome.FA.gz":     ".FA.gz",
		"/g/genome.fa.gz.fai": "",
		"/g/genome.dict":      "",
		"/g/genome.fa.bwt":    "",
		"/g/reads.fq.gz":      "",
		"/g/.fa":              "",
	}
	for in, want := range cases {
		if got := FastaExt(in); got != want {
			t.Errorf("FastaExt(%q) = %q, want %q", in, got, want)
		}
	}
}

func TestIsBgzippedFasta(t *testing.T) {
	for _, in := range []string{"/g/genome.fa.gz", "/g/genome.fasta.gz", "/g/genome.fna.gz", "/g/G.FNA.GZ"} {
		if !IsBgzippedFasta(in) {
			t.Errorf("IsBgzippedFasta(%q) = false, want true", in)
		}
	}
	for _, in := range []string{"/g/genome.fa", "/g/genome.fasta", "/g/genome.dict", "/g/reads.fq.gz"} {
		if IsBgzippedFasta(in) {
			t.Errorf("IsBgzippedFasta(%q) = true, want false", in)
		}
	}
}

// The dictionary is the one sidecar whose name differs between the two forms:
// an uncompressed reference replaces its extension, a compressed one appends.
func TestDictPath(t *testing.T) {
	cases := map[string]string{
		"/g/genome.fa":       "/g/genome.dict",
		"/g/genome.fasta":    "/g/genome.dict",
		"/g/genome.fna":      "/g/genome.dict",
		"/g/genome.fa.gz":    "/g/genome.fa.gz.dict",
		"/g/genome.fasta.gz": "/g/genome.fasta.gz.dict",
		"/g/genome.fna.gz":   "/g/genome.fna.gz.dict",
		"/g/odd.reference":   "/g/odd.dict",
	}
	for in, want := range cases {
		if got := DictPath(in); got != want {
			t.Errorf("DictPath(%q) = %q, want %q", in, got, want)
		}
	}
}

// A dot in the assembly name must not be mistaken for the extension.
func TestDictPathDottedName(t *testing.T) {
	if got, want := DictPath("/g/PN40024.v4.fa.gz"), "/g/PN40024.v4.fa.gz.dict"; got != want {
		t.Errorf("DictPath = %q, want %q", got, want)
	}
	if got, want := DictPath("/g/PN40024.v4.fa"), "/g/PN40024.v4.dict"; got != want {
		t.Errorf("DictPath = %q, want %q", got, want)
	}
}

func TestFaiAndGziPath(t *testing.T) {
	if got, want := FaiPath("/g/genome.fa"), "/g/genome.fa.fai"; got != want {
		t.Errorf("FaiPath = %q, want %q", got, want)
	}
	if got, want := FaiPath("/g/genome.fa.gz"), "/g/genome.fa.gz.fai"; got != want {
		t.Errorf("FaiPath = %q, want %q", got, want)
	}
	if got := GziPath("/g/genome.fa"); got != "" {
		t.Errorf("GziPath(uncompressed) = %q, want \"\"", got)
	}
	if got, want := GziPath("/g/genome.fa.gz"), "/g/genome.fa.gz.gzi"; got != want {
		t.Errorf("GziPath = %q, want %q", got, want)
	}
}

// Discovery on a real directory tree: a bgzipped assembly is found, and where
// both forms of one assembly sit side by side the uncompressed one is used.
func TestGetValidGenomesFromDiskCompressed(t *testing.T) {
	root := t.TempDir()

	write := func(names ...string) {
		for _, name := range names {
			full := filepath.Join(root, name)
			if err := os.MkdirAll(filepath.Dir(full), 0o755); err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(full, nil, 0o644); err != nil {
				t.Fatal(err)
			}
		}
	}

	// bgzipped only, with the sidecars samtools and GATK would write for it
	write(
		"vitis/PN40024/assembly/genome.fa.gz",
		"vitis/PN40024/assembly/genome.fa.gz.dict",
		"vitis/PN40024/assembly/genome.fa.gz.fai",
		"vitis/PN40024/assembly/genome.fa.gz.gzi",
	)
	// both forms present
	write(
		"vitis/V2/assembly/genome.fa",
		"vitis/V2/assembly/genome.fa.gz",
		"vitis/V2/assembly/genome.dict",
	)
	// uncompressed only
	write(
		"malus/GDDH13/assembly/ref.fasta",
		"malus/GDDH13/assembly/ref.dict",
	)

	genomes, err := GetValidGenomesFromDisk(root)
	if err != nil {
		t.Fatal(err)
	}

	find := func(species, refVer string) GenomeRef {
		for _, r := range genomes[species] {
			if r.RefVer == refVer {
				return r
			}
		}
		t.Fatalf("%s/%s not discovered; got %+v", species, refVer, genomes)
		return GenomeRef{}
	}

	if got, want := find("VITIS", "PN40024").FastaPath, filepath.Join(root, "vitis/PN40024/assembly/genome.fa.gz"); got != want {
		t.Errorf("bgzipped assembly: fasta = %q, want %q", got, want)
	}
	if got, want := find("VITIS", "V2").FastaPath, filepath.Join(root, "vitis/V2/assembly/genome.fa"); got != want {
		t.Errorf("both forms present: fasta = %q, want the uncompressed %q", got, want)
	}
	if got, want := find("MALUS", "GDDH13").FastaPath, filepath.Join(root, "malus/GDDH13/assembly/ref.fasta"); got != want {
		t.Errorf("uncompressed assembly: fasta = %q, want %q", got, want)
	}

	// The dictionary each discovered reference resolves to must be the file
	// actually on disk beside it.
	for _, key := range []string{"VITIS", "MALUS"} {
		for _, r := range genomes[key] {
			if _, err := os.Stat(DictPath(r.FastaPath)); err != nil {
				t.Errorf("%s/%s: derived dict %s does not exist", key, r.RefVer, DictPath(r.FastaPath))
			}
		}
	}
}
