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

// The dictionary name is the one GATK derives itself, which means the FASTA
// extension is replaced for compressed and uncompressed alike. genome.fa and
// genome.fa.gz share genome.dict; the engine tools look for nothing else.
func TestDictPath(t *testing.T) {
	cases := map[string]string{
		"/g/genome.fa":       "/g/genome.dict",
		"/g/genome.fasta":    "/g/genome.dict",
		"/g/genome.fna":      "/g/genome.dict",
		"/g/genome.fa.gz":    "/g/genome.dict",
		"/g/genome.fasta.gz": "/g/genome.dict",
		"/g/genome.fna.gz":   "/g/genome.dict",
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
	if got, want := DictPath("/g/PN40024.v4.fa.gz"), "/g/PN40024.v4.dict"; got != want {
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
	// actually on disk beside it. PN40024 here carries the pre-correction
	// genome.fa.gz.dict, so this also covers an unmigrated genome: readers go
	// through ResolveDictPath and still find it.
	for _, key := range []string{"VITIS", "MALUS"} {
		for _, r := range genomes[key] {
			if _, err := os.Stat(ResolveDictPath(r.FastaPath)); err != nil {
				t.Errorf("%s/%s: resolved dict %s does not exist", key, r.RefVer, ResolveDictPath(r.FastaPath))
			}
		}
	}
}

// A genome prepared before DictPath was corrected carries genome.fa.gz.dict.
// Readers must still find it, while anything writing a dictionary is pointed
// at the name GATK will look for.
func TestLegacyAndResolveDictPath(t *testing.T) {
	if got := LegacyDictPath("/g/genome.fa.gz"); got != "/g/genome.fa.gz.dict" {
		t.Errorf("LegacyDictPath(compressed) = %q", got)
	}
	if got := LegacyDictPath("/g/genome.fa"); got != "" {
		t.Errorf("LegacyDictPath(uncompressed) = %q, want empty", got)
	}

	dir := t.TempDir()
	ref := filepath.Join(dir, "genome.fa.gz")

	// Nothing on disk yet: callers are steered to the canonical name so a
	// freshly prepared genome lands where GATK looks.
	if got, want := ResolveDictPath(ref), filepath.Join(dir, "genome.dict"); got != want {
		t.Errorf("ResolveDictPath(no dict) = %q, want %q", got, want)
	}

	// Only the legacy name present: use it rather than pretending the
	// canonical one exists.
	legacy := filepath.Join(dir, "genome.fa.gz.dict")
	if err := os.WriteFile(legacy, []byte("@HD\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if got := ResolveDictPath(ref); got != legacy {
		t.Errorf("ResolveDictPath(legacy only) = %q, want %q", got, legacy)
	}

	// Both present: the canonical name wins, so a migrated genome stops
	// depending on the old file.
	canonical := filepath.Join(dir, "genome.dict")
	if err := os.WriteFile(canonical, []byte("@HD\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if got := ResolveDictPath(ref); got != canonical {
		t.Errorf("ResolveDictPath(both) = %q, want %q", got, canonical)
	}
}

// A BLAST database is recognised under the reference's own full name, and in
// both of the shapes makeblastdb writes: one volume leaves a .nsq, several
// leave only a .nal alias beside the per-volume files.
func TestHasBlastDB(t *testing.T) {
	dir := t.TempDir()
	ref := filepath.Join(dir, "genome.fa.gz")

	if HasBlastDB(ref) {
		t.Error("HasBlastDB(nothing on disk) = true")
	}

	nsq := ref + ".nsq"
	if err := os.WriteFile(nsq, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if !HasBlastDB(ref) {
		t.Error("HasBlastDB(single volume) = false")
	}

	// A volumed database: the plain .nsq is absent and the alias stands in for
	// it, so a check that only knew .nsq would rebuild the whole thing.
	if err := os.Remove(nsq); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(ref+".00.nsq", []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if HasBlastDB(ref) {
		t.Error("HasBlastDB(volume file only, no alias) = true")
	}
	if err := os.WriteFile(ref+".nal", []byte("TITLE genome\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if !HasBlastDB(ref) {
		t.Error("HasBlastDB(volumed) = false")
	}
}

// The protein counterpart, in the same two shapes makeblastdb writes: one
// volume leaves a .psq, several leave only a .pal alias.
func TestHasProtBlastDB(t *testing.T) {
	dir := t.TempDir()
	db := filepath.Join(dir, "AllResistanceGenes.fasta")

	if HasProtBlastDB(db) {
		t.Error("HasProtBlastDB(nothing on disk) = true")
	}

	psq := db + ".psq"
	if err := os.WriteFile(psq, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if !HasProtBlastDB(db) {
		t.Error("HasProtBlastDB(single volume) = false")
	}

	// A volumed database: the plain .psq is absent and the alias stands in.
	if err := os.Remove(psq); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(db+".00.psq", []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if HasProtBlastDB(db) {
		t.Error("HasProtBlastDB(volume file only, no alias) = true")
	}
	if err := os.WriteFile(db+".pal", []byte("TITLE prg\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if !HasProtBlastDB(db) {
		t.Error("HasProtBlastDB(volumed) = false")
	}

	// A nucleotide database must not satisfy the protein check, or a caller
	// would hand blastp a db it cannot read and find out much later.
	nuclOnly := filepath.Join(dir, "genome.fa")
	if err := os.WriteFile(nuclOnly+".nsq", []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if HasProtBlastDB(nuclOnly) {
		t.Error("HasProtBlastDB(nucleotide db) = true")
	}
}

// A bare database name is resolved by blastp through $BLASTDB and cannot be
// checked on disk, so it must not be mistaken for a path.
func TestLooksLikeBlastDBPath(t *testing.T) {
	for _, db := range []string{"/home/godwin/prgdb/AllResistanceGenes.fasta", "./prgdb/x", "prgdb/x"} {
		if !LooksLikeBlastDBPath(db) {
			t.Errorf("LooksLikeBlastDBPath(%q) = false", db)
		}
	}
	for _, db := range []string{"swissprot", "nr"} {
		if LooksLikeBlastDBPath(db) {
			t.Errorf("LooksLikeBlastDBPath(%q) = true", db)
		}
	}
}
