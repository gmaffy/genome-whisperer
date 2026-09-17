package utils

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// fastaExts are the reference FASTA extensions recognised on disk. The
// bgzipped forms come first so that trimming ".fa.gz" never stops short at
// ".gz" and leaves a stray ".fa" behind.
var fastaExts = []string{".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"}

// FastaExt returns the recognised FASTA extension at the end of name, exactly
// as it is spelled on disk, or "" when name is not a FASTA. The comparison
// itself is case-insensitive, so Genome.FA.gz is recognised.
func FastaExt(name string) string {
	// The length guard is on the base name, so a dotfile called ".fa" is not
	// mistaken for an assembly named "" with a .fa extension.
	base := strings.ToLower(filepath.Base(name))
	for _, ext := range fastaExts {
		if len(base) > len(ext) && strings.HasSuffix(base, ext) {
			return name[len(name)-len(ext):]
		}
	}
	return ""
}

// IsFasta reports whether name ends in a recognised FASTA extension,
// compressed or not.
func IsFasta(name string) bool { return FastaExt(name) != "" }

// IsBgzippedFasta reports whether name is one of the block-compressed forms.
func IsBgzippedFasta(name string) bool {
	return strings.HasSuffix(strings.ToLower(FastaExt(name)), ".gz")
}

// DictPath returns the sequence dictionary belonging to a reference FASTA.
//
// The name is the one GATK derives itself: the recognised FASTA extension is
// replaced by ".dict", for compressed and uncompressed alike, so genome.fa and
// genome.fa.gz both map to genome.dict. GATK's engine tools — HaplotypeCaller,
// BaseRecalibrator, ApplyBQSR, GenotypeGVCFs — accept no other spelling. They
// derive the name strictly and never fall back to genome.fa.gz.dict, even
// though the Picard-derived tools in the same jar do accept it. A reference
// carrying only the longer name fails at engine start-up with "Fasta dict file
// ... does not exist".
//
// Both forms of one assembly therefore share a single dictionary. That is
// correct rather than a collision: a bgzipped reference holds the same
// sequences under the same names as the flat file it was compressed from, so
// the dictionary describing one describes the other.
//
// A name that is not a recognised FASTA falls back to replacing its last
// extension, which is what every call site did before this existed.
func DictPath(refFasta string) string {
	if ext := FastaExt(refFasta); ext != "" {
		return strings.TrimSuffix(refFasta, ext) + ".dict"
	}
	return strings.TrimSuffix(refFasta, filepath.Ext(refFasta)) + ".dict"
}

// LegacyDictPath returns the dictionary name earlier versions wrote beside a
// bgzipped reference: the whole file name plus ".dict", the shape samtools
// uses for the .fai and .gzi. Nothing writes it any more, but genomes prepared
// before DictPath was corrected still carry it, and re-deriving a dictionary
// for a multi-gigabyte assembly is slow enough to be worth avoiding. Empty for
// an uncompressed reference, which never had a second spelling.
func LegacyDictPath(refFasta string) string {
	if !IsBgzippedFasta(refFasta) {
		return ""
	}
	return refFasta + ".dict"
}

// ResolveDictPath returns the dictionary that actually exists beside a
// reference: the canonical name when it is there, otherwise the legacy one
// from an older preparation. Readers use it so an unmigrated genome still
// works; anything that writes a dictionary uses DictPath, so what lands on
// disk is always the name GATK itself will look for.
func ResolveDictPath(refFasta string) string {
	canonical := DictPath(refFasta)
	if _, err := os.Stat(canonical); err == nil {
		return canonical
	}
	if legacy := LegacyDictPath(refFasta); legacy != "" {
		if _, err := os.Stat(legacy); err == nil {
			return legacy
		}
	}
	return canonical
}

// EnsureGatkDict reports whether the dictionary GATK itself derives is present
// beside a reference. Stages that go on to invoke GATK with -R call this
// rather than stat'ing ResolveDictPath: a genome carrying only the
// pre-correction genome.fa.gz.dict would satisfy a plain existence check and
// then fail much later, inside the first HaplotypeCaller or BaseRecalibrator,
// under GATK's own less obvious wording. Failing here says which file is
// missing and how to produce it.
func EnsureGatkDict(refFasta string) error {
	canonical := DictPath(refFasta)
	if _, err := os.Stat(canonical); err == nil {
		return nil
	}
	if legacy := LegacyDictPath(refFasta); legacy != "" {
		if _, err := os.Stat(legacy); err == nil {
			return fmt.Errorf("reference dict file %s does not exist: %s is present, but GATK derives the shorter name and will not read it. Run: genome-whisperer PrepareFastaFile -r %s", canonical, legacy, refFasta)
		}
	}
	return fmt.Errorf("reference dict file %s does not exist. Run: genome-whisperer PrepareFastaFile -r %s", canonical, refFasta)
}

// FaiPath returns the samtools index beside a reference. samtools names it
// after the full file name for compressed and uncompressed alike, so unlike
// the dictionary this needs no special case.
func FaiPath(refFasta string) string { return refFasta + ".fai" }

// GziPath returns the bgzip offset index a compressed reference needs before
// anything can seek within it. Empty for an uncompressed reference, which
// needs none.
func GziPath(refFasta string) string {
	if !IsBgzippedFasta(refFasta) {
		return ""
	}
	return refFasta + ".gzi"
}

// HasBlastDB reports whether a nucleotide BLAST database already exists beside
// a reference, under the reference's own full name — genome.fa.gz.nsq beside
// genome.fa.gz, genome.fa.nsq beside genome.fa.
//
// Two shapes count. A database small enough for one volume writes the .nsq
// sequence file directly; one large enough to be split writes per-volume
// .00.nsq, .01.nsq ... plus a .nal alias file naming them, and no plain .nsq
// at all. Checking only for .nsq would therefore rebuild a volumed database on
// every run.
func HasBlastDB(refFasta string) bool {
	for _, ext := range []string{".nsq", ".nal"} {
		if _, err := os.Stat(refFasta + ext); err == nil {
			return true
		}
	}
	return false
}

// HasProtBlastDB reports whether a protein BLAST database exists under a
// prefix — AllResistanceGenes.fasta.psq beside AllResistanceGenes.fasta. As
// with HasBlastDB, a database small enough for one volume writes the .psq
// directly while a split one writes per-volume files plus a .pal alias and no
// plain .psq, so both shapes count.
//
// Unlike HasBlastDB this is not a complete test, and callers must treat a false
// as "cannot confirm" rather than "absent". HasBlastDB appends to a reference
// FASTA because this repo always builds its nucleotide databases as
// -out <the fasta itself>. A blastp -db argument is only a prefix, and a bare
// name is resolved through $BLASTDB, so HasProtBlastDB("swissprot") is false
// for a database blastp finds perfectly. Check it only for a path-shaped db.
func HasProtBlastDB(db string) bool {
	for _, ext := range []string{".psq", ".pal"} {
		if _, err := os.Stat(db + ext); err == nil {
			return true
		}
	}
	return false
}

// LooksLikeBlastDBPath reports whether a -db value names a location on disk
// rather than a bare database name for $BLASTDB to resolve. Only the former can
// be checked with HasProtBlastDB.
func LooksLikeBlastDBPath(db string) bool {
	return strings.ContainsRune(db, os.PathSeparator) || strings.ContainsRune(db, '/')
}
