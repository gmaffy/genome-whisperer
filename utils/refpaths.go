package utils

import (
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
// An uncompressed genome.fa keeps the GATK default of genome.dict — the
// extension is replaced. A bgzipped genome.fa.gz instead gets
// genome.fa.gz.dict: the whole file name plus the suffix, the same shape
// samtools already uses for the .fai and .gzi beside it. Both forms of one
// assembly can then sit in a single assembly directory without sharing, or
// silently fighting over, one dictionary.
//
// A name that is not a recognised FASTA falls back to replacing its last
// extension, which is what every call site did before this existed.
func DictPath(refFasta string) string {
	if IsBgzippedFasta(refFasta) {
		return refFasta + ".dict"
	}
	if ext := FastaExt(refFasta); ext != "" {
		return strings.TrimSuffix(refFasta, ext) + ".dict"
	}
	return strings.TrimSuffix(refFasta, filepath.Ext(refFasta)) + ".dict"
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
