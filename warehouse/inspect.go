package warehouse

import (
	"io/fs"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/gmaffy/genome-whisperer/alignmentdir"
)

// listDir is the only filesystem read in this package, and the unit of work
// handed to the worker pool.
//
// It is one getdents plus one lstat per entry, and nothing else — no
// validation, no subprocess. On a 9p mount that read is the entire cost of a
// scan (measured: 44 s serially across 345 gVCF directories on one crop, 7.3 s
// at 16 workers), which is why the leaf directory rather than the file is what
// gets parallelised.
type dirListing struct {
	Path    string
	Files   []statFile
	SubDirs []string // names only; the walker decides what to descend into
	Exists  bool
	// Err records a directory that could not be listed. It is a string, not an
	// error: this travels into a JSON row, and an error field marshals to {}.
	Err string
}

// statFile is one entry, with only what a stat already told us.
type statFile struct {
	Name    string
	Bytes   int64
	ModTime time.Time
	Mode    fs.FileMode
}

// IsRegular reports whether the entry is a plain file. A symlink or FIFO in a
// data directory must not be reported as sequencing data.
func (f statFile) IsRegular() bool { return f.Mode.IsRegular() }

// MTime formats the modification time for a report row.
func (f statFile) MTime() string { return f.ModTime.UTC().Format(time.RFC3339) }

// listDir reads one directory. A missing directory is not an error — plenty of
// samples have never been aligned — but an unreadable one is recorded and the
// scan continues.
func listDir(path string) dirListing {
	l := dirListing{Path: path}

	entries, err := os.ReadDir(path)
	if err != nil {
		if os.IsNotExist(err) {
			return l
		}
		l.Exists = true
		l.Err = err.Error()
		return l
	}
	l.Exists = true

	for _, e := range entries {
		name := e.Name()
		if e.IsDir() {
			l.SubDirs = append(l.SubDirs, name)
			continue
		}

		f := statFile{Name: name}
		if info, iErr := e.Info(); iErr == nil {
			f.Bytes = info.Size()
			f.ModTime = info.ModTime()
			f.Mode = info.Mode()
		}
		l.Files = append(l.Files, f)
	}

	return l
}

// names indexes the listing by filename for constant-time sibling lookups.
func (l dirListing) names() map[string]statFile {
	m := make(map[string]statFile, len(l.Files))
	for _, f := range l.Files {
		m[f.Name] = f
	}
	return m
}

// TotalBytes is the size of every regular file in the directory.
func (l dirListing) TotalBytes() int64 {
	var total int64
	for _, f := range l.Files {
		if f.IsRegular() {
			total += f.Bytes
		}
	}
	return total
}

// alignmentIndexFor finds a BAM or CRAM's index among names already read.
//
// It is a map lookup over data already paid for, where the pipeline's own
// findIndex would cost up to four os.Stat calls plus a subprocess per file.
// The trade is that a truncated index reads as present here: this answers
// "does an index exist", not "would GATK accept it".
func alignmentIndexFor(byName map[string]statFile, dataFile string) (string, int64, bool) {
	for _, candidate := range alignmentdir.IndexCandidates(dataFile) {
		if f, ok := byName[filepath.Base(candidate)]; ok && f.Bytes > 0 {
			return f.Name, f.Bytes, true
		}
	}
	return "", 0, false
}

// variantIndexFor finds a VCF, gVCF or BCF's index among names already read.
// Tabix is what the pipeline writes; .csi is accepted because bcftools can
// write it and a hand-indexed file will have one.
func variantIndexFor(byName map[string]statFile, dataFile string) (string, int64, bool) {
	for _, suffix := range []string{".tbi", ".csi"} {
		if f, ok := byName[dataFile+suffix]; ok && f.Bytes > 0 {
			return f.Name, f.Bytes, true
		}
	}
	return "", 0, false
}

// unitMatcher resolves a filename's sequence unit.
//
// Both filename spellings have to be tried. GATK and the joint VCF namer write
// the sanitised label (Cp4_1LG01) while DeepVariant writes the dictionary's own
// spelling (Cp4.1LG01), and a reference whose sequence IDs contain dots — pepo
// v4.1's all do — produces both on the same disk.
//
// Matching is by suffix against the known set, longest first, and never by
// splitting on dots: the reference version inside a joint VCF filename
// ("ucd10xv1.1") carries its own dot, so position counting does not survive
// contact with real names.
type unitMatcher struct {
	// spelling -> canonical unit ID
	bySpelling map[string]string
	// spellings longest-first, so chr1 cannot claim a chr12 file
	ordered []string
	// which spelling each entry was, for reporting
	kindOf map[string]string
}

func newUnitMatcher(units []Unit) *unitMatcher {
	m := &unitMatcher{
		bySpelling: make(map[string]string, len(units)*2),
		kindOf:     make(map[string]string, len(units)*2),
	}
	for _, u := range units {
		m.bySpelling[u.Label] = u.ID
		m.kindOf[u.Label] = SpellingSanitised
		if u.ID != u.Label {
			m.bySpelling[u.ID] = u.ID
			m.kindOf[u.ID] = SpellingRaw
		}
	}
	for s := range m.bySpelling {
		m.ordered = append(m.ordered, s)
	}
	// Longest first. Ties broken by name so the order — and therefore the
	// report — is deterministic.
	sortByLengthDesc(m.ordered)
	return m
}

// match reports the unit a stem ends with, the spelling used, and the opaque
// leading part of the name.
//
// The leading part is returned, never interpreted. It is frequently not the
// sample: a sample directory named MO971_LR holds gVCFs stemmed
// MENINA_LONG_READS.aligned.cram. A gVCF's sample comes from its position in
// the tree and from nothing else.
func (m *unitMatcher) match(stem string) (unitID, spelling, prefix string, ok bool) {
	for _, s := range m.ordered {
		if strings.HasSuffix(stem, "."+s) {
			return m.bySpelling[s], m.kindOf[s], strings.TrimSuffix(stem, "."+s), true
		}
	}
	return "", "", "", false
}

// ambiguous reports unit pairs where one unit's raw spelling equals another's
// sanitised spelling, which would make a filename genuinely undecidable.
func ambiguousSpellings(units []Unit) []string {
	sanitised := map[string]string{}
	for _, u := range units {
		sanitised[u.Label] = u.ID
	}
	var clashes []string
	for _, u := range units {
		if other, ok := sanitised[u.ID]; ok && other != u.ID {
			clashes = append(clashes, u.ID+" vs "+other)
		}
	}
	sortStrings(clashes)
	return clashes
}
