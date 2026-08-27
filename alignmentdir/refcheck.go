package alignmentdir

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// sqRecord is one @SQ line reduced to what identifies a contig: the name it is
// known by, and the sequence itself.
type sqRecord struct {
	name   string
	length string
	m5     string
}

// parseSQ pulls the @SQ lines out of SAM-format header text, in file order.
func parseSQ(header string) []sqRecord {
	var records []sqRecord
	for _, line := range strings.Split(header, "\n") {
		if !strings.HasPrefix(line, "@SQ\t") {
			continue
		}
		var rec sqRecord
		for _, field := range strings.Split(line, "\t") {
			switch {
			case strings.HasPrefix(field, "SN:"):
				rec.name = field[3:]
			case strings.HasPrefix(field, "LN:"):
				rec.length = field[3:]
			case strings.HasPrefix(field, "M5:"):
				rec.m5 = field[3:]
			}
		}
		if rec.name != "" {
			records = append(records, rec)
		}
	}
	return records
}

// alignmentSQ reads an alignment file's @SQ lines. Only the header is read, so
// this costs nothing and — unlike decoding — works on a CRAM whose reference
// cannot be found, which is exactly the case this is here to diagnose.
func alignmentSQ(path string) ([]sqRecord, error) {
	out, err := exec.Command("samtools", "view", "-H", path).Output()
	if err != nil {
		return nil, err
	}
	return parseSQ(string(out)), nil
}

// referenceSQ reads the .dict beside a reference FASTA.
func referenceSQ(refFasta string) ([]sqRecord, error) {
	dictPath := utils.DictPath(refFasta)
	data, err := os.ReadFile(dictPath)
	if err != nil {
		return nil, err
	}
	return parseSQ(string(data)), nil
}

// sameSequences reports whether two dictionaries describe the same sequences in
// the same order, whatever they call them. An M5 is the checksum of the
// sequence itself, so matching M5s mean a rename and nothing more.
func sameSequences(a, b []sqRecord) bool {
	if len(a) != len(b) || len(a) == 0 {
		return false
	}
	for i := range a {
		if a[i].m5 == "" || b[i].m5 == "" {
			return false
		}
		if a[i].m5 != b[i].m5 || a[i].length != b[i].length {
			return false
		}
	}
	return true
}

// checkAlignmentContigs reports an error when an alignment file's contigs are
// not the reference's.
//
// A CRAM resolves its reference by contig name. An assembly renamed after the
// alignments were made therefore produces files that samtools cannot decode —
// it falls back to the UR path in the header, which by then is usually gone —
// and that GATK rejects for having an incompatible sequence dictionary. Neither
// says which of the two references is the odd one out, and `samtools quickcheck`
// (all that --quick runs) never looks at contigs at all, so without this the
// mismatch surfaces hours later as a missing-sequence error from inside htslib.
//
// A file this cannot read an answer for is left alone: the point is to name a
// problem that is definitely there, not to invent one.
func checkAlignmentContigs(path, refFasta string) error {
	fileSQ, err := alignmentSQ(path)
	if err != nil || len(fileSQ) == 0 {
		return nil
	}
	refSQ, err := referenceSQ(refFasta)
	if err != nil || len(refSQ) == 0 {
		return nil
	}

	known := make(map[string]struct{}, len(refSQ))
	for _, sq := range refSQ {
		known[sq.name] = struct{}{}
	}

	var unknown []string
	for _, sq := range fileSQ {
		if _, ok := known[sq.name]; !ok {
			unknown = append(unknown, sq.name)
			if len(unknown) == 3 {
				break
			}
		}
	}
	if len(unknown) == 0 {
		return nil
	}

	if sameSequences(fileSQ, refSQ) {
		return fmt.Errorf("contigs are named %s… here but %s… in %s — the same %d sequences under different names "+
			"(every M5 checksum and length matches, in order). Rename this file's header with: "+
			"scripts/reheader-crams.sh -r %s %s",
			strings.Join(unknown, ", "), refSQ[0].name, filepath.Base(refFasta), len(refSQ), refFasta, path)
	}
	return fmt.Errorf("contigs %s… are not in %s — this file was aligned to a different assembly",
		strings.Join(unknown, ", "), filepath.Base(refFasta))
}
