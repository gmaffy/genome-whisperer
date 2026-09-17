package warehouse

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
)

// genomeIndex maps crop and reference names to their verbatim spelling under
// the genomes directory.
//
// It is built here rather than taken from utils.GetValidGenomesFromDisk
// because that function upper-cases both names, which is exactly the
// information a case-drift audit needs: this estate has sample directories
// naming a reference "c39" where the genomes directory spells it "C39". The
// pipeline resolves that fine — it upper-cases too — which is why the drift has
// survived unnoticed, and why reporting it needs both spellings intact.
type genomeIndex struct {
	dir string
	// upper(crop) -> verbatim crop
	crops map[string]string
	// upper(crop) -> upper(ref) -> verbatim ref
	refs map[string]map[string]string
}

func newGenomeIndex(genomesDir string) *genomeIndex {
	idx := &genomeIndex{
		dir:   genomesDir,
		crops: map[string]string{},
		refs:  map[string]map[string]string{},
	}
	if genomesDir == "" {
		return idx
	}

	cropEntries, err := os.ReadDir(genomesDir)
	if err != nil {
		return idx
	}
	for _, c := range cropEntries {
		if !c.IsDir() || IsPruned(c.Name()) {
			continue
		}
		cropKey := strings.ToUpper(c.Name())
		idx.crops[cropKey] = c.Name()
		idx.refs[cropKey] = map[string]string{}

		refEntries, rErr := os.ReadDir(filepath.Join(genomesDir, c.Name()))
		if rErr != nil {
			continue
		}
		for _, r := range refEntries {
			if !r.IsDir() || IsPruned(r.Name()) {
				continue
			}
			idx.refs[cropKey][strings.ToUpper(r.Name())] = r.Name()
		}
	}
	return idx
}

// verbatimRef returns the genomes directory's own spelling of a reference.
func (idx *genomeIndex) verbatimRef(crop, ref string) (cropName, refName string, ok bool) {
	cropKey := strings.ToUpper(crop)
	cropName, ok = idx.crops[cropKey]
	if !ok {
		return "", "", false
	}
	refName, ok = idx.refs[cropKey][strings.ToUpper(ref)]
	if !ok {
		return cropName, "", false
	}
	return cropName, refName, true
}

// dictFor resolves the reference dictionary for a crop and reference as spelled
// in a data directory.
//
// It returns a status rather than an error, because "there is no dictionary" is
// a normal thing to report and a fatal thing to raise: a read-only audit must
// describe an unresolvable reference, not abort on it. The pipeline's
// utils.EnsureGatkDict is deliberately not used for that reason.
func (idx *genomeIndex) dictFor(crop, ref string) (dictPath, fastaPath, verbatim, status string) {
	cropName, refName, ok := idx.verbatimRef(crop, ref)
	if !ok {
		return "", "", refName, DictUnknownRef
	}

	assembly := filepath.Join(idx.dir, cropName, refName, "assembly")
	entries, err := os.ReadDir(assembly)
	if err != nil {
		return "", "", refName, DictMissing
	}

	// Prefer the name GATK derives from the FASTA, which is what
	// utils.ResolveDictPath picks; fall back to any .dict present.
	var fasta, anyDict string
	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		name := e.Name()
		if strings.HasSuffix(strings.ToLower(name), ".dict") && anyDict == "" {
			anyDict = filepath.Join(assembly, name)
		}
		if fasta == "" && utils.IsFasta(name) {
			fasta = filepath.Join(assembly, name)
		}
	}

	if fasta != "" {
		if resolved := utils.ResolveDictPath(fasta); resolved != "" {
			if _, sErr := os.Stat(resolved); sErr == nil {
				return resolved, fasta, refName, DictResolved
			}
		}
	}
	if anyDict != "" {
		return anyDict, fasta, refName, DictResolved
	}
	return "", fasta, refName, DictMissing
}

// referenceRow builds the reference row for a crop and reference, resolving the
// expected units when a dictionary is available.
//
// Units are omitted — not emptied — whenever the dictionary did not resolve.
// The distinction is the whole point: an empty expected set that reads as
// "satisfied" is how a previous inventory of this estate reported complete
// samples as 81,202 gVCFs short.
func (idx *genomeIndex) referenceRow(volume, crop, ref string) (Reference, []Unit, *ReferenceBatch, string, []Nonconformance) {
	var problems []Nonconformance

	dictPath, fastaPath, verbatim, status := idx.dictFor(crop, ref)
	row := Reference{
		Kind:       KindReference,
		Volume:     volume,
		Crop:       crop,
		Reference:  ref,
		Spelling:   ref,
		DictStatus: status,
		DictPath:   dictPath,
		GenomesRef: verbatim,
		Rule:       variants.UnitRule,
	}

	// The reference exists in the genomes directory under a different case.
	// Harmless on a case-insensitive mount, and a broken path the moment this
	// data is copied somewhere case-sensitive.
	if verbatim != "" && verbatim != ref {
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleReferenceCaseDrift,
			Severity: SeverityWarning, Volume: volume, Crop: crop, Reference: ref,
			Found: ref, Expected: verbatim,
			Detail: "data directory and genomes directory disagree on capitalisation",
		})
	}

	switch status {
	case DictUnknownRef:
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleReferenceNotInGenomesDir,
			Severity: SeverityWarning, Volume: volume, Crop: crop, Reference: ref,
			Path:   filepath.Join(idx.dir, crop, ref),
			Detail: "no such reference under the genomes directory; completeness cannot be judged",
		})
		return row, nil, nil, fastaPath, problems
	case DictMissing:
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleDictMissing,
			Severity: SeverityWarning, Volume: volume, Crop: crop, Reference: ref,
			Path:   filepath.Join(idx.dir, crop, verbatim, "assembly"),
			Detail: "no .dict in the assembly directory; completeness cannot be judged",
		})
		return row, nil, nil, fastaPath, problems
	}

	units, err := variants.ExpectedUnits(dictPath)
	if err != nil {
		row.DictStatus = DictUnreadable
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleDictUnreadable,
			Severity: SeverityError, Volume: volume, Crop: crop, Reference: ref,
			Path: dictPath, Detail: err.Error(),
		})
		return row, nil, nil, fastaPath, problems
	}

	rows := make([]Unit, 0, len(units))
	var batch *ReferenceBatch
	seqCount := 0
	for _, u := range units {
		members := len(u.Members)
		if members == 0 {
			members = 1
		}
		seqCount += members
		rows = append(rows, Unit{
			ID:          u.ID,
			Label:       u.Label,
			UnitKind:    u.Kind,
			Length:      u.Length,
			MemberCount: members,
		})

		if u.Kind == variants.UnitContigs {
			seqIDs := make([]string, 0, len(u.Members))
			for _, m := range u.Members {
				seqIDs = append(seqIDs, m.ID)
			}
			sortStrings(seqIDs)
			batch = &ReferenceBatch{
				Kind: KindReferenceBatch, Crop: crop, Reference: ref,
				MemberCount: len(u.Members), TotalLength: u.Length, SeqIDs: seqIDs,
			}
		}
	}

	row.SeqCount = seqCount
	row.Units = rows
	row.ExpectedUnitCount = len(rows)

	// A dictionary holding both "X.1" and "X_1" would make a filename
	// genuinely undecidable between the two spellings.
	if clashes := ambiguousSpellings(rows); len(clashes) > 0 {
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleAmbiguousSeqLabel,
			Severity: SeverityError, Volume: volume, Crop: crop, Reference: ref,
			Path: dictPath, Count: len(clashes),
			Detail: "sequence IDs collide once dots become underscores: " + strings.Join(clashes, ", "),
		})
	}

	return row, rows, batch, fastaPath, problems
}
