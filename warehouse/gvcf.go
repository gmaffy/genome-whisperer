package warehouse

import (
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// variantFileKind is what a file in a gVCF directory actually is.
type variantFileKind int

const (
	fileOther variantFileKind = iota
	fileGvcf
	fileGvcfIndex
	// fileSidecarVcf is DeepVariant's plain VCF, written beside every gVCF.
	// It is not a gVCF, and a *.vcf.gz glob would double every count.
	fileSidecarVcf
	fileSidecarVcfIndex
	// fileSidecarReport is DeepVariant's per-region HTML report.
	fileSidecarReport
)

// classifyVariantFile dispatches on the longest matching suffix.
//
// Order is the whole correctness argument. Every name ending ".g.vcf.gz.tbi"
// also ends ".vcf.gz.tbi", and every ".g.vcf.gz" also ends ".vcf.gz", so the
// gVCF forms have to be tested first or a directory of gVCFs reads as a
// directory of sidecars.
func classifyVariantFile(name string) (variantFileKind, string) {
	switch {
	case strings.HasSuffix(name, SuffixGvcfTbi):
		return fileGvcfIndex, strings.TrimSuffix(name, SuffixGvcfTbi)
	case strings.HasSuffix(name, SuffixGvcfCsi):
		return fileGvcfIndex, strings.TrimSuffix(name, SuffixGvcfCsi)
	case strings.HasSuffix(name, SuffixGvcf):
		return fileGvcf, strings.TrimSuffix(name, SuffixGvcf)
	case strings.HasSuffix(name, SuffixVcfTbi):
		return fileSidecarVcfIndex, strings.TrimSuffix(name, SuffixVcfTbi)
	case strings.HasSuffix(name, SuffixVcfCsi):
		return fileSidecarVcfIndex, strings.TrimSuffix(name, SuffixVcfCsi)
	case strings.HasSuffix(name, SuffixVcf):
		return fileSidecarVcf, strings.TrimSuffix(name, SuffixVcf)
	case strings.HasSuffix(name, SuffixDVReport):
		return fileSidecarReport, strings.TrimSuffix(name, SuffixDVReport)
	}
	return fileOther, name
}

// inspectGvcfDir turns one gVCF directory listing into a report row.
//
// matcher is nil when the reference's dictionary did not resolve. The gVCFs are
// still counted and their bytes still totalled — there is no reason to hide
// files just because the denominator is unknown — but ExpectedUnitCount and
// Missing are left unset, so nothing downstream can compute a fraction against
// a denominator that was never established.
func inspectGvcfDir(l dirListing, dirName, caller string, matcher *unitMatcher, expected []Unit,
	validate, quick bool) (GvcfSet, []Nonconformance) {
	set := GvcfSet{Dir: dirName, Caller: caller}
	var problems []Nonconformance

	byName := l.names()

	present := map[string]bool{}
	indexed := map[string]bool{}
	var stems, spellings, unmatched []string

	for _, f := range l.Files {
		if !f.IsRegular() {
			continue
		}

		kind, _ := classifyVariantFile(f.Name)
		switch kind {
		case fileSidecarVcf:
			set.SidecarVcfCount++
			if validate {
				path := filepath.Join(l.Path, f.Name)
				if err := utils.ValidateGvcf(path, false, quick); err != nil {
					problems = append(problems, Nonconformance{
						Kind: KindNonconformance, Rule: RuleGvcfValidationFailed,
						Severity: SeverityError, Path: path, Detail: err.Error(),
					})
				}
			}
			continue
		case fileSidecarVcfIndex:
			continue
		case fileSidecarReport:
			set.SidecarReportCount++
			continue
		case fileGvcfIndex:
			continue
		case fileOther:
			set.UnrecognisedCount++
			continue
		}

		// A gVCF.
		set.Bytes += f.Bytes
		if validate {
			path := filepath.Join(l.Path, f.Name)
			if err := utils.ValidateGvcf(path, false, quick); err != nil {
				problems = append(problems, Nonconformance{
					Kind: KindNonconformance, Rule: RuleGvcfValidationFailed,
					Severity: SeverityError, Path: path, Detail: err.Error(),
				})
			}
		}
		stem := strings.TrimSuffix(f.Name, SuffixGvcf)

		if matcher == nil {
			set.PresentCount++
			continue
		}

		unitID, spelling, prefix, ok := matcher.match(stem)
		if !ok {
			unmatched = append(unmatched, f.Name)
			set.UnrecognisedCount++
			continue
		}

		present[unitID] = true
		spellings = append(spellings, spelling)
		if prefix != "" {
			stems = append(stems, prefix)
		}

		if _, _, hasIndex := variantIndexFor(byName, f.Name); hasIndex {
			indexed[unitID] = true
		}
	}

	set.Stems = uniqueSorted(stems)
	set.LabelSpellings = uniqueSorted(spellings)

	if matcher != nil {
		set.ExpectedUnitCount = len(expected)
		set.PresentCount = len(present)
		set.IndexedCount = len(indexed)
		for _, u := range expected {
			if !present[u.ID] {
				set.Missing = append(set.Missing, u.Label)
			} else if !indexed[u.ID] {
				set.Unindexed = append(set.Unindexed, u.Label)
			}
		}
		sortStrings(set.Missing)
		sortStrings(set.Unindexed)
	}

	if len(unmatched) > 0 {
		sortStrings(unmatched)
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnrecognisedSeqLabel,
			Severity: SeverityWarning, Path: l.Path, Count: len(unmatched),
			Detail: "gVCFs naming no expected sequence: " + strings.Join(truncate(unmatched, 5), ", "),
		})
	}

	// More than one stem means the directory mixes provenance: gVCFs from
	// different alignments, or leftovers from a renamed sample. The stem is
	// not the sample name and is never read as one, but it is a useful signal
	// that two runs wrote here.
	if len(set.Stems) > 1 {
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleMultipleGvcfStems,
			Severity: SeverityWarning, Path: l.Path, Count: len(set.Stems),
			Detail: "gVCFs from more than one alignment: " + strings.Join(truncate(set.Stems, 5), ", "),
		})
	}

	if len(set.Unindexed) > 0 {
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleGvcfMissingIndex,
			Severity: SeverityWarning, Path: l.Path, Count: len(set.Unindexed),
			Detail: "gVCFs with no .tbi: " + strings.Join(truncate(set.Unindexed, 5), ", "),
		})
	}

	return set, problems
}

// countLegacyGvcfs counts gVCFs in a non-standard directory without
// inventorying them.
//
// Legacy directories are deliberately not inventoried, but they must be
// counted: without the count a crop that was fully called into the old layout
// is indistinguishable from one that was never called at all.
func countLegacyGvcfs(l dirListing) (files int, bytes int64) {
	for _, f := range l.Files {
		if f.IsRegular() && strings.HasSuffix(f.Name, SuffixGvcf) {
			files++
			bytes += f.Bytes
		}
	}
	return files, bytes
}

// truncate caps a list for a human-readable detail string.
func truncate(s []string, n int) []string {
	if len(s) <= n {
		return s
	}
	out := append([]string{}, s[:n]...)
	return append(out, "...")
}
