package warehouse

import (
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// Containers a joint VCF can be stored in.
const (
	ContainerVcfGz = "vcf.gz"
	// ContainerBcf is what GLnexus writes natively. A *.joint.vcf.gz glob
	// finds nothing in such a directory, which is how a whole DeepVariant
	// cohort can appear never to have been merged.
	ContainerBcf = "bcf"
)

// inspectJointDir turns one <caller>_<merger> directory into a report row.
//
// The annotation chain is reported as booleans rather than as file rows,
// because it is a fixed sequence applied to one file: the concatenated
// ".all" VCF is hard-filtered and then annotated, and nothing is annotated per
// chromosome.
func inspectJointDir(l dirListing, volume, crop, reference, tag, caller, merger, layout string,
	matcher *unitMatcher, expected []Unit, validate, quick bool) (JointSet, []Nonconformance) {

	set := JointSet{
		Kind: KindJointSet, Volume: volume, Crop: crop, Reference: reference,
		Tag: tag, Caller: caller, Merger: merger, Layout: layout, Dir: l.Path,
		Container: ContainerVcfGz,
	}
	var problems []Nonconformance

	byName := l.names()
	present := map[string]bool{}
	indexed := map[string]bool{}
	sawBcf := false

	for _, f := range l.Files {
		if !f.IsRegular() {
			continue
		}
		name := f.Name

		// Indexes are accounted for through their data file.
		if strings.HasSuffix(name, ".tbi") || strings.HasSuffix(name, ".csi") {
			continue
		}
		if validate && (strings.HasSuffix(name, SuffixJointVcf) ||
			strings.HasSuffix(name, SuffixJointBcf) ||
			strings.HasSuffix(name, ".all"+SuffixHardFilt) ||
			strings.HasSuffix(name, ".all"+SuffixVcf)) {
			path := filepath.Join(l.Path, name)
			if err := utils.ValidateGvcf(path, false, quick); err != nil {
				problems = append(problems, Nonconformance{
					Kind: KindNonconformance, Rule: RuleJointVcfValidationFailed,
					Severity: SeverityError, Volume: volume, Crop: crop, Reference: reference,
					Path: path, Detail: err.Error(),
				})
			}
		}

		// The annotation chain, longest suffix first: every one of these also
		// ends in ".tsv", and the hard-filtered VCF also ends in ".vcf.gz".
		switch {
		case strings.HasSuffix(name, SuffixSuperTsv):
			set.HasSuperVcfTsv = true
			continue
		case strings.HasSuffix(name, SuffixDescTsv):
			set.HasDescTsv = true
			continue
		case strings.HasSuffix(name, SuffixPrgTsv):
			continue
		case strings.HasSuffix(name, SuffixEffTsv):
			set.HasEffTsv = true
			continue
		case strings.HasSuffix(name, SuffixSnpEffTsv):
			set.HasSnpEffTsv = true
			continue
		case strings.HasSuffix(name, SuffixSnpEffVcf):
			// Never compressed and never indexed by the pipeline.
			set.HasSnpEffVcf = true
			continue
		case strings.HasSuffix(name, ".all"+SuffixHardFilt):
			set.HasHardFiltered = true
			set.JointBytes += f.Bytes
			_, _, set.HasHardFilteredIndex = variantIndexFor(byName, name)
			continue
		case strings.HasSuffix(name, ".all"+SuffixVcf):
			set.HasAllVcf = true
			set.JointBytes += f.Bytes
			_, _, set.HasAllIndex = variantIndexFor(byName, name)
			continue
		}

		// A per-unit joint call.
		var stem string
		switch {
		case strings.HasSuffix(name, SuffixJointVcf):
			stem = strings.TrimSuffix(name, SuffixJointVcf)
		case strings.HasSuffix(name, SuffixJointBcf):
			stem = strings.TrimSuffix(name, SuffixJointBcf)
			sawBcf = true
		default:
			set.OtherFileCount++
			continue
		}

		set.JointBytes += f.Bytes
		if matcher == nil {
			set.PresentCount++
			continue
		}

		unitID, _, _, ok := matcher.match(stem)
		if !ok {
			set.OtherFileCount++
			continue
		}
		present[unitID] = true
		if _, _, hasIndex := variantIndexFor(byName, name); hasIndex {
			indexed[unitID] = true
		}
	}

	if sawBcf {
		set.Container = ContainerBcf
	}

	if matcher != nil {
		set.ExpectedUnitCount = len(expected)
		set.PresentCount = len(present)
		set.IndexedCount = len(indexed)
		for _, u := range expected {
			if !present[u.ID] {
				set.Missing = append(set.Missing, u.Label)
			}
		}
		sortStrings(set.Missing)
	}

	return set, problems
}
