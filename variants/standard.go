package variants

import "strings"

// The directory and filename standard, exported for read-only consumers — the
// warehouse scanner in particular, which inventories the storage arrays and
// must expect exactly what the pipeline writes.
//
// These wrap the unexported implementations rather than restating them, and
// deliberately return the *resolved* answer rather than the raw ingredients: a
// caller handed getChromsAndContigs' two slices still has to re-derive the
// label sanitisation and the "no contigs unit unless there were more than 21
// sequences" rule, which is how an auditor ends up disagreeing with the
// pipeline it is auditing.

// Unit is one thing the pipeline calls variants for: a single chromosome, or
// the synthetic batch of everything too small to get its own job.
//
// ID is what GvcfPath and JointVcfPath are called with. Label is what actually
// lands in the filename, which differs whenever a sequence ID contains a dot.
type Unit struct {
	ID      string    // dict sequence ID, or BatchedLabel
	Label   string    // filename spelling: dots replaced with underscores
	Kind    string    // UnitChrom or UnitContigs
	Length  int       // the sequence's LN, or the sum across a batch
	Members []SeqInfo // nil for a chromosome; every batched sequence otherwise
}

// Unit kinds.
const (
	UnitChrom   = "chrom"
	UnitContigs = "contigs"
)

// BatchedLabel is the single label every sequence too small for its own job is
// called under. It contains no dot, so its ID and Label are the same.
const BatchedLabel = "contigs"

// UnitRule names the rule ExpectedUnits implements, for a report to record
// alongside the units so a later change is visible rather than silent.
const UnitRule = "top21-by-length+MT+Pltd+contigs"

// ExpectedUnits returns every unit a reference produces, in the order
// CreateGvcfs enqueues them: the individually-called chromosomes first, then
// the batched group if there is one.
//
// There is no "contigs" unit when the dict holds 21 sequences or fewer — every
// sequence gets its own job — so the expected count is len(units), never
// len(chroms)+1.
func ExpectedUnits(dictFilePath string) ([]Unit, error) {
	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return nil, err
	}

	units := make([]Unit, 0, len(chroms)+1)
	for _, c := range chroms {
		units = append(units, Unit{
			ID:     c.ID,
			Label:  labelFor(c.ID),
			Kind:   UnitChrom,
			Length: c.Len,
		})
	}

	if len(contigs) > 0 {
		total := 0
		for _, c := range contigs {
			total += c.Len
		}
		units = append(units, Unit{
			ID:      BatchedLabel,
			Label:   BatchedLabel,
			Kind:    UnitContigs,
			Length:  total,
			Members: contigs,
		})
	}

	return units, nil
}

// labelFor is the filename spelling of a sequence ID. It mirrors the
// sanitisation in GvcfPath and JointVcfPath.
//
// Note it is one-way and lossy: DeepVariant gVCFs already on disk keep the
// dotted spelling in their filenames, so anything matching files must try both
// Unit.Label and Unit.ID.
func labelFor(chrom string) string {
	return strings.ReplaceAll(chrom, ".", "_")
}

// GvcfDirNames returns the per-sample gVCF directory names, one per caller, as
// GvcfPath chooses them. A directory named anything else is not standard.
func GvcfDirNames() []string {
	return []string{"gatk_gvcfs", "dv_gvcfs"}
}

// CallerMergerTags returns every caller/merger combination directory name, as
// callerMergerTag builds them. Each also appears inside the filenames within.
func CallerMergerTags() []string {
	return []string{"gatk_gatk", "gatk_glnexus", "dv_glnexus"}
}
