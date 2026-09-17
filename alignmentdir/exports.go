package alignmentdir

import "strings"

// Pure classifiers exported for read-only consumers — the warehouse scanner in
// particular, which must recognise the same files the pipeline does.
//
// Everything here is a string function: no stat, no subprocess. That line is
// deliberate. The inspection helpers in this package are the opposite —
// getAlignmentFileInfo runs samtools per file and findIndex decompresses each
// index — which is right before handing a file to GATK and far too slow for
// walking an estate of tens of thousands of files. A scanner needs the naming
// rules, not the validation, so only the naming rules are exported.

// Alignment roles, as classified by AlignmentRole.
const (
	// RoleIgnore is a file that belongs in a bams directory but is not
	// alignment data: an index, a metrics table, an interval list, a plot.
	RoleIgnore = ""
	// RoleOther is a file the pipeline neither writes nor recognises.
	RoleOther     = "other"
	RoleSortedBam = "sorted_bam"
	RoleRgmdBam   = "rgmd_bam"
	RoleRgmdCram  = "rgmd_cram"
	RoleBqsrBam   = "bqsr_bam"
	RoleBqsrCram  = "bqsr_cram"
)

// AlignmentRole reports what a filename in a bams directory is, using the same
// rules and the same precedence order inspectSampleBamDir applies.
func AlignmentRole(name string) string {
	switch classifyAlignmentFile(name) {
	case slotSortedBam:
		return RoleSortedBam
	case slotRgmdBam:
		return RoleRgmdBam
	case slotRgmdCram:
		return RoleRgmdCram
	case slotBqsrBam:
		return RoleBqsrBam
	case slotBqsrCram:
		return RoleBqsrCram
	case slotOther:
		return RoleOther
	default:
		return RoleIgnore
	}
}

// ReadRole is the role a FASTQ file plays in a clean_reads directory.
type ReadRole int

const (
	// ReadRoleUnknown covers a file matching neither convention. GetReadsPE
	// drops these silently; a scanner has to report them, which is why this
	// classifier is exported separately from the collector.
	ReadRoleUnknown ReadRole = iota
	ReadRoleForward
	ReadRoleReverse
)

func (r ReadRole) String() string {
	switch r {
	case ReadRoleForward:
		return "fwd"
	case ReadRoleReverse:
		return "rev"
	default:
		return "unknown"
	}
}

// ClassifyRead reports whether a clean_reads filename is a forward or reverse
// read, with the same precedence as GetReadsPE: forward is tested first, so a
// name matching both predicates is forward.
//
// A single-end or long-read FASTQ matches neither and returns ReadRoleUnknown;
// that is a real layout, not an error, so a caller must not read "no forward
// and no reverse" as "no reads".
func ClassifyRead(filename string) ReadRole {
	if isFwd(filename) {
		return ReadRoleForward
	}
	if isRev(filename) {
		return ReadRoleReverse
	}
	return ReadRoleUnknown
}

// IsLongReadSample reports whether a sample name marks a long-read sample.
//
// Long reads are identified by the sample directory name alone, never by file
// contents. The same test is spelled three separate ways in this codebase
// (dirAlign.go, and twice in variants/createGvcfs.go), two upper-casing and one
// lower-casing; they agree on every real name, and this is the spelling new
// code should use.
func IsLongReadSample(sample string) bool {
	return strings.HasSuffix(strings.ToUpper(sample), "LR")
}

// IndexCandidates lists every index filename an alignment file might carry,
// current and historical, in the order the pipeline prefers them.
//
// This is the naming knowledge without the cost: findIndex additionally reads
// each candidate to reject one left truncated by an interrupted run. A caller
// that only stats these names learns that an index exists, not that it is
// usable.
func IndexCandidates(path string) []string {
	return indexCandidates(path)
}

// alignmentSlot names the field of SampleBamState a filename belongs to.
type alignmentSlot int

const (
	slotIgnore alignmentSlot = iota
	slotOther
	slotSortedBam
	slotRgmdBam
	slotRgmdCram
	slotBqsrBam
	slotBqsrCram
)

// classifyAlignmentFile maps a filename in a bams directory to its slot.
//
// The case order is load-bearing and must not be reordered: "sorted.bam" is
// tested before "rgmd.bam" so a *.RGMD.sorted.bam counts as the sorted output,
// and "rgmd.bam" before "bqsr.bam" so *.RGMD_bqsr.bam — which does not end in
// "rgmd.bam" — falls through to the BQSR slot.
func classifyAlignmentFile(name string) alignmentSlot {
	lowerName := strings.ToLower(name)
	switch {
	case strings.HasSuffix(lowerName, "sorted.bam"):
		return slotSortedBam
	case strings.HasSuffix(lowerName, "rgmd.bam"):
		return slotRgmdBam
	case strings.HasSuffix(lowerName, "rgmd.cram"):
		return slotRgmdCram
	case strings.HasSuffix(lowerName, "bqsr.bam"):
		return slotBqsrBam
	case strings.HasSuffix(lowerName, "bqsr.cram"):
		return slotBqsrCram
	case strings.HasSuffix(lowerName, ".bai"), strings.HasSuffix(lowerName, ".csi"), strings.HasSuffix(lowerName, ".crai"),
		strings.HasSuffix(lowerName, ".pdf"), strings.HasSuffix(lowerName, ".txt"), strings.HasSuffix(lowerName, ".list"):
		return slotIgnore
	default:
		return slotOther
	}
}

// slotField returns the field a slot writes to, or nil for slots that are not
// a single alignment file.
func (s *SampleBamState) slotField(slot alignmentSlot) *FileInfo {
	switch slot {
	case slotSortedBam:
		return &s.SortedBam
	case slotRgmdBam:
		return &s.RgmdBam
	case slotRgmdCram:
		return &s.RgmdCram
	case slotBqsrBam:
		return &s.BqsrBam
	case slotBqsrCram:
		return &s.BqsrCram
	}
	return nil
}
