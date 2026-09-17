package warehouse

import (
	"regexp"
	"strings"

	"github.com/gmaffy/genome-whisperer/variants"
)

// The layout standard, in one place.
//
// Standard tree:
//
//	<mount>/<data>/<crop>/<year>/<sample>/clean_reads/
//	                                      reference_genomes/<ref>/bams/
//	                                                             /gatk_gvcfs/
//	                                                             /dv_gvcfs/
//	<mount>/<data>/<crop>/MERGED_VCFs/<ref>/<caller>_<merger>/
//
// Only these shapes are inventoried. Older shapes are recognised well enough to
// be counted and reported, never read for contents — so a crop still on an old
// layout reads as "data in the wrong place", not as "no data".

// Directory names at the sample and reference levels.
const (
	DirCleanReads       = "clean_reads"
	DirReferenceGenomes = "reference_genomes"
	DirBams             = "bams"
	DirMergedVcfs       = "MERGED_VCFs"

	// DirLegacyGvcfs predates the per-caller split into gatk_gvcfs/dv_gvcfs.
	DirLegacyGvcfs = "gvcfs"
	// DirLegacyVcfs is the older joint VCF root, both as <crop>/VCFs/<ref>/
	// and as <crop>/<ref>/VCFs/<tag>/.
	DirLegacyVcfs = "VCFs"
)

// yearPattern matches a cohort directory. Everything else at that level is
// something we have to classify rather than assume.
var yearPattern = regexp.MustCompile(`^(19|20)\d{2}$`)

// IsYear reports whether a directory name is a cohort year.
func IsYear(name string) bool { return yearPattern.MatchString(name) }

// GvcfDirCaller maps a standard gVCF directory name to its caller, and reports
// whether the name is standard at all.
func GvcfDirCaller(name string) (string, bool) {
	switch name {
	case "gatk_gvcfs":
		return "gatk", true
	case "dv_gvcfs":
		return "dv", true
	}
	return "", false
}

// SplitCallerMergerTag splits a joint VCF directory name into its caller and
// merger, and reports whether it is one of the standard combinations.
func SplitCallerMergerTag(name string) (caller, merger string, ok bool) {
	for _, tag := range variants.CallerMergerTags() {
		if tag == name {
			parts := strings.SplitN(tag, "_", 2)
			return parts[0], parts[1], true
		}
	}
	return "", "", false
}

// prunedNames are directories never descended into: version-control and
// snapshot machinery, recycle bins, scratch, and the pipeline's own temporary
// output. Without this a walk wanders into a checked-out git repository that
// happens to sit under a crop directory.
var prunedNames = map[string]bool{
	"lost+found": true,
	"#recycle":   true,
	"work":       true,
	"shards":     true,
	"QC":         true,
}

// IsPruned reports whether a directory should not be descended into.
func IsPruned(name string) bool {
	if name == "" {
		return true
	}
	// Dotfiles (.git, .snapshot) and QNAP/Synology volume metadata (@Recycle,
	// @Recently-Snapshot) appear at share roots.
	if name[0] == '.' || name[0] == '@' {
		return true
	}
	// DeepVariant intermediate directories, left behind by a crashed run.
	if strings.HasPrefix(name, "tmp_") {
		return true
	}
	return prunedNames[name]
}

// File suffixes. Order matters wherever these are tried: the longest match has
// to win, or a gVCF index is mistaken for a gVCF.
const (
	SuffixGvcf      = ".g.vcf.gz"
	SuffixGvcfTbi   = ".g.vcf.gz.tbi"
	SuffixGvcfCsi   = ".g.vcf.gz.csi"
	SuffixVcf       = ".vcf.gz"
	SuffixVcfTbi    = ".vcf.gz.tbi"
	SuffixVcfCsi    = ".vcf.gz.csi"
	SuffixBcf       = ".bcf"
	SuffixBcfCsi    = ".bcf.csi"
	SuffixDVReport  = ".visual_report.html"
	SuffixJointVcf  = ".joint.vcf.gz"
	SuffixJointBcf  = ".joint.bcf"
	SuffixHardFilt  = ".hard_filtered.vcf.gz"
	SuffixSnpEffVcf = ".snpEff.vcf"
	SuffixSnpEffTsv = ".snpEff.tsv"
	SuffixEffTsv    = ".snpEff_EFF.tsv"
	SuffixDescTsv   = ".snpEff_EFF_DESC.tsv"
	SuffixPrgTsv    = ".snpEff_EFF_PRG.tsv"
	SuffixSuperTsv  = ".snpEff_EFF_SUPER_VCF.tsv"
)

// fastqSuffixes recognise a read file regardless of pairing convention, so a
// single-end or long-read FASTQ that matches neither the forward nor the
// reverse predicate is still counted as reads.
var fastqSuffixes = []string{".fastq.gz", ".fq.gz", ".fastq", ".fq"}

// IsFastq reports whether a filename is a read file.
func IsFastq(name string) bool {
	lower := strings.ToLower(name)
	for _, s := range fastqSuffixes {
		if strings.HasSuffix(lower, s) {
			return true
		}
	}
	return false
}

// Nonconformance rules. This is a closed enum: a consumer switches on it, so
// new deviations get a new constant rather than a free-text message.
const (
	RuleStrayFileAtCropPosition   = "stray_file_at_crop_position"
	RuleStrayFileAtYearPosition   = "stray_file_at_year_position"
	RuleStrayFileAtSamplePosition = "stray_file_at_sample_position"
	RuleStrayFileInCleanReads     = "stray_file_in_clean_reads"

	RuleForeignDirAtYearPosition = "foreign_dir_at_year_position"
	RuleSampleAtYearPosition     = "sample_at_year_position"

	RuleLegacyJointRoot         = "legacy_joint_root"
	RuleLegacyJointRootRefFirst = "legacy_joint_root_ref_first"
	RuleLegacyGvcfDir           = "legacy_gvcf_dir"

	RuleMissingCleanReads       = "missing_clean_reads"
	RuleMissingReferenceGenomes = "missing_reference_genomes"
	RuleEmptySampleDir          = "empty_sample_dir"

	RulePipelineDirAtSamplePosition = "pipeline_dir_at_sample_position"
	RuleNestedDirInBams             = "nested_dir_in_bams"
	RuleAlignmentOutsideBamsDir     = "alignment_outside_bams_dir"
	RuleUnexpectedDirAtRefPosition  = "unexpected_dir_at_reference_position"

	RuleReferenceCaseDrift        = "reference_case_drift"
	RuleReferenceNotInGenomesDir  = "reference_not_in_genomes_dir"
	RuleDictMissing               = "dict_missing"
	RuleDictUnreadable            = "dict_unreadable"
	RuleUnrecognisedSeqLabel      = "unrecognised_seq_label"
	RuleAmbiguousSeqLabel         = "ambiguous_seq_label_sanitisation"
	RuleMultipleGvcfStems         = "multiple_gvcf_stems"
	RuleGvcfMissingIndex          = "gvcf_missing_index"
	RuleAlignmentMissingIndex     = "alignment_missing_index"
	RuleAlignmentValidationFailed = "alignment_validation_failed"
	RuleGvcfValidationFailed      = "gvcf_validation_failed"
	RuleJointVcfValidationFailed  = "joint_vcf_validation_failed"
	RuleUnreadableDir             = "unreadable_dir"
	RuleCapacityUnavailable       = "capacity_unavailable"
	RuleDataDirMissing            = "data_dir_missing"
)

// Label spellings a filename can use for a sequence unit.
const (
	// SpellingSanitised is the dotted-to-underscore form GATK and the joint
	// VCF namer write.
	SpellingSanitised = "sanitised"
	// SpellingRaw is the dictionary's own spelling, which DeepVariant keeps.
	SpellingRaw = "raw"
)
