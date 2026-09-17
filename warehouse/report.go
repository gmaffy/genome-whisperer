// Package warehouse inventories the sequencing data directories and audits
// them against the genome-whisperer layout standard.
//
// It is strictly read-only: it stats files and reads reference dictionaries,
// and never writes anywhere except the report path it is given. It also never
// grades anything. A row says what is on disk — sizes in bytes, which files are
// present, which expected units are missing — and leaves "is that good?" to
// whatever consumes the report, because thresholds change far more often than
// filesystem layout does.
package warehouse

import (
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"
)

// SchemaVersion is the report format version. Bump it on any field rename or
// removal so a consumer can refuse a report it does not understand.
const SchemaVersion = 1

// StandardVersion is the version of the directory standard being audited.
const StandardVersion = 1

// Row kinds. A consumer switches on these, so they are a closed set.
const (
	KindReport          = "report"
	KindVolume          = "volume"
	KindCrop            = "crop"
	KindReference       = "reference"
	KindReferenceBatch  = "reference_batch"
	KindSample          = "sample"
	KindSampleReference = "sample_reference"
	KindJunkCandidate   = "junk_candidate"
	KindJointSet        = "joint_set"
	KindNonconformance  = "nonconformance"
	KindSummary         = "summary"
)

// Header is the first line of every report.
type Header struct {
	Kind            string   `json:"kind"`
	SchemaVersion   int      `json:"schema_version"`
	StandardVersion int      `json:"standard_version"`
	Tool            string   `json:"tool"`
	ToolVersion     string   `json:"tool_version"`
	Command         string   `json:"command"`
	RunID           string   `json:"run_id"`
	StartedAt       string   `json:"started_at"`
	Host            string   `json:"host"`
	OS              string   `json:"os"`
	RootsRequested  []string `json:"roots_requested"`
	GenomesDir      string   `json:"genomes_dir"`
	Validate        string   `json:"validate"`
	ScanWorkers     int      `json:"scan_workers"`
	UnitRule        string   `json:"unit_rule"`
}

// Volume is one deduplicated filesystem.
//
// Roots lists every requested path that resolved to it, which is how two
// spellings of a case-insensitive mount collapse into one row instead of
// counting their bytes twice.
type Volume struct {
	Kind     string   `json:"kind"`
	Volume   string   `json:"volume"`
	FSDev    uint64   `json:"fs_dev"`
	Roots    []string `json:"roots"`
	DataDir  string   `json:"data_dir"`
	Spelling string   `json:"data_dir_spelling"`
	Present  bool     `json:"present"`

	// CapacityKnown is false when the platform or the mount would not answer.
	// The byte fields are then meaningless and must not be rendered as zero.
	CapacityKnown  bool  `json:"capacity_known"`
	BytesTotal     int64 `json:"bytes_total"`
	BytesFree      int64 `json:"bytes_free"`
	BytesAvailable int64 `json:"bytes_available"`
	BlockSize      int64 `json:"block_size"`

	Crops []string `json:"crops"`
}

// Crop is one crop directory on one volume.
type Crop struct {
	Kind         string   `json:"kind"`
	Volume       string   `json:"volume"`
	Crop         string   `json:"crop"`
	Spelling     string   `json:"crop_spelling"`
	Path         string   `json:"path"`
	Years        []string `json:"years"`
	SampleCount  int      `json:"sample_count"`
	References   []string `json:"references"`
	JointLayouts []string `json:"joint_layouts"`
	HasMerged    bool     `json:"merged_vcfs_present"`
}

// Dict resolution outcomes. Only DictResolved permits a completeness fraction.
const (
	DictResolved   = "resolved"
	DictMissing    = "dict_missing"
	DictUnknownRef = "reference_unknown"
	DictUnreadable = "dict_unreadable"
)

// Unit is one thing the pipeline calls variants for.
//
// ID is the dictionary sequence ID; Label is the spelling that lands in a
// filename. They differ whenever the ID contains a dot, and both are published
// because GATK writes the label while DeepVariant writes the ID.
type Unit struct {
	ID          string `json:"unit"`
	Label       string `json:"label"`
	UnitKind    string `json:"unit_kind"`
	Length      int    `json:"length"`
	MemberCount int    `json:"member_count"`
}

// Reference is one crop/reference pair and the units it is expected to produce.
type Reference struct {
	Kind       string `json:"kind"`
	Volume     string `json:"volume"`
	Crop       string `json:"crop"`
	Reference  string `json:"reference"`
	Spelling   string `json:"reference_spelling"`
	DictStatus string `json:"dict_status"`
	DictPath   string `json:"dict_path,omitempty"`
	GenomesRef string `json:"genomes_reference,omitempty"`
	Rule       string `json:"rule"`
	SeqCount   int    `json:"seq_count"`

	// Units is absent unless DictStatus is DictResolved. Absent is not the
	// same as empty: an empty expected set that reads as "satisfied" is how a
	// previous inventory of this estate declared complete samples incomplete.
	Units             []Unit `json:"units,omitempty"`
	ExpectedUnitCount int    `json:"expected_unit_count,omitempty"`

	SampleReferenceCount int `json:"sample_reference_count"`
}

// ReferenceBatch carries the member sequence IDs of a batched contig group,
// emitted separately because one pepper reference has ~81,000 of them and
// inlining that makes a single multi-megabyte line.
type ReferenceBatch struct {
	Kind        string   `json:"kind"`
	Crop        string   `json:"crop"`
	Reference   string   `json:"reference"`
	MemberCount int      `json:"member_count"`
	TotalLength int      `json:"total_length"`
	SeqIDs      []string `json:"seqids"`
}

// Read layouts.
const (
	LayoutPaired     = "PE"
	LayoutSingle     = "SE"
	LayoutPairedPlus = "PE+SE"
	LayoutNoReads    = "none"
)

// ReadFile is one FASTQ in a clean_reads directory.
type ReadFile struct {
	Name  string `json:"name"`
	Role  string `json:"role"`
	Bytes int64  `json:"bytes"`
	MTime string `json:"mtime"`
}

// Sample is one crop/year/sample directory: its reads and its shape.
//
// A sample row is emitted even when the directory is empty or malformed, so
// "we have never heard of this sample" and "this sample has nothing in it" stay
// distinguishable.
type Sample struct {
	Kind     string `json:"kind"`
	Volume   string `json:"volume"`
	Crop     string `json:"crop"`
	Year     string `json:"year"`
	Sample   string `json:"sample"`
	Path     string `json:"path"`
	LongRead bool   `json:"long_read"`

	HasCleanReads       bool `json:"has_clean_reads"`
	HasReferenceGenomes bool `json:"has_reference_genomes"`

	ReadLayout string     `json:"read_layout"`
	ReadsFwd   int        `json:"reads_fwd_count"`
	ReadsRev   int        `json:"reads_rev_count"`
	ReadsSE    int        `json:"reads_se_count"`
	ReadsBytes int64      `json:"reads_bytes"`
	ReadFiles  []ReadFile `json:"reads_files"`

	References []string `json:"references"`
}

// JunkCandidate is a directory at sample position whose immediate children
// were read successfully and contained neither sample marker directory.
type JunkCandidate struct {
	Kind       string `json:"kind"`
	Volume     string `json:"volume"`
	Crop       string `json:"crop"`
	Year       string `json:"year"`
	Name       string `json:"name"`
	Path       string `json:"path"`
	Reason     string `json:"reason"`
	Bytes      int64  `json:"bytes"`
	BytesKnown bool   `json:"bytes_known"`
	FileCount  int    `json:"file_count"`
}

// AlignmentFile is one BAM or CRAM, with whatever index sits beside it.
type AlignmentFile struct {
	Name       string `json:"name"`
	Role       string `json:"role"`
	Bytes      int64  `json:"bytes"`
	MTime      string `json:"mtime"`
	Index      string `json:"index,omitempty"`
	IndexBytes int64  `json:"index_bytes"`

	// Validated records whether integrity was checked at all. When false,
	// Valid is meaningless and a consumer must render "not checked" — never
	// "invalid".
	Validated     bool   `json:"validated"`
	Valid         bool   `json:"valid,omitempty"`
	ValidateError string `json:"validate_error,omitempty"`
}

// Inventory statuses for a sample/reference pair.
const (
	// InventoryStandard means gVCFs were found in a standard directory.
	InventoryStandard = "standard"
	// InventoryLegacyOnly means gVCFs exist, but only in a legacy directory,
	// so they are counted and not inventoried. This is the value that keeps
	// "never called" separate from "called into the wrong directory".
	InventoryLegacyOnly = "legacy_only"
	InventoryMixed      = "mixed"
	InventoryAbsent     = "absent"
)

// GvcfSet is one per-sample gVCF directory: one caller's output.
type GvcfSet struct {
	Dir    string `json:"dir"`
	Caller string `json:"caller"`

	// ExpectedUnitCount is omitted when the reference's dict did not resolve,
	// so no consumer can divide by a denominator that was never established.
	ExpectedUnitCount int      `json:"expected_unit_count,omitempty"`
	PresentCount      int      `json:"present_count"`
	Missing           []string `json:"missing"`
	IndexedCount      int      `json:"indexed_count"`
	Unindexed         []string `json:"unindexed"`
	Bytes             int64    `json:"bytes"`

	// Stems are the opaque leading parts of the filenames. They are recorded,
	// never parsed: the sample a gVCF belongs to comes from its directory
	// position, because the stem is frequently not the sample name at all.
	// More than one stem in a directory means mixed provenance.
	Stems []string `json:"stems"`

	// LabelSpellings is "sanitised", "raw", or both — which of the two
	// filename conventions this directory actually uses.
	LabelSpellings []string `json:"label_spellings"`

	SidecarVcfCount    int `json:"sidecar_vcf_count"`
	SidecarReportCount int `json:"sidecar_report_count"`
	UnrecognisedCount  int `json:"unrecognised_count"`
}

// SampleReference is the row a warehouse UI renders per sample per reference.
type SampleReference struct {
	Kind      string `json:"kind"`
	Volume    string `json:"volume"`
	Crop      string `json:"crop"`
	Year      string `json:"year"`
	Sample    string `json:"sample"`
	Reference string `json:"reference"`
	Spelling  string `json:"reference_spelling"`
	Path      string `json:"path"`

	AlignmentDirPresent bool            `json:"alignment_dir_present"`
	Alignments          []AlignmentFile `json:"alignments"`
	AlignmentBytes      int64           `json:"alignment_bytes"`
	AlignmentRoles      []string        `json:"alignment_roles"`

	InventoryStatus string    `json:"inventory_status"`
	GvcfSets        []GvcfSet `json:"gvcf_sets"`

	GvcfLegacyDirCount  int `json:"gvcf_legacy_dir_count"`
	GvcfLegacyFileCount int `json:"gvcf_legacy_file_count"`
}

// Joint VCF directory layouts.
const (
	// JointStandard is <crop>/MERGED_VCFs/<ref>/<caller>_<merger>/.
	JointStandard = "standard"
	// JointLegacyRefDir is the older <crop>/VCFs/<ref>/, with no
	// caller/merger directory and no tag in the filenames.
	JointLegacyRefDir = "legacy_vcfs_ref"
	// JointLegacyRefFirst is <crop>/<ref>/VCFs/<caller>_<merger>/ — standard
	// filenames and a standard tag directory under a non-standard parent.
	JointLegacyRefFirst = "legacy_ref_vcfs"
)

// JointSet is one crop/reference/caller-merger combination of joint VCFs.
type JointSet struct {
	Kind      string `json:"kind"`
	Volume    string `json:"volume"`
	Crop      string `json:"crop"`
	Reference string `json:"reference"`
	Tag       string `json:"tag"`
	Caller    string `json:"caller"`
	Merger    string `json:"merger"`
	Layout    string `json:"layout"`
	Dir       string `json:"dir"`

	// Container is "vcf.gz" or "bcf". GLnexus writes BCF natively, and a
	// *.joint.vcf.gz glob would miss such a cohort entirely.
	Container string `json:"container"`

	ExpectedUnitCount int      `json:"expected_unit_count,omitempty"`
	PresentCount      int      `json:"present_count"`
	Missing           []string `json:"missing"`
	IndexedCount      int      `json:"indexed_count"`
	JointBytes        int64    `json:"joint_bytes"`

	HasAllVcf            bool `json:"has_all_vcf"`
	HasAllIndex          bool `json:"has_all_index"`
	HasHardFiltered      bool `json:"has_hard_filtered"`
	HasHardFilteredIndex bool `json:"has_hard_filtered_index"`
	HasSnpEffVcf         bool `json:"has_snpeff_vcf"`
	HasSnpEffTsv         bool `json:"has_snpeff_tsv"`
	HasEffTsv            bool `json:"has_eff_tsv"`
	HasDescTsv           bool `json:"has_desc_tsv"`
	HasSuperVcfTsv       bool `json:"has_super_vcf_tsv"`

	OtherFileCount int `json:"other_file_count"`
}

// Nonconformance severities.
const (
	SeverityInfo    = "info"
	SeverityWarning = "warning"
	SeverityError   = "error"
)

// Nonconformance is one deviation from the standard.
//
// Detail is a string, never an error: an error-typed field marshals to {} and
// silently loses the message.
type Nonconformance struct {
	Kind      string `json:"kind"`
	Rule      string `json:"rule"`
	Severity  string `json:"severity"`
	Volume    string `json:"volume,omitempty"`
	Crop      string `json:"crop,omitempty"`
	Year      string `json:"year,omitempty"`
	Sample    string `json:"sample,omitempty"`
	Reference string `json:"reference,omitempty"`
	Path      string `json:"path"`
	Found     string `json:"found,omitempty"`
	Expected  string `json:"expected,omitempty"`
	Count     int    `json:"count,omitempty"`
	Detail    string `json:"detail,omitempty"`
}

// Summary is the last line of every report.
//
// A report without this line did not finish. Complete alone is not enough to
// mean "clean", so the counts sit beside it.
type Summary struct {
	Kind             string         `json:"kind"`
	Complete         bool           `json:"complete"`
	Truncated        bool           `json:"truncated"`
	FinishedAt       string         `json:"finished_at"`
	DurationMS       int64          `json:"duration_ms"`
	RowCounts        map[string]int `json:"row_counts"`
	BytesInventoried int64          `json:"bytes_inventoried"`
	FilesInventoried int            `json:"files_inventoried"`
	ErrorCount       int            `json:"error_count"`
	UnreadableDirs   int            `json:"unreadable_dir_count"`
}

// Writer streams rows as newline-delimited JSON.
//
// The report is built in a temporary file beside its destination and renamed
// into place on Close, so a consumer never reads a half-written report: it
// either sees the previous one or a complete new one.
type Writer struct {
	dest string
	tmp  *os.File
	enc  *json.Encoder
	out  io.Writer

	counts     map[string]int
	bytes      int64
	files      int
	errors     int
	unreadable int
	started    time.Time
}

// NewWriter creates a report at dest.
func NewWriter(dest string) (*Writer, error) {
	abs, err := filepath.Abs(dest)
	if err != nil {
		return nil, fmt.Errorf("resolving report path %s: %w", dest, err)
	}
	if err := os.MkdirAll(filepath.Dir(abs), 0o755); err != nil {
		return nil, fmt.Errorf("creating report directory: %w", err)
	}

	// Same directory as the destination, so the rename is atomic rather than
	// a cross-device copy.
	tmp, err := os.CreateTemp(filepath.Dir(abs), filepath.Base(abs)+".partial-*")
	if err != nil {
		return nil, fmt.Errorf("creating temporary report: %w", err)
	}

	w := &Writer{
		dest:    abs,
		tmp:     tmp,
		out:     tmp,
		counts:  map[string]int{},
		started: time.Now(),
	}
	w.enc = json.NewEncoder(tmp)
	return w, nil
}

// Dest is the final report path.
func (w *Writer) Dest() string { return w.dest }

// Write emits one row. kind is counted for the summary.
func (w *Writer) Write(kind string, row any) error {
	w.counts[kind]++
	// json.Encoder.Encode appends a newline, which is the NDJSON framing.
	if err := w.enc.Encode(row); err != nil {
		return fmt.Errorf("writing %s row: %w", kind, err)
	}
	return nil
}

// CountFiles records inventoried files and bytes for the summary.
func (w *Writer) CountFiles(n int, bytes int64) {
	w.files += n
	w.bytes += bytes
}

// CountError records a recoverable problem.
func (w *Writer) CountError() { w.errors++ }

// CountUnreadableDir records a directory that could not be listed.
func (w *Writer) CountUnreadableDir() {
	w.unreadable++
	w.errors++
}

// Close writes the summary, flushes, and renames the report into place.
func (w *Writer) Close(complete bool) error {
	summary := Summary{
		Kind:             KindSummary,
		Complete:         complete,
		Truncated:        !complete,
		FinishedAt:       time.Now().UTC().Format(time.RFC3339),
		DurationMS:       time.Since(w.started).Milliseconds(),
		BytesInventoried: w.bytes,
		FilesInventoried: w.files,
		ErrorCount:       w.errors,
		UnreadableDirs:   w.unreadable,
	}
	// Snapshot the counts before adding the summary's own row.
	summary.RowCounts = make(map[string]int, len(w.counts)+1)
	for k, v := range w.counts {
		summary.RowCounts[k] = v
	}
	summary.RowCounts[KindSummary] = 1

	if err := w.enc.Encode(summary); err != nil {
		w.tmp.Close()
		os.Remove(w.tmp.Name())
		return fmt.Errorf("writing summary row: %w", err)
	}

	// Sync before rename: a crash after the rename must not leave a report
	// whose tail never reached the disk.
	if err := w.tmp.Sync(); err != nil {
		w.tmp.Close()
		os.Remove(w.tmp.Name())
		return fmt.Errorf("flushing report: %w", err)
	}
	if err := w.tmp.Close(); err != nil {
		os.Remove(w.tmp.Name())
		return fmt.Errorf("closing report: %w", err)
	}
	if err := os.Rename(w.tmp.Name(), w.dest); err != nil {
		os.Remove(w.tmp.Name())
		return fmt.Errorf("moving report into place: %w", err)
	}
	return nil
}

// Abandon discards a partial report without leaving it behind.
func (w *Writer) Abandon() {
	w.tmp.Close()
	os.Remove(w.tmp.Name())
}
