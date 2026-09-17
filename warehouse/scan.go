package warehouse

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/gmaffy/genome-whisperer/alignmentdir"
	"github.com/gmaffy/genome-whisperer/utils"
	"github.com/gmaffy/genome-whisperer/variants"
	"golang.org/x/sync/errgroup"
)

// Validation scopes.
const (
	ValidateNone       = "none"
	ValidateAlignments = "alignments"
	ValidateAll        = "all"
)

// DefaultWorkers is how many directories are read at once.
//
// This is deliberately not derived from the core count. The cost of a scan is
// network round-trips to the storage array, not computation: measured on one
// crop, reading 345 directories took 44 s serially and 7.3 s at 16 workers, on
// a machine whose core count is irrelevant to either figure.
const DefaultWorkers = 16

// MaxWorkers caps the pool so a typo cannot open thousands of concurrent reads
// against a network mount.
const MaxWorkers = 64

// Options configures a scan.
type Options struct {
	// Roots are mounts or data directories to scan.
	Roots []string
	// GenomesDir holds the reference assemblies.
	GenomesDir string
	// ReportPath is where the NDJSON report is written.
	ReportPath string
	// Validate is one of the Validate* constants.
	Validate string
	// Workers is the directory read concurrency; zero means DefaultWorkers.
	Workers int
	// EmitBatchSeqIDs includes every batched contig sequence ID in the
	// report. One pepper reference has about 81,000 of them.
	EmitBatchSeqIDs bool
	Verbose         bool
	Quick           bool
}

func (o Options) workers() int {
	n := o.Workers
	if n <= 0 {
		n = DefaultWorkers
	}
	if n > MaxWorkers {
		n = MaxWorkers
	}
	if n > runtime.GOMAXPROCS(0)*8 {
		// Still generous: these goroutines are almost always blocked on I/O.
		n = runtime.GOMAXPROCS(0) * 8
	}
	if n < 1 {
		n = 1
	}
	return n
}

// discovered is one sample directory found during discovery.
type discovered struct {
	volume, crop, year, sample, path string
	hasCleanReads                    bool
	hasReferenceGenomes              bool
	references                       []string
}

// collector gathers rows and problems from concurrent workers.
type collector struct {
	mu               sync.Mutex
	problems         []Nonconformance
	samples          []Sample
	sampleReferences []SampleReference
	junk             []JunkCandidate
	events           *EventWriter
	inspectDone      int
	files            int
	bytes            int64
}

func (c *collector) problem(p ...Nonconformance) {
	if len(p) == 0 {
		return
	}
	c.mu.Lock()
	c.problems = append(c.problems, p...)
	c.mu.Unlock()
	for _, problem := range p {
		c.events.Emit(ScanEvent{Event: KindNonconformance,
			Volume: problem.Volume, Path: problem.Path, Payload: problem})
	}
}

func (c *collector) counted(files int, bytes int64) {
	c.mu.Lock()
	c.files += files
	c.bytes += bytes
	totalFiles, totalBytes := c.files, c.bytes
	c.mu.Unlock()
	c.events.Emit(ScanEvent{Event: "inventory_progress",
		Payload: map[string]any{"files_inventoried": totalFiles, "bytes_inventoried": totalBytes}})
}

// Scan inventories every root and writes the report.
//
// It is read-only apart from the report itself, and it does not stop for bad
// data: an unreadable directory or an unresolvable reference becomes a row and
// the walk continues. Only a setup failure — an unusable report path, no
// readable root — returns an error, because those mean the report would be
// worthless rather than incomplete.
func Scan(opts Options) error {
	started := time.Now()
	events := NewEventWriter(os.Stdout)
	stopHeartbeats := make(chan struct{})
	go events.Heartbeats(stopHeartbeats)
	defer close(stopHeartbeats)
	events.Emit(ScanEvent{Event: "phase", Phase: "setup", Message: "Resolving warehouse roots"})

	if opts.ReportPath == "" {
		return fmt.Errorf("no report path given")
	}
	if opts.Validate == "" {
		opts.Validate = ValidateNone
	}
	switch opts.Validate {
	case ValidateNone, ValidateAlignments, ValidateAll:
	default:
		return fmt.Errorf("validate must be %q, %q or %q, got %q",
			ValidateNone, ValidateAlignments, ValidateAll, opts.Validate)
	}

	volumes, volProblems := resolveVolumes(opts.Roots)
	if len(volumes) == 0 {
		return fmt.Errorf("no readable data directory among %v", opts.Roots)
	}

	// A report written inside a scanned tree becomes part of the next scan's
	// input, so two consecutive runs never agree. This estate already has one
	// hand-built inventory sitting in a data directory being counted as a
	// stray file.
	if err := reportPathIsOutside(opts.ReportPath, volumes); err != nil {
		return err
	}

	w, err := NewWriter(opts.ReportPath)
	if err != nil {
		return err
	}

	host, _ := os.Hostname()
	header := Header{
		Kind: KindReport, SchemaVersion: SchemaVersion, StandardVersion: StandardVersion,
		Tool: "genome-whisperer", ToolVersion: "1.0.0", Command: "ScanWarehouse",
		RunID:     started.UTC().Format("20060102T150405Z"),
		StartedAt: started.UTC().Format(time.RFC3339),
		Host:      host, OS: runtime.GOOS,
		RootsRequested: opts.Roots, GenomesDir: opts.GenomesDir,
		Validate: opts.Validate, ScanWorkers: opts.workers(),
		UnitRule: variants.UnitRule,
	}
	if err := w.Write(KindReport, header); err != nil {
		w.Abandon()
		return err
	}

	events.Emit(ScanEvent{Event: "phase", Phase: "discovery", Message: "Discovering warehouse directories"})

	col := &collector{problems: volProblems, events: events}
	genomes := newGenomeIndex(opts.GenomesDir)

	// Phase 1: discover the tree. Cheap and serial down to sample names, so
	// that nonconformance rows come out in a stable order, then parallel for
	// the per-sample reads.
	cropRows, samples := discoverTree(volumes, col, opts)
	events.Emit(ScanEvent{Event: "phase", Phase: "references", Message: "Resolving reference genomes", Total: len(samples)})

	// Phase 2: resolve every reference once. Must precede inspection: without
	// the expected unit list there is no denominator to compare against.
	refIndex := resolveReferences(genomes, samples, col)
	events.Emit(ScanEvent{Event: "phase", Phase: "inspection", Message: "Inspecting samples", Total: len(samples)})

	// Phase 3: inspect. This is the bulk of the I/O and where the pool pays.
	inspectSamples(samples, refIndex, col, opts)

	// Phase 4: the crop-wide joint VCFs.
	events.Emit(ScanEvent{Event: "phase", Phase: "joint", Message: "Inspecting joint variant sets"})
	jointRows := inspectJointSets(volumes, refIndex, col, opts)

	// Phase 5: write, in a fixed order so two scans of an unchanged tree
	// produce identical bytes and a diff shows only real change.
	events.Emit(ScanEvent{Event: "phase", Phase: "report", Message: "Writing final inventory report"})
	if err := writeAll(w, volumes, cropRows, refIndex, col, jointRows, opts); err != nil {
		w.Abandon()
		return err
	}

	w.CountFiles(col.files, col.bytes)
	for _, p := range col.problems {
		if p.Rule == RuleUnreadableDir {
			w.CountUnreadableDir()
		} else if p.Severity == SeverityError {
			w.CountError()
		}
	}

	if err := w.Close(true); err != nil {
		return err
	}
	events.Emit(ScanEvent{Event: "complete", Phase: "report", Message: "Scanner report complete"})
	return nil
}

// reportPathIsOutside refuses a report path inside any scanned tree.
func reportPathIsOutside(report string, volumes []resolvedVolume) error {
	abs, err := filepath.Abs(report)
	if err != nil {
		return fmt.Errorf("resolving report path: %w", err)
	}
	for _, v := range volumes {
		rel, rErr := filepath.Rel(v.DataDir, abs)
		if rErr != nil {
			continue
		}
		if rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
			return fmt.Errorf("report path %s is inside the scanned directory %s; "+
				"write it elsewhere or the next scan will inventory this report", abs, v.DataDir)
		}
	}
	return nil
}

// discoverTree walks volumes down to sample directories.
func discoverTree(volumes []resolvedVolume, col *collector, opts Options) ([]Crop, []*discovered) {
	var cropRows []Crop
	var samples []*discovered

	for _, v := range volumes {
		dataListing := listDir(v.DataDir)
		if dataListing.Err != "" {
			col.problem(Nonconformance{
				Kind: KindNonconformance, Rule: RuleUnreadableDir,
				Severity: SeverityError, Volume: v.Name, Path: v.DataDir,
				Detail: dataListing.Err,
			})
			continue
		}

		// Loose files beside the crops. This estate has three.
		for _, f := range dataListing.Files {
			col.problem(Nonconformance{
				Kind: KindNonconformance, Rule: RuleStrayFileAtCropPosition,
				Severity: SeverityInfo, Volume: v.Name,
				Path: filepath.Join(v.DataDir, f.Name), Found: f.Name,
			})
		}

		crops := append([]string{}, dataListing.SubDirs...)
		sort.Strings(crops)

		for _, crop := range crops {
			if IsPruned(crop) {
				continue
			}
			row, cropSamples := discoverCrop(v, crop, col)
			cropRows = append(cropRows, row)
			samples = append(samples, cropSamples...)
		}
	}

	// Fill in each sample's own contents in parallel: the sample directory and
	// its reference_genomes listing.
	g := new(errgroup.Group)
	g.SetLimit(opts.workers())
	for _, s := range samples {
		s := s
		g.Go(func() error {
			describeSample(s, col)
			return nil
		})
	}
	_ = g.Wait()

	confirmed := samples[:0]
	for _, s := range samples {
		if s.hasCleanReads || s.hasReferenceGenomes {
			confirmed = append(confirmed, s)
		}
	}
	samples = confirmed

	// Now that references are known, complete the crop rows.
	byCrop := map[string]*Crop{}
	for i := range cropRows {
		byCrop[cropRows[i].Volume+"/"+cropRows[i].Crop] = &cropRows[i]
	}
	for _, s := range samples {
		if row, ok := byCrop[s.volume+"/"+s.crop]; ok {
			row.References = append(row.References, s.references...)
		}
	}
	for i := range cropRows {
		cropRows[i].References = uniqueSorted(cropRows[i].References)
		cropRows[i].SampleCount = 0
	}
	for _, s := range samples {
		if row, ok := byCrop[s.volume+"/"+s.crop]; ok {
			row.SampleCount++
		}
	}
	col.events.Emit(ScanEvent{Event: "discovery_complete", Phase: "discovery",
		Message: "Candidate discovery complete", Total: len(samples),
		Payload: map[string]int{"samples": len(samples), "junk_candidates": len(col.junk)}})

	return cropRows, samples
}

// discoverCrop classifies everything at the year position and collects samples.
func discoverCrop(v resolvedVolume, crop string, col *collector) (Crop, []*discovered) {
	cropPath := filepath.Join(v.DataDir, crop)
	row := Crop{
		Kind: KindCrop, Volume: v.Name, Crop: crop, Spelling: crop, Path: cropPath,
	}

	listing := listDir(cropPath)
	if listing.Err != "" {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnreadableDir,
			Severity: SeverityError, Volume: v.Name, Crop: crop, Path: cropPath,
			Detail: listing.Err,
		})
		return row, nil
	}

	for _, f := range listing.Files {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleStrayFileAtYearPosition,
			Severity: SeverityInfo, Volume: v.Name, Crop: crop,
			Path: filepath.Join(cropPath, f.Name), Found: f.Name,
		})
	}

	var samples []*discovered
	entries := append([]string{}, listing.SubDirs...)
	sort.Strings(entries)

	for _, name := range entries {
		if IsPruned(name) {
			continue
		}
		path := filepath.Join(cropPath, name)

		if IsYear(name) {
			row.Years = append(row.Years, name)
			samples = append(samples, discoverYear(v, crop, name, path, col)...)
			continue
		}

		switch name {
		case DirMergedVcfs:
			row.HasMerged = true
			row.JointLayouts = append(row.JointLayouts, JointStandard)
			continue
		case DirLegacyVcfs:
			// Counted, never inventoried: the layout predates the
			// caller/merger directory.
			row.JointLayouts = append(row.JointLayouts, JointLegacyRefDir)
			col.problem(Nonconformance{
				Kind: KindNonconformance, Rule: RuleLegacyJointRoot,
				Severity: SeverityWarning, Volume: v.Name, Crop: crop, Path: path,
				Found:    DirLegacyVcfs,
				Expected: DirMergedVcfs + "/<reference>/<caller>_<merger>",
				Detail:   "joint VCFs here are reported but not inventoried",
			})
			continue
		}

		// A non-year directory. One more read tells us what it actually is,
		// which matters because a real sample sits at this position on this
		// estate and must not be discarded as clutter.
		col.problem(classifyYearPositionDir(v, crop, name, path)...)
		if sample := sampleAtYearPosition(v, crop, name, path); sample != nil {
			samples = append(samples, sample)
		}
		if isRefFirstJointRoot(path) {
			row.JointLayouts = append(row.JointLayouts, JointLegacyRefFirst)
		}
	}

	sort.Strings(row.Years)
	row.JointLayouts = uniqueSorted(row.JointLayouts)
	row.SampleCount = len(samples)
	return row, samples
}

// discoverYear collects the sample directories of one cohort.
func discoverYear(v resolvedVolume, crop, year, yearPath string, col *collector) []*discovered {
	col.events.Emit(ScanEvent{Event: "operation_start", Phase: "discovery", Volume: v.Name, Path: yearPath})
	listing := listDir(yearPath)
	if listing.Err != "" {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnreadableDir,
			Severity: SeverityError, Volume: v.Name, Crop: crop, Year: year,
			Path: yearPath, Detail: listing.Err,
		})
		return nil
	}

	for _, f := range listing.Files {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleStrayFileAtSamplePosition,
			Severity: SeverityInfo, Volume: v.Name, Crop: crop, Year: year,
			Path: filepath.Join(yearPath, f.Name), Found: f.Name,
		})
	}

	names := append([]string{}, listing.SubDirs...)
	sort.Strings(names)

	out := make([]*discovered, 0, len(names))
	for _, name := range names {
		if IsPruned(name) {
			continue
		}
		out = append(out, &discovered{
			volume: v.Name, crop: crop, year: year, sample: name,
			path: filepath.Join(yearPath, name),
		})
	}
	return out
}

// classifyYearPositionDir explains a directory that is not a year.
func classifyYearPositionDir(v resolvedVolume, crop, name, path string) []Nonconformance {
	base := Nonconformance{
		Kind: KindNonconformance, Volume: v.Name, Crop: crop, Path: path, Found: name,
	}

	listing := listDir(path)
	if listing.Err != "" {
		base.Rule, base.Severity, base.Detail = RuleUnreadableDir, SeverityError, listing.Err
		return []Nonconformance{base}
	}

	children := map[string]bool{}
	for _, d := range listing.SubDirs {
		children[d] = true
	}

	switch {
	case children[DirCleanReads] || children[DirReferenceGenomes]:
		base.Rule = RuleSampleAtYearPosition
		base.Severity = SeverityWarning
		base.Expected = "<crop>/<year>/" + name
		base.Detail = "a real sample directly under the crop; inventoried with an empty year"
	case children[DirLegacyVcfs]:
		base.Rule = RuleLegacyJointRootRefFirst
		base.Severity = SeverityWarning
		base.Expected = DirMergedVcfs + "/" + name + "/<caller>_<merger>"
		base.Detail = "reference-first joint VCF layout; reported but not inventoried"
	default:
		base.Rule = RuleForeignDirAtYearPosition
		base.Severity = SeverityInfo
		base.Detail = "not a cohort, a sample or a joint VCF root; not descended into"
	}
	return []Nonconformance{base}
}

// sampleAtYearPosition returns a sample found directly under a crop.
func sampleAtYearPosition(v resolvedVolume, crop, name, path string) *discovered {
	listing := listDir(path)
	for _, d := range listing.SubDirs {
		if d == DirCleanReads || d == DirReferenceGenomes {
			// Year is empty rather than invented: the cohort is genuinely
			// unknown, and a made-up one would be indistinguishable from a
			// real one downstream.
			return &discovered{
				volume: v.Name, crop: crop, year: "", sample: name, path: path,
			}
		}
	}
	return nil
}

// isRefFirstJointRoot reports whether path is <ref>/VCFs/<tag>.
func isRefFirstJointRoot(path string) bool {
	listing := listDir(filepath.Join(path, DirLegacyVcfs))
	return listing.Exists && listing.Err == ""
}

// describeSample reads one sample directory and its reference list.
func describeSample(s *discovered, col *collector) {
	col.events.Emit(ScanEvent{Event: "operation_start", Phase: "discovery", Volume: s.volume, Path: s.path})
	listing := listDir(s.path)
	if listing.Err != "" {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnreadableDir,
			Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
			Sample: s.sample, Path: s.path, Detail: listing.Err,
		})
		return
	}

	for _, d := range listing.SubDirs {
		switch d {
		case DirCleanReads:
			s.hasCleanReads = true
		case DirReferenceGenomes:
			s.hasReferenceGenomes = true
		case DirBams, DirLegacyGvcfs, "gatk_gvcfs", "dv_gvcfs":
			// Pipeline output directly under the sample, skipping the
			// reference_genomes/<ref> levels entirely.
			col.problem(Nonconformance{
				Kind: KindNonconformance, Rule: RulePipelineDirAtSamplePosition,
				Severity: SeverityWarning, Volume: s.volume, Crop: s.crop,
				Year: s.year, Sample: s.sample,
				Path: filepath.Join(s.path, d), Found: d,
				Expected: DirReferenceGenomes + "/<reference>/" + d,
				Detail:   "not inventoried; no reference can be attributed to it",
			})
		}
	}

	if !s.hasCleanReads && !s.hasReferenceGenomes {
		regularFiles := 0
		bytesKnown := len(listing.SubDirs) == 0
		for _, file := range listing.Files {
			if file.IsRegular() {
				regularFiles++
			} else {
				bytesKnown = false
			}
		}
		junk := JunkCandidate{
			Kind: KindJunkCandidate, Volume: s.volume, Crop: s.crop, Year: s.year,
			Name: s.sample, Path: s.path, Reason: "no_sample_markers",
			Bytes: listing.TotalBytes(), BytesKnown: bytesKnown,
			FileCount: regularFiles,
		}
		col.mu.Lock()
		col.junk = append(col.junk, junk)
		col.mu.Unlock()
		col.events.Emit(ScanEvent{Event: KindJunkCandidate, Phase: "discovery",
			Volume: s.volume, Path: s.path, Payload: junk})
		return
	}

	if !s.hasCleanReads {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleMissingCleanReads,
			Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
			Sample: s.sample, Path: s.path, Expected: DirCleanReads,
		})
	}
	if !s.hasReferenceGenomes {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleMissingReferenceGenomes,
			Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
			Sample: s.sample, Path: s.path, Expected: DirReferenceGenomes,
		})
		if !s.hasCleanReads {
			col.problem(Nonconformance{
				Kind: KindNonconformance, Rule: RuleEmptySampleDir,
				Severity: SeverityWarning, Volume: s.volume, Crop: s.crop,
				Year: s.year, Sample: s.sample, Path: s.path,
			})
		}
		return
	}

	refs := listDir(filepath.Join(s.path, DirReferenceGenomes))
	if refs.Err != "" {
		col.problem(Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnreadableDir,
			Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
			Sample: s.sample, Path: refs.Path, Detail: refs.Err,
		})
		return
	}
	for _, name := range refs.SubDirs {
		if IsPruned(name) {
			continue
		}
		s.references = append(s.references, name)
	}
	sort.Strings(s.references)
}

// refKey identifies a crop/reference pair.
type refKey struct{ volume, crop, reference string }

// refInfo is a resolved reference: its row, its expected units and a matcher.
type refInfo struct {
	row     Reference
	units   []Unit
	batch   *ReferenceBatch
	matcher *unitMatcher
	fasta   string
}

// resolveReferences resolves every crop/reference pair seen during discovery.
func resolveReferences(genomes *genomeIndex, samples []*discovered, col *collector) map[refKey]*refInfo {
	index := map[refKey]*refInfo{}
	counts := map[refKey]int{}

	for _, s := range samples {
		for _, ref := range s.references {
			key := refKey{s.volume, s.crop, ref}
			counts[key]++
			if _, seen := index[key]; seen {
				continue
			}
			row, units, batch, fasta, problems := genomes.referenceRow(s.volume, s.crop, ref)
			col.problem(problems...)

			info := &refInfo{row: row, units: units, batch: batch, fasta: fasta}
			if len(units) > 0 {
				info.matcher = newUnitMatcher(units)
			}
			index[key] = info
		}
	}

	for key, n := range counts {
		index[key].row.SampleReferenceCount = n
	}
	return index
}

// inspectSamples visits each discovered sample directory in parallel and collects
// read layouts, alignment files, and gVCF sets.
func inspectSamples(samples []*discovered, refIndex map[refKey]*refInfo, col *collector, opts Options) {
	g := new(errgroup.Group)
	g.SetLimit(opts.workers())

	for _, s := range samples {
		s := s
		g.Go(func() error {
			col.events.Emit(ScanEvent{Event: "operation_start", Phase: "inspection",
				Volume: s.volume, Path: s.path, Total: len(samples)})
			sampleRow, sampleRefs, problems, files, bytes := inspectOneSample(s, refIndex, opts)
			col.mu.Lock()
			col.samples = append(col.samples, sampleRow)
			col.sampleReferences = append(col.sampleReferences, sampleRefs...)
			col.files += files
			col.bytes += bytes
			totalFiles, totalBytes := col.files, col.bytes
			col.inspectDone++
			done := col.inspectDone
			col.mu.Unlock()
			col.problem(problems...)
			col.events.Emit(ScanEvent{Event: "inventory_progress", Phase: "inspection",
				Payload: map[string]any{"files_inventoried": totalFiles, "bytes_inventoried": totalBytes}})
			col.events.Emit(ScanEvent{Event: KindSample, Phase: "inspection",
				Volume: s.volume, Path: s.path, Completed: done, Total: len(samples), Payload: sampleRow})
			for _, row := range sampleRefs {
				col.events.Emit(ScanEvent{Event: KindSampleReference, Phase: "inspection",
					Volume: s.volume, Path: row.Path, Completed: done, Total: len(samples), Payload: row})
			}
			return nil
		})
	}
	_ = g.Wait()
}

func inspectOneSample(s *discovered, refIndex map[refKey]*refInfo, opts Options) (
	Sample, []SampleReference, []Nonconformance, int, int64,
) {
	var problems []Nonconformance
	var totalFiles int
	var totalBytes int64

	isLongRead := alignmentdir.IsLongReadSample(s.sample)
	sampleRow := Sample{
		Kind:                KindSample,
		Volume:              s.volume,
		Crop:                s.crop,
		Year:                s.year,
		Sample:              s.sample,
		Path:                s.path,
		LongRead:            isLongRead,
		HasCleanReads:       s.hasCleanReads,
		HasReferenceGenomes: s.hasReferenceGenomes,
		ReadLayout:          LayoutNoReads,
		References:          append([]string{}, s.references...),
	}
	sortStrings(sampleRow.References)

	if s.hasCleanReads {
		readsDir := filepath.Join(s.path, DirCleanReads)
		listing := listDir(readsDir)
		if listing.Err != "" {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleUnreadableDir,
				Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
				Sample: s.sample, Path: readsDir, Detail: listing.Err,
			})
		} else {
			for _, f := range listing.Files {
				if !f.IsRegular() {
					continue
				}
				if !IsFastq(f.Name) {
					problems = append(problems, Nonconformance{
						Kind: KindNonconformance, Rule: RuleStrayFileInCleanReads,
						Severity: SeverityInfo, Volume: s.volume, Crop: s.crop, Year: s.year,
						Sample: s.sample, Path: filepath.Join(readsDir, f.Name), Found: f.Name,
					})
					continue
				}

				totalFiles++
				totalBytes += f.Bytes
				sampleRow.ReadsBytes += f.Bytes

				role := alignmentdir.ClassifyRead(f.Name)
				roleStr := "se"
				switch role {
				case alignmentdir.ReadRoleForward:
					sampleRow.ReadsFwd++
					roleStr = "fwd"
				case alignmentdir.ReadRoleReverse:
					sampleRow.ReadsRev++
					roleStr = "rev"
				default:
					sampleRow.ReadsSE++
					roleStr = "se"
				}

				sampleRow.ReadFiles = append(sampleRow.ReadFiles, ReadFile{
					Name:  f.Name,
					Role:  roleStr,
					Bytes: f.Bytes,
					MTime: f.MTime(),
				})
			}

			sort.Slice(sampleRow.ReadFiles, func(i, j int) bool {
				return sampleRow.ReadFiles[i].Name < sampleRow.ReadFiles[j].Name
			})

			switch {
			case sampleRow.ReadsFwd > 0 && sampleRow.ReadsRev > 0 && sampleRow.ReadsSE > 0:
				sampleRow.ReadLayout = LayoutPairedPlus
			case sampleRow.ReadsFwd > 0 && sampleRow.ReadsRev > 0:
				sampleRow.ReadLayout = LayoutPaired
			case sampleRow.ReadsSE > 0:
				sampleRow.ReadLayout = LayoutSingle
			default:
				sampleRow.ReadLayout = LayoutNoReads
			}
		}
	}

	var sampleRefs []SampleReference
	for _, ref := range s.references {
		refPath := filepath.Join(s.path, DirReferenceGenomes, ref)
		info := refIndex[refKey{s.volume, s.crop, ref}]

		srRow, srProblems, srFiles, srBytes := inspectSampleReference(s, ref, refPath, info, opts)
		sampleRefs = append(sampleRefs, srRow)
		problems = append(problems, srProblems...)
		totalFiles += srFiles
		totalBytes += srBytes
	}

	return sampleRow, sampleRefs, problems, totalFiles, totalBytes
}

func inspectSampleReference(s *discovered, ref, refPath string, info *refInfo, opts Options) (
	SampleReference, []Nonconformance, int, int64,
) {
	var problems []Nonconformance
	var refFiles int
	var refBytes int64

	spelling := ref
	var matcher *unitMatcher
	var expectedUnits []Unit
	if info != nil {
		spelling = info.row.Spelling
		matcher = info.matcher
		expectedUnits = info.units
	}

	sr := SampleReference{
		Kind:            KindSampleReference,
		Volume:          s.volume,
		Crop:            s.crop,
		Year:            s.year,
		Sample:          s.sample,
		Reference:       ref,
		Spelling:        spelling,
		Path:            refPath,
		InventoryStatus: InventoryAbsent,
	}

	refListing := listDir(refPath)
	if refListing.Err != "" {
		problems = append(problems, Nonconformance{
			Kind: KindNonconformance, Rule: RuleUnreadableDir,
			Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
			Sample: s.sample, Reference: ref, Path: refPath, Detail: refListing.Err,
		})
		return sr, problems, 0, 0
	}

	for _, f := range refListing.Files {
		if f.IsRegular() {
			role := alignmentdir.AlignmentRole(f.Name)
			if role != alignmentdir.RoleIgnore && role != alignmentdir.RoleOther {
				problems = append(problems, Nonconformance{
					Kind: KindNonconformance, Rule: RuleAlignmentOutsideBamsDir,
					Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
					Sample: s.sample, Reference: ref,
					Path:     filepath.Join(refPath, f.Name),
					Found:    f.Name,
					Expected: DirBams + "/" + f.Name,
					Detail:   "alignment file outside bams/ directory",
				})
			}
		}
	}

	bamsPath := filepath.Join(refPath, DirBams)
	bamsListing := listDir(bamsPath)
	if bamsListing.Exists {
		sr.AlignmentDirPresent = true
		if bamsListing.Err != "" {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleUnreadableDir,
				Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
				Sample: s.sample, Reference: ref, Path: bamsPath, Detail: bamsListing.Err,
			})
		} else {
			for _, sub := range bamsListing.SubDirs {
				problems = append(problems, Nonconformance{
					Kind: KindNonconformance, Rule: RuleNestedDirInBams,
					Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
					Sample: s.sample, Reference: ref,
					Path:   filepath.Join(bamsPath, sub),
					Found:  sub,
					Detail: "unexpected nested directory in bams/",
				})
			}

			byName := bamsListing.names()
			for _, f := range bamsListing.Files {
				if !f.IsRegular() {
					continue
				}

				role := alignmentdir.AlignmentRole(f.Name)
				if role == alignmentdir.RoleIgnore {
					continue
				}
				lowerName := strings.ToLower(f.Name)
				isAlignment := strings.HasSuffix(lowerName, ".bam") || strings.HasSuffix(lowerName, ".cram")

				refFiles++
				refBytes += f.Bytes
				sr.AlignmentBytes += f.Bytes
				sr.AlignmentRoles = append(sr.AlignmentRoles, role)

				idxName, idxBytes, hasIdx := alignmentIndexFor(byName, f.Name)
				if hasIdx {
					refFiles++
					refBytes += idxBytes
				} else if isAlignment {
					problems = append(problems, Nonconformance{
						Kind: KindNonconformance, Rule: RuleAlignmentMissingIndex,
						Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
						Sample: s.sample, Reference: ref,
						Path:   filepath.Join(bamsPath, f.Name),
						Detail: "alignment file missing index (.bai / .crai / .csi)",
					})
				}

				af := AlignmentFile{
					Name:       f.Name,
					Role:       role,
					Bytes:      f.Bytes,
					MTime:      f.MTime(),
					Index:      idxName,
					IndexBytes: idxBytes,
					Validated:  false,
				}

				if isAlignment && (opts.Validate == ValidateAlignments || opts.Validate == ValidateAll) && info != nil && info.fasta != "" {
					af.Validated = true
					valErr := utils.ValidateBam(filepath.Join(bamsPath, f.Name), info.fasta, false, opts.Quick)
					af.Valid = (valErr == nil)
					if valErr != nil {
						af.ValidateError = valErr.Error()
						problems = append(problems, Nonconformance{
							Kind: KindNonconformance, Rule: RuleAlignmentValidationFailed,
							Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
							Sample: s.sample, Reference: ref, Path: filepath.Join(bamsPath, f.Name),
							Detail: valErr.Error(),
						})
					}
				}

				sr.Alignments = append(sr.Alignments, af)
			}
			sort.Slice(sr.Alignments, func(i, j int) bool {
				return sr.Alignments[i].Name < sr.Alignments[j].Name
			})
			sr.AlignmentRoles = uniqueSorted(sr.AlignmentRoles)
		}
	}

	for _, sub := range refListing.SubDirs {
		subPath := filepath.Join(refPath, sub)
		subListing := listDir(subPath)
		if subListing.Err != "" {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleUnreadableDir,
				Severity: SeverityError, Volume: s.volume, Crop: s.crop, Year: s.year,
				Sample: s.sample, Reference: ref, Path: subPath, Detail: subListing.Err,
			})
			continue
		}

		if caller, isStandard := GvcfDirCaller(sub); isStandard {
			set, gvcfProblems := inspectGvcfDir(subListing, sub, caller, matcher, expectedUnits,
				opts.Validate == ValidateAll, opts.Quick)
			sr.GvcfSets = append(sr.GvcfSets, set)
			problems = append(problems, gvcfProblems...)
			refFiles += len(subListing.Files)
			refBytes += subListing.TotalBytes()
			continue
		}

		if sub == DirLegacyGvcfs {
			sr.GvcfLegacyDirCount++
			legFiles, legBytes := countLegacyGvcfs(subListing)
			sr.GvcfLegacyFileCount += legFiles
			refFiles += len(subListing.Files)
			refBytes += legBytes
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleLegacyGvcfDir,
				Severity: SeverityWarning, Volume: s.volume, Crop: s.crop, Year: s.year,
				Sample: s.sample, Reference: ref, Path: subPath,
				Found:    DirLegacyGvcfs,
				Expected: "gatk_gvcfs or dv_gvcfs",
				Count:    legFiles,
				Detail:   "legacy gvcfs/ directory; reported but not inventoried",
			})
			continue
		}

		if sub != DirBams {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleUnexpectedDirAtRefPosition,
				Severity: SeverityInfo, Volume: s.volume, Crop: s.crop, Year: s.year,
				Sample: s.sample, Reference: ref, Path: subPath, Found: sub,
				Detail: "unrecognised directory under reference",
			})
		}
	}

	sort.Slice(sr.GvcfSets, func(i, j int) bool {
		return sr.GvcfSets[i].Dir < sr.GvcfSets[j].Dir
	})

	hasStandard := len(sr.GvcfSets) > 0
	hasLegacy := sr.GvcfLegacyDirCount > 0
	switch {
	case hasStandard && !hasLegacy:
		sr.InventoryStatus = InventoryStandard
	case hasStandard && hasLegacy:
		sr.InventoryStatus = InventoryMixed
	case !hasStandard && hasLegacy:
		sr.InventoryStatus = InventoryLegacyOnly
	default:
		sr.InventoryStatus = InventoryAbsent
	}

	return sr, problems, refFiles, refBytes
}

func inspectJointSets(volumes []resolvedVolume, refIndex map[refKey]*refInfo, col *collector, opts Options) []JointSet {
	var jointRows []JointSet

	for _, v := range volumes {
		dataListing := listDir(v.DataDir)
		if dataListing.Err != "" {
			continue
		}

		for _, crop := range dataListing.SubDirs {
			if IsPruned(crop) {
				continue
			}
			mergedPath := filepath.Join(v.DataDir, crop, DirMergedVcfs)
			mergedListing := listDir(mergedPath)
			if !mergedListing.Exists || mergedListing.Err != "" {
				continue
			}

			refs := append([]string{}, mergedListing.SubDirs...)
			sortStrings(refs)

			for _, ref := range refs {
				if IsPruned(ref) {
					continue
				}
				refMergedPath := filepath.Join(mergedPath, ref)
				tagListing := listDir(refMergedPath)
				if tagListing.Err != "" {
					continue
				}

				tags := append([]string{}, tagListing.SubDirs...)
				sortStrings(tags)

				for _, tag := range tags {
					if IsPruned(tag) {
						continue
					}
					caller, merger, ok := SplitCallerMergerTag(tag)
					if !ok {
						continue
					}

					dirPath := filepath.Join(refMergedPath, tag)
					dirList := listDir(dirPath)
					if dirList.Err != "" {
						col.problem(Nonconformance{
							Kind: KindNonconformance, Rule: RuleUnreadableDir,
							Severity: SeverityError, Volume: v.Name, Crop: crop,
							Path: dirPath, Detail: dirList.Err,
						})
						continue
					}

					var matcher *unitMatcher
					var expected []Unit
					if info, hasInfo := refIndex[refKey{v.Name, crop, ref}]; hasInfo {
						matcher = info.matcher
						expected = info.units
					}

					set, problems := inspectJointDir(dirList, v.Name, crop, ref, tag, caller, merger, JointStandard, matcher, expected,
						opts.Validate == ValidateAll, opts.Quick)
					jointRows = append(jointRows, set)
					col.problem(problems...)
					col.counted(len(dirList.Files), dirList.TotalBytes())
				}
			}
		}
	}

	sort.Slice(jointRows, func(i, j int) bool {
		if jointRows[i].Volume != jointRows[j].Volume {
			return jointRows[i].Volume < jointRows[j].Volume
		}
		if jointRows[i].Crop != jointRows[j].Crop {
			return jointRows[i].Crop < jointRows[j].Crop
		}
		if jointRows[i].Reference != jointRows[j].Reference {
			return jointRows[i].Reference < jointRows[j].Reference
		}
		return jointRows[i].Tag < jointRows[j].Tag
	})

	return jointRows
}

func writeAll(w *Writer, volumes []resolvedVolume, cropRows []Crop, refIndex map[refKey]*refInfo, col *collector, jointRows []JointSet, opts Options) error {
	cropsByVol := map[string][]string{}
	for _, c := range cropRows {
		cropsByVol[c.Volume] = append(cropsByVol[c.Volume], c.Crop)
	}

	for _, v := range volumes {
		volRow, problem := capacityFor(v)
		if problem != nil {
			col.problem(*problem)
		}
		volRow.Crops = uniqueSorted(cropsByVol[v.Name])
		if err := w.Write(KindVolume, volRow); err != nil {
			return err
		}
	}

	sort.Slice(cropRows, func(i, j int) bool {
		if cropRows[i].Volume != cropRows[j].Volume {
			return cropRows[i].Volume < cropRows[j].Volume
		}
		return cropRows[i].Crop < cropRows[j].Crop
	})
	for _, c := range cropRows {
		if err := w.Write(KindCrop, c); err != nil {
			return err
		}
	}

	refKeys := make([]refKey, 0, len(refIndex))
	for k := range refIndex {
		refKeys = append(refKeys, k)
	}
	sort.Slice(refKeys, func(i, j int) bool {
		if refKeys[i].volume != refKeys[j].volume {
			return refKeys[i].volume < refKeys[j].volume
		}
		if refKeys[i].crop != refKeys[j].crop {
			return refKeys[i].crop < refKeys[j].crop
		}
		return refKeys[i].reference < refKeys[j].reference
	})

	for _, k := range refKeys {
		info := refIndex[k]
		if err := w.Write(KindReference, info.row); err != nil {
			return err
		}
		if opts.EmitBatchSeqIDs && info.batch != nil {
			if err := w.Write(KindReferenceBatch, *info.batch); err != nil {
				return err
			}
		}
	}

	sort.Slice(col.samples, func(i, j int) bool {
		if col.samples[i].Volume != col.samples[j].Volume {
			return col.samples[i].Volume < col.samples[j].Volume
		}
		if col.samples[i].Crop != col.samples[j].Crop {
			return col.samples[i].Crop < col.samples[j].Crop
		}
		if col.samples[i].Year != col.samples[j].Year {
			return col.samples[i].Year < col.samples[j].Year
		}
		return col.samples[i].Sample < col.samples[j].Sample
	})
	for _, s := range col.samples {
		if err := w.Write(KindSample, s); err != nil {
			return err
		}
	}

	sort.Slice(col.junk, func(i, j int) bool { return col.junk[i].Path < col.junk[j].Path })
	for _, junk := range col.junk {
		if err := w.Write(KindJunkCandidate, junk); err != nil {
			return err
		}
	}

	sort.Slice(col.sampleReferences, func(i, j int) bool {
		if col.sampleReferences[i].Volume != col.sampleReferences[j].Volume {
			return col.sampleReferences[i].Volume < col.sampleReferences[j].Volume
		}
		if col.sampleReferences[i].Crop != col.sampleReferences[j].Crop {
			return col.sampleReferences[i].Crop < col.sampleReferences[j].Crop
		}
		if col.sampleReferences[i].Year != col.sampleReferences[j].Year {
			return col.sampleReferences[i].Year < col.sampleReferences[j].Year
		}
		if col.sampleReferences[i].Sample != col.sampleReferences[j].Sample {
			return col.sampleReferences[i].Sample < col.sampleReferences[j].Sample
		}
		return col.sampleReferences[i].Reference < col.sampleReferences[j].Reference
	})
	for _, sr := range col.sampleReferences {
		if err := w.Write(KindSampleReference, sr); err != nil {
			return err
		}
	}

	for _, js := range jointRows {
		if err := w.Write(KindJointSet, js); err != nil {
			return err
		}
	}

	sort.Slice(col.problems, func(i, j int) bool {
		if col.problems[i].Volume != col.problems[j].Volume {
			return col.problems[i].Volume < col.problems[j].Volume
		}
		if col.problems[i].Crop != col.problems[j].Crop {
			return col.problems[i].Crop < col.problems[j].Crop
		}
		if col.problems[i].Year != col.problems[j].Year {
			return col.problems[i].Year < col.problems[j].Year
		}
		if col.problems[i].Sample != col.problems[j].Sample {
			return col.problems[i].Sample < col.problems[j].Sample
		}
		if col.problems[i].Rule != col.problems[j].Rule {
			return col.problems[i].Rule < col.problems[j].Rule
		}
		return col.problems[i].Path < col.problems[j].Path
	})
	for _, p := range col.problems {
		if err := w.Write(KindNonconformance, p); err != nil {
			return err
		}
	}

	return nil
}
