# `variants` package redesign — plan

Signatures and responsibilities only. No implementations.

## 1. Design goals

The current package implements the same biology twice (`VariantCallingDir` vs
`VariantCalling`) with different directory layouts, different resume strategies and
different concurrency models — and only one of the two reaches a merged VCF. Every
bug found in review traces back to that duplication.

Five principles drive the redesign:

1. **One engine, three front-ends.** Discovery is the *only* thing that differs
   between `data_dir`, `config` and inline mode. All three normalise to a single
   `RunPlan`, and everything downstream is shared code.
2. **One unit of work.** A chromosome and the contigs bundle are the same thing —
   a `SeqUnit`. The current code has four near-copies of the
   call → merge → filter sequence (gatk/dv × chroms/contigs); there will be one.
3. **Paths are computed in exactly one place.** Every filename mismatch found in
   review (`.ALL.joint.vcf.gz` vs `.all.vcf.gz`, `species/VCFs/refVer` vs
   `species/refVer/VCFs`) exists because writers and readers each built paths
   inline. A `Layout` interface owns all of it.
4. **Side effects are injectable.** External tools go through a `Runner`
   interface, so the whole pipeline is testable without GATK, DeepVariant,
   GLnexus, Docker or bcftools installed.
5. **Pure core, thin shell.** Command construction, filter predicates, dict
   parsing, path building and resume logic are pure functions. Only a handful of
   functions touch the filesystem or spawn processes.

### Bugs this design structurally prevents

| Bug (current code) | Prevented by |
|---|---|
| Final VCF built from goroutine-completion order, not reference order | `SeqUnit.Order`, concat sorts on it |
| `LatestVCFDir` returns a chromosome dir instead of a timestamp dir | `IsCohortDirName` regex in `paths.go` |
| `--skip-verification` yields silently empty merges | `Validator` interface; skipping swaps the impl, never the control flow |
| Sample-name comparison order-sensitive (`slices.Equal`) | `SampleIDsMatch` is set-based |
| Resume trusts the journal over the filesystem | `Executor.Do` requires journal **and** outputs **and** validity |
| Resume across incompatible settings | `Fingerprint` guard on journal open |
| Config typo silently ignored | strict YAML decoding |
| Two output-root layouts | single `Layout` per plan |
| `maxWorkers < 1 → 2` oversubscribes small machines | `resolveConcurrency`, unit-tested |
| Every failure exits 0 | all entry points return `error`; `cmd` uses `RunE` |

---

## 2. File layout

One package, many small files. No subpackages initially (avoids export churn);
`vcfio` can be extracted later if it earns its own tests.

```
variants/
  doc.go            package overview, mode summary
  types.go          domain types, no behaviour
  plan.go           the three front-ends → RunPlan
  configfile.go     YAML config schema + loader
  reference.go      fasta/dict resolution, SeqUnit construction
  discover.go       sample / alignment / gVCF discovery
  paths.go          Layout interface + implementations
  runner.go         Command, Runner, ExecRunner, FakeRunner
  journal.go        append-only resume log
  executor.go       skip / log / run wrapper (one place)
  validate.go       integrity checks
  caller.go         gVCF creation (GATK, DeepVariant)
  merger.go         joint genotyping (GATK, GLnexus)
  filter.go         hard filtering
  pipeline.go       orchestration
  cleanup.go        workspace / tmp removal
  testdata/         .dict, VCF headers, configs, golden argv
  *_test.go
```

---

## 3. Core types — `types.go`

Data only; behaviour lives in the files that own it.

```go
type ReadType int                  // ReadShort, ReadLong — replaces the "lr" name-suffix hack
type StageName string              // typed stage constants, no bare strings
type Mode int                      // ModeDataDir, ModeConfig, ModeInline

type Sequence struct               // one .dict entry: ID, Length
type SeqUnit struct                // unit of parallelism: Label, Seqs, IsContigs, Order
type Sample struct                 // SampleID, AlignmentPath, IndexPath, HomeDir, ReadType
type Reference struct              // FastaPath, DictPath, FaiPath, Sequences

type RunPlan struct                // THE normalised run description (see §5)
type UnitResult struct             // Label, Order, Gvcfs, JointVCF, FilteredVCF, Err
type Result struct                 // FinalVCF, Units, Skipped, Failed
```

**`SeqUnit.Order`** carries dict position. It is the single fix for the
final-VCF ordering bug: the concat stage sorts on it instead of trusting the
order goroutines happened to finish in.

**`Sample.HomeDir`** is where per-sample outputs live in `data_dir` mode
(`<CramDir>`), and empty in flat mode. This one field lets a single `Layout`
serve both trees.

Functions:

- `func (u SeqUnit) SeqIDs() []string` — sequence IDs for `-L` arguments.
- `func (u SeqUnit) String() string` — stable display label.
- `func (r Reference) Chromosomes() []Sequence` — sequences above the contig threshold.
- `func (r Reference) Contigs() []Sequence` — the remainder.

---

## 4. The three keystones

### 4.1 `Layout` — all path construction (`paths.go`)

Every path the pipeline reads or writes comes from here. Two implementations
cover the two directory conventions; nothing else in the package joins paths.

```go
type Layout interface {
    GvcfDir(s Sample, u SeqUnit) string
    GvcfPath(s Sample, u SeqUnit) string
    JointVCFPath(u SeqUnit) string
    FilteredVCFPath(u SeqUnit) string
    FinalVCFPath() string
    GenomicsDBPath(u SeqUnit) string
    SampleMapPath(u SeqUnit) string
    TempDir(u SeqUnit) string
    CohortDir() string
    JournalPath() string
}
```

- `func NewSampleTreeLayout(dataDir, species, refVer, callerSubdir string, cohort string) Layout`
  — `data_dir` mode. gVCFs beside each sample at
  `<Sample.HomeDir>/<refVer>/<callerSubdir>/`; cohort outputs under a
  timestamped dir.
- `func NewRunDirLayout(outDir, species, refVer, callerSubdir string) Layout`
  — inline/config mode. Everything under one `--out` root.
- `func SanitizeLabel(label string) string` — the `"." → "_"` rule, defined once.
- `func CohortDirName(t time.Time) string` — timestamp format, defined once.
- `func IsCohortDirName(name string) bool` — strict `^\d{8}_\d{6}$`; stops
  chromosome directories from winning the "latest cohort" sort.
- `func LatestCohortDir(parent string) (string, error)` — newest valid cohort dir.
- `func CreateCohortDir(parent string, now time.Time) (string, error)` — takes
  `now` as a parameter so tests are deterministic; retries on same-second collision.

### 4.2 `Runner` — injectable process execution (`runner.go`)

Commands are argv slices, never interpolated shell strings. This fixes
paths-with-spaces and makes every command assertion a unit test.

```go
type Command struct {
    Name       string
    Args       []string
    StdoutFile string      // for the `> out.bcf` redirections; runner owns the file
    PipeTo     *Command    // for glnexus | bcftools | bgzip
    Dir        string
    Env        []string
}

type Runner interface {
    Run(ctx context.Context, c Command) error
}
```

- `func NewExecRunner(verbose bool, log *slog.Logger) Runner` — real execution.
- `func NewFakeRunner() *FakeRunner` — records calls, returns stubbed results.
- `func (f *FakeRunner) Calls() []Command`
- `func (f *FakeRunner) StubError(nameOrArgSubstring string, err error)`
- `func (f *FakeRunner) StubCreatesFile(match, path string)` — simulates a tool
  producing its output, so `Executor`'s output checks can be exercised.
- `func (c Command) Display() string` — for logs and the journal only; never executed.
- `func CheckTools(names ...string) error` — dependency preflight, moved out of `cmd`.

### 4.3 `Executor.Do` — the one skip/log/run wrapper (`executor.go`)

The current code repeats a `jlog.Info(STARTED)` / `RunCmd` / `jlog.Error(FAILED)` /
`jlog.Info(COMPLETED)` block roughly ten times. It becomes one function.

```go
type Stage struct {
    Name     StageName
    Sample   string          // "ALL" for cohort-level stages
    Unit     string
    Outputs  []string        // artifacts that must exist afterwards
    Validate func() error    // optional deeper integrity check
    Run      func(context.Context) error
}

type Executor struct { /* journal, runner, validator, logger, force */ }
```

- `func NewExecutor(j *Journal, r Runner, v Validator, log *slog.Logger, force bool) *Executor`
- `func (e *Executor) Do(ctx context.Context, s Stage) (skipped bool, err error)`
  — resolves resume by requiring **all three**: journal says COMPLETED, every
  `Outputs` path exists, and `Validate` (if set) passes. Otherwise it logs
  STARTED, runs, records COMPLETED or FAILED. This is where the
  log-says-done-but-file-is-gone class of failure dies.
- `func (e *Executor) Runner() Runner`
- `func (e *Executor) Log() *slog.Logger`

---

## 5. `RunPlan` and the three front-ends — `plan.go`

`RunPlan` is the contract between the front-ends and the engine.

```go
type Inputs struct { /* the raw flag bag from cmd — one struct, not 20 params */ }

type RunPlan struct {
    Mode        Mode
    Species     string
    RefVer      string
    Reference   *Reference
    Samples     []Sample
    Units       []SeqUnit
    Layout      Layout
    Caller      VariantCaller
    Merger      Merger
    Thresholds  FilterThresholds
    Concurrency Concurrency
    Skip        SkipFlags        // NoMerging, NoHardFilter, SkipVerification, Quick
    Issues      []DiscoveryIssue // non-fatal: missing/duplicate/corrupt inputs
}
```

- `func BuildPlan(in Inputs) (*RunPlan, error)` — the dispatcher. **Errors if
  more than one mode is supplied**, instead of today's silent
  config-beats-datadir-beats-inline precedence.
- `func planFromDataDir(in Inputs) (*RunPlan, error)` — discovers samples by
  walking the tree; installs `SampleTreeLayout`.
- `func planFromConfig(in Inputs) (*RunPlan, error)` — loads the config file,
  applies it to `Inputs`, then delegates to the inline builder; installs
  `RunDirLayout`.
- `func planFromInline(in Inputs) (*RunPlan, error)` — builds `[]Sample` from
  `--bam`; installs `RunDirLayout`.
- `func (p *RunPlan) Validate() error` — **every** precondition in one place:
  species/refVer present, reference + dict readable, ≥1 sample, no duplicate
  sample IDs, caller/merger compatible, thresholds sane, output root writable.
  Replaces the ~20 scattered `fmt.Println` + `return` checks.
- `func (p *RunPlan) Describe() string` — dry-run summary: samples, units,
  tools, resolved output paths. Printed before any compute is launched.
- `func (p *RunPlan) Fingerprint() string` — hash of resume-critical settings.
- `func MergeOnlyPlan(in Inputs) (*RunPlan, error)` — entry point for the
  `MergeGvcfs` command: same plan, gVCF-creation stages pre-satisfied, gVCFs
  discovered rather than produced.

### What actually differs between the three modes

| | `data_dir` | `config` | inline |
|---|---|---|---|
| Sample source | tree walk | `bam:` entries | `--bam` flags |
| Reference | `ResolveFastaPath` w/ genomes dir | config field | `--reference` |
| Layout | `SampleTreeLayout` | `RunDirLayout` | `RunDirLayout` |
| Output root | `<data>/<species>/<refVer>/VCFs` | config `OutputDir` | `--out` |
| Everything after `BuildPlan` | identical | identical | identical |

---

## 6. Config file — `configfile.go`

Replaces the hand-rolled `Key: value` parser. Strict YAML fixes the silent-typo
failure and the `SplitN(line, ":", 2)` breakage on `C:\...` paths.

- `type FileConfig struct` — yaml-tagged; reference, species, ref_version,
  output_dir, bams, gvcfs, threads, caller, merger, filter thresholds.
- `func LoadFileConfig(path string) (*FileConfig, error)` — strict decoding;
  unknown fields are an error, not a shrug.
- `func (c *FileConfig) ApplyTo(in *Inputs) error` — config is authoritative,
  explicitly-set flags override, and provenance is recorded per field so
  `Describe()` can show where each value came from.
- `func LoadLegacyConfig(path string) (*FileConfig, error)` — reads the old
  format for one release, emitting a deprecation warning.
- `func WriteTemplateConfig(w io.Writer) error` — backs `CreateTemplate`.

---

## 7. Reference handling — `reference.go`

- `func ResolveFastaPath(explicit, genomesDir, species, refVer string) (string, error)`
  — one resolution rule for all three modes (currently only `data_dir` has one).
- `func LoadReference(fastaPath string) (*Reference, error)` — derives `.dict`
  and `.fai`, requires both, parses sequences.
- `func parseDict(r io.Reader) ([]Sequence, error)` — pure; `@SQ`/`SN`/`LN`.
- `func (r *Reference) SeqUnits(contigMinLength int64) []SeqUnit` — one unit per
  chromosome plus a single contigs unit, each stamped with `Order`.
- `func (r *Reference) UnitByLabel(label string) (SeqUnit, bool)` — for
  merge-only runs that resume from existing per-chromosome outputs.

---

## 8. Discovery — `discover.go`

- `type DiscoveryIssue struct` — Sample, Unit, Reason (`Missing`, `Multiple`, `Corrupt`).
- `func DiscoverSamples(dataDir, species, refVer string, preferBQSR bool) ([]Sample, []DiscoveryIssue, error)`
  — the tree walk; assigns `HomeDir` and infers `ReadType`.
- `func findAlignment(sampleDir, refVer string, preferBQSR bool) (string, []DiscoveryIssue)`
  — cram-then-bam resolution for one sample.
- `func SamplesFromPaths(paths []string) ([]Sample, error)` — inline/config
  construction; reads real sample IDs from headers rather than guessing from filenames.
- `func FindExistingGvcfs(l Layout, samples []Sample, u SeqUnit) (map[string]string, []DiscoveryIssue)`
  — per-unit gVCF inventory keyed by sample ID. **Returns a per-unit map**, which
  is the grouping the merge stage needs and the current flat glob fails to provide.
- `func ReadSampleIDs(vcfPath string) ([]string, error)` — pure-Go `#CHROM`
  header read; replaces the `bcftools query -l` subprocess.
- `func SampleIDsMatch(a, b []string) bool` — set comparison; replaces the
  order-sensitive `slices.Equal`.
- `func DuplicateSampleIDs(ids []string) []string` — catches the duplicate-sample
  condition that makes `GenomicsDBImport` fail late and cryptically.

---

## 9. Resume journal — `journal.go`

- `type Event struct` — Time, Stage, Sample, Unit, Status, Cmd, Outputs, Err.
- `type Journal struct`
- `func OpenJournal(path, fingerprint string) (*Journal, error)` — opens
  append-only JSONL, replays prior events, and **refuses to resume when the
  fingerprint differs** (today a rerun with a different `--caller` or different
  filter thresholds silently resumes on incompatible outputs).
- `func (j *Journal) Completed(stage StageName, sample, unit string) bool`
- `func (j *Journal) Record(ev Event) error`
- `func (j *Journal) Attempts(stage StageName, sample, unit string) int` — enables
  a retry cap instead of unbounded re-attempts.
- `func (j *Journal) Close() error`
- `func parseEvents(r io.Reader) ([]Event, []error)` — pure; tolerant of
  malformed lines and free of the unchecked type assertions that can panic today.

---

## 10. Validation — `validate.go`

- `type Validator interface { ValidateVCF(path string, quick bool) error; ValidateAlignment(path, ref string, quick bool) error }`
- `func NewToolValidator(r Runner) Validator` — shells out for deep checks.
- `func NewNoopValidator() Validator` — what `--skip-verification` swaps in.
  Skipping changes the *implementation*, never the control flow, so it can no
  longer cause the empty-cohort bug.
- `func CheckBGZFTrailer(path string) error` — pure-Go BGZF EOF-block check;
  catches truncated files cheaply, no subprocess.
- `func CheckIndexPresent(vcfPath string) error` — `.tbi` exists and is not older
  than its VCF.

---

## 11. gVCF creation — `caller.go`

Command *construction* is pure and separately testable from command *execution*.

```go
type CallRequest struct {
    Sample    Sample
    Unit      SeqUnit
    Reference *Reference
    OutPath   string
    TempDir   string
    Threads   int
    LogLevel  string
}

type VariantCaller interface {
    Name() string
    Stage() StageName
    GvcfSubdir() string                            // "gatk_gvcfs" | "dv_gvcfs"
    BuildCommand(CallRequest) (Command, error)     // pure
}
```

- `func NewGATKCaller(javaOpts, logLevel string) VariantCaller`
- `func NewDeepVariantCaller(version, defaultModel string) VariantCaller`
- `func (d *deepVariantCaller) modelFor(s Sample) string` — reads
  `Sample.ReadType`, honouring the configured default. Retires the
  `strings.HasSuffix(sample, "lr")` convention that silently overrode `--model-type`.
- `func ResolveCaller(name string, in Inputs) (VariantCaller, error)`
- `func CallGvcf(ctx context.Context, e *Executor, c VariantCaller, req CallRequest) (string, error)`
  — wraps `BuildCommand` in a `Stage` and delegates to `Executor.Do`.

---

## 12. Joint genotyping — `merger.go`

```go
type MergeRequest struct {
    Unit          SeqUnit
    Gvcfs         []string
    Reference     *Reference
    OutPath       string
    WorkspaceDir  string
    SampleMapPath string
    TempDir       string
    LogLevel      string
}

type Merger interface {
    Name() string
    Merge(ctx context.Context, e *Executor, req MergeRequest) (string, error)
}
```

- `func NewGATKMerger(javaOpts string, batchSize int) Merger` — two `Executor`
  stages: `GenomicsDBImport`, then `GenotypeGVCFs`.
- `func NewGLnexusMerger(preset string) Merger` — three stages: `glnexus_cli`,
  `bcftools view | bgzip`, `tabix`.
- `func WriteSampleMap(path string, gvcfs []string) error` — pure-Go via
  `ReadSampleIDs`; **errors** on duplicates or on an input whose header yields no
  sample, rather than silently omitting it as the current suffix check does.
- `func ResolveMerger(name string, in Inputs) (Merger, error)`
- `func ValidateCallerMergerPair(c VariantCaller, m Merger) error` — the
  "DeepVariant requires GLnexus" rule, stated once, checked in `RunPlan.Validate`.

---

## 13. Hard filtering — `filter.go`

The predicate is pure, so it is the cheapest high-value thing in the package to
test exhaustively.

- `type FilterThresholds struct` — SNP and INDEL cutoffs.
- `func DefaultThresholds() FilterThresholds`
- `func LightThresholds() FilterThresholds`
- `func (t FilterThresholds) Validate() error` — rejects contradictory cutoffs.
- `type VariantClass int` — `ClassSNP`, `ClassIndel`, `ClassMNP`, `ClassOther`.
- `func ClassifyVariant(v *vcfgo.Variant) VariantClass` — pure.
- `func (t FilterThresholds) Evaluate(v *vcfgo.Variant) (pass bool, failedOn []string)`
  — pure; returns which annotations failed so the FILTER field is informative
  instead of a bare drop.
- `func HardFilterVCF(ctx context.Context, e *Executor, req FilterRequest) (string, error)`
  — streams input to a bgzipped output, then indexes.
- `type FilterReport struct` — per-class in/out counts, so a run is auditable.
- `func (r FilterReport) String() string`

---

## 14. Orchestration — `pipeline.go`

- `func Run(ctx context.Context, p *RunPlan) (*Result, error)` — the single
  top-level entry. Opens the journal, builds the executor, fans out over
  `p.Units`, barriers, concatenates, returns a populated `Result`.
- `func runUnit(ctx context.Context, e *Executor, p *RunPlan, u SeqUnit) UnitResult`
  — the *one* per-unit path: call gVCFs for every sample → merge → hard filter.
  Chromosomes and contigs both go through it, so the four duplicated blocks in
  today's code collapse to one.
- `func callUnitGvcfs(ctx context.Context, e *Executor, p *RunPlan, u SeqUnit) ([]string, []DiscoveryIssue, error)`
  — per-sample fan-out inside a unit.
- `func concatenate(ctx context.Context, e *Executor, p *RunPlan, units []UnitResult) (string, error)`
  — **sorts by `SeqUnit.Order`** before writing the input list, then runs
  `gatk MergeVcfs` (or copies when there is one input), then indexes. Fixes the
  reference-order bug and the missing final index.
- `func resolveConcurrency(cores, threadsPerJob, jobsFlag int) int` — pure;
  never returns more workers than a small machine can honour.
- `func newWorkerPool(n int) *workerPool` — bounded fan-out with context cancellation.
- `func collectErrors(results []UnitResult) error` — aggregates via
  `errors.Join` so one run surfaces *all* failing units, not just the first.

## 15. Cleanup — `cleanup.go`

- `func Cleanup(ctx context.Context, p *RunPlan, keep bool) error` — removes
  GenomicsDB workspaces and temp dirs once a unit succeeds. These are large and
  currently never reclaimed.

---

## 16. Test plan

No external bioinformatics tool is needed for any test: `FakeRunner` intercepts
execution and `t.TempDir()` provides real filesystem behaviour. Target: the pure
core at high coverage, the shell at smoke level.

| File | Covers | Pins down |
|---|---|---|
| `paths_test.go` | every `Layout` method under both implementations, table-driven | the two output-root conventions and every filename mismatch |
| `reference_test.go` | `parseDict` on golden `.dict`, chrom/contig split, `SeqUnits` ordering | `Order` is dict order |
| `configfile_test.go` | strict YAML, unknown field rejection, Windows paths, comments, `ApplyTo` precedence | silent-typo and `C:\` parse failures |
| `plan_test.go` | all three front-ends produce equivalent plans from equivalent inputs; conflicting modes error; `Validate` rejects each bad precondition | silent mode precedence |
| `journal_test.go` | replay, `Completed`, fingerprint mismatch refusal, malformed lines, FAILED→retry, `Attempts` | resume across incompatible settings; the panicking type assertions |
| `executor_test.go` | skip when journal+outputs+validation agree; **run when journal says done but output is missing**; failure recording | the log-trusting resume bug |
| `validate_test.go` | `CheckBGZFTrailer` on good/truncated fixtures; noop validator | — |
| `caller_test.go` | `BuildCommand` argv golden files for GATK and DeepVariant; `modelFor` by `ReadType` | the `"lr"` suffix override |
| `merger_test.go` | stage sequence and argv for both mergers; `WriteSampleMap` duplicate + no-sample errors | silent gVCF omission |
| `filter_test.go` | `ClassifyVariant` and `Evaluate` across a synthetic variant matrix, boundary values per annotation | the highest-value pure logic in the package |
| `concat_test.go` | **input list is in reference order regardless of completion order**; single-input copy path carries the index | the mis-sorted final VCF |
| `discover_test.go` | tree walk over a synthetic `t.TempDir()` tree; missing/multiple/corrupt classification; `SampleIDsMatch`; `DuplicateSampleIDs` | order-sensitive comparison |
| `pipeline_test.go` | full cohort with `FakeRunner`: happy path, one unit failing while others complete, resume from a partial journal, `--no-merging`, `--skip-verification` | the empty-cohort bug; all-or-nothing merging |
| `concurrency_test.go` | `resolveConcurrency` table incl. cores<threads | `maxWorkers < 1 → 2` |

`testdata/` holds: a small `.dict`, VCF/gVCF headers (single- and multi-sample,
plus sites-only), a truncated `.vcf.gz`, valid and malformed journals, legacy and
YAML configs, and golden argv files.

---

## 17. Sequencing

Each phase ends compiling and green, so it can be reviewed on its own.

| Phase | Deliverable |
|---|---|
| 1 | `types.go`, `paths.go`, `reference.go` + tests. Pure, no I/O. |
| 2 | `runner.go`, `journal.go`, `executor.go`, `validate.go` + tests. The shared machinery. |
| 3 | `caller.go`, `merger.go`, `filter.go` + tests. Argv construction and the filter predicate. |
| 4 | `pipeline.go`, `cleanup.go` + tests. `FakeRunner` end-to-end. |
| 5 | `discover.go`, `configfile.go`, `plan.go` + tests. The three front-ends. |
| 6 | Rewire `cmd/VariantCalling.go` and `cmd/MergeGvcfs.go` to `RunE`; wire `--jobs`; add `--dry-run`. |
| 7 | One real end-to-end run against a small dataset; then delete `*.backup` and `mergeGVCFs.go.1`. |

### Not carried over

- `MergeGvcfsGATK` / `MergeGvcfsGATKDir` and `MergeGvcfsGlnexus` /
  `MergeGvcfsGlnexusDir` — four functions collapsing to two `Merger` implementations.
- `VariantCallingDir` and `VariantCalling` as separate engines — one `Run`.
- `CreateGvcfGATKog`, `CreateGvcfDVog`, `runCmd`, `FindBamOrVcfs`,
  `glnexusSript` — dead or superseded.
- Shell-string command building throughout — replaced by argv.
- `utils.ReadConfig`'s variant-calling fields — superseded by `FileConfig`
  (other commands keep using `utils.Config` until they are migrated too).
