# TODO

## Modular access to the alignment steps (parked)

Context: outside `--data-dir` mode, `AlignReads` is all-or-nothing. There is no way to
run one step — mark duplicates, convert to CRAM, index — against a file of your own.
The question was whether to add `MarkDuplicates` / `BamToCram` / `IndexBam` subcommands
or to tell people to call `gatk`/`samtools` directly.

Conclusion: the wrappers are not worth exposing as three sibling commands. What they
know beyond the raw tool is *where things spill, how files are named, and how writes are
made atomic* — `--TMP_DIR` on the spill disk, `-c` for CSI rather than BAI,
`version=3.0`, `.partial`+rename. Real knowledge, but a rare enough need standalone that
three commands (plus flags, help and the `flags_test.go` conventions) is a poor trade.

### 1. `--print-commands` on AlignReads

Plan the samples exactly as now, print each step's command with the real paths filled
in, and exit without running anything. Gives a person the correct
`gatk MarkDuplicates -R … --TMP_DIR /var/tmp/…` for *their* file, which they can copy,
edit and run — strictly more useful than a wrapper, because it can be modified. One
flag, no new surface, and it makes the conventions visible instead of hiding them.

### 2. `Ingest` — bring a foreign alignment to estate convention

The one operation with no `samtools` one-liner equivalent: verify contigs against the
reference dictionary, rewrite `@RG` to the sample, mark duplicates, convert to CRAM,
index. It is what people already do by hand — `MENINA_LONG_READS.aligned.cram`
(`warehouse/inspect.go:175`) exists because someone did those steps outside the tool and
stopped part-way.

The implementation already exists from the long-read work: `verifyAdoptable` and
`normaliseAdoptedAlignment` in `alignmentdir/longreads.go`, plus the existing
markdup → cram → index chain. Exposing it is mostly flag parsing.

### Rules if either is built

- A subcommand must call the **same exported function the pipeline calls**, and contain
  nothing but flag parsing. This repo already carries three alignment scanners
  (`inspectSampleBamDir`, `ScanAlignments`, the warehouse's) and two unrelated BQSR
  implementations (`alignment/baseRecal.go` vs `alignmentdir`'s `runBQSR`); both started
  as a reasonable small separate thing.
- If more than one lands, nest them — `genome-whisperer bam index|to-cram|ingest`.
  The root help already lists 22 commands and nothing is nested today.

## Known duplication worth collapsing

- **Two BQSR implementations.** `alignment/baseRecal.go` (`Recalibrate`, `BootstrapBqsr`,
  `DbSnpBqsr`, used by the `BQSR` subcommand) and `alignmentdir`'s `runBQSR` /
  `createBootstrapKnownSites`. `alignmentdir` calls none of the former. Only the
  data-dir one scatters over intervals.
- **Three alignment scanners**, as above. `ScanAlignments` also still carries cobra's
  placeholder description.
- **`createBootstrapKnownSitesOld`** in `alignmentdir/dirAlign.go` is dead.
- **Unfinished commands**: `GeneSpace`, `CreateTemplate`, `ScanAlignments` still have
  "A brief description of your command".

## Input-mode inconsistency

`--data-dir` works on `AlignReads`, `CreateGvcfs`, `ScanAlignments`. `--config` works on
eight commands. `MergeGvcfs` and `VariantCalling` take config but not data-dir, so
"estate mode" is not a uniform concept across the CLI.
