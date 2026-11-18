# genome-whisperer

A command-line toolkit for plant genomics, written in Go. It takes sequencing data from raw reads to ranked candidate genes:

```
reads ──► AlignReads ──► VariantCalling ──► VariantAnnotation ──► GoBSAseq ──► GeneSpace / MarkerDesignPipeline
          (+ BQSR)      (gVCFs → joint VCF     (snpEff + PRG +      (QTLs)       (candidate genes / markers)
                         → hard filtering)      descriptions =
                                                "superVCF" table)
```

Each stage is a subcommand. Each can run on its own with your own files, or be chained through the standard data directory.

- **Alignment**: short reads (bwa-mem, bwa-mem2, bowtie2) and long reads (pbmm2), duplicate marking and GATK BQSR. Output is BAM or CRAM.
- **Variant calling**: GATK HaplotypeCaller or DeepVariant, merged with GATK GenomicsDB or GLnexus, then hard-filtered. Interrupted runs resume where they stopped.
- **Annotation**: snpEff effects, plus resistance-gene (PRGdb) BLAST hits and gene descriptions, combined into one tab-separated "superVCF" table.
- **QTL mapping**: bulk segregant analysis (BSA-seq) with one or two bulks, with or without parents.
- **Candidate genes**: GeneSpace ranks the genes in a QTL by genotype–phenotype concordance, variant impact and trait-specific function.
- **Markers**: MarkerDesignPipeline designs markers in QTLs from the superVCF table.
- **Pangenome**: GoPan builds a pangenome by iterative mapping and assembly.
- **Data management**: ScanWarehouse and ScanAlignments take stock of sequencing storage and alignment state.

---

## Installation

Requires Go 1.26+.

```bash
git clone <repo-url> genome-whisperer
cd genome-whisperer
go build -o genome-whisperer .
# or: go install github.com/gmaffy/genome-whisperer@latest
```

### External tools

genome-whisperer calls standard bioinformatics tools. Each command checks for the tools it needs before it starts and tells you which are missing. Install only what you use:

| Used by | Tools |
|---|---|
| Alignment | `samtools`, `gatk`, `java`, plus `bwa`, `bwa-mem2`, `bowtie2` or `pbmm2` (+ `pbmarkdup`) |
| Variant calling | `gatk`, `bcftools`, `tabix`; DeepVariant needs `docker`; GLnexus needs `glnexus_cli`, `bgzip` |
| Annotation | `snpEff`, `java`, `gatk` |
| PRG hits | `blastp` (BLAST+) and a PRGdb protein BLAST database |
| GoPan | `gffread`, `MAC2.0`, `megahit`, `masurca`, `seqtk`, `bedtools`, `bwa`, `bowtie2` |

---

## General usage

```bash
genome-whisperer <Command> [flags]
genome-whisperer <Command> --help
```

### Shared flags

These flags are accepted by every command. Each command reads only the ones that apply to it.

| Flag | Meaning |
|---|---|
| `-r, --reference` | Reference genome FASTA |
| `--ref-version` | Reference genome version (e.g. `AD1.1`) |
| `-s, --species` | Species name |
| `-g, --genomes-dir` | Directory of prepared reference genomes |
| `-c, --config` | Config file (see [Config file](#config-file)) |
| `-d, --data-dir` | Root of the standard data directory |
| `-b, --bam` | BAM/CRAM file (repeatable) |
| `-o, --out-dir` | Output directory |
| `-l, --log` | Run log file |
| `--alignment-fmt` | `bam` or `cram` (default `bam`) |
| `-t, --threads` | Threads per job (default 4) |
| `-v, --verbose` | Stream external tool output |
| `--gatk-log-level` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--quick` / `--skip-verification` | Header-only integrity checks / no checks |
| `--aligner` | `bwa-mem`, `bwa-mem2` (default), `bowtie2`, `pbmm2` |
| `--preset` | pbmm2 preset: `SUBREAD`, `CCS`, `HIFI` (default), `ISOSEQ`, `UNROLLED` |
| `-k, --known-sites` / `--bootstrap` | BQSR known-sites VCF(s) / build them by bootstrapping |
| `--caller` | `gatk` (default) or `deepvariant` |
| `--merger` | `gatk` (default) or `glnexus` |
| `--model-type` | DeepVariant model: `WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA` |
| `--deepvariant-version` | Default `1.10.0` |
| `--merge-jobs`, `--local-work-dir`, `--no-merging` | Joint-genotyping concurrency, scratch dir, stop after gVCFs |
| `--database` | snpEff database name |
| `--bsaseq` | Also export AD/DP/GQ columns (needed for BSA-seq) |
| `--gff`, `--protein`, `--annotation-version` | Genome annotation inputs |
| `--gene-description-tsv` | `gene_id<TAB>description` file (optional wherever it is used, except `AddGeneDescriptions`) |
| `--prg` | PRG top-hits table from `GetProtPrgTopHits` |

### Input modes

Most pipeline commands accept their input in one of three ways. Pick exactly one:

1. **Data directory** (`--data-dir --species --ref-version`). Samples are discovered from the standard layout, and outputs are written beside them:
   ```
   <data-dir>/<species>/<project>/<sample>/reference_genomes/<ref-version>/bams/
   <data-dir>/<species>/<project>/<sample>/reference_genomes/<ref-version>/{gatk_gvcfs|dv_gvcfs}/
   <data-dir>/<species>/MERGED_VCFs/<ref-version>/<caller>_<merger>/
   ```
2. **Config file** (`--config`). Inputs are listed in a text file. Command-line flags win where both are given.
3. **Inline**: every input is given on the command line (`--bam`, `--reference`, `--out-dir`, …).

### Config file

Generate a template with `genome-whisperer CreateTemplate`, which writes `config.txt`:

```
OutputDir: <path to output directory>
Reference: <path to reference fasta file>
gff: <path to gff file>
proteins: <path to protein fasta file>
cds: <path to cds fasta file>
ReadPair: <forward reads> <reverse reads> <sample name> <library name>
SingeEndReads: <reads> <sample name> <library name>
bam: <path to bam/cram file>
BSAseqBam: <path to bam/cram file> <HIGH_PARENT|LOW_PARENT|HIGH_BULK|LOW_BULK>
VCF: <path to vcf file>
gvcf: <path to gvcf file>
```

For BSA-seq, the library name of each read pair is `HIGH_PARENT`, `LOW_PARENT`, `HIGH_BULK` or `LOW_BULK`.

---

## Commands

### Reference and annotation preparation

#### PrepareFastaFile
Indexes a reference with samtools, writes the GATK sequence dictionary and builds the aligner index. The reference can be plain or bgzipped.
```bash
genome-whisperer PrepareFastaFile --reference genome.fa --aligner bwa-mem2
```

#### CreateSnpEffDB
Builds a snpEff database from a genome, GFF3, protein FASTA and CDS FASTA.
```bash
genome-whisperer CreateSnpEffDB --reference genome.fa --gff genes.gff3 \
  --protein proteins.fa --cds cds.fa --species squash --annotation-version v1.1
```
`--reference`, `--gff`, `--protein` and `--cds` can also come from `--config`.

#### GetProtPrgTopHits
BLASTs a protein FASTA against PRGdb and keeps the best hit (highest bitscore) per protein. The full BLAST table is never written. The output, `<SPECIES><VERSION>_PRG_TOP_HITS.txt`, is the file `--prg` expects in every other command.
```bash
genome-whisperer GetProtPrgTopHits --protein proteins.fa --prg-db prgdb.fa \
  --species squash --annotation-version v1.1 --threads 16
```
| Flag | Default | |
|---|---|---|
| `--prg-db` | | blastp `-db` prefix (build with `makeblastdb -in prgdb.fa -dbtype prot -out prgdb.fa`) |
| `--evalue` | `1e-9` | e-value cutoff |
| `--max-target-seqs` | `0` (uncapped) | A cap can hide the best hit |

---

### Alignment

#### AlignReads
Aligns reads, marks duplicates and (optionally) runs BQSR.
```bash
# Inline, paired-end
genome-whisperer AlignReads -1 R1.fq.gz -2 R2.fq.gz --sample S1 --library L1 \
  --reference genome.fa --out-dir aln --aligner bwa-mem2 --bqsr --bootstrap

# Long reads
genome-whisperer AlignReads --single-end hifi.fq.gz --sample S1 --library L1 \
  --reference genome.fa --out-dir aln --aligner pbmm2 --preset HIFI

# Everything listed in a config file, or discovered from a data directory
genome-whisperer AlignReads --config config.txt
genome-whisperer AlignReads --data-dir /data --species squash --ref-version v1.1 --genomes-dir /genomes
```
| Flag | |
|---|---|
| `-1, --forward` / `-2, --reverse` | Paired reads |
| `--single-end` | Single-end or long reads |
| `--sample`, `--library` | Read-group sample and library |
| `--bqsr` | Run base quality score recalibration |

#### BQSR
Runs GATK BaseRecalibrator and ApplyBQSR on duplicate-marked BAMs. Without `--known-sites`, a known-sites set is bootstrapped from the data.
```bash
genome-whisperer BQSR --bam S1.rgmd.bam --reference genome.fa --known-sites known.vcf.gz
```

#### ScanAlignments
Reports the state of every sample's alignment files under a data directory, without changing anything. For each sample it lists the files with size, validity and index state, and the steps `AlignReads` would still run. It uses the same planning code as `AlignReads`, so the two always agree.
```bash
genome-whisperer ScanAlignments --data-dir /data --species squash --ref-version v1.1 --genomes-dir /genomes
```
Each sample gets one status:

| Status | Meaning |
|---|---|
| `COMPLETE` | Finished, nothing to do |
| `CLEANUP` | Outputs done, only intermediates to remove |
| `RESUME` | An alignment exists, later steps remain |
| `ALIGN` | Nothing usable, must align from reads |
| `AMBIGUOUS` | Several external alignments for a long-read sample; needs a person |
| `ERROR` | The bams directory could not be read |

Plans assume `AlignReads --bqsr`; pass `--bqsr=false` for a run without recalibration. Unrecognised files in a short-read sample's bams directory are flagged. `AlignReads` cleanup never removes them, so they keep the sample from counting as complete.

---

### Variant calling

#### VariantCalling
The full pipeline: per-chromosome gVCFs → joint genotyping → concatenation → hard filtering. Re-running skips valid outputs, so an interrupted run resumes.
```bash
# Data directory
genome-whisperer VariantCalling --data-dir /data --species cotton --ref-version AD1.1 \
  --genomes-dir /genomes --threads 8

# Inline BAMs
genome-whisperer VariantCalling --bam s1.cram --bam s2.cram --reference ref.fa \
  --out-dir vcfs --species cotton --ref-version AD1.1

# Long reads with DeepVariant (GLnexus is required as merger)
genome-whisperer VariantCalling --data-dir /data --species cotton --ref-version AD1.1 \
  --caller deepvariant --merger glnexus --model-type PACBIO
```
| Flag | Default | |
|---|---|---|
| `--bqsr` | `true` | Use the BQSR alignment (`--bqsr=false` uses the duplicate-marked one) |
| `--no-merging` / `--no-hard-filter` | | Stop after gVCFs / after merging |
| `--light-filtering` | | Relaxed thresholds |
| `--min-qd-snp`, `--min-qual-snp`, `--max-sor-snp`, `--max-fs-snp`, `--min-mq-snp`, `--min-mqranksum-snp`, `--min-readposranksum-snp` | 2, 30, 3, 60, 40, −12.5, −8 | SNP thresholds |
| `--min-qd-indel`, `--min-qual-indel`, `--max-fs-indel`, `--max-sor-indel`, `--min-readposranksum-indel` | 2, 30, 200, 10, −20 | INDEL thresholds |
| `--min-gq` | 20 | DeepVariant: minimum GQ in at least one sample |

Joint VCFs go into a folder named after the caller and merger combination (`gatk_gatk`, `gatk_glnexus`, `dv_glnexus`), so results from different combinations sit side by side.

#### CreateGvcfs
Only the gVCF stage of `VariantCalling`. Valid existing gVCFs are reused and corrupt ones are re-created.
```bash
genome-whisperer CreateGvcfs --data-dir /data --species cotton --ref-version AD1.1 --genomes-dir /genomes
```

#### MergeGvcfs
Only the joint-genotyping stage: merges gVCFs into one multi-sample VCF.
```bash
genome-whisperer MergeGvcfs --gvcf s1.g.vcf.gz --gvcf s2.g.vcf.gz --reference ref.fa \
  --out-dir vcfs --species cotton --ref-version AD1.1 --merger glnexus
```

#### HardFilterSNPs / HardFilterINDELs / HardFilterVCFs
Applies GATK-recommended hard filters to SNPs, to INDELs, or to both (`HardFilterVCFs` splits, filters and merges the result).
```bash
genome-whisperer HardFilterVCFs -V joint.vcf.gz
```

#### GetVariantType
Keeps only one variant type.
```bash
genome-whisperer GetVariantType -V joint.vcf.gz --variant-type SNP
```

#### CreateSuperVcf
`VariantCalling` with default filters, followed by `VariantAnnotation`, in one command.
```bash
genome-whisperer CreateSuperVcf --data-dir /data --species squash --ref-version v1.1 \
  --genomes-dir /genomes --database squash_v1.1 --prg SQUASHV1.1_PRG_TOP_HITS.txt \
  [--gene-description-tsv desc.tsv] [--bsaseq]
```

---

### Annotation

#### VariantAnnotation
Annotates VCFs with snpEff and writes the **superVCF table**: `<input>…_SUPER_VCF.tsv`.
```bash
genome-whisperer VariantAnnotation -V joint.filtered.vcf.gz --database squash_v1.1 \
  --prg SQUASHV1.1_PRG_TOP_HITS.txt [--gene-description-tsv desc.tsv] [--bsaseq]
```
The table's columns are:
- `CHROM POS REF ALT QUAL TYPE`
- `<sample>.GT` for each sample (plus `.AD .DP .GQ` with `--bsaseq`)
- `SNPEFF_EFFECT SNPEff_IMPACT SNPEFF_FUNCTIONAL_CLASS … SNPEFF_GENE_NAME SNPEFF_TRANSCRIPT_ID …`
- `PercIdent QLEN/MATCH_LEN Eval Description`

Without `--gene-description-tsv`, the description step is skipped and `Description` is `NA`.

#### RunSnpEff
Only the snpEff step: writes an annotated VCF and a TSV with the effect fields split into columns.
```bash
genome-whisperer RunSnpEff -V joint.vcf.gz --database squash_v1.1 [--bsaseq]
```

#### AddPrg / AddGeneDescriptions
Add PRG hit columns or gene descriptions to an existing snpEff table or VCF.
```bash
genome-whisperer AddPrg -V joint.snpEff_EFF.tsv --prg SQUASHV1.1_PRG_TOP_HITS.txt
genome-whisperer AddGeneDescriptions -V joint.snpEff_EFF.tsv --gene-description-tsv desc.tsv
```

---

### QTL mapping and candidate genes

#### GoBSAseq
Bulk segregant analysis to detect QTLs. It takes reads, BAMs (via `--config`) or a VCF. Supported designs:
1. Two bulks only: `-A` and `-B`.
2. Two bulks with one or two parents: `-H -L -A -B`.
3. One bulk with one or two parents: `-H -L` and `-A` or `-B`.

```bash
genome-whisperer GoBSAseq -V cohort.vcf.gz -H ResParent -L SusParent -A ResBulk -B SusBulk \
  -n 20 -m 20 -p F2 -o bsa_out
```
| Flag | Default | |
|---|---|---|
| `-H/-L` | | High (resistant) / low (susceptible) parent names |
| `-A/-B` | | High / low bulk names |
| `--min-{high,low}-parent-depth` | 5 | |
| `--min-{high,low}-bulk-depth` | 40 | |
| `-n/-m` | 20 | Individuals per bulk |
| `-p, --pop-structure` | `F2` | `F2`, `BC` or `RIL` |
| `--rep` | 10000 | Simulation replicates for thresholds |
| `-w, --window-size` / `--window-step` | 2,000,000 / 10,000 | Sliding window |
| `--smooth` | | Smooth the plots |
| `-i, --interactive` | | Prompt for inputs |

#### GeneSpace
Ranks the genes in a QTL window as candidates for a trait. Each gene is scored on two independent axes:

- **Variant score (0–4)**: the most damaging variant whose genotypes follow the resistant/susceptible split. Every SNP is tested under both a dominant and a recessive pattern.
- **Function score (0–3)**: PRG hit quality plus gene-description terms from a built-in trait profile.

The two scores are combined into a tier: **A** lead candidate · **B** strong · **C** plausible · **D** positional only.

```bash
genome-whisperer GeneSpace --gff genes.gff3 -V cohort_SUPER_VCF.tsv \
  --chrom Chr04 --start 1200000 --stop 2400000 \
  --res-lines R1,R2,R3 --sus-lines S1,S2 \
  --prg SQUASHV1.1_PRG_TOP_HITS.txt --gene-description-tsv desc.tsv -o genespace
```
| Flag | Default | |
|---|---|---|
| `-V, --variant` | | superVCF table (from genome-whisperer or GoBSAseq) |
| `--chrom`, `--start`, `-e/--stop` | | QTL window |
| `--res-lines`, `-u/--sus-lines` | | Line names (`R1` or `R1.GT`); unknown names are an error |
| `--trait` | `disease-resistance` | Built-in trait profile |
| `--min-concordance` | 0.9 | Fraction of lines that must fit the pattern |
| `--promoter-length` | 2000 | Upstream bases counted toward a gene (0 disables) |
| `--prg`, `--gene-description-tsv` | | Optional; when absent, that evidence is not scored |

Output: `GeneSpace_<chrom>_<start>-<stop>.xlsx`, with one row per gene, score columns, tier colours and a Legend sheet. See [`gene_space_review.md`](gene_space_review.md) for the scoring rationale.

#### MarkerDesignPipeline
For every QTL (± `--flank-size`), this command:
- slices the superVCF table and names sites from the marker DB,
- keeps sites where resistant and susceptible lines differ,
- adds 50 bp flanks that carry IUPAC codes for every variant in `--gvcf`.

It writes per-QTL TSV and Excel files.
```bash
genome-whisperer MarkerDesignPipeline -V cohort_SUPER_VCF.tsv --qtls qtls.xlsx \
  --res-lines R1,R2 --sus-lines S1,S2 --gvcf joint.vcf.gz --marker-db markers.xlsx \
  --crop Cmos --trait pc --reference genome.fa -o markers
```

---

### Pangenome

#### GoPan
Builds a pangenome by iteratively mapping short reads to a reference and assembling the unmapped reads. Requires `--config`.
```bash
genome-whisperer GoPan --config config.txt --assembler mac   # mac | megahit | masurca
```

---

### Utilities

#### CreateTemplate
Writes a `config.txt` template to the current directory.

#### ScanWarehouse
Inventories sequencing storage and writes an NDJSON report.
```bash
genome-whisperer ScanWarehouse --mounts /mnt/nas1 --mounts /mnt/nas2 --report warehouse.ndjson \
  --validate alignments
```
| Flag | Default | |
|---|---|---|
| `--mounts` | | Roots to scan (repeatable) |
| `--report` | | NDJSON output path |
| `--validate` | `none` | `none`, `alignments` or `all` |
| `--scan-workers` | 16 | Concurrent directory inspections |
| `--emit-batch-seqids` | | Include every sequence ID in batched contig groups |

---

## Typical workflows

**Variants for a panel of lines → candidate genes in a known QTL**
```bash
genome-whisperer PrepareFastaFile --reference genome.fa
genome-whisperer CreateSnpEffDB --config config.txt --species squash --annotation-version v1.1
genome-whisperer GetProtPrgTopHits --config config.txt --prg-db prgdb.fa --species squash --annotation-version v1.1
genome-whisperer AlignReads --config config.txt --bqsr --bootstrap
genome-whisperer CreateSuperVcf --config config.txt --species squash --ref-version v1.1 \
  --database squash_v1.1 --prg SQUASHV1.1_PRG_TOP_HITS.txt
genome-whisperer GeneSpace --gff genes.gff3 -V <..._SUPER_VCF.tsv> --chrom Chr04 \
  --start 1200000 --stop 2400000 --res-lines R1,R2 --sus-lines S1,S2 \
  --prg SQUASHV1.1_PRG_TOP_HITS.txt -o genespace
```

**BSA-seq → QTL → markers**
```bash
genome-whisperer GoBSAseq --config config.txt -H ResParent -L SusParent -A ResBulk -B SusBulk -o bsa
genome-whisperer MarkerDesignPipeline -V <..._SUPER_VCF.tsv> --qtls bsa/<qtls>.xlsx ...
```

---

## Development

```bash
go build ./...
go vet ./...
go test ./...
```

Package layout:

| Package | Contents |
|---|---|
| `cmd/` | Cobra commands and shared flags (`root.go`) |
| `alignment/`, `alignmentdir/` | Read alignment, BQSR, data-directory alignment |
| `variants/` | gVCF creation, merging, hard filtering |
| `annotation/` | snpEff, PRG BLAST, superVCF creation |
| `bsaseq/` | BSA-seq analysis |
| `genespace/` | Candidate gene ranking and trait profiles |
| `pangenome/` | GoPan |
| `warehouse/` | Storage inventory |
| `utils/` | Config parsing, dependency checks, progress bars, logging |

Design notes are in [`docs/architecture.html`](docs/architecture.html).
