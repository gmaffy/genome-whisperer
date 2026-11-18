# GeneSpace Review: Gene Classification & Prioritisation

## 1. What the code does now

**Data collected for each mRNA** in the `chrom:start-stop` window:

| Evidence | Source | How it's used |
|---|---|---|
| `NumSNPs` | variant table | Count only; not used for ranking |
| `NumSigSNPs` | `SNPEFF_EFFECT` not in {UPSTREAM, DOWNSTREAM, SYNONYMOUS_*, INTRON} | Count only; not used for ranking |
| `NumResNotSusSNPs` | All resistant lines are called and differ from one homozygous susceptible genotype | Only as an on/off switch inside tier 4 |
| PRG hit | First BLAST row for that gene in the PRG table | Only `PercIdent > 0` (any hit at all) |
| Description | Free-text TSV | Substring search against about 60 resistance keywords and 50 transcription-factor (TF) keywords |

**Tiers** (`ClassifyPriority`). Each gene gets exactly one tier, and the first matching rule wins:

```
5  TF keyword                              ← checked first, beats everything
4  ResNotSus > 0 AND (PRG OR res keyword)
3  PRG AND res keyword
2  res keyword only
1  PRG only
0  nothing
```

### Main weaknesses

1. **Annotation outweighs genetics.** A gene with a resistance keyword but no segregating variant ranks the same as, or above, genes with real genetic evidence. A TF with no variants at all ranks highest.
2. **The amount of variant evidence is ignored.** One ResNotSus SNP in an intron counts the same as a premature stop codon in all resistant lines.
3. **The PRG hit is accepted at any quality.** A 30% identity hit over 10% of the protein counts. So does an e-value of 1.
4. **Keyword matching on descriptions is unreliable.** You get false hits ("sup**erf**amily" matches `ERF`, "**orga**n" matches `RGA`), descriptions differ between annotation sources, and the keywords only fit disease resistance.
5. **Information is lost.** One tier label hides which evidence was present, so users can't re-sort by what they care about.
6. **No position effect.** A gene at the QTL peak and a gene at the window edge are treated the same.
7. **The genotype rule is all-or-nothing.** One missing call or one mis-phenotyped line drops the SNP entirely.
8. **Rows are per mRNA.** Genes with several isoforms are counted several times.

---

## 2. Proposed approach: separate evidence scores, then a tier

Score each gene on **independent evidence axes**, show each score as its own column, and derive the tier from the combination. That keeps the reasoning visible and lets users filter in Excel.

### Axis 1: Genetic association

Does the gene carry variants whose genotypes follow the phenotype?

- **Concordance score instead of all-or-nothing.** For each SNP: fraction of R/S lines whose genotype matches the expected pattern, with missing calls ignored rather than disqualifying. For example, 9/10 concordant counts as a hit at a threshold you choose.
- **Inheritance mode as a parameter**, because the expected genotype pattern differs:
  - *Dominant*: R lines carry ≥1 alt allele, S lines are homozygous ref (het allowed in R).
  - *Recessive*: R lines are homozygous for an allele that S lines never carry homozygously.
  - *Additive / quantitative*: correlate dosage with a phenotype value.
- **Position:** distance to the QTL peak, or the G′/ΔSNP-index value at the gene if the BSA output is fed in. Genes in the confidence interval core score higher.

### Axis 2: Functional impact of the concordant variants

This uses the impact class already in the variant table (`SNPEff_IMPACT`, written by `annotation/annotation.go`). Only variants that pass Axis 1 count here:

| Concordant variant | Score |
|---|---|
| HIGH (stop gained, frameshift, splice donor/acceptor, start lost) | 4 |
| MODERATE (missense, in-frame indel) | 3 |
| LOW, 5′/3′ UTR, splice region | 2 |
| Upstream ≤2 kb (promoter) | 1 (needs the GFF window extended by 2 kb) |
| Intron / downstream / synonymous | 0 |

Take the maximum per gene, and keep a count column alongside it. Normalise raw SNP counts per kb so that long genes don't float up just because they are long.

### Axis 3: Biological plausibility for the trait

- **Domains over description text.** Run InterProScan/Pfam once per genome and match on domain IDs. For example: NB-ARC PF00931, TIR PF01582, LRR PF00560/PF13855, RPW8 PF05659, LysM PF01476, protein kinase PF00069. This is far more reliable than descriptions.
- **PRG hits with real quality cut-offs**: e.g. identity ≥40%, query coverage (Length/Qlen) ≥60%, e-value ≤1e-10. Also use the PRGdb class (CNL, TNL, RLK, RLP), since CNL/TNL are stronger R-gene candidates than a generic kinase.
- **Keywords as a last resort**, matched on whole words (`\b`), with the lists held in a per-trait profile file rather than in code.
- **TFs as a modifier, not the top tier.** A TF is a plausible regulator for almost any trait, so it adds a point. It shouldn't outrank a gene with direct genetic evidence. A TF that carries a concordant HIGH or MODERATE variant naturally scores well through Axes 1 and 2.

### Axis 4 (optional): Expression

If RNA-seq is available, add a column for differential expression between R and S (ideally after inoculation or stress). This is the strongest independent support for regulatory or promoter candidates.

### Combining the scores into tiers

| Tier | Rule (example) |
|---|---|
| **A: lead candidate** | Concordant variant with impact ≥3 **and** plausibility ≥2 |
| **B: strong** | Concordant variant with impact ≥2 **or** (concordant variant + plausibility ≥2) |
| **C: plausible** | Plausibility ≥2, no concordant functional variant |
| **D: positional only** | In the interval, nothing else |

The rule is that **genetic evidence comes first and annotation adjusts the ranking within that**. Within a tier, sort by total score, then by distance to the peak.

---

## 3. Disease resistance: points that are easy to miss

- **Susceptibility (S) genes.** Recessive resistance often comes from losing a host susceptibility gene. Examples: *mlo* (powdery mildew), *eIF4E/eIFiso4E* (potyviruses), *DMR6*, *SWEET*, *pip* genes. None of these match R-gene keywords. With a recessive mode plus a HIGH-impact variant in the resistant lines, they would surface. They are worth their own keyword/domain list in the resistance profile.
- **Presence/absence variation.** NLRs often exist in clusters that are missing from the reference genome or present only in it. Short reads against one reference can't detect a gene that's absent from the reference. Two cheap signals help: a per-gene coverage ratio between R and S (resistant lines with near-zero coverage suggests a deletion), and a "gene in NLR cluster" flag.
- **Heterozygous resistance in hybrids.** Dominant mode should accept `0/1` in R lines. The current check already does.

---

## 4. Other traits: one engine, different profiles

Axes 1, 2 and 4 don't depend on the trait. Only Axis 3 changes, so a **trait profile** (a YAML or TSV listing domains, gene families, keywords and weights) selected with `--trait` covers all of these:

| Trait | High-plausibility families and domains |
|---|---|
| **Disease resistance** | NLR (NB-ARC, TIR, CC, RPW8), RLK/RLP, LysM, PR proteins, S-genes (MLO, eIF4E, DMR6, SWEET) |
| **Herbicide tolerance** | Target-site enzymes (ALS/AHAS, ACCase, EPSPS, PPO, HPPD, psbA); metabolic: P450s, GSTs, UGTs, ABC transporters |
| **Abiotic stress (drought, heat, salt)** | DREB/CBF, NAC, bZIP-ABF, LEA/dehydrins, HSPs/HSFs, aquaporins, HKT/NHX/SOS transporters, ABA pathway (NCED, PYL, PP2C, SnRK2) |
| **Flowering / maturity** | FT/TFL1 (PEBP), CO/COL, FLC/MADS, GI, PRR/clock genes, photoreceptors (PHY, CRY), SPL |
| **Fruit / flower colour & quality** | R2R3-MYB + bHLH + WD40 (MBW complex), anthocyanin pathway (CHS, CHI, F3H, DFR, ANS), carotenoids (PSY, LCY, CRTISO), SWEET/invertases for sugar |
| **Plant architecture / growth habit** | DELLA/GA pathway, TCP (TB1), auxin (ARF, PIN), CLV/WUS, SP/SFT (determinacy) |
| **Fertility / hybrid systems** | PPR / restorer-of-fertility, mitochondrial ORFs, tapetum TFs |
| **Seed / fruit size, yield** | Cell-number regulators (CNR/FW2.2), GW/GS-type genes, cytokinin oxidases, ubiquitin pathway |

Quantitative traits (yield, size, quality) usually come from **regulatory and dosage changes** more than knockouts. For those profiles, the promoter/UTR score and expression evidence should carry more weight. For qualitative traits (major R-genes, colour, S-genes), coding HIGH/MODERATE variants matter most.

---

## 5. Questions to settle before implementing

1. **Inputs:** is InterProScan/Pfam output, or BSA peak/G′ output, available to feed in? Those decide whether Axis 3 can use domains and whether Axis 1 can include position.
2. **Inheritance mode:** fixed per run, or should a SNP be tested under both dominant and recessive patterns and report which fits?
3. **Output:** keep one sheet with score columns, or add a "Lead candidates" sheet containing tiers A–B?
4. **Profiles:** ship built-in profiles (resistance, herbicide, flowering, …) or only a template file that users fill in?

---

## 6. Implemented (2026-11-18)

Decisions: input is a superVCF table (genome-whisperer or GoBSAseq); every SNP is tested under both dominant and recessive patterns; one data sheet with score columns; built-in profiles, `disease-resistance` only for now.

- **Genes, not mRNAs.** Isoforms are collapsed by `Parent`. Genes are kept if they overlap the window. Malformed GFF attributes are skipped instead of panicking. PRG/description lookups try the gene ID, then each transcript ID.
- **Genotype concordance** (`testSNP`). Every allele, including ref, is tested under DOM and REC. Concordance is the mean of the R and S match fractions. Missing calls are ignored, but a SNP needs ≥1 call per group and ≥50% call rate. `--min-concordance` defaults to 0.9. `1/0` = `0/1`, and `.` is missing.
- **Variant score**: 4 HIGH · 3 MODERATE · 2 UTR/splice region/other LOW · 1 promoter window (`--promoter-length`, default 2000, strand-aware) · 0 otherwise. Synonymous scores 0.
- **Function score**: 3 PRG pass + strong term · 2 either · 1 weak term / weak PRG / TF · 0. PRG pass = identity ≥40%, coverage ≥60%, E ≤1e-10.
- **Tiers** A–D as in section 2. Sort order is tier → variant → function → concordance → # concordant → position (deterministic).
- **Profile** (`genespace/profiles.go`). Word-boundary regexes; S-genes (MLO, eIF4E, DMR6, SWEET) added; generic receptor-like kinases moved to weak; F-box/LEA/ERD/detoxification dropped. `--trait` selects the profile.
- **CLI**: errors exit non-zero; `--prg` and `--gene-description-tsv` are optional; an empty `--out-dir` means `.`; output is `GeneSpace_<chrom>_<start>-<stop>.xlsx`.
- **Tests**: `genespace/genespace_test.go` (unit + end-to-end workbook check).

Not yet done: QTL peak / G′ position scoring, expression evidence, other trait profiles.
