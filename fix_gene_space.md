# Plan: Fix and Optimize `genespace/gene_space.go`

## Objective
Fix critical bugs, improve performance, and add output functionality to the `GeneSpace` tool.

## Problems Identified
1. **Compilation Error**: `prg` is a slice in `genePrgMap`, but accessed as a struct in `GeneSpace`.
2. **Genotype Lookup**: `resLines` and `susLines` names from CLI don't match the `.GT` column names in the VCF table.
3. **Brittle GFF Parsing**: Attribute parsing is too simple and only looks for `mRNA` features by start position.
4. **Performance**: O(N*M) nested loops for genes and SNPs.
5. **Logic Gaps**: `NumResNotSusSNPs` field is never populated.
6. **No Output**: The results are calculated but never saved to a file.

## Proposed Changes

### 1. Robust GFF Parsing
- Update `genePosFromGffMap` to properly parse attributes and find the `ID` or `Name`.
- Consider both `gene` and `mRNA` features (user choice, but usually `mRNA` is used for coordinates).
- Ensure range checking handles gene spans properly.

### 2. Efficient SNP Lookup
- Index `qtlRecords` by position in a map `map[int][]VCFRecord` for O(1) lookup during gene iteration.

### 3. Correct Genotype Mapping
- In `GeneSpace`, map sample names in `resLines` and `susLines` to their corresponding `.GT` column names found in the header.

### 4. Fix PRG Selection
- Update `genePrgMap` to return a `map[string]PrgMatches` by picking the best match (e.g., highest `PercIdent`) for each gene.

### 5. Add Excel Output
- Add `github.com/xuri/excelize/v2` dependency to `go.mod`.
- Implement a `writeExcel` function (inspired by `geneSpace.go.1`) to save the `GeneSpaceRow` slice to an `.xlsx` file.

### 6. Logic Completion
- Ensure all fields in `GeneSpaceRow` are correctly calculated, including `NumResNotSusSNPs`.

## Implementation Steps

### Step 1: Improve Data Structures and Parsers
- Update `PrgMatches` to handle single best match per gene.
- Refactor `genePrgMap` to pick the best match.
- Refactor `genePosFromGffMap` for better attribute extraction.

### Step 2: Optimize `GeneSpace` Function
- Load VCF records and index them by position.
- Fix the genotype column lookup logic.
- Populate all `GeneSpaceRow` fields.

### Step 3: Implement Output
- Add an Excel writer to save results using `excelize`.
- Update `cmd/GeneSpace.go` to handle the output file path (optional, or use a default name).

## Verification Plan
1. **Compilation**: Ensure the package compiles without errors.
2. **Manual Test**: Run `genome-whisperer GeneSpace` with sample data and verify the output `.xlsx` contents.
3. **Performance Check**: Verify that processing large regions is significantly faster with the indexed lookup.
