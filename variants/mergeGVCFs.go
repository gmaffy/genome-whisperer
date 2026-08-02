package variants

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// JointVcfDir returns the directory joint-genotyped VCFs belong in. The last
// component names the caller/merger combination that produced them:
// gatk_gatk, gatk_glnexus or dv_glnexus.
//
// Data-dir mode:  <DataDir>/<species>/<refVer>/VCFs/<caller>_<merger>
// Otherwise:      <OutDir>/VCFs/<caller>_<merger>
//
// Per-chromosome VCFs live directly in here with the chromosome in the filename,
// alongside the final concatenated VCF, so there is one place to look.
//
// Giving each combination its own directory means running a second combination
// cannot overwrite the first, the reuse checks never compare a joint VCF against
// gVCFs from a different caller, and the scratch directories (work/, vcfs.list)
// are per-combination too.
func JointVcfDir(opts Options) string {
	if opts.DataDir != "" {
		return filepath.Join(opts.DataDir, strings.ToLower(opts.Species), opts.RefVer, "VCFs", callerMergerTag(opts))
	}
	return filepath.Join(opts.OutDir, "VCFs", callerMergerTag(opts))
}

// callerMergerTag names the caller/merger combination: gatk_gatk, gatk_glnexus
// or dv_glnexus. It appears both as the directory and inside every filename, so
// a VCF copied out of its directory still says how it was produced.
func callerMergerTag(opts Options) string {
	// "dv" rather than "deepvariant", to match the dv_gvcfs directory name.
	caller := "gatk"
	if strings.ToLower(opts.Caller) != "gatk" {
		caller = "dv"
	}
	return caller + "_" + strings.ToLower(opts.Merger)
}

// JointVcfPath returns the joint VCF path for one chromosome:
//
//	<SPECIES>.<refver>.<caller>_<merger>.<chrom>.joint.vcf.gz
//
// The combination sits before the chromosome so the trailing
// ".<chrom>.joint.vcf.gz" stays intact — that tail is what suffix and glob
// checks match on.
func JointVcfPath(opts Options, chrom string) string {
	label := strings.ReplaceAll(chrom, ".", "_")
	name := fmt.Sprintf("%s.%s.%s.%s.joint.vcf.gz",
		strings.ToUpper(opts.Species), strings.ToLower(opts.RefVer), callerMergerTag(opts), label)
	return filepath.Join(JointVcfDir(opts), name)
}

// FinalVcfPath returns the concatenated multi-sample VCF path:
//
//	<SPECIES>.<refver>.<caller>_<merger>.all.vcf.gz
//
// Same ordering as JointVcfPath, so ".all.vcf.gz" remains the tail.
func FinalVcfPath(opts Options) string {
	name := fmt.Sprintf("%s.%s.%s.all.vcf.gz",
		strings.ToUpper(opts.Species), strings.ToLower(opts.RefVer), callerMergerTag(opts))
	return filepath.Join(JointVcfDir(opts), name)
}

// vcfSampleNames reads the sample columns from a VCF or gVCF header.
func vcfSampleNames(vcf string) ([]string, error) {
	in, cleanup, err := openVCF(vcf)
	if err != nil {
		return nil, fmt.Errorf("open %q: %w", vcf, err)
	}
	defer cleanup()

	scanner := bufio.NewScanner(in)
	// Header lines get long with many samples.
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) <= 9 {
				return nil, fmt.Errorf("%s: no sample columns", vcf)
			}
			return fields[9:], nil
		}
		if !strings.HasPrefix(line, "#") {
			break
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanning %s: %w", vcf, err)
	}
	return nil, fmt.Errorf("%s: no #CHROM header line found", vcf)
}

// sampleNamesMatch compares two sample sets ignoring order.
//
// GATK and GLnexus both emit sorted sample columns while the gVCFs arrive in
// discovery order, so an order-sensitive comparison reports a mismatch on almost
// every run and forces a pointless re-merge.
func sampleNamesMatch(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	x := append([]string(nil), a...)
	y := append([]string(nil), b...)
	sort.Strings(x)
	sort.Strings(y)
	for i := range x {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}

// dictOrder maps each sequence ID to its position in the reference dictionary.
//
// getChromsAndContigs sorts by length so it cannot be used for this: gatk
// MergeVcfs needs its inputs in dictionary order, and feeding it length order
// yields a mis-sorted VCF.
func dictOrder(dictFilePath string) (map[string]int, error) {
	f, err := os.Open(dictFilePath)
	if err != nil {
		return nil, fmt.Errorf("opening reference dict file %s: %w", dictFilePath, err)
	}
	defer f.Close()

	order := make(map[string]int)
	pos := 0
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "@SQ") {
			continue
		}
		for _, part := range strings.Split(line, "\t") {
			if strings.HasPrefix(part, "SN:") {
				order[strings.TrimPrefix(part, "SN:")] = pos
				pos++
				break
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanning dict file: %w", err)
	}
	return order, nil
}

// writeSampleMap writes the GenomicsDBImport sample map (sample<TAB>path).
//
// Sample names come from each gVCF's own header rather than from its filename,
// because that is what GenomicsDBImport keys on. A gVCF with no sample column or
// a duplicated sample name is an error here instead of a confusing failure
// several minutes into the import.
func writeSampleMap(path string, gvcfs []string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("creating sample map %s: %w", path, err)
	}
	defer f.Close()

	seen := make(map[string]string, len(gvcfs))
	for _, gvcf := range gvcfs {
		names, nErr := vcfSampleNames(gvcf)
		if nErr != nil {
			return fmt.Errorf("reading sample name from %s: %w", gvcf, nErr)
		}
		if len(names) != 1 {
			return fmt.Errorf("%s: expected 1 sample, got %d", gvcf, len(names))
		}
		if prev, dup := seen[names[0]]; dup {
			return fmt.Errorf("duplicate sample %q in %s and %s", names[0], prev, gvcf)
		}
		seen[names[0]] = gvcf
		if _, wErr := fmt.Fprintf(f, "%s\t%s\n", names[0], gvcf); wErr != nil {
			return fmt.Errorf("writing sample map: %w", wErr)
		}
	}
	return nil
}

// FindExistingGvcfs inventories the gVCFs already on disk, grouped by chromosome,
// so MergeGvcfs can run standalone without CreateGvcfs having run in-process.
// It uses GvcfPath, so it looks exactly where CreateGvcfs writes.
func FindExistingGvcfs(opts Options) (map[string][]string, error) {
	dictFilePath := opts.RefFasta[:len(opts.RefFasta)-len(filepath.Ext(opts.RefFasta))] + ".dict"
	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return nil, fmt.Errorf("getting chromosomes and contigs: %w", err)
	}

	samples, skipped, err := FindGvcfSamples(opts)
	if err != nil {
		return nil, err
	}
	if len(skipped) > 0 {
		color.Yellow("Skipped %d sample(s) with no usable alignment: %v\n", len(skipped), skipped)
	}
	if len(samples) == 0 {
		return nil, fmt.Errorf("no usable samples found")
	}

	labels := make([]string, 0, len(chroms)+1)
	for _, c := range chroms {
		labels = append(labels, c.ID)
	}
	if len(contigs) > 0 {
		labels = append(labels, "contigs")
	}

	gvcfs := make(map[string][]string)
	for _, s := range samples {
		for _, label := range labels {
			theGVCF := GvcfPath(opts, s, label)
			if _, sErr := os.Stat(theGVCF); sErr != nil {
				color.Red("[%s] gVCF missing for %s: %s\n", s.Sample, label, theGVCF)
				continue
			}
			if !opts.SkipVerification {
				if vErr := utils.ValidateGvcf(theGVCF, opts.Verbose, opts.Quick); vErr != nil {
					color.Red("[%s] gVCF for %s is corrupt: %v\n", s.Sample, label, vErr)
					continue
				}
			}
			gvcfs[label] = append(gvcfs[label], theGVCF)
		}
	}

	for label := range gvcfs {
		sort.Strings(gvcfs[label])
	}
	return gvcfs, nil
}

// gatkIntervalArgs renders the -L value for a group of sequences, plus any extra
// GenomicsDBImport arguments the group needs.
//
// A single chromosome goes straight into -L. The contigs group covers many
// sequences and "contigs" is not a sequence name, so passing the label would make
// GenomicsDBImport fail with "contig contigs not present in the sequence
// dictionary"; an interval list is written into dir instead, the same way
// CreateGvcfGATK does. --merge-input-intervals keeps hundreds of small contigs
// from becoming hundreds of GenomicsDB partitions.
//
// Separate from mergeChromGATK so the interval handling can be tested without
// GATK installed.
func gatkIntervalArgs(dir string, seqs []SeqInfo) (regionArg, extraArgs string, err error) {
	if len(seqs) == 0 {
		return "", "", fmt.Errorf("no sequences given")
	}
	if len(seqs) == 1 {
		return seqs[0].ID, "", nil
	}

	intervals := filepath.Join(dir, "intervals.list")
	f, cErr := os.Create(intervals)
	if cErr != nil {
		return "", "", fmt.Errorf("creating %s: %w", intervals, cErr)
	}
	for _, s := range seqs {
		if _, wErr := fmt.Fprintln(f, s.ID); wErr != nil {
			f.Close()
			return "", "", fmt.Errorf("writing %s: %w", intervals, wErr)
		}
	}
	if clErr := f.Close(); clErr != nil {
		return "", "", fmt.Errorf("closing %s: %w", intervals, clErr)
	}
	return intervals, " --merge-input-intervals", nil
}

// mergeChromGATK joint-genotypes one chromosome with GenomicsDBImport followed by
// GenotypeGVCFs.
func mergeChromGATK(opts Options, chrom string, seqs []SeqInfo, gvcfs []string, jointVCF string) error {
	workDir := filepath.Join(JointVcfDir(opts), "work", strings.ReplaceAll(chrom, ".", "_"))
	theDB := filepath.Join(workDir, "db")
	tmpDir := filepath.Join(workDir, "tmp")

	if err := os.MkdirAll(tmpDir, 0755); err != nil {
		return fmt.Errorf("creating %s: %w", tmpDir, err)
	}
	// GenomicsDBImport refuses to write into an existing workspace.
	if err := os.RemoveAll(theDB); err != nil {
		return fmt.Errorf("removing %s: %w", theDB, err)
	}

	sampleMap := filepath.Join(workDir, "sample_map.txt")
	if err := writeSampleMap(sampleMap, gvcfs); err != nil {
		return err
	}

	regionArg, extraArgs, err := gatkIntervalArgs(workDir, seqs)
	if err != nil {
		return fmt.Errorf("%s: %w", chrom, err)
	}

	gDBCmd := fmt.Sprintf(
		`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport --sample-name-map %s `+
			`--genomicsdb-workspace-path %s --tmp-dir %s -L %s%s `+
			`--genomicsdb-shared-posixfs-optimizations true --batch-size 50 `+
			`--bypass-feature-reader --verbosity %s`,
		sampleMap, theDB, tmpDir, regionArg, extraArgs, opts.GatkLogLevel,
	)
	if err := utils.RunCmd(gDBCmd, opts.Verbose); err != nil {
		return fmt.Errorf("gatk GenomicsDBImport: %w", err)
	}

	genoCmd := fmt.Sprintf(
		`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
		opts.RefFasta, theDB, jointVCF, tmpDir, opts.GatkLogLevel,
	)
	if err := utils.RunCmd(genoCmd, opts.Verbose); err != nil {
		return fmt.Errorf("gatk GenotypeGVCFs: %w", err)
	}

	// The workspace and temp files are large and are not needed once the joint
	// VCF exists.
	if rErr := os.RemoveAll(workDir); rErr != nil {
		color.Yellow("could not clean up %s: %v\n", workDir, rErr)
	}
	return nil
}

// mergeChromGlnexus joint-genotypes one chromosome with glnexus_cli, then
// converts the BCF to a bgzipped, indexed VCF.
func mergeChromGlnexus(opts Options, chrom string, gvcfs []string, jointVCF string) error {
	workDir := filepath.Join(JointVcfDir(opts), "work", strings.ReplaceAll(chrom, ".", "_"))
	if err := os.MkdirAll(workDir, 0755); err != nil {
		return fmt.Errorf("creating %s: %w", workDir, err)
	}
	// glnexus_cli refuses to write into an existing database directory.
	glnexusDB := filepath.Join(workDir, "GLnexus.DB")
	if err := os.RemoveAll(glnexusDB); err != nil {
		return fmt.Errorf("removing %s: %w", glnexusDB, err)
	}

	preset := "DeepVariant"
	if strings.ToLower(opts.Caller) == "gatk" {
		preset = "gatk"
	}

	jointBCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".bcf"
	glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s %s > %s`,
		preset, glnexusDB, strings.Join(gvcfs, " "), jointBCF)
	if err := utils.RunCmd(glnexusCmd, opts.Verbose); err != nil {
		return fmt.Errorf("glnexus_cli: %w", err)
	}

	bcfCmd := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
	if err := utils.RunCmd(bcfCmd, opts.Verbose); err != nil {
		return fmt.Errorf("bcftools view: %w", err)
	}

	if err := utils.RunCmd(fmt.Sprintf(`tabix -f -p vcf %s`, jointVCF), opts.Verbose); err != nil {
		return fmt.Errorf("indexing %s: %w", jointVCF, err)
	}

	os.Remove(jointBCF)
	if rErr := os.RemoveAll(workDir); rErr != nil {
		color.Yellow("could not clean up %s: %v\n", workDir, rErr)
	}
	return nil
}

// MergeGvcfs joint-genotypes each chromosome and concatenates the results into a
// single multi-sample VCF, whose path it returns.
//
// gvcfs is the map CreateGvcfs returns. Pass nil to have MergeGvcfs discover the
// gVCFs itself, which is how the standalone MergeGvcfs subcommand works.
//
// A chromosome whose sample set is short of the largest seen is reported and
// skipped rather than aborting the run, so one bad sample no longer blocks the
// whole cohort. An existing joint VCF is reused when it is valid and already
// holds exactly the expected samples.
func MergeGvcfs(opts Options, gvcfs map[string][]string) (string, error) {

	// ==================================== Validate inputs ===================================== //

	merger := strings.ToLower(opts.Merger)
	if merger != "gatk" && merger != "glnexus" {
		return "", fmt.Errorf("merger must be gatk or glnexus, got %q", opts.Merger)
	}
	if strings.ToLower(opts.Caller) == "deepvariant" && merger != "glnexus" {
		return "", fmt.Errorf("deepvariant gVCFs must be merged with glnexus, not %q", opts.Merger)
	}
	if opts.RefFasta == "" {
		return "", fmt.Errorf("reference fasta must not be empty")
	}

	// Same normalisation as CreateGvcfs, so a standalone MergeGvcfs looks in the
	// same absolute locations the create stage wrote to.
	var err error
	if opts, err = absRoots(opts); err != nil {
		return "", err
	}

	dictFilePath := opts.RefFasta[:len(opts.RefFasta)-len(filepath.Ext(opts.RefFasta))] + ".dict"
	if _, dErr := os.Stat(dictFilePath); dErr != nil {
		return "", fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	if gvcfs == nil {
		color.Cyan("================================== Finding gVCFs ==================================\n\n")
		if gvcfs, err = FindExistingGvcfs(opts); err != nil {
			return "", err
		}
	}
	if len(gvcfs) == 0 {
		return "", fmt.Errorf("no gVCFs to merge")
	}

	vcfDir := JointVcfDir(opts)
	if err := os.MkdirAll(vcfDir, 0755); err != nil {
		return "", fmt.Errorf("creating %s: %w", vcfDir, err)
	}

	// ============================= Decide which chromosomes to merge =========================== //

	// The full cohort size is the largest sample set seen across chromosomes.
	// Anything short of it is missing samples and is not safe to joint-call.
	expected := 0
	for _, paths := range gvcfs {
		if len(paths) > expected {
			expected = len(paths)
		}
	}

	order, err := dictOrder(dictFilePath)
	if err != nil {
		return "", err
	}

	// The sequences behind each label. GATK needs the real sequence list, not the
	// label: "contigs" stands for every small sequence and is not itself a
	// sequence name.
	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return "", fmt.Errorf("getting chromosomes and contigs: %w", err)
	}
	seqsFor := make(map[string][]SeqInfo, len(chroms)+1)
	for _, c := range chroms {
		seqsFor[c.ID] = []SeqInfo{c}
	}
	if len(contigs) > 0 {
		seqsFor["contigs"] = contigs
	}

	labels := make([]string, 0, len(gvcfs))
	for label, paths := range gvcfs {
		if len(paths) < expected {
			color.Yellow("[%s] skipping: only %d/%d samples have a gVCF\n", label, len(paths), expected)
			continue
		}
		labels = append(labels, label)
	}
	if len(labels) == 0 {
		return "", fmt.Errorf("no chromosome has a complete set of %d gVCFs", expected)
	}

	// Dictionary order, contigs last. gatk MergeVcfs needs its inputs in this
	// order; the previous code appended them in whatever order workers finished.
	sort.Slice(labels, func(i, j int) bool {
		pi, oki := order[labels[i]]
		pj, okj := order[labels[j]]
		if !oki {
			pi = len(order) // "contigs"
		}
		if !okj {
			pj = len(order)
		}
		return pi < pj
	})

	// ======================================= Merge ============================================ //

	color.Cyan("=========================== Merging %d chromosome group(s) with %s ===========================\n\n",
		len(labels), merger)

	var jointVcfs []string
	for _, label := range labels {
		jointVCF := JointVcfPath(opts, label)

		// ------------------------- reuse a valid, up-to-date joint VCF ------------------------- //
		if _, sErr := os.Stat(jointVCF); sErr == nil {
			reuse := true

			// --skip-verification skips the integrity check only.
			if !opts.SkipVerification {
				if vErr := utils.ValidateGvcf(jointVCF, opts.Verbose, opts.Quick); vErr != nil {
					color.Yellow("[%s] joint VCF is corrupt, re-merging: %v\n", label, vErr)
					reuse = false
				}
			}

			// The sample-set check always runs, even under --skip-verification: it is
			// only a header read, and skipping it means a newly added sample silently
			// never reaches the output.
			if reuse {
				want, wErr := allGvcfSampleNames(gvcfs[label])
				have, hErr := vcfSampleNames(jointVCF)
				switch {
				case wErr != nil || hErr != nil:
					color.Yellow("[%s] could not compare sample names, re-merging\n", label)
					reuse = false
				case !sampleNamesMatch(want, have):
					color.Yellow("[%s] sample set changed (%d gVCFs vs %d samples in the joint VCF), re-merging\n",
						label, len(want), len(have))
					reuse = false
				}
			}

			if reuse {
				color.Green("[%s] joint VCF is up to date, reusing: %s\n", label, jointVCF)
				jointVcfs = append(jointVcfs, jointVCF)
				continue
			}
			os.Remove(jointVCF)
			os.Remove(jointVCF + ".tbi")
		}

		// An explicit --gvcf list arrives under a synthetic label ("all") that is not
		// a sequence name, so fall back to the whole reference for it.
		seqs := seqsFor[label]
		if len(seqs) == 0 {
			seqs = append(append([]SeqInfo{}, chroms...), contigs...)
		}

		color.Cyan("[%s] merging %d gVCFs ...\n\n", label, len(gvcfs[label]))
		var mErr error
		if merger == "gatk" {
			mErr = mergeChromGATK(opts, label, seqs, gvcfs[label], jointVCF)
		} else {
			mErr = mergeChromGlnexus(opts, label, gvcfs[label], jointVCF)
		}
		if mErr != nil {
			color.Red("[%s] merge FAILED: %v\n", label, mErr)
			continue
		}
		color.Green("[%s] joint VCF created: %s\n\n", label, jointVCF)
		jointVcfs = append(jointVcfs, jointVCF)
	}

	if len(jointVcfs) == 0 {
		return "", fmt.Errorf("no chromosome merged successfully")
	}
	if len(jointVcfs) != len(labels) {
		color.Yellow("%d of %d chromosome group(s) merged; the concatenated VCF will be incomplete\n",
			len(jointVcfs), len(labels))
	}

	// ==================================== Concatenate ========================================= //

	color.Cyan("=========================== Concatenating %d joint VCF(s) ===========================\n\n", len(jointVcfs))

	finalVCF, err := ConcatenateVcfs(jointVcfs, FinalVcfPath(opts), opts.Verbose)
	if err != nil {
		return "", fmt.Errorf("concatenating VCFs: %w", err)
	}

	if _, iErr := os.Stat(finalVCF + ".tbi"); iErr != nil {
		if tErr := utils.RunCmd(fmt.Sprintf(`tabix -f -p vcf %s`, finalVCF), opts.Verbose); tErr != nil {
			return "", fmt.Errorf("indexing %s: %w", finalVCF, tErr)
		}
	}

	color.Green("Multi-sample VCF: %s\n\n", finalVCF)
	return finalVCF, nil
}

// allGvcfSampleNames collects the sample name from each gVCF in order.
func allGvcfSampleNames(gvcfs []string) ([]string, error) {
	var all []string
	for _, gvcf := range gvcfs {
		names, err := vcfSampleNames(gvcf)
		if err != nil {
			return nil, err
		}
		all = append(all, names...)
	}
	return all, nil
}

// ---------------------------------------------------------------------------
// Absorbed from the retired RunVariantCaller.go / RunVariantCallerDir.go
// ---------------------------------------------------------------------------

// ConcatenateVcfs writes vcfs into outVCF with gatk MergeVcfs, copying instead
// when there is only one input. vcfs must already be in reference-dictionary
// order; MergeVcfs will not reorder them.
//
// The destination is passed in rather than derived here, so every VCF name in the
// package comes from JointVcfPath or FinalVcfPath and cannot drift.
func ConcatenateVcfs(vcfs []string, outVCF string, verbose bool) (string, error) {
	if len(vcfs) == 0 {
		return "", fmt.Errorf("no VCFs to concatenate")
	}

	if len(vcfs) == 1 {
		if err := utils.CopyFile(vcfs[0], outVCF); err != nil {
			return "", fmt.Errorf("copying %s to %s: %w", vcfs[0], outVCF, err)
		}
		if _, statErr := os.Stat(vcfs[0] + ".tbi"); statErr == nil {
			if cErr := utils.CopyFile(vcfs[0]+".tbi", outVCF+".tbi"); cErr != nil {
				return "", fmt.Errorf("copying index %s.tbi to %s.tbi: %w", vcfs[0], outVCF, cErr)
			}
		}
		return outVCF, nil
	}

	// The input list lives next to the output, so parallel runs of different
	// caller/merger combinations cannot overwrite each other's list.
	vcfListPath := strings.TrimSuffix(outVCF, ".vcf.gz") + ".vcfs.list"
	f, err := os.Create(vcfListPath)
	if err != nil {
		return "", fmt.Errorf("creating %s: %w", vcfListPath, err)
	}
	for _, vcf := range vcfs {
		if _, wErr := fmt.Fprintf(f, "%s\n", vcf); wErr != nil {
			f.Close()
			return "", fmt.Errorf("writing %s: %w", vcfListPath, wErr)
		}
	}
	if cErr := f.Close(); cErr != nil {
		return "", fmt.Errorf("closing %s: %w", vcfListPath, cErr)
	}

	cmd := fmt.Sprintf(`gatk MergeVcfs -I %s -O %s`, vcfListPath, outVCF)
	if err := utils.RunCmd(cmd, verbose); err != nil {
		return "", fmt.Errorf("gatk MergeVcfs: %w", err)
	}
	return outVCF, nil
}
