package variants

import (
	"bufio"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// JointVcfDir returns the directory joint-genotyped VCFs belong in. The last
// component names the caller/merger combination that produced them:
// gatk_gatk, gatk_glnexus or dv_glnexus.
//
// Data-dir mode:  <DataDir>/<species>/MERGED_VCFs/<refVer>/<caller>_<merger>
// Otherwise:      <OutDir>/VCFs/<caller>_<merger>
//
// Per-chromosome VCFs live directly in here with the chromosome in the filename,
// alongside the final concatenated VCF, so there is one place to look.
//
// Giving each combination its own directory means running a second combination
// cannot overwrite the first, the reuse checks never compare a joint VCF against
// gVCFs from a different caller, and the concatenation list (vcfs.list) is
// per-combination too. The scratch work directory (GenomicsDB workspace,
// GLnexus database) is not under here — see mergeWorkDir — since it does
// small-file I/O unsuited to whatever storage this directory lives on.
func JointVcfDir(opts Options) string {
	if opts.DataDir != "" {
		return filepath.Join(opts.DataDir, strings.ToLower(opts.Species), "MERGED_VCFs", opts.RefVer, callerMergerTag(opts))
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
// It uses GvcfPath, so it looks exactly where CreateGvcfs writes. The returned
// int is the number of samples discovered — the true cohort size MergeGvcfs
// should hold every chromosome to, not an estimate inferred from the gVCFs
// found.
func FindExistingGvcfs(opts Options) (map[string][]string, int, error) {
	dictFilePath := utils.DictPath(opts.RefFasta)
	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return nil, 0, fmt.Errorf("getting chromosomes and contigs: %w", err)
	}

	samples, skipped, err := FindSampleAlignments(opts)
	if err != nil {
		return nil, 0, err
	}
	if len(skipped) > 0 {
		color.Yellow("Skipped %d sample(s) with no usable alignment: %v\n", len(skipped), skipped)
	}
	if len(samples) == 0 {
		return nil, 0, fmt.Errorf("no usable samples found")
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
	return gvcfs, len(samples), nil
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

// gdbImportDoneFile is written next to the GenomicsDB workspace once
// GenomicsDBImport has returned successfully, and holds the fingerprint of the
// inputs that were imported. It sits beside the workspace rather than inside it
// because GenomicsDBImport insists on creating the workspace directory itself.
const gdbImportDoneFile = "db.import.done"

// gdbFingerprint identifies the exact inputs a GenomicsDB workspace was built
// from: every gVCF (path, size and modification time) and the intervals that
// were imported. A workspace whose fingerprint no longer matches was built from
// a different cohort, or from gVCFs that have since been re-created, so
// genotyping it would produce a VCF that does not match the gVCFs on disk.
func gdbFingerprint(gvcfs []string, seqs []SeqInfo) (string, error) {
	h := sha256.New()
	for _, gvcf := range gvcfs {
		info, err := os.Stat(gvcf)
		if err != nil {
			return "", fmt.Errorf("stat %s: %w", gvcf, err)
		}
		fmt.Fprintf(h, "gvcf\t%s\t%d\t%d\n", gvcf, info.Size(), info.ModTime().UnixNano())
	}
	for _, s := range seqs {
		fmt.Fprintf(h, "seq\t%s\n", s.ID)
	}
	return hex.EncodeToString(h.Sum(nil)), nil
}

// gdbReuse decides whether an existing GenomicsDB workspace can be genotyped as
// it stands. discardReason is non-empty only when there is a workspace that has
// to be thrown away, so a first run says nothing.
//
// The marker file is the authority: mergeChromGATK writes it only after
// GenomicsDBImport returns successfully, so its absence means the import was
// killed part-way and the workspace holds an unknown subset of the cohort. The
// metadata check is a second, cheaper line of defence against a workspace that
// was truncated or partly deleted after the import finished. Both fail towards
// re-importing, which costs time but cannot produce a wrong VCF.
func gdbReuse(theDB, marker, fingerprint string) (reuse bool, discardReason string) {
	if _, err := os.Stat(theDB); err != nil {
		return false, ""
	}
	recorded, err := os.ReadFile(marker)
	if err != nil {
		return false, "the import never completed"
	}
	if strings.TrimSpace(string(recorded)) != fingerprint {
		return false, "it was built from different gVCFs or intervals"
	}
	for _, meta := range []string{"callset.json", "vidmap.json"} {
		if _, sErr := os.Stat(filepath.Join(theDB, meta)); sErr != nil {
			return false, "its " + meta + " is missing"
		}
	}
	return true, ""
}

// mergeWorkDir returns the local scratch directory for one chromosome group's
// GenomicsDB workspace or GLnexus database. Unlike JointVcfDir, this is rooted
// at opts.LocalWorkDir (default $HOME/.genome-whisperer/merge-work), not
// DataDir/OutDir, because GenomicsDBImport and glnexus_cli both do heavy
// small-file, random I/O that performs badly on NAS/NFS storage — the final
// joint VCF is what belongs on shared storage, not the workspace that builds
// it.
func mergeWorkDir(opts Options, chrom string) string {
	return filepath.Join(opts.LocalWorkDir, strings.ToLower(opts.Species), opts.RefVer,
		callerMergerTag(opts), "work", strings.ReplaceAll(chrom, ".", "_"))
}

// gatkScratchDir returns the directory passed as GATK's own --tmp-dir. It is
// kept separate from mergeWorkDir, rooted at utils.TmpBase() rather than at
// opts.LocalWorkDir, so it lives on the machine's scratch disk and is cleaned up
// independently of the GenomicsDB workspace.
//
// Unlike the alignment pipeline's large spills this stays on the default base
// even when that is a RAM-backed tmpfs, because GenotypeGVCFs streams and the
// import writes its bulk into the workspace, not here. What GenomicsDBImport
// does put here still grows with cohort size and interval length, so a cohort
// too large for the default base is moved with $GENOME_WHISPERER_TMPDIR, which
// utils.TmpBase() reads.
func gatkScratchDir(opts Options, chrom string) string {
	return filepath.Join(utils.TmpBase(), "genome-whisperer", "merge",
		strings.ToLower(opts.Species), opts.RefVer, callerMergerTag(opts),
		strings.ReplaceAll(chrom, ".", "_"))
}

// mergeChromGATK joint-genotypes one chromosome with GenomicsDBImport followed by
// GenotypeGVCFs.
//
// The import is the expensive half and is skipped when a previous run already
// completed it for exactly these gVCFs — which is what a run interrupted between
// the two steps, or one whose GenotypeGVCFs failed, leaves behind. An incomplete
// or stale workspace is deleted and rebuilt.
//
// The GenomicsDB workspace (workDir/theDB) lives under opts.LocalWorkDir and
// the GATK --tmp-dir scratch under the OS temp dir — both local, not the NAS
// the final joint VCF is written to (JointVcfDir) — since GenomicsDBImport's
// small-file random I/O performs badly over the network.
//
// memGB is this job's share of the machine; both JVM heaps are sized from it, so
// running several chromosomes at once shrinks each one rather than
// oversubscribing memory.
func mergeChromGATK(opts Options, chrom string, seqs []SeqInfo, gvcfs []string, jointVCF, tag string, memGB int) error {
	workDir := mergeWorkDir(opts, chrom)
	theDB := filepath.Join(workDir, "db")
	tmpDir := gatkScratchDir(opts, chrom)
	marker := filepath.Join(workDir, gdbImportDoneFile)

	if err := os.MkdirAll(workDir, 0755); err != nil {
		return fmt.Errorf("creating %s: %w", workDir, err)
	}
	if err := os.MkdirAll(tmpDir, 0755); err != nil {
		return fmt.Errorf("creating %s: %w", tmpDir, err)
	}
	// Every return below this point has to take the scratch with it. Unlike
	// workDir, which an interrupted import resumes from, nothing in here is
	// reusable — and on a tmpfs base a failed chromosome would otherwise hold
	// its scratch as memory until the machine is rebooted.
	defer func() {
		if rErr := os.RemoveAll(tmpDir); rErr != nil {
			color.Yellow("could not clean up %s: %v\n", tmpDir, rErr)
		}
	}()

	fingerprint, err := gdbFingerprint(gvcfs, seqs)
	if err != nil {
		return err
	}

	importHeapGB, genotypeHeapGB := gatkMergeHeaps(memGB)

	reuseDB, discardReason := gdbReuse(theDB, marker, fingerprint)
	if reuseDB {
		color.Green("%s[%s] GenomicsDB workspace is complete, skipping the import\n", tag, chrom)
	} else {
		if discardReason != "" {
			color.Yellow("%s[%s] discarding the GenomicsDB workspace: %s\n", tag, chrom, discardReason)
		}
		// The marker goes first: if the removal below is interrupted, what is left
		// must not look importable. GenomicsDBImport also refuses to write into an
		// existing workspace, so the directory has to go either way.
		if rErr := os.Remove(marker); rErr != nil && !os.IsNotExist(rErr) {
			return fmt.Errorf("removing %s: %w", marker, rErr)
		}
		if rErr := os.RemoveAll(theDB); rErr != nil {
			return fmt.Errorf("removing %s: %w", theDB, rErr)
		}

		sampleMap := filepath.Join(workDir, "sample_map.txt")
		if sErr := writeSampleMap(sampleMap, gvcfs); sErr != nil {
			return sErr
		}

		regionArg, extraArgs, iErr := gatkIntervalArgs(workDir, seqs)
		if iErr != nil {
			return fmt.Errorf("%s: %w", chrom, iErr)
		}

		gDBCmd := fmt.Sprintf(
			`gatk --java-options "-Xmx%dg -Xms%dg" GenomicsDBImport --sample-name-map %s `+
				`--genomicsdb-workspace-path %s --tmp-dir %s -L %s%s `+
				`--genomicsdb-shared-posixfs-optimizations true --batch-size 50 `+
				`--bypass-feature-reader --verbosity %s`,
			importHeapGB, importHeapGB, sampleMap, theDB, tmpDir, regionArg, extraArgs, opts.GatkLogLevel,
		)
		if rErr := utils.RunCmd(gDBCmd, opts.Verbose); rErr != nil {
			return fmt.Errorf("gatk GenomicsDBImport: %w", rErr)
		}

		if wErr := os.WriteFile(marker, []byte(fingerprint), 0644); wErr != nil {
			return fmt.Errorf("writing %s: %w", marker, wErr)
		}
	}

	genoCmd := fmt.Sprintf(
		`gatk --java-options "-Xmx%dg" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
		genotypeHeapGB, opts.RefFasta, theDB, jointVCF, tmpDir, opts.GatkLogLevel,
	)
	if err := utils.RunCmd(genoCmd, opts.Verbose); err != nil {
		return fmt.Errorf("gatk GenotypeGVCFs: %w", err)
	}

	// The workspace is large and is not needed once the joint VCF exists. The
	// scratch directory is handled by the deferred cleanup above, which also
	// covers the failure paths.
	if rErr := os.RemoveAll(workDir); rErr != nil {
		color.Yellow("could not clean up %s: %v\n", workDir, rErr)
	}
	return nil
}

// mergeChromGlnexus joint-genotypes one chromosome with glnexus_cli, then
// converts the BCF to a bgzipped, indexed VCF.
//
// memGB and threads are this job's share of the machine. Left to itself
// glnexus_cli claims most of the system's memory and every core, which is
// correct for one job and ruinous for several running side by side, so both are
// passed explicitly.
//
// The GLnexus database lives under opts.LocalWorkDir, not the NAS JointVcfDir
// is on — like GenomicsDB, it does heavy small-file random I/O that performs
// badly over the network.
func mergeChromGlnexus(opts Options, chrom string, gvcfs []string, jointVCF, tag string, memGB, threads int) error {
	workDir := mergeWorkDir(opts, chrom)
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
	glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s --threads %d --mem-gbytes %d %s > %s`,
		preset, glnexusDB, threads, memGB, strings.Join(gvcfs, " "), jointBCF)
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

	dictFilePath := utils.DictPath(opts.RefFasta)
	if _, dErr := os.Stat(dictFilePath); dErr != nil {
		return "", fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	// expected is the true cohort size, known from stage 1 (CreateGvcfs sets
	// opts.ExpectedSamples) or discovered here when running standalone. It is
	// only inferred from the gVCFs themselves as a last resort, since that
	// can't tell a sample missing from every chromosome apart from a full
	// cohort — see the fallback below.
	expected := opts.ExpectedSamples

	if gvcfs == nil {
		color.Cyan("================================== Finding gVCFs ==================================\n\n")
		var sampleCount int
		if gvcfs, sampleCount, err = FindExistingGvcfs(opts); err != nil {
			return "", err
		}
		if expected == 0 {
			expected = sampleCount
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

	// Fallback for callers with no real sample count available (e.g. a
	// standalone run given explicit --gvcf paths): infer it from the largest
	// gVCF set seen across chromosomes. Anything short of it is missing
	// samples and is not safe to joint-call.
	if expected == 0 {
		for _, paths := range gvcfs {
			if len(paths) > expected {
				expected = len(paths)
			}
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

	// Chromosome groups are independent — separate gVCFs in, separate work
	// directory, separate output — so they run concurrently. How many at once is
	// bounded by memory rather than cores: a GATK job asks the JVM for a 12 GB
	// heap and a GLnexus job will use whatever it is given.
	jobs, memPerJobGB := mergeJobs(opts, merger, len(labels))

	color.Cyan("=========================== Merging %d chromosome group(s) with %s ===========================\n\n",
		len(labels), merger)
	if merger == "gatk" {
		importHeapGB, genotypeHeapGB := gatkMergeHeaps(memPerJobGB)
		color.Cyan("Running %d merge job(s) at a time, %d GB each (import heap %dg, genotype heap %dg)\n\n",
			jobs, memPerJobGB, importHeapGB, genotypeHeapGB)
	} else {
		color.Cyan("Running %d merge job(s) at a time, %d GB and %d thread(s) each\n\n",
			jobs, memPerJobGB, jobThreads(opts))
	}

	// Results are indexed by label, not appended, so the dictionary order
	// established above survives workers finishing out of order — gatk MergeVcfs
	// will not reorder its inputs. Each worker writes to its own index, so no
	// lock is needed. An empty string means that group did not produce a VCF.
	results := make([]string, len(labels))

	taskCh := make(chan int, len(labels))
	for i := range labels {
		taskCh <- i
	}
	close(taskCh)

	var wg sync.WaitGroup
	for w := 1; w <= jobs; w++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()

			// Only tag lines when there is more than one worker to tell apart.
			tag := ""
			if jobs > 1 {
				tag = fmt.Sprintf("[Worker %d] ", workerID)
			}

			for i := range taskCh {
				label := labels[i]

				// An explicit --gvcf list arrives under a synthetic label ("all") that is
				// not a sequence name, so fall back to the whole reference for it.
				seqs := seqsFor[label]
				if len(seqs) == 0 {
					seqs = append(append([]SeqInfo{}, chroms...), contigs...)
				}

				results[i] = mergeOneGroup(opts, merger, label, seqs, gvcfs[label], tag, memPerJobGB)
			}
		}(w)
	}
	wg.Wait()

	var jointVcfs []string
	for _, vcf := range results {
		if vcf != "" {
			jointVcfs = append(jointVcfs, vcf)
		}
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

// mergeOneGroup produces the joint VCF for one chromosome group and returns its
// path, or "" if the group could not be merged. A failure is reported and
// returned as "" rather than aborting, so one bad chromosome does not cost the
// whole cohort.
//
// An existing joint VCF is reused when it is valid and already holds exactly the
// samples these gVCFs carry.
func mergeOneGroup(opts Options, merger, label string, seqs []SeqInfo, gvcfs []string, tag string, memPerJobGB int) string {
	jointVCF := JointVcfPath(opts, label)

	// ------------------------- reuse a valid, up-to-date joint VCF ------------------------- //
	if _, sErr := os.Stat(jointVCF); sErr == nil {
		reuse := true

		// --skip-verification skips the integrity check only.
		if !opts.SkipVerification {
			if vErr := utils.ValidateGvcf(jointVCF, opts.Verbose, opts.Quick); vErr != nil {
				color.Yellow("%s[%s] joint VCF is corrupt, re-merging: %v\n", tag, label, vErr)
				reuse = false
			}
		}

		// The sample-set check always runs, even under --skip-verification: it is
		// only a header read, and skipping it means a newly added sample silently
		// never reaches the output.
		if reuse {
			want, wErr := allGvcfSampleNames(gvcfs)
			have, hErr := vcfSampleNames(jointVCF)
			switch {
			case wErr != nil || hErr != nil:
				color.Yellow("%s[%s] could not compare sample names, re-merging\n", tag, label)
				reuse = false
			case !sampleNamesMatch(want, have):
				color.Yellow("%s[%s] sample set changed (%d gVCFs vs %d samples in the joint VCF), re-merging\n",
					tag, label, len(want), len(have))
				reuse = false
			}
		}

		if reuse {
			color.Green("%s[%s] joint VCF is up to date, reusing: %s\n", tag, label, jointVCF)
			return jointVCF
		}
		os.Remove(jointVCF)
		os.Remove(jointVCF + ".tbi")
	}

	color.Cyan("%s[%s] merging %d gVCFs ...\n\n", tag, label, len(gvcfs))

	var mErr error
	if merger == "gatk" {
		mErr = mergeChromGATK(opts, label, seqs, gvcfs, jointVCF, tag, memPerJobGB)
	} else {
		mErr = mergeChromGlnexus(opts, label, gvcfs, jointVCF, tag, memPerJobGB, jobThreads(opts))
	}
	if mErr != nil {
		color.Red("%s[%s] merge FAILED: %v\n", tag, label, mErr)
		return ""
	}

	color.Green("%s[%s] joint VCF created: %s\n\n", tag, label, jointVCF)
	return jointVCF
}

// Neither merger has a fixed appetite, so these are floors: the least memory a
// job of each kind is worth starting with. Each job is then handed an equal
// share of what is actually free, and the GATK heaps are sized from that share
// (see gatkMergeHeaps) rather than being fixed in the command line.
//
// 8 GB for GATK is what a GenotypeGVCFs heap of 6 GB plus JVM metaspace, thread
// stacks and GenomicsDB's off-heap buffers comes to. Budgeting the old fixed
// 12 GB heap instead cost real parallelism: on a 32 GB machine it allowed one
// merge job, where two chromosomes could have run side by side.
const (
	gatkMergeJobMemGB  = 8
	glnexusJobMinMemGB = 8

	// unknownMemJobCap bounds parallelism when free memory cannot be read at
	// all. Merging is the most memory-hungry stage in the pipeline and an
	// over-subscribed machine will OOM-kill a job hours in, so the blind default
	// stays small; --merge-jobs raises it.
	unknownMemJobCap = 2
)

// gatkMergeHeaps sizes the two JVM heaps in a GATK merge job from that job's
// memory share.
//
// GenotypeGVCFs holds the most, so it gets three quarters of the share, leaving
// the rest for metaspace, thread stacks and GC structures. GenomicsDBImport
// streams in batches of 50 samples and needs less, but its buffers live off the
// heap, so it gets half.
//
// Both are clamped. The floors keep a small share from producing a heap too
// small to start; the ceilings are the values this pipeline used to hardcode,
// which were sized for a large cohort and are ample — handing GenotypeGVCFs a
// 25 GB heap because a machine happens to be idle buys nothing and makes GC
// pauses worse.
func gatkMergeHeaps(memPerJobGB int) (importGB, genotypeGB int) {
	return clampGB(memPerJobGB/2, 2, 8), clampGB(memPerJobGB*3/4, 4, 12)
}

func clampGB(v, lo, hi int) int {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}

// jobThreads is the per-job thread count, guarding against a zero or negative
// --threads reaching a tool that would reject it.
func jobThreads(opts Options) int {
	if opts.Threads < 1 {
		return 1
	}
	return opts.Threads
}

// availableMemGB reports free RAM in GB from /proc/meminfo's MemAvailable, which
// unlike MemFree counts reclaimable page cache. It returns 0 when that cannot be
// read — a non-Linux host, or a restricted container — and callers fall back to
// a fixed cap.
func availableMemGB() int {
	data, err := os.ReadFile("/proc/meminfo")
	if err != nil {
		return 0
	}
	for _, line := range strings.Split(string(data), "\n") {
		if !strings.HasPrefix(line, "MemAvailable:") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return 0
		}
		kb, cErr := strconv.Atoi(fields[1])
		if cErr != nil {
			return 0
		}
		return kb / (1024 * 1024)
	}
	return 0
}

// mergeJobs decides how many chromosome groups to merge at once, and how much
// memory each job may use.
//
// Memory is the binding constraint, not cores: the cores/--threads figure that
// bounds gVCF creation would start a dozen JVMs here. An explicit --merge-jobs
// is honoured as given, on the basis that someone setting it knows their
// machine, but it is warned about when it oversubscribes what is free.
//
// memPerJobGB is each job's equal share of free memory. GLnexus is handed it
// directly; the GATK path sizes its heaps from it.
func mergeJobs(opts Options, merger string, groups int) (jobs, memPerJobGB int) {
	if groups < 1 {
		return 1, glnexusJobMinMemGB
	}

	avail := availableMemGB()
	minPerJob := gatkMergeJobMemGB
	if merger == "glnexus" {
		minPerJob = glnexusJobMinMemGB
	}

	if opts.MergeJobs > 0 {
		jobs = opts.MergeJobs
		if avail > 0 && jobs > groups {
			// Capped to groups below, so only warn about what will actually run.
			jobs = groups
		}
		if avail > 0 && jobs*minPerJob > avail {
			color.Yellow("--merge-jobs %d wants about %d GB but only %d GB is free; jobs may be killed\n",
				jobs, jobs*minPerJob, avail)
		}
	} else {
		jobs = runtime.NumCPU() / jobThreads(opts)

		if avail <= 0 {
			if jobs > unknownMemJobCap {
				jobs = unknownMemJobCap
			}
		} else if avail/minPerJob < jobs {
			jobs = avail / minPerJob
		}
	}

	if jobs > groups {
		jobs = groups
	}
	if jobs < 1 {
		jobs = 1
	}

	// The floor keeps an oversubscribed --merge-jobs, or a host whose free memory
	// could not be read, from producing a share too small to size a heap from.
	memPerJobGB = avail / jobs
	if memPerJobGB < minPerJob {
		memPerJobGB = minPerJob
	}
	return jobs, memPerJobGB
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

	// The list is scratch, so it lives in the scratch directory — keyed by the
	// output's own directory, and named after the output, so parallel runs of
	// different caller/merger combinations cannot overwrite each other's list.
	tmpDir := utils.WorkTmpDir(outVCF)
	vcfListPath := filepath.Join(tmpDir,
		strings.TrimSuffix(filepath.Base(outVCF), ".vcf.gz")+".vcfs.list")
	defer os.Remove(vcfListPath)
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

	cmd := fmt.Sprintf(`gatk MergeVcfs -I %s -O %s --tmp-dir %s`, vcfListPath, outVCF, tmpDir)
	if err := utils.RunCmd(cmd, verbose); err != nil {
		return "", fmt.Errorf("gatk MergeVcfs: %w", err)
	}
	return outVCF, nil
}
