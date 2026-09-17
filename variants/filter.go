package variants

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// FilterVcf hard filters a joint VCF and returns the path of the filtered file
// (<input>.hard_filtered.vcf.gz, indexed). An existing filtered VCF that is
// still up to date is reused rather than rebuilt — see staleFilteredVcf.
//
// The threshold profile follows the caller, because the two produce different
// annotations:
//
//   - gatk:        QualByDepth, FisherStrand, StrandOddsRatio, MQ and the rank-sum
//     tests, per GATK best practices.
//   - deepvariant: QUAL and genotype quality only. DeepVariant emits none of the
//     GATK annotations, so applying the GATK profile to its output silently
//     filtered on QUAL alone and made every other threshold a no-op.
func FilterVcf(opts Options, vcf string) (string, error) {
	var filteredVcf string
	switch {
	case strings.HasSuffix(vcf, ".vcf.gz"):
		filteredVcf = strings.TrimSuffix(vcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
	case strings.HasSuffix(vcf, ".vcf"):
		filteredVcf = strings.TrimSuffix(vcf, ".vcf") + ".hard_filtered.vcf.gz"
	default:
		return "", fmt.Errorf("vcf file %q does not end with .vcf or .vcf.gz", vcf)
	}

	cfg := opts.HardFilter
	if cfg.LightFilter {
		// Light filtering keeps everything that is not obviously bad: QUAL only.
		cfg = utils.HardFilterConfig{
			LightFilter:    true,
			SNP_QUAL_Min:   cfg.SNP_QUAL_Min,
			INDEL_QUAL_Min: cfg.INDEL_QUAL_Min,
		}
	}

	minGQ := opts.MinGQ
	if minGQ <= 0 {
		minGQ = 20
	}

	keep := func(v *vcfgo.Variant) bool { return PassesHardFilter(v, cfg) }
	profile := "GATK best practices"
	profileKey := "gatk"
	needSamples := strings.ToLower(opts.Caller) == "deepvariant"
	if needSamples {
		keep = func(v *vcfgo.Variant) bool { return passesDeepVariant(v, cfg, minGQ) }
		profile = fmt.Sprintf("DeepVariant (QUAL + GQ >= %d)", minGQ)
		profileKey = "deepvariant"
	}

	stamp := filterStamp(cfg, profileKey, minGQ)

	// ------------------------ reuse a valid, up-to-date filtered VCF ------------------------ //

	if _, sErr := os.Stat(filteredVcf); sErr == nil {
		why := staleFilteredVcf(opts, vcf, filteredVcf, stamp)
		if why == "" {
			color.Green("Filtered VCF is up to date, reusing: %s\n\n", filteredVcf)
			return filteredVcf, nil
		}
		color.Yellow("Existing filtered VCF %s, re-filtering: %s\n\n", why, filteredVcf)
		os.Remove(filteredVcf)
		os.Remove(filteredVcf + ".tbi")
	}

	color.Cyan("Hard filtering %s using %s\n\n", vcf, profile)

	in, cleanup, err := openVCF(vcf)
	if err != nil {
		return "", fmt.Errorf("open %q: %w", vcf, err)
	}
	defer cleanup()

	// The second argument is lazySamples: when true, vcfgo leaves Variant.Samples
	// empty until the caller parses them. The GATK profile only reads INFO fields,
	// so lazy is right there, but the DeepVariant profile needs per-sample GQ and
	// would otherwise see zero samples and fall back to filtering on QUAL alone.
	rdr, err := vcfgo.NewReader(in, !needSamples)
	if err != nil {
		return "", fmt.Errorf("VCF header %q: %w", vcf, err)
	}

	// Record what is being applied in the output's own header, so the next run
	// can tell "already filtered" from "already filtered, with other thresholds".
	// Any stamp inherited from the input is dropped first: filtering a filtered
	// VCF would otherwise leave two, and the reader takes the first it finds.
	extras := rdr.Header.Extras[:0]
	for _, line := range rdr.Header.Extras {
		if !strings.HasPrefix(line, "##"+filterStampKey+"=") {
			extras = append(extras, line)
		}
	}
	rdr.Header.Extras = append(extras, "##"+filterStampKey+"="+stamp)

	outFile, err := os.Create(filteredVcf)
	if err != nil {
		return "", fmt.Errorf("create %q: %w", filteredVcf, err)
	}
	defer outFile.Close()

	bgzfW := bgzf.NewWriter(outFile, runtime.GOMAXPROCS(0))
	w, err := vcfgo.NewWriter(bgzfW, rdr.Header)
	if err != nil {
		bgzfW.Close()
		return "", fmt.Errorf("VCF writer: %w", err)
	}

	var read, written int
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			continue
		}
		read++
		if keep(v) {
			w.WriteVariant(v)
			written++
		}
	}
	if rErr := rdr.Error(); rErr != nil {
		bgzfW.Close()
		return "", fmt.Errorf("reading variants from %s: %w", vcf, rErr)
	}
	if cErr := bgzfW.Close(); cErr != nil {
		return "", fmt.Errorf("close bgzf: %w", cErr)
	}

	if out, tErr := exec.Command("tabix", "-f", "-p", "vcf", filteredVcf).CombinedOutput(); tErr != nil {
		return "", fmt.Errorf("tabix %q: %w\n%s", filteredVcf, tErr, out)
	}

	color.Green("Kept %d of %d variants -> %s\n\n", written, read, filteredVcf)
	return filteredVcf, nil
}

// filterStampKey names the header line FilterVcf writes into every VCF it
// produces, recording the thresholds that produced it.
const filterStampKey = "GenomeWhispererHardFilter"

// filterStamp renders a threshold set as one VCF header value. Nothing else in
// a filtered VCF says how it was filtered, so this is what lets a later run
// tell an output it can reuse from one built to different thresholds. Writing
// the values out rather than a hash of them means the header doubles as
// provenance for anyone reading the VCF by hand.
//
// cfg is the effective config — after the --light-filter reduction — so a light
// run and a full run with the same QUAL floors still stamp differently.
func filterStamp(cfg utils.HardFilterConfig, profile string, minGQ int) string {
	num := func(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }

	return strings.Join([]string{
		"profile=" + profile,
		"light=" + strconv.FormatBool(cfg.LightFilter),
		"snpQD=" + num(cfg.SNP_QD_Min),
		"snpQUAL=" + num(cfg.SNP_QUAL_Min),
		"snpSOR=" + num(cfg.SNP_SOR_Max),
		"snpFS=" + num(cfg.SNP_FS_Max),
		"snpMQ=" + num(cfg.SNP_MQ_Min),
		"snpMQRankSum=" + num(cfg.SNP_MQRankSum_Min),
		"snpReadPosRankSum=" + num(cfg.SNP_ReadPosRankSum_Min),
		"indelQD=" + num(cfg.INDEL_QD_Min),
		"indelQUAL=" + num(cfg.INDEL_QUAL_Min),
		"indelFS=" + num(cfg.INDEL_FS_Max),
		"indelReadPosRankSum=" + num(cfg.INDEL_ReadPosRankSum_Min),
		"indelSOR=" + num(cfg.INDEL_SOR_Max),
		"minGQ=" + strconv.Itoa(minGQ),
	}, ",")
}

// staleFilteredVcf reports why the filtered VCF beside a joint VCF cannot be
// reused, or "" when it can. It mirrors the reuse check in mergeOneGroup: an
// intact output that holds the right samples and was built the same way is
// worth hours of re-filtering on a large cohort.
//
// Freshness is judged from content, never from timestamps. MergeGvcfs
// re-concatenates the joint VCF on every run, so its mtime is always newer than
// a filtered VCF from a previous run and an mtime comparison would reuse
// nothing.
//
// The gap this leaves: a joint VCF whose records changed while its sample set
// did not — a single chromosome re-merged, say — is not detected, and the
// previous filtered VCF is reused. Delete it to force a rebuild.
func staleFilteredVcf(opts Options, vcf, filteredVcf, stamp string) string {
	// --skip-verification skips the integrity check only, as in mergeOneGroup.
	if !opts.SkipVerification {
		if vErr := utils.ValidateGvcf(filteredVcf, opts.Verbose, opts.Quick); vErr != nil {
			return fmt.Sprintf("is corrupt (%v)", vErr)
		}
	}

	// The remaining checks run even under --skip-verification: they are header
	// reads, and skipping them means a changed threshold or a re-merged cohort
	// silently keeps the previous run's variants.
	had, hErr := vcfHeaderValue(filteredVcf, filterStampKey)
	switch {
	case hErr != nil:
		return fmt.Sprintf("has an unreadable header (%v)", hErr)
	case had == "":
		// Written before FilterVcf stamped its output, so what produced it is
		// unknowable. Rebuilt once, after which the stamp is there.
		return "records no filter settings"
	case had != stamp:
		return fmt.Sprintf("was filtered with different settings (%s)", had)
	}

	// A sites-only input has no sample columns for vcfSampleNames to return, so
	// it lands in the branch below and is re-filtered every run. Nothing in this
	// pipeline produces one — GATK and GLnexus both emit sample columns — so the
	// lost reuse is theoretical.
	want, wErr := vcfSampleNames(vcf)
	if wErr != nil {
		return fmt.Sprintf("cannot be checked against %s (%v)", filepath.Base(vcf), wErr)
	}
	have, sErr := vcfSampleNames(filteredVcf)
	if sErr != nil {
		return fmt.Sprintf("has no readable sample columns (%v)", sErr)
	}
	if !sampleNamesMatch(want, have) {
		return fmt.Sprintf("holds %d samples, %s now holds %d", len(have), filepath.Base(vcf), len(want))
	}

	return ""
}

// vcfHeaderValue returns the value of the first ##<key>=<value> line in a VCF
// header, or "" when the header has none. It stops at #CHROM, so it never
// touches the records.
func vcfHeaderValue(vcf, key string) (string, error) {
	in, cleanup, err := openVCF(vcf)
	if err != nil {
		return "", fmt.Errorf("open %q: %w", vcf, err)
	}
	defer cleanup()

	prefix := "##" + key + "="

	scanner := bufio.NewScanner(in)
	// Header lines get long with many samples.
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, prefix) {
			return strings.TrimPrefix(line, prefix), nil
		}
		if !strings.HasPrefix(line, "##") {
			break
		}
	}
	if err := scanner.Err(); err != nil {
		return "", fmt.Errorf("scanning %s: %w", vcf, err)
	}
	return "", nil
}

// passesDeepVariant applies the DeepVariant threshold profile: site QUAL plus a
// genotype-quality floor.
//
// A site is kept when at least one sample is genotyped confidently (GQ >= minGQ).
// A site where no sample reaches that carries no usable genotype, whatever its
// QUAL.
func passesDeepVariant(v *vcfgo.Variant, cfg utils.HardFilterConfig, minGQ int) bool {
	isSNP, isIndel, isMNP := classifyVariant(v)

	qualMin := cfg.SNP_QUAL_Min
	if isIndel || isMNP {
		qualMin = cfg.INDEL_QUAL_Min
	}
	if !isSNP && !isIndel && !isMNP {
		// Structural or otherwise unclassified: keep, as the GATK profile does.
		return true
	}
	if float64(v.Quality) < qualMin {
		return false
	}

	// Sites-only VCFs have no genotypes to judge; QUAL is all there is.
	if len(v.Samples) == 0 {
		return true
	}
	for _, s := range v.Samples {
		if s != nil && s.GQ >= minGQ {
			return true
		}
	}
	return false
}

// ---------------------------------------------------------------------------
// Absorbed from the retired RunVariantCaller.go / RunVariantCallerDir.go
// ---------------------------------------------------------------------------

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	cleanup := func() { f.Close() }

	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		cleanup = func() { gz.Close(); f.Close() }
		return gz, cleanup, nil
	}
	return f, cleanup, nil
}

func classifyVariant(v *vcfgo.Variant) (isSNP bool, isIndel bool, isMNP bool) {
	refLen := len(v.Ref())
	alts := v.Alt()

	if len(alts) == 0 {
		return false, false, false
	}

	isSNP = true
	isIndel = false
	isMNP = true

	if refLen != 1 {
		isSNP = false
	}

	for _, alt := range alts {
		altLen := len(alt)
		// If lengths differ, it's an Indel
		if altLen != refLen {
			isIndel = true
			isMNP = false
			isSNP = false
		}
	}

	if isIndel {
		isMNP = false
		isSNP = false
	}

	if refLen == 1 {
		isMNP = false
	}

	return isSNP, isIndel, isMNP
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	isSNP, isIndel, isMNP := classifyVariant(v)

	// Pre-fetch all commonly used INFO fields
	qd, hasQD := utils.GetFloat(v, "QD")
	fs, hasFS := utils.GetFloat(v, "FS")
	sor, hasSOR := utils.GetFloat(v, "SOR")
	mq, hasMQ := utils.GetFloat(v, "MQ")
	mqRankSum, hasMQRankSum := utils.GetFloat(v, "MQRankSum")
	readPosRankSum, hasReadPosRankSum := utils.GetFloat(v, "ReadPosRankSum")

	// vcfgo stores Quality directly on the struct
	qual := float64(v.Quality)

	switch {
	case isSNP:
		if qual < hfcfg.SNP_QUAL_Min {
			return false
		}
		if hasQD && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if hasFS && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if hasSOR && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if hasMQ && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if hasMQRankSum && mqRankSum < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if hasReadPosRankSum && readPosRankSum < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel, isMNP:
		// We group MNPs with Indels here to hold them to the same robust standards.
		if qual < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if hasQD && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if hasFS && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if hasSOR && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if hasReadPosRankSum && readPosRankSum < hfcfg.INDEL_ReadPosRankSum_Min {
			return false
		}
		return true

	default:
		// these may be SVs. so we keep
		return true
	}
}
