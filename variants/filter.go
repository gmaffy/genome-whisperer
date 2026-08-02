package variants

import (
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strings"

	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// FilterVcf hard filters a joint VCF and returns the path of the filtered file
// (<input>.hard_filtered.vcf.gz, indexed).
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
	if strings.ToLower(opts.Caller) == "deepvariant" {
		keep = func(v *vcfgo.Variant) bool { return passesDeepVariant(v, cfg, minGQ) }
		profile = fmt.Sprintf("DeepVariant (QUAL + GQ >= %d)", minGQ)
	}
	color.Cyan("Hard filtering %s using %s\n\n", vcf, profile)

	in, cleanup, err := openVCF(vcf)
	if err != nil {
		return "", fmt.Errorf("open %q: %w", vcf, err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(in, true)
	if err != nil {
		return "", fmt.Errorf("VCF header %q: %w", vcf, err)
	}

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
