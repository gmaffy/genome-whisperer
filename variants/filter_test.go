package variants

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/genome-whisperer/utils"
)

// gatkDefaults mirrors the flag defaults on VariantCalling.
func gatkDefaults() utils.HardFilterConfig {
	return utils.HardFilterConfig{
		SNP_QD_Min: 2.0, SNP_QUAL_Min: 30.0, SNP_SOR_Max: 3.0, SNP_FS_Max: 60.0,
		SNP_MQ_Min: 40.0, SNP_MQRankSum_Min: -12.5, SNP_ReadPosRankSum_Min: -8.0,
		INDEL_QD_Min: 2.0, INDEL_QUAL_Min: 30.0, INDEL_FS_Max: 200.0,
		INDEL_ReadPosRankSum_Min: -20.0, INDEL_SOR_Max: 10.0,
	}
}

// readVariants parses records through vcfgo so the tests exercise the same typing
// path as production rather than hand-built Info maps. parseSamples mirrors
// FilterVcf's lazySamples choice.
func readVariants(t *testing.T, parseSamples bool, samples []string, records ...string) []*vcfgo.Variant {
	t.Helper()
	path := writeVCF(t, filepath.Join(t.TempDir(), "in.vcf"), samples, records...)

	in, cleanup, err := openVCF(path)
	if err != nil {
		t.Fatal(err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(in, !parseSamples)
	if err != nil {
		t.Fatal(err)
	}
	var out []*vcfgo.Variant
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		out = append(out, v)
	}
	if err := rdr.Error(); err != nil {
		t.Fatalf("parsing test VCF: %v", err)
	}
	if len(out) != len(records) {
		t.Fatalf("parsed %d variants, expected %d", len(out), len(records))
	}
	return out
}

const goodInfo = "QD=20;FS=1.0;SOR=0.5;MQ=60;MQRankSum=0.0;ReadPosRankSum=0.0"

// ---------------------------------------------------------------------------
// classifyVariant
// ---------------------------------------------------------------------------

func TestClassifyVariant(t *testing.T) {
	vs := readVariants(t, false, []string{"S1"},
		"A01\t100\t.\tA\tG\t500\t.\t"+goodInfo+"\tGT:GQ\t0/1:60",   // SNP
		"A01\t200\t.\tAT\tA\t500\t.\t"+goodInfo+"\tGT:GQ\t0/1:60",  // deletion
		"A01\t300\t.\tA\tATG\t500\t.\t"+goodInfo+"\tGT:GQ\t0/1:60", // insertion
		"A01\t400\t.\tAT\tGC\t500\t.\t"+goodInfo+"\tGT:GQ\t0/1:60", // MNP
	)

	for i, want := range []struct{ snp, indel, mnp bool }{
		{true, false, false},
		{false, true, false},
		{false, true, false},
		{false, false, true},
	} {
		snp, indel, mnp := classifyVariant(vs[i])
		if snp != want.snp || indel != want.indel || mnp != want.mnp {
			t.Errorf("variant %d: got snp=%v indel=%v mnp=%v, want %+v", i, snp, indel, mnp, want)
		}
	}
}

// ---------------------------------------------------------------------------
// GATK profile
// ---------------------------------------------------------------------------

func TestPassesHardFilterGATKProfile(t *testing.T) {
	cfg := gatkDefaults()

	cases := []struct {
		name string
		info string
		qual string
		want bool
	}{
		{"all annotations good", goodInfo, "500", true},
		{"QUAL below minimum", goodInfo, "10", false},
		{"QD below minimum", "QD=1.0;FS=1.0;SOR=0.5;MQ=60", "500", false},
		{"FS above maximum", "QD=20;FS=99;SOR=0.5;MQ=60", "500", false},
		{"SOR above maximum", "QD=20;FS=1.0;SOR=9.9;MQ=60", "500", false},
		{"MQ below minimum", "QD=20;FS=1.0;SOR=0.5;MQ=10", "500", false},
		{"MQRankSum below minimum", "QD=20;FS=1.0;SOR=0.5;MQ=60;MQRankSum=-20.0", "500", false},
		{"ReadPosRankSum below minimum", "QD=20;FS=1.0;SOR=0.5;MQ=60;ReadPosRankSum=-30.0", "500", false},
		// Each threshold is guarded by a presence check, so a missing annotation is
		// not a failure. That is deliberate for GATK output but is exactly why the
		// GATK profile is useless on DeepVariant output, which has none of them.
		{"no annotations at all", ".", "500", true},
	}

	for _, c := range cases {
		v := readVariants(t, false, []string{"S1"},
			"A01\t100\t.\tA\tG\t"+c.qual+"\t.\t"+c.info+"\tGT:GQ\t0/1:60")[0]
		if got := PassesHardFilter(v, cfg); got != c.want {
			t.Errorf("%s: PassesHardFilter = %v, want %v", c.name, got, c.want)
		}
	}
}

// A zero-valued HardFilterConfig sets every maximum to 0, so any variant with
// FS or SOR above zero is discarded. bsaseq passed exactly that with filtering
// enabled, which gutted its VCF.
func TestZeroHardFilterConfigDiscardsAlmostEverything(t *testing.T) {
	v := readVariants(t, false, []string{"S1"},
		"A01\t100\t.\tA\tG\t500\t.\t"+goodInfo+"\tGT:GQ\t0/1:60")[0]

	if PassesHardFilter(v, utils.HardFilterConfig{}) {
		t.Error("a zero config should reject a variant with SOR>0 and FS>0; " +
			"if this passes, the bsaseq default no longer matters")
	}
	if !PassesHardFilter(v, gatkDefaults()) {
		t.Error("the same variant must pass with real defaults")
	}
}

// ---------------------------------------------------------------------------
// DeepVariant profile
// ---------------------------------------------------------------------------

func TestPassesDeepVariantQualAndGQ(t *testing.T) {
	cfg := gatkDefaults() // supplies SNP_QUAL_Min / INDEL_QUAL_Min

	cases := []struct {
		name    string
		record  string
		minGQ   int
		want    bool
		comment string
	}{
		{
			name:   "QUAL ok and one sample confident",
			record: "A01\t100\t.\tA\tG\t500\t.\t.\tGT:GQ\t0/1:60\t0/0:5",
			minGQ:  20, want: true,
		},
		{
			name:   "QUAL below minimum",
			record: "A01\t100\t.\tA\tG\t5\t.\t.\tGT:GQ\t0/1:60\t0/1:60",
			minGQ:  20, want: false,
		},
		{
			name:   "no sample reaches minGQ",
			record: "A01\t100\t.\tA\tG\t500\t.\t.\tGT:GQ\t0/1:5\t0/0:9",
			minGQ:  20, want: false,
		},
		{
			name:   "exactly at minGQ is kept",
			record: "A01\t100\t.\tA\tG\t500\t.\t.\tGT:GQ\t0/1:20\t0/0:1",
			minGQ:  20, want: true,
		},
		{
			name:   "indel uses the indel QUAL threshold",
			record: "A01\t100\t.\tAT\tA\t20\t.\t.\tGT:GQ\t0/1:60\t0/0:60",
			minGQ:  20, want: false,
		},
	}

	for _, c := range cases {
		v := readVariants(t, true, []string{"S1", "S2"}, c.record)[0]
		if got := passesDeepVariant(v, cfg, c.minGQ); got != c.want {
			t.Errorf("%s: passesDeepVariant = %v, want %v", c.name, got, c.want)
		}
	}
}

// The GATK annotations DeepVariant never emits must not cause a rejection: the
// DeepVariant profile judges QUAL and GQ only.
func TestPassesDeepVariantIgnoresMissingGATKAnnotations(t *testing.T) {
	v := readVariants(t, true, []string{"S1"},
		"A01\t100\t.\tA\tG\t500\t.\t.\tGT:GQ\t0/1:60")[0]

	if !passesDeepVariant(v, gatkDefaults(), 20) {
		t.Error("a DeepVariant record with no QD/FS/SOR/MQ should still pass on QUAL+GQ")
	}
}

// A sites-only VCF has no genotypes to judge, so QUAL is all there is.
func TestPassesDeepVariantSitesOnly(t *testing.T) {
	v := readVariants(t, true, nil, "A01\t100\t.\tA\tG\t500\t.\t.")[0]
	if !passesDeepVariant(v, gatkDefaults(), 20) {
		t.Error("sites-only record above the QUAL threshold should be kept")
	}

	low := readVariants(t, true, nil, "A01\t100\t.\tA\tG\t5\t.\t.")[0]
	if passesDeepVariant(low, gatkDefaults(), 20) {
		t.Error("sites-only record below the QUAL threshold should be dropped")
	}
}

// This is the bug that writing these tests exposed. FilterVcf called
// vcfgo.NewReader(in, true), i.e. lazySamples, so Variant.Samples was empty and
// passesDeepVariant took its len(v.Samples)==0 branch and returned true for every
// record: the GQ filter never ran. The reader must parse samples for this profile.
func TestDeepVariantProfileNeedsParsedSamples(t *testing.T) {
	record := "A01\t100\t.\tA\tG\t500\t.\t.\tGT:GQ\t0/1:5\t0/0:9" // both below minGQ

	withSamples := readVariants(t, true, []string{"S1", "S2"}, record)[0]
	if len(withSamples.Samples) == 0 {
		t.Fatal("expected parsed samples when lazySamples is false")
	}
	if passesDeepVariant(withSamples, gatkDefaults(), 20) {
		t.Error("with samples parsed, a record where no sample reaches minGQ must be dropped")
	}

	lazy := readVariants(t, false, []string{"S1", "S2"}, record)[0]
	if len(lazy.Samples) != 0 {
		t.Skip("this vcfgo build populates Samples even when lazy; the guard below is moot")
	}
	if !passesDeepVariant(lazy, gatkDefaults(), 20) {
		t.Fatal("unexpected: lazy variant was filtered")
	}
	t.Log("confirmed: with lazySamples the GQ check silently passes everything, " +
		"which is why FilterVcf now parses samples for the DeepVariant profile")
}

// ---------------------------------------------------------------------------
// Reuse of an existing filtered VCF
// ---------------------------------------------------------------------------

// writeStampedVCF writes a VCF carrying a GenomeWhispererHardFilter header line,
// the shape FilterVcf leaves behind.
func writeStampedVCF(t *testing.T, path, stamp string, samples []string, records ...string) string {
	t.Helper()
	writeVCF(t, path, samples, records...)

	body, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	stamped := strings.Replace(string(body), "#CHROM",
		"##"+filterStampKey+"="+stamp+"\n#CHROM", 1)
	if err := os.WriteFile(path, []byte(stamped), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestFilterStampSeparatesSettings(t *testing.T) {
	base := filterStamp(gatkDefaults(), "gatk", 20)

	changed := gatkDefaults()
	changed.SNP_QD_Min = 4.0

	light := gatkDefaults()
	light.LightFilter = true

	cases := []struct {
		name  string
		stamp string
	}{
		{"one threshold changed", filterStamp(changed, "gatk", 20)},
		{"light filtering", filterStamp(light, "gatk", 20)},
		{"other caller profile", filterStamp(gatkDefaults(), "deepvariant", 20)},
		{"other minGQ", filterStamp(gatkDefaults(), "gatk", 30)},
	}
	for _, c := range cases {
		if c.stamp == base {
			t.Errorf("%s: stamp is unchanged (%s), so a rerun would reuse the wrong VCF", c.name, c.stamp)
		}
	}

	if filterStamp(gatkDefaults(), "gatk", 20) != base {
		t.Error("the same settings must stamp identically, or nothing is ever reused")
	}
}

func TestStaleFilteredVcf(t *testing.T) {
	stamp := filterStamp(gatkDefaults(), "gatk", 20)
	record := "A01\t100\t.\tA\tG\t500\t.\t" + goodInfo + "\tGT\t0/1\t0/0"

	// --skip-verification keeps the check to header reads: the integrity step
	// shells out to bcftools, which these tests do not depend on.
	opts := Options{SkipVerification: true}

	cases := []struct {
		name     string
		joint    []string // samples in the joint VCF
		filtered []string // samples in the filtered VCF
		stamp    string
		reuse    bool
	}{
		{"same samples and settings", []string{"S1", "S2"}, []string{"S1", "S2"}, stamp, true},
		{"samples in a different order", []string{"S1", "S2"}, []string{"S2", "S1"}, stamp, true},
		{"a sample was added", []string{"S1", "S2", "S3"}, []string{"S1", "S2"}, stamp, false},
		{"a sample was dropped", []string{"S1"}, []string{"S1", "S2"}, stamp, false},
		{"thresholds changed", []string{"S1", "S2"}, []string{"S1", "S2"},
			filterStamp(gatkDefaults(), "deepvariant", 20), false},
		{"no stamp at all", []string{"S1", "S2"}, []string{"S1", "S2"}, "", false},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dir := t.TempDir()

			jointRecord := "A01\t100\t.\tA\tG\t500\t.\t" + goodInfo + "\tGT" +
				strings.Repeat("\t0/1", len(c.joint))
			joint := writeVCF(t, filepath.Join(dir, "cohort.vcf"), c.joint, jointRecord)

			filtered := filepath.Join(dir, "cohort.hard_filtered.vcf")
			if c.stamp == "" {
				writeVCF(t, filtered, c.filtered, record)
			} else {
				writeStampedVCF(t, filtered, c.stamp, c.filtered, record)
			}

			why := staleFilteredVcf(opts, joint, filtered, stamp)
			if c.reuse && why != "" {
				t.Errorf("expected reuse, got refused: %s", why)
			}
			if !c.reuse && why == "" {
				t.Error("expected the filtered VCF to be rebuilt, but it was accepted for reuse")
			}
		})
	}
}

// The stamp is only useful if it survives the writer. vcfgo keeps unrecognised
// ##key=value lines in Header.Extras, and this pins that: append a stamp on the
// way in, read it back off the file that comes out.
func TestFilterStampRoundTripsThroughVcfgo(t *testing.T) {
	dir := t.TempDir()
	in := writeVCF(t, filepath.Join(dir, "in.vcf"), []string{"S1"},
		"A01\t100\t.\tA\tG\t500\t.\t"+goodInfo+"\tGT\t0/1")

	r, cleanup, err := openVCF(in)
	if err != nil {
		t.Fatal(err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(r, true)
	if err != nil {
		t.Fatal(err)
	}

	stamp := filterStamp(gatkDefaults(), "gatk", 20)
	rdr.Header.Extras = append(rdr.Header.Extras, "##"+filterStampKey+"="+stamp)

	out := filepath.Join(dir, "out.vcf")
	f, err := os.Create(out)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := vcfgo.NewWriter(f, rdr.Header); err != nil {
		t.Fatal(err)
	}
	f.Close()

	got, err := vcfHeaderValue(out, filterStampKey)
	if err != nil {
		t.Fatal(err)
	}
	if got != stamp {
		t.Errorf("stamp did not survive the writer:\n got %q\nwant %q", got, stamp)
	}
}

func TestVcfHeaderValueMissingKey(t *testing.T) {
	path := writeVCF(t, filepath.Join(t.TempDir(), "plain.vcf"), []string{"S1"},
		"A01\t100\t.\tA\tG\t500\t.\t"+goodInfo+"\tGT\t0/1")

	got, err := vcfHeaderValue(path, filterStampKey)
	if err != nil {
		t.Fatal(err)
	}
	if got != "" {
		t.Errorf("expected no stamp, got %q", got)
	}
}

// End to end: the reuse check and the stamp FilterVcf writes have to agree, and
// only a run against the real bgzf/tabix output proves they do. Skipped where
// the tools are absent, since nothing else in this package needs them.
func TestFilterVcfReusesItsOwnOutput(t *testing.T) {
	for _, tool := range []string{"bgzip", "tabix", "bcftools"} {
		if _, err := exec.LookPath(tool); err != nil {
			t.Skipf("%s is not installed", tool)
		}
	}

	dir := t.TempDir()
	plain := writeVCF(t, filepath.Join(dir, "cohort.vcf"), []string{"S1", "S2"},
		"A01\t100\t.\tA\tG\t500\t.\t"+goodInfo+"\tGT\t0/1\t0/0",
		"A01\t200\t.\tA\tG\t5\t.\t"+goodInfo+"\tGT\t0/1\t0/0")
	if out, err := exec.Command("bgzip", "-f", plain).CombinedOutput(); err != nil {
		t.Fatalf("bgzip: %v\n%s", err, out)
	}
	joint := plain + ".gz"
	if out, err := exec.Command("tabix", "-f", "-p", "vcf", joint).CombinedOutput(); err != nil {
		t.Fatalf("tabix: %v\n%s", err, out)
	}

	opts := Options{Caller: "gatk", HardFilter: gatkDefaults(), MinGQ: 20}

	first, err := FilterVcf(opts, joint)
	if err != nil {
		t.Fatalf("first run: %v", err)
	}
	before := modTime(t, first)

	stamp, err := vcfHeaderValue(first, filterStampKey)
	if err != nil {
		t.Fatal(err)
	}
	if stamp == "" {
		t.Fatal("FilterVcf wrote no stamp, so no later run can reuse its output")
	}

	if _, err := FilterVcf(opts, joint); err != nil {
		t.Fatalf("second run: %v", err)
	}
	if got := modTime(t, first); !got.Equal(before) {
		t.Error("unchanged settings rewrote the filtered VCF instead of reusing it")
	}

	tighter := opts
	tighter.HardFilter.SNP_QUAL_Min = 1000
	if _, err := FilterVcf(tighter, joint); err != nil {
		t.Fatalf("third run: %v", err)
	}
	if got := modTime(t, first); got.Equal(before) {
		t.Error("a changed threshold reused the previous filtered VCF")
	}
}

func modTime(t *testing.T, path string) time.Time {
	t.Helper()
	fi, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	return fi.ModTime()
}
