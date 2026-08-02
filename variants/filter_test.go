package variants

import (
	"path/filepath"
	"testing"

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
