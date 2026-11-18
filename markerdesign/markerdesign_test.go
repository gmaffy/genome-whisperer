package markerdesign

import (
	"slices"
	"testing"
)

func TestSimplify(t *testing.T) {
	header := []string{"CHROM", "POS", "REF", "ALT", "A.GT", "B.GT", "C.GT"}
	rows := [][]string{
		{"c", "1", "A", "G", "A/A", "A|G", "./."},        // SNP
		{"c", "2", "AT", "A", "AT/AT", "AT/A", "A/A"},    // deletion
		{"c", "3", "A", "G,AT", "A/AT", "G/G", "A/G"},    // MIXED: SNP alleles keep their base
		{"c", "4", "A", "T,*", "*/*", "A/*", "T/T"},      // spanning deletion
		{"c", "5", "ATT", "A,AT", "A/A", "AT/AT", "./."}, // three indel lengths
	}
	want := [][]string{
		{"A", "R", "N", "A", "False"},
		{"+", "0", "-", "+", "False"},
		{"N", "G", "R", "A", "False"},
		{"-", "N", "T", "A", "False"},
		{"-", "+", "N", "+", "True"},
	}
	h, got := simplify(header, rows, []int{4, 5, 6}, 2, 3)
	if !slices.Equal(h[7:], []string{"REF_CODE", "INDEL_MULTILEN"}) {
		t.Fatalf("header = %v", h)
	}
	for i := range rows {
		if !slices.Equal(got[i][4:], want[i]) {
			t.Errorf("row %d: got %v, want %v", i, got[i][4:], want[i])
		}
	}
}

func TestSplits(t *testing.T) {
	cases := []struct {
		a, b []string
		want bool
	}{
		{[]string{"A", "A", "N"}, []string{"G", "G"}, true},
		{[]string{"A", "A"}, []string{"R"}, true},  // het in b is allowed
		{[]string{"R", "R"}, []string{"G"}, false}, // a must be homozygous
		{[]string{"A", "G"}, []string{"T"}, false}, // a must agree
		{[]string{"A"}, []string{"G", "N"}, false}, // b must be fully called
		{[]string{"A"}, []string{"G", "A"}, false}, // b never carries a's call
		{[]string{"N", "N"}, []string{"G"}, false}, // a has no call
	}
	for _, c := range cases {
		if got := splits(c.a, c.b); got != c.want {
			t.Errorf("splits(%v, %v) = %v, want %v", c.a, c.b, got, c.want)
		}
	}
}

func TestAddPrimers(t *testing.T) {
	header := []string{"CHROM", "POS", "REF", "ALT"}
	col := map[string]int{"CHROM": 0, "POS": 1, "REF": 2, "ALT": 3}
	seq := "ACGTAGGC"
	vars := []variant{{2, "CG", "TA"}, {4, "TA", "T"}, {7, "G", "A"}} // MNP, the target, SNP
	// TA/T at 4 -> T[A/-]: padding T moves into the left flank
	h, got := addPrimers(header, [][]string{{"c", "4", "TA", "T"}}, col, seq, vars, 3)
	want := []string{"CGT", "GGC", "YRT", "1", "NO", "GRC", "1", "NO", "CGT[A/-]GGC", "YRT[A/-]GRC"}
	if !slices.Equal(h[4:], primerCols) || !slices.Equal(got[0][4:], want) {
		t.Errorf("got %v, want %v", got[0][4:], want)
	}
}
