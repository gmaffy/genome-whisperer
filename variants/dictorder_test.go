package variants

import (
	"os"
	"strings"
	"testing"
)

// The reference this pipeline is run against numbers its sequences, so
// dictionary order and lexicographic order disagree from scaffold_100 onwards.
// That disagreement is the whole bug: a multi-interval GenomicsDB workspace is
// genotyped in lexicographic order, which reads as backwards to everything
// downstream.
func TestContigOrderViolation(t *testing.T) {
	// Chr01..Chr11 then scaffold_12.. — the shape of a real .dict.
	order := map[string]int{
		"Chr01": 0, "Chr11": 10,
		"scaffold_12": 11, "scaffold_21": 20, "scaffold_22": 21,
		"scaffold_100": 66, "scaffold_104": 69, "scaffold_1037": 323,
	}

	tests := []struct {
		name    string
		contigs []string
		want    string
	}{
		{
			name:    "dictionary order passes",
			contigs: []string{"Chr01", "Chr11", "scaffold_12", "scaffold_22", "scaffold_100", "scaffold_1037"},
			want:    "",
		},
		{
			// Exactly what GenotypeGVCFs produced for the contigs group.
			name:    "lexicographic order is caught",
			contigs: []string{"scaffold_100", "scaffold_1037", "scaffold_104"},
			want:    "scaffold_104 (dict position 69) comes after scaffold_1037 (dict position 323)",
		},
		{
			name:    "a contig missing from the dictionary is caught",
			contigs: []string{"Chr01", "scaffold_9999"},
			want:    "scaffold_9999 is not in the reference dictionary",
		},
		{
			// A group that only ever holds one sequence cannot be misordered.
			name:    "single contig passes",
			contigs: []string{"scaffold_1037"},
			want:    "",
		},
		{
			name:    "gaps are fine, only going backwards is not",
			contigs: []string{"Chr01", "scaffold_1037"},
			want:    "",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if got := contigOrderViolation(tc.contigs, order); got != tc.want {
				t.Errorf("contigOrderViolation() = %q, want %q", got, tc.want)
			}
		})
	}
}

// dictOrder is what feeds the check, so the two have to agree on what "position"
// means: the order sequences appear in the dict, not their names sorted.
func TestDictOrderMatchesFileOrder(t *testing.T) {
	dict := t.TempDir() + "/ref.dict"
	body := "@HD\tVN:1.6\n" +
		"@SQ\tSN:Chr01\tLN:100\n" +
		"@SQ\tSN:scaffold_22\tLN:50\n" +
		"@SQ\tSN:scaffold_100\tLN:10\n"
	if err := os.WriteFile(dict, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}

	order, err := dictOrder(dict)
	if err != nil {
		t.Fatal(err)
	}
	if order["scaffold_22"] >= order["scaffold_100"] {
		t.Fatalf("dictOrder did not preserve file order: %v", order)
	}

	// And the two compose: file order passes, lexicographic order does not.
	if got := contigOrderViolation([]string{"Chr01", "scaffold_22", "scaffold_100"}, order); got != "" {
		t.Errorf("dictionary order rejected: %s", got)
	}
	if got := contigOrderViolation([]string{"scaffold_100", "scaffold_22"}, order); !strings.Contains(got, "scaffold_22") {
		t.Errorf("lexicographic order accepted, got %q", got)
	}
}
