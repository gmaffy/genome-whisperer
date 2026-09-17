package annotation

import (
	"os"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// A hand-built table with the two things the reducer has to get right: more
// than one query, and a tie at the top bitscore. Only the first and last fields
// are ever inspected, so the middle is filled with plausible values.
func row(qseqID string, bitscore string, marker string) string {
	return strings.Join([]string{
		qseqID, "558", marker, "481", "61.538", "143", "18", "1",
		"118", "223", "291", "433", "3.22e-40", bitscore,
	}, "\t")
}

func TestTopHitsPicksHighestBitscorePerQuery(t *testing.T) {
	in := strings.Join([]string{
		row("geneB", "100", "b-low"),
		row("geneA", "156", "a-mid"),
		row("geneA", "400", "a-best"),
		row("geneB", "250", "b-best"),
		row("geneA", "12", "a-low"),
	}, "\n") + "\n"

	rows, malformed, err := topHits(strings.NewReader(in))
	if err != nil {
		t.Fatalf("topHits: %v", err)
	}
	if malformed != 0 {
		t.Errorf("malformed = %d, want 0", malformed)
	}
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want one per query (2)", len(rows))
	}
	// Sorted by query id, not by the order they appeared.
	if got := field(rows[0], 0); got != "geneA" {
		t.Errorf("first row query = %q, want geneA (rows should be sorted)", got)
	}
	if got := field(rows[0], 2); got != "a-best" {
		t.Errorf("geneA kept %q, want a-best", got)
	}
	if got := field(rows[1], 2); got != "b-best" {
		t.Errorf("geneB kept %q, want b-best", got)
	}
}

// The tie rule is load-bearing: 8.4% of queries in the real Maxima table have
// more than one hit at the maximum bitscore, and pandas' .iloc[0] kept the
// first. A >= comparison here would silently report a different subject.
func TestTopHitsTieKeepsFirstRow(t *testing.T) {
	in := strings.Join([]string{
		row("geneA", "300", "first"),
		row("geneA", "300", "second"),
		row("geneA", "300", "third"),
	}, "\n") + "\n"

	rows, _, err := topHits(strings.NewReader(in))
	if err != nil {
		t.Fatalf("topHits: %v", err)
	}
	if len(rows) != 1 {
		t.Fatalf("got %d rows, want 1", len(rows))
	}
	if got := field(rows[0], 2); got != "first" {
		t.Errorf("tie kept %q, want the earliest row (first)", got)
	}
}

// Rows of the wrong shape are counted and skipped, never returned early on —
// the caller is draining a pipe that blastp is still writing into.
func TestTopHitsCountsMalformedRows(t *testing.T) {
	in := strings.Join([]string{
		row("geneA", "300", "good"),
		"not\ta\ttable",
		row("geneB", "notanumber", "bad-score"),
		"",
		row("geneB", "50", "good-b"),
	}, "\n") + "\n"

	rows, malformed, err := topHits(strings.NewReader(in))
	if err != nil {
		t.Fatalf("topHits: %v", err)
	}
	if malformed != 2 {
		t.Errorf("malformed = %d, want 2", malformed)
	}
	if len(rows) != 2 {
		t.Errorf("got %d rows, want 2", len(rows))
	}
}

func TestTopHitsEmptyInput(t *testing.T) {
	rows, malformed, err := topHits(strings.NewReader(""))
	if err != nil {
		t.Fatalf("topHits: %v", err)
	}
	if len(rows) != 0 || malformed != 0 {
		t.Errorf("got %d rows / %d malformed, want 0/0", len(rows), malformed)
	}
}

// The real thing: 1500 rows taken verbatim from Maximav1.1.PRGblast.txt, the
// table the Python pipeline produced for the squash genome.
func TestTopHitsAgainstRealBlastTable(t *testing.T) {
	f, err := os.Open("testdata/prgblast_sample.tsv")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	rows, malformed, err := topHits(f)
	if err != nil {
		t.Fatalf("topHits: %v", err)
	}
	if malformed != 0 {
		t.Errorf("malformed = %d, want 0 on a real blastp table", malformed)
	}

	// Recompute the expected answer independently, straight from the file.
	raw, err := os.ReadFile("testdata/prgblast_sample.tsv")
	if err != nil {
		t.Fatal(err)
	}
	wantBest := map[string]float64{}
	var order []string
	for _, line := range strings.Split(strings.TrimRight(string(raw), "\n"), "\n") {
		parts := strings.Split(line, "\t")
		score, convErr := strconv.ParseFloat(parts[13], 64)
		if convErr != nil {
			t.Fatalf("fixture row has an unparseable bitscore: %q", line)
		}
		if prev, seen := wantBest[parts[0]]; !seen {
			wantBest[parts[0]] = score
			order = append(order, parts[0])
		} else if score > prev {
			wantBest[parts[0]] = score
		}
	}

	if len(rows) != len(wantBest) {
		t.Fatalf("got %d rows, want one per distinct query (%d)", len(rows), len(wantBest))
	}

	sort.Strings(order)
	for i, r := range rows {
		q := field(r, 0)
		if q != order[i] {
			t.Fatalf("row %d query = %q, want %q (output must be sorted by QseqID)", i, q, order[i])
		}
		got, convErr := strconv.ParseFloat(field(r, 13), 64)
		if convErr != nil {
			t.Fatalf("row %d bitscore %q: %v", i, field(r, 13), convErr)
		}
		if got != wantBest[q] {
			t.Errorf("%s kept bitscore %v, want the maximum %v", q, got, wantBest[q])
		}
	}
}

func TestPrgTopHitsName(t *testing.T) {
	for _, tc := range []struct{ species, version, want string }{
		{"Maxima", "v1.1", "MAXIMAV1.1_PRG_TOP_HITS.txt"},
		{"Tomato", "v4.1", "TOMATOV4.1_PRG_TOP_HITS.txt"},
		{"Moschata", "", "MOSCHATA_PRG_TOP_HITS.txt"},
	} {
		if got := PrgTopHitsName(tc.species, tc.version); got != tc.want {
			t.Errorf("PrgTopHitsName(%q, %q) = %q, want %q", tc.species, tc.version, got, tc.want)
		}
	}
}

// The header must stay in step with the -outfmt column list, since every reader
// downstream looks columns up by name and then indexes by position.
func TestPrgHeaderMatchesOutfmt(t *testing.T) {
	cols := strings.Fields(prgColumns)
	if cols[0] != "6" {
		t.Fatalf("prgColumns should start with the tabular format 6, got %q", cols[0])
	}
	if got, want := len(cols)-1, prgFields; got != want {
		t.Errorf("prgColumns names %d fields, prgFields says %d", got, want)
	}
	if got := len(strings.Split(prgHeader, "\t")); got != prgFields {
		t.Errorf("prgHeader has %d columns, want %d", got, prgFields)
	}
	// The names AddPrg and genespace.genePrgMap require.
	names := strings.Split(prgHeader, "\t")
	for _, required := range []string{"QseqID", "Qlen", "PercIdent", "Length", "Eval"} {
		found := false
		for _, n := range names {
			if n == required {
				found = true
			}
		}
		if !found {
			t.Errorf("prgHeader is missing the required column %q", required)
		}
	}
	if names[bitscoreIdx] != "Bitscore" {
		t.Errorf("bitscoreIdx points at %q, want Bitscore", names[bitscoreIdx])
	}
}

func field(row string, i int) string {
	parts := strings.Split(row, "\t")
	if i >= len(parts) {
		return ""
	}
	return parts[i]
}

// The banner is meant to be pasted into a shell, which means -outfmt's value —
// the one argument containing spaces — has to come back quoted.
func TestShellQuoted(t *testing.T) {
	got := shellQuoted([]string{"blastp", "-outfmt", "6 qseqid bitscore", "-query", "/a/b.fa"})
	want := `blastp -outfmt '6 qseqid bitscore' -query /a/b.fa`
	if got != want {
		t.Errorf("shellQuoted() = %q, want %q", got, want)
	}
	if got := shellQuoted([]string{"x", "it's"}); got != `x 'it'\''s'` {
		t.Errorf("shellQuoted(embedded quote) = %q", got)
	}
}
