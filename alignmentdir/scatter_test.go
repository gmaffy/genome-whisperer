package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"testing"
)

var spanPattern = regexp.MustCompile(`^[0-9]+-[0-9]+$`)

// parseShards reads back the interval files written by writeShardIntervals and
// returns, per contig, the intervals claimed by any shard.
func parseShards(t *testing.T, paths []string, lengths map[string]int64) map[string][][2]int64 {
	t.Helper()
	covered := make(map[string][][2]int64)

	for _, path := range paths {
		data, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("reading shard %s: %v", path, err)
		}
		for _, line := range strings.Split(strings.TrimSpace(string(data)), "\n") {
			if line == "" {
				continue
			}
			name, span := line, ""
			// Only a trailing start-end counts as a span: a contig name may
			// itself contain a colon.
			if idx := strings.LastIndex(line, ":"); idx >= 0 && spanPattern.MatchString(line[idx+1:]) {
				name, span = line[:idx], line[idx+1:]
			}
			length, ok := lengths[name]
			if !ok {
				t.Fatalf("shard %s names unknown contig %q", path, name)
			}
			if span == "" {
				covered[name] = append(covered[name], [2]int64{1, length})
				continue
			}
			bounds := strings.SplitN(span, "-", 2)
			if len(bounds) != 2 {
				t.Fatalf("shard %s has malformed interval %q", path, line)
			}
			start, sErr := strconv.ParseInt(bounds[0], 10, 64)
			end, eErr := strconv.ParseInt(bounds[1], 10, 64)
			if sErr != nil || eErr != nil {
				t.Fatalf("shard %s has malformed interval %q", path, line)
			}
			if start < 1 || end > length || start > end {
				t.Fatalf("shard %s interval %q outside 1-%d", path, line, length)
			}
			covered[name] = append(covered[name], [2]int64{start, end})
		}
	}
	return covered
}

// assertFullCoverage checks that the shards tile every contig exactly once: a
// gap silently drops part of the genome from the call set, and an overlap has
// the same region called twice.
func assertFullCoverage(t *testing.T, covered map[string][][2]int64, lengths map[string]int64) {
	t.Helper()
	if len(covered) != len(lengths) {
		t.Fatalf("shards cover %d contigs, reference has %d", len(covered), len(lengths))
	}
	for name, length := range lengths {
		spans := covered[name]
		if len(spans) == 0 {
			t.Fatalf("contig %s is not covered by any shard", name)
		}
		// Intervals are emitted in ascending order per contig, so a single
		// forward walk is enough.
		var next int64 = 1
		for _, span := range spans {
			if span[0] != next {
				t.Fatalf("contig %s: expected next interval to start at %d, got %d", name, next, span[0])
			}
			next = span[1] + 1
		}
		if next-1 != length {
			t.Fatalf("contig %s: shards cover %d of %d bases", name, next-1, length)
		}
	}
}

func writeFai(t *testing.T, dir string, entries []faiEntry) (string, map[string]int64) {
	t.Helper()
	refFasta := filepath.Join(dir, "ref.fa")
	var sb strings.Builder
	lengths := make(map[string]int64, len(entries))
	var offset int64
	for _, e := range entries {
		fmt.Fprintf(&sb, "%s\t%d\t%d\t60\t61\n", e.name, e.length, offset)
		offset += e.length
		lengths[e.name] = e.length
	}
	if err := os.WriteFile(refFasta+".fai", []byte(sb.String()), 0o644); err != nil {
		t.Fatalf("writing fai: %v", err)
	}
	return refFasta, lengths
}

func TestWriteShardIntervalsTilesTheGenome(t *testing.T) {
	dir := t.TempDir()
	refFasta, lengths := writeFai(t, dir, []faiEntry{
		{name: "big1", length: 1_000_000_000},
		{name: "big2", length: 1_000_000_000},
		{name: "small1", length: 10_000_000},
		{name: "small2", length: 10_000_000},
		{name: "small3", length: 10_000_000},
	})

	shards, err := writeShardIntervals(refFasta, filepath.Join(dir, "shards"), 4)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	if len(shards) < 4 {
		t.Fatalf("expected at least 4 shards for 4 requested, got %d", len(shards))
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)
}

func TestWriteShardIntervalsSingleShard(t *testing.T) {
	dir := t.TempDir()
	refFasta, lengths := writeFai(t, dir, []faiEntry{
		{name: "c1", length: 500_000_000},
		{name: "c2", length: 500_000_000},
	})

	shards, err := writeShardIntervals(refFasta, filepath.Join(dir, "shards"), 1)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)
}

// TestWriteShardIntervalsRealReference runs the splitter over a reference on this
// machine when one is present. It is skipped elsewhere; the synthetic cases above
// are what CI relies on.
func TestWriteShardIntervalsRealReference(t *testing.T) {
	refFasta := os.Getenv("GW_TEST_REF")
	if refFasta == "" {
		t.Skip("GW_TEST_REF not set")
	}
	entries, err := readFai(refFasta)
	if err != nil {
		t.Skipf("cannot read fai: %v", err)
	}

	lengths := make(map[string]int64, len(entries))
	var total int64
	for _, e := range entries {
		lengths[e.name] = e.length
		total += e.length
	}

	shards, err := writeShardIntervals(refFasta, t.TempDir(), 24)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)
	t.Logf("%d contigs, %.2f Gb, %d shards, ~%.0f Mb per shard",
		len(entries), float64(total)/1e9, len(shards), float64(total)/float64(len(shards))/1e6)
}

// A genome too small to be worth scattering must come back as a single shard:
// below minShardBases, JVM start-up costs more than the shard's own work.
func TestWriteShardIntervalsSmallGenomeIsNotScattered(t *testing.T) {
	dir := t.TempDir()
	refFasta, lengths := writeFai(t, dir, []faiEntry{
		{name: "c1", length: 1_000_000},
		{name: "c2", length: 1_000_000},
		{name: "c3", length: 1_000_000},
	})

	shards, err := writeShardIntervals(refFasta, filepath.Join(dir, "shards"), 24)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	if len(shards) != 1 {
		t.Fatalf("a 3 Mb genome should not be scattered, got %d shards", len(shards))
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)
}

// The shard count scales with the genome rather than being fixed, so a
// mid-sized reference gets proportionally fewer shards than a 16 Gb one.
func TestWriteShardIntervalsScalesWithGenomeSize(t *testing.T) {
	dir := t.TempDir()
	refFasta, lengths := writeFai(t, dir, []faiEntry{
		{name: "c1", length: 100_000_000},
	})

	shards, err := writeShardIntervals(refFasta, filepath.Join(dir, "shards"), 24)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	if len(shards) != 5 {
		t.Fatalf("100 Mb at %d bases per shard should give 5 shards, got %d", minShardBases, len(shards))
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)
}

// A contig whose name contains a colon cannot be expressed as name:start-end,
// so it has to be passed whole rather than split.
func TestWriteShardIntervalsKeepsColonNamedContigsWhole(t *testing.T) {
	dir := t.TempDir()
	refFasta, lengths := writeFai(t, dir, []faiEntry{
		{name: "chr1:v2", length: 900_000_000},
		{name: "chr2", length: 100_000_000},
	})

	shards, err := writeShardIntervals(refFasta, filepath.Join(dir, "shards"), 4)
	if err != nil {
		t.Fatalf("writeShardIntervals: %v", err)
	}
	assertFullCoverage(t, parseShards(t, shards, lengths), lengths)

	for _, path := range shards {
		data, rErr := os.ReadFile(path)
		if rErr != nil {
			t.Fatalf("reading %s: %v", path, rErr)
		}
		for _, line := range strings.Split(strings.TrimSpace(string(data)), "\n") {
			if strings.HasPrefix(line, "chr1:v2") && line != "chr1:v2" {
				t.Fatalf("colon-named contig was split into %q", line)
			}
		}
	}
}

// The route a reference takes is decided by its total length, so a small genome
// keeps the single-process path it has always used.
func TestUseScatterPicksRouteBySize(t *testing.T) {
	cases := []struct {
		name    string
		entries []faiEntry
		scatter bool
	}{
		{"tomato-sized", []faiEntry{{name: "c1", length: 800_000_000}}, false},
		{"just under threshold", []faiEntry{{name: "c1", length: scatterGenomeBases - 1}}, false},
		{"exactly at threshold", []faiEntry{{name: "c1", length: scatterGenomeBases}}, true},
		{"pepper-sized", []faiEntry{{name: "c1", length: 3_300_000_000}}, true},
		{"summed across contigs", []faiEntry{
			{name: "c1", length: 600_000_000},
			{name: "c2", length: 600_000_000},
		}, true},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			refFasta, _ := writeFai(t, t.TempDir(), tc.entries)
			if got := useScatter(refFasta); got != tc.scatter {
				t.Fatalf("useScatter = %v, want %v", got, tc.scatter)
			}
		})
	}
}

// With no readable .fai, the FASTA's own size stands in, so a missing index
// cannot silently route a large genome down the single-process path.
func TestGenomeBasesFallsBackToFastaSize(t *testing.T) {
	dir := t.TempDir()
	refFasta := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refFasta, make([]byte, 4096), 0o644); err != nil {
		t.Fatalf("writing fasta: %v", err)
	}

	if got := genomeBases(refFasta); got != 4096 {
		t.Fatalf("genomeBases = %d, want the file size 4096", got)
	}
	if genomeBases(filepath.Join(dir, "absent.fa")) != 0 {
		t.Fatal("a reference that is not there should measure 0")
	}
}
