package warehouse

import (
	"os"
	"path/filepath"
	"testing"
)

// mountedOrSkip skips a test when the estate is not attached.
func mountedOrSkip(t *testing.T, paths ...string) {
	t.Helper()
	for _, p := range paths {
		if _, err := os.Stat(p); err != nil {
			t.Skipf("%s not mounted", p)
		}
	}
}

// TestResolveVolumesCollapsesCaseSpellings is the reason identity comes from
// the filesystem rather than the string. These mounts are case-insensitive, so
// three spellings name one directory; treating them as three would treble
// every byte total in the report.
func TestResolveVolumesCollapsesCaseSpellings(t *testing.T) {
	mountedOrSkip(t, "/mnt/u")

	vols, problems := resolveVolumes([]string{"/mnt/u", "/mnt/u/data", "/mnt/u/DATA", "/mnt/u/Data"})
	for _, p := range problems {
		t.Errorf("unexpected problem: %s %s %s", p.Rule, p.Path, p.Detail)
	}
	if len(vols) != 1 {
		t.Fatalf("got %d volumes, want 1", len(vols))
	}
	if len(vols[0].Roots) != 4 {
		t.Errorf("roots = %v, want all four spellings recorded", vols[0].Roots)
	}
	if vols[0].Spelling == "" {
		t.Error("spelling is empty; the verbatim on-disk name is the audit answer")
	}
	t.Logf("one volume: name=%q data_dir=%q spelling=%q", vols[0].Name, vols[0].DataDir, vols[0].Spelling)
}

// TestResolveVolumesKeepsDistinctDevices guards the trap that deduplicating on
// inode alone would spring: on this estate /mnt/u and /mnt/v report the same
// inode number on different devices, so an inode-keyed map loses a volume and
// every crop on it.
func TestResolveVolumesKeepsDistinctDevices(t *testing.T) {
	mountedOrSkip(t, "/mnt/u", "/mnt/v", "/mnt/y", "/mnt/z")

	vols, _ := resolveVolumes([]string{"/mnt/u", "/mnt/v", "/mnt/y", "/mnt/z"})
	if len(vols) != 4 {
		t.Fatalf("got %d volumes, want 4 — distinct devices must not merge", len(vols))
	}

	seenIno := map[uint64][]string{}
	for _, v := range vols {
		seenIno[v.ino] = append(seenIno[v.ino], v.Name)
	}
	for ino, names := range seenIno {
		if len(names) > 1 {
			t.Logf("inode %d is shared by %v — confirms the pair must be the key", ino, names)
		}
	}
}

func TestCapacityIsKnownAndUsesTheMountsBlockSize(t *testing.T) {
	mountedOrSkip(t, "/mnt/u", "/mnt/y")

	vols, _ := resolveVolumes([]string{"/mnt/u", "/mnt/y"})
	if len(vols) != 2 {
		t.Fatalf("got %d volumes, want 2", len(vols))
	}

	for _, v := range vols {
		row, problem := capacityFor(v)
		if problem != nil {
			t.Fatalf("%s: capacity unavailable: %s", v.Name, problem.Detail)
		}
		if !row.CapacityKnown {
			t.Errorf("%s: capacity not known", v.Name)
		}
		if row.BytesTotal <= 0 || row.BlockSize <= 0 {
			t.Errorf("%s: total=%d block=%d, both must be positive", v.Name, row.BytesTotal, row.BlockSize)
		}
		if row.BytesFree > row.BytesTotal {
			t.Errorf("%s: free %d exceeds total %d", v.Name, row.BytesFree, row.BytesTotal)
		}
		t.Logf("%s: %.2f TiB total, %.2f TiB free, %d-byte blocks",
			v.Name,
			float64(row.BytesTotal)/(1<<40),
			float64(row.BytesFree)/(1<<40),
			row.BlockSize)
	}
}

// TestResolveVolumesFindsDataDirPortably covers the resolution logic without
// needing the estate: a synthetic mount holding a data directory with one crop.
func TestResolveVolumesFindsDataDirPortably(t *testing.T) {
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, "DATA", "pepper", "2025", "S1", DirCleanReads), 0o755); err != nil {
		t.Fatal(err)
	}

	vols, problems := resolveVolumes([]string{root})
	for _, p := range problems {
		t.Errorf("unexpected problem: %s %s", p.Rule, p.Detail)
	}
	if len(vols) != 1 {
		t.Fatalf("got %d volumes, want 1", len(vols))
	}
	if vols[0].Spelling != "DATA" {
		t.Errorf("spelling = %q, want the verbatim %q", vols[0].Spelling, "DATA")
	}
}

// TestResolveVolumesAcceptsADataDirDirectly lets a caller point at the data
// directory rather than the mount.
func TestResolveVolumesAcceptsADataDirDirectly(t *testing.T) {
	root := t.TempDir()
	data := filepath.Join(root, "data")
	if err := os.MkdirAll(filepath.Join(data, "pepo", "2026", "S1", DirReferenceGenomes), 0o755); err != nil {
		t.Fatal(err)
	}

	vols, _ := resolveVolumes([]string{data})
	if len(vols) != 1 {
		t.Fatalf("got %d volumes, want 1", len(vols))
	}
	if vols[0].DataDir != data {
		t.Errorf("data dir = %q, want %q", vols[0].DataDir, data)
	}
}

func TestResolveVolumesReportsAMissingRoot(t *testing.T) {
	_, problems := resolveVolumes([]string{filepath.Join(t.TempDir(), "absent")})
	if len(problems) != 1 || problems[0].Rule != RuleDataDirMissing {
		t.Fatalf("problems = %+v, want one %s", problems, RuleDataDirMissing)
	}
}
