package utils

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestWorkTmpDirIsUnderOSTempDir(t *testing.T) {
	got := WorkTmpDirFor("/mnt/x/data/beans/2016/Edmund/reference_genomes/v2.1/bams")
	root := filepath.Join(os.TempDir(), "genome-whisperer")
	if !strings.HasPrefix(got, root+string(os.PathSeparator)) {
		t.Errorf("scratch dir %s is not under %s", got, root)
	}
	t.Logf("scratch dir: %s", got)
}

// Two samples' output directories must never share a scratch directory: their
// spills would be visible to each other, and one sample's cleanup would delete
// the other's live temporary files.
func TestWorkTmpDirForSeparatesDirectories(t *testing.T) {
	a := WorkTmpDirFor("/data/Edmund/reference_genomes/v2.1/bams")
	b := WorkTmpDirFor("/data/Kranskop/reference_genomes/v2.1/bams")
	if a == b {
		t.Errorf("two samples share a scratch dir: %s", a)
	}

	// A path whose separators flatten to the same label as another must still
	// get its own directory — that is what the hash in the name is for.
	c := WorkTmpDirFor("/data/x_y/bams")
	d := WorkTmpDirFor("/data/x/y/bams")
	if c == d {
		t.Errorf("flattened labels collided: %s", c)
	}
}

// Cleanup calls WorkTmpDirFor to find what WorkTmpDir created, so the two must
// agree on the path.
func TestWorkTmpDirMatchesWorkTmpDirFor(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "Edmund.RGMD_bqsr.bam")
	created := WorkTmpDir(out)
	if want := WorkTmpDirFor(dir); created != want {
		t.Fatalf("WorkTmpDir = %s, WorkTmpDirFor = %s", created, want)
	}
	if info, err := os.Stat(created); err != nil || !info.IsDir() {
		t.Fatalf("WorkTmpDir did not create %s: %v", created, err)
	}
	_ = os.RemoveAll(created)
}

// The whole point of SpillTmpDir: a whole-genome sort must not be sent into RAM.
func TestSpillTmpDirAvoidsRAMBackedTemp(t *testing.T) {
	got := SpillTmpDirFor("/mnt/x/data/beans/2016/Edmund/reference_genomes/v2.1/bams")
	if isRAMBacked(os.TempDir()) && strings.HasPrefix(got, os.TempDir()+string(os.PathSeparator)) {
		t.Errorf("spill dir %s is on the RAM-backed %s", got, os.TempDir())
	}
	if isRAMBacked(got) {
		t.Errorf("spill dir %s is RAM-backed", got)
	}
	t.Logf("spill dir: %s (os.TempDir()=%s, RAM-backed=%v)", got, os.TempDir(), isRAMBacked(os.TempDir()))
}

func TestSpillTmpDirHonoursExplicitEnv(t *testing.T) {
	base := t.TempDir()
	t.Setenv("GENOME_WHISPERER_TMPDIR", base)
	got := SpillTmpDirFor("/data/Edmund/reference_genomes/v2.1/bams")
	if !strings.HasPrefix(got, filepath.Join(base, "genome-whisperer")+string(os.PathSeparator)) {
		t.Errorf("GENOME_WHISPERER_TMPDIR ignored: %s", got)
	}
}

// Cleanup removes both scratch locations, so they must not be the same path and
// each must be findable without having been created.
func TestSpillAndWorkScratchAreDistinct(t *testing.T) {
	dir := "/mnt/x/data/beans/2016/Edmund/reference_genomes/v2.1/bams"
	if SpillTmpDirFor(dir) == WorkTmpDirFor(dir) {
		t.Errorf("spill and work scratch are the same dir: %s", WorkTmpDirFor(dir))
	}
}

func TestIsRAMBackedReadsProcMounts(t *testing.T) {
	if _, err := os.Stat("/proc/mounts"); err != nil {
		t.Skip("no /proc/mounts")
	}
	if isRAMBacked("/proc/self/no/such/path/under/root") {
		t.Error("/ reported as RAM-backed")
	}
}

func TestTmpBaseHonoursOverride(t *testing.T) {
	if got, want := TmpBase(), os.TempDir(); got != want {
		t.Errorf("TmpBase with no override = %s, want %s", got, want)
	}
	base := t.TempDir()
	t.Setenv("GENOME_WHISPERER_TMPDIR", base)
	if got := TmpBase(); got != base {
		t.Errorf("TmpBase = %s, want %s", got, base)
	}
	if got := WorkTmpDirFor("/data/x/bams"); !strings.HasPrefix(got, base+string(os.PathSeparator)) {
		t.Errorf("work scratch ignored the override: %s", got)
	}
}
