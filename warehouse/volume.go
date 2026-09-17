package warehouse

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"syscall"
)

// dataDirName is the directory holding crops, matched case-insensitively.
// Spelling varies across this estate and the mounts themselves are
// case-insensitive, so the name is discovered rather than assumed.
const dataDirName = "data"

// resolvedVolume is one distinct filesystem to scan.
type resolvedVolume struct {
	// Name is a short label for the volume, used to join rows together.
	Name string
	// Roots is every requested path that resolved here, sorted.
	Roots []string
	// DataDir is the crop-holding directory, in its verbatim on-disk spelling.
	DataDir string
	// Spelling is that directory's name as the filesystem spells it.
	Spelling string

	dev uint64
	ino uint64
}

// resolveVolumes turns requested roots into distinct filesystems to scan.
//
// Two problems make this more than a loop. The mounts are case-insensitive, so
// /mnt/u/data, /mnt/u/DATA and /mnt/u/Data are one directory under three
// names — globbing the spellings would scan the same tree three times and
// treble every byte count. And identity has to key on the device and inode
// *pair*: on this estate /mnt/u and /mnt/v happen to share an inode number,
// so deduplicating on inode alone silently discards a whole volume.
//
// Any root that cannot be resolved is returned as a problem rather than an
// error, so one unmounted volume does not abandon the scan.
func resolveVolumes(roots []string) ([]resolvedVolume, []Nonconformance) {
	var problems []Nonconformance
	// Keyed by device and inode together.
	type fsKey struct{ dev, ino uint64 }
	index := map[fsKey]*resolvedVolume{}
	var order []fsKey

	for _, root := range roots {
		root = strings.TrimSpace(root)
		if root == "" {
			continue
		}
		abs, err := filepath.Abs(root)
		if err != nil {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleDataDirMissing,
				Severity: SeverityError, Path: root,
				Detail: fmt.Sprintf("resolving path: %v", err),
			})
			continue
		}

		dataDirs, dErr := findDataDirs(abs)
		if dErr != "" {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleDataDirMissing,
				Severity: SeverityError, Path: abs, Detail: dErr,
			})
			continue
		}
		if len(dataDirs) == 0 {
			problems = append(problems, Nonconformance{
				Kind: KindNonconformance, Rule: RuleDataDirMissing,
				Severity: SeverityError, Path: abs,
				Expected: dataDirName,
				Detail:   "no data directory here and this path holds no crops",
			})
			continue
		}

		for _, dataDir := range dataDirs {
			info, sErr := os.Stat(dataDir)
			if sErr != nil {
				problems = append(problems, Nonconformance{
					Kind: KindNonconformance, Rule: RuleDataDirMissing,
					Severity: SeverityError, Path: dataDir,
					Detail: sErr.Error(),
				})
				continue
			}

			dev, ino, ok := fsIdentity(info)
			if !ok {
				// Without an identity we cannot deduplicate, so treat the
				// path itself as the identity and carry on.
				dev, ino = 0, uint64(len(order)+1)
			}

			key := fsKey{dev, ino}
			if existing, seen := index[key]; seen {
				existing.Roots = append(existing.Roots, abs)
				continue
			}

			index[key] = &resolvedVolume{
				Name:     volumeName(dataDir),
				Roots:    []string{abs},
				DataDir:  dataDir,
				Spelling: filepath.Base(dataDir),
				dev:      dev,
				ino:      ino,
			}
			order = append(order, key)
		}
	}

	out := make([]resolvedVolume, 0, len(order))
	for _, key := range order {
		v := index[key]
		v.Roots = uniqueSorted(v.Roots)
		out = append(out, *v)
	}
	return out, problems
}

// findDataDirs returns the crop-holding directories under root.
//
// Normally there is exactly one, found by a case-insensitive name match over a
// single ReadDir — which yields the filesystem's own spelling and makes a
// {data,DATA} glob unnecessary. On a genuinely case-sensitive filesystem two
// differently-cased directories can both exist and both are returned; the
// device/inode check upstream collapses them if they turn out to be the same.
//
// If root has no data child but does look like a data directory itself, root is
// returned, so a caller may pass either a mount or the data directory.
func findDataDirs(root string) ([]string, string) {
	entries, err := os.ReadDir(root)
	if err != nil {
		return nil, err.Error()
	}

	var found []string
	for _, e := range entries {
		if !e.IsDir() {
			continue
		}
		if strings.EqualFold(e.Name(), dataDirName) {
			found = append(found, filepath.Join(root, e.Name()))
		}
	}
	if len(found) > 0 {
		sortStrings(found)
		return found, ""
	}

	// No data child. Treat root as the data directory if anything under it
	// looks like a crop — a directory holding a year or a joint VCF root.
	if looksLikeDataDir(root, entries) {
		return []string{root}, ""
	}
	return nil, ""
}

// looksLikeDataDir reports whether a directory holds crops, judged by finding
// one child that in turn holds a year directory or a joint VCF root.
func looksLikeDataDir(root string, entries []os.DirEntry) bool {
	for _, e := range entries {
		if !e.IsDir() || IsPruned(e.Name()) {
			continue
		}
		children, err := os.ReadDir(filepath.Join(root, e.Name()))
		if err != nil {
			continue
		}
		for _, c := range children {
			if !c.IsDir() {
				continue
			}
			if IsYear(c.Name()) || c.Name() == DirMergedVcfs || c.Name() == DirLegacyVcfs {
				return true
			}
		}
	}
	return false
}

// fsIdentity extracts the device and inode of a directory.
func fsIdentity(info os.FileInfo) (dev, ino uint64, ok bool) {
	st, cast := info.Sys().(*syscall.Stat_t)
	if !cast {
		return 0, 0, false
	}
	return uint64(st.Dev), uint64(st.Ino), true
}

// volumeName is a short label for a volume, taken from the mount point so rows
// read as "u" and "y" rather than repeating a full path.
//
// Mount letters are run-scoped: this estate has already moved a crop from one
// letter to another, so a consumer must key on the logical crop/year/sample
// tuple and treat this as a label.
func volumeName(dataDir string) string {
	parent := filepath.Dir(dataDir)
	if name := filepath.Base(parent); name != "" && name != "." && name != string(filepath.Separator) {
		return name
	}
	return filepath.Base(dataDir)
}

// capacityFor reads a volume's capacity, returning a problem row when the
// filesystem will not say.
func capacityFor(v resolvedVolume) (Volume, *Nonconformance) {
	row := Volume{
		Kind:     KindVolume,
		Volume:   v.Name,
		FSDev:    v.dev,
		Roots:    v.Roots,
		DataDir:  v.DataDir,
		Spelling: v.Spelling,
		Present:  true,
	}

	total, free, avail, blockSize, ok := volumeCapacity(v.DataDir)
	if !ok {
		return row, &Nonconformance{
			Kind: KindNonconformance, Rule: RuleCapacityUnavailable,
			Severity: SeverityWarning, Volume: v.Name, Path: v.DataDir,
			Detail: "statfs did not answer; capacity is unknown, not zero",
		}
	}

	row.CapacityKnown = true
	row.BytesTotal = total
	row.BytesFree = free
	row.BytesAvailable = avail
	row.BlockSize = blockSize
	return row, nil
}
