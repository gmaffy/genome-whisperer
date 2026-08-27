package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

// indexCandidates lists every index name an alignment file might carry, current
// and historical: .csi is what BamIndex writes now, .bai what earlier runs left
// behind, and both the "beside the file" and "instead of the extension" spellings
// have been in use.
func indexCandidates(path string) []string {
	stem := strings.TrimSuffix(path, filepath.Ext(path))
	if strings.HasSuffix(strings.ToLower(path), ".cram") {
		return []string{path + ".crai", stem + ".crai"}
	}
	return []string{path + ".csi", stem + ".csi", path + ".bai", stem + ".bai"}
}

// findIndex returns the first usable index for path, and whether one was found.
//
// Existence is not enough. An index left behind by an interrupted run is present
// and non-empty but truncated, and GATK given a truncated index does not
// obviously fail — it silently sees only part of the file. So the index is read
// here: .crai and .csi are BGZF, where a truncated stream is detectable
// directly, and a legacy .bai is checked through samtools, which reads the index
// without touching the alignment data.
func findIndex(path string, verbose bool) (string, int64, bool) {
	for _, candidate := range indexCandidates(path) {
		info, err := os.Stat(candidate)
		if err != nil || info.Size() == 0 {
			continue
		}

		var check string
		if strings.HasSuffix(candidate, ".bai") {
			check = fmt.Sprintf("samtools idxstats %s > /dev/null", shQuote(path))
		} else {
			check = fmt.Sprintf("gzip -t %s", shQuote(candidate))
		}
		if err := utils.RunCmd(check, verbose); err != nil {
			color.Yellow("Index %s is unreadable (%v) — treating it as missing\n", candidate, err)
			continue
		}

		// An index older than its data file may belong to a previous version of
		// it, but copying a directory can reorder timestamps on its own, and
		// rebuilding the index of a large CRAM costs an hour. So this is said
		// out loud and not acted on; the integrity check above is what decides.
		if dataInfo, dErr := os.Stat(path); dErr == nil && info.ModTime().Before(dataInfo.ModTime()) {
			color.Yellow("Index %s is older than %s — reusing it anyway, delete it to force a rebuild\n", candidate, filepath.Base(path))
		}

		return candidate, info.Size(), true
	}
	return "", 0, false
}

// setIndexInfo fills in the index fields of info from whatever index its file
// actually carries, so a scan reports the same thing the pipeline decided.
func setIndexInfo(info *FileInfo, verbose bool) {
	_, size, found := findIndex(info.Path, verbose)
	info.IndexPresent = found
	info.IndexSize = size
}
