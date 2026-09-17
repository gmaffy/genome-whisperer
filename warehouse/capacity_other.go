//go:build !unix

package warehouse

// volumeCapacity has no portable implementation off unix. Reporting ok=false
// keeps the scan working and the capacity fields explicitly unknown.
func volumeCapacity(string) (total, free, avail, blockSize int64, ok bool) {
	return 0, 0, 0, 0, false
}
