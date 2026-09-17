//go:build unix

package warehouse

import "golang.org/x/sys/unix"

// volumeCapacity reports a filesystem's size via statfs.
//
// The block size must come from the mount, not a constant: on this estate the
// data volumes report 4096-byte blocks and the others 1024, so a hardcoded
// 4096 overstates the smaller ones fourfold.
//
// A failure returns ok=false rather than zeros, because a zero-byte disk is a
// plausible-looking lie and "we could not ask" is the truth.
func volumeCapacity(path string) (total, free, avail, blockSize int64, ok bool) {
	var st unix.Statfs_t
	if err := unix.Statfs(path, &st); err != nil {
		return 0, 0, 0, 0, false
	}

	// Frsize is the fragment size and the one statvfs multiplies by; Bsize is
	// the preferred I/O size. They agree on every mount here, but Frsize is
	// the correct choice and Bsize is the fallback when a filesystem leaves it
	// unset.
	blockSize = int64(st.Frsize)
	if blockSize == 0 {
		blockSize = int64(st.Bsize)
	}
	if blockSize == 0 {
		return 0, 0, 0, 0, false
	}

	return int64(st.Blocks) * blockSize,
		int64(st.Bfree) * blockSize,
		int64(st.Bavail) * blockSize,
		blockSize,
		true
}
