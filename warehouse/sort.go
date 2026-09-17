package warehouse

import "sort"

// Every slice that reaches the report is sorted, so two scans of an unchanged
// tree produce byte-identical output and a diff between runs shows only real
// change.

func sortStrings(s []string) { sort.Strings(s) }

func sortByLengthDesc(s []string) {
	sort.Slice(s, func(i, j int) bool {
		if len(s[i]) != len(s[j]) {
			return len(s[i]) > len(s[j])
		}
		return s[i] < s[j]
	})
}

// uniqueSorted returns the distinct values of s, sorted.
func uniqueSorted(s []string) []string {
	if len(s) == 0 {
		return nil
	}
	seen := make(map[string]bool, len(s))
	out := make([]string, 0, len(s))
	for _, v := range s {
		if !seen[v] {
			seen[v] = true
			out = append(out, v)
		}
	}
	sort.Strings(out)
	return out
}
