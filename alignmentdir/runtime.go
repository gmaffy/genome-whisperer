package alignmentdir

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/fatih/color"

	"github.com/gmaffy/genome-whisperer/utils"
)

// runtimeOpts is the resource budget for one run, worked out once in
// RunAlignReadsDir so that every external tool is sized from the same figures
// instead of each call site guessing.
type runtimeOpts struct {
	threads    int    // per-sample thread budget: aligner -t, samtools -@
	javaOpts   string // --java-options for the one-JVM-per-sample GATK steps
	shardJava  string // --java-options for the scattered GATK steps
	shardCount int    // interval shards a whole-genome GATK step is split into
	pairHmm    int    // --native-pair-hmm-threads per HaplotypeCaller shard
}

// gatkSlots bounds how many GATK JVMs run at once across the whole process.
//
// Sample-level parallelism and shard-level fan-out multiply: five samples each
// scattering a HaplotypeCaller over twenty intervals is a hundred JVMs, and
// each one costs a few GB. Only leaf commands take a slot — a function that
// spawns shards must not hold one itself, or the shards it waits on can never
// acquire theirs.
var gatkSlots chan struct{}

func initGatkSlots(n int) {
	if n < 1 {
		n = 1
	}
	gatkSlots = make(chan struct{}, n)
}

// runGatk runs one GATK command, waiting for a JVM slot first.
//
// Without --verbose the command's output is captured rather than discarded and
// the tail replayed if it fails. A scattered step that dies otherwise reports
// nothing but an exit status, once per shard, which is no use at all when the
// failure arrives hours into a run.
func runGatk(cmdStr string, verbose bool) error {
	if gatkSlots != nil {
		gatkSlots <- struct{}{}
		defer func() { <-gatkSlots }()
	}
	fmt.Printf("\n-------------------------------------------------------------------\nRunning: %s\n------------------------------------------------------------------\n\n", cmdStr)

	if verbose {
		return runBash(cmdStr, true)
	}

	out, err := exec.Command("bash", "-c", cmdStr).CombinedOutput()
	if err != nil {
		color.Red("Command failed: %v\n%s\n", err, tailLines(string(out), 25))
	}
	return err
}

// tailLines returns the last n non-blank lines of s.
func tailLines(s string, n int) string {
	var lines []string
	for _, line := range strings.Split(s, "\n") {
		if strings.TrimSpace(line) != "" {
			lines = append(lines, line)
		}
	}
	if len(lines) > n {
		lines = lines[len(lines)-n:]
	}
	return strings.Join(lines, "\n")
}

// availableMemGB reports free RAM in GB from /proc/meminfo's MemAvailable, which
// unlike MemFree counts reclaimable page cache. It returns 0 when that cannot be
// read — a non-Linux host, or a restricted container — and callers fall back to
// a fixed figure.
func availableMemGB() int {
	data, err := os.ReadFile("/proc/meminfo")
	if err != nil {
		return 0
	}
	for _, line := range strings.Split(string(data), "\n") {
		if !strings.HasPrefix(line, "MemAvailable:") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return 0
		}
		kb, cErr := strconv.Atoi(fields[1])
		if cErr != nil {
			return 0
		}
		return kb / (1024 * 1024)
	}
	return 0
}

func clamp(v, lo, hi int) int {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}

// newRuntimeOpts sizes the run from the cores and free memory actually present.
//
// Heaps are deliberately two sizes: MarkDuplicates and ApplyBQSR get one JVM per
// sample and need room to spill, while a HaplotypeCaller shard is a streaming
// job that runs many-at-once and needs far less. Both are bounded by the JVM
// slot count rather than by the sample count, which is what keeps the two
// dimensions of parallelism from multiplying into an out-of-memory kill.
func newRuntimeOpts(threads, parallelJobs int) runtimeOpts {
	if threads < 1 {
		threads = 1
	}
	if parallelJobs < 1 {
		parallelJobs = 1
	}

	usable := availableMemGB() - 4
	if usable < 4 {
		usable = 4
	}

	slots := clamp(runtime.NumCPU()/2, 1, 8)
	initGatkSlots(slots)

	sampleHeap := clamp(usable/parallelJobs, 4, 24)
	shardHeap := clamp(usable/slots, 3, 8)

	// One shard per ~500 Mb of reference: fine enough to spread over the cores
	// and to make a killed run resume cheaply, coarse enough that JVM start-up
	// stays noise.
	shards := clamp(slots*3, 8, 64)

	return runtimeOpts{
		threads:    threads,
		javaOpts:   fmt.Sprintf("-Xms2g -Xmx%dg", sampleHeap),
		shardJava:  fmt.Sprintf("-Xms1g -Xmx%dg", shardHeap),
		shardCount: shards,
		pairHmm:    clamp(runtime.NumCPU()/slots, 1, 4),
	}
}

// scratchDirs lists the scratch directories a sample's scatter steps create, so
// cleanup can remove them once the final CRAMs exist.
func scratchDirs(bamDir string) []string {
	return []string{
		filepath.Join(bamDir, "shards"),
		utils.WorkTmpDirFor(bamDir),
		// The sort/MarkDuplicates spill lives on a different disk from the rest
		// of the scratch, so it is a second directory to remove, not the same
		// one under another name.
		utils.SpillTmpDirFor(bamDir),
	}
}
