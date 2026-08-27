package alignmentdir

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/alignment"
	"github.com/gmaffy/genome-whisperer/utils"
)

// scatterGenomeBases is the reference size at which the whole-genome GATK steps
// switch from one process each to an interval scatter.
//
// Below it, a single HaplotypeCaller or BaseRecalibrator pass finishes in a time
// that the scatter's overhead — a JVM and a sequence-dictionary read per shard,
// then a merge — would be a noticeable fraction of. Above it, that pass is the
// longest thing in the run by a wide margin and splitting it is the only way to
// use more than one core on it. A gigabase sits between the reference sizes this
// pipeline has always handled a process at a time and the ones where a single
// pass runs for days.
const scatterGenomeBases = 1_000_000_000

// genomeBases reports a reference's total length, preferring the .fai because
// summing one column of it costs nothing. If the index cannot be read it falls
// back to the size of the FASTA on disk, which overstates the true length by the
// headers and line breaks — close enough for a decision with this much slack in
// it, and it keeps a missing .fai from deciding the route by accident.
func genomeBases(refFasta string) int64 {
	if entries, err := readFai(refFasta); err == nil {
		var total int64
		for _, e := range entries {
			total += e.length
		}
		return total
	}
	if info, err := os.Stat(refFasta); err == nil {
		return info.Size()
	}
	return 0
}

// useScatter reports whether a reference is large enough that the scattered
// implementations in this file should be used in place of the single-process
// ones in dirAlign.go.
func useScatter(refFasta string) bool {
	return genomeBases(refFasta) >= scatterGenomeBases
}

// baiMaxPosition is the highest coordinate a BAI index can address. Its linear
// index is a fixed 2^15 windows of 2^14 bases, so nothing past 2^29-1 fits.
const baiMaxPosition = 1<<29 - 1

// gatkReadsNeedBam reports whether GATK has to be handed a BAM rather than a
// CRAM for this reference.
//
// htsjdk has no CRAM index reader. Whenever a CRAM carries a .crai it converts
// it into an in-memory BAI as it opens the file
// (CRAIIndex.openCraiFileAsBaiStream), and a BAI cannot address a position past
// baiMaxPosition. So on a reference with a contig longer than that, every read
// beyond 512 Mb falls outside the linear index being built and the open dies
// with
//
//	ArrayIndexOutOfBoundsException: Index 36621 out of bounds for length 32770
//
// before the tool has read a single record — which is every shard of a
// scattered step, not an unlucky one. The .crai itself is sound and samtools
// uses it happily; the ceiling belongs to the index htsjdk builds from it. A
// BAM indexed as CSI has no such limit, so on these references the pipeline
// keeps an rgmd.bam alongside the cram and gives GATK that instead.
//
// A reference whose .fai cannot be read reports false: the scatter path needs
// the .fai for its intervals anyway and fails with a better message than this
// one would.
func gatkReadsNeedBam(refFasta string) bool {
	return refExceedsBaiLimit(refFasta)
}

// refExceedsBaiLimit reports whether any contig runs past what the BAI-style
// binning indexes can address.
//
// A reference whose .fai cannot be read reports false: the scatter path needs
// the .fai for its intervals anyway and fails with a better message than this
// one would.
func refExceedsBaiLimit(refFasta string) bool {
	entries, err := readFai(refFasta)
	if err != nil {
		return false
	}
	for _, e := range entries {
		if e.length > baiMaxPosition {
			return true
		}
	}
	return false
}

// variantSuffixes returns the extension a VCF written for this reference should
// carry, and the index that belongs beside it.
//
// A .tbi is the same binning index as a BAI and stops at the same ceiling.
// bcftools refuses to write one past it ("cannot be stored in a tbi index. Try
// using a csi index") and GATK gets all the way through a call before dying in
// close(), building one:
//
//	TabixIndexCreator.finalizeIndex → BinningIndexBuilder.processFeature
//	ArrayIndexOutOfBoundsException: Index 36621 out of bounds for length 32770
//
// A .csi has no such limit, but GATK cannot read one for a VCF at all — it
// reports "An index is required but was not found" and stops. That leaves the
// uncompressed VCF, which GATK indexes with a Tribble .idx that addresses the
// whole contig, so that is what the bootstrap chain writes on these references.
// It costs disk over bgzip; the shard VCFs and everything derived from them are
// intermediates that the sample cleanup removes.
func variantSuffixes(refFasta string) (vcfExt, idxExt string) {
	if refExceedsBaiLimit(refFasta) {
		return ".vcf", ".idx"
	}
	return ".vcf.gz", ".tbi"
}

// describeRoute names the route taken for a reference, for the run log.
func describeRoute(refFasta string) string {
	bases := genomeBases(refFasta)
	if bases >= scatterGenomeBases {
		return fmt.Sprintf("scattered (reference is %.2f Gb)", float64(bases)/1e9)
	}
	return fmt.Sprintf("single-process (reference is %.2f Gb)", float64(bases)/1e9)
}

// faiEntry is one line of a .fai: a contig and its length.
type faiEntry struct {
	name   string
	length int64
}

func readFai(refFasta string) ([]faiEntry, error) {
	faiPath := refFasta + ".fai"
	data, err := os.ReadFile(faiPath)
	if err != nil {
		return nil, fmt.Errorf("reading %s: %w", faiPath, err)
	}

	var entries []faiEntry
	for _, line := range strings.Split(string(data), "\n") {
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		length, cErr := strconv.ParseInt(fields[1], 10, 64)
		if cErr != nil {
			return nil, fmt.Errorf("bad length %q for contig %s in %s", fields[1], fields[0], faiPath)
		}
		entries = append(entries, faiEntry{name: fields[0], length: length})
	}
	if len(entries) == 0 {
		return nil, fmt.Errorf("%s lists no contigs", faiPath)
	}
	return entries, nil
}

// minShardBases is the smallest amount of reference a shard is worth starting a
// JVM for. Below roughly this much sequence, the tens of seconds GATK spends
// starting up and reading the sequence dictionary outweigh the work the shard
// does, so a small genome is scattered less — or, if it is small enough, not at
// all — rather than being cut into pieces that cost more to launch than to run.
const minShardBases = 20_000_000

// writeShardIntervals splits the reference into at most wanted roughly equal
// interval lists under dir and returns their paths in reference order.
//
// Contigs longer than the target are cut into sub-intervals rather than left
// whole: on a genome delivered as ~1 Gb pieces, one-shard-per-contig would leave
// most cores idle waiting for the biggest piece. The small-scaffold tail is
// packed together so it does not become thousands of one-contig jobs.
func writeShardIntervals(refFasta, dir string, wanted int) ([]string, error) {
	entries, err := readFai(refFasta)
	if err != nil {
		return nil, err
	}
	if wanted < 1 {
		wanted = 1
	}

	var total int64
	for _, e := range entries {
		total += e.length
	}

	// Scale the shard count down to the size of the genome in hand, so the same
	// figure works for a 16 Gb reference and a 100 Mb one.
	if maxUseful := int(total / minShardBases); wanted > maxUseful {
		wanted = maxUseful
	}
	if wanted < 1 {
		wanted = 1
	}

	target := total / int64(wanted)
	if target < 1 {
		target = total
	}

	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}

	var (
		shards  [][]string
		current []string
		size    int64
	)
	flush := func() {
		if len(current) > 0 {
			shards = append(shards, current)
			current = nil
			size = 0
		}
	}

	for _, e := range entries {
		// A contig whose name contains a colon cannot be written as
		// name:start-end without the interval parser mistaking part of the name
		// for coordinates, so it is kept whole however long it is.
		if e.length > target && !strings.Contains(e.name, ":") {
			flush()
			for start := int64(1); start <= e.length; start += target {
				end := start + target - 1
				if end > e.length {
					end = e.length
				}
				shards = append(shards, []string{fmt.Sprintf("%s:%d-%d", e.name, start, end)})
			}
			continue
		}
		if e.length >= target {
			// Too long to pack with others, but not splittable — give it a
			// shard of its own.
			flush()
			shards = append(shards, []string{e.name})
			continue
		}
		current = append(current, e.name)
		size += e.length
		if size >= target {
			flush()
		}
	}
	flush()

	paths := make([]string, 0, len(shards))
	for i, intervals := range shards {
		path := filepath.Join(dir, fmt.Sprintf("shard_%03d.list", i))
		if wErr := os.WriteFile(path, []byte(strings.Join(intervals, "\n")+"\n"), 0o644); wErr != nil {
			return nil, wErr
		}
		paths = append(paths, path)
	}
	return paths, nil
}

// runShards runs one command per shard concurrently and returns the first error.
// Concurrency comes from the JVM slot pool inside runGatk, so every shard is
// started here and the pool decides how many actually run.
func runShards(label string, cmds []string, verbose bool) error {
	var (
		wg       sync.WaitGroup
		mu       sync.Mutex
		firstErr error
		done     int
	)

	for i, cmdStr := range cmds {
		wg.Add(1)
		go func(idx int, cmdStr string) {
			defer wg.Done()
			err := runGatk(cmdStr, verbose)

			mu.Lock()
			defer mu.Unlock()
			done++
			if err != nil {
				color.Red("[%s] shard %d/%d failed: %v\n", label, idx, len(cmds), err)
				if firstErr == nil {
					firstErr = fmt.Errorf("%s shard %d: %w", label, idx, err)
				}
				return
			}
			color.Green("[%s] shard %d done (%d/%d)\n", label, idx, done, len(cmds))
		}(i, cmdStr)
	}

	wg.Wait()
	return firstErr
}

// scatterHaplotypeCaller calls variants over the whole genome as a scatter of
// per-interval jobs merged back into one VCF.
//
// HaplotypeCaller has no whole-file threading beyond its pair-HMM, so one
// process over a 15 Gb reference is the single longest step in a bootstrap BQSR
// run by a wide margin. Splitting by interval is how GATK's own pipelines
// parallelise it. Variant calls immediately at a shard boundary can differ from
// an unscattered run, since an active region is assembled from the reads inside
// its own shard.
func scatterHaplotypeCaller(refFasta, input, rawVCF, shardDir string, opts runtimeOpts, verbose bool) error {
	shards, err := writeShardIntervals(refFasta, shardDir, opts.shardCount)
	if err != nil {
		return err
	}
	vcfExt, idxExt := variantSuffixes(refFasta)
	color.Cyan("HaplotypeCaller: %d interval shards for %s\n", len(shards), filepath.Base(input))

	var (
		cmds     []string
		outputs  []string
		reusable int
	)
	for i, shard := range shards {
		shardVCF := filepath.Join(shardDir, fmt.Sprintf("shard_%03d.raw%s", i, vcfExt))
		outputs = append(outputs, shardVCF)

		// A shard already on disk is reused, so a run killed part-way through
		// resumes at the shard it died on rather than at the start.
		if info, sErr := os.Stat(shardVCF); sErr == nil && info.Size() > 0 {
			if vErr := utils.ValidateGvcf(shardVCF, false, true); vErr == nil {
				reusable++
				continue
			}
			_ = removeIfExists(shardVCF, shardVCF+idxExt)
		}

		cmds = append(cmds, fmt.Sprintf(`gatk --java-options "%s" HaplotypeCaller -R %s -I %s -L %s -O %s --native-pair-hmm-threads %d --tmp-dir %s`,
			opts.shardJava, shQuote(refFasta), shQuote(input), shQuote(shard), shQuote(shardVCF), opts.pairHmm, shQuote(alignment.WorkTmpDir(shardVCF))))
	}
	if reusable > 0 {
		color.Green("HaplotypeCaller: reusing %d/%d shards already on disk\n", reusable, len(shards))
	}

	if len(cmds) > 0 {
		if err := runShards("HaplotypeCaller", cmds, verbose); err != nil {
			return err
		}
	}

	var inputArgs []string
	for _, out := range outputs {
		inputArgs = append(inputArgs, "-I "+shQuote(out))
	}
	mergeCmd := fmt.Sprintf(`gatk --java-options "%s" MergeVcfs %s -O %s --TMP_DIR %s`,
		opts.javaOpts, strings.Join(inputArgs, " "), shQuote(rawVCF), shQuote(alignment.WorkTmpDir(rawVCF)))
	if err := runGatk(mergeCmd, verbose); err != nil {
		return fmt.Errorf("merging HaplotypeCaller shards: %w", err)
	}

	// The bootstrap steps downstream treat a missing index as an invalid VCF and
	// redo the call, so the index is made sure of here rather than assumed.
	if _, sErr := os.Stat(rawVCF + idxExt); sErr != nil {
		var iErr error
		if idxExt == ".tbi" {
			iErr = utils.RunCmd(fmt.Sprintf("bcftools index --tbi --force %s", shQuote(rawVCF)), verbose)
		} else {
			// bcftools cannot index an uncompressed VCF; the Tribble .idx that
			// GATK reads is GATK's own to write.
			iErr = runGatk(fmt.Sprintf(`gatk --java-options "%s" IndexFeatureFile -I %s`, opts.javaOpts, shQuote(rawVCF)), verbose)
		}
		if iErr != nil {
			return fmt.Errorf("indexing %s: %w", rawVCF, iErr)
		}
	}
	return nil
}

// scatterBaseRecalibrator builds a recalibration table as a scatter of
// per-interval tables gathered into one, which is exactly how GATK's own
// pipelines run this step. GatherBQSRReports sums the covariate counts, so the
// gathered table is identical to the unscattered one.
func scatterBaseRecalibrator(refFasta, input, recalTable, knownSitesArgs, shardDir, label string, opts runtimeOpts, verbose bool) error {
	shards, err := writeShardIntervals(refFasta, shardDir, opts.shardCount)
	if err != nil {
		return err
	}

	var (
		cmds    []string
		outputs []string
	)
	for i, shard := range shards {
		shardTable := filepath.Join(shardDir, fmt.Sprintf("%s_shard_%03d.recal.txt", label, i))
		_ = removeIfExists(shardTable)
		outputs = append(outputs, shardTable)
		cmds = append(cmds, fmt.Sprintf(`gatk --java-options "%s" BaseRecalibrator -R %s -I %s %s -L %s -O %s --maximum-cycle-value %d --tmp-dir %s`,
			opts.shardJava, shQuote(refFasta), shQuote(input), knownSitesArgs, shQuote(shard), shQuote(shardTable), alignment.MaxCycleValue, shQuote(alignment.WorkTmpDir(shardTable))))
	}

	if err := runShards("BaseRecalibrator/"+label, cmds, verbose); err != nil {
		return err
	}

	var inputArgs []string
	for _, out := range outputs {
		inputArgs = append(inputArgs, "-I "+shQuote(out))
	}
	gatherCmd := fmt.Sprintf(`gatk --java-options "%s" GatherBQSRReports %s -O %s --tmp-dir %s`,
		opts.javaOpts, strings.Join(inputArgs, " "), shQuote(recalTable), shQuote(alignment.WorkTmpDir(recalTable)))
	if err := runGatk(gatherCmd, verbose); err != nil {
		return fmt.Errorf("gathering %s recalibration tables: %w", label, err)
	}

	for _, out := range outputs {
		_ = removeIfExists(out)
	}
	return nil
}
