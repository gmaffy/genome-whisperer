package annotation

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// PRG top hits: one best BLAST hit per query protein.
//
// This is the Go port of the three-layer Python stack — blast_prots_to_prgdb.sh
// running blastp, make_prg_top_hits.py orchestrating, top_prg_hits.py reducing
// with pandas. It produces the file AddPrg, CreateSuperVcf and
// genespace.genePrgMap all read, and which nothing in the Go CLI could write.
//
// The raw BLAST table is never written. The Python version kept it — 127 MB for
// one squash genome — and then discarded all but one row per query. Here blastp
// writes to a pipe and the reduction happens as the rows stream past, so the
// only thing that reaches disk is the answer.

// prgColumns is the -outfmt 6 field list, and prgHeader is the header written
// above the reduced table. The two must stay in step with each other and with
// the column names the readers look up by name — QseqID, PercIdent, Length,
// Qlen and Eval are required by both genespace.genePrgMap and AddPrg. Reordering
// either breaks every consumer silently, since they index by header position.
const (
	prgColumns = "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
	prgHeader  = "QseqID\tQlen\tSseqID\tSlen\tPercIdent\tLength\tMISMATCH\tGAPOPEN\tQstart\tQend\tSstart\tSend\tEval\tBitscore"
	prgFields  = 14
	// bitscoreIdx is the last field, which is why the reducer can find it with
	// LastIndexByte instead of splitting every row.
	bitscoreIdx = 13
)

// DefaultPrgEvalue is the cutoff the shell script used, spelled 0.000000001.
const DefaultPrgEvalue = "1e-9"

type PrgBlastOptions struct {
	Proteins          string // protein FASTA, blastp -query
	PrgDb             string // blastp -db prefix
	Species           string // with AnnotationVersion, names the output file
	AnnotationVersion string
	OutDir            string // default: the directory holding Proteins
	Evalue            string // default DefaultPrgEvalue
	Threads           int    // 0 omits -num_threads
	MaxTargetSeqs     int    // 0 omits -max_target_seqs; see the note below
	Verbose           bool
}

// PrgTopHitsName is the output file name for a species and annotation version:
// the two concatenated without a separator, upper-cased. The concatenation is
// the convention the Python raw file already used (Maximav1.1.PRGblast.txt) and
// the upper-casing matches the existing MOSCHATA_PRG_TOP_HITS.txt files.
func PrgTopHitsName(species, annotationVersion string) string {
	return strings.ToUpper(species+annotationVersion) + "_PRG_TOP_HITS.txt"
}

// BlastProtsToPrgDb blasts a protein FASTA against the PRG database and writes
// one best hit per protein beside the protein file. It returns the path written.
func BlastProtsToPrgDb(opts PrgBlastOptions) (string, error) {
	if err := utils.CheckDeps([]string{"blastp"}); err != nil {
		return "", fmt.Errorf("dependency check failed: %w", err)
	}

	// Everything cheap is checked first. blastp against the PRG database takes
	// long enough that finding out about an empty --species at the end of it is
	// a wasted afternoon.
	if opts.Proteins == "" {
		return "", fmt.Errorf("no protein FASTA given")
	}
	if _, err := os.Stat(opts.Proteins); err != nil {
		return "", fmt.Errorf("protein FASTA %s is not readable: %w", opts.Proteins, err)
	}
	if opts.PrgDb == "" {
		return "", fmt.Errorf("no PRG BLAST database given")
	}
	if strings.TrimSpace(opts.Species) == "" {
		return "", fmt.Errorf("no species given; it names the output file")
	}

	evalue := opts.Evalue
	if evalue == "" {
		evalue = DefaultPrgEvalue
	}
	if _, err := strconv.ParseFloat(evalue, 64); err != nil {
		return "", fmt.Errorf("evalue %q is not a number: %w", evalue, err)
	}

	// Only a path-shaped -db can be checked on disk; a bare name is resolved by
	// blastp itself through $BLASTDB and is passed through untouched.
	if utils.LooksLikeBlastDBPath(opts.PrgDb) && !utils.HasProtBlastDB(opts.PrgDb) {
		return "", fmt.Errorf("no protein BLAST database beside %s (looked for .psq and .pal). Run: makeblastdb -in %s -dbtype prot -out %s",
			opts.PrgDb, opts.PrgDb, opts.PrgDb)
	}

	// The output belongs with the genome it describes, which by default means
	// the directory the protein FASTA came from.
	outDir := opts.OutDir
	if outDir == "" {
		outDir = filepath.Dir(opts.Proteins)
	}
	if err := utils.EnsureWritableDir(outDir); err != nil {
		return "", fmt.Errorf("output directory: %w", err)
	}
	outPath := filepath.Join(outDir, PrgTopHitsName(opts.Species, opts.AnnotationVersion))

	// blastp has no gzip support: handed a .gz it parses the compressed bytes as
	// FASTA, writes "CFastaReader: ... doesn't look like plausible data" to
	// stderr and exits 0 with no hits. Decompress to a temporary file first so a
	// .gz protein FASTA behaves like every other .gz input the tool accepts.
	queryPath := opts.Proteins
	if strings.HasSuffix(queryPath, ".gz") {
		tmp, dErr := decompressToTemp(queryPath)
		if dErr != nil {
			return "", dErr
		}
		defer os.Remove(tmp)
		utils.Printf("Decompressed %s to %s for blastp\n", opts.Proteins, tmp)
		queryPath = tmp
	}

	args := []string{
		"-query", queryPath,
		"-db", opts.PrgDb,
		"-evalue", evalue,
		// -max_hsps caps HSPs per query-subject pair. A subject's best HSP is
		// reported first, so this cannot change which subject wins on bitscore
		// — unlike -max_target_seqs, which caps subjects during the preliminary
		// stage, before the final gapped ranking, and so can drop the true best
		// hit (Shah et al. 2018). Since this command reports exactly one row per
		// query, that substitution would be invisible in the output. The real
		// table averages 94.5 HSPs per query and all but one are discarded, so
		// this is where the saving is anyway.
		"-max_hsps", "1",
		"-outfmt", prgColumns,
	}
	if opts.Threads > 0 {
		args = append(args, "-num_threads", strconv.Itoa(opts.Threads))
	}
	if opts.MaxTargetSeqs > 0 {
		args = append(args, "-max_target_seqs", strconv.Itoa(opts.MaxTargetSeqs))
	}

	cmd := exec.Command("blastp", args...)
	printable := shellQuoted(cmd.Args)
	utils.Printf("\n-------------------------------------------------------------------\n%s\n-------------------------------------------------------------------\n\n", printable)

	// stderr gets its own pipe rather than sharing stdout's. blastp writes
	// warnings there routinely — "Selenocysteine (U) at position N replaced by
	// X" — and a merged warning line would land in the middle of the table.
	// Non-empty stderr is therefore not a failure; only a non-zero exit is.
	stderrTail := &ringBuffer{limit: 40}
	if opts.Verbose {
		out := utils.LogWriter()
		defer out.Close()
		cmd.Stderr = io.MultiWriter(out, stderrTail)
	} else {
		cmd.Stderr = stderrTail
	}

	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return "", fmt.Errorf("blastp stdout pipe: %w", err)
	}
	if err := cmd.Start(); err != nil {
		return "", fmt.Errorf("starting blastp: %w", err)
	}

	// Order matters here. The pipe must be drained to EOF before Wait, which
	// closes it, and the reducer must never return early: blastp would block
	// writing into a full pipe and neither side would move again.
	rows, malformed, scanErr := topHits(stdout)
	waitErr := cmd.Wait()

	if waitErr != nil {
		return "", fmt.Errorf("blastp failed: %w\n  command: %s\n%s", waitErr, printable, stderrTail.report())
	}
	if scanErr != nil {
		return "", fmt.Errorf("reading blastp output: %w", scanErr)
	}
	if malformed > 0 {
		utils.Printf("Skipped %d blastp row(s) that did not have %d tab-separated fields\n", malformed, prgFields)
	}
	// blastp exiting 0 with nothing to show is what a wrong --db looks like.
	// Writing a header-only file here would push the emptiness downstream,
	// where it reads as "no gene matched the PRG" rather than "nothing ran".
	if len(rows) == 0 {
		// The stderr tail is the difference between "wrong database" and "blastp
		// could not read the query at all", which look identical from exit code
		// and stdout alone.
		return "", fmt.Errorf("blastp returned no hits for %s against %s; check that --prg-db is the right protein database and that the query is readable FASTA\n  command: %s\n%s",
			opts.Proteins, opts.PrgDb, printable, stderrTail.report())
	}

	if err := writeTopHits(outPath, rows); err != nil {
		return "", err
	}

	utils.Printf("Wrote %d top hits to %s\n", len(rows), outPath)
	return outPath, nil
}

// topHits reduces raw -outfmt 6 rows to the single best-scoring hit per query.
//
// It returns the surviving rows sorted by query id, how many rows were the wrong
// shape, and any read error. It never returns early on a bad row, because the
// caller is draining a pipe a live blastp is still writing into.
func topHits(r io.Reader) (rows []string, malformed int, err error) {
	// One line per query, not per HSP: ~32k proteins for a squash genome, so a
	// few MB against the 127 MB the Python version put on disk.
	best := make(map[string]string)
	bestScore := make(map[string]float64)

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 1<<20), 10<<20)

	for scanner.Scan() {
		// Text, not Bytes: Bytes aliases the scanner's buffer and would be
		// overwritten by the next Scan, long before these rows are written out.
		line := strings.TrimRight(scanner.Text(), "\r")
		if line == "" {
			continue
		}

		// A correct row has exactly 13 tabs. Checking the count is both the
		// shape test and what makes the two index lookups below safe.
		if strings.Count(line, "\t") != prgFields-1 {
			malformed++
			continue
		}
		// The query id is the first field and the bitscore the last, so both can
		// be taken by index without splitting the row into 14 strings — worth it
		// at ~94 HSPs per query across 32k queries.
		tab := strings.IndexByte(line, '\t')
		lastTab := strings.LastIndexByte(line, '\t')

		score, convErr := strconv.ParseFloat(line[lastTab+1:], 64)
		if convErr != nil {
			malformed++
			continue
		}

		qseqID := line[:tab]
		// Strictly greater, so a tie keeps the row seen first — what pandas'
		// .iloc[0] did. This is not a nicety: on the real Maxima table 8.4% of
		// queries have more than one hit at the maximum bitscore.
		if prev, seen := bestScore[qseqID]; !seen || score > prev {
			bestScore[qseqID] = score
			best[qseqID] = line
		}
	}
	if err := scanner.Err(); err != nil {
		// Swallowing this would truncate the table and still report success.
		return nil, malformed, err
	}

	ids := make([]string, 0, len(best))
	for id := range best {
		ids = append(ids, id)
	}
	// Sorted output, where the Python iterated a set and so produced a different
	// row order on every run.
	sort.Strings(ids)

	rows = make([]string, 0, len(ids))
	for _, id := range ids {
		rows = append(rows, best[id])
	}
	return rows, malformed, nil
}

// writeTopHits writes the header and rows through a .partial file, so an
// interrupted run cannot leave a half-written table under the real name.
func writeTopHits(outPath string, rows []string) error {
	partial := outPath + ".partial"
	f, err := os.Create(partial)
	if err != nil {
		return fmt.Errorf("creating %s: %w", partial, err)
	}

	w := bufio.NewWriter(f)
	writeErr := func() error {
		if _, err := w.WriteString(prgHeader + "\n"); err != nil {
			return err
		}
		for _, row := range rows {
			if _, err := w.WriteString(row + "\n"); err != nil {
				return err
			}
		}
		return w.Flush()
	}()

	if writeErr == nil {
		writeErr = f.Close()
	} else {
		f.Close()
	}
	if writeErr != nil {
		os.Remove(partial)
		return fmt.Errorf("writing %s: %w", partial, writeErr)
	}

	if err := os.Rename(partial, outPath); err != nil {
		os.Remove(partial)
		return fmt.Errorf("renaming %s to %s: %w", partial, outPath, err)
	}
	return nil
}

// shellQuoted renders an argv the way it would have to be typed at a prompt.
// Nothing here goes through a shell, but the banner and the failure message are
// most useful when they can be pasted into one, and -outfmt's value contains
// spaces.
func shellQuoted(argv []string) string {
	quoted := make([]string, 0, len(argv))
	for _, a := range argv {
		if strings.ContainsAny(a, " \t\"'\\") {
			quoted = append(quoted, "'"+strings.ReplaceAll(a, "'", `'\''`)+"'")
			continue
		}
		quoted = append(quoted, a)
	}
	return strings.Join(quoted, " ")
}

// ringBuffer keeps the last few lines a command wrote to stderr, for the error
// message. utils has one of these in tailBuf, but it is unexported and this
// package cannot reach it.
type ringBuffer struct {
	limit int
	lines []string
	rest  string
}

func (b *ringBuffer) Write(p []byte) (int, error) {
	b.rest += string(p)
	for {
		i := strings.IndexByte(b.rest, '\n')
		if i < 0 {
			break
		}
		b.add(b.rest[:i])
		b.rest = b.rest[i+1:]
	}
	return len(p), nil
}

func (b *ringBuffer) add(line string) {
	if strings.TrimSpace(line) == "" {
		return
	}
	b.lines = append(b.lines, line)
	if len(b.lines) > b.limit {
		b.lines = b.lines[len(b.lines)-b.limit:]
	}
}

func (b *ringBuffer) report() string {
	b.add(strings.TrimRight(b.rest, "\r\n"))
	b.rest = ""
	if len(b.lines) == 0 {
		return ""
	}
	var sb strings.Builder
	fmt.Fprintf(&sb, "  last %d line(s) of blastp output:\n", len(b.lines))
	for _, line := range b.lines {
		fmt.Fprintf(&sb, "    %s\n", line)
	}
	return sb.String()
}

// decompressToTemp writes a gzipped file out uncompressed and returns the
// temporary path, which the caller owns and must remove.
func decompressToTemp(src string) (string, error) {
	rc, err := openMaybeGz(src)
	if err != nil {
		return "", fmt.Errorf("opening %s: %w", src, err)
	}
	defer rc.Close()

	base := strings.TrimSuffix(filepath.Base(src), ".gz")
	tmp, err := os.CreateTemp("", "gw-*-"+base)
	if err != nil {
		return "", fmt.Errorf("creating temp file for %s: %w", src, err)
	}
	if _, err := io.Copy(tmp, rc); err != nil {
		tmp.Close()
		os.Remove(tmp.Name())
		return "", fmt.Errorf("decompressing %s: %w", src, err)
	}
	if err := tmp.Close(); err != nil {
		os.Remove(tmp.Name())
		return "", fmt.Errorf("closing %s: %w", tmp.Name(), err)
	}
	return tmp.Name(), nil
}
