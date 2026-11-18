package annotation

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/gmaffy/genome-whisperer/utils"
)

// snpEffPaths locates the snpEff install from the snpEff launcher on PATH and
// returns the install root and its config file.
func snpEffPaths() (snpEffDir, configPath string, err error) {
	snpEffPath, err := exec.LookPath("snpEff")
	if err != nil {
		return "", "", fmt.Errorf("snpEff not found on PATH: %w", err)
	}
	snpEffDir = filepath.Dir(filepath.Dir(snpEffPath))
	configPath = filepath.Join(snpEffDir, "snpEff.config")
	if _, sErr := os.Stat(configPath); sErr != nil {
		return "", "", fmt.Errorf("snpEff config not found at %s: %w", configPath, sErr)
	}
	return snpEffDir, configPath, nil
}

// dbIsBuilt reports whether a database has actually been built, rather than
// merely declared in snpEff.config. A failed build leaves the config entry
// behind, so trusting "snpEff databases" alone would make a retry report
// success without building anything.
func dbIsBuilt(snpEffDir, db string) bool {
	_, err := os.Stat(filepath.Join(snpEffDir, "data", db, "snpEffectPredictor.bin"))
	return err == nil
}

// openMaybeGz opens a file and transparently decompresses it when gzipped.
func openMaybeGz(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if !strings.HasSuffix(path, ".gz") {
		return f, nil
	}
	zr, err := gzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("reading gzip %s: %w", path, err)
	}
	return struct {
		io.Reader
		io.Closer
	}{zr, multiCloser{zr, f}}, nil
}

type multiCloser []io.Closer

func (m multiCloser) Close() error {
	var first error
	for _, c := range m {
		if err := c.Close(); err != nil && first == nil {
			first = err
		}
	}
	return first
}

// transcriptIDs maps every name a GFF3 transcript is known by to its ID. The ID
// is what snpEff uses as the transcript identifier, so it is what the CDS and
// protein FASTA names have to agree with. GWH (CNCB Genome Warehouse) GFF3s
// also carry the accession their FASTA files are named by in Accession=, so
// that is mapped back to the ID as well.
func transcriptIDs(gff string) (map[string]string, error) {
	rc, err := openMaybeGz(gff)
	if err != nil {
		return nil, err
	}
	defer rc.Close()

	ids := make(map[string]string)
	sc := bufio.NewScanner(rc)
	sc.Buffer(make([]byte, 0, 1024*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 9 {
			continue
		}
		switch cols[2] {
		case "mRNA", "transcript":
		default:
			continue
		}
		var id, acc string
		for _, attr := range strings.Split(cols[8], ";") {
			attr = strings.TrimSpace(attr)
			if v, ok := strings.CutPrefix(attr, "ID="); ok {
				id = v
			} else if v, ok := strings.CutPrefix(attr, "Accession="); ok {
				acc = v
			}
		}
		if id == "" {
			continue
		}
		ids[id] = id
		if _, taken := ids[acc]; acc != "" && !taken {
			ids[acc] = id
		}
	}
	return ids, sc.Err()
}

// stageFasta copies a CDS or protein FASTA into the snpEff data directory,
// reconciling its record names with the GFF3 transcript IDs.
//
// Phytozome ships GFF3 IDs like Phvul.001G000400.1.v2.1 but names the matching
// FASTA record Phvul.001G000400.1, with the full ID only in an ID= key on the
// header. GWH uses its own accessions in the GFF3 Accession= and as FASTA
// record names, with the GFF3 ID in OriID= (CDS) or OriTrascriptID= (protein),
// and names protein records by protein accession. snpEff matches on the record
// name, so left alone every transcript fails its CDS and protein check and the
// build aborts. When a record name is already a GFF3 ID this is a byte-for-byte
// copy.
func stageFasta(src, dst string, txIDs map[string]string) (renamed int, err error) {
	rc, err := openMaybeGz(src)
	if err != nil {
		return 0, fmt.Errorf("opening %s: %w", src, err)
	}
	defer rc.Close()

	out, err := os.Create(dst)
	if err != nil {
		return 0, fmt.Errorf("creating %s: %w", dst, err)
	}
	defer out.Close()

	var w io.Writer = out
	if strings.HasSuffix(dst, ".gz") {
		zw := gzip.NewWriter(out)
		defer zw.Close()
		w = zw
	}
	bw := bufio.NewWriterSize(w, 1<<20)

	sc := bufio.NewScanner(rc)
	sc.Buffer(make([]byte, 0, 1024*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 || line[0] != '>' {
			if _, wErr := bw.WriteString(line + "\n"); wErr != nil {
				return renamed, wErr
			}
			continue
		}
		if id, ok := headerTranscriptID(line[1:], txIDs); ok {
			line = ">" + id
			renamed++
		}
		if _, wErr := bw.WriteString(line + "\n"); wErr != nil {
			return renamed, wErr
		}
	}
	if sErr := sc.Err(); sErr != nil {
		return renamed, sErr
	}
	return renamed, bw.Flush()
}

// headerTranscriptID resolves a FASTA header to the GFF3 transcript ID it
// belongs to. ok is false when the record name already is a GFF3 ID, or when
// nothing in the header resolves, and the header is left alone either way.
// The record name is tried first, then the value of each key=value field.
func headerTranscriptID(header string, txIDs map[string]string) (id string, ok bool) {
	fields := strings.Fields(header)
	if len(fields) == 0 {
		return "", false
	}
	if id, known := txIDs[fields[0]]; known {
		if id == fields[0] {
			return "", false
		}
		return id, true
	}
	for _, f := range fields[1:] {
		if _, v, found := strings.Cut(f, "="); found {
			if id, known := txIDs[v]; known {
				return id, true
			}
		}
	}
	return "", false
}

func CreateCustomDb(ref, prot, cds, species, gff, version string) error {
	snpEffDir, configPath, err := snpEffPaths()
	if err != nil {
		return err
	}

	db := fmt.Sprintf("%s%s", species, version)

	// ------------------------------ Already built? ------------------------------ //
	if dbIsBuilt(snpEffDir, db) {
		fmt.Printf("Database %s is already built. Nothing to do.\n\n", db)
		return nil
	}

	// --------------------------- Register in snpEff.config ----------------------- //
	// Skip when already present: a previous failed build leaves the entry behind
	// and appending again would give snpEff a duplicate genome.
	cfgBytes, err := os.ReadFile(configPath)
	if err != nil {
		return fmt.Errorf("reading snpEff config %s: %w", configPath, err)
	}
	entry := fmt.Sprintf("%s.genome", db)
	if !strings.Contains(string(cfgBytes), entry) {
		fmt.Printf("Adding %s to %s ...\n\n", entry, configPath)
		f, oErr := os.OpenFile(configPath, os.O_APPEND|os.O_WRONLY, 0644)
		if oErr != nil {
			return fmt.Errorf("opening snpEff config for append: %w", oErr)
		}
		_, wErr := fmt.Fprintf(f, "\n# %s genome, version %s\n%s : %s\n", species, db, entry, species)
		cErr := f.Close()
		if wErr != nil {
			return fmt.Errorf("writing snpEff config: %w", wErr)
		}
		if cErr != nil {
			return fmt.Errorf("closing snpEff config: %w", cErr)
		}
	} else {
		fmt.Printf("%s already declared in %s\n\n", entry, configPath)
	}

	// -------------------------------- Data directory ----------------------------- //
	dbDir := filepath.Join(snpEffDir, "data", db)
	if mErr := os.MkdirAll(dbDir, 0755); mErr != nil {
		return fmt.Errorf("creating %s: %w", dbDir, mErr)
	}

	// ---------------------------------- Reference -------------------------------- //
	seqs := filepath.Join(dbDir, "sequences.fa")
	if strings.HasSuffix(ref, ".gz") {
		seqs += ".gz"
	}
	fmt.Printf("Copying %s -> %s ...\n", ref, seqs)
	if cErr := CopyFile(ref, seqs); cErr != nil {
		return cErr
	}

	// ---------------------------------- Annotation ------------------------------- //
	// snpEff reads GFF3 natively. Converting to GTF with gffread first drops the
	// gene records and the gene_name attribute, which makes snpEff invent one
	// synthetic gene per transcript named null.N and destroys gene-level
	// identity in every downstream report.
	var buildFmt string
	switch {
	case strings.HasSuffix(gff, ".gff3.gz"), strings.HasSuffix(gff, ".gff.gz"):
		buildFmt = "-gff3"
		if cErr := CopyFile(gff, filepath.Join(dbDir, "genes.gff.gz")); cErr != nil {
			return cErr
		}
	case strings.HasSuffix(gff, ".gff3"), strings.HasSuffix(gff, ".gff"):
		buildFmt = "-gff3"
		if cErr := CopyFile(gff, filepath.Join(dbDir, "genes.gff")); cErr != nil {
			return cErr
		}
	case strings.HasSuffix(gff, ".gtf.gz"):
		buildFmt = "-gtf22"
		if cErr := CopyFile(gff, filepath.Join(dbDir, "genes.gtf.gz")); cErr != nil {
			return cErr
		}
	case strings.HasSuffix(gff, ".gtf"):
		buildFmt = "-gtf22"
		if cErr := CopyFile(gff, filepath.Join(dbDir, "genes.gtf")); cErr != nil {
			return cErr
		}
	default:
		return fmt.Errorf("unrecognised annotation format %q: expected .gff3, .gff, .gtf (optionally .gz)", gff)
	}
	fmt.Printf("Copied %s -> %s (building with %s)\n\n", gff, dbDir, buildFmt)

	// ------------------------- CDS and protein, ID-reconciled -------------------- //
	var txIDs map[string]string
	if buildFmt == "-gff3" {
		txIDs, err = transcriptIDs(gff)
		if err != nil {
			return fmt.Errorf("reading transcript IDs from %s: %w", gff, err)
		}
		n := 0
		for name, id := range txIDs {
			if name == id {
				n++
			}
		}
		fmt.Printf("Found %d transcript IDs in %s\n", n, gff)
	}

	for _, f := range []struct{ src, name string }{{prot, "protein.fa"}, {cds, "cds.fa"}} {
		dst := filepath.Join(dbDir, f.name)
		if strings.HasSuffix(f.src, ".gz") {
			dst += ".gz"
		}
		renamed, sErr := stageFasta(f.src, dst, txIDs)
		if sErr != nil {
			return fmt.Errorf("staging %s: %w", f.src, sErr)
		}
		if renamed > 0 {
			fmt.Printf("Copied %s -> %s (%d records renamed to their GFF3 ID)\n", f.src, dst, renamed)
		} else {
			fmt.Printf("Copied %s -> %s\n", f.src, dst)
		}
	}

	// ----------------------------------- Build ----------------------------------- //
	cmdStr := fmt.Sprintf("snpEff build %s -v %s", buildFmt, db)
	fmt.Printf("\nRunning: %s\n\n", cmdStr)
	if bErr := utils.RunBashCmdVerbose(cmdStr); bErr != nil {
		return fmt.Errorf("snpEff build failed for %s: %w", db, bErr)
	}

	if !dbIsBuilt(snpEffDir, db) {
		return fmt.Errorf("snpEff build reported success but %s was not written",
			filepath.Join(dbDir, "snpEffectPredictor.bin"))
	}

	fmt.Printf("\nDatabase %s built successfully.\n", db)
	return nil
}

func CopyFile(src, dst string) error {
	sourceFile, sErr := os.Open(src)
	if sErr != nil {
		return fmt.Errorf("couldn't open source file %s: %w", src, sErr)
	}
	defer sourceFile.Close()

	dstFile, dErr := os.Create(dst)
	if dErr != nil {
		return fmt.Errorf("couldn't create destination file %s: %w", dst, dErr)
	}
	defer dstFile.Close()

	if _, err := io.Copy(dstFile, sourceFile); err != nil {
		return fmt.Errorf("failed to copy file contents: %w", err)
	}
	return dstFile.Close()
}

func CreateCustomDbFromConfig(configFile, species, version string) error {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		return fmt.Errorf("reading config %s: %w", configFile, err)
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Proteins:", cfg.Proteins)
	fmt.Println("CDS:", cfg.CDS)
	fmt.Println("GFF:", cfg.GFF)

	for _, f := range []struct{ label, path string }{
		{"reference", cfg.Reference},
		{"protein", cfg.Proteins},
		{"CDS", cfg.CDS},
		{"GFF", cfg.GFF},
	} {
		if f.path == "" {
			return fmt.Errorf("config %s does not set a %s file", configFile, f.label)
		}
		if _, sErr := os.Stat(f.path); sErr != nil {
			return fmt.Errorf("%s file %s is not readable: %w", f.label, f.path, sErr)
		}
	}
	if species == "" {
		return fmt.Errorf("please provide a species name with --species")
	}
	if version == "" {
		return fmt.Errorf("please provide the annotation version with --annotation-version")
	}

	return CreateCustomDb(cfg.Reference, cfg.Proteins, cfg.CDS, species, cfg.GFF, version)
}
