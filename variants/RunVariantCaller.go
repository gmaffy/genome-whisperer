package variants

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/genome-whisperer/utils"
)

type SeqInfo struct {
	ID  string
	Len int
}

func CreateGvcfGATK(bam string, refFile string, chroms []SeqInfo, species string, refVer string, logLevel string, verbose bool, out string) (string, error) {

	var regions string
	var gVCF string
	if len(chroms) == 1 {
		regions = chroms[0].ID
		chromDirName := strings.ReplaceAll(regions, ".", "_")
		chromDir := filepath.Join(out, chromDirName)
		gvcfDir := filepath.Join(chromDir, "gvcfs")
		gVCF = filepath.Join(gvcfDir, species+refVer+"."+regions+".g.vcf.gz")

		// ---------------------- Create Dirs ------------------------------- //
		for _, dir := range []string{chromDir, gvcfDir} {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				if cErr := os.MkdirAll(dir, 0755); cErr != nil {
					return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
				}
				//fmt.Printf("Created directory: %s\n", dir)
			}
		}
	} else {
		f, err := os.CreateTemp("", "gatk_intervals_contigs_*.list")
		if err != nil {
			return "", fmt.Errorf("creating GATK interval list: %w", err)
		}
		defer os.Remove(f.Name())
		defer f.Close()
		for _, c := range chroms {
			fmt.Fprintln(f, c.ID)
		}
		regions = f.Name()

		chromDir := filepath.Join(out, "contigs")
		gvcfDir := filepath.Join(chromDir, "gvcfs")
		// -------------------------- Create Directories ------------------------------------- //

		for _, dir := range []string{chromDir, gvcfDir} {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				if cErr := os.MkdirAll(dir, 0755); cErr != nil {
					return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
				}
			}
		}
		gVCF = filepath.Join(gvcfDir, species+refVer+"Contigs.g.vcf.gz")

	}

	hapCmdStr := fmt.Sprintf(
		`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`, refFile, bam, regions, gVCF, logLevel,
	)
	fmt.Printf("\n%s\n\n", hapCmdStr)
	return gVCF, runCmd(hapCmdStr, verbose)

}

func CreateGvcfDV(bam string, refFile string, chroms []SeqInfo, species, refVer, dvVer string, modelType string, verbose bool, out string) (string, error) {

	var regions string
	var gVCF string

	if len(chroms) == 1 {
		regions = chroms[0].ID
		chromDirName := strings.ReplaceAll(regions, ".", "_")
		chromDir := filepath.Join(out, chromDirName)
		gvcfDir := filepath.Join(chromDir, "gvcfs")
		// --------------------------------- Create directories -------------------------------- //
		for _, dir := range []string{chromDir, gvcfDir} {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				if cErr := os.MkdirAll(dir, 0755); cErr != nil {
					return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
				}
				//fmt.Printf("Created directory: %s\n", dir)
			}
		}

		gVCF = filepath.Join(gvcfDir, species+refVer+"."+regions+".g.vcf.gz")
	} else {

		chromDir := filepath.Join(out, "contigs")
		gvcfDir := filepath.Join(chromDir, "gvcfs")

		for _, dir := range []string{chromDir, gvcfDir} {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				if cErr := os.MkdirAll(dir, 0755); cErr != nil {
					return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
				}
				//fmt.Printf("Created directory: %s\n", dir)
			}
		}

		gVCF = filepath.Join(gvcfDir, species+refVer+"."+"Contigs.g.vcf.gz")

		f, err := os.CreateTemp(gvcfDir, "deepvariant_intervals_*.bed")
		if err != nil {
			return "", fmt.Errorf("creating DeepVariant interval list: %w", err)
		}
		defer os.Remove(f.Name())
		defer f.Close()
		for _, c := range chroms {
			fmt.Fprintf(f, "%s\t0\t%d\n", c.ID, c.Len)
		}
		regions = f.Name()

	}

	bamDir := filepath.Dir(bam)
	bamName := filepath.Base(bam)
	refDir := filepath.Dir(refFile)
	refName := filepath.Base(refFile)
	gvcfName := filepath.Base(gVCF)
	vcfName := strings.Replace(gvcfName, ".g.vcf.gz", ".vcf.gz", 1)

	safeRegion := strings.NewReplacer(string(os.PathSeparator), "_", ".", "_", "/", "_").Replace(regions)
	intermediateName := fmt.Sprintf("tmp_%s_%s", strings.TrimSuffix(bamName, filepath.Ext(bamName)), safeRegion)

	dvCmdStr := fmt.Sprintf(
		`docker run -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
			`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
			`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
			`--intermediate_results_dir /output/%s`,
		bamDir, refDir, out, dvVer,
		modelType, refName, bamName,
		regions, vcfName, gvcfName,
		intermediateName,
	)
	fmt.Printf("\n%s\n\n", dvCmdStr)
	return gVCF, runCmd(dvCmdStr, verbose)
}

func CreateSampleMap(gvcfs []string, outDir string) (string, error) {
	f, err := os.Create(filepath.Join(outDir, "sample_map.txt"))
	if err != nil {
		return "", fmt.Errorf("creating sample map: %w", err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	for _, sample := range gvcfs {
		if strings.HasSuffix(sample, ".g.vcf.gz") {

			out, err := exec.Command("bcftools", "query", "-l", sample).Output()
			if err != nil {
				return "", fmt.Errorf("%s: bcftools query -l: %w", sample, err)
			}
			names := strings.Fields(string(out))
			if len(names) != 1 {
				return "", fmt.Errorf("%s: expected 1 sample, got %d", sample, len(names))
			}
			_, err = w.WriteString(fmt.Sprintf("%s\t%s\n", names[0], sample))
			if err != nil {
				return "", fmt.Errorf("writing sample map: %w", err)
			}
		}
	}
	return f.Name(), nil
}

func MergeGvcfsGATK(gvcfs []string, refFile string, chrom string, species string, refVer string, verbose bool, logLevel string, outDir string) (string, error) {
	//outDir := filepath.Dir(gvcfDir)
	theDB := filepath.Join(outDir, chrom+"DB")
	tmpDir := filepath.Join(outDir, "tmp")
	tmpDir2 := filepath.Join(outDir, "tmp2")
	vcfDir := filepath.Join(outDir, "VCFs")

	for _, dir := range []string{tmpDir, tmpDir2, vcfDir} {
		if _, err := os.Stat(dir); os.IsNotExist(err) {
			if cErr := os.MkdirAll(dir, 0755); cErr != nil {
				return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
			}
			//fmt.Printf("Created directory: %s\n", dir)
		}
	}

	if rErr := os.RemoveAll(theDB); rErr != nil {
		return "", fmt.Errorf("removing %s: %w", theDB, rErr)
	}
	// --------------------------------------------- GenomicsDBImport -------------------------------------------- //
	sampleMap, err := CreateSampleMap(gvcfs, outDir)
	if err != nil {
		return "", fmt.Errorf("creating sample map: %w", err)
	}

	gDBCmd := fmt.Sprintf(
		`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport --sample-name-map %s `+
			`--genomicsdb-workspace-path %s --tmp-dir %s -L %s `+
			`--genomicsdb-shared-posixfs-optimizations true --batch-size 50 `+
			`--bypass-feature-reader --verbosity %s`,
		sampleMap, theDB, tmpDir, chrom, logLevel,
	)
	err = runCmd(gDBCmd, verbose)
	if err != nil {
		return "", fmt.Errorf("gatk GenomicsDBImport: %w", err)
	}
	// ---------------------------------------- GenotypeGvcfs ------------------------------------------------------- //
	jointVCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.vcf.gz")
	genoCmd := fmt.Sprintf(
		`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
		refFile, theDB, jointVCF, tmpDir2, logLevel,
	)
	err = runCmd(genoCmd, verbose)
	if err != nil {
		return "", fmt.Errorf("gatk GenotypeGVCFs: %w", err)
	}

	return jointVCF, err
}

func MergeGvcfsGlnexus(gvcfs []string, chrom string, species string, refVer string, caller string, verbose bool, outDir string) (string, error) {
	//outDir := filepath.Dir(gvcfPath)
	vcfDir := filepath.Join(outDir, "VCFs")
	for _, dir := range []string{vcfDir} {
		if _, err := os.Stat(dir); os.IsNotExist(err) {
			if cErr := os.MkdirAll(dir, 0755); cErr != nil {
				return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
			}
			//fmt.Printf("Created directory: %s\n", dir)
		}
	}
	var glnexusPreset string
	if caller == "gatk" {
		glnexusPreset = "gatk"
	} else {
		glnexusPreset = "DeepVariant"
	}
	jointVCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.vcf.gz")
	jointBCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.bcf")
	gvcfFiles := strings.Join(gvcfs, " ")
	glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s %s > %s`, glnexusPreset, filepath.Join(vcfDir, "GLnexus.DB"), gvcfFiles, jointBCF)
	err := runCmd(glnexusCmd, verbose)
	if err != nil {
		return "", fmt.Errorf("glnexus_cli: %w", err)
	}
	bcfCmd := fmt.Sprintf("bcftools view %s | bgzip -@ 4 -c > %s", jointBCF, jointVCF)
	err = runCmd(bcfCmd, verbose)
	if err != nil {
		return "", fmt.Errorf("bcftools view: %w", err)
	}

	indexCmd := fmt.Sprintf(`tabix -p vcf %s`, jointVCF)
	err = runCmd(indexCmd, verbose)
	if err != nil {
		return "", fmt.Errorf("indexing joint VCF with tabix: %w", err)
	}

	return jointVCF, nil
}

// ----------------------------------------- Filter --------------------------------------------------------------- //

func classifyVariant(v *vcfgo.Variant) (isSNP, isIndel bool) {
	refLen := len(v.Ref())
	isSNP = refLen == 1
	for _, alt := range v.Alt() {
		if alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<') {
			continue
		}
		if len(alt) != 1 {
			isSNP = false
		}
		if len(alt) != refLen {
			isIndel = true
		}
	}
	return
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	isSNP, isIndel := classifyVariant(v)

	switch {
	case isSNP:
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if mq, ok := utils.GetFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if mqrs, ok := utils.GetFloat(v, "MQRankSum"); ok && mqrs < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.INDEL_ReadPosRankSum_Min {
			return false
		}
		return true

	default:
		return false
	}
}

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	cleanup := func() { f.Close() }

	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		cleanup = func() { gz.Close(); f.Close() }
		return gz, cleanup, nil
	}
	return f, cleanup, nil
}

func HardFilterVcf(vcf string, cfg utils.HardFilterConfig) (string, error) {

	var hardFilteredVcfPath string
	if strings.HasSuffix(vcf, ".vcf") {
		hardFilteredVcfPath = strings.TrimSuffix(vcf, ".vcf") + ".hard_filtered.vcf.gz"
	} else if strings.HasSuffix(vcf, ".vcf.gz") {
		hardFilteredVcfPath = strings.TrimSuffix(vcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
	} else {
		return "", fmt.Errorf("vcf file %q does not end with .vcf or .vcf.gz", vcf)
	}

	in, cleanup, err := openVCF(vcf)
	if err != nil {
		return "", fmt.Errorf("open %q: %w", vcf, err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(in, false)
	if err != nil {
		return "", fmt.Errorf("VCF header %q: %w", vcf, err)
	}

	outFile, err := os.Create(hardFilteredVcfPath)
	if err != nil {
		return "", fmt.Errorf("create %q: %w", hardFilteredVcfPath, err)
	}
	defer outFile.Close()

	bgzfW := bgzf.NewWriter(outFile, runtime.GOMAXPROCS(0))

	w, err := vcfgo.NewWriter(bgzfW, rdr.Header)
	if err != nil {
		bgzfW.Close()
		return "", fmt.Errorf("VCF writer: %w", err)
	}

	ch := make(chan *vcfgo.Variant, 512)
	readErr := make(chan error, 1)

	go func() {
		defer close(ch)
		for {
			v := rdr.Read()
			if v == nil {
				readErr <- rdr.Error() // nil on clean EOF, error otherwise
				return
			}
			ch <- v
		}
	}()

	for v := range ch {
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			continue
		}
		if PassesHardFilter(v, cfg) {
			w.WriteVariant(v)
		}
	}

	if err := <-readErr; err != nil {
		return "", fmt.Errorf("reading variants: %w", err)
	}
	if err := bgzfW.Close(); err != nil {
		return "", fmt.Errorf("close bgzf: %w", err)
	}
	out, err := exec.Command("tabix", "-p", "vcf", hardFilteredVcfPath).CombinedOutput()
	if err != nil {
		return "", fmt.Errorf("tabix %q: %w\n%s", hardFilteredVcfPath, err, out)
	}

	return hardFilteredVcfPath, nil
}

// --------------------------------------- VC Pipeline ------------------------------------------------------------ //

func getChromsAndContigs(dictFilePath string) ([]SeqInfo, []SeqInfo, error) {
	dictFile, err := os.Open(dictFilePath)
	if err != nil {
		return nil, nil, fmt.Errorf("opening reference dict file %s: %w", dictFilePath, err)
	}
	defer dictFile.Close()

	scanner := bufio.NewScanner(dictFile)
	var seqs []SeqInfo

	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "@SQ") {
			continue
		}
		parts := strings.Split(line, "\t")
		seqID := strings.TrimPrefix(parts[1], "SN:")
		seqLenStr := strings.TrimPrefix(parts[2], "LN:")
		seqLen, aErr := strconv.Atoi(seqLenStr)
		if aErr != nil {
			return nil, nil, fmt.Errorf("parsing LN field for sequence %q in dict file: %w", seqID, aErr)
		}
		seqs = append(seqs, SeqInfo{ID: seqID, Len: seqLen})
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, fmt.Errorf("scanning dict file: %w", err)
	}

	// Sort by length descending so the biggest sequences become "chroms".
	sort.Slice(seqs, func(i, j int) bool { return seqs[i].Len > seqs[j].Len })

	var chroms, contigs []SeqInfo
	if len(seqs) > 21 {
		chroms = append(chroms, seqs[:21]...)
		contigs = append(contigs, seqs[21:]...)
	} else {
		chroms = append(chroms, seqs...)
	}

	// Always promote MT and Pltd into the chroms group.
	for i := 0; i < len(contigs); {
		if contigs[i].ID == "MT" || contigs[i].ID == "Pltd" {
			chroms = append(chroms, contigs[i])
			contigs = append(contigs[:i], contigs[i+1:]...)
		} else {
			i++
		}
	}

	return chroms, contigs, nil
}

func ConcatenateVcfs(vcfs []string, species string, refVer string, outDir string, verbose bool) (string, error) {
	vcfListPath := filepath.Join(outDir, "vcfs.list")
	f, err := os.Create(vcfListPath)
	if err != nil {
		return "", fmt.Errorf("creating %s: %w", vcfListPath, err)
	}
	defer f.Close()
	for _, vcf := range vcfs {
		fmt.Fprintf(f, "%s\n", vcf)
	}

	concatVcfName := fmt.Sprintf("%s.%s.all.vcf.gz", strings.ToUpper(species), strings.ToLower(refVer))
	if len(vcfs) == 1 {
		err = utils.CopyFile(vcfs[0], filepath.Join(outDir, concatVcfName))
		if err != nil {
			return "", fmt.Errorf("copying %s to %s: %w", vcfs[0], filepath.Join(outDir, concatVcfName), err)
		}

	} else {
		cmd := fmt.Sprintf(`gatk MergeVcfs -I %s -O %s`, vcfListPath, filepath.Join(outDir, concatVcfName))
		fmt.Println(cmd)
		err = runCmd(cmd, verbose)

		if err != nil {
			return "", fmt.Errorf("gatk MergeVcfs error: %w\n%s", err)
		}
	}
	return filepath.Join(outDir, concatVcfName), nil
}

func VariantCalling(refFile string, bams []string, out string, species string, refVer string, threadsPerSample int, gatkLogLevel string, caller string, merger string, logFilePath string, dvVer string, modelType string, verbose bool, noMerging bool, cfg utils.HardFilterConfig, noHardFilter bool) (string, error) {

	dictFilePath := refFile[:len(refFile)-len(filepath.Ext(refFile))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		return "", fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return "", fmt.Errorf("getting chroms and contigs from dict file: %w", err)
	}

	finalVCF := ""
	switch strings.ToLower(caller) {
	case "gatk":
		fmt.Printf("Using GATK HaplotypeCaller with log level %s\n", gatkLogLevel)

		jointVcfs := []string{}
		for _, chrom := range chroms {
			// ------------------------------------- Create Gvcfs --------------------------------------------------------- //
			var gvcfs []string
			for _, bam := range bams {
				chromGVCF, err := CreateGvcfGATK(bam, refFile, []SeqInfo{chrom}, species, refVer, gatkLogLevel, verbose, out)
				if err != nil {
					return "", fmt.Errorf("creating gVCF for %s with GATK HaplotypeCaller: %w", chrom.ID, err)
				}
				fmt.Printf("Created gVCF: %s\n", chromGVCF)
				gvcfs = append(gvcfs, chromGVCF)
			}

			if noMerging {
				fmt.Println("Skipping merging of gVCFs")
				continue
			}

			// ------------------------------------- Merge gvcfs ---------------------------------------------------------- //
			var finalChromVcf string
			var mergedChromVcf string
			switch merger {
			case "gatk":
				mergedChromVcf, err = MergeGvcfsGATK(gvcfs, refFile, chrom.ID, species, refVer, verbose, gatkLogLevel, out)
				if err != nil {
					return "", fmt.Errorf("merging gVCFs for chrom %s with GATK: %w", chrom.ID, err)
				}
				fmt.Printf("Merged gVCFs for chromosome %s into joint VCF: %s\n", chrom.ID, mergedChromVcf)
				finalChromVcf = mergedChromVcf
			case "glnexus":
				mergedChromVcf, err = MergeGvcfsGlnexus(gvcfs, chrom.ID, species, refVer, caller, verbose, out)
				finalChromVcf = mergedChromVcf
				if err != nil {
					return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
				}
				fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedChromVcf)
			default:
				return "", fmt.Errorf("unsupported merger: %s", merger)
			}

			// ----------------------------------------- Hard filter vcf ------------------------------------------- //
			if !noHardFilter {
				hardFilteredVcf, err := HardFilterVcf(mergedChromVcf, cfg)
				if err != nil {
					return "", fmt.Errorf("hard filtering VCF: %w", err)
				}
				finalChromVcf = hardFilteredVcf
			}
			jointVcfs = append(jointVcfs, finalChromVcf)

		}
		if len(contigs) > 0 {
			// ------------------------------------- Create Gvcfs --------------------------------------------------- //
			contGvcfs := []string{}
			for _, bam := range bams {
				contigsGVCF, err := CreateGvcfGATK(bam, refFile, contigs, species, refVer, gatkLogLevel, verbose, out)
				if err != nil {
					return "", fmt.Errorf("creating gVCF for contigs with GATK HaplotypeCaller: %w", err)
				}
				fmt.Printf("Created gVCF: %s\n", contigsGVCF)
				contGvcfs = append(contGvcfs, contigsGVCF)
			}
			// ------------------------------------- Merge gvcfs ----------------------------------------------------- //
			var finalContigsVcf string
			var mergedContigsVcf string
			if !noMerging {
				switch merger {
				case "gatk":
					mergedContigsVcf, err = MergeGvcfsGATK(contGvcfs, refFile, "contigs", species, refVer, verbose, gatkLogLevel, out)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GATK: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				case "glnexus":
					mergedContigsVcf, err = MergeGvcfsGlnexus(contGvcfs, "contigs", species, refVer, caller, verbose, out)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				default:
					return "", fmt.Errorf("unsupported merger: %s", merger)
				}
			}

			// ----------------------------------------- Hard filter vcf ------------------------------------------- //
			if !noHardFilter {
				hardFilteredVcf, err := HardFilterVcf(mergedContigsVcf, cfg)
				if err != nil {
					return "", fmt.Errorf("hard filtering VCF: %w", err)
				}
				finalContigsVcf = hardFilteredVcf
			}
			jointVcfs = append(jointVcfs, finalContigsVcf)
		}

		finalVCF, err = ConcatenateVcfs(jointVcfs, species, refVer, out, verbose)
		if err != nil {
			return "", fmt.Errorf("concatenating VCFs: %w", err)
		}

	case "deepvariant":
		fmt.Printf("Using DeepVariant %s with model type %s\n", dvVer, modelType)
		//
		for _, chrom := range chroms {
			// ------------------------------------- Create Gvcfs --------------------------------------------------------- //
			var gvcfs []string
			for _, bam := range bams {
				chromGVCF, err := CreateGvcfDV(bam, refFile, []SeqInfo{chrom}, species, refVer, dvVer, modelType, verbose, out)
				if err != nil {
					return "", fmt.Errorf("creating gVCF for %s with GATK HaplotypeCaller: %w", chrom.ID, err)
				}
				fmt.Printf("Created gVCF: %s\n", chromGVCF)
				gvcfs = append(gvcfs, chromGVCF)
			}

			if noMerging {
				fmt.Println("Skipping merging of gVCFs")
				continue
			}

			// ------------------------------------- Merge gvcfs ---------------------------------------------------------- //

			switch merger {
			case "gatk":
				theVCF, err := MergeGvcfsGATK(gvcfs, refFile, chrom.ID, species, refVer, verbose, gatkLogLevel, out)
				if err != nil {
					return "", fmt.Errorf("merging gVCFs with GATK: %w", err)
				}
				fmt.Printf("Merged gVCFs into joint VCF: %s\n", theVCF)
			case "glnexus":
				theVCF, err := MergeGvcfsGlnexus(gvcfs, chrom.ID, species, refVer, caller, verbose, out)
				if err != nil {
					return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
				}
				fmt.Printf("Merged gVCFs into joint VCF: %s\n", theVCF)
			default:
				return "", fmt.Errorf("unsupported merger: %s", merger)
			}
		}
		if len(contigs) > 0 {
			// ------------------------------------- Create Gvcfs --------------------------------------------------- //
			contGvcfs := []string{}
			for _, bam := range bams {
				contigsGVCF, err := CreateGvcfDV(bam, refFile, contigs, species, refVer, dvVer, modelType, verbose, out)
				if err != nil {
					return "", fmt.Errorf("creating gVCF for contigs with GATK HaplotypeCaller: %w", err)
				}
				fmt.Printf("Created gVCF: %s\n", contigsGVCF)
				contGvcfs = append(contGvcfs, contigsGVCF)
			}
			// ------------------------------------- Merge gvcfs ----------------------------------------------------- //
			if !noMerging {
				switch merger {
				case "gatk":
					theVCF, err := MergeGvcfsGATK(contGvcfs, refFile, "contigs", species, refVer, verbose, gatkLogLevel, out)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GATK: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", theVCF)
				case "glnexus":
					theVCF, err := MergeGvcfsGlnexus(contGvcfs, "contigs", species, refVer, caller, verbose, out)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", theVCF)
				default:
					return "", fmt.Errorf("unsupported merger: %s", merger)
				}
			}
		}

	default:
		return "", fmt.Errorf("unsupported caller: %s", caller)

	}
	return finalVCF, nil
}

	// ------------------------------------- Merge gvcfs ---------------------------------------------------------- //

	// ----------------------------------------- Hard Filter vcfs ------------------------------------------------- //

	// ----------------------------------------- open / create log file --------------------------------------------- //
	fmt.Println("Opening log file ...", logFilePath)
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		return "", fmt.Errorf("opening log file %s: %w", logFilePath, err)
	}
	defer logFile.Close()

	jlog := slog.New(slog.NewJSONHandler(logFile, nil))

	finalVcf := filepath.Join(out, species+".joint_hard_filtered.vcf.gz")

	jlog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

	// ── resume check ─────────────────────────────────────────────────────────
	logged := utils.ParseLogFile(logFilePath)
	if utils.StageHasCompleted(logged, stageMergeVcfs, "ALL", "ALL") {
		fmt.Println("All stages completed. Skipping.")
		return finalVcf, nil
	}

	// ── scan FASTA sequences ─────────────────────────────────────────────────
	r := fasta.NewReader(reader, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)

	totalCores := runtime.NumCPU()
	fmt.Printf("Available CPU cores: %d\n", totalCores)

	maxParallelJobs := totalCores / threadsPerSample
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
	}

	//var (
	//	wg          sync.WaitGroup
	//	mu          sync.Mutex // protects jointvSlice
	//	sem         = make(chan struct{}, maxParallelJobs)
	//	errCh       = make(chan error, 64) // FIX: collect goroutine errors instead of log.Fatalf
	//	jointvSlice []string
	//)

	var (
		wg          sync.WaitGroup
		mu          sync.Mutex // protects jointvSlice
		sem         = make(chan struct{}, maxParallelJobs)
		callerSem   = make(chan struct{}, maxParallelJobs) // FIX: Limits GATK global concurrency
		errCh       = make(chan error, 64)                 // FIX: collect goroutine errors instead of log.Fatalf
		jointvSlice []string
	)

	// ── per-chromosome goroutines ─────────────────────────────────────────────
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)

		chromDir := strings.ReplaceAll(seq.ID, ".", "_")
		chromDirPath := filepath.Join(out, chromDir)
		gvcfPath := filepath.Join(chromDirPath, "gvcfs")
		tmpPath := filepath.Join(chromDirPath, "tmp")
		tmp2Path := filepath.Join(chromDirPath, "tmp2")
		vcfPath := filepath.Join(chromDirPath, "VCFs")

		jointVCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.vcf.gz")
		jointBCF := filepath.Join(vcfPath, species+"_"+chromDir+".joint.bcf")
		snpVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".SNP.vcf.gz"
		indelVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".INDEL.vcf.gz"
		hardFilteredSNPs := strings.TrimSuffix(snpVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		hardFilteredINDELs := strings.TrimSuffix(indelVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		hardFilteredVCF := strings.TrimSuffix(jointVCF, ".vcf.gz") + ".hard_filtered.vcf.gz"
		theDB := filepath.Join(chromDirPath, chromDir+"DB")

		// Create required directories before spawning the goroutine.
		for _, dir := range []string{chromDirPath, gvcfPath, tmpPath, tmp2Path, vcfPath} {
			if _, err := os.Stat(dir); os.IsNotExist(err) {
				if cErr := os.MkdirAll(dir, 0755); cErr != nil {
					slog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "STATUS", fmt.Sprintf("FAILED: %v", cErr))
					jlog.Error("VARIANT CALLING", "PROGRAM", "mkdir", "STATUS", fmt.Sprintf("FAILED: %v", cErr))
					return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
				}
				//fmt.Printf("Created directory: %s\n", dir)
			}
		}

		sem <- struct{}{}
		wg.Add(1)

		switch strings.ToLower(caller) {

		// ══════════════════════════════════════════════════════════════════════
		// GATK branch
		// ══════════════════════════════════════════════════════════════════════
		case "gatk":
			go func(seq *linear.Seq, gvcfPath, tmpPath, vcfPath, theDB string, jointVCF, jointBCF, snpVCF, indelVCF string, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF string) {
				defer func() { wg.Done(); <-sem }()

				slog.Info(fmt.Sprintf("Working on chromosome %s …\n", seq.ID))

				var (
					callerWg  sync.WaitGroup
					vSlice    []string
					callerErr error
					callerMu  sync.Mutex
				)

				for _, bam := range bams {
					bamName := filepath.Base(bam)

					theGVCF := gvcfForBam(bamName, gvcfPath, chromDir)
					if theGVCF == "" {
						callerMu.Lock()
						callerErr = errors.Join(callerErr,
							fmt.Errorf("BAM %s has neither .bam nor .cram extension", bamName))
						callerMu.Unlock()
						continue
					}

					vSlice = append(vSlice, "-V "+theGVCF)

					callerSem <- struct{}{}
					callerWg.Add(1)
					go func(bam, bamName, theGVCF string) {
						defer func() { callerWg.Done(); <-callerSem }()

						if utils.StageHasCompleted(logged, stageHaplotypeCaller, bamName, seq.ID) {
							slog.Info(fmt.Sprintf("%s already completed for %s / %s. Skipping.", stageHaplotypeCaller, bamName, seq.ID))
							return
						}

						hapCmdStr := fmt.Sprintf(
							`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`,
							refFile, bam, seq.ID, theGVCF, gatkLogLevel,
						)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", hapCmdStr)
						slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						if hapErr := runCmd(hapCmdStr, verbose); hapErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", hapErr))
							// FIX: propagate error instead of log.Fatalf
							callerMu.Lock()
							callerErr = errors.Join(callerErr, hapErr)
							callerMu.Unlock()
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					}(bam, bamName, theGVCF)
				}
				callerWg.Wait()

				if callerErr != nil {
					errCh <- fmt.Errorf("chromosome %s HaplotypeCaller: %w", seq.ID, callerErr)
					return
				}
				fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

				if noMerging {
					return
				}

				// ── joint calling ─────────────────────────────────────────────
				switch merger {
				case "gatk":
					if !utils.StageHasCompleted(logged, stageGenomicsDBImport, "ALL", seq.ID) {
						vArgs := strings.Join(vSlice, " ")
						if rErr := os.RemoveAll(theDB); rErr != nil {
							errCh <- fmt.Errorf("removing %s: %w", theDB, rErr)
							return
						}
						gDBCmd := fmt.Sprintf(
							`gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport %s `+
								`--genomicsdb-workspace-path %s --tmp-dir %s -L %s `+
								`--genomicsdb-shared-posixfs-optimizations true --batch-size 50 `+
								`--bypass-feature-reader --verbosity %s`,
							vArgs, theDB, tmpPath, seq.ID, gatkLogLevel,
						)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						if gErr := runCmd(gDBCmd, verbose); gErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gErr))
							errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGenomicsDBImport, gErr)
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenomicsDBImport, seq.ID))
					}

					if !utils.StageHasCompleted(logged, stageGenotypeGVCFs, "ALL", seq.ID) {
						genoCmd := fmt.Sprintf(
							`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
							refFile, theDB, jointVCF, tmpPath, gatkLogLevel,
						)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						if gtErr := runCmd(genoCmd, verbose); gtErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", gtErr))
							errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGenotypeGVCFs, gtErr)
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenotypeGVCFs, seq.ID))
					}

				case "glnexus":
					if !utils.StageHasCompleted(logged, stageGLNEXUS, "ALL", seq.ID) {
						glnexusCmd := fmt.Sprintf(`glnexus_cli --config gatk --dir %s %s/*.g.vcf.gz > %s`, filepath.Join(tmpPath, "GLnexus.DB"), gvcfPath, jointBCF)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", glnexusCmd)
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						if glErr := runCmd(glnexusCmd, verbose); glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGLNEXUS, glErr)
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUS, seq.ID))
					}

					if !utils.StageHasCompleted(logged, stageGLNEXUSBCFTOOLS, "ALL", seq.ID) {
						bcfCmd := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", bcfCmd)
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						if glErr := runCmd(bcfCmd, verbose); glErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
							errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGLNEXUSBCFTOOLS, glErr)
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSBCFTOOLS, seq.ID))
					}

					if !utils.StageHasCompleted(logged, stageGLNEXUSIndex, "ALL", seq.ID) {
						indexCmd := fmt.Sprintf(`gatk IndexFeatureFile -I %s`, jointVCF)
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", indexCmd)
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						if indErr := runCmd(indexCmd, verbose); indErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", indErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", indErr))
							// Indexing failure is non-fatal; log and continue.
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSIndex, seq.ID))
					}

				default:
					errCh <- fmt.Errorf("chromosome %s: unknown merger %q (must be \"gatk\" or \"glnexus\")", seq.ID, merger)
					return
				}

				// ------------------------------------- hard filtering --------------------------------------------- //
				_, _, hardFilteredSNPs, hardFilteredINDELs, err := runHardFiltering(seq.ID, jointVCF, snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, gatkLogLevel, verbose, logged, jlog)
				if err != nil {
					errCh <- fmt.Errorf("chromosome %s hard filtering: %w", seq.ID, err)
					return
				}

				// --------------------------------------- Per Chrom Merging ---------------------------------------- //
				if err := mergeHardFiltered(seq.ID, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF, gatkLogLevel, verbose, logged, jlog); err != nil {
					errCh <- fmt.Errorf("chromosome %s MergeVcfs: %w", seq.ID, err)
					return
				}

				mu.Lock()
				jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)
				mu.Unlock()

			}(seq, gvcfPath, tmpPath, vcfPath, theDB, jointVCF, jointBCF, snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF)

		// ══════════════════════════════════════════════════════════════════════
		// DeepVariant branch
		// ══════════════════════════════════════════════════════════════════════
		case "deepvariant":
			fmt.Printf("Variant Calling using DeepVariant on chromosome %s …\n\n", seq.ID)

			go func(seq *linear.Seq, gvcfPath, vcfPath string, jointVCF, jointBCF, snpVCF, indelVCF string, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF string) {
				defer func() { wg.Done(); <-sem }()
				var (
					callerWg  sync.WaitGroup
					callerErr error
					callerMu  sync.Mutex
				)

				for _, bam := range bams {
					bamName := filepath.Base(bam)

					gvcfName, vcfName := gvcfAndVcfNamesForBam(bamName, chromDir)
					if gvcfName == "" {
						callerMu.Lock()
						callerErr = errors.Join(callerErr,
							fmt.Errorf("BAM %s has neither .bam nor .cram extension", bamName))
						callerMu.Unlock()
						continue
					}

					callerSem <- struct{}{}
					callerWg.Add(1)
					go func(bam, bamName, vcfName, gvcfName string) {
						defer func() { callerWg.Done(); <-callerSem }()

						if utils.StageHasCompleted(logged, stageDeepVariant, bamName, seq.ID) {
							slog.Info(fmt.Sprintf("%s already completed for %s / %s. Skipping.", stageDeepVariant, bamName, seq.ID))
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "STARTED")

						if dvErr := runDeepVariantContainer(bam, refFile, gvcfPath, gvcfName, vcfName, seq.ID, dvVer, modelType, verbose); dvErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", dvErr))
							callerMu.Lock()
							callerErr = errors.Join(callerErr, dvErr)
							callerMu.Unlock()
							return
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					}(bam, bamName, vcfName, gvcfName)
				}
				callerWg.Wait()

				if callerErr != nil {
					errCh <- fmt.Errorf("chromosome %s DeepVariant: %w", seq.ID, callerErr)
					return
				}

				if noMerging {
					return
				}

				if merger != "glnexus" {
					errCh <- fmt.Errorf("chromosome %s: only \"glnexus\" merger is supported with DeepVariant (got %q)", seq.ID, merger)
					return
				}

				// ── GLnexus ───────────────────────────────────────────────────
				if !utils.StageHasCompleted(logged, stageGLNEXUS, "ALL", seq.ID) {
					// Use DeepVariantWGS preset (or DeepVariantWES for exome).
					// Do NOT pass the caller string — GLnexus preset names are
					// case-sensitive and differ from our internal caller names.
					glnexusPreset := "DeepVariant"
					if strings.EqualFold(modelType, "WES") {
						glnexusPreset = "DeepVariantWES"
					}

					glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s %s/*.g.vcf.gz > %s`, glnexusPreset, filepath.Join(vcfPath, "GLnexus.DB"), gvcfPath, jointBCF)

					jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", glnexusCmd)
					slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					if glErr := runCmd(glnexusCmd, verbose); glErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
						slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
						errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGLNEXUS, glErr)
						return
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				} else {
					slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUS, seq.ID))
				}

				if !utils.StageHasCompleted(logged, stageGLNEXUSBCFTOOLS, "ALL", seq.ID) {
					bcfCmd := fmt.Sprintf(`bcftools view %s | bgzip -@ 4 -c > %s`, jointBCF, jointVCF)
					jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED", "CMD", bcfCmd)
					slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "STARTED")
					if glErr := runCmd(bcfCmd, verbose); glErr != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
						slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", fmt.Sprintf("FAILED: %v", glErr))
						errCh <- fmt.Errorf("chromosome %s %s: %w", seq.ID, stageGLNEXUSBCFTOOLS, glErr)
						return
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", seq.ID, "STATUS", "COMPLETED")
				} else {
					slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSBCFTOOLS, seq.ID))
				}

				// ── hard filtering ────────────────────────────────────────────
				_, _, hardFilteredSNPs, hardFilteredINDELs, err :=
					runHardFiltering(seq.ID, jointVCF, snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, gatkLogLevel, verbose, logged, jlog)
				if err != nil {
					errCh <- fmt.Errorf("chromosome %s hard filtering: %w", seq.ID, err)
					return
				}

				// ── per-chrom merge ───────────────────────────────────────────
				if err := mergeHardFiltered(seq.ID, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF, gatkLogLevel, verbose, logged, jlog); err != nil {
					errCh <- fmt.Errorf("chromosome %s MergeVcfs: %w", seq.ID, err)
					return
				}

				mu.Lock()
				jointvSlice = append(jointvSlice, "-I "+hardFilteredVCF)
				mu.Unlock()

			}(seq, gvcfPath, vcfPath, jointVCF, jointBCF, snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF)

		default:
			wg.Done()
			<-sem
			return "", fmt.Errorf("unknown caller %q: must be \"gatk\" or \"deepvariant\"", caller)
		}
	}

	// Wait for all chromosome goroutines, then close the error channel.
	go func() {
		wg.Wait()
		close(errCh)
	}()

	// Collect any errors reported by goroutines (This will now run immediately and continuously drain errors)
	var combinedErr error
	for e := range errCh {
		combinedErr = errors.Join(combinedErr, e)
		// Optional: you can print them as they happen so you don't have to wait for the end
		// jlog.Error("Pipeline Error", "details", e.Error())
	}

	if combinedErr != nil {
		return "", combinedErr
	}

	// ── final multi-sample merge ──────────────────────────────────────────────
	if !noMerging {
		fmt.Println("Merging ALL hard-filtered VCFs …")
		mergeCmdStr2 := fmt.Sprintf(`gatk MergeVcfs %s -O %s --VERBOSITY %s`,
			strings.Join(jointvSlice, " "), finalVcf, gatkLogLevel)

		jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED", "CMD", mergeCmdStr2)
		slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

		if mErr := runCmd(mergeCmdStr2, verbose); mErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
			slog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", mErr))
			return "", fmt.Errorf("final MergeVcfs: %w", mErr)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		fmt.Printf("%s\n\n", strings.Repeat("-", 90))
	}

	return finalVcf, nil
}

// FailedTask tracks which sample/chromosome combinations failed.
type FailedTask struct {
	Sample string
	Chrom  string
	Reason error
}

// SampleWork bundles per-sample data needed by the worker goroutine.
type SampleWork struct {
	Sample  string
	Cram    string
	CramDir string // parent of the "bams" directory; used to derive gvcf output path
}

// ─────────────────────────────────────────────────────────────────────────────
// runCmd is a small helper that respects the verbose flag without repeating
// the if/else everywhere.
// ─────────────────────────────────────────────────────────────────────────────

func runCmd(cmdStr string, verbose bool) error {
	if verbose {
		return utils.RunBashCmdVerbose(cmdStr)
	}
	return utils.RunBashCmd(cmdStr)
}

// ─────────────────────────────────────────────────────────────────────────────
// HapCaller – thin wrapper kept for backwards-compat with other callers
// ─────────────────────────────────────────────────────────────────────────────

func HapCaller(ref string, bam string, verbose bool, vcfType string, vcf string) error {
	cmdStrHap := fmt.Sprintf(`gatk HaplotypeCaller -R %s -I %s -O %s`, ref, bam, vcf)
	fmt.Println(cmdStrHap)
	slog.Info("HaplotypeCaller", "cmd", cmdStrHap)
	return runCmd(cmdStrHap, verbose)
}

// ─────────────────────────────────────────────────────────────────────────────
// getChromsAndContigs reads a .dict file and splits sequences into the 21
// longest chromosomes (plus MT/Pltd) and everything else (contigs).
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// FindBamOrVcfs locates BAM/CRAM or gVCF files for a sample using glob.
// ─────────────────────────────────────────────────────────────────────────────

func FindBamOrVcfs(dataDirAbs string, species string, sample string, refVer string, bamOrVcf string, chrom string) ([]string, error) {
	var pattern string
	switch bamOrVcf {
	case "bams":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*bqsr.cram",
			dataDirAbs, strings.ToLower(species), sample, refVer)
	case "gvcfs":
		pattern = fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/gvcfs/*%s.g.vcf.gz",
			dataDirAbs, strings.ToLower(species), sample, refVer, chrom)
	}

	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, fmt.Errorf("glob error: %w", err)
	}
	if len(matches) == 0 {
		return nil, fmt.Errorf("no %s files found in %s", bamOrVcf, sample)
	}
	if len(matches) > 1 {
		return matches, fmt.Errorf("multiple %s files found in %s", bamOrVcf, sample)
	}
	return matches, nil
}

func runDeepVariantContainer(bam, refFile, outputDir, gvcfName, vcfName, regions, dvVer, modelType string, verbose bool) error {
	bamDir := filepath.Dir(bam)
	bamName := filepath.Base(bam)
	refDir := filepath.Dir(refFile)
	refName := filepath.Base(refFile)

	// Build the --regions argument.
	// If `regions` is a file path it must live inside outputDir (mounted as
	// /output) so the container can reach it.  We pass the container-side path.
	var containerRegions string
	if filepath.IsAbs(regions) || strings.ContainsRune(regions, '/') {
		// It's a file path — translate host path to container path.
		rel, err := filepath.Rel(outputDir, regions)
		if err != nil {
			return fmt.Errorf("interval list %s is not inside outputDir %s: %w", regions, outputDir, err)
		}
		containerRegions = "/output/" + rel
	} else {
		// Plain chromosome/contig ID — pass as-is.
		containerRegions = regions
	}

	// Unique intermediate dir per (bam, region) prevents FileExistsError races
	// when multiple containers write to the same /output mount.
	safeRegion := strings.NewReplacer(
		string(os.PathSeparator), "_",
		".", "_",
		"/", "_",
	).Replace(regions)
	intermediateName := fmt.Sprintf("tmp_%s_%s",
		strings.TrimSuffix(bamName, filepath.Ext(bamName)),
		safeRegion,
	)

	// Note: image tag is NOT quoted — docker does not accept quoted tags.
	dvCmdStr := fmt.Sprintf(
		`docker run -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
			`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
			`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
			`--intermediate_results_dir /output/%s`,
		bamDir, refDir, outputDir, dvVer,
		modelType, refName, bamName,
		containerRegions, vcfName, gvcfName,
		intermediateName,
	)
	fmt.Printf("\n%s\n\n", dvCmdStr)
	return runCmd(dvCmdStr, verbose)
}

// ─────────────────────────────────────────────────────────────────────────────
// CreateGvcf runs either GATK HaplotypeCaller or DeepVariant (via docker) to
// produce a single per-sample/per-region gVCF.
// ─────────────────────────────────────────────────────────────────────────────

func CreateGvcf(bam string, refFile string, chroms []SeqInfo, theGVCF string, gatkLogLevel string, caller string, dvVer string, modelType string, verbose bool) (string, error) {
	var regionArg string

	switch strings.ToLower(caller) {
	case "gatk":
		if len(chroms) == 1 {
			regionArg = chroms[0].ID
		} else {
			// Use a unique temp file per call to avoid races when multiple
			// goroutines run CreateGvcf concurrently.
			f, err := os.CreateTemp("", "gatk_intervals_contigs_*.list")
			if err != nil {
				return "", fmt.Errorf("creating GATK interval list: %w", err)
			}
			defer os.Remove(f.Name())
			defer f.Close()
			for _, c := range chroms {
				fmt.Fprintln(f, c.ID)
			}
			regionArg = f.Name()
		}

		hapCmdStr := fmt.Sprintf(
			`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s`,
			refFile, bam, regionArg, theGVCF, gatkLogLevel,
		)
		fmt.Printf("\n%s\n\n", hapCmdStr)
		return theGVCF, runCmd(hapCmdStr, verbose)

	case "deepvariant":
		gvcfDir := filepath.Dir(theGVCF)
		gvcfName := filepath.Base(theGVCF)
		vcfName := strings.Replace(gvcfName, ".g.vcf.gz", ".vcf.gz", 1)

		var regions string
		if len(chroms) == 1 {
			regions = chroms[0].ID
		} else {
			// Write a BED file into gvcfDir (which is mounted as /output) so
			// the container can read it via its /output mount.
			f, err := os.CreateTemp(gvcfDir, "deepvariant_intervals_*.bed")
			if err != nil {
				return "", fmt.Errorf("creating DeepVariant interval list: %w", err)
			}
			defer os.Remove(f.Name())
			defer f.Close()
			for _, c := range chroms {
				fmt.Fprintf(f, "%s\t0\t%d\n", c.ID, c.Len)
			}
			regions = f.Name() // absolute host path inside gvcfDir → translatable to /output/...
		}

		return theGVCF, runDeepVariantContainer(bam, refFile, gvcfDir, gvcfName, vcfName, regions, dvVer, modelType, verbose)

	default:
		return "", fmt.Errorf("unknown caller %q: must be \"gatk\" or \"deepvariant\"", caller)
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// runHardFiltering selects SNPs and INDELs and hard-filters each.
// Returns updated paths for the caller to record.
// ─────────────────────────────────────────────────────────────────────────────

func runHardFiltering(chromID string, jointVCF, snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs string, gatkLogLevel string, verbose bool, logged []utils.LogEntry, jlog *slog.Logger) (string, string, string, string, error) {

	fmt.Println("Hard filtering joint VCF …")

	if !utils.StageHasCompleted(logged, stageSelectSNPs, "ALL", chromID) {
		jlog.Info("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		newSnpVCF, sErr := GetVariantType(jointVCF, "SNP", gatkLogLevel, verbose)
		if sErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
			slog.Error("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", sErr))
			return snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, sErr
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageSelectSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		snpVCF = newSnpVCF
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageSelectSNPs, chromID))
	}
	fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

	if !utils.StageHasCompleted(logged, stageSelectINDELs, "ALL", chromID) {
		jlog.Info("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		newIndelVCF, iErr := GetVariantType(jointVCF, "INDEL", gatkLogLevel, verbose)
		if iErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
			slog.Error("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", iErr))
			return snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, iErr
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageSelectINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		indelVCF = newIndelVCF
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageSelectINDELs, chromID))
	}
	fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

	if !utils.StageHasCompleted(logged, stageHardFilterSNPs, "ALL", chromID) {
		jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		newHFS, hsErr := HardFilterSNPs(snpVCF, gatkLogLevel, verbose)
		if hsErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
			slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", hsErr))
			return snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, hsErr
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterSNPs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		hardFilteredSNPs = newHFS
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageHardFilterSNPs, chromID))
	}
	fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

	if !utils.StageHasCompleted(logged, stageHardFilterINDELs, "ALL", chromID) {
		jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
		newHFI, hiErr := HardFilterINDELs(indelVCF, gatkLogLevel, verbose)
		if hiErr != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
			slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", hiErr))
			return snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, hiErr
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilterINDELs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
		hardFilteredINDELs = newHFI
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageHardFilterINDELs, chromID))
	}
	fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))

	return snpVCF, indelVCF, hardFilteredSNPs, hardFilteredINDELs, nil
}

// mergeHardFiltered merges per-type hard-filtered VCFs into a single per-chrom VCF.
func mergeHardFiltered(chromID, hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF string, gatkLogLevel string, verbose bool, logged []utils.LogEntry, jlog *slog.Logger) error {
	if utils.StageHasCompleted(logged, stageMergeVcfs, "ALL", chromID) {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageMergeVcfs, chromID))
		return nil
	}
	jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")
	slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "STARTED")

	mergeCmd := fmt.Sprintf(`gatk MergeVcfs -I %s -I %s -O %s --VERBOSITY %s`,
		hardFilteredSNPs, hardFilteredINDELs, hardFilteredVCF, gatkLogLevel)
	if mErr := runCmd(mergeCmd, verbose); mErr != nil {
		jlog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		slog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", fmt.Sprintf("FAILED: %v", mErr))
		return mErr
	}
	jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
	slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", chromID, "STATUS", "COMPLETED")
	fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	return nil
}

// ─────────────────────────────────────────────────────────────────────────────
// gvcfForBam returns the gVCF output path for a given BAM name and chromosome
// directory.  Returns "" if the extension is unrecognised.
// ─────────────────────────────────────────────────────────────────────────────

func gvcfForBam(bamName, gvcfPath, chromDir string) string {
	switch {
	case strings.HasSuffix(bamName, ".bam"):
		return filepath.Join(gvcfPath, strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
	case strings.HasSuffix(bamName, ".cram"):
		return filepath.Join(gvcfPath, strings.Replace(bamName, "cram", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1))
	default:
		return ""
	}
}

// gvcfAndVcfNamesForBam returns (gvcfName, vcfName) base-names for DeepVariant output.
// Returns ("", "") if the extension is unrecognised.
func gvcfAndVcfNamesForBam(bamName, chromDir string) (string, string) {
	switch {
	case strings.HasSuffix(bamName, ".bam"):
		return strings.Replace(bamName, "bam", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1),
			strings.Replace(bamName, "bam", fmt.Sprintf("%s.vcf.gz", chromDir), 1)
	case strings.HasSuffix(bamName, ".cram"):
		return strings.Replace(bamName, "cram", fmt.Sprintf("%s.g.vcf.gz", chromDir), 1),
			strings.Replace(bamName, "cram", fmt.Sprintf("%s.vcf.gz", chromDir), 1)
	default:
		return "", ""
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// VariantCallingConfig reads a config file and delegates to VariantCalling.
// ─────────────────────────────────────────────────────────────────────────────

func VariantCallingConfig(configFile string, species string, maxParallelJobs int, gatkLogLevel string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool) {
	fmt.Println("Reading config file ...")
	cfg, err := utils.ReadConfig(configFile)
	if err != nil {
		fmt.Printf("Error reading config: %v\n", err)
		return
	}
	fmt.Println("Reference:", cfg.Reference)
	fmt.Println("Bams:", cfg.Bams)
	fmt.Println("Output:", cfg.OutputDir)

	refFile := cfg.Reference
	if _, rErr := os.Stat(refFile); rErr != nil {
		fmt.Printf("Reference file %s does not exist\n", refFile)
		return
	}

	bams := cfg.Bams
	if len(bams) == 0 {
		fmt.Println("You must provide at least one BAM file")
		return
	}
	for _, b := range bams {
		if _, err := os.Stat(b); err != nil {
			fmt.Printf("BAM file %s is not a valid file path: %v\n", b, err)
			return
		}
	}

	outputDir := cfg.OutputDir
	outInfo, outErr := os.Stat(outputDir)
	if outErr != nil {
		if os.IsNotExist(outErr) {
			fmt.Printf("Output directory %s does not exist. Attempting to create it.\n", outputDir)
			if createErr := os.MkdirAll(outputDir, 0755); createErr != nil {
				fmt.Printf("Failed to create output directory %s: %v\n", outputDir, createErr)
				return
			}
			fmt.Printf("Output directory %s created successfully.\n", outputDir)
		} else {
			fmt.Printf("Error accessing output directory %s: %v\n", outputDir, outErr)
			return
		}
	} else if !outInfo.IsDir() {
		fmt.Printf("Output path %s is not a directory\n", outputDir)
		return
	}

	logFilePath := filepath.Join(outputDir, "variant_calling.log")
	finalVcf, fErr := VariantCalling(refFile, bams, outputDir, species, maxParallelJobs, gatkLogLevel, caller, merger, logFilePath, dvVer, modelType, verbose, noMerging)
	if fErr != nil {
		fmt.Printf("Error calling variants: %v\nNo multisample VCF: %s\n", fErr, finalVcf)
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// VariantCallingDir performs variant calling on all discovered samples in a
// directory tree.
// ─────────────────────────────────────────────────────────────────────────────

func VariantCallingDir(dataDir string, species string, refVer string, genomesDir string, refFasta string, caller string, merger string, dvVer string, modelType string, verbose bool, noMerging bool, gatkLogLevel string, quick bool, threadsPerJob int) error {

	// ── validate paths ────────────────────────────────────────────────────────
	fmt.Println("Checking paths ...")

	dInfo, err := os.Stat(dataDir)
	if err != nil {
		return fmt.Errorf("accessing data directory %s: %w", dataDir, err)
	}
	if !dInfo.IsDir() {
		return fmt.Errorf("data directory %s is not a directory", dataDir)
	}
	dataDirAbs, err := filepath.Abs(dataDir)
	if err != nil {
		return fmt.Errorf("getting absolute path for data directory %s: %w", dataDir, err)
	}

	if species == "" {
		return fmt.Errorf("species name must not be empty")
	}
	if refVer == "" {
		return fmt.Errorf("reference version must not be empty")
	}

	resolvedFasta, resolveErr := utils.ResolveRefFasta(refFasta, genomesDir, species, refVer)
	if resolveErr != nil {
		return fmt.Errorf("could not resolve reference fasta: %w", resolveErr)
	}
	fastaInfo, err := os.Stat(resolvedFasta)
	if err != nil {
		return fmt.Errorf("accessing reference fasta file %s: %w", resolvedFasta, err)
	}
	if !fastaInfo.Mode().IsRegular() {
		return fmt.Errorf("reference fasta %s is not a regular file", resolvedFasta)
	}

	dictFilePath := resolvedFasta[:len(resolvedFasta)-len(filepath.Ext(resolvedFasta))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		return fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return fmt.Errorf("getting chromosomes and IDs: %w", err)
	}

	color.Cyan("All file paths valid\n…\n\n")

	// ── discover samples ──────────────────────────────────────────────────────
	color.Cyan("================================== Discovering samples =================================\n\n")
	pattern := filepath.Join(dataDir, species, "*", "*", "reference_genomes")
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return fmt.Errorf("finding samples in %s: %w", pattern, err)
	}

	color.Green("SAMPLES FOUND:\n%s\n\n", strings.Repeat("-", 69))
	seen := make(map[string]struct{}, len(matches))
	var samples []string
	for _, match := range matches {
		s := filepath.Base(filepath.Dir(match))
		if _, ok := seen[s]; !ok {
			seen[s] = struct{}{}
			samples = append(samples, s)
			fmt.Println(s)
		}
	}
	color.Cyan("\nFound %d sample(s) for %s\n%s\n\n", len(samples), species, strings.Repeat("=", 34))

	// ── resolve valid samples ─────────────────────────────────────────────────
	fmt.Println("Looking for BAM files …")

	var (
		mu           sync.Mutex
		missingBams  []string
		multipleBams []string
		validSamples []SampleWork
		wg1          sync.WaitGroup
	)

	for _, sample := range samples {
		wg1.Add(1)
		go func(sample string) {
			defer wg1.Done()

			pat := fmt.Sprintf("%s/%s/*/%s/reference_genomes/%s/bams/*bqsr.cram",
				dataDirAbs, strings.ToLower(species), sample, refVer)
			cramFiles, cErr := filepath.Glob(pat)
			if cErr != nil {
				color.Red("Error finding BAM files: %v\n", cErr)
				mu.Lock()
				missingBams = append(missingBams, sample)
				mu.Unlock()
				return
			}

			switch len(cramFiles) {
			case 0:
				color.Red("[%s] bqsr.cram MISSING ❌\n", sample)
				mu.Lock()
				missingBams = append(missingBams, sample)
				mu.Unlock()
			case 1:
				color.Green("[%s] bqsr.cram FOUND ✅: %s\n", sample, cramFiles[0])
				color.Cyan("Validating BAM file …\n")
				if valErr := utils.ValidateBam(cramFiles[0], resolvedFasta, verbose, quick); valErr != nil {
					color.Red("[%s] bqsr.cram not valid: %v\n", sample, valErr)
					mu.Lock()
					missingBams = append(missingBams, sample)
					mu.Unlock()
				} else {
					color.Green("[%s] bqsr.cram valid ✅\n", sample)
					mu.Lock()
					validSamples = append(validSamples, SampleWork{
						Sample:  sample,
						Cram:    cramFiles[0],
						CramDir: filepath.Dir(filepath.Dir(filepath.Dir(cramFiles[0]))),
					})
					mu.Unlock()
				}
			default:
				color.Red("[%s] Multiple bqsr.cram files found — skipping ❌\n", sample)
				mu.Lock()
				multipleBams = append(multipleBams, sample)
				mu.Unlock()
			}
		}(sample)
	}
	wg1.Wait()

	color.Green("\nValid:    %d\n", len(validSamples))
	color.Red("Missing:  %d\n", len(missingBams))
	color.Yellow("Multiple: %d\n", len(multipleBams))
	fmt.Printf("\n%s\n\n", strings.Repeat("=", 89))

	if len(validSamples) == 0 {
		return fmt.Errorf("no valid samples found (missing: %d, multiple: %d)", len(missingBams), len(multipleBams))
	}

	color.Cyan("========= Creating gVCFs for %d valid samples =========\n\n", len(validSamples))

	// ── resource-aware concurrent gVCF creation ───────────────────────────────
	var (
		failedTasks []FailedTask
		failedMu    sync.Mutex
	)

	totalCores := runtime.NumCPU()
	if threadsPerJob <= 0 {
		threadsPerJob = 4
	}
	maxParallelJobs := totalCores / threadsPerJob
	if maxParallelJobs < 1 {
		maxParallelJobs = 1
	}
	color.Cyan("Machine has %d cores. %d threads per job → max %d parallel jobs\n\n",
		totalCores, threadsPerJob, maxParallelJobs)

	type workItem struct {
		sw       SampleWork
		chroms   []SeqInfo
		label    string
		isContig bool
	}

	totalWork := len(validSamples) * (len(chroms) + 1)
	workCh := make(chan workItem, totalWork)
	var wg2 sync.WaitGroup

	for i := range maxParallelJobs {
		wg2.Add(1)
		go func(workerID int) {
			defer wg2.Done()
			for item := range workCh {
				sw := item.sw
				label := item.label

				cramName := filepath.Base(sw.Cram)
				var gvcfName string
				if item.isContig {
					gvcfName = strings.Replace(cramName, ".cram", ".contigs.g.vcf.gz", 1)
					gvcfName = strings.Replace(gvcfName, ".bam", ".contigs.g.vcf.gz", 1)
				} else {
					gvcfName = strings.Replace(cramName, ".cram", fmt.Sprintf(".%s.g.vcf.gz", label), 1)
					gvcfName = strings.Replace(gvcfName, ".bam", fmt.Sprintf(".%s.g.vcf.gz", label), 1)
				}
				gvcfPath := filepath.Join(sw.CramDir, refVer, "gvcfs", gvcfName)

				var vcfFiles []string
				if item.isContig {
					if _, statErr := os.Stat(gvcfPath); statErr == nil {
						vcfFiles = []string{gvcfPath}
					}
				} else {
					vcfFiles, _ = FindBamOrVcfs(dataDirAbs, species, sw.Sample, refVer, "gvcfs", label)
				}

				switch len(vcfFiles) {
				case 0:
					color.Cyan("[Worker %d] [%s] gVCF not found for %s — creating …\n\n", workerID, sw.Sample, label)
					if _, err := CreateGvcf(sw.Cram, resolvedFasta, item.chroms, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
						color.Red("[Worker %d] [%s] Error creating gVCF for %s: %v\n\n", workerID, sw.Sample, label, err)
						failedMu.Lock()
						failedTasks = append(failedTasks, FailedTask{Sample: sw.Sample, Chrom: label, Reason: err})
						failedMu.Unlock()
					} else {
						color.Green("[Worker %d] [%s] gVCF for %s created successfully\n\n", workerID, sw.Sample, label)
					}

				case 1:
					vcf := vcfFiles[0]
					color.Green("[Worker %d] [%s] gVCF for %s exists: %s\n\n", workerID, sw.Sample, label, vcf)
					fmt.Printf("[Worker %d] [%s] checking integrity of %s …\n", workerID, sw.Sample, color.BlueString(vcf))
					if vErr := utils.ValidateGvcf(vcf, verbose, quick); vErr != nil {
						color.Red("[Worker %d] [%s] gVCF %s corrupted (%v) — re-creating\n", workerID, sw.Sample, color.BlueString(vcf), vErr)
						if _, err := CreateGvcf(sw.Cram, resolvedFasta, item.chroms, gvcfPath, gatkLogLevel, caller, dvVer, modelType, verbose); err != nil {
							color.Red("[Worker %d] [%s] Error re-creating gVCF %s: %v\n", workerID, sw.Sample, color.BlueString(vcf), err)
							failedMu.Lock()
							failedTasks = append(failedTasks, FailedTask{Sample: sw.Sample, Chrom: label, Reason: err})
							failedMu.Unlock()
						} else {
							color.Green("[Worker %d] [%s] gVCF %s re-created successfully\n", workerID, sw.Sample, color.BlueString(vcf))
						}
					} else {
						color.Green("[Worker %d] [%s] gVCF %s is valid ✅\n\n", workerID, sw.Sample, vcf)
					}

				default:
					color.Red("[Worker %d] [%s] Multiple gVCF files found for %s — please remove extras.\n\n", workerID, sw.Sample, label)
					failedMu.Lock()
					failedTasks = append(failedTasks, FailedTask{
						Sample: sw.Sample, Chrom: label,
						Reason: fmt.Errorf("multiple gVCF files found"),
					})
					failedMu.Unlock()
				}
			}
		}(i)
	}

	for _, sw := range validSamples {
		for _, chrom := range chroms {
			workCh <- workItem{sw: sw, chroms: []SeqInfo{chrom}, label: chrom.ID}
		}
	}
	if len(contigs) > 0 {
		for _, sw := range validSamples {
			workCh <- workItem{sw: sw, chroms: contigs, label: "contigs", isContig: true}
		}
	}
	close(workCh)
	wg2.Wait()

	// ── summary ───────────────────────────────────────────────────────────────
	fmt.Printf("\n%s FINAL SUMMARY %s\n", strings.Repeat("=", 29), strings.Repeat("=", 29))
	fmt.Printf("Machine cores:     %d\n", totalCores)
	fmt.Printf("Threads per job:   %d\n", threadsPerJob)
	fmt.Printf("Max parallel jobs: %d\n", maxParallelJobs)
	fmt.Printf("Samples processed: %d\n", len(validSamples))
	fmt.Printf("Missing BAMs:      %d %v\n", len(missingBams), missingBams)

	if len(failedTasks) > 0 {
		color.Red("Failed gVCF tasks: %d\n", len(failedTasks))
		for _, f := range failedTasks {
			fmt.Printf("  - Sample: %s, Chrom/Label: %s, Reason: %v\n", f.Sample, f.Chrom, f.Reason)
		}
	} else {
		color.Green("All gVCF tasks completed successfully!\n")
	}
	fmt.Printf("%s\n\n", strings.Repeat("=", 69))

	if !noMerging {
		if len(missingBams) == 0 && len(failedTasks) == 0 {
			color.Cyan("\nNo missing BAMs and no failed samples. Ready for merging.\n")
			// TODO: implement merge call here
		} else {
			color.Yellow("Skipping merge due to missing BAMs (%d) or failed tasks (%d)\n", len(missingBams), len(failedTasks))
		}
	}

	if len(failedTasks) > 0 {
		return fmt.Errorf("%d gVCF task(s) failed", len(failedTasks))
	}
	return nil
}
