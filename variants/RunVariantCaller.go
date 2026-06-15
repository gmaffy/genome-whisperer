package variants

import (
	"bufio"
	"compress/gzip"
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

	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
	"github.com/fatih/color"

	"github.com/gmaffy/genome-whisperer/utils"
)

// ─────────────────────────────────────────────────────────────────────────────
// Stage-name constants (defined once, used everywhere to avoid typos)
// ─────────────────────────────────────────────────────────────────────────────

const (
	stageHaplotypeCaller  = "HaplotypeCaller"
	stageDeepVariant      = "DEEPVARIANT"
	stageGenomicsDBImport = "GenomicsDBImport"
	stageGenotypeGVCFs    = "GenotypeGVCFs"
	stageGLNEXUS          = "GLNEXUS"
	stageGLNEXUSBCFTOOLS  = "GLNEXUS_BCFTOOLS"
	stageGLNEXUSIndex     = "GLNEXUS_INDEX"
	stageHardFilter       = "HardFilterVcf"
	stageMergeVcfs        = "MergeVcfs"
)

type SeqInfo struct {
	ID  string
	Len int
}

func CreateGvcfGATK(bam string, refFile string, chroms []SeqInfo, species string, refVer string, logLevel string, verbose bool, out string, threadsPerJob int) (string, error) {

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
		gVCF = filepath.Join(gvcfDir, species+refVer+".Contigs.g.vcf.gz")

	}

	hapCmdStr := fmt.Sprintf(
		`gatk HaplotypeCaller -R %s -I %s -L %s -O %s -ERC GVCF --verbosity %s --native-pair-hmm-threads %d`, refFile, bam, regions, gVCF, logLevel, threadsPerJob,
	)
	fmt.Printf("\n%s\n\n", hapCmdStr)
	return gVCF, utils.RunCmd(hapCmdStr, verbose)

}

func CreateGvcfDV(bam string, refFile string, chroms []SeqInfo, species, refVer, dvVer string, modelType string, verbose bool, out string, threadsPerJob int) (string, error) {

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

	//dvCmdStr := fmt.Sprintf(
	//	`docker run -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
	//		`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
	//		`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
	//		`--intermediate_results_dir /output/%s`,
	//	bamDir, refDir, out, dvVer,
	//	modelType, refName, bamName,
	//	regions, vcfName, gvcfName,
	//	intermediateName,
	//)
	dvCmdStr := fmt.Sprintf(
		`docker run -v "%s":/bam -v "%s":/ref -v "%s":/output google/deepvariant:%s `+
			`/opt/deepvariant/bin/run_deepvariant --model_type=%s --ref=/ref/%s --reads=/bam/%s `+
			`--regions "%s" --output_vcf=/output/%s --output_gvcf=/output/%s `+
			`--num_shards=%d --intermediate_results_dir /output/%s`,
		bamDir, refDir, out, dvVer,
		modelType, refName, bamName,
		regions, vcfName, gvcfName,
		threadsPerJob, intermediateName,
	)
	fmt.Printf("\n%s\n\n", dvCmdStr)
	return gVCF, utils.RunCmd(dvCmdStr, verbose)
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

func MergeGvcfsGATK(gvcfs []string, refFile string, chrom string, species string, refVer string, verbose bool, logLevel string, outDir string, logged []utils.LogEntry, jlog *slog.Logger) (string, error) {

	theDB := filepath.Join(outDir, chrom+"DB")
	tmpDir := filepath.Join(outDir, "tmp")
	tmpDir2 := filepath.Join(outDir, "tmp2")
	vcfDir := filepath.Join(outDir, "VCFs")

	for _, dir := range []string{tmpDir, tmpDir2, vcfDir} {
		if _, err := os.Stat(dir); os.IsNotExist(err) {
			if cErr := os.MkdirAll(dir, 0755); cErr != nil {
				return "", fmt.Errorf("creating directory %s: %w", dir, cErr)
			}
		}
	}

	jointVCF := filepath.Join(vcfDir, species+refVer+"."+chrom+".joint.vcf.gz")

	// ── GenomicsDBImport ────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGenomicsDBImport, "ALL", chrom) {
		if rErr := os.RemoveAll(theDB); rErr != nil {
			return "", fmt.Errorf("removing %s: %w", theDB, rErr)
		}

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
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", gDBCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(gDBCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("gatk GenomicsDBImport: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenomicsDBImport, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenomicsDBImport, chrom))
	}

	// ── GenotypeGVCFs ───────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGenotypeGVCFs, "ALL", chrom) {
		genoCmd := fmt.Sprintf(
			`gatk --java-options "-Xmx12g" GenotypeGVCFs -R %s -V gendb://%s -O %s --tmp-dir %s --verbosity %s`,
			refFile, theDB, jointVCF, tmpDir2, logLevel,
		)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", genoCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(genoCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("gatk GenotypeGVCFs: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGenotypeGVCFs, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGenotypeGVCFs, chrom))
	}

	return jointVCF, nil
}

func MergeGvcfsGlnexus(gvcfs []string, chrom string, species string, refVer string, caller string, verbose bool, outDir string, logged []utils.LogEntry, jlog *slog.Logger) (string, error) {
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

	// ── GLnexus ──────────────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGLNEXUS, "ALL", chrom) {
		glnexusCmd := fmt.Sprintf(`glnexus_cli --config %s --dir %s %s > %s`, glnexusPreset, filepath.Join(vcfDir, "GLnexus.DB"), gvcfFiles, jointBCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", glnexusCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(glnexusCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("glnexus_cli: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUS, chrom))
	}

	// ── bcftools view | bgzip ──────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGLNEXUSBCFTOOLS, "ALL", chrom) {
		bcfCmd := fmt.Sprintf("bcftools view %s | bgzip -@ 4 -c > %s", jointBCF, jointVCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", bcfCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(bcfCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("bcftools view: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSBCFTOOLS, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSBCFTOOLS, chrom))
	}

	// ── tabix index ──────────────────────────────────────────────────────────
	if !utils.StageHasCompleted(logged, stageGLNEXUSIndex, "ALL", chrom) {
		indexCmd := fmt.Sprintf(`tabix -p vcf %s`, jointVCF)
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED", "CMD", indexCmd)
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "STARTED")

		if err := utils.RunCmd(indexCmd, verbose); err != nil {
			jlog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			slog.Error("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", fmt.Sprintf("FAILED: %v", err))
			return "", fmt.Errorf("indexing joint VCF with tabix: %w", err)
		}
		jlog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		slog.Info("VARIANT CALLING", "PROGRAM", stageGLNEXUSIndex, "SAMPLE", "ALL", "CHROMOSOME", chrom, "STATUS", "COMPLETED")
		fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
	} else {
		slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageGLNEXUSIndex, chrom))
	}

	return jointVCF, nil
}

// ----------------------------------------- Filter --------------------------------------------------------------- //

func classifyVariant(v *vcfgo.Variant) (isSNP bool, isIndel bool, isMNP bool) {
	refLen := len(v.Ref())
	isSNP = true
	isIndel = false
	isMNP = true

	if refLen != 1 {
		isSNP = false
	}

	for _, alt := range v.Alt() {
		altLen := len(alt)
		// If lengths differ, it's an Indel
		if altLen != refLen {
			isIndel = true
			isMNP = false
			isSNP = false
		}
	}

	if isIndel {
		isMNP = false
		isSNP = false
	}

	if refLen == 1 {
		isMNP = false
	}

	return isSNP, isIndel, isMNP
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	isSNP, isIndel, isMNP := classifyVariant(v)

	// Pre-fetch all commonly used INFO fields
	qd, hasQD := utils.GetFloat(v, "QD")
	fs, hasFS := utils.GetFloat(v, "FS")
	sor, hasSOR := utils.GetFloat(v, "SOR")
	mq, hasMQ := utils.GetFloat(v, "MQ")
	mqRankSum, hasMQRankSum := utils.GetFloat(v, "MQRankSum")
	readPosRankSum, hasReadPosRankSum := utils.GetFloat(v, "ReadPosRankSum")

	// vcfgo stores Quality directly on the struct
	qual := float64(v.Quality)

	switch {
	case isSNP:
		if qual < hfcfg.SNP_QUAL_Min {
			return false
		}
		if hasQD && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if hasFS && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if hasSOR && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if hasMQ && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if hasMQRankSum && mqRankSum < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if hasReadPosRankSum && readPosRankSum < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel, isMNP:
		// We group MNPs with Indels here to hold them to the same robust standards.
		if qual < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if hasQD && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if hasFS && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if hasSOR && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if hasReadPosRankSum && readPosRankSum < hfcfg.INDEL_ReadPosRankSum_Min {
			return false
		}
		return true

	default:
		// these may be SVs. so we keep
		return true
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
		err = utils.RunCmd(cmd, verbose)

		if err != nil {
			return "", fmt.Errorf("gatk MergeVcfs error: %w", err)
		}
	}
	return filepath.Join(outDir, concatVcfName), nil
}

func VariantCalling(bams []string, refFile string, species string, refVer string, out string, caller string, merger string, noMerging bool, noHardFilter bool, cfg utils.HardFilterConfig, verbose bool, gatkLogLevel string, threads int, logFilePath string, dvVer string, modelType string) (string, error) {

	dictFilePath := refFile[:len(refFile)-len(filepath.Ext(refFile))] + ".dict"
	if _, dicfErr := os.Stat(dictFilePath); dicfErr != nil {
		return "", fmt.Errorf("reference dict file %s does not exist", dictFilePath)
	}

	chroms, contigs, err := getChromsAndContigs(dictFilePath)
	if err != nil {
		return "", fmt.Errorf("getting chroms and contigs from dict file: %w", err)
	}

	// ----------------------------------------- open / create log file --------------------------------------------- //
	fmt.Println("Opening log file ...", logFilePath)
	logFile, err := os.OpenFile(logFilePath, os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		return "", fmt.Errorf("opening log file %s: %w", logFilePath, err)
	}
	defer logFile.Close()

	jlog := slog.New(slog.NewJSONHandler(logFile, nil))

	jlog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")
	slog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "STARTED")

	// ── resume check ─────────────────────────────────────────────────────────
	logged := utils.ParseLogFile(logFilePath)
	if !noMerging && utils.StageHasCompleted(logged, stageMergeVcfs, "ALL", "ALL") {
		fmt.Println("All stages completed. Skipping.")
		return filepath.Join(out, fmt.Sprintf("%s.%s.all.vcf.gz", strings.ToUpper(species), strings.ToLower(refVer))), nil
	}

	finalVCF := ""
	switch strings.ToLower(caller) {
	case "gatk":
		fmt.Printf("Using GATK HaplotypeCaller with log level %s\n", gatkLogLevel)
		totalThreads := runtime.NumCPU()
		color.Cyan("# Avalilable threads: %v\n# Jobs run (with %v threads per job): %v\n", totalThreads, threads, totalThreads/threads)
		maxWorkers := runtime.NumCPU() / threads
		if maxWorkers < 1 {
			maxWorkers = 2 // Fallback default
		}

		semaphore := make(chan struct{}, maxWorkers)
		var wg sync.WaitGroup
		errChan := make(chan error, len(chroms))
		var jointVcfs []string
		var mu sync.Mutex

		for _, chrom := range chroms {
			wg.Add(1)
			semaphore <- struct{}{}
			go func(c SeqInfo) {
				defer wg.Done()
				defer func() { <-semaphore }()
				// ------------------------------------- Create Gvcfs --------------------------------------------------------- //
				var gvcfs []string
				for _, bam := range bams {
					bamName := filepath.Base(bam)

					if utils.StageHasCompleted(logged, stageHaplotypeCaller, bamName, c.ID) {
						slog.Info(fmt.Sprintf("%s already completed for %s / %s. Skipping.", stageHaplotypeCaller, bamName, chrom.ID))
						chromDirName := strings.ReplaceAll(c.ID, ".", "_")
						chromDir := filepath.Join(out, chromDirName)
						gvcfDir := filepath.Join(chromDir, "gvcfs")
						chromGVCF := filepath.Join(gvcfDir, species+refVer+"."+c.ID+".g.vcf.gz")
						gvcfs = append(gvcfs, chromGVCF)
						continue
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", "STARTED")

					chromGVCF, err := CreateGvcfGATK(bam, refFile, []SeqInfo{chrom}, species, refVer, gatkLogLevel, verbose, out, threads)
					if err != nil {
						//errChan <- fmt.Errorf("failed on contig %s: %w", chrom.ID, err)
						jlog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", fmt.Sprintf("FAILED: %v", err))
						slog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", fmt.Sprintf("FAILED: %v", err))
						errChan <- fmt.Errorf("creating gVCF for %s with GATK HaplotypeCaller: %w", chrom.ID, err)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", c.ID, "STATUS", "COMPLETED")
					fmt.Printf("Created gVCF: %s\n", chromGVCF)
					fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
					gvcfs = append(gvcfs, chromGVCF)
				}

				if noMerging {
					fmt.Println("Skipping merging of gVCFs")

				}

				// ------------------------------------- Merge gvcfs ---------------------------------------------------------- //
				var finalChromVcf string
				var mergedChromVcf string
				switch merger {
				case "gatk":
					mergedChromVcf, err = MergeGvcfsGATK(gvcfs, refFile, chrom.ID, species, refVer, verbose, gatkLogLevel, out, logged, jlog)
					if err != nil {
						errChan <- fmt.Errorf("merging gVCFs for chrom %s with GATK: %w", chrom.ID, err)
					}
					fmt.Printf("Merged gVCFs for chromosome %s into joint VCF: %s\n", chrom.ID, mergedChromVcf)
					finalChromVcf = mergedChromVcf
				case "glnexus":
					mergedChromVcf, err = MergeGvcfsGlnexus(gvcfs, chrom.ID, species, refVer, caller, verbose, out, logged, jlog)
					finalChromVcf = mergedChromVcf
					if err != nil {
						errChan <- fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedChromVcf)
				default:
					errChan <- fmt.Errorf("unsupported merger: %s", merger)
				}

				// ----------------------------------------- Hard filter vcf ------------------------------------------- //
				if !noHardFilter {
					hardFilteredVcf := strings.TrimSuffix(mergedChromVcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
					if !utils.StageHasCompleted(logged, stageHardFilter, "ALL", chrom.ID) {
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "STARTED")

						hfVcf, hfErr := HardFilterVcf(mergedChromVcf, cfg)
						if hfErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							errChan <- fmt.Errorf("hard filtering VCF: %w", hfErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
						hardFilteredVcf = hfVcf
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageHardFilter, chrom.ID))
					}
					finalChromVcf = hardFilteredVcf
				}
				mu.Lock()
				jointVcfs = append(jointVcfs, finalChromVcf)
				mu.Unlock()
			}(chrom)

		}
		wg.Wait()
		close(errChan)

		for err := range errChan {
			if err != nil {
				return "", fmt.Errorf("pipeline aborted due to worker error: %w", err)
			}
		}

		if len(contigs) > 0 {
			// ------------------------------------- Create Gvcfs --------------------------------------------------- //
			contGvcfs := []string{}
			for _, bam := range bams {
				bamName := filepath.Base(bam)

				if utils.StageHasCompleted(logged, stageHaplotypeCaller, bamName, "contigs") {
					slog.Info(fmt.Sprintf("%s already completed for %s / contigs. Skipping.", stageHaplotypeCaller, bamName))
					contigsGVCF, err := CreateGvcfGATK(bam, refFile, contigs, species, refVer, gatkLogLevel, verbose, out, threads)
					if err != nil {
						return "", fmt.Errorf("resolving existing gVCF path for %s / contigs: %w", bamName, err)
					}
					contGvcfs = append(contGvcfs, contigsGVCF)
					continue
				}

				jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "STARTED")
				slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "STARTED")

				contigsGVCF, err := CreateGvcfGATK(bam, refFile, contigs, species, refVer, gatkLogLevel, verbose, out, threads)
				if err != nil {
					jlog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", err))
					slog.Error("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", err))
					return "", fmt.Errorf("creating gVCF for contigs with GATK HaplotypeCaller: %w", err)
				}
				jlog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
				slog.Info("VARIANT CALLING", "PROGRAM", stageHaplotypeCaller, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
				fmt.Printf("Created gVCF: %s\n", contigsGVCF)
				fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
				contGvcfs = append(contGvcfs, contigsGVCF)
			}
			// ------------------------------------- Merge gvcfs ----------------------------------------------------- //
			var finalContigsVcf string
			var mergedContigsVcf string
			if !noMerging {
				switch merger {
				case "gatk":
					mergedContigsVcf, err = MergeGvcfsGATK(contGvcfs, refFile, "contigs", species, refVer, verbose, gatkLogLevel, out, logged, jlog)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GATK: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				case "glnexus":
					mergedContigsVcf, err = MergeGvcfsGlnexus(contGvcfs, "contigs", species, refVer, caller, verbose, out, logged, jlog)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				default:
					return "", fmt.Errorf("unsupported merger: %s", merger)
				}

				// ----------------------------------------- Hard filter vcf ------------------------------------------- //
				if !noHardFilter {
					hardFilteredVcf := strings.TrimSuffix(mergedContigsVcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
					if !utils.StageHasCompleted(logged, stageHardFilter, "ALL", "contigs") {
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "STARTED")

						hfVcf, hfErr := HardFilterVcf(mergedContigsVcf, cfg)
						if hfErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							return "", fmt.Errorf("hard filtering VCF: %w", hfErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
						hardFilteredVcf = hfVcf
					} else {
						slog.Info(fmt.Sprintf("%s already completed for contigs. Skipping.", stageHardFilter))
					}
					finalContigsVcf = hardFilteredVcf
				}
				jointVcfs = append(jointVcfs, finalContigsVcf)
			}
		}

		if !noMerging {
			finalVCF, err = ConcatenateVcfs(jointVcfs, species, refVer, out, verbose)
			if err != nil {
				jlog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
				slog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
				return "", fmt.Errorf("concatenating VCFs: %w", err)
			}
			jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

	case "deepvariant":
		fmt.Printf("Using DeepVariant %s with model type %s\n", dvVer, modelType)
		maxWorkers := threads
		if maxWorkers < 1 {
			maxWorkers = 2 // Fallback default
		}

		semaphore := make(chan struct{}, maxWorkers)
		var wg sync.WaitGroup
		errChan := make(chan error, len(chroms))
		var jointVcfs []string
		var mu sync.Mutex

		//jointVcfs := []string{}
		for _, chrom := range chroms {
			wg.Add(1)
			semaphore <- struct{}{}
			go func(c SeqInfo) {
				defer wg.Done()
				defer func() { <-semaphore }()
				// ------------------------------------- Create Gvcfs --------------------------------------------------------- //
				var gvcfs []string
				for _, bam := range bams {
					bamName := filepath.Base(bam)

					if utils.StageHasCompleted(logged, stageDeepVariant, bamName, chrom.ID) {
						slog.Info(fmt.Sprintf("%s already completed for %s / %s. Skipping.", stageDeepVariant, bamName, chrom.ID))
						chromGVCF, err := CreateGvcfDV(bam, refFile, []SeqInfo{chrom}, species, refVer, dvVer, modelType, verbose, out, threads)
						if err != nil {
							errChan <- fmt.Errorf("resolving existing gVCF path for %s / %s: %w", bamName, chrom.ID, err)
						}
						gvcfs = append(gvcfs, chromGVCF)
						continue
					}

					jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", "STARTED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", "STARTED")

					chromGVCF, err := CreateGvcfDV(bam, refFile, []SeqInfo{chrom}, species, refVer, dvVer, modelType, verbose, out, threads)
					if err != nil {
						jlog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", err))
						slog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", err))
						errChan <- fmt.Errorf("creating gVCF for %s with DeepVariant: %w", chrom.ID, err)
					}
					jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
					slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
					fmt.Printf("Created gVCF: %s\n", chromGVCF)
					fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
					gvcfs = append(gvcfs, chromGVCF)
				}

				if noMerging {
					fmt.Println("Skipping merging of gVCFs")
					//continue
				}

				// ------------------------------------- Merge gvcfs ---------------------------------------------------------- //
				var finalChromVcf string
				var mergedChromVcf string
				switch merger {
				case "gatk":
					mergedChromVcf, err = MergeGvcfsGATK(gvcfs, refFile, chrom.ID, species, refVer, verbose, gatkLogLevel, out, logged, jlog)
					if err != nil {
						errChan <- fmt.Errorf("merging gVCFs for chrom %s with GATK: %w", chrom.ID, err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedChromVcf)
					finalChromVcf = mergedChromVcf
				case "glnexus":
					mergedChromVcf, err = MergeGvcfsGlnexus(gvcfs, chrom.ID, species, refVer, caller, verbose, out, logged, jlog)
					if err != nil {
						errChan <- fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedChromVcf)
					finalChromVcf = mergedChromVcf
				default:
					errChan <- fmt.Errorf("unsupported merger: %s", merger)
				}

				// ----------------------------------------- Hard filter vcf ------------------------------------------- //
				if !noHardFilter {
					hardFilteredVcf := strings.TrimSuffix(mergedChromVcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
					if !utils.StageHasCompleted(logged, stageHardFilter, "ALL", chrom.ID) {
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "STARTED")

						hfVcf, hfErr := HardFilterVcf(mergedChromVcf, cfg)
						if hfErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							errChan <- fmt.Errorf("hard filtering VCF: %w", hfErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", chrom.ID, "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
						hardFilteredVcf = hfVcf
					} else {
						slog.Info(fmt.Sprintf("%s already completed for %s. Skipping.", stageHardFilter, chrom.ID))
					}
					finalChromVcf = hardFilteredVcf
				}
				mu.Lock()
				jointVcfs = append(jointVcfs, finalChromVcf)
				mu.Unlock()
			}(chrom)
		}
		wg.Wait()
		close(errChan)

		for err := range errChan {
			if err != nil {
				return "", fmt.Errorf("pipeline aborted due to worker error: %w", err)
			}
		}
		if len(contigs) > 0 {
			// ------------------------------------- Create Gvcfs --------------------------------------------------- //
			contGvcfs := []string{}
			for _, bam := range bams {
				bamName := filepath.Base(bam)

				if utils.StageHasCompleted(logged, stageDeepVariant, bamName, "contigs") {
					slog.Info(fmt.Sprintf("%s already completed for %s / contigs. Skipping.", stageDeepVariant, bamName))
					contigsGVCF, err := CreateGvcfDV(bam, refFile, contigs, species, refVer, dvVer, modelType, verbose, out, threads)
					if err != nil {
						return "", fmt.Errorf("resolving existing gVCF path for %s / contigs: %w", bamName, err)
					}
					contGvcfs = append(contGvcfs, contigsGVCF)
					continue
				}

				jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "STARTED")
				slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "STARTED")

				contigsGVCF, err := CreateGvcfDV(bam, refFile, contigs, species, refVer, dvVer, modelType, verbose, out, threads)
				if err != nil {
					jlog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", err))
					slog.Error("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", err))
					return "", fmt.Errorf("creating gVCF for contigs with DeepVariant: %w", err)
				}
				jlog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
				slog.Info("VARIANT CALLING", "PROGRAM", stageDeepVariant, "SAMPLE", bamName, "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
				fmt.Printf("Created gVCF: %s\n", contigsGVCF)
				fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
				contGvcfs = append(contGvcfs, contigsGVCF)
			}
			// ------------------------------------- Merge gvcfs ----------------------------------------------------- //
			var finalContigsVcf string
			var mergedContigsVcf string
			if !noMerging {
				switch merger {
				case "gatk":
					mergedContigsVcf, err = MergeGvcfsGATK(contGvcfs, refFile, "contigs", species, refVer, verbose, gatkLogLevel, out, logged, jlog)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GATK: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				case "glnexus":
					mergedContigsVcf, err = MergeGvcfsGlnexus(contGvcfs, "contigs", species, refVer, caller, verbose, out, logged, jlog)
					if err != nil {
						return "", fmt.Errorf("merging gVCFs with GLnexus: %w", err)
					}
					fmt.Printf("Merged gVCFs into joint VCF: %s\n", mergedContigsVcf)
					finalContigsVcf = mergedContigsVcf
				default:
					return "", fmt.Errorf("unsupported merger: %s", merger)
				}

				// ----------------------------------------- Hard filter vcf ------------------------------------------- //
				if !noHardFilter {
					hardFilteredVcf := strings.TrimSuffix(mergedContigsVcf, ".vcf.gz") + ".hard_filtered.vcf.gz"
					if !utils.StageHasCompleted(logged, stageHardFilter, "ALL", "contigs") {
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "STARTED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "STARTED")

						hfVcf, hfErr := HardFilterVcf(mergedContigsVcf, cfg)
						if hfErr != nil {
							jlog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							slog.Error("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", fmt.Sprintf("FAILED: %v", hfErr))
							return "", fmt.Errorf("hard filtering VCF: %w", hfErr)
						}
						jlog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
						slog.Info("VARIANT CALLING", "PROGRAM", stageHardFilter, "SAMPLE", "ALL", "CHROMOSOME", "contigs", "STATUS", "COMPLETED")
						fmt.Printf("\n\n%s\n\n", strings.Repeat("-", 90))
						hardFilteredVcf = hfVcf
					} else {
						slog.Info(fmt.Sprintf("%s already completed for contigs. Skipping.", stageHardFilter))
					}
					finalContigsVcf = hardFilteredVcf
				}
				jointVcfs = append(jointVcfs, finalContigsVcf)
			}
		}

		if !noMerging {
			finalVCF, err = ConcatenateVcfs(jointVcfs, species, refVer, out, verbose)
			if err != nil {
				jlog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
				slog.Error("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", fmt.Sprintf("FAILED: %v", err))
				return "", fmt.Errorf("concatenating VCFs: %w", err)
			}
			jlog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
			slog.Info("VARIANT CALLING", "PROGRAM", stageMergeVcfs, "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
		}

	default:
		return "", fmt.Errorf("unsupported caller: %s", caller)

	}

	jlog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")
	slog.Info("VARIANT CALLING", "PROGRAM", "INITIALISE", "SAMPLE", "ALL", "CHROMOSOME", "ALL", "STATUS", "COMPLETED")

	return finalVCF, nil
}
