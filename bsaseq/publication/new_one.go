package main

import (
    "bufio"
    "fmt"
    "math"
    "os"
    "sort"
    "strconv"
    "strings"
)

/* ============================
   DATA STRUCTURES
============================ */

type SNP struct {
    Chr string
    Pos int

    AH, RH int
    AL, RL int

    SNPiH float64
    SNPiL float64
    Delta float64
    ED    float64
    G     float64

    DeltaK float64
    EDK    float64
    Gprime float64
}

type Block struct {
    Chr   string
    Start int
    End   int
    AFD   float64
    Z     float64
    Sig   bool
}

/* ============================
   BASIC STATISTICS
============================ */

func snpIndex(a, r int) float64 {
    if a+r == 0 {
        return math.NaN()
    }
    return float64(a) / float64(a+r)
}

func euclideanDelta(pH, pL float64) float64 {
    return math.Sqrt2 * math.Abs(pH-pL)
}

func gStatistic(aH, rH, aL, rL int) float64 {
    N := float64(aH + rH + aL + rL)
    if N == 0 {
        return 0
    }

    rowH := float64(aH + rH)
    rowL := float64(aL + rL)
    colA := float64(aH + aL)
    colR := float64(rH + rL)

    exp := func(r, c float64) float64 {
        return r * c / N
    }

    var G float64
    cells := []struct{ o, e float64 }{
        {float64(aH), exp(rowH, colA)},
        {float64(rH), exp(rowH, colR)},
        {float64(aL), exp(rowL, colA)},
        {float64(rL), exp(rowL, colR)},
    }

    for _, c := range cells {
        if c.o > 0 {
            G += 2 * c.o * math.Log(c.o/c.e)
        }
    }
    return G
}

/* ============================
   KERNELS
============================ */

func biweight(x float64) float64 {
    if x >= 1 {
        return 0
    }
    t := 1 - x*x
    return t * t
}

func tricube(x float64) float64 {
    if x >= 1 {
        return 0
    }
    return math.Pow(1-math.Pow(x, 3), 3)
}

/* ============================
   VCF PARSING
============================ */

func parseVCF(path string) ([]SNP, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, err
    }
    defer f.Close()

    var snps []SNP
    sc := bufio.NewScanner(f)

    for sc.Scan() {
        l := sc.Text()
        if strings.HasPrefix(l, "#") {
            continue
        }

        flds := strings.Split(l, "\t")
        chr := flds[0]
        pos, _ := strconv.Atoi(flds[1])
        format := strings.Split(flds[8], ":")
        hi := strings.Split(flds[9], ":")
        lo := strings.Split(flds[10], ":")

        adIdx := -1
        for i, f := range format {
            if f == "AD" {
                adIdx = i
                break
            }
        }
        if adIdx == -1 {
            continue
        }

        parseAD := func(s string) (a, r int) {
            p := strings.Split(s, ",")
            r, _ = strconv.Atoi(p[0])
            a, _ = strconv.Atoi(p[1])
            return
        }

        aH, rH := parseAD(hi[adIdx])
        aL, rL := parseAD(lo[adIdx])

        pH := snpIndex(aH, rH)
        pL := snpIndex(aL, rL)

        snps = append(snps, SNP{
            Chr:   chr,
            Pos:   pos,
            AH:    aH, RH: rH,
            AL:    aL, RL: rL,
            SNPiH: pH,
            SNPiL: pL,
            Delta: pH - pL,
            ED:    euclideanDelta(pH, pL),
            G:     gStatistic(aH, rH, aL, rL),
        })
    }

    return snps, sc.Err()
}

/* ============================
   KERNEL SMOOTHING (PER CHR)
============================ */

func smoothChromosome(snps []SNP, bandwidth int) {
    n := len(snps)

    for i := 0; i < n; i++ {
        center := snps[i].Pos
        var numD, numE, numG float64
        var denD, denE, denG float64

        for j := i; j >= 0; j-- {
            d := float64(center - snps[j].Pos)
            if d > float64(bandwidth) {
                break
            }
            x := d / float64(bandwidth)
            wB := biweight(x)
            wT := tricube(x)

            numD += wB * snps[j].Delta
            numE += wB * snps[j].ED
            numG += wT * snps[j].G
            denD += wB
            denE += wB
            denG += wT
        }

        for j := i + 1; j < n; j++ {
            d := float64(snps[j].Pos - center)
            if d > float64(bandwidth) {
                break
            }
            x := d / float64(bandwidth)
            wB := biweight(x)
            wT := tricube(x)

            numD += wB * snps[j].Delta
            numE += wB * snps[j].ED
            numG += wT * snps[j].G
            denD += wB
            denE += wB
            denG += wT
        }

        if denD > 0 {
            snps[i].DeltaK = numD / denD
            snps[i].EDK = numE / denE
        }
        if denG > 0 {
            snps[i].Gprime = numG / denG
        }
    }
}

/* ============================
   BRM BLOCK REGRESSION
============================ */

func brmBlocks(snps []SNP, blockSize int, zThresh float64) []Block {
    byChr := map[string][]SNP{}
    for _, s := range snps {
        byChr[s.Chr] = append(byChr[s.Chr], s)
    }

    var blocks []Block

    for chr, list := range byChr {
        sort.Slice(list, func(i, j int) bool {
            return list[i].Pos < list[j].Pos
        })

        for i := 0; i < len(list); {
            start := list[i].Pos
            end := start + blockSize
            var sum float64
            var n int

            for i < len(list) && list[i].Pos < end {
                sum += list[i].Delta
                n++
                i++
            }
            if n == 0 {
                continue
            }

            afd := sum / float64(n)
            z := afd / math.Sqrt(1.0/float64(n))

            blocks = append(blocks, Block{
                Chr:   chr,
                Start: start,
                End:   end,
                AFD:   afd,
                Z:     z,
                Sig:   math.Abs(z) >= zThresh,
            })
        }
    }
    return blocks
}

/* ============================
   MAIN
============================ */

func main() {
    if len(os.Args) < 2 {
        fmt.Println("usage: bsa_stats input.vcf")
        os.Exit(1)
    }

    // PARAMETERS (easy to promote to CLI flags)
    smoothWindow := 50000  // 50 kb smoothing
    blockSize := 100000    // 100 kb BRM blocks
    zThreshold := 3.0      // conservative significance

    snps, err := parseVCF(os.Args[1])
    if err != nil {
        panic(err)
    }

    // group and smooth per chromosome
    byChr := map[string][]SNP{}
    for _, s := range snps {
        byChr[s.Chr] = append(byChr[s.Chr], s)
    }
    for chr := range byChr {
        sort.Slice(byChr[chr], func(i, j int) bool {
            return byChr[chr][i].Pos < byChr[chr][j].Pos
        })
        smoothChromosome(byChr[chr], smoothWindow)
    }

    // flatten SNPs
    var allSnps []SNP
    for _, s := range byChr {
        allSnps = append(allSnps, s...)
    }

    // write SNP table
    snpsOut, _ := os.Create("snps.tsv")
    defer snpsOut.Close()
    fmt.Fprintln(snpsOut, "CHR\tPOS\tSNPiH\tSNPiL\tDELTA\tDELTAK\tED\tEDK\tG\tGPRIME")
    for _, s := range allSnps {
        fmt.Fprintf(snpsOut, "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\n",
            s.Chr, s.Pos, s.SNPiH, s.SNPiL,
            s.Delta, s.DeltaK, s.ED, s.EDK, s.G, s.Gprime)
    }

    // BRM blocks
    blocks := brmBlocks(allSnps, blockSize, zThreshold)
    blkOut, _ := os.Create("brm_blocks.tsv")
    defer blkOut.Close()
    fmt.Fprintln(blkOut, "CHR\tSTART\tEND\tAFD\tZ\tSIGNIFICANT")
    for _, b := range blocks {
        fmt.Fprintf(blkOut, "%s\t%d\t%d\t%.4f\t%.2f\t%v\n",
            b.Chr, b.Start, b.End, b.AFD, b.Z, b.Sig)
    }
}