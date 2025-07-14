package bsaseq

import (
	"fmt"
	"golang.org/x/exp/rand"
	"log"
	"os"
	"path/filepath"
	"time"
)

func createResultsDir(outputDir string) (string, error) {

	baseDir := filepath.Join(outputDir, "goBSAseqResults")
	bErr := os.MkdirAll(filepath.Join(outputDir, "goBSAseqResults"), 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
		return "", bErr
	}

	now := time.Now()
	resultsDir := filepath.Join(baseDir, fmt.Sprintf("%02d_%02d_%04d_%02d_%02d_%02d", now.Day(), now.Month(), now.Year(), now.Hour(), now.Minute(), now.Second()))

	err := os.MkdirAll(resultsDir, 0755)
	if err != nil {
		log.Fatalf("Error creating results directory: %s\n", err)
		return "", err

	}
	fmt.Printf("Created results directory at %s ..\n\n", resultsDir)

	return resultsDir, nil
}

func simulateAF(popStruc string, bulkSize float64, rep int) float64 {

	var prob []float64

	switch popStruc {
	case "F2":
		prob = []float64{0.25, 0.5, 0.25}
	case "RIL":
		prob = []float64{0.5, 0.0, 0.5}
	case "BC":
		prob = []float64{0.5, 0.5, 0.0}
	default:
		fmt.Println("Invalid population structure")
		return 0.0
	}

	var totalFreq float64
	for i := 0; i < rep; i++ {
		var sumFreq float64
		for j := 0; j < int(bulkSize); j++ {
			// Use a simple weighted random choice
			r := rand.Float64()
			var allele float64
			if r < prob[0] {
				allele = 0.0
			} else if r < prob[0]+prob[1] {
				allele = 0.5
			} else {
				allele = 1.0
			}
			sumFreq += allele
		}
		totalFreq += sumFreq / bulkSize
	}
	return totalFreq / float64(rep)

}

func detectQtlPeaks(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {

	abovePoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] > y9Data[i] {
			abovePoints = append(abovePoints, i)
		}
	}

	if len(abovePoints) == 0 {
		return 0, 0, 0, 0, false
	}

	maxY := yData[abovePoints[0]]
	peakIndex := abovePoints[0]
	for _, idx := range abovePoints {
		if yData[idx] > maxY {
			maxY = yData[idx]
			peakIndex = idx
		}
	}
	peakX = x[peakIndex]
	peakY = maxY

	var intersections []int
	for i := 1; i < len(yData); i++ {

		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) || (yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}
	leftIntersect := -1
	rightIntersect := -1

	for _, idx := range intersections {
		if idx < peakIndex && (leftIntersect == -1 || idx > leftIntersect) {
			leftIntersect = idx
		} else if idx > peakIndex && (rightIntersect == -1 || idx < rightIntersect) {
			rightIntersect = idx
		}
	}

	if leftIntersect == -1 || rightIntersect == -1 {
		return peakX, peakY, 0, 0, false
	}

	return peakX, peakY, x[leftIntersect], x[rightIntersect], true
}

func detectQtlValleys(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {

	belowPoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] < y9Data[i] {
			belowPoints = append(belowPoints, i)
		}
	}

	if len(belowPoints) == 0 {
		return 0, 0, 0, 0, false
	}

	minY := yData[belowPoints[0]]
	peakIndex := belowPoints[0]
	for _, idx := range belowPoints {
		if yData[idx] < minY {
			minY = yData[idx]
			peakIndex = idx
		}
	}
	peakX = x[peakIndex]
	peakY = minY

	var intersections []int
	for i := 1; i < len(yData); i++ {
		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) ||
			(yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}

	leftIntersect := -1
	rightIntersect := -1

	for _, idx := range intersections {
		if idx < peakIndex && (leftIntersect == -1 || idx > leftIntersect) {
			leftIntersect = idx
		} else if idx > peakIndex && (rightIntersect == -1 || idx < rightIntersect) {
			rightIntersect = idx
		}
	}

	if leftIntersect == -1 || rightIntersect == -1 {
		return peakX, peakY, 0, 0, false
	}

	return peakX, peakY, x[leftIntersect], x[rightIntersect], true
}
