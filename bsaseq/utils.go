package bsaseq

import (
	"fmt"
	"golang.org/x/exp/rand"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"time"
)

func CheckDeps() error {
	deps := []string{"gatk", "samtools", "bwa"}

	for _, dep := range deps {
		if _, err := exec.LookPath(dep); err != nil {
			return fmt.Errorf("%s not found: %w", dep, err)
		}
		fmt.Printf("%s OK\n", dep)
	}

	return nil
}

func createResultsDir() string {
	bErr := os.MkdirAll("goBSAseqResults", 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
	}

	baseDir := "goBSAseqResults"
	now := time.Now()
	resultsDir := filepath.Join(baseDir, fmt.Sprintf("%02d_%02d_%04d_%02d_%02d_%02d", now.Day(), now.Month(), now.Year(), now.Hour(), now.Minute(), now.Second()))

	err := os.MkdirAll(resultsDir, 0755)
	if err != nil {
		log.Fatalf("Error creating results directory: %s\n", err)

	}
	fmt.Printf("Created results directory at %s ..\n\n", resultsDir)

	return resultsDir
}

func simulateAF(popStruc string, bulkSize float64, rep int) float64 {

	var prob []float64

	// Define probabilities based on population structure
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

func FindIntersections(x []int, y1 []float64, y2 []float64) []int {
	intersections := []int{}

	for i := 0; i < len(x)-1; i++ {
		diff1 := y1[i] - y2[i]
		diff2 := y1[i+1] - y2[i+1]

		// Check if there's a sign change (i.e., crossing)
		if diff1*diff2 <= 0 {
			// Linear interpolation to estimate crossing point
			dx := float64(x[i+1] - x[i])
			if dx == 0 {
				continue // skip vertical lines
			}

			dy1 := diff1
			dy2 := diff2
			t := dy1 / (dy1 - dy2) // fraction along the segment
			xIntersect := float64(x[i]) + t*dx
			intersections = append(intersections, int(xIntersect))
		}
	}

	return intersections
}

/*
	func MaxIndex(y []float64) (int, float64) {
		maxIdx := 0
		maxVal := y[0]
		for i := 1; i < len(y); i++ {
			if y[i] > maxVal {
				maxVal = y[i]
				maxIdx = i
			}
		}
		return maxIdx, maxVal
	}

	func MinIndex(y []float64) (int, float64) {
		minIdx := 0
		minVal := y[0]
		for i := 1; i < len(y); i++ {
			if y[i] < minVal {
				minVal = y[i]
				minIdx = i
			}
		}
		return minIdx, minVal
	}

	func FlankingIntersections(intersections []int, x []int, peakIdx int) (int, int) {
		peakPos := x[peakIdx]

		left := -1
		right := -1

		for i := 0; i < len(intersections); i++ {
			if intersections[i] < peakPos {
				left = intersections[i]
			} else if intersections[i] > peakPos {
				right = intersections[i]
				break
			}
		}

		return left, right
	}

func detectQtlPeaks(x []int, y []float64, qtlThresh []float64) (float64, int, int, bool) {

	deltaWith99 := FindIntersections(x, y, qtlThresh)
	peakX, peakY := MaxIndex(y)
	qtlStart, qtlEnd := FlankingIntersections(deltaWith99, x, peakX)

	if qtlStart == -1 && qtlEnd == -1 {
		return 0, 0, 0.0, false
	} else {
		return peakY, qtlStart, qtlEnd, true

	}

}

	func detectQtlValleys(x []int, y []float64, qtlThresh []float64) (float64, int, int, bool) {
		deltaWith99 := FindIntersections(x, y, qtlThresh)
		peakX, peakY := MinIndex(y)
		qtlStart, qtlEnd := FlankingIntersections(deltaWith99, x, peakX)

		if qtlStart == -1 && qtlEnd == -1 {
			return 0, 0, 0.0, false
		} else {
			return peakY, qtlStart, qtlEnd, true

		}
	}
*/

func detectQtlPeaks(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {
	// First find all points where yData is above y9Data
	abovePoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] > y9Data[i] {
			abovePoints = append(abovePoints, i)
		}
	}

	if len(abovePoints) == 0 {
		return 0, 0, 0, 0, false
	}

	// Find the highest point among these
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
		// Check if the lines crossed between i-1 and i
		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) || (yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			// Approximate the intersection by taking the current point

			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}

	// Find the two intersections closest to the peak (one on each side)
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

/*
	func detectQtlPeaks(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX, rightIntersectX float64, found bool) {
		// First find all points where yData is above y9Data
		abovePoints := make([]int, 0)
		for i := 0; i < len(yData); i++ {
			if yData[i] > y9Data[i] {
				abovePoints = append(abovePoints, i)
			}
		}

		if len(abovePoints) == 0 {
			return 0, 0, 0, 0, false
		}

		// Find the highest point among these
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

		// Find intersections (where yData crosses y9Data)
		var intersections []float64
		for i := 1; i < len(yData); i++ {
			// Check if the lines crossed between i-1 and i
			if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) || (yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
				// Calculate the intersection point using linear interpolation
				x1, y1 := float64(x[i-1]), yData[i-1]-y9Data[i-1]
				x2, y2 := float64(x[i]), yData[i]-y9Data[i]
				intersectionX := x1 - y1*(x2-x1)/(y2-y1)
				intersections = append(intersections, intersectionX)
			}
		}

		if len(intersections) < 2 {
			return peakX, peakY, 0, 0, false
		}

		// Find the two intersections closest to the peak (one on each side)
		leftIntersect := math.Inf(1)
		rightIntersect := math.Inf(-1)

		for _, intersectX := range intersections {
			if intersectX < float64(peakX) && intersectX > leftIntersect {
				leftIntersect = intersectX
			} else if intersectX > float64(peakX) && intersectX < rightIntersect {
				rightIntersect = intersectX
			}
		}

		if leftIntersect == math.Inf(1) || rightIntersect == math.Inf(-1) {
			return peakX, peakY, 0, 0, false
		}

		return peakX, peakY, leftIntersect, rightIntersect, true
	}
*/
func detectQtlValleys(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {
	// First find all points where yData is below y9Data
	belowPoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] < y9Data[i] {
			belowPoints = append(belowPoints, i)
		}
	}

	if len(belowPoints) == 0 {
		return 0, 0, 0, 0, false
	}

	// Find the lowest point among these
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

	// Find intersections (where yData crosses y9Data)
	var intersections []int
	for i := 1; i < len(yData); i++ {
		// Check if the lines crossed between i-1 and i
		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) ||
			(yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			// Approximate the intersection by taking the current point
			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}

	// Find the two intersections closest to the peak (one on each side)
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
