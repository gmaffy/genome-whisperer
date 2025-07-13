package main

import (
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/gmaffy/genome-whisperer/utils"
)

func TestParseLog(t *testing.T) {
	// Create a temporary log file with the example entries
	logContent := `{"time":"2025-06-18T21:11:02.572267197+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"INITIALISE","SAMPLE":"ALL","CHROMOSOME":"ALL","STATUS":"STARTED","CMD":"ALL"}
{"time":"2025-06-18T21:11:03.397122518+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr01","STATUS":"STARTED"}
{"time":"2025-06-18T21:11:04.124962114+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr02","STATUS":"STARTED"}
{"time":"2025-06-18T21:11:05.01947693+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr03","STATUS":"STARTED"}
{"time":"2025-06-18T21:11:06.687393372+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr04","STATUS":"STARTED"}
{"time":"2025-06-18T21:20:17.308876904+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr02","STATUS":"COMPLETED"}
{"time":"2025-06-18T21:20:17.310433516+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"SEOL_TAGADI.RGMD.bam","CHROMOSOME":"Cmo_Chr02","STATUS":"STARTED"}
{"time":"2025-06-18T21:23:58.626151562+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"INITIALISE","SAMPLE":"ALL","CHROMOSOME":"ALL","STATUS":"STARTED","CMD":"ALL"}
{"time":"2025-06-18T21:23:58.952009702+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr01","STATUS":"STARTED"}
{"time":"2025-06-18T21:23:59.23049438+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr02","STATUS":"STARTED"}
{"time":"2025-06-18T21:23:59.644223773+02:00","level":"INFO","msg":"VARIANT CALLING","PROGRAM":"HaplotypeCaller","SAMPLE":"NIGERIAN_LOCAL.RGMD.bam","CHROMOSOME":"Cmo_Chr03","STATUS":"STARTED"}`

	// Create a temporary directory for the test
	tempDir := filepath.Join(os.TempDir(), "genome-whisperer-test")
	err := os.MkdirAll(tempDir, 0755)
	if err != nil {
		fmt.Printf("Error creating temp directory: %v\n", err)
		return
	}
	defer os.RemoveAll(tempDir)

	// Create the log file
	logFilePath := filepath.Join(tempDir, "test.log")
	err = os.WriteFile(logFilePath, []byte(logContent), 0644)
	if err != nil {
		fmt.Printf("Error writing log file: %v\n", err)
		return
	}

	// Parse the log file
	logEntries := utils.ParseLogFile(logFilePath)

	// Print the results
	fmt.Printf("Found %d log entries\n", len(logEntries))
	for i, entry := range logEntries {
		fmt.Printf("Entry %d:\n", i+1)
		fmt.Printf("  Timestamp: %s\n", entry.Timestamp)
		fmt.Printf("  Tool: %s\n", entry.Tool)
		fmt.Printf("  Program: %s\n", entry.Program)
		fmt.Printf("  Sample: %s\n", entry.Sample)
		fmt.Printf("  Chromosome: %s\n", entry.Chromosome)
		fmt.Printf("  Status: %s\n", entry.Status)
		fmt.Printf("  Cmd: %s\n", entry.Cmd)
		fmt.Println()
	}

	// Test the StageHasCompleted function
	completed := utils.StageHasCompleted(logEntries, "HaplotypeCaller", "NIGERIAN_LOCAL.RGMD.bam", "Cmo_Chr02")
	fmt.Printf("HaplotypeCaller for NIGERIAN_LOCAL.RGMD.bam on Cmo_Chr02 completed: %v\n", completed)

	notCompleted := utils.StageHasCompleted(logEntries, "HaplotypeCaller", "NIGERIAN_LOCAL.RGMD.bam", "Cmo_Chr03")
	fmt.Printf("HaplotypeCaller for NIGERIAN_LOCAL.RGMD.bam on Cmo_Chr03 completed: %v\n", notCompleted)
}