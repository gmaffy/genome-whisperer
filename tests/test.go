package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"log/slog"
	"os"
)

func hasCompletedStatus(logFile string, targetChrom int) (bool, error) {
	file, err := os.Open(logFile)
	if err != nil {
		return false, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		var entry map[string]interface{}
		if err := json.Unmarshal(scanner.Bytes(), &entry); err != nil {
			continue // skip malformed line
		}

		chromVal, chromOk := entry["CHROM"]
		statusVal, statusOk := entry["STATUS"]

		// Check CHROM and STATUS
		if chromOk && statusOk {
			if intVal, ok := toInt(chromVal); ok && intVal == targetChrom {
				if statusStr, ok := statusVal.(string); ok && statusStr == "COMPLETED" {
					return true, nil
				}
			}
		}
	}

	return false, scanner.Err()
}

// Helper to convert float64 or int JSON value to int
func toInt(v interface{}) (int, bool) {
	switch val := v.(type) {
	case float64:
		return int(val), true
	case int:
		return val, true
	default:
		return 0, false
	}
}

func main() {

	logFile, err := os.OpenFile("log.json", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("failed to open log file: %v", err)
	}
	defer logFile.Close()

	// Create a JSON handler for the file
	jsonHandler := slog.NewJSONHandler(logFile, nil)
	fileLogger := slog.New(jsonHandler)

	// Log some entries
	fmt.Println("Logging entries...")
	fileLogger.Info("Hello", "CHROM", 1, "PROGRAM", "HAPCALLER", "STATUS", "STATED")
	fileLogger.Info("Hello", "CHROM", 2, "PROGRAM", "HAPCALLER", "STATUS", "STATED")
	fileLogger.Info("Hello", "CHROM", 1, "PROGRAM", "HAPCALLER", "STATUS", "COMPLETED") // This is the one we'll find
	fileLogger.Info("Hello", "CHROM", 3, "PROGRAM", "HAPCALLER", "STATUS", "STATED")
	fileLogger.Info("Hello", "CHROM", 2, "PROGRAM", "HAPCALLER", "STATUS", "COMPLETED")
	fileLogger.Info("Goodbye", "ExtraKey", "ExtraValue") // Example of a log with different fields

	completed, err := hasCompletedStatus("log.json", 3)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}
	if completed {
		fmt.Println("CHROM 1 is COMPLETED")
	} else {
		fmt.Println("CHROM 1 is NOT completed")
	}
}
