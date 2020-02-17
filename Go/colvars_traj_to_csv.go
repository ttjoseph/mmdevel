// Converts colvars variables trajectory (*.colvars.traj) from its weird fixed space
// format where headers are repeated every so often into CSV so as to make further processing
// and analysis easier.
//
// Reads from stdin and prints to stdout.
//
// Tom Joseph - thomas dot joseph at pennmedicine dot upenn dot edu
package main

import (
	"os"
	"fmt"
	"bufio"
	"strings"
	"encoding/csv"
)


func main() {
	scanner := bufio.NewScanner(os.Stdin)
	writer := csv.NewWriter(os.Stdout)
	headerWasPrinted := false
	headerString := ""
	for scanner.Scan() {
		// Trim leading '#' and leading/trailing whitespace.
		// The header has a '#' to start for some reason. If we get rid of that
		// we can treat it like any other line, unless we've already printed headers.
		s := strings.TrimSpace(scanner.Text())
		if s[:1] == "#" {
			s = strings.TrimSpace(strings.TrimPrefix(s, "#"))
			if headerWasPrinted == false {
				headerString = s
				headerWasPrinted = true
			} else {
				if headerString != s {
					fmt.Fprintln(os.Stderr, "\n=========== First header string ===========")
					fmt.Fprintln(os.Stderr, headerString)
					fmt.Fprintln(os.Stderr, "\n----------- Latest header string -----------")
					fmt.Fprintln(os.Stderr, s)
					panic("Encountered a header string that doesn't match the first one. Did you mix output from different colvars configs?")
				}
				continue
			}
		}

		if err := scanner.Err(); err != nil {
			fmt.Fprintln(os.Stderr, "reading standard input: ", err)
		}

		// Tokenize text into record
		record := strings.Fields(s)
		// Write the record in CSV format
		writer.Write(record)
	}
	writer.Flush()
}