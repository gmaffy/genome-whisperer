package utils

import (
	"bytes"
	"fmt"
	"os"
	"strings"
	"sync"
	"sync/atomic"
	"testing"

	"github.com/fatih/color"
	"github.com/schollz/progressbar/v3"
)

// attachTestBar attaches a bar writing into buf, bypassing the terminal check
// so the test does not need a TTY. It returns the release func.
//
// AttachBar itself refuses to attach without a terminal, which is the behaviour
// under test elsewhere; here the goal is the routing, so the state is set up
// the way AttachBar would.
func attachTestBar(t *testing.T, buf *bytes.Buffer, max int64) func() {
	t.Helper()

	bar := progressbar.NewOptions64(max,
		progressbar.OptionSetWriter(buf),
		progressbar.OptionThrottle(0),
		progressbar.OptionShowCount(),
		progressbar.OptionShowIts(),
		progressbar.OptionSetWidth(10),
	)

	attachMu.Lock()
	if activeBar.Load() != nil {
		attachMu.Unlock()
		t.Fatal("a bar is already attached")
	}
	savedColorOutput, savedColorError = color.Output, color.Error
	sink := barSink{}
	color.Output, color.Error = sink, sink
	activeBar.Store(bar)
	attachMu.Unlock()

	return func() {
		attachMu.Lock()
		defer attachMu.Unlock()
		activeBar.Store(nil)
		color.Output, color.Error = savedColorOutput, savedColorError
	}
}

// The bar renders with carriage returns and no newlines, so stripping every \r
// segment that is not a line of ours leaves the log behind. Log lines are the
// ones the test wrote, so match them by prefix instead of trying to un-render
// the bar.
func loggedLines(buf *bytes.Buffer, prefix string) []string {
	var out []string
	for _, chunk := range strings.Split(buf.String(), "\n") {
		for _, part := range strings.Split(chunk, "\r") {
			if strings.HasPrefix(part, prefix) {
				out = append(out, part)
			}
		}
	}
	return out
}

func TestPrintfGoesThroughBarUnsplit(t *testing.T) {
	var buf bytes.Buffer
	release := attachTestBar(t, &buf, 200)
	defer release()

	const workers, each = 25, 8

	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			for i := 0; i < each; i++ {
				Printf("LOG worker=%d line=%d\n", w, i)
			}
		}(w)
	}
	wg.Wait()

	got := loggedLines(&buf, "LOG ")
	if len(got) != workers*each {
		t.Fatalf("got %d log lines, want %d", len(got), workers*each)
	}

	seen := make(map[string]int, len(got))
	for _, line := range got {
		seen[line]++
	}
	for w := 0; w < workers; w++ {
		for i := 0; i < each; i++ {
			want := fmt.Sprintf("LOG worker=%d line=%d", w, i)
			if seen[want] != 1 {
				t.Errorf("line %q appeared %d times, want 1", want, seen[want])
			}
		}
	}
}

func TestLogWriterEmitsWholeLinesOnly(t *testing.T) {
	var buf bytes.Buffer
	release := attachTestBar(t, &buf, 100)
	defer release()

	// Each child gets its own writer, which is what RunBashCmdVerbose does, so
	// a line arriving in fragments cannot be interleaved with another's.
	const writers = 10

	var wg sync.WaitGroup
	for n := 0; n < writers; n++ {
		wg.Add(1)
		go func(n int) {
			defer wg.Done()
			w := LogWriter()
			defer w.Close()
			line := fmt.Sprintf("LOG child=%d payload=abcdefghij\n", n)
			// One byte at a time: only the newline may release it.
			for i := 0; i < len(line); i++ {
				if _, err := w.Write([]byte{line[i]}); err != nil {
					t.Error(err)
					return
				}
			}
		}(n)
	}
	wg.Wait()

	got := loggedLines(&buf, "LOG ")
	if len(got) != writers {
		t.Fatalf("got %d lines %q, want %d whole lines", len(got), got, writers)
	}
	for _, line := range got {
		if !strings.HasSuffix(line, "payload=abcdefghij") {
			t.Errorf("line was split or interleaved: %q", line)
		}
	}
}

func TestLogWriterTerminators(t *testing.T) {
	tests := []struct {
		name  string
		write string
		want  string
	}{
		// A carriage return overwrites the pending line rather than emitting
		// it, so a tool redrawing one line settles on its final text instead
		// of filling the scrollback with a frame per update. Close emits what
		// it settled on.
		{"carriage return overwrites", "LOG stale 10%\rLOG progress 90%", "LOG progress 90%"},
		{"newline", "LOG done\n", "LOG done"},
		// A final line with no terminator at all is flushed by Close.
		{"unterminated flushed on close", "LOG truncated", "LOG truncated"},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			release := attachTestBar(t, &buf, 10)
			defer release()

			w := LogWriter()
			if _, err := w.Write([]byte(tc.write)); err != nil {
				t.Fatal(err)
			}
			if err := w.Close(); err != nil {
				t.Fatal(err)
			}

			got := loggedLines(&buf, "LOG ")
			if len(got) != 1 || got[0] != tc.want {
				t.Fatalf("got %q, want exactly [%q]", got, tc.want)
			}
		})
	}
}

func TestColorIsRoutedThroughBarAndRestored(t *testing.T) {
	before, beforeErr := color.Output, color.Error

	var buf bytes.Buffer
	release := attachTestBar(t, &buf, 10)

	if color.Output == before {
		t.Error("color.Output was not redirected while a bar is attached")
	}
	if color.Error == beforeErr {
		t.Error("color.Error was not redirected while a bar is attached")
	}

	// NoColor is set from stdout being a terminal, which it is not under `go
	// test`, so the text arrives unescaped and is straightforward to assert.
	color.Red("LOG via color\n")

	release()

	if color.Output != before || color.Error != beforeErr {
		t.Error("color writers were not restored on release")
	}
	if got := loggedLines(&buf, "LOG "); len(got) != 1 {
		t.Fatalf("got %q, want the color.Red line routed through the bar", got)
	}
}

func TestAttachIsNotReentrant(t *testing.T) {
	var outer bytes.Buffer
	release := attachTestBar(t, &outer, 10)
	defer release()

	// The real entry point must refuse rather than swap the owner out.
	inner := progressbar.NewOptions64(10, progressbar.OptionSetWriter(&bytes.Buffer{}))
	innerRelease := AttachBar(inner)
	if activeBar.Load() == inner {
		t.Fatal("an inner bar took the terminal from the attached one")
	}
	innerRelease()

	if activeBar.Load() == nil {
		t.Fatal("releasing the inner no-op detached the outer bar")
	}
}

func TestWriteFallsBackToStdoutWithNoBar(t *testing.T) {
	if activeBar.Load() != nil {
		t.Fatal("a bar is attached; tests are not isolated")
	}

	// LogWriter must hand back os.Stdout itself, so that assigning it to an
	// exec.Cmd passes the descriptor straight through as it did before.
	w := LogWriter()
	if got, ok := w.(nopCloser); !ok || got.Writer != os.Stdout {
		t.Fatalf("LogWriter() = %T, want nopCloser wrapping os.Stdout", w)
	}

	// And a write must not panic or land anywhere unexpected.
	var reached atomic.Bool
	func() {
		defer func() {
			if r := recover(); r != nil {
				t.Errorf("write panicked with no bar attached: %v", r)
			}
		}()
		Printf("")
		reached.Store(true)
	}()
	if !reached.Load() {
		t.Error("Printf did not complete with no bar attached")
	}
}

// fatih/color emits the SGR prefix, the text and the reset as three separate
// writes, so a sink that forwards each write on its own hands the bar a lone
// "\x1b[31m" and the bar repaints itself red for every remaining line.
func TestColorEscapesDoNotLeakOntoTheBar(t *testing.T) {
	var buf bytes.Buffer
	release := attachTestBar(t, &buf, 10)
	defer release()

	// color decides on escapes from stdout, which is not a terminal under `go
	// test`; force them on so this exercises the real terminal case.
	noColor := color.NoColor
	color.NoColor = false
	defer func() { color.NoColor = noColor }()

	color.Red("LOG failure\n")
	_ = activeBar.Load().Add(1)

	out := buf.String()
	for _, chunk := range strings.Split(out, "\n") {
		for _, part := range strings.Split(chunk, "\r") {
			opened := strings.Count(part, "\x1b[31m")
			if opened == 0 {
				continue
			}
			if !strings.Contains(part, "\x1b[0m") {
				t.Errorf("a line opened a colour it never reset, so the bar inherits it: %q", part)
			}
		}
	}
	// The bar's own repaint must never carry a colour it did not set.
	if i := strings.LastIndex(out, "\x1b[31m"); i >= 0 {
		rest := out[i:]
		if !strings.Contains(rest[:min(len(rest), 64)], "\x1b[0m") {
			t.Errorf("colour left open before the trailing bar render: %q", rest)
		}
	}
}

// A run that finds nothing builds a bar over zero items. progressbar rejects
// every Add on a zero-max bar ("max must be greater than 0"), so the nudge that
// flushes buffered lines fails and anything printed would be swallowed.
func TestPrintfSurvivesAnEmptyBar(t *testing.T) {
	if activeBar.Load() != nil {
		t.Fatal("a bar is attached; tests are not isolated")
	}

	// NewBar makes a zero-length bar invisible, and AttachBar refuses one, so
	// the message has to fall through to stdout rather than into a bar that
	// cannot tick.
	bar := NewBar(0, "nothing to do")
	release := AttachBar(bar)
	defer release()

	if activeBar.Load() != nil {
		t.Error("a bar that cannot tick was attached")
	}

	stdout := captureStdout(t, func() {
		Printf("LOG still visible\n")
	})
	if !strings.Contains(stdout, "LOG still visible") {
		t.Errorf("message was swallowed by a zero-length bar; stdout = %q", stdout)
	}

	// And the same when the guard is bypassed the way attachTestBar does, which
	// is the shape a future caller could reach by building a bar by hand.
	var buf bytes.Buffer
	zero := progressbar.NewOptions64(0, progressbar.OptionSetWriter(&buf), progressbar.OptionThrottle(0))
	attachMu.Lock()
	activeBar.Store(zero)
	attachMu.Unlock()
	defer func() {
		attachMu.Lock()
		activeBar.Store(nil)
		attachMu.Unlock()
	}()

	// write falls back to a direct print, because a zero-max bar rejects the
	// Add that would flush it and so can never render what it was handed.
	stdout = captureStdout(t, func() {
		Printf("LOG through a zero-max bar\n")
	})
	if !strings.Contains(stdout, "LOG through a zero-max bar") {
		t.Errorf("line lost to a zero-max bar; stdout = %q, bar buffer = %q", stdout, buf.String())
	}
}

// A failing external tool reports why in its last lines. Losing them was the
// difference between "gatk MergeVcfs: exit status 3" and a message naming a
// read-only output volume, so the tail has to survive a chatty tool that wrote
// thousands of progress lines before it died.
func TestRunCmdScanReportsTailOnFailure(t *testing.T) {
	script := `for i in $(seq 1 5000); do echo "Processed $i records"; done; echo "" >&2; echo "Caused by: Read-only file system" >&2; exit 3`

	stdout := captureStdout(t, func() {
		if err := RunCmdScan(script, false, func(string) {}); err == nil {
			t.Fatal("expected an error from a command that exits 3")
		}
	})

	if !strings.Contains(stdout, "Caused by: Read-only file system") {
		t.Errorf("failure report dropped the reason:\n%s", stdout)
	}
	if !strings.Contains(stdout, "command: "+script) {
		t.Errorf("failure report did not name the command:\n%s", stdout)
	}
	// Blank lines are dropped, so a tool that pads its output cannot push the
	// reason out of the tail.
	if strings.Contains(stdout, "Processed 4959 records") {
		t.Errorf("tail kept more than %d lines:\n%s", failTail, stdout)
	}
}

// Under --verbose every line has already been printed, so repeating the tail
// would duplicate what is on screen.
func TestRunCmdScanVerboseReportsCommandOnly(t *testing.T) {
	stdout := captureStdout(t, func() {
		if err := RunCmdScan(`printf 'bo%s\\n' om >&2; exit 1`, true, func(string) {}); err == nil {
			t.Fatal("expected an error")
		}
	})

	if strings.Count(stdout, "boom") != 1 {
		t.Errorf("verbose failure repeated the output:\n%s", stdout)
	}
	if !strings.Contains(stdout, "command: ") {
		t.Errorf("verbose failure did not name the command:\n%s", stdout)
	}
}

// captureStdout redirects os.Stdout for the duration of fn.
func captureStdout(t *testing.T, fn func()) string {
	t.Helper()

	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	saved := os.Stdout
	os.Stdout = w

	done := make(chan string, 1)
	go func() {
		var b bytes.Buffer
		_, _ = b.ReadFrom(r)
		done <- b.String()
	}()

	fn()

	os.Stdout = saved
	_ = w.Close()
	out := <-done
	_ = r.Close()
	return out
}
