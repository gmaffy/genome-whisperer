package utils

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/fatih/color"
	"github.com/mattn/go-isatty"
	"github.com/schollz/progressbar/v3"
)

// A progress bar owns one line of the terminal and repaints that line with a
// carriage return, so anything else that reaches the terminal mid-repaint
// strands the bar's own characters there — the fragmented bar seen under
// --verbose, where every worker prints a command banner while its samtools or
// GATK child writes straight to the stdout and stderr it inherited.
//
// While a bar is attached here it is the only writer to the terminal.
// Everything else is handed to the bar, which buffers it and flushes it above
// itself on its next repaint, so the bar stays the bottom line and the log
// scrolls past it.
//
// Two properties keep this to a handful of edits rather than a repo-wide
// refactor: fatih/color's Output and Error are package variables every
// color.Red/Green/Cyan call writes through, so attaching redirects all of them
// at once; and RunBashCmdVerbose is the single place a verbose child process is
// spawned, so pointing its stdio at LogWriter covers every caller.
var (
	// attachMu serializes attach against release. The bar itself is held in an
	// atomic pointer so that write, which runs on every log line from every
	// worker, never contends for it.
	attachMu  sync.Mutex
	activeBar atomic.Pointer[progressbar.ProgressBar]

	// lineMu is held across one whole logical line, so a line buffered by one
	// worker cannot be split by another's.
	lineMu sync.Mutex

	savedColorOutput io.Writer
	savedColorError  io.Writer
)

// RunCmdScan runs cmdStr and hands every line it writes to scan, forwarding
// those lines to the console only when verbose is set.
//
// It exists for the steps whose progress only the external tool knows. gatk
// MergeVcfs spends minutes inside java and reports how far it has got in its own
// log, so that log is the only thing a progress bar can be driven from — and it
// has to be read whether or not the user asked to see it. Without this, the
// choice would be a bar or --verbose, not both.
//
// scan runs on the goroutine os/exec uses to drain the pipe, so it must not
// block; feeding a progress bar is exactly the right amount of work for it.
func RunCmdScan(cmdStr string, verbose bool, scan func(line string)) error {
	cmd := exec.Command("bash", "-c", cmdStr)

	var tail tailBuf
	out := &lineWriter{emit: func(line string) {
		scan(line)
		tail.add(line)
		if verbose {
			write(line + "\n")
		}
	}}
	// Same writer for both, so os/exec gives the child one pipe and keeps the
	// two streams in the order the tool wrote them.
	cmd.Stdout = out
	cmd.Stderr = out

	err := cmd.Run()
	_ = out.Close()
	if err != nil {
		Printf("CMD error: %v\n%s\n", err, tail.report(cmdStr, verbose))
		return err
	}
	return nil
}

// failTail is how many of a failed command's last lines are reported with the
// error. A tool that dies says why in its last few lines — for a java tool that
// is the exception message plus its stack trace — so the tail is what makes an
// "exit status 3" diagnosable. GATK's is around twenty lines, so this leaves
// room for the "Caused by" underneath it.
const failTail = 40

// tailBuf keeps the last failTail non-blank lines a command wrote.
//
// It is filled from the goroutine os/exec drains the pipe on and read only
// after cmd.Run has returned, which has already waited for that goroutine, so
// it needs no lock of its own.
type tailBuf struct{ lines []string }

func (t *tailBuf) add(line string) {
	if strings.TrimSpace(line) == "" {
		return
	}
	t.lines = append(t.lines, line)
	if len(t.lines) > failTail {
		t.lines = t.lines[len(t.lines)-failTail:]
	}
}

// report renders the failing command and its last lines. Under --verbose those
// lines have already scrolled past, so only the command is repeated — the tail
// would be an immediate duplicate of what is on screen.
func (t *tailBuf) report(cmdStr string, verbose bool) string {
	var b strings.Builder
	fmt.Fprintf(&b, "  command: %s\n", cmdStr)
	if verbose || len(t.lines) == 0 {
		return b.String()
	}
	fmt.Fprintf(&b, "  last %d line(s) of output:\n", len(t.lines))
	for _, line := range t.lines {
		fmt.Fprintf(&b, "    %s\n", line)
	}
	return b.String()
}

// StderrIsTerminal reports whether the progress bars have a terminal to draw
// on. When they do not — a redirected run, a CI log — a bar is not worth
// drawing and its ANSI repaint frames are noise in the file.
func StderrIsTerminal() bool {
	fd := os.Stderr.Fd()
	return isatty.IsTerminal(fd) || isatty.IsCygwinTerminal(fd)
}

// NewBar is progressbar.Default with three pipeline defaults added:
//
//   - it is unthrottled, so a log line handed to it never waits on the next
//     tick to become visible — a render inside the throttle window returns
//     before it copies the buffered lines out, and a bar of known length has no
//     timer of its own, so a throttled bar could hold the last line back for as
//     long as a HaplotypeCaller call takes;
//   - it is invisible when stderr is not a terminal, because a carriage-return
//     repaint is line noise in a redirected log;
//   - it is invisible for max 0, the shape a run that found nothing produces.
//     progressbar rejects every Add on such a bar, so a visible one would paint
//     a blank bar it could never advance or terminate. A max of -1, meaning
//     indeterminate, is a real bar and stays visible.
//
// Extra options are appended, so a caller can override any of them.
func NewBar(max int64, desc string, extra ...progressbar.Option) *progressbar.ProgressBar {
	opts := []progressbar.Option{
		progressbar.OptionSetDescription(desc),
		progressbar.OptionSetWriter(os.Stderr),
		progressbar.OptionSetWidth(10),
		progressbar.OptionShowTotalBytes(true),
		progressbar.OptionThrottle(0),
		progressbar.OptionShowCount(),
		progressbar.OptionShowIts(),
		progressbar.OptionOnCompletion(func() { fmt.Fprint(os.Stderr, "\n") }),
		progressbar.OptionSpinnerType(14),
		progressbar.OptionFullWidth(),
		progressbar.OptionSetRenderBlankState(true),
		progressbar.OptionSetVisibility(StderrIsTerminal() && max != 0),
	}
	return progressbar.NewOptions64(max, append(opts, extra...)...)
}

// AttachBar makes bar the terminal's owner until the returned func is called.
//
// It is a no-op in two cases, both of which return a release that does nothing.
// When stderr is not a terminal there is no bar on screen to protect, and
// routing through an invisible bar would quietly move messages from stdout onto
// stderr — an invisible bar reports its output as uncacheable, so Bprintf goes
// straight to the bar's own writer. And when another bar is already attached the
// first owner is left in place, so a bar opened inside another's scope cannot
// take the terminal away from it.
func AttachBar(bar *progressbar.ProgressBar) (release func()) {
	// A bar over an empty work list has max 0, and progressbar rejects every
	// Add on one ("max must be greater than 0"). Since the nudge that flushes
	// buffered lines is an Add, attaching such a bar would swallow every
	// message printed under it. A run that found nothing has no bar worth
	// protecting anyway. An indeterminate bar, built with max -1, reports its
	// width here instead and is attachable.
	if bar == nil || !StderrIsTerminal() || bar.GetMax64() <= 0 {
		return func() {}
	}

	attachMu.Lock()
	defer attachMu.Unlock()

	if activeBar.Load() != nil {
		return func() {}
	}

	savedColorOutput, savedColorError = color.Output, color.Error
	sink := barSink{}
	color.Output, color.Error = sink, sink
	activeBar.Store(bar)

	return func() {
		attachMu.Lock()
		defer attachMu.Unlock()

		if activeBar.Load() != bar {
			return
		}
		activeBar.Store(nil)
		color.Output, color.Error = savedColorOutput, savedColorError
	}
}

// Printf writes one message without breaking the active bar.
func Printf(format string, a ...interface{}) {
	write(fmt.Sprintf(format, a...))
}

// Println writes one message, newline-terminated, without breaking the active
// bar.
func Println(a ...interface{}) {
	write(fmt.Sprintln(a...))
}

// write is the single path everything in this file funnels into. With a bar
// attached the text is buffered by the bar and the Add(0) that follows forces
// the repaint that flushes it: Add reaches render for these bars because
// ShowCount and ShowIts are set, and render only bails early on the throttle,
// which NewBar sets to zero. So the buffer never holds more than one line, even
// when a chatty child is streaming through it.
func write(s string) {
	if s == "" {
		return
	}
	if bar := activeBar.Load(); bar != nil {
		lineMu.Lock()
		defer lineMu.Unlock()
		_, _ = progressbar.Bprintf(bar, "%s", s)
		if err := bar.Add(0); err == nil {
			return
		}
		// The bar cannot tick, so it will never flush what was just handed to
		// it. AttachBar rejects the one bar that behaves this way, so reaching
		// here means a bar was attached by some other route; print directly
		// rather than lose the line. The buffered copy stays stranded in a bar
		// that can never render, so this cannot double up.
	}
	fmt.Print(s)
}

// barSink adapts write to io.Writer so fatih/color can be pointed at it.
//
// It has to buffer, because color does not write a line at a time. color.Red
// makes three writes to color.Output — the SGR prefix, the text, then the reset
// — so forwarding each one as it arrives hands the bar a lone "\x1b[31m" and
// every repaint after it comes out red until the reset finally lands. Holding
// the bytes until the newline keeps prefix, text and reset in one line, and any
// line that still opened a colour it did not close gets a reset appended, so a
// colour can never escape onto the bar's own line.
type barSink struct{}

var (
	sinkMu  sync.Mutex
	sinkBuf bytes.Buffer
)

func (barSink) Write(p []byte) (int, error) {
	sinkMu.Lock()
	pending := func() []string {
		var lines []string
		sinkBuf.Write(p)
		for {
			line, err := sinkBuf.ReadString('\n')
			if err != nil {
				// Not a whole line yet. An escape-only tail is the reset color
				// writes after the newline; there is nothing more coming for
				// it, so drop it rather than prepend it to the next message.
				if line != "" && !hasPrintable(line) {
					line = ""
				}
				sinkBuf.Reset()
				sinkBuf.WriteString(line)
				return lines
			}
			lines = append(lines, closeColor(line))
		}
	}()
	sinkMu.Unlock()

	for _, line := range pending {
		write(line)
	}
	return len(p), nil
}

// closeColor appends a reset to a line that opened an SGR sequence and never
// closed it, so the colour cannot bleed into whatever is drawn next.
func closeColor(line string) string {
	if !strings.Contains(line, "\x1b[") {
		return line
	}
	body := strings.TrimSuffix(line, "\n")
	if strings.HasSuffix(body, colorReset) {
		return line
	}
	return body + colorReset + "\n"
}

const colorReset = "\x1b[0m"

// hasPrintable reports whether s contains anything other than ANSI escape
// sequences and whitespace.
func hasPrintable(s string) bool {
	for {
		i := strings.Index(s, "\x1b[")
		if i < 0 {
			break
		}
		end := i + 2
		for end < len(s) && !(s[end] >= 0x40 && s[end] <= 0x7e) {
			end++
		}
		if end < len(s) {
			end++
		}
		s = s[:i] + s[end:]
	}
	return strings.TrimSpace(s) != ""
}

// LogWriter returns the sink an external tool's streamed output should be
// written to.
//
// Each call returns a fresh writer with its own partial-line buffer, so one
// child's half-written line cannot land in the middle of another's; only whole
// lines take the shared lock. Close it once the child has exited, to flush a
// final line that arrived without a terminator.
//
// With no bar attached it hands back os.Stdout unchanged. Assigning that to
// both Stdout and Stderr of an exec.Cmd passes the file descriptor straight
// through, so streaming outside a bar stays byte for byte what it was before
// this file existed — including a tool that draws its own progress with
// carriage returns and no newline at all.
func LogWriter() io.WriteCloser {
	if activeBar.Load() == nil {
		return nopCloser{os.Stdout}
	}
	return &lineWriter{emit: func(line string) { write(line + "\n") }}
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

// lineWriter accumulates a child process's bytes and emits whole lines.
//
// Only a newline emits. A carriage return discards what precedes it on the
// pending line instead, which is what the terminal it replaces would have done
// with it: a tool that redraws one line — the docker pull in CreateGvcfDV,
// bcftools' record counters — settles on its final text rather than filling the
// scrollback with one line per frame, each of which would repaint the bar.
//
// Holding a partial line back is a correctness requirement, not tidiness.
// os/exec delivers whatever chunk it read, not lines; a chunk ending mid-line
// makes the bar repaint on the same line as that text, and the bar's next
// repaint clears its line with spaces and erases the text along with it.
type lineWriter struct {
	buf bytes.Buffer

	// emit receives one complete line, without its terminator.
	emit func(line string)
}

// lineCap bounds the pending line, so a tool that never emits a newline at all
// cannot grow this without limit.
const lineCap = 64 << 10

func (w *lineWriter) Write(p []byte) (int, error) {
	w.buf.Write(p)

	for {
		s := w.buf.String()
		i := strings.IndexAny(s, "\n\r")
		if i < 0 {
			break
		}
		if s[i] == '\r' {
			// Overwrite: drop the pending text and the \r itself.
			w.buf.Next(i + 1)
			continue
		}
		line := s[:i]
		w.buf.Next(i + 1)
		w.emit(line)
	}

	if w.buf.Len() > lineCap {
		w.flush()
	}
	return len(p), nil
}

// Close emits a last line the tool left without a newline. cmd.Run has already
// waited for os/exec's copying goroutine by the time it returns, so closing
// after Run cannot race a write.
func (w *lineWriter) Close() error {
	w.flush()
	return nil
}

func (w *lineWriter) flush() {
	if w.buf.Len() > 0 {
		line := w.buf.String()
		w.buf.Reset()
		w.emit(line)
	}
}
