package warehouse

import (
	"encoding/json"
	"io"
	"sync"
	"time"
)

// ScanEvent is the live, append-only event contract consumed by data_warehouse.
// The final report remains the authoritative atomic inventory snapshot.
type ScanEvent struct {
	Kind      string `json:"kind"`
	Sequence  int64  `json:"sequence"`
	EmittedAt string `json:"emitted_at"`
	Event     string `json:"event"`
	Phase     string `json:"phase,omitempty"`
	Message   string `json:"message,omitempty"`
	Volume    string `json:"volume,omitempty"`
	Path      string `json:"path,omitempty"`
	Completed int    `json:"completed,omitempty"`
	Total     int    `json:"total,omitempty"`
	Payload   any    `json:"payload,omitempty"`
}

// EventWriter serializes events from concurrent inspection workers.
type EventWriter struct {
	mu  sync.Mutex
	enc *json.Encoder
	seq int64
}

func NewEventWriter(out io.Writer) *EventWriter {
	if out == nil {
		return nil
	}
	return &EventWriter{enc: json.NewEncoder(out)}
}

func (w *EventWriter) Emit(event ScanEvent) {
	if w == nil {
		return
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	w.seq++
	event.Kind = "scan_event"
	event.Sequence = w.seq
	event.EmittedAt = time.Now().UTC().Format(time.RFC3339Nano)
	_ = w.enc.Encode(event)
}

func (w *EventWriter) Heartbeats(stop <-chan struct{}) {
	if w == nil {
		return
	}
	ticker := time.NewTicker(2 * time.Second)
	defer ticker.Stop()
	for {
		select {
		case <-stop:
			return
		case <-ticker.C:
			w.Emit(ScanEvent{Event: "heartbeat"})
		}
	}
}
