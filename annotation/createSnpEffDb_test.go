package annotation

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestHeaderTranscriptID(t *testing.T) {
	txIDs := map[string]string{
		"Cpep01G0000100.1":        "Cpep01G0000100.1",
		"GWHTHOJH000001.1":        "Cpep01G0000100.1",
		"Phvul.001G000400.1.v2.1": "Phvul.001G000400.1.v2.1",
	}
	cases := []struct {
		name, header, want string
		ok                 bool
	}{
		{"already an ID", "Cpep01G0000100.1 some description", "", false},
		{"GWH CDS named by accession", "GWHTHOJH000001.1\tProtein=GWHPHOJH000001.1\tGene=GWHGHOJH000001.1\tPosition=GWHHOJH00000001.1: 14315-14928: -\tOriID=Cpep01G0000100.1", "Cpep01G0000100.1", true},
		{"GWH protein named by protein accession", "GWHPHOJH000001.1\tmRNA=GWHTHOJH000001.1\tGene=GWHGHOJH000001.1\tOriTrascriptID=Cpep01G0000100.1", "Cpep01G0000100.1", true},
		{"Phytozome ID= key", "Phvul.001G000400.1 pacid=37102 transcript=Phvul.001G000400.1 locus=Phvul.001G000400 ID=Phvul.001G000400.1.v2.1 annot-version=v2.1", "Phvul.001G000400.1.v2.1", true},
		{"nothing resolves", "unknown.1 locus=nope", "", false},
		{"empty header", "", "", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, ok := headerTranscriptID(c.header, txIDs)
			if got != c.want || ok != c.ok {
				t.Errorf("headerTranscriptID(%q) = %q, %v; want %q, %v", c.header, got, ok, c.want, c.ok)
			}
		})
	}
}

func TestTranscriptIDsMapsGWHAccession(t *testing.T) {
	dir := t.TempDir()
	gff := filepath.Join(dir, "genes.gff.gz")
	f, err := os.Create(gff)
	if err != nil {
		t.Fatal(err)
	}
	zw := gzip.NewWriter(f)
	zw.Write([]byte(strings.Join([]string{
		"##gff-version 3",
		"GWHHOJH00000001.1\tGWH\tgene\t14315\t16355\t.\t-\t.\tID=Cpep01G0000100;Accession=GWHGHOJH000001.1",
		"GWHHOJH00000001.1\tGWH\tmRNA\t14315\t16355\t.\t-\t.\tID=Cpep01G0000100.1;Accession=GWHTHOJH000001.1;Parent=Cpep01G0000100",
		"GWHHOJH00000001.1\tGWH\tCDS\t14315\t14928\t.\t-\t0\tID=Cpep01G0000100.1.cds1;Parent=Cpep01G0000100.1",
	}, "\n") + "\n"))
	zw.Close()
	f.Close()

	ids, err := transcriptIDs(gff)
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]string{
		"Cpep01G0000100.1": "Cpep01G0000100.1",
		"GWHTHOJH000001.1": "Cpep01G0000100.1",
	}
	if len(ids) != len(want) {
		t.Fatalf("got %v, want %v", ids, want)
	}
	for k, v := range want {
		if ids[k] != v {
			t.Errorf("ids[%q] = %q, want %q", k, ids[k], v)
		}
	}
}
