package genespace

import (
	"fmt"
	"regexp"
	"sort"
	"strings"
)

// Term is one gene-description pattern in a trait profile. Label is what the
// "Matched terms" column shows, so a user can see why a gene scored.
type Term struct {
	Label string
	Re    *regexp.Regexp
}

// Profile holds the trait-specific knowledge GeneSpace uses to judge how
// plausible a gene is for the trait. Everything else in the scoring (genotype
// concordance, variant impact) is trait-independent.
//
//	Strong – gene families with a direct, well-established role in the trait.
//	Weak   – families that are often involved but are large or indirect.
//	TF     – transcription-factor families that commonly regulate the trait.
type Profile struct {
	Name        string
	Description string
	Strong      []Term
	Weak        []Term
	TF          []Term
}

// term compiles a case-insensitive pattern. Patterns use \b word boundaries so
// short family names don't match inside unrelated words ("ERF" in
// "superfamily", "RGA" in "organ", "NAC" in "spinach").
func term(label, pattern string) Term {
	return Term{Label: label, Re: regexp.MustCompile(`(?i)` + pattern)}
}

// DefaultProfile is the trait used when none is given.
const DefaultProfile = "disease-resistance"

var profiles = map[string]Profile{
	"disease-resistance": {
		Name:        "disease-resistance",
		Description: "Disease resistance: NLRs, PRR receptor kinases and susceptibility (S) genes",
		Strong: []Term{
			// NLRs (intracellular immune receptors)
			term("NLR", `\bNBS-LRR\b|\bNB-ARC\b|\bNLRs?\b|\bNBS\b|\bTIR-NBS|\bCC-NBS|\bRPW8\b`),
			term("disease resistance protein", `disease resistance (protein|gene)`),
			term("RGA", `\bRGA\d*\b`),
			term("RPM1", `\bRPM1\b`),
			term("ZAR1", `\bZAR1\b`),
			// Pattern-recognition receptors
			term("LRR-RLK", `\bLRR-RLKs?\b|leucine-rich repeat receptor-like (serine/threonine-)?(protein )?kinase|\bLRR receptor-like`),
			term("RLP", `\bRLP\d*\b|leucine-rich repeat receptor-like protein|receptor-like protein \d+`),
			term("LysM receptor", `lysm (domain[- ])?(containing )?receptor-like kinase|\bLYK\d+\b|\bCERK1\b|\bLYM\d+\b`),
			// Susceptibility genes: loss of function gives (usually recessive) resistance
			term("S-gene:MLO", `\bMLO\d*\b`),
			term("S-gene:eIF4E", `\beIF\(?iso\)?4E\b|initiation factor (\(?iso\)?)?4E\b`),
			term("S-gene:DMR6", `\bDMR6\b|downy mildew resistance 6`),
			term("S-gene:SWEET", `\bSWEET\d*\b`),
		},
		Weak: []Term{
			term("PR protein", `pathogenesis-related|\bPR-?\d+\b`),
			term("thaumatin", `thaumatin`),
			term("chitinase", `chitinase`),
			term("glucanase", `glucanase|glucan endo-1,3-beta|\(1->3\)-beta-glucan endohydrolase`),
			term("peroxidase", `\bperoxidase`),
			term("lipoxygenase", `lipoxygenase`),
			term("defensin", `defensin`),
			term("RBOH", `respiratory burst oxidase`),
			term("CDR1", `\bCDR1\b`),
			term("subtilisin", `subtilisin`),
			term("PGIP", `polygalacturonase inhibitor|\bPGIP\d*\b`),
			term("dirigent", `dirigent|disease resistance response`),
			term("immune signalling", `\bEDS1\b|\bPAD4\b|\bNDR1\b|\bNPR\d\b`),
			term("receptor-like kinase", `receptor-like (serine/threonine-)?(protein )?kinase|wall-associated receptor kinase`),
		},
		TF: []Term{
			term("WRKY", `\bWRKY\d*\b`),
			term("MYB", `\bMYB`),
			term("NAC", `\bNAC\d*\b|\bNAC domain`),
			term("bHLH", `\bbHLH\d*\b`),
			term("bZIP", `\bbZIP\d*\b|basic leucine zipper`),
			term("AP2/ERF", `\bAP2\b|\bERF\d*\b|ethylene-responsive transcription factor`),
			term("transcription factor", `transcription factor`),
		},
	},
}

// GetProfile returns the built-in profile with the given name.
func GetProfile(name string) (Profile, error) {
	if p, ok := profiles[name]; ok {
		return p, nil
	}
	return Profile{}, fmt.Errorf("unknown trait profile %q (available: %s)", name, strings.Join(ProfileNames(), ", "))
}

// ProfileNames lists the built-in profiles, sorted.
func ProfileNames() []string {
	names := make([]string, 0, len(profiles))
	for n := range profiles {
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// matchTerms returns the labels of every term whose pattern occurs in desc.
func matchTerms(terms []Term, desc string) []string {
	var hits []string
	for _, t := range terms {
		if t.Re.MatchString(desc) {
			hits = append(hits, t.Label)
		}
	}
	return hits
}
