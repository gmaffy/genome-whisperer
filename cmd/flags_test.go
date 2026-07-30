package cmd

// Flag consistency tests.
//
// These work purely by inspecting the commands cobra has actually registered,
// so they need no helper machinery in the package itself: flags are declared
// with plain pflag calls (shared ones on rootCmd.PersistentFlags(),
// command-local ones in each command's init) and these tests hold the result
// to a consistent shape.

import (
	"go/ast"
	"go/parser"
	"go/token"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// kebabCase is the naming convention every long flag name must follow.
var kebabCase = regexp.MustCompile(`^[a-z0-9]+(-[a-z0-9]+)*$`)

// flagsLocalByDesign are shared by name but cannot be persistent flags, so each
// command declares its own. The reason is recorded here so the "same flag,
// same shape" test can allow the specific difference and nothing else.
var flagsLocalByDesign = map[string]string{
	"variant": "a list in the annotation commands, a single file in the filter commands",
	"bqsr":    "defaults false for AlignReads/GoBSAseq, true for VariantCalling",
}

func allCommands() []*cobra.Command { return rootCmd.Commands() }

// eachLocalFlag visits every command-local flag (excludes inherited persistent
// flags, which are visited once via rootCmd.PersistentFlags()).
func eachLocalFlag(fn func(cmdName string, f *pflag.Flag)) {
	for _, c := range allCommands() {
		name := c.Name()
		c.LocalNonPersistentFlags().VisitAll(func(f *pflag.Flag) { fn(name, f) })
	}
}

// ---------------------------------------------------------------------------
// Naming convention
// ---------------------------------------------------------------------------

func TestPersistentFlagNamesAreKebabCase(t *testing.T) {
	rootCmd.PersistentFlags().VisitAll(func(f *pflag.Flag) {
		if !kebabCase.MatchString(f.Name) {
			t.Errorf("persistent flag %q is not lowercase kebab-case", f.Name)
		}
	})
}

func TestLocalFlagNamesAreKebabCase(t *testing.T) {
	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Name == "help" {
			return
		}
		if !kebabCase.MatchString(f.Name) {
			t.Errorf("%s: flag %q is not lowercase kebab-case", cmdName, f.Name)
		}
	})
}

// ---------------------------------------------------------------------------
// Shorthands
// ---------------------------------------------------------------------------

// TestShorthandsAreGloballyUnique makes sure a shorthand letter never means two
// different things across the CLI. Before the flag cleanup, -v meant verbose,
// version and vcf-table depending on the subcommand, and -g meant four
// different flags.
func TestShorthandsAreGloballyUnique(t *testing.T) {
	owner := make(map[string]string) // shorthand -> long name
	where := make(map[string]string) // shorthand -> first command seen

	rootCmd.PersistentFlags().VisitAll(func(f *pflag.Flag) {
		if f.Shorthand != "" {
			owner[f.Shorthand] = f.Name
			where[f.Shorthand] = "<root persistent>"
		}
	})

	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Shorthand == "" || f.Name == "help" {
			return
		}
		if prev, seen := owner[f.Shorthand]; seen && prev != f.Name {
			t.Errorf("shorthand -%s means %q in %s but %q in %s",
				f.Shorthand, prev, where[f.Shorthand], f.Name, cmdName)
			return
		}
		owner[f.Shorthand] = f.Name
		where[f.Shorthand] = cmdName
	})
}

// TestLocalFlagsDoNotShadowPersistentOnes catches a command redeclaring a flag
// that already exists on root, which would silently override the shared
// definition for that one command.
func TestLocalFlagsDoNotShadowPersistentOnes(t *testing.T) {
	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Name == "help" {
			return
		}
		if rootCmd.PersistentFlags().Lookup(f.Name) == nil {
			return
		}
		if _, ok := flagsLocalByDesign[f.Name]; ok {
			t.Errorf("%s declares %q locally AND it exists on root; pick one", cmdName, f.Name)
			return
		}
		t.Errorf("%s: local flag %q shadows the persistent flag of the same name", cmdName, f.Name)
	})
}

// ---------------------------------------------------------------------------
// Shared-but-local flags must still agree with each other
// ---------------------------------------------------------------------------

// TestLocalByDesignFlagsAgree checks that --variant and --bqsr, which are
// declared per command rather than on root, still use the same shorthand and
// description everywhere. Only the difference each is exempted for (type for
// --variant, default for --bqsr) is allowed to vary.
func TestLocalByDesignFlagsAgree(t *testing.T) {
	type occurrence struct{ cmd, shorthand, usage, def, typ string }
	seen := make(map[string][]occurrence)

	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if _, ok := flagsLocalByDesign[f.Name]; !ok {
			return
		}
		seen[f.Name] = append(seen[f.Name], occurrence{
			cmd: cmdName, shorthand: f.Shorthand, usage: f.Usage,
			def: f.DefValue, typ: f.Value.Type(),
		})
	})

	if len(seen) != len(flagsLocalByDesign) {
		t.Errorf("expected every flagsLocalByDesign entry to appear; got %d of %d",
			len(seen), len(flagsLocalByDesign))
	}

	for name, occs := range seen {
		first := occs[0]
		for _, o := range occs[1:] {
			if o.shorthand != first.shorthand {
				t.Errorf("%q: shorthand %q in %s vs %q in %s",
					name, first.shorthand, first.cmd, o.shorthand, o.cmd)
			}
			switch name {
			case "variant":
				// Type may differ (string vs stringSlice); the default follows
				// from the type, so it is exempt too. Usage must still match.
				if !strings.HasPrefix(o.usage, "Path to a VCF/variant file") {
					t.Errorf("%q usage in %s should describe a VCF path: %q", name, o.cmd, o.usage)
				}
			case "bqsr":
				// Default may differ; type and usage must not.
				if o.typ != first.typ {
					t.Errorf("%q: type %s in %s vs %s in %s", name, first.typ, first.cmd, o.typ, o.cmd)
				}
				if o.usage != first.usage {
					t.Errorf("%q: usage differs between %s and %s", name, first.cmd, o.cmd)
				}
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Every flag is consumed
// ---------------------------------------------------------------------------

// TestEveryFlagIsRead guards against flags that are declared and then never
// consumed. --jobs and --bulk_size were both in that state before the cleanup:
// genome-whisperer accepted them and silently ignored them.
func TestEveryFlagIsRead(t *testing.T) {
	read := flagNamesReadInPackage(t)

	rootCmd.PersistentFlags().VisitAll(func(f *pflag.Flag) {
		if !read[f.Name] {
			t.Errorf("persistent flag %q is never read via Flags().Get*", f.Name)
		}
	})
	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Name == "help" {
			return
		}
		if !read[f.Name] {
			t.Errorf("%s: flag %q is never read via Flags().Get*", cmdName, f.Name)
		}
	})
}

// ---------------------------------------------------------------------------
// Help text
// ---------------------------------------------------------------------------

// TestUsageStringsAreMeaningful catches the copy-paste help strings that used to
// describe --species as "number of threads", --genomes-dir as "Skip
// verification" and --database as "Species name".
func TestUsageStringsAreMeaningful(t *testing.T) {
	banned := map[string][]string{
		"species":     {"thread"},
		"genomes-dir": {"skip verification"},
		"database":    {"species name"},
	}

	check := func(cmdName string, f *pflag.Flag) {
		if strings.TrimSpace(f.Usage) == "" {
			t.Errorf("%s: flag %q has an empty usage string", cmdName, f.Name)
			return
		}
		lower := strings.ToLower(f.Usage)
		for _, b := range banned[f.Name] {
			if strings.Contains(lower, b) {
				t.Errorf("%s: %q usage looks copy-pasted (contains %q): %q",
					cmdName, f.Name, b, f.Usage)
			}
		}
	}

	rootCmd.PersistentFlags().VisitAll(func(f *pflag.Flag) { check("<root persistent>", f) })
	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Name != "help" {
			check(cmdName, f)
		}
	})
}

// TestThresholdPrefixesMatchDirection pins the rename that fixed --min-FS-SNP
// and --min-SOR-SNP, which were maximums wearing a "min-" prefix.
func TestThresholdPrefixesMatchDirection(t *testing.T) {
	fs := VariantCallingCmd.LocalNonPersistentFlags()

	for _, name := range []string{"max-fs-snp", "max-sor-snp", "max-fs-indel", "max-sor-indel"} {
		f := fs.Lookup(name)
		if f == nil {
			t.Errorf("%q should be registered on VariantCalling", name)
			continue
		}
		if !strings.Contains(strings.ToLower(f.Usage), "maximum") {
			t.Errorf("%q is a ceiling; usage should say \"maximum\": %q", name, f.Usage)
		}
	}
	for _, name := range []string{
		"min-qd-snp", "min-qual-snp", "min-mq-snp", "min-mqranksum-snp",
		"min-readposranksum-snp", "min-qd-indel", "min-qual-indel",
		"min-readposranksum-indel",
	} {
		f := fs.Lookup(name)
		if f == nil {
			t.Errorf("%q should be registered on VariantCalling", name)
			continue
		}
		if !strings.Contains(strings.ToLower(f.Usage), "minimum") {
			t.Errorf("%q is a floor; usage should say \"minimum\": %q", name, f.Usage)
		}
	}
}

// ---------------------------------------------------------------------------
// Help ordering
// ---------------------------------------------------------------------------

// TestHelpUsesRegistrationOrder confirms SortFlags is off everywhere, so --help
// lists flags in the order they were declared rather than alphabetically.
func TestHelpUsesRegistrationOrder(t *testing.T) {
	if rootCmd.PersistentFlags().SortFlags {
		t.Error("rootCmd.PersistentFlags().SortFlags should be false")
	}
	for _, c := range allCommands() {
		if c.Name() == "CreateTemplate" || c.Name() == "help" {
			continue // no flags of its own
		}
		if c.Flags().SortFlags {
			t.Errorf("%s: Flags().SortFlags should be false", c.Name())
		}
	}
}

// ---------------------------------------------------------------------------
// Inventory (informational, printed by `go test -v`)
// ---------------------------------------------------------------------------

func TestShorthandInventory(t *testing.T) {
	byLetter := make(map[string]string)

	rootCmd.PersistentFlags().VisitAll(func(f *pflag.Flag) {
		if f.Shorthand != "" {
			byLetter[f.Shorthand] = f.Name + "  (global)"
		}
	})
	eachLocalFlag(func(cmdName string, f *pflag.Flag) {
		if f.Shorthand != "" && f.Name != "help" {
			byLetter[f.Shorthand] = f.Name
		}
	})

	letters := make([]string, 0, len(byLetter))
	for l := range byLetter {
		letters = append(letters, l)
	}
	sort.Strings(letters)
	for _, l := range letters {
		t.Logf("-%s  %s", l, byLetter[l])
	}
}

// ---------------------------------------------------------------------------
// Static analysis helper
// ---------------------------------------------------------------------------

// flagNamesReadInPackage parses the cmd package and returns the set of flag
// names passed to a Flags().Get* call.
func flagNamesReadInPackage(t *testing.T) map[string]bool {
	t.Helper()

	fset := token.NewFileSet()
	pkgs, err := parser.ParseDir(fset, ".", nil, 0)
	if err != nil {
		t.Fatalf("parsing cmd package: %v", err)
	}

	read := make(map[string]bool)
	for _, pkg := range pkgs {
		for _, f := range pkg.Files {
			ast.Inspect(f, func(n ast.Node) bool {
				call, ok := n.(*ast.CallExpr)
				if !ok || len(call.Args) == 0 {
					return true
				}
				sel, ok := call.Fun.(*ast.SelectorExpr)
				if !ok || !strings.HasPrefix(sel.Sel.Name, "Get") {
					return true
				}
				lit, ok := call.Args[0].(*ast.BasicLit)
				if !ok || lit.Kind != token.STRING {
					return true
				}
				if v, err := strconv.Unquote(lit.Value); err == nil {
					read[v] = true
				}
				return true
			})
		}
	}
	return read
}
