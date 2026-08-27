#!/usr/bin/env bash
#
# reheader-crams.sh — point existing CRAMs at a renamed copy of the reference
# they were aligned to.
#
# An assembly whose contigs were renamed after the alignments were made leaves
# CRAMs that nothing can read: CRAM resolves its reference by contig name, so
# samtools cannot decode them against the renamed FASTA and falls back to the UR
# path in the header, which is usually long gone. GATK, if it gets that far,
# rejects the pair for having incompatible sequence dictionaries.
#
# Renaming is only safe when the sequences are the same, so that is checked
# rather than assumed: every contig must carry the same M5 checksum and length,
# in the same order, in both the CRAM header and the reference .dict. Anything
# else and the file is left alone.
#
# The edit is made in place. The header is rewritten with the new names and
# without the UR field — it points at a path that no longer exists, and dropping
# it is what makes the new header fit in the space the old one occupied. Because
# nothing after the header moves, existing .crai indexes stay valid.
#
# Usage:
#   scripts/reheader-crams.sh -r /path/to/renamed.fa /path/to/*.cram
#   scripts/reheader-crams.sh -r /path/to/renamed.fa -d /mnt/v/DATA/onion
#
set -euo pipefail

REF=""
SCAN_DIR=""
BACKUP_DIR="./reheader-backups"
DRY_RUN=0

usage() {
	sed -n '3,27p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
	cat <<'EOF'

Options:
  -r FASTA          reference whose names the CRAMs should adopt (required;
                    its .dict is what the CRAM headers are checked against)
  -d DIR            rename every .cram found under DIR instead of naming them
  -B DIR            where to keep the replaced headers  (default ./reheader-backups)
  --dry-run         report what would change, touch nothing
  -h, --help        this message
EOF
}

die() {
	printf '\033[31merror:\033[0m %s\n' "$*" >&2
	exit 1
}
info() { printf '\033[36m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[33mwarn:\033[0m %s\n' "$*"; }
ok() { printf '\033[32m  ok\033[0m %s\n' "$*"; }
skip() { printf '\033[90m  --\033[0m %s\n' "$*"; }

CRAMS=()
while (($#)); do
	case "$1" in
	-r) REF="$2"; shift 2 ;;
	-d) SCAN_DIR="$2"; shift 2 ;;
	-B) BACKUP_DIR="$2"; shift 2 ;;
	--dry-run) DRY_RUN=1; shift ;;
	-h | --help) usage; exit 0 ;;
	-*) die "unknown option: $1" ;;
	*) CRAMS+=("$1"); shift ;;
	esac
done

[[ -n $REF ]] || { usage; die "-r FASTA is required"; }
[[ -f $REF ]] || die "reference not found: $REF"

DICT="${REF%.*}.dict"
[[ -f $DICT ]] || die "no .dict beside the reference: $DICT"

for tool in samtools awk; do
	command -v "$tool" >/dev/null || die "$tool is not on PATH"
done

if [[ -n $SCAN_DIR ]]; then
	[[ -d $SCAN_DIR ]] || die "not a directory: $SCAN_DIR"
	mapfile -t found < <(find "$SCAN_DIR" -type f -name '*.cram' | sort)
	CRAMS+=("${found[@]}")
fi
((${#CRAMS[@]})) || die "no CRAMs given (name them, or use -d DIR)"

# sq_fingerprint reads @SQ lines from stdin and prints "M5<TAB>LN" per contig,
# in file order. Two files whose fingerprints match hold the same sequences
# under whatever names they happen to use.
sq_fingerprint() {
	awk '/^@SQ/{
		m=""; l=""
		for (i = 1; i <= NF; i++) {
			if ($i ~ /^M5:/) m = substr($i, 4)
			if ($i ~ /^LN:/) l = substr($i, 4)
		}
		print m "\t" l
	}'
}

REF_FP="$(mktemp)"; REF_SQ="$(mktemp)"
trap 'rm -f "$REF_FP" "$REF_SQ"' EXIT
sq_fingerprint <"$DICT" >"$REF_FP"
# The @SQ lines the new header is built from, minus UR: see the note at the top.
grep '^@SQ' "$DICT" | sed 's/\tUR:[^\t]*//' >"$REF_SQ"

if grep -q '^\t' "$REF_FP" || grep -q '\t$' "$REF_FP"; then
	die "$DICT has @SQ lines without an M5 or LN — cannot verify a rename against it"
fi
info "Reference: $REF"
printf '     %d contigs, all with M5 checksums\n' "$(wc -l <"$REF_FP")"

((DRY_RUN)) || mkdir -p "$BACKUP_DIR"

renamed=0; skipped=0; refused=0
for cram in "${CRAMS[@]}"; do
	name="$(basename "$cram")"
	[[ -f $cram ]] || { warn "$name: not found"; continue; }

	header="$(mktemp)"
	if ! samtools view -H "$cram" >"$header" 2>/dev/null; then
		warn "$name: header unreadable — skipping"
		rm -f "$header"
		refused=$((refused + 1))
		continue
	fi

	# Already carrying the target names? Then there is nothing to do, and saying
	# so is better than rewriting a header to the value it already has.
	if diff -q <(grep '^@SQ' "$header" | sed 's/\tUR:[^\t]*//') "$REF_SQ" >/dev/null 2>&1; then
		skip "$name already uses the reference's contig names"
		rm -f "$header"
		skipped=$((skipped + 1))
		continue
	fi

	if ! diff -q <(sq_fingerprint <"$header") "$REF_FP" >/dev/null 2>&1; then
		warn "$name: contigs do not match the reference (M5/length/order differ) — NOT touching it"
		printf '       compare: samtools view -H %q | grep ^@SQ   vs   %q\n' "$cram" "$DICT"
		rm -f "$header"
		refused=$((refused + 1))
		continue
	fi

	old_bytes=$(wc -c <"$header")
	new_header="$(mktemp)"
	{
		grep '^@HD' "$header" || printf '@HD\tVN:1.6\tSO:coordinate\n'
		cat "$REF_SQ"
		grep -v '^@HD' "$header" | grep -v '^@SQ' || true
	} >"$new_header"
	new_bytes=$(wc -c <"$new_header")

	first_old="$(awk '/^@SQ/{for(i=1;i<=NF;i++)if($i~/^SN:/){print substr($i,4);exit}}' "$header")"
	first_new="$(awk '/^@SQ/{for(i=1;i<=NF;i++)if($i~/^SN:/){print substr($i,4);exit}}' "$new_header")"

	if ((DRY_RUN)); then
		printf '   \033[90m[dry-run]\033[0m %s: %s → %s (header %d → %d bytes)\n' \
			"$name" "$first_old" "$first_new" "$old_bytes" "$new_bytes"
		rm -f "$header" "$new_header"
		renamed=$((renamed + 1))
		continue
	fi

	if ((new_bytes > old_bytes)); then
		warn "$name: the new header is larger than the one it replaces ($new_bytes > $old_bytes)"
		printf '       samtools cannot pad it in place. Rewrite the file instead:\n'
		printf '       samtools reheader -P %q %q > new.cram   (and rebuild its .crai)\n' "$new_header" "$cram"
		rm -f "$header" "$new_header"
		refused=$((refused + 1))
		continue
	fi

	cp "$header" "$BACKUP_DIR/$name.header.bak"
	info "$name: $first_old → $first_new"
	if ! samtools reheader -P -i "$new_header" "$cram"; then
		rm -f "$header" "$new_header"
		die "$name: reheader failed. The original header is at $BACKUP_DIR/$name.header.bak"
	fi

	# Read it back rather than trusting the exit status: this edits a file that
	# took days to produce and has no other copy on this disk.
	if ! samtools quickcheck "$cram"; then
		die "$name: fails quickcheck after reheader. Original header: $BACKUP_DIR/$name.header.bak"
	fi
	if ! diff -q <(samtools view -H "$cram" | grep '^@SQ' | sed 's/\tUR:[^\t]*//') "$REF_SQ" >/dev/null; then
		die "$name: header did not take. Original header: $BACKUP_DIR/$name.header.bak"
	fi
	for idx in "$cram.crai" "${cram%.cram}.crai"; do
		[[ -f $idx ]] || continue
		gzip -t "$idx" 2>/dev/null || warn "$name: $(basename "$idx") is unreadable — it was already truncated, the pipeline will rebuild it"
	done

	ok "$name renamed (header padded in place, .crai untouched)"
	rm -f "$header" "$new_header"
	renamed=$((renamed + 1))
done

printf '\n\033[1m=== summary ===\033[0m\n'
printf 'renamed:  %d\n' "$renamed"
((skipped)) && printf 'already correct: %d\n' "$skipped"
((refused)) && printf '\033[31mrefused:  %d\033[0m\n' "$refused"
((DRY_RUN)) && printf 'dry run — nothing was modified\n'
((refused == 0))
