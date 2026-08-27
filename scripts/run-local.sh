#!/usr/bin/env bash
#
# run-local.sh — stage a genome-whisperer AlignReads run on local disk, run it,
# publish the results back to the NAS.
#
# The NAS shares are 9p mounts under WSL2 and read at roughly 20-35 MB/s, while
# the ext4 root disk does well over 1 GB/s. A BQSR run reads each CRAM several
# times over and writes an uncompressed BAM the size of two of them, so running
# it in place on the NAS spends most of its wall clock waiting on the mount.
# Copying the inputs down once, working locally, and copying the finished CRAMs
# back moves far fewer bytes over 9p than running in place does.
#
# Samples are processed in batches so that peak local disk stays bounded: a
# batch is staged, run, published, and removed before the next one starts.
#
# Usage:
#   scripts/run-local.sh -s onion -R GCA_030765085.1_ASM3076508v1 -b 4 \
#       -- --bqsr --bootstrap --quick -t 20
#
# Everything after -- is passed to genome-whisperer verbatim.

set -euo pipefail

# ----------------------------------------------------------------- defaults --

SPECIES=""
REFVER=""
NAS_DATA="/mnt/v/DATA"
NAS_GENOMES="/mnt/z/genomes"
WORK="${HOME}/gw-work"
SAMPLES=""
BATCH=0
WITH_READS=0
KEEP_LOCAL=0
DRY_RUN=0
DO_INSTALL=0
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GW_ARGS=()

# Local space to keep free on top of what a batch is expected to need, so the
# root filesystem never fills to the point of breaking other things.
RESERVE_GB=50

usage() {
	sed -n '3,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
	cat <<'EOF'

Options:
  -s SPECIES        species directory name under the data and genomes roots (required)
  -R REFVER         reference version directory name (required)
  -d DIR            NAS data root                   (default /mnt/v/DATA)
  -g DIR            NAS genomes root                (default /mnt/z/genomes)
  -w DIR            local work directory            (default $HOME/gw-work)
  -S LIST           comma-separated samples to run  (default: all found)
  -b N              samples per batch, 0 = all      (default 0)
  --with-reads      also stage the clean_reads FASTQs (needed to align from scratch)
  --keep-local      keep the local copy after publishing
  --install         run "go install ." in the repo before starting
  --dry-run         print what would happen, touch nothing
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

run() {
	if ((DRY_RUN)); then
		printf '   \033[90m[dry-run]\033[0m %s\n' "$*"
		return 0
	fi
	"$@"
}

# -------------------------------------------------------------------- flags --

while (($#)); do
	case "$1" in
	-s) SPECIES="$2"; shift 2 ;;
	-R) REFVER="$2"; shift 2 ;;
	-d) NAS_DATA="$2"; shift 2 ;;
	-g) NAS_GENOMES="$2"; shift 2 ;;
	-w) WORK="$2"; shift 2 ;;
	-S) SAMPLES="$2"; shift 2 ;;
	-b) BATCH="$2"; shift 2 ;;
	--with-reads) WITH_READS=1; shift ;;
	--keep-local) KEEP_LOCAL=1; shift ;;
	--install) DO_INSTALL=1; shift ;;
	--dry-run) DRY_RUN=1; shift ;;
	-h | --help) usage; exit 0 ;;
	--) shift; GW_ARGS=("$@"); break ;;
	*) die "unknown option: $1 (use -- before genome-whisperer's own flags)" ;;
	esac
done

[[ -n $SPECIES ]] || { usage; die "-s SPECIES is required"; }
[[ -n $REFVER ]] || { usage; die "-R REFVER is required"; }
[[ $BATCH =~ ^[0-9]+$ ]] || die "-b takes a number"

# Whether BQSR was asked for decides which outputs a sample must have produced
# before its results are published.
WANT_BQSR=0
for arg in "${GW_ARGS[@]}"; do
	[[ $arg == "--bqsr" ]] && WANT_BQSR=1
done

# ---------------------------------------------------------------- preflight --

for tool in rsync samtools numfmt flock; do
	command -v "$tool" >/dev/null || die "$tool is not on PATH"
done

if ((DO_INSTALL)); then
	info "Installing genome-whisperer from $REPO"
	run env -C "$REPO" go install .
fi

GW="$(command -v genome-whisperer || true)"
[[ -n $GW ]] || die "genome-whisperer is not on PATH (run with --install, or go install . in $REPO)"

# A binary older than the source it was built from is the classic way to spend a
# day watching the bug you already fixed.
if [[ -d $REPO ]]; then
	newest_src="$(find "$REPO" -name '*.go' -newer "$GW" -not -path '*/old/*' -print -quit 2>/dev/null || true)"
	[[ -n $newest_src ]] && warn "$GW is older than $newest_src — consider --install"
fi

[[ -d $NAS_DATA ]] || die "data root not found: $NAS_DATA"
[[ -d $NAS_GENOMES ]] || die "genomes root not found: $NAS_GENOMES"

# ------------------------------------------------------------ the reference --

ASSEMBLY_DIR="$NAS_GENOMES/$SPECIES/$REFVER/assembly"
[[ -d $ASSEMBLY_DIR ]] || die "assembly directory not found: $ASSEMBLY_DIR"

NAS_FASTA=""
for candidate in "$ASSEMBLY_DIR"/*.fa "$ASSEMBLY_DIR"/*.fasta "$ASSEMBLY_DIR"/*.fna; do
	[[ -f $candidate ]] || continue
	# A reference is only usable here if it carries both the .fai the shard
	# splitter reads and the .dict GATK insists on.
	[[ -f ${candidate}.fai && -f ${candidate%.*}.dict ]] || continue
	NAS_FASTA="$candidate"
	break
done
[[ -n $NAS_FASTA ]] || die "no FASTA with both a .fai and a .dict in $ASSEMBLY_DIR"

REF_STEM="$(basename "${NAS_FASTA%.*}")"
REF_EXT="${NAS_FASTA##*.}"
LOCAL_REF_DIR="$WORK/ref/$SPECIES/$REFVER"
LOCAL_FASTA="$LOCAL_REF_DIR/$REF_STEM.$REF_EXT"

REF_BASES="$(awk '{s+=$2} END {printf "%d", s}' "${NAS_FASTA}.fai")"
info "Reference: $NAS_FASTA"
printf '     %s bases, %s contigs\n' \
	"$(numfmt --to=si "$REF_BASES")" "$(wc -l <"${NAS_FASTA}.fai")"

# ----------------------------------------------------------- find the samples --

mapfile -t READ_DIRS < <(find "$NAS_DATA/$SPECIES" -mindepth 3 -maxdepth 3 -type d -name clean_reads 2>/dev/null | sort)
((${#READ_DIRS[@]})) || die "no <year>/<sample>/clean_reads directories under $NAS_DATA/$SPECIES"

declare -a WANTED=()
for read_dir in "${READ_DIRS[@]}"; do
	sample_dir="$(dirname "$read_dir")"
	sample="$(basename "$sample_dir")"

	# Long-read samples are skipped by the aligner itself, so there is no point
	# moving their data across the network.
	if [[ ${sample^^} == *LR ]]; then
		warn "skipping $sample (long-read sample, unsupported by AlignReads)"
		continue
	fi
	if [[ -n $SAMPLES && ",$SAMPLES," != *",$sample,"* ]]; then
		continue
	fi
	WANTED+=("$sample_dir")
done
((${#WANTED[@]})) || die "no samples selected"

info "Samples selected: ${#WANTED[@]}"

((BATCH > 0)) || BATCH=${#WANTED[@]}

# ------------------------------------------------------------------ helpers --

# copy_in SRC DST — one file, skipped when the destination already matches.
copy_in() {
	local src="$1" dst="$2"
	[[ -f $src ]] || return 0
	if [[ -f $dst ]] && [[ "$(stat -c%s "$src")" == "$(stat -c%s "$dst")" ]]; then
		ok "already staged: $(basename "$dst")"
		return 0
	fi
	run rsync -h --partial --inplace --info=progress2 "$src" "$dst"
}

# index_is_sound PATH — true when an index can actually be read. Both .crai and
# .csi are BGZF, so a stream truncated by an interrupted run fails here; there
# is no point copying one across the network for the pipeline to discard.
index_is_sound() {
	local idx="$1"
	[[ -s $idx ]] || return 1
	case "$idx" in
	*.crai | *.csi) gzip -t "$idx" 2>/dev/null ;;
	*) samtools idxstats "${idx%.*}" >/dev/null 2>&1 ;;
	esac
}

# dir_bytes DIR [find-args...] — total size of matching files, 0 if none.
dir_bytes() {
	local dir="$1"; shift
	[[ -d $dir ]] || { echo 0; return; }
	find "$dir" -maxdepth 1 -type f "$@" -printf '%s\n' 2>/dev/null |
		awk '{s+=$1} END {printf "%d", s+0}'
}

# --------------------------------------------------------------- the batches --

mkdir -p "$WORK" "$WORK/logs"
exec 9>"$WORK/.lock"
flock -n 9 || die "another run-local.sh is using $WORK"

STAMP="$(date +%Y%m%d-%H%M%S)"
RUN_LOG="$WORK/logs/run-$STAMP.log"
declare -a PUBLISHED=() FAILED=() SKIPPED=()

total_batches=$(((${#WANTED[@]} + BATCH - 1) / BATCH))
batch_no=0

LOCAL_DATA="$WORK/batch/data"

# The reference is staged once and reused by every batch.
info "Staging reference to $LOCAL_REF_DIR"
run mkdir -p "$LOCAL_REF_DIR"
copy_in "$NAS_FASTA" "$LOCAL_FASTA"
copy_in "${NAS_FASTA}.fai" "${LOCAL_FASTA}.fai"
copy_in "${NAS_FASTA%.*}.dict" "${LOCAL_FASTA%.*}.dict"

for ((offset = 0; offset < ${#WANTED[@]}; offset += BATCH)); do
	batch=("${WANTED[@]:offset:BATCH}")
	batch_no=$((batch_no + 1))

	printf '\n\033[1m=== batch %d/%d: %s ===\033[0m\n' \
		"$batch_no" "$total_batches" "$(basename -a "${batch[@]}" | paste -sd' ')"

	# ---- work out what this batch will cost on local disk ----
	staged_bytes=0
	for sample_dir in "${batch[@]}"; do
		nas_bams="$sample_dir/reference_genomes/$REFVER/bams"
		staged_bytes=$((staged_bytes + $(dir_bytes "$nas_bams" -name '*.cram')))
		((WITH_READS)) && staged_bytes=$((staged_bytes + $(dir_bytes "$sample_dir/clean_reads")))
	done

	# On a reference whose contigs run past the BAI limit, GATK cannot read the
	# staged CRAM at all, so the pipeline decodes it back to a BAM first (see
	# gatkReadsNeedBam). At the peak a sample then holds its CRAM, that decoded
	# BAM, the uncompressed BAM ApplyBQSR writes, and the recalibrated CRAM, each
	# BAM roughly twice its CRAM — call it six times what was staged, plus the
	# bootstrap shard VCFs. Deliberately generous: a wrong guess here costs a
	# run, not a rerun.
	need_gb=$(((staged_bytes * 6) / 1000000000 + RESERVE_GB))
	avail_gb=$(($(df --output=avail -B1G "$WORK" | tail -1)))
	printf '     staging %s, expecting to need ~%d GB, %d GB free\n' \
		"$(numfmt --to=si "$staged_bytes")" "$need_gb" "$avail_gb"
	if ((need_gb > avail_gb)) && ((!DRY_RUN)); then
		die "not enough local space for this batch: need ~${need_gb} GB, have ${avail_gb} GB. Use a smaller -b."
	fi

	# ---- stage each sample ----
	declare -a staged=()
	for sample_dir in "${batch[@]}"; do
		sample="$(basename "$sample_dir")"
		year="$(basename "$(dirname "$sample_dir")")"
		nas_bams="$sample_dir/reference_genomes/$REFVER/bams"
		local_sample="$LOCAL_DATA/$SPECIES/$year/$sample"
		local_bams="$local_sample/reference_genomes/$REFVER/bams"

		info "Staging $sample"
		run mkdir -p "$local_bams" "$local_sample/clean_reads"

		# CRAMs the pipeline can resume from, with their indexes when sound.
		shopt -s nullglob
		for cram in "$nas_bams"/*.cram; do
			copy_in "$cram" "$local_bams/$(basename "$cram")"
			for idx in "$cram.crai" "${cram%.cram}.crai"; do
				[[ -f $idx ]] || continue
				if index_is_sound "$idx"; then
					copy_in "$idx" "$local_bams/$(basename "$idx")"
				else
					warn "$(basename "$idx") is truncated — leaving it behind, it will be rebuilt locally"
				fi
			done
		done

		if ((WITH_READS)); then
			for fq in "$sample_dir"/clean_reads/*.f*q.gz "$sample_dir"/clean_reads/*.f*q; do
				copy_in "$fq" "$local_sample/clean_reads/$(basename "$fq")"
			done
		fi
		shopt -u nullglob

		# Nothing to work from is a setup mistake worth naming now rather than
		# letting the aligner report it per sample later.
		if ((!WITH_READS)) && ! compgen -G "$local_bams/*.cram" >/dev/null && ((!DRY_RUN)); then
			warn "$sample has no CRAM staged and --with-reads was not given — it will be skipped"
			SKIPPED+=("$sample")
			continue
		fi
		staged+=("$sample_dir")
	done

	((${#staged[@]})) || { warn "nothing staged for this batch"; continue; }

	# ---- run ----
	cmd=("$GW" AlignReads -d "$LOCAL_DATA" -r "$LOCAL_FASTA" -s "$SPECIES" --ref-version "$REFVER")
	cmd+=("${GW_ARGS[@]}")
	info "Running: ${cmd[*]}"
	printf '     log: %s\n' "$RUN_LOG"
	if ((DRY_RUN)); then
		printf '   \033[90m[dry-run]\033[0m would run the above\n'
	else
		set +e
		"${cmd[@]}" 2>&1 | tee -a "$RUN_LOG"
		rc=${PIPESTATUS[0]}
		set -e
		((rc == 0)) || warn "genome-whisperer exited $rc — per-sample results are still checked below"
	fi

	# ---- publish whatever finished ----
	for sample_dir in "${staged[@]}"; do
		sample="$(basename "$sample_dir")"
		year="$(basename "$(dirname "$sample_dir")")"
		local_bams="$LOCAL_DATA/$SPECIES/$year/$sample/reference_genomes/$REFVER/bams"
		nas_bams="$sample_dir/reference_genomes/$REFVER/bams"

		# A sample is only published once the CRAMs it was supposed to produce
		# are there and readable. Publishing a half-written CRAM over a good one
		# on the NAS is the one outcome worth going out of the way to prevent.
		missing=""
		for expected in "$local_bams/$sample.RGMD.cram"; do
			[[ -f $expected ]] || missing="$missing $(basename "$expected")"
		done
		if ((WANT_BQSR)); then
			[[ -f $local_bams/${sample}.RGMD_bqsr.cram ]] || missing="$missing ${sample}.RGMD_bqsr.cram"
		fi
		if [[ -n $missing ]] && ((!DRY_RUN)); then
			warn "$sample did not finish (missing:$missing) — not publishing"
			FAILED+=("$sample")
			continue
		fi

		for cram in "$local_bams"/*.cram; do
			[[ -f $cram ]] || continue
			if ! samtools quickcheck "$cram" 2>/dev/null; then
				warn "$sample: $(basename "$cram") fails quickcheck — not publishing"
				missing="$missing $(basename "$cram")"
			fi
		done
		if [[ -n $missing ]] && ((!DRY_RUN)); then
			FAILED+=("$sample")
			continue
		fi

		info "Publishing $sample to $nas_bams"
		run mkdir -p "$nas_bams"
		# Only the durable artefacts go back: the BAM intermediates are large,
		# already deleted by the pipeline's own cleanup on success, and of no use
		# on the NAS. --partial lets an interrupted transfer resume.
		# -r to descend into the directory at all, -t to keep modification
		# times so an index stays newer than the file it indexes.
		run rsync -rt -h --partial --info=progress2 \
			--exclude='shards/' --exclude='tmp/' \
			--exclude='*.bam' --exclude='*.bam.bai' --exclude='*.bam.csi' \
			--include='*.cram' --include='*.crai' --include='*.txt' --include='*.pdf' \
			--exclude='*' \
			"$local_bams"/ "$nas_bams"/

		# Confirm the copy landed: sizes match and the CRAM still reads on the
		# far side. Reading a header over 9p is cheap; a silent short write is not.
		published_ok=1
		for cram in "$local_bams"/*.cram; do
			[[ -f $cram ]] || continue
			remote="$nas_bams/$(basename "$cram")"
			if ((DRY_RUN)); then continue; fi
			if [[ ! -f $remote ]] || [[ "$(stat -c%s "$cram")" != "$(stat -c%s "$remote")" ]]; then
				warn "$sample: $(basename "$cram") did not copy back intact"
				published_ok=0
				continue
			fi
			samtools quickcheck "$remote" 2>/dev/null || {
				warn "$sample: published $(basename "$cram") fails quickcheck on the NAS"
				published_ok=0
			}
		done

		if ((DRY_RUN)); then
			ok "$sample would be published"
			PUBLISHED+=("$sample")
		elif ((published_ok)); then
			ok "$sample published"
			PUBLISHED+=("$sample")
			if ((!KEEP_LOCAL)); then
				info "Removing local copy of $sample"
				run rm -rf "$LOCAL_DATA/$SPECIES/$year/$sample"
			fi
		else
			FAILED+=("$sample")
		fi
	done
done

# ------------------------------------------------------------------ summary --

printf '\n\033[1m=== summary ===\033[0m\n'
printf 'published: %d%s\n' "${#PUBLISHED[@]}" \
	"$( ((${#PUBLISHED[@]})) && printf ' — %s' "$(printf '%s ' "${PUBLISHED[@]}")")"
if ((${#FAILED[@]})); then
	printf '\033[31mfailed:    %d — %s\033[0m\n' "${#FAILED[@]}" "$(printf '%s ' "${FAILED[@]}")"
fi
if ((${#SKIPPED[@]})); then
	printf '\033[33mskipped:   %d — %s\033[0m\n' "${#SKIPPED[@]}" "$(printf '%s ' "${SKIPPED[@]}")"
fi
printf 'log:       %s\n' "$RUN_LOG"
((KEEP_LOCAL)) && printf 'local copy kept under %s\n' "$WORK/batch"

((${#FAILED[@]} == 0))
