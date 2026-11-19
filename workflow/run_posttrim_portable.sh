#!/usr/bin/env bash
set -euo pipefail

PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RUN="${1:?Usage: run_posttrim_portable.sh SAMPLE RUN_ROOT [CONFIG]}"
RUN_ROOT="${2:?Usage: run_posttrim_portable.sh SAMPLE RUN_ROOT [CONFIG]}"
CONFIG="${3:-$PKG/config/default.env}"
PREDEPLETED_R1="${4:-}"
PREDEPLETED_R2="${5:-}"

if [ -n "$PREDEPLETED_R1$PREDEPLETED_R2" ] && { [ -z "$PREDEPLETED_R1" ] || [ -z "$PREDEPLETED_R2" ]; }; then
    echo "ERROR: both pre-depleted FASTQ files are required" >&2
    exit 2
fi

mkdir -p "$RUN_ROOT"
RUN_ROOT="$(readlink -f "$RUN_ROOT")"

export MYCOWES_CONFIG="$CONFIG"
source "$PKG/workflow/runtime.sh"

MM2="$MINIMAP2_BIN"
SAMTOOLS="$SAMTOOLS_BIN"

TRIMDIR="$RUN_ROOT/01_trim"
HOSTDIR="$RUN_ROOT/02_host_depletion"
JOINTDIR="$RUN_ROOT/03_joint_alignment"
LOGDIR="$RUN_ROOT/logs"
LOCKROOT="$RUN_ROOT/.locks"

TRIM_R1="$TRIMDIR/$RUN.trimmed_R1.fastq.gz"
TRIM_R2="$TRIMDIR/$RUN.trimmed_R2.fastq.gz"

UNMAPPED_BAM="$HOSTDIR/$RUN.trimmed.hg38_both_unmapped.qname.bam"
UNMAP_R1="$HOSTDIR/$RUN.trimmed.hg38_both_unmapped_R1.fastq.gz"
UNMAP_R2="$HOSTDIR/$RUN.trimmed.hg38_both_unmapped_R2.fastq.gz"

JOINT_BAM="$JOINTDIR/$RUN.trimmed.AB_CDE_joint.qname.bam"

MANIFEST="$RUN_ROOT/posttrim_manifest.txt"

mkdir -p \
    "$TRIMDIR" \
    "$HOSTDIR" \
    "$JOINTDIR" \
    "$LOGDIR" \
    "$LOCKROOT" \
    "$RUN_ROOT/tmp"

# Prevent duplicate execution for the same sample.
LOCKDIR="$LOCKROOT/$RUN.posttrim.lock"

if ! mkdir "$LOCKDIR" 2>/dev/null; then
    echo "ERROR: another optimized run for $RUN is active, or a stale lock exists:"
    echo "$LOCKDIR"
    exit 1
fi

TMPROOT="$RUN_ROOT/tmp/${RUN}.optimized.$$.$(date +%s)"
mkdir -p "$TMPROOT"

UNMAPPED_BAM_TMP="$UNMAPPED_BAM.tmp.$$"
JOINT_BAM_TMP="$JOINT_BAM.tmp.$$"

cleanup() {
    rm -f "$UNMAPPED_BAM_TMP" "$JOINT_BAM_TMP"
    rm -rf "$TMPROOT"
    rmdir "$LOCKDIR" 2>/dev/null || true
}
trap cleanup EXIT

required_files=("$MM2" "$SAMTOOLS" "$JOINT_MMI")
if [ -n "$PREDEPLETED_R1" ]; then
    PREDEPLETED_R1="$(readlink -f "$PREDEPLETED_R1")"
    PREDEPLETED_R2="$(readlink -f "$PREDEPLETED_R2")"
    required_files+=("$PREDEPLETED_R1" "$PREDEPLETED_R2")
else
    required_files+=("$HG38_MMI" "$TRIM_R1" "$TRIM_R2")
fi

for file in "${required_files[@]}"
do
    if [ ! -s "$file" ]; then
        echo "ERROR: missing or empty file: $file"
        exit 1
    fi
done

START_EPOCH=$(date +%s)

{
    echo "sample=$RUN"
    echo "start=$(date --iso-8601=seconds)"
    echo "pipeline=posttrim_optimized_v0.1"
    echo "input_stage=$([ -n "$PREDEPLETED_R1" ] && echo GRCh38_PAIRED_BOTH_UNMAPPED || echo TRIMMED_FASTQ)"
    echo "trimmed_r1=$TRIM_R1"
    echo "trimmed_r2=$TRIM_R2"
    echo "hg38_index=$([ -n "$PREDEPLETED_R1" ] && echo SKIPPED_PREDEPLETED_INPUT || echo "$HG38_MMI")"
    echo "joint_index=$JOINT_MMI"
    echo "host_strategy=filter_before_name_sort"
    echo "host_flags=-f12_-F2304"
    echo "host_secondary=minimap2_default_filtered_by_-F2304"
    echo "temporary_directory=$TMPROOT"
} > "$MANIFEST"

if [ -n "$PREDEPLETED_R1" ]; then
    echo "===== STEP 1: HG38 DEPLETION ALREADY COMPLETE -> SKIP ====="
    gzip -t "$PREDEPLETED_R1" "$PREDEPLETED_R2"
    for pair in "$PREDEPLETED_R1:$UNMAP_R1" "$PREDEPLETED_R2:$UNMAP_R2"; do
        source_file="${pair%%:*}"
        target_file="${pair#*:}"
        tmp_link="$target_file.tmp.$$"
        rm -f "$tmp_link"
        ln -s "$source_file" "$tmp_link"
        mv -Tf "$tmp_link" "$target_file"
    done
    echo "host_depletion=SKIPPED_PREDEPLETED_INPUT" >> "$MANIFEST"
    UNMAP_PAIRS=$(gzip -cd "$UNMAP_R1" | awk 'END {print int(NR/4)}')
else
    echo "===== STEP 1: HG38 ALIGNMENT AND EARLY FILTERING ====="

# Critical optimization:
# Filter to primary read pairs where read and mate are both unmapped
# before name sorting.
"$MM2" \
    -ax sr \
    -t "$HOST_THREADS" \
    "$HG38_MMI" \
    "$TRIM_R1" \
    "$TRIM_R2" \
    2> "$LOGDIR/$RUN.optimized.hg38.minimap2.log" \
| "$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -u \
    -f 12 \
    -F 2304 \
    - \
| "$SAMTOOLS" sort \
    -@ "$SAMTOOLS_THREADS" \
    -n \
    -m 512M \
    -T "$TMPROOT/hg38_unmapped_sort" \
    -o "$UNMAPPED_BAM_TMP" \
    -

"$SAMTOOLS" quickcheck -v "$UNMAPPED_BAM_TMP"

mv -f "$UNMAPPED_BAM_TMP" "$UNMAPPED_BAM"

"$SAMTOOLS" flagstat \
    -@ "$SAMTOOLS_THREADS" \
    "$UNMAPPED_BAM" \
    > "$LOGDIR/$RUN.optimized.hg38_unmapped.flagstat.txt"

UNMAP_PAIRS=$(
    "$SAMTOOLS" view \
        -c \
        -f 64 \
        "$UNMAPPED_BAM"
)

echo "===== STEP 2: OUTPUT HUMAN-UNMAPPED FASTQ ====="

rm -f "$UNMAP_R1" "$UNMAP_R2"

"$SAMTOOLS" fastq \
    -@ "$SAMTOOLS_THREADS" \
    -n \
    -1 "$UNMAP_R1" \
    -2 "$UNMAP_R2" \
    -0 /dev/null \
    -s /dev/null \
    "$UNMAPPED_BAM" \
    2> "$LOGDIR/$RUN.optimized.samtools_fastq.log"

gzip -t "$UNMAP_R1"
gzip -t "$UNMAP_R2"
fi

echo "human_unmapped_pairs=$UNMAP_PAIRS" | tee -a "$MANIFEST"

echo "===== STEP 3: AB+CDE JOINT ALIGNMENT ====="

if [ -n "$PREDEPLETED_R1" ] && [ -s "$JOINT_BAM" ] && "$SAMTOOLS" quickcheck "$JOINT_BAM" 2>/dev/null; then
    echo "JOINT ALIGNMENT already complete -> SKIP"
else
"$MM2" \
    -ax sr \
    --secondary=yes \
    -N 100 \
    --MD \
    -R "@RG\tID:${RUN}_trimmed\tSM:${RUN}\tPL:ILLUMINA" \
    -t "$MICROBIAL_THREADS" \
    "$JOINT_MMI" \
    "$UNMAP_R1" \
    "$UNMAP_R2" \
    2> "$LOGDIR/$RUN.optimized.AB_CDE.minimap2.log" \
| "$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -b \
    - \
| "$SAMTOOLS" sort \
    -@ "$SAMTOOLS_THREADS" \
    -n \
    -m 512M \
    -T "$TMPROOT/joint_sort" \
    -o "$JOINT_BAM_TMP" \
    -

"$SAMTOOLS" quickcheck -v "$JOINT_BAM_TMP"

mv -f "$JOINT_BAM_TMP" "$JOINT_BAM"

"$SAMTOOLS" flagstat \
    -@ "$SAMTOOLS_THREADS" \
    "$JOINT_BAM" \
    > "$LOGDIR/$RUN.optimized.AB_CDE.flagstat.txt"
fi

END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))

{
    echo "human_unmapped_bam=$UNMAPPED_BAM"
    echo "human_unmapped_r1=$UNMAP_R1"
    echo "human_unmapped_r2=$UNMAP_R2"
    echo "joint_bam=$JOINT_BAM"
    echo "elapsed_seconds=$ELAPSED"
    echo "finish=$(date --iso-8601=seconds)"
    echo "status=COMPLETE"
} >> "$MANIFEST"

rm -rf "$TMPROOT"

echo
echo "===== OPTIMIZED POST-TRIM PIPELINE COMPLETE ====="
cat "$MANIFEST"
