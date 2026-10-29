#!/usr/bin/env bash
set -euo pipefail

PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
    echo "Usage:"
    echo "  $0 --sample SAMPLE --r1 R1.fastq.gz --r2 R2.fastq.gz --run-root DIR [--config FILE] [--nominal-pe-length BP] [--host-depleted]"
}

SAMPLE=""
R1=""
R2=""
RUN_ROOT=""
CONFIG="$PKG/config/default.env"
HOST_DEPLETED=0
NOMINAL_PE_LENGTH=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --sample)   SAMPLE="$2"; shift 2 ;;
        --r1)       R1="$2"; shift 2 ;;
        --r2)       R2="$2"; shift 2 ;;
        --run-root) RUN_ROOT="$2"; shift 2 ;;
        --config)   CONFIG="$2"; shift 2 ;;
        --nominal-pe-length) NOMINAL_PE_LENGTH="$2"; shift 2 ;;
        --host-depleted) HOST_DEPLETED=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

[ -n "$SAMPLE" ] || { echo "ERROR: --sample required" >&2; exit 2; }
[ -n "$R1" ] || { echo "ERROR: --r1 required" >&2; exit 2; }
[ -n "$R2" ] || { echo "ERROR: --r2 required" >&2; exit 2; }
[ -n "$RUN_ROOT" ] || { echo "ERROR: --run-root required" >&2; exit 2; }

[ -s "$R1" ] || { echo "ERROR: R1 missing: $R1" >&2; exit 3; }
[ -s "$R2" ] || { echo "ERROR: R2 missing: $R2" >&2; exit 3; }
[ -s "$CONFIG" ] || { echo "ERROR: config missing: $CONFIG" >&2; exit 3; }

R1="$(readlink -f "$R1")"
R2="$(readlink -f "$R2")"

mkdir -p "$RUN_ROOT"
RUN_ROOT="$(readlink -f "$RUN_ROOT")"
CONFIG="$(readlink -f "$CONFIG")"

export MYCOWES_CONFIG="$CONFIG"
source "$PKG/workflow/runtime.sh"

[[ "${CORE_MIN_ALIGNED_BP:-}" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: CORE_MIN_ALIGNED_BP must be a positive integer" >&2
    exit 3
}

if [ -z "$NOMINAL_PE_LENGTH" ]; then
    NOMINAL_PE_LENGTH=100
    NOMINAL_PE_LENGTH_SOURCE="LEGACY_DEFAULT"
elif [[ "$NOMINAL_PE_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
    NOMINAL_PE_LENGTH_SOURCE="EXPLICIT"
else
    echo "ERROR: --nominal-pe-length must be a positive integer" >&2
    exit 2
fi

if [ "$CORE_MIN_ALIGNED_BP" -lt "$NOMINAL_PE_LENGTH" ]; then
    EFFECTIVE_MIN_ALIGNED_BP="$CORE_MIN_ALIGNED_BP"
else
    EFFECTIVE_MIN_ALIGNED_BP="$NOMINAL_PE_LENGTH"
fi

PARSER="$PKG/workflow/parse_fragment_evidence_v0.1.py"
AUDIT="$PKG/workflow/audit_target_candidates_v0.1.py"
RISKER="$PKG/workflow/annotate_fragment_risk_v0.1.py"
FILTERPY="$PKG/workflow/build_primary_strict_qnames_v0.2.py"

POST="$PKG/workflow/run_posttrim_portable.sh"
BREADTH="$PKG/workflow/run_dedup_dispersion_portable.sh"

for F in \
    "$PARSER" \
    "$AUDIT" \
    "$RISKER" \
    "$FILTERPY" \
    "$POST" \
    "$BREADTH"
do
    [ -s "$F" ] || {
        echo "ERROR: required workflow file missing: $F" >&2
        exit 4
    }
done

TRIMDIR="$RUN_ROOT/01_trim"
HOSTDIR="$RUN_ROOT/02_host_depletion"
JOINTDIR="$RUN_ROOT/03_joint_alignment"
EVIDDIR="$RUN_ROOT/04_fragment_evidence"
RISKDIR="$RUN_ROOT/05_risk"
FILTERDIR="$RUN_ROOT/06_filter"
BREADTHDIR="$RUN_ROOT/07_breadth"
LOGDIR="$RUN_ROOT/logs"

mkdir -p \
    "$TRIMDIR" \
    "$HOSTDIR" \
    "$JOINTDIR" \
    "$EVIDDIR" \
    "$RISKDIR" \
    "$FILTERDIR" \
    "$BREADTHDIR" \
    "$LOGDIR"

LOCK="$RUN_ROOT/.sample_running"

if ! mkdir "$LOCK" 2>/dev/null; then
    echo "ERROR: sample already running or stale lock:"
    echo "$LOCK"
    exit 5
fi

trap 'rm -rf "$LOCK"' EXIT

READ_LENGTH_PARAMETERS="$RUN_ROOT/read_length_parameters.tsv"
EXPECTED_READ_LENGTH_PARAMETERS=$(printf 'key\tvalue\nnominal_pe_length\t%s\nnominal_pe_length_source\t%s\ncore_min_aligned_bp\t%s\neffective_min_aligned_bp\t%s\n' \
    "$NOMINAL_PE_LENGTH" "$NOMINAL_PE_LENGTH_SOURCE" "$CORE_MIN_ALIGNED_BP" "$EFFECTIVE_MIN_ALIGNED_BP")

if [ -s "$READ_LENGTH_PARAMETERS" ]; then
    [ "$(cat "$READ_LENGTH_PARAMETERS")" = "$EXPECTED_READ_LENGTH_PARAMETERS" ] || {
        echo "ERROR: existing run-root has different PE read-length parameters; use a new run-root" >&2
        exit 6
    }
elif [ -s "$EVIDDIR/$SAMPLE.target_candidates.audit.tsv" ] || [ -s "$RISKDIR/$SAMPLE.fragment_risk_annotated.tsv" ]; then
    echo "ERROR: existing threshold-dependent outputs have no read-length parameter record; use a new run-root" >&2
    exit 6
else
    printf '%s\n' "$EXPECTED_READ_LENGTH_PARAMETERS" > "$READ_LENGTH_PARAMETERS.tmp"
    mv "$READ_LENGTH_PARAMETERS.tmp" "$READ_LENGTH_PARAMETERS"
fi

echo "============================================================"
echo "MycoWES v0.2 portable single-sample engine"
echo "============================================================"
echo "SAMPLE=$SAMPLE"
echo "R1=$R1"
echo "R2=$R2"
echo "RUN_ROOT=$RUN_ROOT"
echo "CONFIG=$CONFIG"
echo "NOMINAL_PE_LENGTH=$NOMINAL_PE_LENGTH"
echo "NOMINAL_PE_LENGTH_SOURCE=$NOMINAL_PE_LENGTH_SOURCE"
echo "EFFECTIVE_MIN_ALIGNED_BP=$EFFECTIVE_MIN_ALIGNED_BP"
echo "START=$(date --iso-8601=seconds)"
echo "============================================================"

# ============================================================
# STEP 1: adapter trimming
# ============================================================

echo
echo "===== STEP 1: CUTADAPT ====="

TRIM_R1="$TRIMDIR/$SAMPLE.trimmed_R1.fastq.gz"
TRIM_R2="$TRIMDIR/$SAMPLE.trimmed_R2.fastq.gz"

if [ "$HOST_DEPLETED" -eq 1 ]; then
    echo "INPUT_STAGE=GRCh38_PAIRED_BOTH_UNMAPPED -> SKIP CUTADAPT"
elif (
    [ -s "$TRIM_R1" ] &&
    [ -s "$TRIM_R2" ] &&
    gzip -t "$TRIM_R1" 2>/dev/null &&
    gzip -t "$TRIM_R2" 2>/dev/null
)
then
    echo "TRIM already complete -> SKIP"
else
    TMP_R1="$TRIMDIR/$SAMPLE.tmp.trimmed_R1.fastq.gz"
    TMP_R2="$TRIMDIR/$SAMPLE.tmp.trimmed_R2.fastq.gz"

    rm -f "$TMP_R1" "$TMP_R2"

    "$CUTADAPT_BIN" \
        -j "$THREADS" \
        -e 0.10 \
        -O 5 \
        -n 2 \
        -m 50 \
        --pair-filter=any \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        -o "$TMP_R1" \
        -p "$TMP_R2" \
        "$R1" \
        "$R2" \
        > "$LOGDIR/$SAMPLE.cutadapt.log" \
        2>&1

    gzip -t "$TMP_R1"
    gzip -t "$TMP_R2"

    mv "$TMP_R1" "$TRIM_R1"
    mv "$TMP_R2" "$TRIM_R2"

    echo "TRIM COMPLETE"
fi

# ============================================================
# STEP 2: human depletion + joint alignment
# ============================================================

echo
echo "===== STEP 2: HUMAN DEPLETION + JOINT ALIGNMENT ====="

if [ "$HOST_DEPLETED" -eq 1 ]; then
    bash "$POST" "$SAMPLE" "$RUN_ROOT" "$CONFIG" "$R1" "$R2"
else
    bash "$POST" "$SAMPLE" "$RUN_ROOT" "$CONFIG"
fi

JOINT_BAM="$JOINTDIR/$SAMPLE.trimmed.AB_CDE_joint.qname.bam"

UNMAP_R1="$HOSTDIR/$SAMPLE.trimmed.hg38_both_unmapped_R1.fastq.gz"
UNMAP_R2="$HOSTDIR/$SAMPLE.trimmed.hg38_both_unmapped_R2.fastq.gz"

"$SAMTOOLS_BIN" quickcheck "$JOINT_BAM"

# ============================================================
# STEP 3: human-unmapped count
# ============================================================

echo
echo "===== STEP 3: HUMAN-UNMAPPED COUNT ====="

EXPECTED=$(gzip -cd "$UNMAP_R1" | awk 'END {print int(NR/4)}')
EXPECTED2=$(gzip -cd "$UNMAP_R2" | awk 'END {print int(NR/4)}')

echo "human_unmapped_R1=$EXPECTED"
echo "human_unmapped_R2=$EXPECTED2"

[ "$EXPECTED" -eq "$EXPECTED2" ] || {
    echo "ERROR: paired FASTQ count mismatch" >&2
    exit 31
}

# ============================================================
# STEP 4: fragment evidence
# ============================================================

echo
echo "===== STEP 4: FRAGMENT EVIDENCE ====="

FRAG="$EVIDDIR/$SAMPLE.fragment_evidence.raw.tsv"
FRAG_SUM="$EVIDDIR/$SAMPLE.fragment_assignment_summary.tsv"

if [ -s "$FRAG" ]; then
    echo "FRAGMENT TABLE exists -> SKIP"
else
    "$PYTHON_BIN" "$PARSER" \
        --sample "$SAMPLE" \
        --bam "$JOINT_BAM" \
        --taxonomy "$TAXONOMY" \
        --output "$FRAG" \
        --summary "$FRAG_SUM" \
        --expected "$EXPECTED" \
        --samtools "$SAMTOOLS_BIN"
fi

# ============================================================
# STEP 5: target audit
# ============================================================

echo
echo "===== STEP 5: TARGET AUDIT ====="

AUDIT_PREFIX="$EVIDDIR/$SAMPLE"
AUDIT_TABLE="$AUDIT_PREFIX.target_candidates.audit.tsv"

ADAPTER_AUDIT_ARGS=()
adapter_audit_values=(
    "${ADAPTER_REF_DIR:-}"
    "${ADAPTER_REF_MIN_OVERLAP:-}"
    "${ADAPTER_REF_MIN_IDENTITY:-}"
)
if [ -n "${ADAPTER_REF_DIR:-}${ADAPTER_REF_MIN_OVERLAP:-}${ADAPTER_REF_MIN_IDENTITY:-}" ]; then
    for value in "${adapter_audit_values[@]}"; do
        [ -n "$value" ] || {
            echo "ERROR: adapter reference audit requires ADAPTER_REF_DIR, ADAPTER_REF_MIN_OVERLAP, and ADAPTER_REF_MIN_IDENTITY" >&2
            exit 4
        }
    done
    ADAPTER_AUDIT_ARGS=(
        --adapter-ref-dir "$ADAPTER_REF_DIR"
        --adapter-min-overlap "$ADAPTER_REF_MIN_OVERLAP"
        --adapter-min-identity "$ADAPTER_REF_MIN_IDENTITY"
    )
fi

if [ -s "$AUDIT_TABLE" ]; then
    echo "TARGET AUDIT exists -> SKIP"
else
    "$PYTHON_BIN" "$AUDIT" \
        --sample "$SAMPLE" \
        --table "$FRAG" \
        --bam "$JOINT_BAM" \
        --output-prefix "$AUDIT_PREFIX" \
        --samtools "$SAMTOOLS_BIN" \
        --min-aligned-bp "$EFFECTIVE_MIN_ALIGNED_BP" \
        "${ADAPTER_AUDIT_ARGS[@]}"
fi

# ============================================================
# STEP 6: reference risk annotation
# ============================================================

echo
echo "===== STEP 6: RISK ANNOTATION ====="

RISK_PREFIX="$RISKDIR/$SAMPLE"
RISK_OUT="$RISK_PREFIX.fragment_risk_annotated.tsv"

if [ -s "$RISK_OUT" ]; then
    echo "RISK annotation exists -> SKIP"
else
    "$PYTHON_BIN" "$RISKER" \
        --sample "$SAMPLE" \
        --audit-table "$AUDIT_TABLE" \
        --risk-bed "$RISK_BED" \
        --min-aligned-bp "$EFFECTIVE_MIN_ALIGNED_BP" \
        --min-mapq "$CORE_MIN_MAPQ" \
        --max-nm-rate "$CORE_MAX_NM_RATE" \
        --output-prefix "$RISK_PREFIX"
fi

# ============================================================
# STEP 7: PRIMARY / STRICT
# ============================================================

echo
echo "===== STEP 7: PRIMARY / STRICT FILTER ====="

"$PYTHON_BIN" "$FILTERPY" \
    "$RISK_OUT" \
    "$GLOBAL_CAPTURE_DF2" \
    "$FILTERDIR"

PRIMARY_Q="$FILTERDIR/PRIMARY.keep_qnames.txt"
STRICT_Q="$FILTERDIR/STRICT_DF2.keep_qnames.txt"

[ -f "$PRIMARY_Q" ] || {
    echo "ERROR: PRIMARY qname file missing" >&2
    exit 71
}

[ -f "$STRICT_Q" ] || {
    echo "ERROR: STRICT qname file missing" >&2
    exit 72
}

# ============================================================
# STEP 8: PRIMARY breadth
# ============================================================

echo
echo "===== STEP 8: PRIMARY DEDUP + DISPERSION ====="

PN=$(wc -l < "$PRIMARY_Q")
echo "PRIMARY_FRAGMENT_COUNT=$PN"

PRIMARY_OUT="$BREADTHDIR/PRIMARY"

if [ "$PN" -gt 0 ]; then
    rm -rf "$PRIMARY_OUT"

    bash "$BREADTH" \
        "$SAMPLE" \
        "$JOINT_BAM" \
        "$FRAG" \
        "$PRIMARY_Q" \
        "$PRIMARY_OUT" \
        "$CONFIG"
else
    echo "PRIMARY=0 -> NO BREADTH RUN"
    rm -rf "$PRIMARY_OUT"
fi

# ============================================================
# STEP 9: STRICT breadth
# ============================================================

echo
echo "===== STEP 9: STRICT DEDUP + DISPERSION ====="

SN=$(wc -l < "$STRICT_Q")
echo "STRICT_FRAGMENT_COUNT=$SN"

STRICT_OUT="$BREADTHDIR/STRICT_DF2"

if [ "$SN" -gt 0 ]; then
    rm -rf "$STRICT_OUT"

    bash "$BREADTH" \
        "$SAMPLE" \
        "$JOINT_BAM" \
        "$FRAG" \
        "$STRICT_Q" \
        "$STRICT_OUT" \
        "$CONFIG"
else
    echo "STRICT=0 -> NO BREADTH RUN"
    rm -rf "$STRICT_OUT"
fi

# ============================================================
# run manifest
# ============================================================

cat > "$RUN_ROOT/run_manifest.tsv" <<EOF
key	value
software_version	$SOFTWARE_VERSION
reference_release	$REFERENCE_RELEASE
reference_profile	$REFERENCE_PROFILE
excluded_reference_assemblies	${EXCLUDED_REFERENCE_ASSEMBLIES:-NONE}
reference_sha256	$REFERENCE_SHA256
sample_id	$SAMPLE
input_stage	$([ "$HOST_DEPLETED" -eq 1 ] && echo GRCh38_PAIRED_BOTH_UNMAPPED || echo RAW_FASTQ)
r1	$R1
r2	$R2
config	$CONFIG
nominal_pe_length	$NOMINAL_PE_LENGTH
nominal_pe_length_source	$NOMINAL_PE_LENGTH_SOURCE
core_min_aligned_bp	$CORE_MIN_ALIGNED_BP
effective_min_aligned_bp	$EFFECTIVE_MIN_ALIGNED_BP
cutadapt	$CUTADAPT_BIN
minimap2	$MINIMAP2_BIN
samtools	$SAMTOOLS_BIN
python	$PYTHON_BIN
hg38_mmi	$HG38_MMI
joint_mmi	$JOINT_MMI
target_reference	$AB_TARGET_FASTA
taxonomy	$TAXONOMY
risk_bed	$RISK_BED
global_capture_df2	$GLOBAL_CAPTURE_DF2
status	COMPLETE
EOF

echo
echo "============================================================"
echo "MYCOWES_SINGLE_SAMPLE_COMPLETE=$SAMPLE"
echo "FINISH=$(date --iso-8601=seconds)"
echo "OUTPUT=$RUN_ROOT"
echo "============================================================"
