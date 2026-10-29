#!/usr/bin/env bash
set -euo pipefail

PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SAMPLE="${1:?Usage: run_dedup_dispersion_portable.sh SAMPLE BAM TABLE KEEP_QNAMES OUT [CONFIG]}"
BAM="${2:?missing BAM}"
TABLE="${3:?missing fragment table}"
KEEP_QNAMES="${4:?missing keep-qnames file}"
OUT="${5:?missing output directory}"
CONFIG="${6:-$PKG/config/default.env}"

export MYCOWES_CONFIG="$CONFIG"
source "$PKG/workflow/runtime.sh"

MM2="$MINIMAP2_BIN"
SAMTOOLS="$SAMTOOLS_BIN"
AB="$AB_TARGET_FASTA"

mkdir -p "$OUT" "$OUT/logs"

[ -f "$BAM" ] || {
    echo "ERROR: missing BAM: $BAM" >&2
    exit 1
}

[ -f "$TABLE" ] || {
    echo "ERROR: missing fragment table: $TABLE" >&2
    exit 1
}

[ -f "$KEEP_QNAMES" ] || {
    echo "ERROR: missing keep-qnames: $KEEP_QNAMES" >&2
    exit 1
}

[ -f "$AB" ] || {
    echo "ERROR: missing AB FASTA: $AB" >&2
    exit 1
}

LOCK="$OUT/.running"

if ! mkdir "$LOCK" 2>/dev/null; then
    echo "ERROR: sample already running or stale lock:"
    echo "$LOCK"
    exit 1
fi

TMP=$(mktemp -d "$OUT/tmp.XXXXXX")

cleanup() {
    rm -rf "$TMP"
    rm -rf "$LOCK"
}

trap cleanup EXIT

die () {
    echo "ERROR: $1" >&2
    exit 1
}

echo "=================================================="
echo "MycoWES portable dedup + dispersion"
echo "SAMPLE=$SAMPLE"
echo "OUT=$OUT"
date
echo "=================================================="

# ------------------------------------------------------------
# 1. candidate fragment IDs
# ------------------------------------------------------------

QN="$OUT/00_input_fragment_qnames.txt"

if [ -n "$KEEP_QNAMES" ]; then

    [ -f "$KEEP_QNAMES" ] || die "KEEP_QNAMES file missing: $KEEP_QNAMES"

    sort -u "$KEEP_QNAMES" > "$QN"

    echo "Input mode: external filtered qnames"
    echo "KEEP_QNAMES=$KEEP_QNAMES"

else

"$PYTHON_BIN" - "$TABLE" "$QN" <<'PY'
import csv,sys

table,out=sys.argv[1:]

with open(table) as f:

    r=csv.DictReader(f,delimiter="\t")
    fields=r.fieldnames or []

    qcol=next(
        (x for x in [
            "fragment_id",
            "qname",
            "query_name",
            "read_name",
            "name"
        ] if x in fields),
        None
    )

    acol=next(
        (x for x in [
            "initial_assignment",
            "assignment",
            "final_assignment"
        ] if x in fields),
        None
    )

    if qcol is None or acol is None:
        raise SystemExit(
            "Cannot find qname/assignment columns: "
            + ",".join(fields)
        )

    q=set()

    for x in r:

        if (x.get(acol) or "").strip()=="target_candidate":

            name=(x.get(qcol) or "").strip()

            if name:
                q.add(name)

with open(out,"w") as o:

    for x in sorted(q):
        o.write(x+"\n")

print("target_candidate_fragments =",len(q))
PY

fi

N=$(wc -l < "$QN")

echo "input_fragments=$N"

[ "$N" -gt 0 ] || die "zero input fragments"


# ------------------------------------------------------------
# 2. extract original paired FASTQ
# ------------------------------------------------------------

SUB="$TMP/subset.qname.bam"

"$SAMTOOLS" view \
    -u \
    -N "$QN" \
    "$BAM" \
| "$SAMTOOLS" view \
    -u \
    -F 2304 \
    - \
| "$SAMTOOLS" sort \
    -n \
    -@ "$SAMTOOLS_THREADS" \
    -o "$SUB" \
    -

RC=("${PIPESTATUS[@]}")

if [ "${RC[0]}" -ne 0 ] || \
   [ "${RC[1]}" -ne 0 ] || \
   [ "${RC[2]}" -ne 0 ]; then
    die "BAM extraction failed"
fi

"$SAMTOOLS" quickcheck "$SUB" || die "subset BAM quickcheck failed"

R1="$TMP/raw.R1.fastq"
R2="$TMP/raw.R2.fastq"

"$SAMTOOLS" fastq \
    -@ "$SAMTOOLS_THREADS" \
    -n \
    -1 "$R1" \
    -2 "$R2" \
    -0 /dev/null \
    -s "$TMP/singleton.fastq" \
    "$SUB" \
    > "$OUT/logs/samtools_fastq.stdout.log" \
    2> "$OUT/logs/samtools_fastq.stderr.log"

[ "$?" -eq 0 ] || die "samtools fastq failed"


# ------------------------------------------------------------
# 3. exact sequence dedup
# ------------------------------------------------------------

DR1="$OUT/01_exact_dedup.R1.fastq"
DR2="$OUT/01_exact_dedup.R2.fastq"
GROUP="$OUT/01_exact_dedup.groups.tsv"

"$PYTHON_BIN" \
"$PKG/workflow/exact_pair_dedup_v0.2.py" \
"$R1" \
"$R2" \
"$DR1" \
"$DR2" \
"$GROUP" \
| tee "$OUT/01_exact_dedup.summary.txt"

[ "$?" -eq 0 ] || die "exact dedup failed"

REP="$OUT/01_exact_dedup.representatives.txt"

awk 'NR>1 {print $1}' "$GROUP" \
| sort -u \
> "$REP"


# ------------------------------------------------------------
# 4. evidence-supported assembly + coordinate dedup + dispersion
# ------------------------------------------------------------

"$PYTHON_BIN" \
"$PKG/workflow/evidence_anchor_breadth_v0.2.py" \
--bam "$BAM" \
--ab "$AB" \
--r1 "$DR1" \
--r2 "$DR2" \
--groups "$GROUP" \
--reps "$REP" \
--out "$OUT" \
--minimap2 "$MM2" \
--samtools "$SAMTOOLS" \
--retention 0.90 \
| tee "$OUT/07_console_summary.txt"

[ "$?" -eq 0 ] || die "evidence-supported breadth failed"


# ------------------------------------------------------------
# 5. method note
# ------------------------------------------------------------

cat > "$OUT/METHOD_v0.2.txt" <<'EOF'
MycoWES secondary dedup + genomic dispersion v0.2

1. Evidence unit:
   paired-end fragment.

2. Exact dedup:
   identical R1 + identical R2 sequences are collapsed.

3. Species assignment:
   determined from AB target alignments in the original
   joint AB+CDE alignment.

4. Breadth/dispersion assembly:
   all AB assemblies belonging to that species are tested.

5. Assembly selection:
   highest paired-fragment retention;
   tie-break by total alignment score,
   fewer contigs, larger assembly, accession.

6. Anchor QC:
   retained exact-dedup fragments / species exact-dedup fragments
   must be >= 0.90.

   If <0.90:
   breadth is not called and status = ANCHOR_FAIL.

7. Coordinate dedup:
   same selected assembly
   + same contig
   + same fragment start
   + same fragment end
   + same paired orientation
   => one coordinate family.

8. Dispersion metrics:
   coordinate-unique fragments
   covered bp
   occupied 250-bp bins
   largest-bin fraction
   distinct contigs.

9. Development tiers:

   supportive:
   >=3 coordinate-unique fragments
   AND >=2 occupied 250-bp bins
   AND largest-bin fraction <0.80

   strong:
   >=5 coordinate-unique fragments
   AND >=3 occupied 250-bp bins
   AND covered bp >=300
   AND largest-bin fraction <=0.50

These are development dispersion tiers,
NOT final sample positivity thresholds.
EOF


echo
echo "=================================================="
echo "FINAL SUMMARY"
echo "=================================================="

column -t -s $'\t' \
"$OUT/07_species_dedup_dispersion_summary.tsv"

echo
echo "OUTPUT=$OUT"
echo "===== COMPLETE ====="
date
