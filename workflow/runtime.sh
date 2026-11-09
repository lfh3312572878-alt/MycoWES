#!/usr/bin/env bash

# ============================================================
# MycoWES v0.3 shared runtime
# ============================================================

MYCOWES_PKG="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.." &&
    pwd
)"

MYCOWES_CONFIG="${MYCOWES_CONFIG:-$MYCOWES_PKG/config/default.env}"

MYCOWES_REF_ROOT="${MYCOWES_REF_ROOT:-$MYCOWES_PKG/reference}"

MYCOWES_MICROBIAL_REF_ROOT="${MYCOWES_MICROBIAL_REF_ROOT:-$MYCOWES_REF_ROOT/microbial}"

MYCOWES_HOST_REF_ROOT="${MYCOWES_HOST_REF_ROOT:-$MYCOWES_REF_ROOT/host}"

MYCOWES_CAPTURE_ROOT="${MYCOWES_CAPTURE_ROOT:-$MYCOWES_REF_ROOT/capture}"

MYCOWES_INDEX_ROOT="${MYCOWES_INDEX_ROOT:-$MYCOWES_PKG/cache/indexes}"

[ -f "$MYCOWES_CONFIG" ] || {
    echo "ERROR: config missing: $MYCOWES_CONFIG" >&2
    return 1 2>/dev/null || exit 1
}

# shellcheck disable=SC1090
source "$MYCOWES_CONFIG"


# ============================================================
# Programs
# ============================================================

resolve_program() {

    local name="$1"
    local configured="$2"

    if [[ "$configured" == */* ]]; then

        [ -x "$configured" ] || {
            echo "ERROR: $name executable not found: $configured" >&2
            return 1
        }

        readlink -f "$configured"
        return 0
    fi

    command -v "$configured" 2>/dev/null || {
        echo "ERROR: cannot find $name in PATH: $configured" >&2
        return 1
    }
}

CUTADAPT_BIN="$(resolve_program cutadapt "$CUTADAPT")" || exit 1
MINIMAP2_BIN="$(resolve_program minimap2 "$MINIMAP2")" || exit 1
SAMTOOLS_BIN="$(resolve_program samtools "$SAMTOOLS")" || exit 1
PYTHON_BIN="$(resolve_program python "$PYTHON")" || exit 1


# ============================================================
# Microbial reference
# ============================================================

JOINT_FASTA="${MYCOWES_JOINT_FASTA:-$MYCOWES_MICROBIAL_REF_ROOT/AB_CDE_JOINT.fa}"
if [ ! -s "$JOINT_FASTA" ] && [ -s "$MYCOWES_MICROBIAL_REF_ROOT/AB_CDE_JOINT.fa.gz" ]; then
    JOINT_FASTA="$MYCOWES_MICROBIAL_REF_ROOT/AB_CDE_JOINT.fa.gz"
fi

JOINT_MMI="${MYCOWES_JOINT_MMI:-$MYCOWES_INDEX_ROOT/microbial/AB_CDE_JOINT.sr.mmi}"

TAXONOMY="${MYCOWES_TAXONOMY:-$MYCOWES_MICROBIAL_REF_ROOT/AB_CDE_taxonomy_lineage_map.tsv}"

RISK_BED="${MYCOWES_RISK_BED:-$MYCOWES_MICROBIAL_REF_ROOT/AB_reference_risk.bed}"

AB_TARGET_FASTA="${MYCOWES_TARGET_FASTA:-${MYCOWES_AB_TARGET_FASTA:-$MYCOWES_MICROBIAL_REF_ROOT/AB_TARGET.fa}}"
if [ ! -s "$AB_TARGET_FASTA" ] && [ -s "$MYCOWES_MICROBIAL_REF_ROOT/AB_TARGET.fa.gz" ]; then
    AB_TARGET_FASTA="$MYCOWES_MICROBIAL_REF_ROOT/AB_TARGET.fa.gz"
fi

GLOBAL_CAPTURE_DF2="${MYCOWES_GLOBAL_CAPTURE_DF2:-$MYCOWES_CAPTURE_ROOT/global/panWES.design_family_support_ge2.bed}"


# ============================================================
# Host reference
# ============================================================

HOST_REFERENCE_PROFILE="${HOST_REFERENCE_PROFILE:-MycoWES_GRCh38_v1}"

HOST_FASTA="${MYCOWES_HOST_FASTA:-${HOST_FASTA:-}}"
HG38_MMI="${MYCOWES_HOST_MMI:-${HOST_MMI:-}}"


# If user did not explicitly specify host files,
# look for an installed standard host bundle.

if [ -z "$HOST_FASTA" ]; then

    CANDIDATE="$MYCOWES_HOST_REF_ROOT/$HOST_REFERENCE_PROFILE/hg38.fa"

    if [ -s "$CANDIDATE" ]; then
        HOST_FASTA="$CANDIDATE"
    fi

    if [ -z "$HOST_FASTA" ]; then
        CANDIDATE="$MYCOWES_HOST_REF_ROOT/$HOST_REFERENCE_PROFILE/hg38.fa.gz"
        if [ -s "$CANDIDATE" ]; then
            HOST_FASTA="$CANDIDATE"
        fi
    fi
fi


if [ -z "$HG38_MMI" ]; then

    CANDIDATE="$MYCOWES_INDEX_ROOT/host/$HOST_REFERENCE_PROFILE/hg38.mmi"

    if [ -s "$CANDIDATE" ]; then
        HG38_MMI="$CANDIDATE"
    fi
fi


# ============================================================
# Reference identity metadata
# ============================================================

MICROBIAL_IDENTITY="${MYCOWES_MICROBIAL_IDENTITY:-$MYCOWES_MICROBIAL_REF_ROOT/identity.tsv}"
if [ ! -s "$MICROBIAL_IDENTITY" ]; then
    MICROBIAL_IDENTITY="$MYCOWES_REF_ROOT/core_profiles/AB_CDE_v0.25.identity.tsv"
fi

HOST_IDENTITY="$MYCOWES_HOST_REF_ROOT/$HOST_REFERENCE_PROFILE/identity.tsv"
if [ ! -s "$HOST_IDENTITY" ]; then
    HOST_IDENTITY="$MYCOWES_REF_ROOT/host_profiles/$HOST_REFERENCE_PROFILE.identity.tsv"
fi

identity_value() {
    local identity_file="$1"
    local identity_key="$2"
    awk -F '\t' -v key="$identity_key" '$1 == key { print $2; exit }' "$identity_file"
}

[ -s "$MICROBIAL_IDENTITY" ] || {
    echo "ERROR: microbial reference identity metadata is missing." >&2
    return 1 2>/dev/null || exit 1
}

[ "$(identity_value "$MICROBIAL_IDENTITY" profile_id)" = "${MYCOWES_MICROBIAL_PROFILE_ID:-MycoWES_AB_CDE_v0.25}" ] || {
    echo "ERROR: microbial reference profile identity mismatch." >&2
    return 1 2>/dev/null || exit 1
}

SOFTWARE_VERSION="$(tr -d '\r\n' < "$MYCOWES_PKG/VERSION")"
REFERENCE_RELEASE="${MYCOWES_REFERENCE_RELEASE:-$(identity_value "$MICROBIAL_IDENTITY" reference_release)}"
REFERENCE_PROFILE="${MYCOWES_REFERENCE_PROFILE:-$(identity_value "$MICROBIAL_IDENTITY" reference_profile)}"
HISTORICAL_SOURCE_MMI_SHA256="$(identity_value "$MICROBIAL_IDENTITY" joint_mmi_sha256)"
# Actual local index identity is recorded by the managed CLI, not copied from a source snapshot hash.
REFERENCE_SHA256="NOT_COMPUTED_SEE_LOCAL_FILE_IDENTITIES"
REFERENCE_EXCLUSION_MANIFEST="${MYCOWES_REFERENCE_EXCLUSION_MANIFEST:-}"
EXCLUDED_REFERENCE_ASSEMBLIES=""
if [ -s "$REFERENCE_EXCLUSION_MANIFEST" ]; then
    EXCLUDED_REFERENCE_ASSEMBLIES="$(awk -F '\t' 'NR>1 {print $1}' "$REFERENCE_EXCLUSION_MANIFEST" | paste -sd, -)"
fi

HOST_REFERENCE_STATUS=CUSTOM
if [ -s "$HOST_IDENTITY" ]; then
    if [ "$(identity_value "$HOST_IDENTITY" profile_id)" = "MycoWES_GRCh38_v1" ] && \
       [ "$(identity_value "$HOST_IDENTITY" sequence_identity_sha256)" = "283229c8ad1bb1844f9419a9f5667827002405a8d557bf04d3c9615d4c08f66c" ]; then
        HOST_REFERENCE_STATUS=STANDARD_MATCH
    elif [ "$HOST_REFERENCE_PROFILE" = "MycoWES_GRCh38_v1" ]; then
        echo "ERROR: standard host profile identity mismatch." >&2
        return 1 2>/dev/null || exit 1
    fi
fi


# ============================================================
# Required runtime assets
# ============================================================

if [ -z "$HG38_MMI" ]; then
    echo "ERROR: human reference index is not configured." >&2
    echo "Install a host index or set MYCOWES_HOST_MMI." >&2
    return 1 2>/dev/null || exit 1
fi


for F in \
    "$HG38_MMI" \
    "$JOINT_MMI" \
    "$TAXONOMY" \
    "$RISK_BED" \
    "$AB_TARGET_FASTA" \
    "$GLOBAL_CAPTURE_DF2"
do

    [ -s "$F" ] || {
        echo "ERROR: required runtime asset missing: $F" >&2
        return 1 2>/dev/null || exit 1
    }

done
