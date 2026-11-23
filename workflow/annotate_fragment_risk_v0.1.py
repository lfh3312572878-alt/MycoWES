#!/usr/bin/env python3

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def to_int(value):
    if value in ("", None):
        return None
    return int(float(value))


def to_float(value):
    if value in ("", None):
        return None
    return float(value)


def reference_interval_0based(position_1based, cigar):
    """
    SAM POS: 1-based inclusive.
    Return BED-style 0-based half-open interval.
    """

    if position_1based is None or not cigar or cigar == "*":
        return None

    reference_length = 0

    for number, operation in CIGAR_RE.findall(cigar):
        if operation in "MDN=X":
            reference_length += int(number)

    if reference_length == 0:
        return None

    start = position_1based - 1
    end = start + reference_length

    return start, end


def interval_overlap(start1, end1, start2, end2):
    return max(0, min(end1, end2) - max(start1, start2))


def normalize_risk_label(label):
    text = label.lower()

    if "human_homology" in text:
        return "human_homology"

    if "low_complexity" in text:
        return "low_complexity"

    if "rrna" in text:
        return "rRNA"

    if "cross_species" in text:
        return "cross_species"

    if "ambiguous" in text:
        return "reference_ambiguous"

    return label.split("|")[0] or "other_risk"


def load_risk_bed(path):
    risks = defaultdict(list)

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue

            if line.startswith(("#", "track", "browser")):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue

            label = fields[3] if len(fields) >= 4 else "risk"

            risks[fields[0]].append({
                "start": start,
                "end": end,
                "label": label,
                "category": normalize_risk_label(label),
            })

    for reference in risks:
        risks[reference].sort(
            key=lambda record: (record["start"], record["end"])
        )

    return risks


def find_overlaps(reference, start, end, risks):
    hits = []

    for record in risks.get(reference, []):
        if record["start"] >= end:
            break

        if record["end"] <= start:
            continue

        overlap = interval_overlap(
            start,
            end,
            record["start"],
            record["end"],
        )

        if overlap > 0:
            hits.append({
                **record,
                "overlap_bp": overlap,
            })

    return hits


def core_filter_reason(row, min_aligned_bp, min_mapq=20, max_nm=0.08, mapq_policy='native'):
    if row["pair_support_class"] != "paired_target_consistent":
        return "PAIR_SUPPORT_NOT_PAIRED_TARGET_CONSISTENT"

    r1_aligned = to_int(row["R1_primary_aligned_bases"])
    r2_aligned = to_int(row["R2_primary_aligned_bases"])

    r1_mapq = to_int(row["R1_primary_MAPQ"])
    r2_mapq = to_int(row["R2_primary_MAPQ"])

    max_nm_rate = to_float(row["maximum_primary_NM_rate"])

    if r1_aligned is None or r2_aligned is None:
        return "ALIGNED_BP_MISSING"

    if r1_aligned < min_aligned_bp or r2_aligned < min_aligned_bp:
        return "ALIGNED_BP_LT_EFFECTIVE_MIN"

    if r1_mapq is None or r2_mapq is None:
        return "MAPQ_MISSING"

    species = row.get('R1_primary_species', '')
    same_species_pair = (bool(species) and ';' not in species and
        row.get('R1_primary_class') == row.get('R2_primary_class') == 'TARGET' and
        species == row.get('R2_primary_species') == row.get('R1_best_target_species') == row.get('R2_best_target_species'))
    exception_field = {'assembly-ties-only': '_assembly_tie_only',
                       'assembly-top-ties-gap10': '_assembly_top_tie_gap10'}.get(mapq_policy)
    if any(q < min_mapq and not (exception_field and same_species_pair and
           str(row.get(end+exception_field, '0')) == '1')
           for end, q in [('R1', r1_mapq), ('R2', r2_mapq)]):
        return f"MAPQ_LT{min_mapq}"

    if max_nm_rate is None or max_nm_rate > max_nm:
        return f"NM_RATE_MISSING_OR_GT{max_nm}"

    return "PASS"


def species_label(row):
    label = row.get("target_species_label", "")

    if label:
        return label

    r1 = row.get("R1_best_target_species", "")
    r2 = row.get("R2_best_target_species", "")

    if r1 and r2 and r1 == r2:
        return r1

    return r1 or r2 or "unresolved_target"


def exact_key(row):
    return (
        row.get("exact_pair_sequence_sha256", "")
        or row["fragment_id"]
    )


def deduplicated_count(rows):
    return len({
        exact_key(row)
        for row in rows
    })


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample", required=True)
    parser.add_argument("--audit-table", required=True)
    parser.add_argument("--risk-bed", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--min-aligned-bp", type=int, default=100)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--max-nm-rate", type=float, default=0.08)
    parser.add_argument('--mapq-policy', choices=('native', 'assembly-ties-only', 'assembly-top-ties-gap10'), default='native')

    args = parser.parse_args()

    if args.min_aligned_bp < 1:
        parser.error("--min-aligned-bp must be positive")
    if not 0 <= args.min_mapq <= 255 or not 0 <= args.max_nm_rate <= 1:
        parser.error("quality thresholds out of range")

    risks = load_risk_bed(args.risk_bed)

    with open(args.audit_table, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    annotated_rows = []
    risk_category_counts = Counter()
    core_reason_counts = Counter()

    for row in rows:
        all_hits = []
        end_results = {}

        for read_end in ("R1", "R2"):
            reference = row.get(
                f"{read_end}_primary_reference",
                "",
            )

            primary_class = row.get(
                f"{read_end}_primary_class",
                "",
            )

            position = to_int(
                row.get(f"{read_end}_primary_position")
            )

            cigar = row.get(
                f"{read_end}_primary_CIGAR",
                "",
            )

            hits = []

            if primary_class == "TARGET" and reference:
                interval = reference_interval_0based(
                    position,
                    cigar,
                )

                if interval is not None:
                    start, end = interval

                    hits = find_overlaps(
                        reference,
                        start,
                        end,
                        risks,
                    )

                    end_results[read_end] = {
                        "start": start,
                        "end": end,
                        "hits": hits,
                    }
                else:
                    end_results[read_end] = {
                        "start": "",
                        "end": "",
                        "hits": [],
                    }
            else:
                end_results[read_end] = {
                    "start": "",
                    "end": "",
                    "hits": [],
                }

            all_hits.extend(hits)

        categories = sorted({
            hit["category"]
            for hit in all_hits
        })

        labels = sorted({
            hit["label"]
            for hit in all_hits
        })

        for category in categories:
            risk_category_counts[category] += 1

        annotated = dict(row)

        recorded_minimum = row.get("effective_min_aligned_bp", "")
        if recorded_minimum and int(recorded_minimum) != args.min_aligned_bp:
            raise SystemExit(
                "audit-table effective_min_aligned_bp does not match --min-aligned-bp"
            )

        reason = core_filter_reason(row, args.min_aligned_bp, args.min_mapq, args.max_nm_rate, args.mapq_policy)
        core_reason_counts[reason] += 1

        annotated.update({
            "effective_min_aligned_bp": args.min_aligned_bp,
            "exploratory_core_alignment": int(reason == "PASS"),
            "core_filter_reason_code": reason,
            "risk_any": int(bool(all_hits)),
            "risk_categories": ";".join(categories),
            "risk_labels": ";".join(labels),
            "risk_overlap_record_count": len(all_hits),
            "risk_overlap_bp_sum_across_read_alignments": sum(
                hit["overlap_bp"]
                for hit in all_hits
            ),

            "R1_reference_start_0based": (
                end_results["R1"]["start"]
            ),
            "R1_reference_end_0based_exclusive": (
                end_results["R1"]["end"]
            ),
            "R1_risk_any": int(
                bool(end_results["R1"]["hits"])
            ),

            "R2_reference_start_0based": (
                end_results["R2"]["start"]
            ),
            "R2_reference_end_0based_exclusive": (
                end_results["R2"]["end"]
            ),
            "R2_risk_any": int(
                bool(end_results["R2"]["hits"])
            ),
        })

        annotated_rows.append(annotated)

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    details_path = Path(
        f"{args.output_prefix}.fragment_risk_annotated.tsv"
    )

    with open(
        details_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        annotated_fieldnames = input_fieldnames + [
            *([] if "effective_min_aligned_bp" in input_fieldnames else ["effective_min_aligned_bp"]),
            "exploratory_core_alignment", "core_filter_reason_code",
            "risk_any", "risk_categories",
            "risk_labels", "risk_overlap_record_count",
            "risk_overlap_bp_sum_across_read_alignments",
            "R1_reference_start_0based", "R1_reference_end_0based_exclusive",
            "R1_risk_any", "R2_reference_start_0based",
            "R2_reference_end_0based_exclusive", "R2_risk_any",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=annotated_fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(annotated_rows)

    all_core = [
        row
        for row in annotated_rows
        if int(row["exploratory_core_alignment"]) == 1
    ]

    risk_core = [
        row
        for row in all_core
        if int(row["risk_any"]) == 1
    ]

    risk_free_core = [
        row
        for row in all_core
        if int(row["risk_any"]) == 0
    ]

    summary_path = Path(
        f"{args.output_prefix}.fragment_risk_summary.tsv"
    )

    with open(
        summary_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "raw", "exact_sequence_deduplicated"])

        writer.writerow([
            "effective_min_aligned_bp",
            args.min_aligned_bp,
            "NA",
        ])

        def write_metric(label, selected_rows):
            writer.writerow([
                label,
                len(selected_rows),
                deduplicated_count(selected_rows),
            ])

        write_metric(
            "all_target_candidates",
            annotated_rows,
        )

        for reason, count in sorted(core_reason_counts.items()):
            writer.writerow([f"core_filter_reason:{reason}", count, "NA"])

        write_metric(
            "target_candidates_with_any_risk",
            [
                row for row in annotated_rows
                if int(row["risk_any"]) == 1
            ],
        )

        write_metric(
            "target_candidates_risk_free",
            [
                row for row in annotated_rows
                if int(row["risk_any"]) == 0
            ],
        )

        write_metric(
            "exploratory_core_all",
            all_core,
        )

        write_metric(
            "exploratory_core_with_any_risk",
            risk_core,
        )

        write_metric(
            "exploratory_core_risk_free",
            risk_free_core,
        )

    category_path = Path(
        f"{args.output_prefix}.fragment_risk_category_counts.tsv"
    )

    with open(
        category_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["risk_category", "candidate_fragments"])

        for category, count in sorted(
            risk_category_counts.items()
        ):
            writer.writerow([category, count])

    print(f"summary={summary_path}")
    print(f"categories={category_path}")
    print(f"details={details_path}")


if __name__ == "__main__":
    main()
