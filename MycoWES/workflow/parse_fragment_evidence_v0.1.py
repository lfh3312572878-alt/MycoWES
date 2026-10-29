#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_stats(cigar, pos):
    query_aligned = 0
    reference_aligned = 0
    soft_clipped = 0

    if cigar != "*":
        for number, operation in CIGAR_RE.findall(cigar):
            number = int(number)

            if operation in "MI=X":
                query_aligned += number

            if operation in "MDN=X":
                reference_aligned += number

            if operation == "S":
                soft_clipped += number

    reference_end = (
        pos + reference_aligned - 1
        if reference_aligned
        else 0
    )

    return (
        query_aligned,
        reference_aligned,
        soft_clipped,
        reference_end,
    )


def load_taxonomy(path):
    with open(path, encoding="utf-8") as handle:
        taxonomy = {
            row["reference_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }

    return taxonomy


def parse_record(line, taxonomy):
    fields = line.rstrip().split("\t")

    flag = int(fields[1])
    mapped = not bool(flag & 4)

    tags = {}

    for item in fields[11:]:
        parts = item.split(":", 2)

        if len(parts) == 3:
            tags[parts[0]] = parts[2]

    query_aligned, reference_aligned, soft_clipped, reference_end = (
        cigar_stats(fields[5], int(fields[3]))
    )

    tax = taxonomy.get(fields[2], {}) if mapped else {}

    return {
        "qname": fields[0],
        "flag": flag,
        "read_end": (
            "R1"
            if flag & 64
            else "R2"
            if flag & 128
            else "UNPAIRED"
        ),
        "mapped": mapped,
        "secondary": bool(flag & 256),
        "supplementary": bool(flag & 2048),
        "primary": not bool(flag & 256) and not bool(flag & 2048),
        "proper_pair": bool(flag & 2),
        "reverse": bool(flag & 16),
        "reference": fields[2],
        "position": int(fields[3]),
        "MAPQ": int(fields[4]),
        "CIGAR": fields[5],
        "TLEN": int(fields[8]),
        "AS": int(tags["AS"]) if "AS" in tags else None,
        "NM": int(tags["NM"]) if "NM" in tags else None,
        "query_aligned_bases": query_aligned,
        "reference_aligned_bases": reference_aligned,
        "soft_clipped_bases": soft_clipped,
        "reference_end": reference_end,
        "species": tax.get("reference_species", ""),
        "species_taxid": tax.get("species_taxid", ""),
        "genus": tax.get("genus", ""),
        "genus_taxid": tax.get("genus_taxid", ""),
        "family": tax.get("family", ""),
        "family_taxid": tax.get("family_taxid", ""),
        "order": tax.get("order", ""),
        "order_taxid": tax.get("order_taxid", ""),
        "reference_class": tax.get("reference_class", ""),
        "role": tax.get("role", ""),
    }


def get_primary(records, read_end):
    selected = [
        record
        for record in records
        if record["read_end"] == read_end
        and record["primary"]
    ]

    if len(selected) != 1:
        raise RuntimeError(
            f"{records[0]['qname']}: "
            f"{read_end} primary count={len(selected)}"
        )

    return selected[0]


def collapse_species(records, reference_class):
    best_by_species = {}

    for record in records:
        if not record["mapped"]:
            continue

        if record["supplementary"]:
            continue

        if record["reference_class"] != reference_class:
            continue

        if record["AS"] is None:
            continue

        species_key = (
            record["species_taxid"]
            or record["species"]
        )

        current = best_by_species.get(species_key)

        candidate_rank = (
            record["AS"],
            record["MAPQ"],
            record["query_aligned_bases"],
            -record["soft_clipped_bases"],
        )

        if current is None:
            best_by_species[species_key] = record
            continue

        current_rank = (
            current["AS"],
            current["MAPQ"],
            current["query_aligned_bases"],
            -current["soft_clipped_bases"],
        )

        if candidate_rank > current_rank:
            best_by_species[species_key] = record

    if not best_by_species:
        return None, []

    best_score = max(
        record["AS"]
        for record in best_by_species.values()
    )

    tied_records = sorted(
        [
            record
            for record in best_by_species.values()
            if record["AS"] == best_score
        ],
        key=lambda record: (
            record["species_taxid"],
            record["species"],
        ),
    )

    return best_score, tied_records


def lca_assignment(records):
    if not records:
        return "", "", ""

    ranks = [
        ("species", "species_taxid", "species"),
        ("genus", "genus_taxid", "genus"),
        ("family", "family_taxid", "family"),
        ("order", "order_taxid", "order"),
    ]

    for rank, taxid_field, name_field in ranks:
        taxids = {
            record[taxid_field]
            for record in records
            if record[taxid_field]
        }

        names = {
            record[name_field]
            for record in records
            if record[name_field]
        }

        if len(taxids) == 1 and len(names) == 1:
            return (
                rank,
                next(iter(taxids)),
                next(iter(names)),
            )

    return "ambiguous", "", ""


def summarize_read_end(records, read_end):
    end_records = [
        record
        for record in records
        if record["read_end"] == read_end
    ]

    target_as, target_records = collapse_species(
        end_records,
        "TARGET",
    )

    competitor_as, competitor_records = collapse_species(
        end_records,
        "COMPETITOR",
    )

    if target_as is None and competitor_as is None:
        end_class = "uninformative"

    elif competitor_as is None:
        end_class = "target"

    elif target_as is None:
        end_class = "competitor"

    elif target_as > competitor_as:
        end_class = "target"

    elif competitor_as > target_as:
        end_class = "competitor"

    else:
        end_class = "target_competitor_tie"

    target_rank, target_taxid, target_label = (
        lca_assignment(target_records)
    )

    competitor_rank, competitor_taxid, competitor_label = (
        lca_assignment(competitor_records)
    )

    delta_as = None

    if target_as is not None and competitor_as is not None:
        delta_as = target_as - competitor_as

    return {
        "end_class": end_class,
        "target_AS": target_as,
        "competitor_AS": competitor_as,
        "delta_AS": delta_as,
        "target_species": ";".join(
            sorted({
                record["species"]
                for record in target_records
            })
        ),
        "competitor_species": ";".join(
            sorted({
                record["species"]
                for record in competitor_records
            })
        ),
        "target_rank": target_rank,
        "target_taxid": target_taxid,
        "target_label": target_label,
        "competitor_rank": competitor_rank,
        "competitor_taxid": competitor_taxid,
        "competitor_label": competitor_label,
        "target_records": target_records,
        "competitor_records": competitor_records,
    }


def classify_pair(r1_summary, r2_summary):
    r1_class = r1_summary["end_class"]
    r2_class = r2_summary["end_class"]

    if r1_class == "uninformative" and r2_class == "uninformative":
        return (
            "no_microbial_alignment",
            "no_microbial_alignment",
        )

    if "target_competitor_tie" in (r1_class, r2_class):
        return (
            "fragment_ambiguous",
            "ambiguous",
        )

    if r1_class == "target" and r2_class == "target":
        rank, _, _ = lca_assignment(
            r1_summary["target_records"]
            + r2_summary["target_records"]
        )

        if rank == "ambiguous":
            return (
                "paired_target_discordant",
                "ambiguous",
            )

        return (
            "paired_target_consistent",
            "target_candidate",
        )

    if {r1_class, r2_class} == {"target", "uninformative"}:
        return (
            "single_end_target_candidate",
            "target_candidate",
        )

    if r1_class == "competitor" and r2_class == "competitor":
        rank, _, _ = lca_assignment(
            r1_summary["competitor_records"]
            + r2_summary["competitor_records"]
        )

        if rank == "ambiguous":
            return (
                "paired_competitor_discordant",
                "ambiguous",
            )

        return (
            "paired_competitor_consistent",
            "competitor_not_target",
        )

    if {r1_class, r2_class} == {"competitor", "uninformative"}:
        return (
            "single_end_competitor_candidate",
            "competitor_not_target",
        )

    if {r1_class, r2_class} == {"target", "competitor"}:
        return (
            "target_competitor_discordant",
            "ambiguous",
        )

    return (
        "fragment_ambiguous",
        "ambiguous",
    )


def empty_if_none(value):
    return "" if value is None else value


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--expected", required=True, type=int)
    parser.add_argument(
        "--samtools",
        default="/usr/local/bin/samtools",
    )

    args = parser.parse_args()

    taxonomy = load_taxonomy(args.taxonomy)

    Path(args.output).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    process = subprocess.Popen(
        [args.samtools, "view", args.bam],
        stdout=subprocess.PIPE,
        text=True,
    )

    current_qname = None
    current_records = []

    fragment_rows = []
    summary_counts = Counter()

    def flush_fragment(records):
        if not records:
            return

        r1_primary = get_primary(records, "R1")
        r2_primary = get_primary(records, "R2")

        r1_summary = summarize_read_end(records, "R1")
        r2_summary = summarize_read_end(records, "R2")

        pair_support, initial_assignment = classify_pair(
            r1_summary,
            r2_summary,
        )

        both_primary_mapped = (
            r1_primary["mapped"]
            and r2_primary["mapped"]
        )

        same_contig = int(
            both_primary_mapped
            and r1_primary["reference"]
            == r2_primary["reference"]
        )

        orientation = ""

        if both_primary_mapped:
            orientation = (
                ("R" if r1_primary["reverse"] else "F")
                + ("R" if r2_primary["reverse"] else "F")
            )

        fragment_span = ""

        if same_contig:
            left = min(
                r1_primary["position"],
                r2_primary["position"],
            )

            right = max(
                r1_primary["reference_end"],
                r2_primary["reference_end"],
            )

            fragment_span = right - left + 1

        row = {
            "sample_id": args.sample,
            "fragment_id": records[0]["qname"],
            "alignment_records": len(records),
            "secondary_records": sum(
                record["secondary"]
                for record in records
            ),
            "supplementary_records": sum(
                record["supplementary"]
                for record in records
            ),

            "R1_primary_mapped": int(r1_primary["mapped"]),
            "R1_primary_reference": (
                r1_primary["reference"]
                if r1_primary["mapped"]
                else ""
            ),
            "R1_primary_class": r1_primary["reference_class"],
            "R1_primary_species": r1_primary["species"],
            "R1_primary_AS": empty_if_none(r1_primary["AS"]),
            "R1_primary_MAPQ": (
                r1_primary["MAPQ"]
                if r1_primary["mapped"]
                else ""
            ),
            "R1_primary_NM": empty_if_none(r1_primary["NM"]),
            "R1_primary_CIGAR": (
                r1_primary["CIGAR"]
                if r1_primary["mapped"]
                else ""
            ),
            "R1_primary_aligned_bases": (
                r1_primary["query_aligned_bases"]
                if r1_primary["mapped"]
                else ""
            ),
            "R1_primary_softclip": (
                r1_primary["soft_clipped_bases"]
                if r1_primary["mapped"]
                else ""
            ),

            "R2_primary_mapped": int(r2_primary["mapped"]),
            "R2_primary_reference": (
                r2_primary["reference"]
                if r2_primary["mapped"]
                else ""
            ),
            "R2_primary_class": r2_primary["reference_class"],
            "R2_primary_species": r2_primary["species"],
            "R2_primary_AS": empty_if_none(r2_primary["AS"]),
            "R2_primary_MAPQ": (
                r2_primary["MAPQ"]
                if r2_primary["mapped"]
                else ""
            ),
            "R2_primary_NM": empty_if_none(r2_primary["NM"]),
            "R2_primary_CIGAR": (
                r2_primary["CIGAR"]
                if r2_primary["mapped"]
                else ""
            ),
            "R2_primary_aligned_bases": (
                r2_primary["query_aligned_bases"]
                if r2_primary["mapped"]
                else ""
            ),
            "R2_primary_softclip": (
                r2_primary["soft_clipped_bases"]
                if r2_primary["mapped"]
                else ""
            ),

            "R1_best_target_AS": empty_if_none(
                r1_summary["target_AS"]
            ),
            "R1_best_competitor_AS": empty_if_none(
                r1_summary["competitor_AS"]
            ),
            "R1_delta_AS": empty_if_none(
                r1_summary["delta_AS"]
            ),
            "R1_best_target_species": (
                r1_summary["target_species"]
            ),
            "R1_best_competitor_species": (
                r1_summary["competitor_species"]
            ),
            "R1_end_class": r1_summary["end_class"],

            "R2_best_target_AS": empty_if_none(
                r2_summary["target_AS"]
            ),
            "R2_best_competitor_AS": empty_if_none(
                r2_summary["competitor_AS"]
            ),
            "R2_delta_AS": empty_if_none(
                r2_summary["delta_AS"]
            ),
            "R2_best_target_species": (
                r2_summary["target_species"]
            ),
            "R2_best_competitor_species": (
                r2_summary["competitor_species"]
            ),
            "R2_end_class": r2_summary["end_class"],

            "primary_same_contig": same_contig,
            "primary_orientation": orientation,
            "primary_proper_pair": int(
                r1_primary["proper_pair"]
                and r2_primary["proper_pair"]
            ),
            "primary_fragment_span": fragment_span,
            "pair_support_class": pair_support,
            "initial_assignment": initial_assignment,
        }

        fragment_rows.append(row)

        summary_counts[
            f"assignment:{initial_assignment}"
        ] += 1

        summary_counts[
            f"support:{pair_support}"
        ] += 1

    assert process.stdout is not None

    for line in process.stdout:
        record = parse_record(line, taxonomy)

        if current_qname is None:
            current_qname = record["qname"]

        if record["qname"] != current_qname:
            flush_fragment(current_records)
            current_records = []
            current_qname = record["qname"]

        current_records.append(record)

    flush_fragment(current_records)

    if process.wait() != 0:
        raise SystemExit("samtools view failed")

    if len(fragment_rows) != args.expected:
        raise SystemExit(
            f"Fragment count mismatch: "
            f"{len(fragment_rows)} != {args.expected}"
        )

    with open(
        args.output,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(fragment_rows[0].keys() if fragment_rows else """sample_id	fragment_id	alignment_records	secondary_records	supplementary_records	R1_primary_mapped	R1_primary_reference	R1_primary_class	R1_primary_species	R1_primary_AS	R1_primary_MAPQ	R1_primary_NM	R1_primary_CIGAR	R1_primary_aligned_bases	R1_primary_softclip	R2_primary_mapped	R2_primary_reference	R2_primary_class	R2_primary_species	R2_primary_AS	R2_primary_MAPQ	R2_primary_NM	R2_primary_CIGAR	R2_primary_aligned_bases	R2_primary_softclip	R1_best_target_AS	R1_best_competitor_AS	R1_delta_AS	R1_best_target_species	R1_best_competitor_species	R1_end_class	R2_best_target_AS	R2_best_competitor_AS	R2_delta_AS	R2_best_target_species	R2_best_competitor_species	R2_end_class	primary_same_contig	primary_orientation	primary_proper_pair	primary_fragment_span	pair_support_class	initial_assignment""".split("\t")),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(fragment_rows)

    with open(
        args.summary,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(["metric", "value"])
        writer.writerow(["fragments", len(fragment_rows)])

        for metric, count in sorted(summary_counts.items()):
            writer.writerow([metric, count])

    print(f"fragments={len(fragment_rows)}")
    print(f"fragment_table={args.output}")
    print(f"summary={args.summary}")


if __name__ == "__main__":
    main()
