#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import re
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True)
    parser.add_argument("--window-map", required=True)
    parser.add_argument("--read-length", type=int, required=True)
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--near-full-fraction", type=float, default=0.925)
    parser.add_argument("--high-identity", type=float, default=0.95)
    parser.add_argument("--moderate-identity", type=float, default=0.90)
    return parser.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_tags(fields):
    tags = {}

    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]

    return tags


def cigar_metrics(cigar):
    operations = [
        (int(length), operation)
        for length, operation in CIGAR_RE.findall(cigar)
    ]

    query_aligned = sum(
        length
        for length, operation in operations
        if operation in {"M", "I", "=", "X"}
    )

    reference_aligned = sum(
        length
        for length, operation in operations
        if operation in {"M", "D", "N", "=", "X"}
    )

    softclip = sum(
        length
        for length, operation in operations
        if operation == "S"
    )

    return query_aligned, reference_aligned, softclip


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = f"AB_TARGET.{args.profile_id}.human_homology"

    details_out = outdir / f"{prefix}.details.tsv"
    bed_out = outdir / f"{prefix}.risk.bed"
    summary_out = outdir / f"{prefix}.summary.tsv"

    near_full_bp = math.ceil(
        args.read_length * args.near_full_fraction
    )

    candidate_hits = {}
    mapq_bins = Counter()

    mapped_records = 0
    duplicate_primary_qnames = 0

    with open_text(args.hits) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                continue

            qname = fields[0]
            human_reference = fields[2]
            human_position = fields[3]
            mapq = int(fields[4])
            cigar = fields[5]

            tags = parse_tags(fields)

            nm = int(tags["NM"]) if "NM" in tags else -1
            alignment_score = int(tags["AS"]) if "AS" in tags else -1

            query_aligned, reference_aligned, softclip = (
                cigar_metrics(cigar)
            )

            identity_proxy = (
                max(query_aligned - nm, 0) / query_aligned
                if query_aligned > 0 and nm >= 0
                else 0.0
            )

            aligned_fraction = (
                query_aligned / args.read_length
                if args.read_length > 0
                else 0.0
            )

            mapped_records += 1

            if mapq == 0:
                mapq_bins["MAPQ_0"] += 1
            elif mapq < 20:
                mapq_bins["MAPQ_1_19"] += 1
            elif mapq < 30:
                mapq_bins["MAPQ_20_29"] += 1
            else:
                mapq_bins["MAPQ_GE30"] += 1

            if mapq < args.min_mapq:
                continue

            hit = {
                "MAPQ": mapq,
                "NM": nm,
                "AS": alignment_score,
                "human_hit": f"{human_reference}:{human_position}",
                "CIGAR": cigar,
                "query_aligned_bp": query_aligned,
                "reference_aligned_bp": reference_aligned,
                "softclip_bp": softclip,
                "identity_proxy": identity_proxy,
                "aligned_fraction": aligned_fraction,
            }

            old = candidate_hits.get(qname)

            if old is not None:
                duplicate_primary_qnames += 1

                new_key = (
                    alignment_score,
                    query_aligned,
                    -nm,
                    mapq,
                )
                old_key = (
                    old["AS"],
                    old["query_aligned_bp"],
                    -old["NM"],
                    old["MAPQ"],
                )

                if new_key <= old_key:
                    continue

            candidate_hits[qname] = hit

    columns = [
        "window_id",
        "seq_id",
        "start",
        "end",
        "MAPQ",
        "NM",
        "AS",
        "human_hit",
        "CIGAR",
        "query_aligned_bp",
        "reference_aligned_bp",
        "softclip_bp",
        "identity_proxy",
        "aligned_fraction",
        "risk_tier",
    ]

    tier_counts = Counter()
    total_windows = 0
    written = 0

    with open(args.window_map) as map_handle, \
         details_out.open("w", newline="") as details_handle, \
         bed_out.open("w", newline="") as bed_handle:

        reader = csv.DictReader(map_handle, delimiter="\t")

        writer = csv.DictWriter(
            details_handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            total_windows += 1

            window_id = row["window_id"]
            hit = candidate_hits.get(window_id)

            if hit is None:
                continue

            if (
                hit["query_aligned_bp"] >= near_full_bp
                and hit["identity_proxy"] >= args.high_identity
            ):
                risk_tier = "human_homology_high_candidate"

            elif (
                hit["query_aligned_bp"] >= near_full_bp
                and hit["identity_proxy"] >= args.moderate_identity
            ):
                risk_tier = "human_homology_moderate_candidate"

            else:
                risk_tier = "human_homology_local_weak"

            output_row = {
                "window_id": window_id,
                "seq_id": row["seq_id"],
                "start": row["start"],
                "end": row["end"],
                "MAPQ": hit["MAPQ"],
                "NM": hit["NM"],
                "AS": hit["AS"],
                "human_hit": hit["human_hit"],
                "CIGAR": hit["CIGAR"],
                "query_aligned_bp": hit["query_aligned_bp"],
                "reference_aligned_bp": hit["reference_aligned_bp"],
                "softclip_bp": hit["softclip_bp"],
                "identity_proxy": f'{hit["identity_proxy"]:.6f}',
                "aligned_fraction": f'{hit["aligned_fraction"]:.6f}',
                "risk_tier": risk_tier,
            }

            writer.writerow(output_row)

            bed_handle.write(
                "\t".join([
                    row["seq_id"],
                    row["start"],
                    row["end"],
                    f"{risk_tier}|{window_id}",
                    str(hit["MAPQ"]),
                    ".",
                ]) + "\n"
            )

            tier_counts[risk_tier] += 1
            written += 1

    with summary_out.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"profile_id\t{args.profile_id}\n")
        handle.write(f"read_length\t{args.read_length}\n")
        handle.write(
            f"near_full_fraction\t{args.near_full_fraction}\n"
        )
        handle.write(f"near_full_min_bp\t{near_full_bp}\n")
        handle.write(f"min_MAPQ\t{args.min_mapq}\n")
        handle.write(f"high_identity\t{args.high_identity}\n")
        handle.write(
            f"moderate_identity\t{args.moderate_identity}\n"
        )
        handle.write(f"total_windows\t{total_windows}\n")
        handle.write(f"primary_mapped_records\t{mapped_records}\n")

        for name in [
            "MAPQ_0",
            "MAPQ_1_19",
            "MAPQ_20_29",
            "MAPQ_GE30",
        ]:
            handle.write(f"{name}\t{mapq_bins[name]}\n")

        handle.write(f"MAPQ_GE20_candidates\t{written}\n")

        for name in [
            "human_homology_high_candidate",
            "human_homology_moderate_candidate",
            "human_homology_local_weak",
        ]:
            handle.write(f"{name}\t{tier_counts[name]}\n")

        handle.write(
            f"duplicate_primary_qnames\t"
            f"{duplicate_primary_qnames}\n"
        )

    print(f"total_windows={total_windows}")
    print(f"primary_mapped_records={mapped_records}")
    print(f"MAPQ_GE20_candidates={written}")
    print(
        "human_homology_high_candidate="
        f"{tier_counts['human_homology_high_candidate']}"
    )
    print(
        "human_homology_moderate_candidate="
        f"{tier_counts['human_homology_moderate_candidate']}"
    )
    print(
        "human_homology_local_weak="
        f"{tier_counts['human_homology_local_weak']}"
    )
    print(f"near_full_min_bp={near_full_bp}")
    print(f"details={details_out}")
    print(f"risk_bed={bed_out}")
    print(f"summary={summary_out}")


if __name__ == "__main__":
    main()
