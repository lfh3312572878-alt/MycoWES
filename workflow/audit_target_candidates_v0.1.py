#!/usr/bin/env python3

import argparse
import bisect
import csv
import hashlib
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ADAPTER_CORE = "AGATCGGAAGAGC"
ADAPTER_SHORT = "AGATCGGA"


def to_int(value):
    if value in ("", None):
        return None
    return int(value)


def cigar_query_aligned(cigar):
    if not cigar or cigar == "*":
        return 0

    total = 0

    for number, operation in CIGAR_RE.findall(cigar):
        if operation in "MI=X":
            total += int(number)

    return total


def reverse_complement(sequence):
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return sequence.translate(table)[::-1]


def load_adapter_references(directory):
    """Load, uppercase, and deduplicate all FASTA sequences in a directory."""
    sequences = defaultdict(list)
    for path in sorted(Path(directory).iterdir()):
        if path.suffix.lower() not in {".fa", ".fasta", ".fna"}:
            continue
        name = None
        chunks = []
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if name and chunks:
                        sequences["".join(chunks).upper()].append(
                            f"{path.name}:{name}"
                        )
                    name, chunks = line[1:].split()[0], []
                elif line:
                    chunks.append(line)
        if name and chunks:
            sequences["".join(chunks).upper()].append(f"{path.name}:{name}")

    oriented = defaultdict(list)
    for sequence, names in sequences.items():
        oriented[sequence].extend(f"{name}:forward" for name in names)
        oriented[reverse_complement(sequence)].extend(
            f"{name}:reverse_complement" for name in names
        )
    return [
        (";".join(sorted(set(names))), sequence)
        for sequence, names in sorted(oriented.items())
    ]


def best_adapter_overlap(read, references, minimum_overlap, minimum_identity):
    """Return the best qualifying ungapped read/reference overlap."""
    read = read.upper()
    best = None
    for name, adapter in references:
        for offset in range(-len(adapter) + 1, len(read)):
            read_start = max(0, offset)
            adapter_start = max(0, -offset)
            diagonal = min(len(read) - read_start, len(adapter) - adapter_start)
            if diagonal < minimum_overlap:
                continue

            match_prefix = [0]
            score_prefix = [0.0]
            for index in range(diagonal):
                matched = read[read_start + index] == adapter[adapter_start + index]
                match_prefix.append(match_prefix[-1] + matched)
                score_prefix.append(score_prefix[-1] + matched - minimum_identity)

            coordinates = sorted(score_prefix)
            earliest = [diagonal + 1] * (len(coordinates) + 1)

            def update(position, value):
                position += 1
                while position < len(earliest):
                    earliest[position] = min(earliest[position], value)
                    position += position & -position

            def query(position):
                result = diagonal + 1
                while position:
                    result = min(result, earliest[position])
                    position -= position & -position
                return result

            for end in range(minimum_overlap, diagonal + 1):
                start_to_add = end - minimum_overlap
                update(bisect.bisect_left(coordinates, score_prefix[start_to_add]), start_to_add)
                rank = bisect.bisect_right(coordinates, score_prefix[end] + 1e-12)
                start = query(rank)
                if start > diagonal:
                    continue
                overlap = end - start
                matches = match_prefix[end] - match_prefix[start]
                identity = matches / overlap
                candidate = (overlap, identity, name, matches, offset)
                if best is None or candidate[:2] > best[:2] or (
                    candidate[:2] == best[:2] and name < best[2]
                ):
                    best = candidate
    return best


def restored_original_sequence(sequence, flag):
    if sequence == "*":
        return ""

    if flag & 16:
        return reverse_complement(sequence)

    return sequence


def quantile(values, probability):
    if not values:
        return ""

    values = sorted(values)

    if len(values) == 1:
        return values[0]

    position = (len(values) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    fraction = position - lower

    return values[lower] + (
        values[upper] - values[lower]
    ) * fraction


def summarize_numeric(writer, prefix, values):
    values = [value for value in values if value is not None]

    if not values:
        return

    writer.writerow([f"{prefix}_n", len(values)])
    writer.writerow([f"{prefix}_min", min(values)])
    writer.writerow([f"{prefix}_q25", quantile(values, 0.25)])
    writer.writerow([f"{prefix}_median", statistics.median(values)])
    writer.writerow([f"{prefix}_q75", quantile(values, 0.75)])
    writer.writerow([f"{prefix}_max", max(values)])


def both_primary_aligned_at_least(rows, minimum):
    return sum(
        1
        for row in rows
        if int(row["R1_primary_aligned_bases"] or 0) >= minimum
        and int(row["R2_primary_aligned_bases"] or 0) >= minimum
    )


def species_label(row):
    r1 = row["R1_best_target_species"]
    r2 = row["R2_best_target_species"]

    if r1 and r2 and r1 == r2:
        return r1

    if r1 and r2:
        return f"{r1} || {r2}"

    return r1 or r2 or "unresolved_target"


def parse_primary_bam(bam, samtools, target_ids):
    records = defaultdict(dict)

    process = subprocess.Popen(
        [samtools, "view", "-F", "2304", bam],
        stdout=subprocess.PIPE,
        text=True,
    )

    assert process.stdout is not None

    for line in process.stdout:
        fields = line.rstrip().split("\t")
        qname = fields[0]

        if qname not in target_ids:
            continue

        flag = int(fields[1])

        if flag & 64:
            read_end = "R1"
        elif flag & 128:
            read_end = "R2"
        else:
            continue

        sequence = restored_original_sequence(
            fields[9],
            flag,
        )

        records[qname][read_end] = {
            "flag": flag,
            "reference": fields[2],
            "position": int(fields[3]),
            "MAPQ": int(fields[4]),
            "CIGAR": fields[5],
            "TLEN": int(fields[8]),
            "sequence": sequence,
            "query_aligned": cigar_query_aligned(fields[5]),
            "strand": "-" if flag & 16 else "+",
        }

    return_code = process.wait()

    if return_code != 0:
        raise SystemExit("samtools view failed")

    return records


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument(
        "--samtools",
        default="samtools",
    )
    parser.add_argument("--adapter-ref-dir")
    parser.add_argument("--adapter-min-overlap", type=int)
    parser.add_argument("--adapter-min-identity", type=float)
    parser.add_argument("--min-aligned-bp", type=int, default=100)

    args = parser.parse_args()

    if args.min_aligned_bp < 1:
        parser.error("--min-aligned-bp must be positive")

    adapter_options = (
        args.adapter_ref_dir,
        args.adapter_min_overlap,
        args.adapter_min_identity,
    )
    if any(value is not None for value in adapter_options) and not all(
        value is not None for value in adapter_options
    ):
        parser.error("adapter reference audit requires all three adapter options")
    if args.adapter_min_overlap is not None and args.adapter_min_overlap < 1:
        parser.error("--adapter-min-overlap must be positive")
    if args.adapter_min_identity is not None and not 0 <= args.adapter_min_identity <= 1:
        parser.error("--adapter-min-identity must be between 0 and 1")

    adapter_references = []
    adapter_audit_status = "NOT_CONFIGURED"
    if args.adapter_ref_dir:
        adapter_references = load_adapter_references(args.adapter_ref_dir)
        if not adapter_references:
            parser.error("--adapter-ref-dir contains no FASTA sequences")
        adapter_audit_status = "EVALUATED"

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    with open(args.table, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fieldnames = list(reader.fieldnames or [])
        all_rows = list(reader)

    target_rows = [
        row
        for row in all_rows
        if row["initial_assignment"] == "target_candidate"
    ]

    target_ids = {
        row["fragment_id"]
        for row in target_rows
    }

    primary_records = parse_primary_bam(
        args.bam,
        args.samtools,
        target_ids,
    )

    species_counts = Counter()
    support_counts = Counter()
    exact_sequence_groups = defaultdict(list)
    coordinate_groups = defaultdict(list)

    derived_rows = []

    for row in target_rows:
        fragment_id = row["fragment_id"]
        label = species_label(row)

        species_counts[label] += 1
        support_counts[row["pair_support_class"]] += 1

        r1 = primary_records.get(fragment_id, {}).get("R1", {})
        r2 = primary_records.get(fragment_id, {}).get("R2", {})

        r1_sequence = r1.get("sequence", "")
        r2_sequence = r2.get("sequence", "")

        sequence_key = ""

        if r1_sequence or r2_sequence:
            sequence_key = hashlib.sha256(
                f"{r1_sequence}|{r2_sequence}".encode()
            ).hexdigest()

            exact_sequence_groups[sequence_key].append(fragment_id)

        coordinate_key = (
            r1.get("reference", ""),
            r1.get("position", ""),
            r1.get("strand", ""),
            r1.get("CIGAR", ""),
            r2.get("reference", ""),
            r2.get("position", ""),
            r2.get("strand", ""),
            r2.get("CIGAR", ""),
        )

        coordinate_groups[coordinate_key].append(fragment_id)

        r1_aligned = to_int(row["R1_primary_aligned_bases"])
        r2_aligned = to_int(row["R2_primary_aligned_bases"])

        r1_mapq = to_int(row["R1_primary_MAPQ"])
        r2_mapq = to_int(row["R2_primary_MAPQ"])

        r1_nm = to_int(row["R1_primary_NM"])
        r2_nm = to_int(row["R2_primary_NM"])

        r1_softclip = to_int(row["R1_primary_softclip"])
        r2_softclip = to_int(row["R2_primary_softclip"])

        aligned_values = [
            value
            for value in (r1_aligned, r2_aligned)
            if value is not None
        ]

        mapq_values = [
            value
            for value in (r1_mapq, r2_mapq)
            if value is not None
        ]

        softclip_values = [
            value
            for value in (r1_softclip, r2_softclip)
            if value is not None
        ]

        nm_rates = []

        if r1_nm is not None and r1_aligned:
            nm_rates.append(r1_nm / r1_aligned)

        if r2_nm is not None and r2_aligned:
            nm_rates.append(r2_nm / r2_aligned)

        adapter_r1 = int(
            ADAPTER_CORE in r1_sequence
            or ADAPTER_SHORT in r1_sequence
        )

        adapter_r2 = int(
            ADAPTER_CORE in r2_sequence
            or ADAPTER_SHORT in r2_sequence
        )

        ref_r1 = ref_r2 = None
        if adapter_references:
            ref_r1 = best_adapter_overlap(
                r1_sequence, adapter_references,
                args.adapter_min_overlap, args.adapter_min_identity,
            )
            ref_r2 = best_adapter_overlap(
                r2_sequence, adapter_references,
                args.adapter_min_overlap, args.adapter_min_identity,
            )

        competitor_observed = int(
            bool(row["R1_best_competitor_AS"])
            or bool(row["R2_best_competitor_AS"])
        )

        derived = dict(row)

        derived.update({
            "target_species_label": label,
            "effective_min_aligned_bp": args.min_aligned_bp,
            "R1_primary_position": r1.get("position", ""),
            "R1_primary_strand": r1.get("strand", ""),
            "R2_primary_position": r2.get("position", ""),
            "R2_primary_strand": r2.get("strand", ""),
            "minimum_primary_aligned_bases": (
                min(aligned_values) if aligned_values else ""
            ),
            "minimum_primary_MAPQ": (
                min(mapq_values) if mapq_values else ""
            ),
            "maximum_primary_softclip": (
                max(softclip_values) if softclip_values else ""
            ),
            "maximum_primary_NM_rate": (
                max(nm_rates) if nm_rates else ""
            ),
            "competitor_alignment_observed": competitor_observed,
            "adapter_core_in_R1": adapter_r1,
            "adapter_core_in_R2": adapter_r2,
            "adapter_core_in_either": int(adapter_r1 or adapter_r2),
            "adapter_ref_audit_status": adapter_audit_status,
            "adapter_ref_hit_R1": int(bool(ref_r1)) if adapter_references else "",
            "adapter_ref_hit_R2": int(bool(ref_r2)) if adapter_references else "",
            "adapter_ref_hit_either": int(bool(ref_r1 or ref_r2)) if adapter_references else "",
            "adapter_ref_name_R1": ref_r1[2] if ref_r1 else "",
            "adapter_ref_name_R2": ref_r2[2] if ref_r2 else "",
            "adapter_ref_overlap_R1": ref_r1[0] if ref_r1 else "",
            "adapter_ref_overlap_R2": ref_r2[0] if ref_r2 else "",
            "adapter_ref_identity_R1": ref_r1[1] if ref_r1 else "",
            "adapter_ref_identity_R2": ref_r2[1] if ref_r2 else "",
            "adapter_ref_matches_R1": ref_r1[3] if ref_r1 else "",
            "adapter_ref_matches_R2": ref_r2[3] if ref_r2 else "",
            "adapter_ref_offset_R1": ref_r1[4] if ref_r1 else "",
            "adapter_ref_offset_R2": ref_r2[4] if ref_r2 else "",
            "exact_pair_sequence_sha256": sequence_key,
            "coordinate_signature": "|".join(map(str, coordinate_key)),
        })

        derived_rows.append(derived)

    exact_group_sizes = {
        key: len(fragment_ids)
        for key, fragment_ids in exact_sequence_groups.items()
    }

    coordinate_group_sizes = {
        key: len(fragment_ids)
        for key, fragment_ids in coordinate_groups.items()
    }

    for row in derived_rows:
        row["exact_pair_sequence_group_size"] = exact_group_sizes.get(
            row["exact_pair_sequence_sha256"],
            0,
        )

        coordinate_tuple = tuple(
            row["coordinate_signature"].split("|")
        )

        row["coordinate_group_size"] = coordinate_group_sizes.get(
            coordinate_tuple,
            0,
        )

    derived_rows.sort(
        key=lambda row: (
            0 if row["pair_support_class"] == "paired_target_consistent" else 1,
            -int(row["minimum_primary_aligned_bases"] or 0),
            -int(row["minimum_primary_MAPQ"] or 0),
            float(row["maximum_primary_NM_rate"] or 999),
            row["fragment_id"],
        )
    )

    detailed_path = Path(
        f"{args.output_prefix}.target_candidates.audit.tsv"
    )

    with open(
        detailed_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        derived_fieldnames = input_fieldnames + [
            "target_species_label",
            *([] if "effective_min_aligned_bp" in input_fieldnames else ["effective_min_aligned_bp"]),
            "R1_primary_position", "R1_primary_strand",
            "R2_primary_position", "R2_primary_strand",
            "minimum_primary_aligned_bases", "minimum_primary_MAPQ",
            "maximum_primary_softclip", "maximum_primary_NM_rate",
            "competitor_alignment_observed",
            "adapter_core_in_R1", "adapter_core_in_R2",
            "adapter_core_in_either", "adapter_ref_audit_status",
            "adapter_ref_hit_R1", "adapter_ref_hit_R2",
            "adapter_ref_hit_either", "adapter_ref_name_R1",
            "adapter_ref_name_R2", "adapter_ref_overlap_R1",
            "adapter_ref_overlap_R2", "adapter_ref_identity_R1",
            "adapter_ref_identity_R2", "adapter_ref_matches_R1",
            "adapter_ref_matches_R2", "adapter_ref_offset_R1",
            "adapter_ref_offset_R2", "exact_pair_sequence_sha256",
            "coordinate_signature", "exact_pair_sequence_group_size",
            "coordinate_group_size",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=derived_fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(derived_rows)

    species_path = Path(
        f"{args.output_prefix}.target_species_counts.tsv"
    )

    with open(
        species_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_species_label", "fragments"])

        for label, count in species_counts.most_common():
            writer.writerow([label, count])

    exact_duplicates_path = Path(
        f"{args.output_prefix}.exact_sequence_duplicates.tsv"
    )

    with open(
        exact_duplicates_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "group_size",
            "pair_sequence_sha256",
            "fragment_ids",
        ])

        duplicate_groups = sorted(
            (
                (len(fragment_ids), key, fragment_ids)
                for key, fragment_ids in exact_sequence_groups.items()
                if len(fragment_ids) >= 2
            ),
            reverse=True,
        )

        for size, key, fragment_ids in duplicate_groups:
            writer.writerow([
                size,
                key,
                ";".join(fragment_ids[:50]),
            ])

    coordinate_duplicates_path = Path(
        f"{args.output_prefix}.coordinate_duplicates.tsv"
    )

    with open(
        coordinate_duplicates_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "group_size",
            "coordinate_signature",
            "fragment_ids",
        ])

        duplicate_groups = sorted(
            (
                (len(fragment_ids), key, fragment_ids)
                for key, fragment_ids in coordinate_groups.items()
                if len(fragment_ids) >= 2
            ),
            reverse=True,
        )

        for size, key, fragment_ids in duplicate_groups:
            writer.writerow([
                size,
                "|".join(map(str, key)),
                ";".join(fragment_ids[:50]),
            ])

    paired_rows = [
        row for row in derived_rows
        if row["pair_support_class"] == "paired_target_consistent"
    ]

    min_aligned_values = [
        int(row["minimum_primary_aligned_bases"])
        for row in derived_rows
        if row["minimum_primary_aligned_bases"] != ""
    ]

    min_mapq_values = [
        int(row["minimum_primary_MAPQ"])
        for row in derived_rows
        if row["minimum_primary_MAPQ"] != ""
    ]

    max_nm_rate_values = [
        float(row["maximum_primary_NM_rate"])
        for row in derived_rows
        if row["maximum_primary_NM_rate"] != ""
    ]

    span_values = [
        int(row["primary_fragment_span"])
        for row in paired_rows
        if row["primary_fragment_span"] != ""
    ]

    summary_path = Path(
        f"{args.output_prefix}.target_candidate_summary.tsv"
    )

    with open(
        summary_path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])

        writer.writerow(["sample_id", args.sample])
        writer.writerow(["effective_min_aligned_bp", args.min_aligned_bp])
        writer.writerow(["target_candidates", len(derived_rows)])
        writer.writerow([
            "paired_target_consistent",
            support_counts["paired_target_consistent"],
        ])
        writer.writerow([
            "single_end_target_candidate",
            support_counts["single_end_target_candidate"],
        ])
        writer.writerow([
            "distinct_target_species_labels",
            len(species_counts),
        ])
        writer.writerow([
            "competitor_alignment_observed",
            sum(
                int(row["competitor_alignment_observed"])
                for row in derived_rows
            ),
        ])
        writer.writerow([
            "adapter_core_in_either_read",
            sum(
                int(row["adapter_core_in_either"])
                for row in derived_rows
            ),
        ])
        writer.writerow([
            "primary_proper_pair",
            sum(
                int(row["primary_proper_pair"])
                for row in derived_rows
            ),
        ])
        writer.writerow([
            "primary_same_contig",
            sum(
                int(row["primary_same_contig"])
                for row in derived_rows
            ),
        ])
        writer.writerow([
            f"both_primary_aligned_ge{args.min_aligned_bp}",
            both_primary_aligned_at_least(paired_rows, args.min_aligned_bp),
        ])
        writer.writerow([
            "both_primary_aligned_ge120",
            sum(
                1
                for row in paired_rows
                if int(row["R1_primary_aligned_bases"] or 0) >= 120
                and int(row["R2_primary_aligned_bases"] or 0) >= 120
            ),
        ])
        writer.writerow([
            "both_primary_MAPQ_ge20",
            sum(
                1
                for row in paired_rows
                if int(row["R1_primary_MAPQ"] or 0) >= 20
                and int(row["R2_primary_MAPQ"] or 0) >= 20
            ),
        ])
        writer.writerow([
            "both_primary_MAPQ_ge30",
            sum(
                1
                for row in paired_rows
                if int(row["R1_primary_MAPQ"] or 0) >= 30
                and int(row["R2_primary_MAPQ"] or 0) >= 30
            ),
        ])

        writer.writerow([
            "unique_exact_pair_sequences",
            len(exact_sequence_groups),
        ])
        writer.writerow([
            "exact_sequence_duplicate_excess",
            len(derived_rows) - len(exact_sequence_groups),
        ])
        writer.writerow([
            "largest_exact_sequence_group",
            max(exact_group_sizes.values(), default=0),
        ])
        writer.writerow([
            "unique_coordinate_signatures",
            len(coordinate_groups),
        ])
        writer.writerow([
            "coordinate_duplicate_excess",
            len(derived_rows) - len(coordinate_groups),
        ])
        writer.writerow([
            "largest_coordinate_group",
            max(coordinate_group_sizes.values(), default=0),
        ])

        summarize_numeric(
            writer,
            "minimum_primary_aligned_bases",
            min_aligned_values,
        )
        summarize_numeric(
            writer,
            "minimum_primary_MAPQ",
            min_mapq_values,
        )
        summarize_numeric(
            writer,
            "maximum_primary_NM_rate",
            max_nm_rate_values,
        )
        summarize_numeric(
            writer,
            "paired_fragment_span",
            span_values,
        )

    print(f"summary={summary_path}")
    print(f"species={species_path}")
    print(f"details={detailed_path}")
    print(f"exact_duplicates={exact_duplicates_path}")
    print(f"coordinate_duplicates={coordinate_duplicates_path}")


if __name__ == "__main__":
    main()
