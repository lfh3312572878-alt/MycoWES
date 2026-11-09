#!/usr/bin/env python3
"""Validate supplied historical summaries without running an analysis pipeline."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import defaultdict


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tier_summary(data):
    by_dose = {}
    for row in data:
        dose = int(row["dose"])
        count = int(row["conditions"])
        denominator = int(row["denominator"])
        tier = row["tier"]
        if tier not in {"None", "Weak", "Supportive", "Strong"} or count < 0 or denominator <= 0:
            raise ValueError("invalid tier summary row")
        group = by_dose.setdefault(dose, {"denominator": denominator, "tiers": {}})
        if denominator != group["denominator"] or tier in group["tiers"]:
            raise ValueError("conflicting denominator or repeated tier")
        if abs(float(row["fraction"]) - count / denominator) > 1e-9:
            raise ValueError("fraction/count mismatch")
        group["tiers"][tier] = count
    for group in by_dose.values():
        if sum(group["tiers"].values()) != group["denominator"]:
            raise ValueError("tier counts do not exhaust denominator")
        group["supportive_or_strong"] = sum(group["tiers"].get(t, 0) for t in ["Supportive", "Strong"])
        group["any_non_none"] = group["denominator"] - group["tiers"].get("None", 0)
    return by_dose


def interval_union(intervals):
    total, previous_contig, end = 0, None, -1
    for contig, start, stop in sorted(intervals):
        if contig != previous_contig or start > end:
            total += stop - start
            end = stop
        elif stop > end:
            total += stop - end
            end = stop
        previous_contig = contig
    return total


def bed_summary(lines):
    groups = defaultdict(list)
    for number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3:
            raise ValueError(f"BED line {number}: fewer than three fields")
        start, stop = int(fields[1]), int(fields[2])
        if not fields[0] or start < 0 or stop <= start:
            raise ValueError(f"BED line {number}: invalid half-open interval")
        category = fields[3].split("|", 1)[0] if len(fields) > 3 else "interval"
        groups[category].append((fields[0], start, stop))
    all_intervals = [interval for group in groups.values() for interval in group]
    return {"records": len(all_intervals), "union_bp": interval_union(all_intervals),
            "categories": {key: {"records": len(group), "union_bp": interval_union(group)}
                           for key, group in sorted(groups.items())}}


def self_test():
    source = [dict(dose="20", tier="Weak", conditions="2", denominator="3", fraction=str(2/3)),
              dict(dose="20", tier="Supportive", conditions="1", denominator="3", fraction=str(1/3))]
    result = tier_summary(source)[20]
    assert result["supportive_or_strong"] == 1 and result["any_non_none"] == 3
    try:
        tier_summary(source + [source[0]])
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate tier was accepted")
    bed = bed_summary(["chr1\t0\t10\tlow|a\n", "chr1\t5\t15\tlow|b\n", "chr2\t0\t5\thigh\n"])
    assert bed["union_bp"] == 20 and bed["categories"]["low"]["union_bp"] == 15
    try:
        bed_summary(["chr1\t4\t3\n"])
    except ValueError:
        pass
    else:
        raise AssertionError("invalid interval was accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parent / "historical")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--risk-bed", type=Path, help="optional existing shared risk BED to summarize")
    parser.add_argument("--capture-bed", type=Path, help="optional existing DF2 capture BED to summarize")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    tiers = tier_summary(rows(args.source_root / "saved_summary/Fig2E_source.tsv"))
    if set(tiers) != {0, 1, 2, 3, 5, 10, 20, 50} or tiers[0]["denominator"] != 3:
        raise ValueError("not the recorded seven-dose plus three-blank design")
    if any(tiers[d]["denominator"] != 60 for d in tiers if d):
        raise ValueError("nonzero dose must contain 60 technical conditions")
    manifest = args.source_root / "metadata/REFERENCE_ARM_ASSEMBLY_MANIFEST.tsv"
    expected_sha = "88a5f0ea9a7d7b225505107bc44968e29648cd23e8d96f6b3babccf7ab6819a6"
    if hashlib.sha256(manifest.read_bytes()).hexdigest() != expected_sha:
        raise ValueError("selection manifest differs from the recorded profile identity")
    targets = [r for r in rows(manifest) if r["reference_class"] == "TARGET"]
    chosen = [r for r in targets if r["include_one_rep_target"] == "YES"]
    if len(targets) != 67 or len(chosen) != 18 or len({r["species_taxid"] for r in chosen}) != 18:
        raise ValueError("unexpected target/representative membership")
    result = {"scope": "historical saved summaries, not current final plot reproduction",
                      "nonzero_conditions": 420, "independent_blanks": 3,
                      "target_assemblies": len(targets), "representatives": len(chosen),
                      "dose_tiers": tiers, "native_mapper_run": False}
    for name, path in [("shared_risk", args.risk_bed), ("DF2_capture", args.capture_bed)]:
        if path is not None:
            with path.open(encoding="utf-8-sig") as handle:
                result[name] = bed_summary(handle)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
