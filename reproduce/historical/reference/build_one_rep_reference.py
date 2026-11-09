#!/usr/bin/env python3
import csv
import hashlib
import json
import sys
from pathlib import Path


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fasta_records(path):
    header = None
    seq = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line, []
            else:
                if header is None:
                    raise ValueError("sequence before first FASTA header")
                seq.append(line)
    if header is not None:
        yield header, seq


if len(sys.argv) != 5:
    raise SystemExit("usage: build_one_rep_reference.py ARM_MANIFEST FORMAL_TAXONOMY FORMAL_FASTA OUTDIR")
manifest_path, taxonomy_path, fasta_path, outdir = map(Path, sys.argv[1:])
outdir.mkdir(parents=True, exist_ok=True)
out_fasta = outdir / "ONE_REP_TARGET.AB_CDE_JOINT.fa"
out_tax = outdir / "ONE_REP_TARGET.taxonomy.tsv"
out_excl = outdir / "ONE_REP_TARGET.exclusions.tsv"
out_qc = outdir / "BUILD_QC.json"
for path in (out_fasta, out_tax, out_excl, out_qc):
    if path.exists():
        raise SystemExit(f"refusing to overwrite {path}")

manifest = read_tsv(manifest_path)
included_assemblies = {r["assembly_accession"] for r in manifest if r["include_one_rep_target"] == "YES"}
removed_targets = [r for r in manifest if r["reference_class"] == "TARGET" and r["include_one_rep_target"] == "NO"]
taxonomy = read_tsv(taxonomy_path)
selected_taxonomy = [r for r in taxonomy if r["assembly_accession"] in included_assemblies]
selected_ids = {r["reference_id"] for r in selected_taxonomy}
if len(selected_ids) != len(selected_taxonomy):
    raise SystemExit("duplicate reference_id in selected taxonomy")

seen = set()
bp = 0
with out_fasta.open("w", encoding="ascii", newline="") as out:
    for header, seq in fasta_records(fasta_path):
        reference_id = header[1:].split()[0]
        if reference_id not in selected_ids:
            continue
        if reference_id in seen:
            raise SystemExit(f"duplicate FASTA reference id: {reference_id}")
        seen.add(reference_id)
        out.write(header)
        out.writelines(seq)
        bp += sum(len(s.strip()) for s in seq)
if seen != selected_ids:
    missing = sorted(selected_ids - seen)
    raise SystemExit(f"selected taxonomy ids missing from FASTA: {missing[:10]} ({len(missing)})")

with taxonomy_path.open(encoding="utf-8", newline="") as handle:
    fields = next(csv.reader(handle, delimiter="\t"))
with out_tax.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected_taxonomy)

with out_excl.open("w", encoding="utf-8", newline="") as handle:
    fields = ["assembly_accession", "species", "n_contigs_removed", "reason", "evidence_status", "software_version"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    by_assembly = {r["assembly_accession"]: r for r in removed_targets}
    for accession, row in sorted(by_assembly.items()):
        writer.writerow({
            "assembly_accession": accession,
            "species": row["species_name"],
            "n_contigs_removed": row["contig_records"],
            "reason": "DEVELOPMENT_ONE_REP_TARGET_ARM",
            "evidence_status": "PRESPECIFIED_FROZEN_QUALITY_RANK",
            "software_version": "MycoWES_reference_structure_research_20261015",
        })

counts = {}
for cls in ("TARGET", "COMPETITOR"):
    rows = [r for r in selected_taxonomy if r["reference_class"] == cls]
    counts[cls] = {
        "contigs": len(rows),
        "assemblies": len({r["assembly_accession"] for r in rows}),
        "species": len({r["species_taxid"] for r in rows}),
    }
qc = {
    "status": "PASS",
    "records": len(seen),
    "bases": bp,
    "fasta_sha256": sha256(out_fasta),
    "taxonomy_sha256": sha256(out_tax),
    "excluded_target_assemblies": len({r["assembly_accession"] for r in removed_targets}),
    "counts": counts,
}
assert qc["records"] == 7941
assert counts["TARGET"] == {"contigs": 118, "assemblies": 18, "species": 18}
assert counts["COMPETITOR"] == {"contigs": 7823, "assemblies": 560, "species": 553}
out_qc.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
