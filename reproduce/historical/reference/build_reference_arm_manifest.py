#!/usr/bin/env python3
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


FORMAL_ROOT = Path("/data/h011/project/ref/MycoWES_Ref_v0.1/25_reference_qc_beta_v0.28")
FORMAL_TAX = FORMAL_ROOT / "06_metadata/AB_CDE_reference_taxonomy_lineage_map.reference_qc_beta_v0.28.tsv"
FORMAL_ID = FORMAL_ROOT / "06_metadata/identity.tsv"
FORMAL_FASTA = FORMAL_ROOT / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa"
FORMAL_MMI = FORMAL_ROOT / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi"
SUP = Path("/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/FIG345_MISSING_MATERIALS_SUPPLEMENT_20261014T182500/SF1")
RANKING = SUP / "reference_chain/frozen_preANI_selection_v0.8/MycoWES_preANI_candidate_ranking_v0.8.tsv"
RETAINED = SUP / "reference_chain/frozen_final_AD_assembly_decisions_v0.21/MycoWES_AD_retained_fasta_manifest_v0.21.tsv"
CLUSTERS = SUP / "reference_chain/frozen_near_clone_clusters_v0.11/MycoWES_near_clone_cluster_membership_v0.11.tsv"
SOURCE_META = Path("/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Fig3_Fig4_TXT_COMPLETION_20261010_v1/SPIKE_SOURCE_METADATA.tsv")
SOURCE_PROV = Path("/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Fig3_Fig4_TXT_COMPLETION_20261010_v1/SPIKE_SOURCE_PROVENANCE_UPDATE.tsv")
STAGE3 = Path("/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/MYCOWES_METHOD_RESEARCH_20261015T023520/phase3")
RUNS = {"ERR14129335", "SRR1655810", "SRR15402804", "SRR35637520"}


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def as_int(value, default=10**9):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def main(out_arg):
    out = Path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "STATUS.json").exists():
        raise SystemExit(f"refusing to overwrite existing result: {out}")
    required = [FORMAL_TAX, FORMAL_ID, FORMAL_FASTA, FORMAL_MMI, RANKING, RETAINED, CLUSTERS, SOURCE_META, SOURCE_PROV]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise SystemExit("missing required sources: " + ";".join(missing))

    formal_rows = read_tsv(FORMAL_TAX)
    rank_rows = read_tsv(RANKING)
    retained_rows = read_tsv(RETAINED)
    cluster_rows = read_tsv(CLUSTERS)
    with FORMAL_ID.open(encoding="utf-8") as handle:
        formal_identity = dict(line.rstrip("\n").split("\t", 1) for line in handle if line.strip())

    assemblies = {}
    for row in formal_rows:
        accession = row["assembly_accession"]
        key = (accession, row["species_taxid"], row["reference_class"])
        if key not in assemblies:
            assemblies[key] = {
                "assembly_accession": accession,
                "species_taxid": row["species_taxid"],
                "species_name": row["ncbi_species"],
                "reference_class": row["reference_class"],
                "role": row["role"],
                "group": row["group"],
                "contig_records": 0,
            }
        assemblies[key]["contig_records"] += 1

    rank = {r["Assembly Accession"]: r for r in rank_rows}
    retained = {r["assembly_accession"]: r for r in retained_rows}
    cluster = {r["assembly_accession"]: r for r in cluster_rows}
    target_by_species = defaultdict(list)
    for row in assemblies.values():
        if row["reference_class"] == "TARGET":
            row["quality_rank"] = rank.get(row["assembly_accession"], {}).get("quality_rank_within_species_v0.8", "")
            row["selection_order"] = rank.get(row["assembly_accession"], {}).get("selection_order_within_species_v0.8", "")
            row["cluster_id"] = cluster.get(row["assembly_accession"], {}).get("cluster_id", "")
            row["formal_retained_v021"] = "YES" if row["assembly_accession"] in retained else "NO"
            target_by_species[row["species_taxid"]].append(row)

    one_rep = set()
    selection_reason = {}
    unranked_species = []
    for taxid, rows in target_by_species.items():
        ranked = sorted(rows, key=lambda r: (as_int(r["quality_rank"]), as_int(r["selection_order"]), r["assembly_accession"]))
        if not ranked[0]["quality_rank"]:
            unranked_species.append(taxid)
        chosen = ranked[0]
        one_rep.add(chosen["assembly_accession"])
        selection_reason[chosen["assembly_accession"]] = (
            "minimum frozen quality_rank_within_species_v0.8; then selection_order; then accession"
        )

    multi_rep = {r["assembly_accession"] for rows in target_by_species.values() for r in rows if r["assembly_accession"] in retained}
    multi_missing = sorted(t for t, rows in target_by_species.items() if not any(r["assembly_accession"] in multi_rep for r in rows))
    full_target = {r["assembly_accession"] for rows in target_by_species.values() for r in rows}
    if unranked_species or multi_missing:
        raise SystemExit(f"objective metadata incomplete: unranked={unranked_species}; multi_missing={multi_missing}")
    assert len(one_rep) == len(target_by_species)
    assert multi_rep.issubset(full_target)

    inventory = []
    manifest = []
    for row in sorted(assemblies.values(), key=lambda r: (r["reference_class"], r["species_taxid"], r["assembly_accession"])):
        accession = row["assembly_accession"]
        is_target = row["reference_class"] == "TARGET"
        q = rank.get(accession, {})
        c = cluster.get(accession, {})
        inventory.append({
            **row,
            "quality_rank_within_species_v0.8": q.get("quality_rank_within_species_v0.8", ""),
            "assembly_level": q.get("Assembly Level", ""),
            "refseq_category": q.get("Assembly Refseq Category", ""),
            "checkm_completeness": q.get("CheckM completeness", ""),
            "checkm_contamination": q.get("CheckM contamination", ""),
            "contig_n50": q.get("Assembly Stats Contig N50", ""),
            "n_contigs": q.get("Assembly Stats Number of Contigs", ""),
            "cluster_id_v0.11": c.get("cluster_id", ""),
            "formal_retained_v0.21": "YES" if accession in retained else ("NOT_APPLICABLE_COMPETITOR" if not is_target else "NO"),
        })
        manifest.append({
            "assembly_accession": accession,
            "species_taxid": row["species_taxid"],
            "species_name": row["species_name"],
            "reference_class": row["reference_class"],
            "role": row["role"],
            "contig_records": row["contig_records"],
            "include_full_v028": "YES",
            "include_one_rep_target": "YES" if (not is_target or accession in one_rep) else "NO",
            "include_multi_nonredundant_target": "YES" if (not is_target or accession in multi_rep) else "NO",
            "include_species_aggregated_full": "YES",
            "quality_rank_v0.8": q.get("quality_rank_within_species_v0.8", ""),
            "cluster_id_v0.11": c.get("cluster_id", ""),
            "selection_reason_one_rep": selection_reason.get(accession, "competitor unchanged" if not is_target else "not selected"),
            "selection_reason_multi": retained.get(accession, {}).get("final_AD_decision_reason_v0.21", "competitor unchanged" if not is_target else "not retained by frozen v0.21 decision"),
        })

    arms = [
        ("FULL_V028", "ACTUAL_FLOW_REFERENCE", "reuse formal v0.28 FASTA/MMI", "READY_REUSE"),
        ("ONE_REP_TARGET", "ACTUAL_FLOW_REFERENCE", "one frozen-quality representative per TARGET species; COMPETITOR unchanged", "MANIFEST_READY_INDEX_MISSING"),
        ("MULTI_NONREDUNDANT_TARGET", "ACTUAL_FLOW_REFERENCE", "frozen v0.21 retained TARGET assemblies; COMPETITOR unchanged", "ALIAS_FULL_NO_DUPLICATE_INDEX" if multi_rep == full_target else "MANIFEST_READY_INDEX_MISSING"),
        ("SPECIES_AGGREGATED_FULL", "DERIVED_SPECIES_AGGREGATION", "reuse full v0.28 alignment; aggregate competition support by species without changing mapping", "READY_REUSE_FULL_ALIGNMENT"),
    ]
    summary = []
    for arm, mode, definition, status in arms:
        field = {
            "FULL_V028": "include_full_v028",
            "ONE_REP_TARGET": "include_one_rep_target",
            "MULTI_NONREDUNDANT_TARGET": "include_multi_nonredundant_target",
            "SPECIES_AGGREGATED_FULL": "include_species_aggregated_full",
        }[arm]
        inc = [r for r in manifest if r[field] == "YES"]
        summary.append({
            "arm_id": arm,
            "execution_class": mode,
            "definition": definition,
            "status": status,
            "assemblies_total": len(inc),
            "assemblies_target": sum(r["reference_class"] == "TARGET" for r in inc),
            "assemblies_competitor": sum(r["reference_class"] == "COMPETITOR" for r in inc),
            "species_total": len({r["species_taxid"] for r in inc}),
            "target_species": len({r["species_taxid"] for r in inc if r["reference_class"] == "TARGET"}),
            "contig_records": sum(int(r["contig_records"]) for r in inc),
            "requires_new_mapping": "YES" if mode == "ACTUAL_FLOW_REFERENCE" and arm != "FULL_V028" and status != "ALIAS_FULL_NO_DUPLICATE_INDEX" else "NO",
        })

    source_meta = {r["run_accession"]: r for r in read_tsv(SOURCE_META) if r["run_accession"] in RUNS}
    source_prov = {r["run_accession"]: r for r in read_tsv(SOURCE_PROV) if r["run_accession"] in RUNS}
    source_identity = []
    for run in sorted(RUNS):
        m, p = source_meta[run], source_prov[run]
        source_identity.append({
            "source_run": run,
            "taxid": m["tax_id"],
            "species": m["scientific_name"],
            "assembly_accession": m["assembly_accessions"],
            "assembly_source_type": m["assembly_source_type"],
            "local_assembly_file": m["assembly_file"],
            "local_assembly_md5_status": m["assembly_md5_status"],
            "assembly_in_formal_target": "YES" if m["assembly_accessions"] in full_target else "NO",
            "identity_status": "RESOLVED_ACCESSION",
            "remaining_gap": "source assembly sequence linkage unavailable locally" if not m["assembly_file"] else "none for accession/file identity; divergence still unevaluated",
            "provenance_interpretation": p["interpretation"],
        })
    assert len(source_identity) == 4

    stage3_files = [p for p in STAGE3.rglob("*") if p.is_file()]
    gate_ready = any("freeze" in p.name.lower() and "status" in p.name.lower() for p in stage3_files)
    gate = [{
        "gate": "STAGE3_SPLIT_FREEZE_BEFORE_CANDIDATE_PERFORMANCE",
        "status": "PASS" if gate_ready else "BLOCKED",
        "evidence": ";".join(str(p) for p in stage3_files) if stage3_files else "phase3 has no files",
        "action": "candidate pilot/performance forbidden until explicit split-freeze status exists",
        "sealed_membership_read": "NO",
    }]

    write_tsv(out / "REFERENCE_ASSEMBLY_INVENTORY.tsv", inventory, list(inventory[0]))
    write_tsv(out / "REFERENCE_ARM_ASSEMBLY_MANIFEST.tsv", manifest, list(manifest[0]))
    write_tsv(out / "ARM_SUMMARY.tsv", summary, list(summary[0]))
    write_tsv(out / "SOURCE_ASSEMBLY_IDENTITY.tsv", source_identity, list(source_identity[0]))
    write_tsv(out / "STAGE3_GATE.tsv", gate, list(gate[0]))
    index_rows = [
        {"arm_id": "FULL_V028", "fasta": str(FORMAL_FASTA), "index": str(FORMAL_MMI), "reuse_status": "VALID_RELEASE_INDEX_REUSE", "identity_sha256": formal_identity.get("joint_mmi_sha256", "")},
        {"arm_id": "ONE_REP_TARGET", "fasta": "NOT_BUILT", "index": "NOT_BUILT", "reuse_status": "NO_MATCHING_EXISTING_INDEX_FOUND", "identity_sha256": ""},
        {"arm_id": "MULTI_NONREDUNDANT_TARGET", "fasta": str(FORMAL_FASTA) if multi_rep == full_target else "NOT_BUILT", "index": str(FORMAL_MMI) if multi_rep == full_target else "NOT_BUILT", "reuse_status": "ALIAS_FULL_IDENTICAL_ASSEMBLY_SET" if multi_rep == full_target else "NO_MATCHING_EXISTING_INDEX_FOUND", "identity_sha256": formal_identity.get("joint_mmi_sha256", "") if multi_rep == full_target else ""},
        {"arm_id": "SPECIES_AGGREGATED_FULL", "fasta": str(FORMAL_FASTA), "index": str(FORMAL_MMI), "reuse_status": "REUSE_FULL_MAPPING_DERIVED_AGGREGATION", "identity_sha256": formal_identity.get("joint_mmi_sha256", "")},
    ]
    write_tsv(out / "INDEX_REUSE_AUDIT.tsv", index_rows, list(index_rows[0]))

    criteria = """# Reference-structure arm freeze (pre-performance)\n\n+## Fixed scope\n+\n+- Development data are the historical 300 low-dose conditions only. No sealed-test membership or score was read.\n+- Competitor C/D/E assemblies remain unchanged in all physical-reference arms; the intervention isolates within-target reference multiplicity.\n+- The formal v0.28 library is immutable and reused as the baseline.\n+\n+## One representative per target species\n+\n+Choose the assembly with the smallest frozen `quality_rank_within_species_v0.8`; break ties by frozen `selection_order_within_species_v0.8`, then accession. These ranks were generated before this Stage 2 experiment from assembly level, RefSeq category, type-material flag, CheckM completeness/contamination, N50, contig count, total length, provenance and prespecified human relevance. Benchmark outcomes are not inputs.\n+\n+## Multiple nonredundant target representatives\n+\n+Use the frozen v0.21 retained decisions derived from v0.10 within-species ANI, v0.11 near-clone clusters and v0.12 representative review. Keep cluster representatives and retained diversity assemblies; do not reselect by development performance. If this retained target set equals formal v0.28, the arm is an explicit no-op alias and no duplicate index is built.\n+\n+## Species-aggregated competition evidence\n+\n+Reuse the full v0.28 mapping. Aggregate support across assemblies sharing `species_taxid` while keeping cross-species ambiguity separate. This is a derived scoring/routing experiment, not a changed reference and not evidence for remapping behavior until an actual-flow rule implementation is validated.\n+\n+## Gates\n+\n+Any physical reference change requires actual remapping/routing, a new FASTA/index hash and a correctness pilot. Candidate performance is blocked until an explicit Stage 3 split-freeze status exists. The sealed membership must remain unread.\n+"""
    (out / "SELECTION_CRITERIA.md").write_text(criteria, encoding="utf-8")

    sources = []
    for path in required:
        sources.append({
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path) if path.stat().st_size < 10_000_000 else "RECORDED_IN_FORMAL_RELEASE_OR_SOURCE_MANIFEST",
            "purpose": {
                FORMAL_TAX: "formal contig/assembly/species/class inventory",
                FORMAL_ID: "formal reference hashes",
                FORMAL_FASTA: "immutable full reference FASTA",
                FORMAL_MMI: "immutable full minimap2 index",
                RANKING: "pre-performance objective assembly ranking",
                RETAINED: "frozen nonredundant retained decisions",
                CLUSTERS: "frozen near-clone cluster membership",
                SOURCE_META: "source-run assembly identity",
                SOURCE_PROV: "source identity provenance and remaining gap",
            }[path],
        })
    write_tsv(out / "SOURCE_INDEX.tsv", sources, list(sources[0]))

    status = {
        "status": "MANIFEST_FROZEN_PERFORMANCE_BLOCKED_STAGE3_SPLIT_FREEZE" if not gate_ready else "MANIFEST_FROZEN_READY_FOR_PILOT",
        "development_conditions": 300,
        "sealed_membership_read": False,
        "formal_target_species": len(target_by_species),
        "formal_target_assemblies": len(full_target),
        "one_rep_target_assemblies": len(one_rep),
        "multi_target_assemblies": len(multi_rep),
        "multi_is_full_alias": multi_rep == full_target,
        "downloads": 0,
        "indexes_built": 0,
        "performance_runs": 0,
    }
    (out / "STATUS.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    assert status["one_rep_target_assemblies"] == status["formal_target_species"]
    checks = []
    for path in sorted(p for p in out.iterdir() if p.is_file() and p.name != "CHECKSUMS.sha256"):
        checks.append(f"{sha256(path)}  {path.name}")
    (out / "CHECKSUMS.sha256").write_text("\n".join(checks) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_reference_arm_manifest.py OUTPUT_DIR")
    main(sys.argv[1])
