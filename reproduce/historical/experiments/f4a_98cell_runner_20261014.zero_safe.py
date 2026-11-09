#!/usr/bin/env python3
"""Complete the frozen 14-input x 7-architecture Fig4A representative matrix."""
import collections
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path("/data/h011/project")
PAPER = PROJECT / "analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1"
EXEC = PAPER / "Fig3_Fig4_TXT_COMPLETION_EXEC_20261010_v2"
RUNNER = EXEC / "scripts/confirmed_branches.py"
OUT = Path(sys.argv[1]) if len(sys.argv) == 2 else None
WORK_LIMIT = 35 * 1024**3
ARCH = ["FULL", "REMOVE_ONLY_COMPETITOR", "REMOVE_ONLY_CAPTURE_RISK", "REMOVE_ONLY_DEDUP",
        "REMOVE_ONLY_DISPERSION", "REMOVE_ONLY_ENHANCED_HOST", "CUMULATIVE_ALL_REMOVED"]


def now():
    return datetime.now(timezone.utc).isoformat()


def load_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def atomic_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run(cmd, log, record):
    log = Path(log); log.parent.mkdir(parents=True, exist_ok=True)
    started = now()
    with log.open("a") as handle:
        handle.write("COMMAND=" + json.dumps(list(map(str, cmd))) + "\n"); handle.flush()
        rc = subprocess.run(list(map(str, cmd)), stdout=handle, stderr=subprocess.STDOUT).returncode
    atomic_json(record, {"command": list(map(str, cmd)), "started_utc": started, "ended_utc": now(), "exit_code": rc})
    return rc


def tree_bytes(path):
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def ensure_budget(cb, working):
    used = tree_bytes(working) if Path(working).exists() else 0
    if used > WORK_LIMIT:
        raise RuntimeError(f"working cap exceeded: {used} > {WORK_LIMIT}")
    if cb.quota_headroom_gib() < 10:
        raise RuntimeError("STOPPED_LOW_QUOTA")


def qset(path):
    path = Path(path)
    return {x.strip() for x in path.read_text().splitlines() if x.strip()} if path.is_file() else set()


def tier(n, bins, covered, frac, no_dispersion=False):
    if n >= 5 and covered >= 300 and (no_dispersion or (bins >= 3 and frac <= .50)):
        return "strong"
    if n >= 3 and (no_dispersion or (bins >= 2 and frac < .80)):
        return "supportive"
    return "weak" if n else "none"


def metrics(folder, weighted=False, no_dispersion=False):
    folder = Path(folder)
    rows = load_tsv(folder / "07_species_dedup_dispersion_summary.tsv") if (folder / "07_species_dedup_dispersion_summary.tsv").is_file() else []
    groups = load_tsv(folder / "04_coordinate_duplicate_groups.tsv") if (folder / "04_coordinate_duplicate_groups.tsv").is_file() else []
    out = []
    for row in rows:
        n = int(row.get("coordinate_dedup_fragments", 0) or 0)
        bins = int(row.get("occupied_250bp_bins", 0) or 0)
        covered = int(row.get("covered_bp", 0) or 0)
        frac = float(row.get("largest_bin_fraction", 0) or 0)
        if weighted:
            selected = [g for g in groups if g.get("species") == row.get("species") and g.get("assembly") == row.get("selected_assembly")]
            if len(selected) != n:
                raise RuntimeError("coordinate groups do not match formal summary")
            counts = collections.Counter()
            for g in selected:
                midpoint = (int(g["start"]) + int(g["end"])) // 2
                counts[(g["rname"], midpoint // 250)] += int(g.get("coordinate_group_size", 1) or 1)
            n = sum(counts.values()); bins = len(counts); frac = max(counts.values()) / n if n else 0
        out.append({"species": row.get("species", ""), "selected_assembly": row.get("selected_assembly", ""),
                    "evidence_count": n, "occupied_bins": bins, "covered_bp": covered,
                    "largest_bin_fraction": frac, "tier": tier(n, bins, covered, frac, no_dispersion)})
    return {"STRICT_qnames": None, "coordinate_unique": sum(x["evidence_count"] for x in out),
            "tier": ";".join(sorted({x["tier"] for x in out})) or "none", "species_results": out}


def historical_canonical(row):
    species = json.loads(row.get("species_results_json") or "[]")
    strict = sum(int(x.get("raw_equivalent_fragments", 0) or 0) for x in species)
    coords = sum(int(x.get("coordinate_dedup_fragments", 0) or 0) for x in species)
    tiers = ";".join(sorted({x.get("dispersion_tier", "none") for x in species})) or "none"
    return {"STRICT_qnames": strict, "coordinate_unique": coords, "tier": tiers,
            "species": sorted((x.get("species", ""), x.get("selected_assembly", ""),
                               int(x.get("raw_equivalent_fragments", 0) or 0),
                               int(x.get("coordinate_dedup_fragments", 0) or 0),
                               int(x.get("covered_bp", 0) or 0), x.get("dispersion_tier", "none")) for x in species)}


def observed_canonical(runroot):
    b = Path(runroot) / "07_breadth/STRICT_DF2"
    rows = load_tsv(b / "07_species_dedup_dispersion_summary.tsv") if (b / "07_species_dedup_dispersion_summary.tsv").is_file() else []
    strict = len(qset(Path(runroot) / "06_filter/STRICT_DF2.keep_qnames.txt"))
    coords = sum(int(x.get("coordinate_dedup_fragments", 0) or 0) for x in rows)
    tiers = ";".join(sorted({x.get("dispersion_tier", "none") for x in rows})) or "none"
    return {"STRICT_qnames": strict, "coordinate_unique": coords, "tier": tiers,
            "species": sorted((x.get("species", ""), x.get("selected_assembly", ""),
                               int(x.get("raw_equivalent_fragments", 0) or 0),
                               int(x.get("coordinate_dedup_fragments", 0) or 0),
                               int(x.get("covered_bp", 0) or 0), x.get("dispersion_tier", "none")) for x in rows)}


def preserve_small(runroot, dest):
    runroot = Path(runroot); dest = Path(dest)
    wanted = ["run_manifest.tsv", "02_host_depletion/*both_unmapped*.fastq.gz", "03_joint_alignment/*.bam",
              "03_joint_alignment/*.bai", "04_fragment_evidence/*.tsv", "05_risk/*.tsv",
              "06_filter/*", "07_breadth/STRICT_DF2/*"]
    for pattern in wanted:
        for src in runroot.glob(pattern):
            if src.is_file() and src.stat().st_size <= 512 * 1024**2:
                rel = src.relative_to(runroot); target = dest / rel
                target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(src, target)


def run_pipeline(cb, sample, r1, r2, cfg, root, log):
    if not cb.is_complete_run(root):
        rc = run(["bash", cb.ENGINE, "--sample", sample, "--r1", r1, "--r2", r2,
                  "--run-root", root, "--config", cfg], log, Path(log).with_suffix(".command.json"))
        if rc or not cb.is_complete_run(root):
            raise RuntimeError(f"pipeline failed {sample}; rc={rc}")
    return observed_canonical(root)


def rescore_capture(cb, sample, runroot, cfg, out, log):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    keep = Path(runroot) / "06_filter/PRIMARY.keep_qnames.txt"
    if not qset(keep):
        return {"STRICT_qnames": 0, "coordinate_unique": 0, "tier": "none", "species_results": []}
    bam = next((Path(runroot) / "03_joint_alignment").glob("*.qname.bam"))
    raw = next((Path(runroot) / "04_fragment_evidence").glob("*.fragment_evidence.raw.tsv"))
    rc = run(["bash", cb.RUNTIME / "workflow/run_dedup_dispersion_portable.sh", sample, bam, raw, keep, out, cfg],
             log, Path(log).with_suffix(".command.json"))
    if rc:
        raise RuntimeError(f"capture-risk rescore failed {sample}")
    m = metrics(out); m["STRICT_qnames"] = len(qset(keep)); return m


def configure(cb, out):
    target = EXEC / "analyses/TRUTH_CHALLENGE/target_only_reference/target_only.env"
    if not target.is_file():
        raise RuntimeError("frozen target-only config is missing")
    overrides = []
    for line in Path(target).read_text().splitlines():
        if line.startswith(("MYCOWES_JOINT_FASTA=", "MYCOWES_JOINT_MMI=", "MYCOWES_MICROBIAL_IDENTITY=",
                            "MYCOWES_MICROBIAL_PROFILE_ID=", "MYCOWES_REFERENCE_PROFILE=")):
            overrides.append(line)
    mmi = PROJECT / "ref/Fig4_external_resources_20261008_v1/T2T_CHM13v2.0/GRCh38_plus_T2T_CHM13v2.0_maskedY_rCRS.sr.mmi"
    if not mmi.is_file() or mmi.stat().st_size < 10_000_000_000:
        raise RuntimeError("T2T host index failed QC")
    configs = out / "configs"; configs.mkdir(exist_ok=True)
    full = configs / "full_t2t.env"
    full.write_text(Path(cb.CFG).read_text() + f"\nMYCOWES_HOST_MMI={mmi}\nHOST_REFERENCE_PROFILE=GRCh38_plus_T2T_CHM13v2_maskedY_rCRS\nTHREADS=4\nSAMTOOLS_THREADS=2\n")
    competitor = configs / "target_only_t2t.env"
    competitor.write_text(full.read_text() + "\n" + "\n".join(overrides) + "\n")
    cumulative = configs / "target_only_grch38.env"
    cumulative.write_text(Path(cb.CFG).read_text() + "\n" + "\n".join(overrides) + "\nTHREADS=4\nSAMTOOLS_THREADS=2\n")
    audit = [{"path": str(p), "sha256": sha256(p)} for p in (cb.ENGINE, cb.CFG,
             cb.RUNTIME / "workflow/build_primary_strict_qnames_v0.2.py",
             cb.RUNTIME / "workflow/run_dedup_dispersion_portable.sh", full, competitor, cumulative)]
    write_tsv(out / "RUNTIME_IDENTITY.tsv", audit)
    return full, competitor, cumulative


def main():
    if OUT is None:
        raise SystemExit(f"usage: {sys.argv[0]} F4A_OUTPUT_DIR")
    while True:
        state = json.loads((OUT / "STATUS.json").read_text())
        if state.get("status") == "RESTORE_COMPLETE":
            break
        if state.get("status") == "FAILED" or (OUT / "FAILED.json").exists():
            raise RuntimeError("subset restoration failed")
        time.sleep(30)
    lock = (PAPER / ".f4a_98cell_exec.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX)
    spec = importlib.util.spec_from_file_location("confirmed_branches", RUNNER)
    cb = importlib.util.module_from_spec(spec); spec.loader.exec_module(cb)
    working = OUT / "working"; working.mkdir(exist_ok=True)
    full_cfg, competitor_cfg, cumulative_cfg = configure(cb, OUT)
    missing = load_tsv(OUT / "MISSING_70_CELLS.tsv")
    reused = load_tsv(OUT / "REUSED_28_ACTUAL_ENDPOINTS.tsv")
    historical = {r["pair_id"]: r for r in reused if r["architecture"] == "FULL"}
    cells = []
    for r in reused:
        cells.append({"matrix_id": r["matrix_id"], "pair_id": r["pair_id"], "architecture": r["architecture"],
                      "status": "COMPLETE_REUSED", "execution_class": "HISTORICAL_REUSE", "result_path": str(OUT / "REUSED_28_ACTUAL_ENDPOINTS.tsv"), "failure_reason": ""})
    for r in missing:
        cells.append({"matrix_id": r["matrix_id"], "pair_id": r["pair_id"], "architecture": r["architecture"],
                      "status": "PENDING", "execution_class": "", "result_path": "", "failure_reason": ""})
    byid = {x["matrix_id"]: x for x in cells}
    fields = list(cells[0])
    write_tsv(OUT / "CELL_STATUS.tsv", cells, fields)
    atomic_json(OUT / "STATUS.json", {"phase": "ENDPOINT_EQUIVALENCE_AND_70_CELLS", "status": "RUNNING",
                "completed_cells": 28, "failed_cells": 0, "started_utc": now()})
    manifest = load_tsv(OUT / "F4D_28_RUN_MANIFEST_LOCAL.tsv")
    t2t_rows = {r["pair_id"]: r for r in manifest if r["reference"].startswith("GRCh38_T2T")}
    results = []
    for pair_id in sorted(t2t_rows):
        row = t2t_rows[pair_id]; ensure_budget(cb, working)
        pairwork = working / pair_id; pairwork.mkdir(parents=True, exist_ok=True)
        m1 = pairwork / f"{pair_id}.R1.fastq.gz"; m2 = pairwork / f"{pair_id}.R2.fastq.gz"
        base = OUT / "analyses/F4D_28_RUNS/inputs" / f"{row['biological_background']}.base"
        base1 = Path(str(base) + ".R1.fastq.gz"); base2 = Path(str(base) + ".R2.fastq.gz")
        dose = int(row["dose_pairs"])
        truth_meta = []
        if dose:
            truth_root = EXEC / "analyses/TRUTH_CHALLENGE/source_pairs"
            prefix = f"{row['spike_assembly']}.S{row['seed']}"
            tr1 = truth_root / f"{prefix}.PE150.R1.fastq.gz"
            tr2 = truth_root / f"{prefix}.PE150.R2.fastq.gz"
            meta = truth_root / f"{prefix}.truth.tsv"
            if not all(p.is_file() for p in (tr1, tr2, meta)):
                raise RuntimeError(f"frozen truth pair source is missing: {prefix}")
            truth_meta = load_tsv(meta)
            c1 = pairwork / "truth.R1.fastq.gz"; c2 = pairwork / "truth.R2.fastq.gz"
            with gzip.open(c1, "wt") as a, gzip.open(c2, "wt") as b:
                for i, (x, y) in enumerate(zip(cb.gz_records(tr1), cb.gz_records(tr2))):
                    if i >= dose: break
                    a.write("\n".join(x) + "\n"); b.write("\n".join(y) + "\n")
            cb.concat_gz(m1, base1, c1); cb.concat_gz(m2, base2, c2)
        else:
            cb.concat_gz(m1, base1); cb.concat_gz(m2, base2)
        identity = {"pair_id": pair_id, "background": row["biological_background"], "condition": row["condition"],
                    "dose_pairs": dose, "seed": row["seed"], "assembly": row["spike_assembly"],
                    "base_R1_sha256": sha256(base1), "base_R2_sha256": sha256(base2),
                    "mixed_R1_sha256": sha256(m1), "mixed_R2_sha256": sha256(m2),
                    "truth_qnames": [x["qname"] for x in truth_meta[:dose]]}
        atomic_json(OUT / "input_identity" / f"{pair_id}.json", identity)
        ensure_budget(cb, working)
        try:
            fullroot = pairwork / "FULL_RESTORE"
            observed = run_pipeline(cb, pair_id + "__FULL_RESTORE", m1, m2, full_cfg, fullroot, OUT / "logs" / f"{pair_id}.FULL_RESTORE.log")
            expected = historical_canonical(historical[pair_id])
            equivalence = {"pair_id": pair_id, "expected": expected, "observed": observed,
                           "strict_match": expected["STRICT_qnames"] == observed["STRICT_qnames"],
                           "coordinate_match": expected["coordinate_unique"] == observed["coordinate_unique"],
                           "tier_match": expected["tier"] == observed["tier"],
                           "species_terminal_match": expected["species"] == observed["species"]}
            equivalence["status"] = "PASS" if all(equivalence[k] for k in ("strict_match", "coordinate_match", "tier_match", "species_terminal_match")) else "FAIL"
            atomic_json(OUT / "equivalence" / f"{pair_id}.json", equivalence)
            preserve_small(fullroot, OUT / "evidence" / pair_id / "FULL_RESTORE")
            if equivalence["status"] != "PASS":
                raise RuntimeError("historical endpoint equivalence failed")
            cap = rescore_capture(cb, pair_id + "__NO_CAPTURE", fullroot, full_cfg,
                                  pairwork / "capture_rescore", OUT / "logs" / f"{pair_id}.NO_CAPTURE.log")
            cap_id = pair_id + "__REMOVE_ONLY_CAPTURE_RISK"
            atomic_json(OUT / "results" / f"{cap_id}.json", {"matrix_id": cap_id, "pair_id": pair_id,
                        "architecture": "REMOVE_ONLY_CAPTURE_RISK", "execution_class": "SAVED_EVIDENCE_RESCORE", **cap})
            shutil.copytree(pairwork / "capture_rescore", OUT / "evidence" / pair_id / "REMOVE_ONLY_CAPTURE_RISK", dirs_exist_ok=True)
            byid[cap_id].update(status="COMPLETE", execution_class="SAVED_EVIDENCE_RESCORE",
                                result_path=str(OUT / "results" / f"{cap_id}.json"), failure_reason="")
            formal_folder = fullroot / "07_breadth/STRICT_DF2"
            dedup = metrics(formal_folder, weighted=True); dedup["STRICT_qnames"] = observed["STRICT_qnames"]
            dedup_id = pair_id + "__REMOVE_ONLY_DEDUP"
            atomic_json(OUT / "results" / f"{dedup_id}.json", {"matrix_id": dedup_id, "pair_id": pair_id,
                        "architecture": "REMOVE_ONLY_DEDUP", "execution_class": "SAVED_EVIDENCE_COUNTERFACTUAL", **dedup})
            byid[dedup_id].update(status="COMPLETE", execution_class="SAVED_EVIDENCE_COUNTERFACTUAL",
                                  result_path=str(OUT / "results" / f"{dedup_id}.json"), failure_reason="")
            disp = metrics(formal_folder, no_dispersion=True); disp["STRICT_qnames"] = observed["STRICT_qnames"]
            disp_id = pair_id + "__REMOVE_ONLY_DISPERSION"
            atomic_json(OUT / "results" / f"{disp_id}.json", {"matrix_id": disp_id, "pair_id": pair_id,
                        "architecture": "REMOVE_ONLY_DISPERSION", "execution_class": "SAVED_EVIDENCE_COUNTERFACTUAL", **disp})
            byid[disp_id].update(status="COMPLETE", execution_class="SAVED_EVIDENCE_COUNTERFACTUAL",
                                 result_path=str(OUT / "results" / f"{disp_id}.json"), failure_reason="")
            comproot = pairwork / "REMOVE_ONLY_COMPETITOR"
            comp = run_pipeline(cb, pair_id + "__REMOVE_ONLY_COMPETITOR", m1, m2, competitor_cfg, comproot,
                                OUT / "logs" / f"{pair_id}.REMOVE_ONLY_COMPETITOR.log")
            comp_id = pair_id + "__REMOVE_ONLY_COMPETITOR"
            atomic_json(OUT / "results" / f"{comp_id}.json", {"matrix_id": comp_id, "pair_id": pair_id,
                        "architecture": "REMOVE_ONLY_COMPETITOR", "execution_class": "ACTUAL_TARGET_ONLY_T2T_RUN", **comp})
            preserve_small(comproot, OUT / "evidence" / pair_id / "REMOVE_ONLY_COMPETITOR")
            byid[comp_id].update(status="COMPLETE", execution_class="ACTUAL_TARGET_ONLY_T2T_RUN",
                                 result_path=str(OUT / "results" / f"{comp_id}.json"), failure_reason="")
            cumroot = pairwork / "CUMULATIVE_ACTUAL"
            run_pipeline(cb, pair_id + "__CUMULATIVE_ACTUAL", m1, m2, cumulative_cfg, cumroot,
                         OUT / "logs" / f"{pair_id}.CUMULATIVE_ACTUAL.log")
            cum_breadth = pairwork / "cumulative_rescore"
            rescore_capture(cb, pair_id + "__CUMULATIVE_RESCORE", cumroot, cumulative_cfg, cum_breadth,
                            OUT / "logs" / f"{pair_id}.CUMULATIVE_RESCORE.log")
            cumulative = metrics(cum_breadth, weighted=True, no_dispersion=True)
            cumulative["STRICT_qnames"] = len(qset(cumroot / "06_filter/PRIMARY.keep_qnames.txt"))
            cum_id = pair_id + "__CUMULATIVE_ALL_REMOVED"
            atomic_json(OUT / "results" / f"{cum_id}.json", {"matrix_id": cum_id, "pair_id": pair_id,
                        "architecture": "CUMULATIVE_ALL_REMOVED",
                        "execution_class": "ACTUAL_TARGET_ONLY_GRCH38_PLUS_SAVED_EVIDENCE_RESCORE", **cumulative})
            preserve_small(cumroot, OUT / "evidence" / pair_id / "CUMULATIVE_ACTUAL")
            shutil.copytree(cum_breadth, OUT / "evidence" / pair_id / "CUMULATIVE_ALL_REMOVED", dirs_exist_ok=True)
            byid[cum_id].update(status="COMPLETE", execution_class="ACTUAL_TARGET_ONLY_GRCH38_PLUS_SAVED_EVIDENCE_RESCORE",
                                result_path=str(OUT / "results" / f"{cum_id}.json"), failure_reason="")
        except Exception as exc:
            atomic_json(OUT / "failures" / f"{pair_id}.json", {"pair_id": pair_id, "error": repr(exc), "traceback": traceback.format_exc(), "failed_utc": now()})
            for architecture in ARCH:
                mid = pair_id + "__" + architecture
                if mid in byid and byid[mid]["status"] == "PENDING":
                    byid[mid].update(status="BLOCKED", failure_reason=repr(exc))
        finally:
            write_tsv(OUT / "CELL_STATUS.tsv", [byid[k] for k in sorted(byid)], fields)
            completed = sum(x["status"].startswith("COMPLETE") for x in byid.values())
            failed = sum(x["status"] == "BLOCKED" for x in byid.values())
            atomic_json(OUT / "STATUS.json", {"phase": "ENDPOINT_EQUIVALENCE_AND_70_CELLS", "status": "RUNNING",
                        "current_pair": pair_id, "completed_cells": completed, "failed_cells": failed,
                        "updated_utc": now()})
            shutil.rmtree(pairwork, ignore_errors=True)
    result_rows = []
    for path in sorted((OUT / "results").glob("*.json")):
        row = json.loads(path.read_text()); row["species_results_json"] = json.dumps(row.pop("species_results", []), separators=(",", ":")); result_rows.append(row)
    write_tsv(OUT / "MISSING_70_RESULTS.tsv", result_rows)
    complete = sum(x["status"].startswith("COMPLETE") for x in byid.values())
    blocked = sum(x["status"] == "BLOCKED" for x in byid.values())
    final = {"phase": "COMPLETE" if complete == 98 else "PARTIAL", "status": "COMPLETE" if complete == 98 else "PARTIAL",
             "completed_cells": complete, "failed_cells": blocked, "expected_cells": 98, "completed_utc": now()}
    atomic_json(OUT / "STATUS.json", final)
    if complete == 98:
        atomic_json(OUT / "SUCCESS.json", final)
    else:
        atomic_json(OUT / "PARTIAL.json", final)


if __name__ == "__main__":
    main()
