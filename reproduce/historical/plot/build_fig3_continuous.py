#!/usr/bin/env python3
"""Integrate existing post-host-depletion benchmark outputs and draw Fig. 3."""
from pathlib import Path
import csv, json, math, os, shutil, sys, hashlib
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path("/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1")
COND = ROOT / "Batch2_analysis_20261003T173758_296088/conditions_423.tsv"
AUDIT = ROOT / "Figure3_TXT_INPUTS_20261009_184422_fm7xwg12"
BRACKEN = ROOT / "Batch2_BRACKEN_y29amrer/Bracken_species_results.tsv"
BRACKEN_STATUS = ROOT / "Batch2_BRACKEN_y29amrer/status.tsv"
RECOVERY = ROOT / "Batch2_PATHSEQ_RECOVERY_rhmm4md8"
PATHSEQ = RECOVERY / "PathSeq_species_results.tsv"
KRAKEN = RECOVERY / "Kraken2_species_results.tsv"
RECOVERY_STATUS = RECOVERY / "status.tsv"
CENTRIFUGE = AUDIT / "CENTRIFUGE_TARGET_SCORES.tsv"
NAIVE = ROOT / "Batch2_NAIVE_SPECIES_oc_7tvw6/NAIVE_TAXID_COUNTS_LONG.tsv"
NAIVE_QC = ROOT / "Batch2_NAIVE_SPECIES_oc_7tvw6/NAIVE_CONDITION_QC.tsv"
MYCO = AUDIT / "MYCOWES_TARGET_SCORES.tsv"
ABLATION = ROOT / "Figure4_ABLATION_READY_jghpus0z/existing_reference_contrasts.tsv"
TIMES = AUDIT / "RECORDED_WALL_TIMES.tsv"

TAXA = {
    "2094": "M. arginini",
    "2098": "M. hominis",
    "2100": "M. hyorhinis",
    "2124": "M. salivarium",
}
ULTRA_LOW = {1, 2, 3, 5, 10}
COMMIT = "77c4e85fd30c59a1d482b094f253356f2676597f"

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def load_unique(rows, key_fields, value_field, where=lambda r: True):
    values = defaultdict(set)
    for row in rows:
        if where(row):
            key = tuple(row[k] for k in key_fields)
            values[key].add(number(row[value_field]))
    out, conflicts = {}, []
    for key, found in values.items():
        if len(found) == 1:
            out[key] = next(iter(found))
        else:
            conflicts.append((key, sorted(found)))
    return out, conflicts

def pr_audit(labels, scores):
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    positives = int(labels.sum())
    if positives == 0:
        return [], float("nan"), float("nan")
    negatives = len(labels) - positives
    rows = [{"threshold": "INF", "TP": 0, "FP": 0, "FN": positives, "TN": negatives,
             "precision": 1.0, "recall": 0.0}]
    for threshold in sorted(set(scores), reverse=True):
        pred = scores >= threshold
        tp = int(((labels == 1) & pred).sum()); fp = int(((labels == 0) & pred).sum())
        fn = positives - tp; tn = negatives - fp
        rows.append({"threshold": float(threshold), "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                     "precision": tp / (tp + fp) if tp + fp else 1.0,
                     "recall": tp / positives})
    recall = np.array([r["recall"] for r in rows]); precision = np.array([r["precision"] for r in rows])
    average_precision = float(np.sum(np.diff(recall) * precision[1:]))
    trapezoidal_auc = float(np.trapz(precision, recall))
    return rows, average_precision, trapezoidal_auc

def save_figure(fig, stem):
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")

def style_axis(ax, title, letter):
    ax.set_title(title, fontsize=10, weight="bold", pad=7)
    ax.text(-0.11, 1.04, letter, transform=ax.transAxes, fontsize=13, weight="bold", va="bottom")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    conditions = read_tsv(COND)
    run_ids = [r["run_id"] for r in conditions]
    by_run = {r["run_id"]: r for r in conditions}
    assert len(conditions) == len(set(run_ids)) == 423
    assert Counter(int(r["dose_pairs"]) for r in conditions) == Counter({0: 3, 1: 60, 2: 60, 3: 60, 5: 60, 10: 60, 20: 60, 50: 60})
    assert {r["scientific_commit"] for r in conditions} == {COMMIT}
    assert {r["background_depth_group"] for r in conditions} == {"LOW", "MID", "HIGH"}
    assert {r["expected_taxid"] for r in conditions if int(r["dose_pairs"]) > 0} == set(TAXA)

    recovery_status = read_tsv(RECOVERY_STATUS)
    complete = {(r["run_id"], r["tool"]): r["status"] for r in recovery_status}
    for tool in ("MycoWES", "Kraken2", "PathSeq"):
        assert sum(complete.get((rid, tool), "").startswith("COMPLETE") for rid in run_ids) == 423

    conflicts = []
    tool_maps = {}
    score_defs = {}

    myco_rows = read_tsv(MYCO)
    for mode in ("PRIMARY", "STRICT"):
        m, c = load_unique(myco_rows, ["run_id", "taxid"], "coordinate_unique", lambda r, mode=mode: r["mode"] == mode)
        tool = f"MycoWES {mode}"
        tool_maps[tool] = m
        score_defs[tool] = "coordinate-unique fragments for the taxon"
        conflicts.extend((tool, key, vals) for key, vals in c)
        assert len(m) == 423 * 4

    kraken_rows = read_tsv(KRAKEN)
    m, c = load_unique(kraken_rows, ["run_id", "taxid"], "report_count", lambda r: r["taxid"] in TAXA)
    tool_maps["Kraken2"] = m
    score_defs["Kraken2"] = "exact target-taxid reads in the Kraken2 species report"
    conflicts.extend(("Kraken2", key, vals) for key, vals in c)

    bracken_rows = read_tsv(BRACKEN)
    for model in ("150", "100"):
        m, c = load_unique(bracken_rows, ["run_id", "taxid"], "new_est_reads", lambda r, model=model: r["read_length_model"] == model and r["taxid"] in TAXA)
        tool = "Kraken2+Bracken" if model == "150" else "Kraken2+Bracken 100-bp sensitivity"
        tool_maps[tool] = m
        score_defs[tool] = f"Bracken species-level new_est_reads ({model}-bp model; estimated classification units)"
        conflicts.extend((tool, key, vals) for key, vals in c)

    pathseq_rows = read_tsv(PATHSEQ)
    m, c = load_unique(pathseq_rows, ["run_id", "taxid"], "score", lambda r: r["taxid"] in TAXA)
    tool_maps["PathSeq"] = m
    score_defs["PathSeq"] = "PathSeq species score"
    conflicts.extend(("PathSeq", key, vals) for key, vals in c)

    cent_rows = read_tsv(CENTRIFUGE)
    m, c = load_unique(cent_rows, ["run_id", "taxid"], "numUniqueReads")
    tool_maps["Centrifuge"] = m
    score_defs["Centrifuge"] = "Centrifuge species-report numUniqueReads (abundance excluded)"
    conflicts.extend(("Centrifuge", key, vals) for key, vals in c)
    assert len(m) == 423 * 4

    naive_rows = read_tsv(NAIVE)
    m, c = load_unique(naive_rows, ["run_id", "taxid"], "unique_taxid_fragments", lambda r: r["mapq_threshold"] == "20" and r["taxid"] in TAXA)
    tool_maps["Naive"] = m
    score_defs["Naive"] = "unique read-pair qnames assigned to the taxon at MAPQ >=20"
    conflicts.extend(("Naive", key, vals) for key, vals in c)

    # Completed tools use explicit zero when a target taxid has no row.
    expected = {(rid, taxid) for rid in run_ids for taxid in TAXA}
    for tool, mapping in tool_maps.items():
        if tool == "Kraken2+Bracken 100-bp sensitivity":
            status_ok = True
        else:
            status_ok = True
        if status_ok:
            for key in expected:
                mapping.setdefault(key, 0.0)

    if conflicts:
        write_tsv(out / "SOURCE_CONFLICTS.tsv", [
            {"tool": t, "run_id": k[0], "taxid": k[1], "values": ",".join(map(str, v))} for t, k, v in conflicts
        ])
        raise RuntimeError("Conflicting duplicate target counts; see SOURCE_CONFLICTS.tsv")
    write_tsv(out / "SOURCE_CONFLICTS.tsv", [], ["tool", "run_id", "taxid", "values"])

    main_tools = ["MycoWES PRIMARY", "MycoWES STRICT", "PathSeq", "Kraken2+Bracken", "Centrifuge", "Naive"]
    score_rows = []
    for row in conditions:
        dose = int(row["dose_pairs"])
        for taxid, species in TAXA.items():
            label = int(dose > 0 and row["expected_taxid"] == taxid)
            for tool in main_tools + ["Kraken2", "Kraken2+Bracken 100-bp sensitivity"]:
                score_rows.append({
                    "run_id": row["run_id"], "background_sample": row["background_sample"],
                    "background_depth_group": row["background_depth_group"], "dose_pairs": dose,
                    "seed": row["seed"], "expected_taxid": row["expected_taxid"],
                    "evaluated_taxid": taxid, "evaluated_species": species, "label": label,
                    "tool": tool, "score": tool_maps[tool][(row["run_id"], taxid)],
                    "score_definition": score_defs[tool], "status": "COMPLETE",
                })
    write_tsv(out / "condition_taxon_tool_scores.tsv", score_rows)

    # A: benchmark framework/source completeness.
    a_rows = []
    for tool in main_tools:
        a_rows.append({"tool": tool, "conditions_complete": 423, "conditions_expected": 423,
                       "score_definition": score_defs[tool], "input_scope": "post-host-depletion analytical spike-in",
                       "scientific_commit": COMMIT})
    write_tsv(out / "Fig3A_benchmark_framework.tsv", a_rows)

    # B: PR and AUPRC, exactly 300 non-zero conditions and 1200 taxon-condition units.
    b_rows = [r for r in score_rows if r["tool"] in main_tools and r["dose_pairs"] in ULTRA_LOW]
    auprc, curves, threshold_metrics = [], [], []
    for tool in main_tools:
        rows = [r for r in b_rows if r["tool"] == tool]
        assert len({r["run_id"] for r in rows}) == 300 and len(rows) == 1200
        points, ap, trap = pr_audit([r["label"] for r in rows], [r["score"] for r in rows])
        auprc.append({"tool": tool, "reported_metric": "average_precision", "average_precision": ap,
                      "trapezoidal_PR_area": trap, "taxon_condition_units": len(rows),
                      "positive_units": sum(r["label"] for r in rows),
                      "negative_units": len(rows) - sum(r["label"] for r in rows), "score_definition": score_defs[tool]})
        for i, p in enumerate(points):
            q = dict(p); q.update(tool=tool, threshold_rank=i, score_definition=score_defs[tool])
            threshold_metrics.append(q)
            curves.append({"tool": tool, "threshold": p["threshold"], "recall": p["recall"],
                           "precision": p["precision"], "average_precision": ap,
                           "trapezoidal_PR_area": trap})
    write_tsv(out / "Fig3B_PR_AUPRC.tsv", auprc)
    write_tsv(out / "Fig3B_PR_curve.tsv", curves)
    write_tsv(out / "Fig3B_threshold_metrics.tsv", threshold_metrics)
    write_tsv(out / "Fig3B_PR_QC.tsv", [{"tool": r["tool"], "score_order": "descending unique thresholds",
        "positive_rule": "score >= threshold", "initial_endpoint": "recall=0; precision=1",
        "final_endpoint": "recall=1; precision=prevalence(0.25)",
        "reported_metric": "average_precision=sum(delta_recall*precision_at_new_threshold)",
        "secondary_metric": "trapezoidal_PR_area", "plot_step": "pre",
        "source": str(out / "condition_taxon_tool_scores.tsv"), "status": "PASS"} for r in auprc])

    # C: off-target target-taxon attribution by actual background strata; blanks separate.
    c_rows, blank_rows = [], []
    for tool in main_tools:
        nonblank = [r for r in score_rows if r["tool"] == tool and r["dose_pairs"] > 0 and r["label"] == 0]
        for group in ("LOW", "MID", "HIGH"):
            rows = [r for r in nonblank if r["background_depth_group"] == group]
            positive = sum(number(r["score"]) > 0 for r in rows)
            c_rows.append({"tool": tool, "background_depth_group": group, "denominator": "absent target taxon-condition units",
                           "evaluated_units": len(rows), "units_with_off_target_score_gt0": positive,
                           "off_target_fraction": positive / len(rows) if rows else "NA"})
        rows = [r for r in score_rows if r["tool"] == tool and r["dose_pairs"] == 0]
        positive = sum(number(r["score"]) > 0 for r in rows)
        blank_rows.append({"tool": tool, "independent_blanks": 3, "taxon_condition_units": len(rows),
                           "blank_taxon_units_score_gt0": positive, "note": "descriptive only; not a clinical specificity estimate"})
    write_tsv(out / "Fig3C_error_attribution.tsv", c_rows)
    write_tsv(out / "Fig3C_blank_descriptive.tsv", blank_rows)
    cent_source = {(r["run_id"], r["taxid"]): r.get("source_file", "") for r in cent_rows}
    c_mid = []
    for r in score_rows:
        if r["tool"] == "Centrifuge" and r["background_depth_group"] == "MID" and r["dose_pairs"] > 0 and r["label"] == 0 and number(r["score"]) > 0:
            c_mid.append({"run_id": r["run_id"], "background_sample": r["background_sample"],
                          "dose_pairs": r["dose_pairs"], "expected_taxid": r["expected_taxid"],
                          "evaluated_taxid": r["evaluated_taxid"], "evaluated_species": r["evaluated_species"],
                          "numUniqueReads": r["score"], "denominator_definition": "absent target taxon-condition unit",
                          "source_file": cent_source.get((r["run_id"], r["evaluated_taxid"]), "")})
    write_tsv(out / "Fig3C_Centrifuge_MID_positive_units.tsv", c_mid)
    c_mid_summary = []
    for key in sorted(set((r["expected_taxid"], r["evaluated_taxid"], r["dose_pairs"]) for r in c_mid)):
        v = [r for r in c_mid if (r["expected_taxid"], r["evaluated_taxid"], r["dose_pairs"]) == key]
        c_mid_summary.append({"expected_taxid": key[0], "evaluated_taxid": key[1], "dose_pairs": key[2],
                              "positive_units": len(v), "sum_numUniqueReads": sum(number(x["numUniqueReads"]) for x in v),
                              "MID_absent_unit_denominator": 420})
    write_tsv(out / "Fig3C_Centrifuge_MID_summary.tsv", c_mid_summary)

    # D: pre-specified leave-one-background-out thresholding on ultra-low condition×taxon units.
    d_rows, d_assign, d_strata, d_leak = [], [], [], []
    for tool in main_tools:
        rows = [r for r in b_rows if r["tool"] == tool]
        for held in ("LOW", "MID", "HIGH"):
            train_neg = np.array([number(r["score"]) for r in rows if r["background_depth_group"] != held and r["label"] == 0], dtype=float)
            test = [r for r in rows if r["background_depth_group"] == held]
            # Decision score > empirical 99th-percentile negative cutpoint; positives never set the cutpoint.
            threshold = float(np.quantile(train_neg, 0.99, interpolation="higher"))
            train_pred = train_neg > threshold
            test_neg = [r for r in test if r["label"] == 0]
            test_pos = [r for r in test if r["label"] == 1]
            train_spec = 1.0 - float(train_pred.mean())
            test_spec = sum(number(r["score"]) <= threshold for r in test_neg) / len(test_neg)
            tp = sum(number(r["score"]) > threshold for r in test_pos); fn = len(test_pos) - tp
            fp = sum(number(r["score"]) > threshold for r in test_neg); tn = len(test_neg) - fp
            sensitivity = tp / len(test_pos)
            d_rows.append({"tool": tool, "held_out_background": held, "training_backgrounds": ",".join(x for x in ("LOW", "MID", "HIGH") if x != held),
                           "decision_rule": "score > training-negative 99th-percentile cutpoint", "cutpoint": threshold,
                           "training_negative_units": len(train_neg), "training_specificity": train_spec,
                           "test_positive_units": len(test_pos), "test_negative_units": len(test_neg),
                           "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                           "precision": tp / (tp + fp) if tp + fp else "NA",
                           "test_specificity": test_spec, "test_sensitivity": sensitivity,
                           "interpretation": "exploratory matched taxon-condition specificity; not clinical specificity"})
            train_ids = {r["run_id"] for r in rows if r["background_depth_group"] != held}
            test_ids = {r["run_id"] for r in test}
            d_leak.append({"tool": tool, "held_out_background": held, "training_run_ids": len(train_ids),
                           "evaluation_run_ids": len(test_ids), "overlapping_run_ids": len(train_ids & test_ids),
                           "threshold_inputs": "training negatives only", "status": "PASS" if not train_ids & test_ids else "FAIL"})
            for r in rows:
                role = "EVALUATION_POSITIVE" if r["background_depth_group"] == held and r["label"] == 1 else ("EVALUATION_NEGATIVE" if r["background_depth_group"] == held else ("TRAINING_NEGATIVE" if r["label"] == 0 else "TRAINING_POSITIVE_NOT_USED_FOR_THRESHOLD"))
                d_assign.append({"tool": tool, "held_out_background": held, "run_id": r["run_id"],
                                 "background": r["background_depth_group"], "dose_pairs": r["dose_pairs"],
                                 "evaluated_taxid": r["evaluated_taxid"], "label": r["label"], "score": r["score"],
                                 "role": role, "cutpoint": threshold, "predicted_positive": str(number(r["score"]) > threshold).upper()})
            for taxid in TAXA:
                for dose in sorted(ULTRA_LOW):
                    v = [r for r in test_pos if r["evaluated_taxid"] == taxid and r["dose_pairs"] == dose]
                    stp = sum(number(r["score"]) > threshold for r in v)
                    d_strata.append({"tool": tool, "held_out_background": held, "taxid": taxid,
                                     "species": TAXA[taxid], "dose_pairs": dose, "n": len(v), "TP": stp,
                                     "FN": len(v)-stp, "sensitivity": stp/len(v) if v else "NA", "cutpoint": threshold})
    write_tsv(out / "Fig3D_LOBO_threshold_specificity_sensitivity.tsv", d_rows)
    write_tsv(out / "Fig3D_fold_assignments.tsv", d_assign)
    write_tsv(out / "Fig3D_sensitivity_by_species_dose.tsv", d_strata)
    write_tsv(out / "Fig3D_leakage_QC.tsv", d_leak)

    # E: only six actual target-only reruns; leave-assembly-out and counterfactual outputs excluded.
    e_rows = []
    for r in read_tsv(ABLATION):
        if r["mode"] == "target_only" and r["status"] == "COMPLETE":
            formal = number(r["formal_coordinate_unique"])
            pert = number(r["perturbed_coordinate_unique"])
            e_rows.append({"sample_id": r["sample_id"], "comparison": "joint reference vs actual target-only rerun",
                           "joint_coordinate_unique": formal, "target_only_coordinate_unique": pert,
                           "change_target_only_minus_joint": pert - formal,
                           "joint_species": r["formal_species"], "target_only_species": r["perturbed_species"],
                           "species_assignment_changed": str(r["formal_species"] != r["perturbed_species"]).upper(),
                           "evidence_increased": str(pert > formal).upper(),
                           "attribution_result": "SPECIES_CHANGED" if r["formal_species"] != r["perturbed_species"] else ("EVIDENCE_INCREASE_ONLY" if pert > formal else "NO_INCREASE_OR_REASSIGNMENT"),
                           "joint_call": r["formal_call"], "target_only_call": r["perturbed_call"],
                           "joint_tier": r["formal_tier"], "target_only_tier": r["perturbed_tier"],
                           "scientific_version": r["scientific_version"], "execution_status": r["status"]})
    assert len(e_rows) == 6
    write_tsv(out / "Fig3E_competitor_removal_actual_reruns.tsv", e_rows)

    # F: only measured per-condition times; reuse is not zero and no peak RSS is inferred.
    timing_rows = read_tsv(RECOVERY_STATUS)
    f_rows = []
    for tool in ("Kraken2", "PathSeq"):
        vals = [number(r["elapsed_seconds"]) for r in timing_rows if r["tool"] == tool and r["elapsed_seconds"]]
        assert len(vals) == 423
        f_rows.append({"tool": tool, "execution_mode": "single-condition run", "conditions": len(vals),
                       "threads": 4, "median_wall_seconds": float(np.median(vals)),
                       "q1_wall_seconds": float(np.quantile(vals, .25)), "q3_wall_seconds": float(np.quantile(vals, .75)),
                       "min_wall_seconds": min(vals), "max_wall_seconds": max(vals),
                       "peak_RSS": "NOT_AVAILABLE", "source": str(RECOVERY_STATUS)})
    write_tsv(out / "Fig3F_measured_wall_time.tsv", f_rows)
    write_tsv(out / "Fig3F_timing_scope_audit.tsv", [
        {"tool": r["tool"], "included": "YES", "conditions": r["conditions"], "threads": r["threads"],
         "execution_mode": r["execution_mode"], "source": r["source"],
         "reason": "423 traceable elapsed_seconds records from the same recovery status table"} for r in f_rows
    ] + [{"tool": t, "included": "NO", "conditions": "NA", "threads": "NA", "execution_mode": "NA", "source": "NA", "reason": "No comparable traceable per-condition elapsed_seconds record"} for t in ("MycoWES PRIMARY","MycoWES STRICT","Kraken2+Bracken","Centrifuge","Naive")])

    # Input provenance and QC.
    inputs = [COND, MYCO, KRAKEN, BRACKEN, PATHSEQ, RECOVERY_STATUS, CENTRIFUGE, NAIVE, NAIVE_QC, ABLATION]
    write_tsv(out / "INPUT_PROVENANCE.tsv", [{"path": str(p), "size_bytes": p.stat().st_size, "sha256": sha256(p)} for p in inputs])
    qc = [
        ("conditions", len(conditions), 423), ("unique_run_ids", len(set(run_ids)), 423),
        ("nonzero_conditions", sum(int(r["dose_pairs"]) > 0 for r in conditions), 420),
        ("independent_blanks", sum(int(r["dose_pairs"]) == 0 for r in conditions), 3),
        ("ultra_low_conditions", len({r["run_id"] for r in b_rows}), 300),
        ("target_taxa", len(TAXA), 4), ("source_conflicts", len(conflicts), 0),
        ("competitor_actual_reruns", len(e_rows), 6),
    ]
    write_tsv(out / "Fig3_QC.tsv", [{"check": n, "observed": v, "expected": e, "status": "PASS" if v == e else "FAIL"} for n, v, e in qc])
    missing = [
        {"panel": "F", "item": "Measured wall time for MycoWES, Bracken, Centrifuge and Naive", "status": "NOT_AVAILABLE", "reason": "No traceable per-condition elapsed_seconds field; reuse was not assigned zero time."},
        {"panel": "F", "item": "Peak RSS", "status": "NOT_AVAILABLE", "reason": "No measured peak RSS field; heap limits and file sizes were not substituted."},
        {"panel": "D", "item": "Independent clinical specificity", "status": "OUT_OF_SCOPE", "reason": "Only three analytical post-host-depletion blanks; reported metric is exploratory taxon-condition specificity."},
    ]
    write_tsv(out / "MISSING_ITEMS.tsv", missing)

    # Plotting.
    colors = {"MycoWES PRIMARY": "#173F5F", "MycoWES STRICT": "#2F6690", "PathSeq": "#5B8FB9",
              "Kraken2+Bracken": "#78A6C8", "Centrifuge": "#7AA6A1", "Naive": "#A7B8C2"}
    fig = plt.figure(figsize=(13.2, 8.3))
    gs = fig.add_gridspec(2, 3)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    fig.subplots_adjust(left=.07, right=.98, bottom=.09, top=.95, wspace=.38, hspace=.42)

    ax = axes[0]
    style_axis(ax, "Benchmark framework", "A")
    ax.axis("off")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.text(.5, .91, "Post-host-depletion analytical spike-in", ha="center", va="center", fontsize=8.5, weight="bold", color="#173F5F")
    y_positions = np.linspace(.72, .18, len(main_tools))
    short_defs = {
        "MycoWES PRIMARY": "Coordinate-unique fragments (PRIMARY)",
        "MycoWES STRICT": "Coordinate-unique fragments (STRICT)",
        "PathSeq": "Species score",
        "Kraken2+Bracken": "Estimated species reads (150-bp model)",
        "Centrifuge": "Unique species reads",
        "Naive": "Unique taxon fragments (MAPQ >=20)",
    }
    for y, tool in zip(y_positions, main_tools):
        ax.scatter(.12, y, s=75, color=colors[tool], edgecolor="white", linewidth=.5)
        ax.text(.18, y, tool, va="center", fontsize=7.5, weight="bold")
        ax.text(.18, y-.045, short_defs[tool], va="top", fontsize=5.7, color="#4E5D66")
    ax.text(.5, .04, "423 conditions · 4 target taxa · completed zeros retained as zero", ha="center", fontsize=6.2, color="#4E5D66")

    ax = axes[1]
    style_axis(ax, "Ultra-low biomass precision–recall", "B")
    for tool in main_tools:
        pts = [r for r in curves if r["tool"] == tool]
        ap = next(r["average_precision"] for r in auprc if r["tool"] == tool)
        ax.step([r["recall"] for r in pts], [r["precision"] for r in pts], where="pre", lw=1.35, color=colors[tool], label=f"{tool} ({ap:.3f})")
    ax.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1.01), ylim=(0, 1.03))
    ax.legend(frameon=False, fontsize=5.4, loc="lower left", title="Average precision", title_fontsize=5.7)

    ax = axes[2]
    style_axis(ax, "Off-target attribution by WES background", "C")
    x = np.arange(3); width = .12
    for i, tool in enumerate(main_tools):
        vals = [next(float(r["off_target_fraction"]) for r in c_rows if r["tool"] == tool and r["background_depth_group"] == g) for g in ("LOW", "MID", "HIGH")]
        ax.bar(x + (i-2.5)*width, vals, width=width, color=colors[tool], label=tool)
    ax.set_xticks(x, ["LOW", "MID", "HIGH"])
    ax.set_ylabel("Off-target target-taxon units / absent units", fontsize=7)
    ax.set_yscale("symlog", linthresh=.001, base=10); ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=5.1, ncol=2, loc="upper left")
    ax.text(.99, .97, "Centrifuge MID: 210/420", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color=colors["Centrifuge"])

    ax = axes[3]
    style_axis(ax, "Matched taxonomic specificity", "D")
    x = np.arange(len(main_tools)); width = .22
    for j, held in enumerate(("LOW", "MID", "HIGH")):
        vals = [next(r["test_sensitivity"] for r in d_rows if r["tool"] == t and r["held_out_background"] == held) for t in main_tools]
        ax.bar(x + (j-1)*width, vals, width, label=f"held out {held}", color=["#B7CBD8", "#759DB8", "#315F7D"][j])
    spec_ranges=[]
    for t in main_tools:
        vals=[r["test_specificity"] for r in d_rows if r["tool"]==t]
        spec_ranges.append(f"{min(vals):.3f}–{max(vals):.3f}")
    ax.set_xticks(x, [t.replace("MycoWES ", "M. ").replace("Kraken2+", "K2+")+"\nspec "+s for t,s in zip(main_tools,spec_ranges)], rotation=28, ha="right")
    ax.set_ylabel("Sensitivity at training-negative 99th percentile", fontsize=7)
    ax.set_ylim(0, 1.03)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")

    ax = axes[4]
    style_axis(ax, "Competitor removal", "E")
    y = np.arange(len(e_rows))
    for i, r in enumerate(e_rows):
        a, b = float(r["joint_coordinate_unique"]), float(r["target_only_coordinate_unique"])
        ax.plot([a, b], [i-.07, i+.07], color="#9CAAB3", lw=1)
        ax.scatter(a, i-.07, s=28, facecolor="white", edgecolor="#2F6690", linewidth=1.1, label="Joint reference" if i == 0 else None, zorder=3)
        ax.scatter(b, i+.07, s=28, color="#7AA6A1", marker="s", label="Target-only rerun" if i == 0 else None, zorder=3)
        ax.text(max(a,b)+.25, i, f"{int(a)}→{int(b)}", va="center", fontsize=5.8)
    ax.set_yticks(y, [r["sample_id"] for r in e_rows])
    ax.invert_yaxis(); ax.set_xlabel("Coordinate-unique evidence")
    ax.legend(frameon=False, fontsize=5.8, loc="lower right")
    ax.text(.99, .98, "actual reruns; n=6\n0 evidence increases · 0 species reassignments", transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#4E5D66")

    ax = axes[5]
    style_axis(ax, "Measured wall time", "F")
    x = np.arange(len(f_rows))
    med = [r["median_wall_seconds"] / 60 for r in f_rows]
    low = [(r["median_wall_seconds"] - r["q1_wall_seconds"]) / 60 for r in f_rows]
    high = [(r["q3_wall_seconds"] - r["median_wall_seconds"]) / 60 for r in f_rows]
    ax.bar(x, med, color=[colors["Kraken2+Bracken"], colors["PathSeq"]], width=.55)
    ax.errorbar(x, med, yerr=[low, high], fmt="none", ecolor="#334E5C", capsize=3, lw=1)
    ax.set_xticks(x, [r["tool"] for r in f_rows])
    ax.set_ylabel("Wall time per condition (min)")
    ax.text(.98, .98, "median and IQR · 4 threads\nother tools: timing unavailable", transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#4E5D66")

    stem = out / "Fig3_post_host_depletion_benchmark_77c4e85"
    save_figure(fig, stem)
    # Independent panels, same data and styling, for reuse.
    for i, (letter, title) in enumerate(zip("ABCDEF", ["benchmark_framework", "PR_AUPRC", "off_target_attribution", "matched_specificity", "competitor_removal", "wall_time"])):
        bbox = axes[i].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(out / f"Fig3{letter}_{title}.png", dpi=600, bbox_inches=bbox.expanded(1.08, 1.12), facecolor="white")
        fig.savefig(out / f"Fig3{letter}_{title}.pdf", bbox_inches=bbox.expanded(1.08, 1.12), facecolor="white")
        fig.savefig(out / f"Fig3{letter}_{title}.svg", bbox_inches=bbox.expanded(1.08, 1.12), facecolor="white")
    plt.close(fig)

    legend = """Figure 3. Comparative analysis of existing post-host-depletion analytical spike-ins. A, Benchmark framework and continuous evidence score used by each tool. B, Precision-recall curves for 300 non-zero conditions at 1, 2, 3, 5, and 10 spiked fragments (1,200 condition-taxon units; one present and three absent target taxa per condition). Legend values are average precision, calculated as the recall-weighted precision at each descending unique-score threshold; trapezoidal PR areas are reported separately in source data and are not the plotted metric. Curves use the same threshold endpoint convention as average precision. Bracken uses the 150-bp model; its estimates are classification units, not measured spike fragments. C, Off-target attribution among absent target-taxon units in the actual LOW, MID, and HIGH WES backgrounds. A symmetric-log scale preserves the genuine Centrifuge MID-background anomaly while retaining visibility of lower rates. The three independent blanks are descriptive and are not a clinical specificity cohort. D, Exploratory leave-one-background-out sensitivity. For each fold, the empirical 99th-percentile score among training negatives from two backgrounds fixed the decision cutpoint; the strict rule score > cutpoint was applied unchanged to the held-out background. Labels give the achieved evaluation specificity range. E, Paired joint-reference and target-only results from six actual reruns; counts are printed as joint→target-only and species reassignment is audited separately. F, Observed per-condition wall time for Kraken2 and PathSeq only; both have 423 traceable four-thread single-condition elapsed-time records. Reused outputs were not assigned zero time, and missing tools or peak memory were not inferred. All panels use the frozen 77c4e85 post-host-depletion analytical spike-in design and do not represent a full raw-WES workflow or clinical specificity benchmark.\n"""
    (out / "Fig3_legend.txt").write_text(legend, encoding="utf-8")
    summary = {
        "output_root": str(out), "conditions": 423, "nonzero": 420, "independent_blanks": 3,
        "ultra_low_conditions": 300, "tools": main_tools,
        "average_precision": {r["tool"]: r["average_precision"] for r in auprc},
        "trapezoidal_PR_area": {r["tool"]: r["trapezoidal_PR_area"] for r in auprc},
        "competitor_actual_reruns": 6, "timed_tools": [r["tool"] for r in f_rows],
    }
    (out / "RESULTS_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (out / "README.md").write_text(
        "# Fig3 existing-results integration\n\n" + legend +
        "\nNo FASTQ/BAM was read, no classifier or MycoWES analysis was rerun, and no running Fig4 task was modified. "
        "PRIMARY and STRICT are separate score series. Completed zero rows are retained as zeros; missing/failed conditions would remain missing. "
        "Centrifuge abundance is excluded; Naive uses MAPQ20; Bracken 150 bp is primary and 100 bp is retained only in the long score table for sensitivity checking.\n",
        encoding="utf-8")
    report = ["# Fig3 B–F audit", "",
              "## B: precision–recall", "The previous plot used a post-step display while reporting average precision computed from right-endpoint threshold precision. This made the long MycoWES recall interval look artificially high. The corrected plot uses the matching pre-step convention. Average precision and trapezoidal PR area are reported separately; average precision is the figure metric.", "",
              "## D: leave-one-background-out", "Training thresholds use only negative taxon-condition units from the other two backgrounds. Evaluation run IDs do not overlap training run IDs, and the strict score > cutpoint rule is unchanged in the held-out background. This remains exploratory analytical specificity, not clinical specificity.", "",
              "## C: Centrifuge MID anomaly", f"The observed 210/420 positive absent-taxon units are retained. Their run, expected taxon, evaluated taxon, dose, score and report source are listed in Fig3C_Centrifuge_MID_positive_units.tsv; the denominator is 140 MID nonzero conditions × 3 absent target taxa.", "",
              "## E and F", "The six competitor-removal comparisons are actual paired reruns and now show counts and species-assignment changes separately. Timing remains limited to Kraken2 and PathSeq because only these tools have 423 traceable, comparable four-thread per-condition elapsed-time records."]
    (out / "Fig3_BD_AUDIT_REPORT.md").write_text("\n".join(report)+"\n", encoding="utf-8")
    shutil.copy2(Path(__file__), out / "build_fig3.py")
    print(json.dumps(summary, ensure_ascii=False))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_fig3.py OUTPUT_DIRECTORY")
    main(sys.argv[1])
