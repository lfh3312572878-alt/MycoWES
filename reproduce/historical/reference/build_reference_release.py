#!/usr/bin/env python3
"""Create the independent v0.28 reference-QC beta profile without touching frozen inputs."""
import csv, hashlib, shutil, sys
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths

root = Path(str(_paths.ref_root / "MycoWES_Ref_v0.1"))
out = _paths.out_dir
analysis = Path(str(_paths.project_root / "analysis/MycoWES_BAM29_hostdepleted/ablation_GCF_900118515/reference_ablation"))
excluded = "GCF_900118515.1"
if out.exists():
    raise SystemExit(f"REF_RELEASE_EXISTS={out}")
for d in ["01_microbe", "05_annotations", "06_metadata", "06_capture_bias", "07_scripts"]:
    (out / d).mkdir(parents=True)

def copy(src, dst):
    shutil.copy2(src, dst)

def filter_text(src, dst):
    with open(src, errors="strict") as i, open(dst, "w") as o:
        for line in i:
            if excluded not in line:
                o.write(line)

def filter_fasta(src, dst):
    kept = removed = 0
    write = False
    with open(src) as i, open(dst, "w") as o:
        for line in i:
            if line.startswith(">"):
                write = excluded not in line
                kept += write; removed += not write
            if write: o.write(line)
    return kept, removed

copy(analysis / "AB_CDE_Joint_Competition_v0.25_without_GCF_900118515.fa", out / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa")
copy(analysis / "AB_CDE_Joint_Competition_v0.25_without_GCF_900118515.sr.mmi", out / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi")
target_kept, target_removed = filter_fasta(root / "16_MycoWES_reference_release_v0.27/01_microbe/AB_TARGET.fa", out / "01_microbe/AB_TARGET.reference_qc_beta_v0.28.fa")

old_meta = root / "16_MycoWES_reference_release_v0.27/06_metadata"
for name in ["AB_CDE_reference_taxonomy_map.v0.27.tsv", "AB_CDE_reference_taxonomy_lineage_map.v0.27.tsv"]:
    filter_text(old_meta / name, out / "06_metadata" / name.replace("v0.27", "reference_qc_beta_v0.28"))

old_ann = root / "16_MycoWES_reference_release_v0.27/05_annotations"
for src in old_ann.iterdir():
    if src.is_file():
        filter_text(src, out / "05_annotations" / src.name.replace("v0.27", "reference_qc_beta_v0.28"))

capture = root / "21_WES_capture_risk_release_v0.1_frozen_20260807/01_panWES/panWES.design_family_support_ge2.bed"
filter_text(capture, out / "06_capture_bias/panWES.design_family_support_ge2.reference_qc_beta_v0.28.bed")
shutil.copytree(root / "14_reference_core_v0.25/04_adapters", out / "05_annotations/adapters")

with open(out / "06_metadata/reference_exclusion_manifest.tsv", "w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["assembly_accession", "species", "n_contigs_removed", "reason", "evidence_status", "date", "software_version"])
    w.writerow([excluded, "Mycoplasmopsis_arginini", 22, "SYSTEMATIC_REFERENCE_ATTRACTOR", "VALIDATED_BY_ABLATION", "2026-08-27", "MycoWES_v0.3-beta1"])

readme = """# MycoWES reference QC beta v0.28

Profile: REFERENCE_QC_BETA

This independent profile excludes assembly GCF_900118515.1 from TARGET and JOINT
assets after validated cohort-wide reference-attractor behavior. The operational
label is REFERENCE_ATTRACTOR_EXCLUDED with reason SYSTEMATIC_REFERENCE_ATTRACTOR.
This does not establish that the source assembly is biologically contaminated.
All v0.27 and earlier frozen reference assets remain unchanged.
"""
(out / "README.md").write_text(readme)
(out / "METHOD.txt").write_text(
    "SOURCE_REFERENCE=16_MycoWES_reference_release_v0.27\n"
    "SOURCE_ABLATION=ablation_GCF_900118515/reference_ablation\n"
    "FILTER=remove every FASTA/taxonomy/annotation record containing assembly accession from exclusion manifest\n"
    "JOINT_MMI=validated prebuilt ablation index; no in-place rebuild or overwrite\n"
    f"TARGET_RECORDS_KEPT={target_kept}\nTARGET_CONTIGS_REMOVED={target_removed}\n"
)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(block)
    return h.hexdigest()

joint = out / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa"
mmi = out / "01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi"
target = out / "01_microbe/AB_TARGET.reference_qc_beta_v0.28.fa"
identity = [
    ("profile_id", "MycoWES_AB_CDE_v0.28_reference_qc_beta"),
    ("reference_release", "MycoWES_Ref_v0.28_REFERENCE_QC_BETA"),
    ("reference_profile", "REFERENCE_QC_BETA"),
    ("excluded_reference_assemblies", excluded),
    ("joint_fasta_sha256", sha(joint)), ("joint_mmi_sha256", sha(mmi)),
    ("target_fasta_sha256", sha(target)),
]
with open(out / "06_metadata/identity.tsv", "w") as h:
    for k, v in identity: h.write(f"{k}\t{v}\n")

with open(out / "SHA256SUMS", "w") as h:
    for path in sorted(x for x in out.rglob("*") if x.is_file() and x.name != "SHA256SUMS"):
        h.write(f"{sha(path)}  {path.relative_to(out).as_posix()}\n")
print(f"REFERENCE_RELEASE={out}")
print(f"TARGET_RECORDS_KEPT={target_kept}")
print(f"TARGET_CONTIGS_REMOVED={target_removed}")
