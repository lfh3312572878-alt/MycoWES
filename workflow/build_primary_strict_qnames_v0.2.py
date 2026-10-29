#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict, Counter

riskfn, df2bed, outdir = sys.argv[1:]


def is_true(v):
    return str(v).strip().lower() in {
        "1","true","yes","y","t"
    }


def is_false(v):
    return str(v).strip().lower() in {
        "0","false","no","n","f",""
    }


def iv(v):
    try:
        return int(v)
    except:
        return None


# ------------------------------------------------------------
# ≥2 independent design-family WES capture regions
# ------------------------------------------------------------

db = defaultdict(list)

with open(df2bed) as f:
    for line in f:

        if not line.strip() or line.startswith("#"):
            continue

        x = line.rstrip().split("\t")

        if len(x) < 3:
            continue

        try:
            s = int(x[1])
            e = int(x[2])
        except:
            continue

        db[x[0]].append((s,e))


for k in db:
    db[k].sort()


def overlap(ref, s, e):

    if not ref or s is None or e is None:
        return False

    if e <= s:
        return False

    for a,b in db.get(ref, []):

        if b <= s:
            continue

        if a >= e:
            break

        if min(b,e) > max(a,s):
            return True

    return False


primary = []
strict = []

counts = Counter()
species = defaultdict(Counter)
core_reasons = Counter()

audit = []


with open(riskfn) as f:

    r = csv.DictReader(f, delimiter="\t")

    for x in r:

        q = (x.get("fragment_id") or "").strip()

        if not q:
            continue

        sp = (
            x.get("target_species_label")
            or "NA"
        ).strip()

        target = (
            (x.get("initial_assignment") or "").strip()
            == "target_candidate"
        )

        core = is_true(
            x.get("exploratory_core_alignment", "")
        )

        riskfree = is_false(
            x.get("risk_any", "")
        )

        adapterfree = is_false(
            x.get("adapter_core_in_either", "")
        )

        if target:
            counts["target_candidate"] += 1

        if target and core:
            counts["core"] += 1

        if target:
            core_reasons[
                (x.get("core_filter_reason_code") or ("PASS" if core else "UNSPECIFIED"))
            ] += 1

        if target and core and riskfree:
            counts["core_riskfree"] += 1

        primary_pass = (
            target
            and core
            and riskfree
            and adapterfree
        )


        h1 = overlap(
            x.get("R1_primary_reference",""),
            iv(x.get("R1_reference_start_0based")),
            iv(x.get("R1_reference_end_0based_exclusive"))
        )

        h2 = overlap(
            x.get("R2_primary_reference",""),
            iv(x.get("R2_reference_start_0based")),
            iv(x.get("R2_reference_end_0based_exclusive"))
        )

        high_capture = h1 or h2

        strict_pass = (
            primary_pass
            and not high_capture
        )


        if primary_pass:

            primary.append(q)

            counts["PRIMARY"] += 1
            species[sp]["PRIMARY"] += 1

            if high_capture:

                counts["HIGH_CAPTURE_DF2"] += 1
                species[sp]["HIGH_CAPTURE_DF2"] += 1


        if strict_pass:

            strict.append(q)

            counts["STRICT"] += 1
            species[sp]["STRICT"] += 1


        if primary_pass:

            audit.append({
                "fragment_id": q,
                "species": sp,
                "high_capture_DF2": int(high_capture),
                "R1_hit": int(h1),
                "R2_hit": int(h2),
                "STRICT_keep": int(strict_pass),
                "exact_pair_sequence_sha256":
                    x.get("exact_pair_sequence_sha256",""),
                "coordinate_signature":
                    x.get("coordinate_signature",""),
            })


with open(outdir + "/PRIMARY.keep_qnames.txt","w") as o:

    for q in sorted(set(primary)):
        o.write(q + "\n")


with open(outdir + "/STRICT_DF2.keep_qnames.txt","w") as o:

    for q in sorted(set(strict)):
        o.write(q + "\n")


with open(outdir + "/filter_summary.tsv","w") as o:

    o.write("metric\tfragments\n")

    for k in [
        "target_candidate",
        "core",
        "core_riskfree",
        "PRIMARY",
        "HIGH_CAPTURE_DF2",
        "STRICT",
    ]:
        o.write(f"{k}\t{counts[k]}\n")


with open(outdir + "/species_filter_summary.tsv","w") as o:

    o.write(
        "species\tPRIMARY\tHIGH_CAPTURE_DF2\tSTRICT\n"
    )

    for sp in sorted(species):

        o.write(
            f"{sp}\t"
            f"{species[sp]['PRIMARY']}\t"
            f"{species[sp]['HIGH_CAPTURE_DF2']}\t"
            f"{species[sp]['STRICT']}\n"
        )


with open(outdir + "/core_filter_reason_summary.tsv","w") as o:

    o.write("stage\treason_code\tfragments\n")

    for reason in sorted(core_reasons):
        o.write(f"CORE_ALIGNMENT\t{reason}\t{core_reasons[reason]}\n")


with open(outdir + "/PRIMARY_fragment_audit.tsv","w") as o:

    cols = [
        "fragment_id",
        "species",
        "high_capture_DF2",
        "R1_hit",
        "R2_hit",
        "STRICT_keep",
        "exact_pair_sequence_sha256",
        "coordinate_signature",
    ]

    o.write("\t".join(cols) + "\n")

    for x in audit:
        o.write(
            "\t".join(str(x[c]) for c in cols) + "\n"
        )


print("target_candidate =", counts["target_candidate"])
print("core =", counts["core"])
print("core_riskfree =", counts["core_riskfree"])
print("PRIMARY =", counts["PRIMARY"])
print("HIGH_CAPTURE_DF2 =", counts["HIGH_CAPTURE_DF2"])
print("STRICT =", counts["STRICT"])
