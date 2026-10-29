#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
from collections import defaultdict, Counter

BIN=250

parser=argparse.ArgumentParser()

parser.add_argument("--bam",required=True)
parser.add_argument("--ab",required=True)
parser.add_argument("--r1",required=True)
parser.add_argument("--r2",required=True)
parser.add_argument("--groups",required=True)
parser.add_argument("--reps",required=True)
parser.add_argument("--out",required=True)
parser.add_argument("--minimap2",required=True)
parser.add_argument("--samtools",required=True)
parser.add_argument("--retention",type=float,default=0.90)

args=parser.parse_args()

OUT=Path(args.out)
OUT.mkdir(parents=True,exist_ok=True)

TMP=OUT/"tmp_assembly_selection"
TMP.mkdir(exist_ok=True)


# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

def tags_from_header(h):
    d={}
    for x in h.split("|"):
        if "=" in x:
            k,v=x.split("=",1)
            d[k]=v
    return d


def cigar_ref_span(cigar):
    n=0
    for length,op in re.findall(r"(\d+)([MIDNSHP=X])",cigar):
        if op in "MDN=X":
            n += int(length)
    return n


def get_AS(fields):
    for x in fields[11:]:
        if x.startswith("AS:i:"):
            try:
                return int(x[5:])
            except:
                return -10**9
    return -10**9


def fq_dict(fn):
    d={}
    with open(fn) as f:
        while True:
            h=f.readline()
            if not h:
                break
            s=f.readline().rstrip()
            f.readline()
            q=f.readline().rstrip()
            name=h[1:].strip().split()[0]
            d[name]=(s,q)
    return d


def fasta_iter(fn):
    h=None
    seq=[]
    with open(fn) as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield h,"".join(seq)
                h=line[1:].strip().split()[0]
                seq=[]
            else:
                seq.append(line.strip())
        if h is not None:
            yield h,"".join(seq)


# ------------------------------------------------------------
# exact-group sizes
# ------------------------------------------------------------

exact_size={}

with open(args.groups) as f:
    r=csv.DictReader(f,delimiter="\t")
    for x in r:
        exact_size[x["representative_qname"]]=int(x["exact_group_size"])


reps=set()

with open(args.reps) as f:
    for line in f:
        q=line.strip()
        if q:
            reps.add(q)


R1=fq_dict(args.r1)
R2=fq_dict(args.r2)


# ------------------------------------------------------------
# 1. assign each exact-dedup fragment to TARGET SPECIES
#    using original AB+CDE BAM
# ------------------------------------------------------------

# q -> mate -> species -> best AS
best_species_mate=defaultdict(
    lambda: defaultdict(
        lambda: defaultdict(lambda:-10**9)
    )
)

cmd=[
    args.samtools,
    "view",
    "-N",args.reps,
    args.bam
]

p=subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    text=True
)

for line in p.stdout:

    x=line.rstrip().split("\t")

    if len(x)<11:
        continue

    q=x[0]

    if q not in reps:
        continue

    flag=int(x[1])

    if flag & 0x4:
        continue

    # exclude supplementary, keep secondary
    if flag & 0x800:
        continue

    mate=1 if flag & 0x40 else 2 if flag & 0x80 else 0

    if mate==0:
        continue

    rname=x[2]
    meta=tags_from_header(rname)

    role=meta.get("role","")

    if not (role.startswith("A") or role.startswith("B")):
        continue

    species=meta.get("species")

    if not species:
        continue

    AS=get_AS(x)

    if AS > best_species_mate[q][mate][species]:
        best_species_mate[q][mate][species]=AS

p.wait()

if p.returncode != 0:
    raise SystemExit("samtools view failed during species assignment")


species_of={}
ambiguous=[]

species_rows=[]

for q in sorted(reps):

    species=set()

    for mate in (1,2):
        species.update(best_species_mate[q][mate].keys())

    ranked=[]

    for sp in species:

        scores=[]

        for mate in (1,2):
            v=best_species_mate[q][mate].get(sp,-10**9)
            if v > -10**8:
                scores.append(v)

        mapped_mates=len(scores)
        total_AS=sum(scores)

        ranked.append((mapped_mates,total_AS,sp))

    if not ranked:
        species_rows.append((q,"NO_TARGET_ALIGNMENT","NA",0,0))
        continue

    ranked.sort(reverse=True)

    best_mates,best_AS,best_sp=ranked[0]

    tied=[
        x for x in ranked
        if x[0]==best_mates and x[1]==best_AS
    ]

    if len(tied)>1:
        ambiguous.append(q)
        species_rows.append(
            (q,"AMBIGUOUS_SPECIES",";".join(sorted(x[2] for x in tied)),
             best_mates,best_AS)
        )
        continue

    species_of[q]=best_sp

    species_rows.append(
        (q,"ASSIGNED",best_sp,best_mates,best_AS)
    )


with open(OUT/"01_species_assignment.tsv","w") as o:

    o.write(
        "qname\tstatus\tspecies\t"
        "target_mapped_mates\ttarget_AS_sum\n"
    )

    for x in species_rows:
        o.write("\t".join(map(str,x))+"\n")


# ------------------------------------------------------------
# 2. parse AB FASTA into species -> assemblies
# ------------------------------------------------------------

assemblies=defaultdict(
    lambda: defaultdict(list)
)

assembly_bp=Counter()
assembly_contigs=Counter()

for h,seq in fasta_iter(args.ab):

    meta=tags_from_header(h)

    species=meta.get("species")
    asm=meta.get("asm")

    if not species or not asm:
        continue

    assemblies[species][asm].append((h,seq))

    assembly_bp[(species,asm)] += len(seq)
    assembly_contigs[(species,asm)] += 1


by_species=defaultdict(list)

for q,sp in species_of.items():
    if q in R1 and q in R2:
        by_species[sp].append(q)


# ------------------------------------------------------------
# 3. assembly testing
# ------------------------------------------------------------

def run_one_assembly(species,asm,records,qnames):

    safe=re.sub(r"[^A-Za-z0-9_.-]+","_",species)

    d=TMP/safe
    d.mkdir(exist_ok=True)

    fa=d/f"{asm}.fa"

    if not fa.exists():

        with open(fa,"w") as o:

            for h,seq in records:

                o.write(f">{h}\n")

                for i in range(0,len(seq),80):
                    o.write(seq[i:i+80]+"\n")


    r1=d/"reads.R1.fastq"
    r2=d/"reads.R2.fastq"

    if not r1.exists():

        with open(r1,"w") as a, open(r2,"w") as b:

            for q in qnames:

                s1,qual1=R1[q]
                s2,qual2=R2[q]

                a.write(f"@{q}\n{s1}\n+\n{qual1}\n")
                b.write(f"@{q}\n{s2}\n+\n{qual2}\n")


    cmd=[
        args.minimap2,
        "-ax","sr",
        "--secondary=yes",
        "-N","100",
        "--MD",
        "-t","2",
        str(fa),
        str(r1),
        str(r2)
    ]

    proc=subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True
    )


    # q -> mate -> contig -> best record
    aln=defaultdict(
        lambda: defaultdict(dict)
    )

    for line in proc.stdout:

        if line.startswith("@"):
            continue

        x=line.rstrip().split("\t")

        if len(x)<11:
            continue

        q=x[0]
        flag=int(x[1])

        if flag & 0x4:
            continue

        if flag & 0x800:
            continue

        mate=1 if flag & 0x40 else 2 if flag & 0x80 else 0

        if mate==0:
            continue

        rname=x[2]

        start=int(x[3])-1
        end=start+cigar_ref_span(x[5])

        AS=get_AS(x)

        rec={
            "qname":q,
            "rname":rname,
            "start":start,
            "end":end,
            "strand":"-" if flag & 0x10 else "+",
            "AS":AS,
            "mapq":int(x[4])
        }

        old=aln[q][mate].get(rname)

        if old is None or AS > old["AS"]:
            aln[q][mate][rname]=rec


    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {species} {asm}")


    retained={}

    total_AS=0

    for q in qnames:

        a=aln[q].get(1,{})
        b=aln[q].get(2,{})

        common=set(a).intersection(b)

        if not common:
            continue

        candidates=[]

        for contig in common:

            r1x=a[contig]
            r2x=b[contig]

            score=r1x["AS"]+r2x["AS"]

            candidates.append(
                (score,contig,r1x,r2x)
            )

        candidates.sort(
            key=lambda x:(-x[0],x[1])
        )

        score,contig,r1x,r2x=candidates[0]

        total_AS += score

        retained[q]={
            "qname":q,
            "species":species,
            "assembly":asm,
            "rname":contig,
            "start":min(r1x["start"],r2x["start"]),
            "end":max(r1x["end"],r2x["end"]),
            "orientation":f'R1{r1x["strand"]}_R2{r2x["strand"]}',
            "mapq1":r1x["mapq"],
            "mapq2":r2x["mapq"],
            "pair_AS":score
        }


    return retained,total_AS


selection_rows=[]
selected={}


for species in sorted(by_species):

    qnames=sorted(by_species[species])

    candidates=[]

    for asm,records in assemblies.get(species,{}).items():

        retained,total_AS=run_one_assembly(
            species,
            asm,
            records,
            qnames
        )

        n=len(retained)
        retention=n/len(qnames) if qnames else 0

        row={
            "species":species,
            "assembly":asm,
            "species_exact_fragments":len(qnames),
            "retained_fragments":n,
            "retention_fraction":retention,
            "total_pair_AS":total_AS,
            "assembly_contigs":assembly_contigs[(species,asm)],
            "assembly_bp":assembly_bp[(species,asm)]
        }

        selection_rows.append(row)

        candidates.append(
            (row,retained)
        )


    if not candidates:
        continue


    candidates.sort(
        key=lambda z:(
            -z[0]["retained_fragments"],
            -z[0]["total_pair_AS"],
            z[0]["assembly_contigs"],
            -z[0]["assembly_bp"],
            z[0]["assembly"]
        )
    )

    bestrow,bestmap=candidates[0]

    selected[species]={
        "row":bestrow,
        "mapping":bestmap
    }


with open(OUT/"02_assembly_selection_all.tsv","w") as o:

    cols=[
        "species",
        "assembly",
        "species_exact_fragments",
        "retained_fragments",
        "retention_fraction",
        "total_pair_AS",
        "assembly_contigs",
        "assembly_bp"
    ]

    o.write("\t".join(cols)+"\n")

    for x in sorted(
        selection_rows,
        key=lambda z:(
            z["species"],
            -z["retained_fragments"],
            -z["total_pair_AS"]
        )
    ):

        vals=[]

        for c in cols:

            v=x[c]

            if isinstance(v,float):
                vals.append(f"{v:.6f}")
            else:
                vals.append(str(v))

        o.write("\t".join(vals)+"\n")


# ------------------------------------------------------------
# 4. coordinate dedup + dispersion metrics
# ------------------------------------------------------------

coverage_rows=[]
bin_rows=[]
coord_rows=[]
fragment_rows=[]
summary_rows=[]


def merged_intervals(ivs):

    if not ivs:
        return []

    ivs=sorted(ivs)

    out=[]

    s,e=ivs[0]

    for ns,ne in ivs[1:]:

        if ns <= e:
            e=max(e,ne)
        else:
            out.append((s,e))
            s,e=ns,ne

    out.append((s,e))

    return out


for species in sorted(selected):

    info=selected[species]

    srow=info["row"]
    mappings=info["mapping"]

    retention=srow["retention_fraction"]

    anchor_status=(
        "PASS"
        if retention >= args.retention
        else "FAIL"
    )


    raw_equiv=sum(
        exact_size.get(q,1)
        for q in by_species[species]
    )

    exact_n=len(by_species[species])
    retained_n=len(mappings)


    if anchor_status=="FAIL":

        summary_rows.append({
            "species":species,
            "selected_assembly":srow["assembly"],
            "raw_equivalent_fragments":raw_equiv,
            "exact_dedup_fragments":exact_n,
            "anchor_retained_fragments":retained_n,
            "anchor_retention_fraction":retention,
            "anchor_status":"FAIL",
            "coordinate_dedup_fragments":"NA",
            "coordinate_duplicate_excess":"NA",
            "covered_bp":"NA",
            "occupied_250bp_bins":"NA",
            "largest_bin_fragments":"NA",
            "largest_bin_fraction":"NA",
            "distinct_contigs":"NA",
            "dispersion_tier":"ANCHOR_FAIL"
        })

        continue


    groups=defaultdict(list)

    for q,x in mappings.items():

        key=(
            x["rname"],
            x["start"],
            x["end"],
            x["orientation"]
        )

        groups[key].append(x)


    unique=[]

    for key,members in groups.items():

        members=sorted(
            members,
            key=lambda x:x["qname"]
        )

        rep=members[0].copy()

        rep["coordinate_group_size"]=len(members)

        rep["raw_equiv_group_size"]=sum(
            exact_size.get(x["qname"],1)
            for x in members
        )

        unique.append(rep)

        coord_rows.append({
            "species":species,
            "assembly":srow["assembly"],
            "rname":key[0],
            "start":key[1],
            "end":key[2],
            "orientation":key[3],
            "representative_qname":rep["qname"],
            "coordinate_group_size":len(members),
            "member_qnames":",".join(x["qname"] for x in members)
        })


    for x in unique:

        fragment_rows.append({
            **x
        })


    coord_n=len(unique)

    by_contig=defaultdict(list)
    bins=Counter()

    for x in unique:

        by_contig[x["rname"]].append(
            (x["start"],x["end"])
        )

        midpoint=(x["start"]+x["end"])//2

        bstart=(midpoint//BIN)*BIN

        bins[(x["rname"],bstart)] += 1


    covered_bp=0

    for rname,ivs in by_contig.items():

        merged=merged_intervals(ivs)

        covered_bp += sum(e-s for s,e in merged)

        for s,e in merged:

            coverage_rows.append(
                (rname,s,e,species,srow["assembly"])
            )


    occupied_bins=len(bins)

    largest=max(bins.values()) if bins else 0

    largest_fraction=(
        largest/coord_n
        if coord_n
        else 0
    )


    for (rname,bstart),n in sorted(bins.items()):

        bin_rows.append({
            "species":species,
            "assembly":srow["assembly"],
            "rname":rname,
            "bin_start":bstart,
            "bin_end":bstart+BIN,
            "fragment_count":n,
            "fraction":n/coord_n if coord_n else 0
        })


    # development dispersion tiers
    if (
        coord_n >= 5
        and
        occupied_bins >= 3
        and
        covered_bp >= 300
        and
        largest_fraction <= 0.50
    ):

        tier="strong"

    elif (
        coord_n >= 3
        and
        occupied_bins >= 2
        and
        largest_fraction < 0.80
    ):

        tier="supportive"

    else:

        tier="weak"


    summary_rows.append({
        "species":species,
        "selected_assembly":srow["assembly"],
        "raw_equivalent_fragments":raw_equiv,
        "exact_dedup_fragments":exact_n,
        "anchor_retained_fragments":retained_n,
        "anchor_retention_fraction":retention,
        "anchor_status":"PASS",
        "coordinate_dedup_fragments":coord_n,
        "coordinate_duplicate_excess":retained_n-coord_n,
        "covered_bp":covered_bp,
        "occupied_250bp_bins":occupied_bins,
        "largest_bin_fragments":largest,
        "largest_bin_fraction":largest_fraction,
        "distinct_contigs":len(by_contig),
        "dispersion_tier":tier
    })


# selected assembly fragment table

with open(OUT/"03_selected_anchor_fragments.tsv","w") as o:

    cols=[
        "species",
        "assembly",
        "qname",
        "rname",
        "start",
        "end",
        "orientation",
        "mapq1",
        "mapq2",
        "pair_AS",
        "coordinate_group_size",
        "raw_equiv_group_size"
    ]

    o.write("\t".join(cols)+"\n")

    for x in sorted(
        fragment_rows,
        key=lambda z:(
            z["species"],
            z["rname"],
            z["start"]
        )
    ):

        o.write(
            "\t".join(str(x[c]) for c in cols)+"\n"
        )


# coordinate groups

with open(OUT/"04_coordinate_duplicate_groups.tsv","w") as o:

    cols=[
        "species",
        "assembly",
        "rname",
        "start",
        "end",
        "orientation",
        "representative_qname",
        "coordinate_group_size",
        "member_qnames"
    ]

    o.write("\t".join(cols)+"\n")

    for x in sorted(
        coord_rows,
        key=lambda z:(
            z["species"],
            z["rname"],
            z["start"]
        )
    ):

        o.write(
            "\t".join(str(x[c]) for c in cols)+"\n"
        )


# merged coverage

with open(OUT/"05_coverage_regions.bed","w") as o:

    for rname,s,e,sp,asm in sorted(
        coverage_rows,
        key=lambda x:(x[0],x[1],x[2])
    ):

        o.write(
            f"{rname}\t{s}\t{e}\t{sp}|{asm}\n"
        )


# bins

with open(OUT/"06_bin_distribution_250bp.tsv","w") as o:

    cols=[
        "species",
        "assembly",
        "rname",
        "bin_start",
        "bin_end",
        "fragment_count",
        "fraction"
    ]

    o.write("\t".join(cols)+"\n")

    for x in bin_rows:

        vals=[]

        for c in cols:
            v=x[c]
            if isinstance(v,float):
                vals.append(f"{v:.6f}")
            else:
                vals.append(str(v))

        o.write("\t".join(vals)+"\n")


# final summary

cols=[
    "species",
    "selected_assembly",
    "raw_equivalent_fragments",
    "exact_dedup_fragments",
    "anchor_retained_fragments",
    "anchor_retention_fraction",
    "anchor_status",
    "coordinate_dedup_fragments",
    "coordinate_duplicate_excess",
    "covered_bp",
    "occupied_250bp_bins",
    "largest_bin_fragments",
    "largest_bin_fraction",
    "distinct_contigs",
    "dispersion_tier"
]

with open(OUT/"07_species_dedup_dispersion_summary.tsv","w") as o:

    o.write("\t".join(cols)+"\n")

    for x in summary_rows:

        vals=[]

        for c in cols:
            v=x[c]

            if isinstance(v,float):
                vals.append(f"{v:.6f}")
            else:
                vals.append(str(v))

        o.write("\t".join(vals)+"\n")


# concise console output

print()
print("===== SPECIES ASSIGNMENT =====")
print("exact representatives =",len(reps))
print("species assigned =",len(species_of))
print("ambiguous species =",len(ambiguous))

print()
print("===== FINAL DEDUP + DISPERSION =====")

for x in summary_rows:

    print(
        x["species"],
        "asm=",x["selected_assembly"],
        "raw=",x["raw_equivalent_fragments"],
        "exact=",x["exact_dedup_fragments"],
        "retained=",x["anchor_retained_fragments"],
        "retention=",x["anchor_retention_fraction"],
        "anchor=",x["anchor_status"],
        "coord=",x["coordinate_dedup_fragments"],
        "covered_bp=",x["covered_bp"],
        "bins=",x["occupied_250bp_bins"],
        "largest_bin_frac=",x["largest_bin_fraction"],
        "tier=",x["dispersion_tier"]
    )
