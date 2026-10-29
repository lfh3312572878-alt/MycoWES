#!/usr/bin/env python3

import sys

r1fn,r2fn,o1fn,o2fn,groupfn=sys.argv[1:]

def fq(fn):
    with open(fn) as f:
        while True:
            h=f.readline()
            if not h:
                break
            s=f.readline().rstrip()
            f.readline()
            q=f.readline().rstrip()
            name=h[1:].strip().split()[0]
            yield name,s,q

groups={}
order=[]
raw=0

for a,b in zip(fq(r1fn),fq(r2fn)):
    q1,s1,qual1=a
    q2,s2,qual2=b

    if q1 != q2:
        raise SystemExit(f"FASTQ mismatch: {q1} != {q2}")

    raw += 1
    key=(s1.upper(),s2.upper())

    if key not in groups:
        groups[key]={
            "rep":q1,
            "members":[],
            "r1":(s1,qual1),
            "r2":(s2,qual2)
        }
        order.append(key)

    groups[key]["members"].append(q1)

with open(o1fn,"w") as o1, \
     open(o2fn,"w") as o2, \
     open(groupfn,"w") as g:

    g.write("representative_qname\texact_group_size\tmember_qnames\n")

    for key in order:
        x=groups[key]

        rep=x["rep"]
        s1,q1=x["r1"]
        s2,q2=x["r2"]

        o1.write(f"@{rep}\n{s1}\n+\n{q1}\n")
        o2.write(f"@{rep}\n{s2}\n+\n{q2}\n")

        g.write(
            f'{rep}\t{len(x["members"])}\t'
            f'{",".join(x["members"])}\n'
        )

print("raw_paired_fragments =",raw)
print("exact_dedup_fragments =",len(groups))
print("exact_duplicate_excess =",raw-len(groups))
