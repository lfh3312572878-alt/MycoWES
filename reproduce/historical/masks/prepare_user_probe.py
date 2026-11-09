#!/usr/bin/env python3
"""Build a USER_PROBE AB candidate BED using the frozen capture rule."""
import argparse, csv, json, os, re, subprocess
from pathlib import Path

MIN_ALIGNED_BP = 40
MIN_IDENTITY = 0.90
DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.I)
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

def parse_fasta(path):
    records, name, chunks = [], None, []
    with open(path, encoding="utf-8-sig") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0] if line[1:].strip() else "", []
                if not name: raise ValueError(f"line {number}: empty FASTA identifier")
            else:
                if name is None: raise ValueError(f"line {number}: sequence before FASTA header")
                if not DNA.fullmatch(line): raise ValueError(f"line {number}: invalid nucleotide symbol")
                chunks.append(line)
    if name is not None: records.append((name, "".join(chunks).upper()))
    if not records: raise ValueError("no FASTA records")
    if any(not seq for _, seq in records): raise ValueError("empty probe sequence")
    if len(records) != len({name for name, _ in records}): raise ValueError("duplicate FASTA identifier")
    return records

def target_contigs(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return {r["reference_id"] for r in csv.DictReader(handle, delimiter="\t") if r.get("reference_class") == "TARGET" or r.get("group") in {"A", "B"}}

def candidate_intervals(lines, allowed):
    intervals = set()
    for line in lines:
        if not line or line.startswith("@"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11: continue
        flag, ref, pos, cigar = int(f[1]), f[2], int(f[3]), f[5]
        if flag & 4 or flag & 2048 or ref not in allowed: continue
        tags = {p[0]: p[2] for item in f[11:] if len(p := item.split(":", 2)) == 3}
        if "NM" not in tags: continue
        ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
        aligned = sum(n for n, op in ops if op in "MI=X")
        span = sum(n for n, op in ops if op in "MDN=X")
        identity = max(aligned - int(tags["NM"]), 0) / aligned if aligned else 0
        if aligned >= MIN_ALIGNED_BP and identity >= MIN_IDENTITY:
            intervals.add((ref, pos - 1, pos - 1 + span, f[0]))
    return sorted(intervals)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--probe-fasta",required=True); p.add_argument("--microbial-index"); p.add_argument("--taxonomy",required=True); p.add_argument("--candidate-bed",required=True); p.add_argument("--summary",required=True); p.add_argument("--alignment-sam"); p.add_argument("--threads",type=int,default=1); args=p.parse_args()
    probe=Path(args.probe_fasta).resolve(); records=parse_fasta(probe); allowed=target_contigs(args.taxonomy)
    if args.alignment_sam:
        with open(args.alignment_sam,encoding="utf-8") as handle: intervals=candidate_intervals(handle,allowed)
    else:
        index=args.microbial_index or str(Path(os.environ["MYCOWES_INDEX_ROOT"])/"microbial/AB_CDE_JOINT.sr.mmi")
        run=subprocess.Popen(["minimap2","-a","-x","sr","--secondary=yes","-N","100","--MD","-t",str(args.threads),index,str(probe)],stdout=subprocess.PIPE,text=True)
        assert run.stdout is not None; intervals=candidate_intervals(run.stdout,allowed)
        if run.wait(): raise SystemExit("USER_PROBE_BUILD_FAILED: minimap2")
    bed=Path(args.candidate_bed); bed.parent.mkdir(parents=True,exist_ok=True); bed.write_text("".join(f"{r}\t{s}\t{e}\t{q}\n" for r,s,e,q in intervals),encoding="utf-8")
    state="READY" if intervals else "READY_ZERO_CANDIDATE"
    summary=Path(args.summary); summary.parent.mkdir(parents=True,exist_ok=True); summary.write_text(json.dumps({"status":state,"probe_count":len(records),"candidate_regions":len(intervals),"min_aligned_bp":MIN_ALIGNED_BP,"min_identity":MIN_IDENTITY},indent=2)+"\n",encoding="utf-8")
    print(f"USER_PROBE_STATUS={state}"); print(f"CANDIDATE_REGIONS={len(intervals)}")

if __name__=="__main__": main()
