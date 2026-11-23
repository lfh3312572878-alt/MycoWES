import csv,json,zipfile,hashlib
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths
root=_paths.out_dir;root.mkdir(parents=True,exist_ok=True)
fai=Path(str(_paths.ref_root / "human/hg38.fa.fai"));lengths={x.split('\t')[0]:int(x.split('\t')[1]) for x in fai.read_text().splitlines()}
specs=[('Agilent_V6',str(_paths.ref_root / "capture_bed/Agilent/S07604514_hg38.zip"),'S07604514_Regions.bed','target_regions'),('Agilent_V8',str(_paths.ref_root / "capture_bed/Agilent/S33266436_hg38.zip"),'S33266436_Regions.bed','target_regions'),('IDT_xGen_v2',str(_paths.ref_root / "Fig4_external_resources_20261008_v1/normalized/IDT_xGen_Exome_v2.targets.hg38.normalized.sorted.bed"),'', 'target_regions'),('Twist_Exome_2',str(_paths.ref_root / "Fig4_external_resources_20261008_v1/normalized/Twist_Exome_2.0_v2.0.2.covered_targets.hg38.normalized.sorted.bed"),'', 'covered_targets_separate_stratum')]
rows=[];bad=[]
for kit,source,member,kind in specs:
    p=Path(source);text=zipfile.ZipFile(p).read(member).decode() if member else p.read_text();intervals=[]
    for n,line in enumerate(text.splitlines(),1):
        if not line.strip() or line.startswith(('#','track','browser')):continue
        f=line.split('\t');c,a,b=f[0],int(f[1]),int(f[2])
        if c not in lengths or a<0 or b<=a or b>lengths.get(c,0):bad.append(dict(kit=kit,line=n,contig=c,start=a,end=b,error='MISSING_CONTIG_OR_INVALID_COORDINATE'));continue
        intervals.append((c,a,b))
    dest=root/(kit+'.targets.hg38.bed');dest.write_text(''.join(f'{c}\t{a}\t{b}\n' for c,a,b in sorted(set(intervals))))
    rows.append(dict(kit=kit,type=kind,build='GRCh38',source=source,zip_member=member,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),normalized_bed=str(dest),regions=len(set(intervals)),invalid_or_missing_contigs=sum(x['kit']==kit for x in bad),status='PASS' if not any(x['kit']==kit for x in bad) else 'BLOCKED_CONTIG_QC'))
for name,data,fields in [('BED_RESOURCE_QC.tsv',rows,list(rows[0])),('MISSING_CONTIG_QC.tsv',bad,['kit','line','contig','start','end','error'])]:
    with (root/name).open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(data)
(root/'SCOPE.txt').write_text('Approved BED-region sequence similarity, not commercial oligo/hybridization. Agilent Regions and IDT targets are primary target-region stratum. Twist covered_targets retained separately, no silent equivalence to unexpanded targets. Human assemblyhg38; everycontigvalidated. No scientific threshold change. Nextsteps: reuse matching profiles or missing region mapping, actualreadassociation and maskcomparison.\n')
print(json.dumps(rows))
