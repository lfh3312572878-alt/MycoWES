import csv,json,subprocess,os,fcntl,hashlib,time
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths
OUT=_paths.out_dir
ENV=str(_paths.project_root / "envs/mycowes_beta3_runtime/bin");env=dict(os.environ,PATH=ENV+':/usr/local/bin:/usr/bin:/bin')
def save(x): (OUT/'REGION_MAPPING_STATUS.json').write_text(json.dumps(x,indent=2)+'\n')
def headroom():
    lines=subprocess.check_output(['quota','-v'],text=True).splitlines()
    for line in reversed(lines):
        f=line.split()
        if len(f)>=3 and all(v.isdigit() for v in f[:3]):return (int(f[2])-int(f[0]))/1024**2
    raise RuntimeError('quota unavailable')
lock=(OUT/'region_worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
    with (OUT/'BED_RESOURCE_QC.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert all(r['status']=='PASS' for r in rows)
    (OUT/'TXT_MAPPING_SCOPE_CHECK.txt').write_text('User-approved BED-region sequence similarity alternative. Existing frozen builder uses minimap2 sr,secondaryyesN100,>=40alignedbp,>=90%identity. Same hg38 and frozen microbial index. Product covered-target range separate from target-regions. This mapping is region sequence similarity only, no baitTm or hybridization conclusion.\n')
    for r in rows:
        kit=r['kit'];dest=OUT/kit;dest.mkdir(exist_ok=True);bed=dest/'AB_region_candidates.bed';summary=dest/'REGION_SIMILARITY.json'
        if bed.is_file() and summary.is_file():continue
        if headroom()<10:save(dict(status='STOPPED_LOW_QUOTA',kit=kit));break
        save(dict(status='RUNNING',kit=kit,threads=2,pid=os.getpid(),loadavg=os.getloadavg()))
        named=dest/'regions.named.bed';fa=dest/'regions.hg38.fa'
        with Path(r['normalized_bed']).open() as src,named.open('w') as dst:
            for i,line in enumerate(src):dst.write(line.rstrip()+f'\t{kit}_REGION_{i+1}\n')
        start=time.time()
        with (dest/'mapping.log').open('a') as log:
            subprocess.run(['bedtools','getfasta','-fi',str(_paths.ref_root / "human/hg38.fa"),'-bed',str(named),'-nameOnly','-fo',str(fa)],env=env,stdout=log,stderr=log,check=True)
            subprocess.run(['nice','-n','15','ionice','-c2','-n7','python3',str(Path(__file__).with_name('prepare_user_probe.py')),'--probe-fasta',str(fa),'--microbial-index',str(_paths.ref_root / "MycoWES_Ref_v0.1/25_reference_qc_beta_v0.28/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi"),'--taxonomy',str(_paths.ref_root / "MycoWES_Ref_v0.1/25_reference_qc_beta_v0.28/06_metadata/AB_CDE_reference_taxonomy_lineage_map.reference_qc_beta_v0.28.tsv"),'--candidate-bed',str(bed),'--summary',str(summary),'--threads','2'],env=env,stdout=log,stderr=log,check=True)
        (dest/'PROVENANCE.json').write_text(json.dumps(dict(resource=r,region_fasta_sha256=hashlib.sha256(fa.read_bytes()).hexdigest(),wall_seconds=time.time()-start,loadavg=os.getloadavg(),mapper_threads=2,interpretation='BEDregion_sequence_similarity_not_commercial_bait'),indent=2)+'\n')
    else:save(dict(status='REGION_MAPPING_COMPLETE',products=4,next='actualreadsassociation_and_maskcomparison_pending'))
except Exception as e:save(dict(status='FAILED',error=str(e)));raise
