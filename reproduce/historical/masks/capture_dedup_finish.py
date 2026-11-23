import csv,json,subprocess,os,fcntl,shutil,hashlib
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths
BASE=Path(str(_paths.project_root / "analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1"));OUT=_paths.out_dir;OLD=Path(str(_paths.project_root / "analysis/MycoWES_Fig3_Fig4_PRE77c4e85_FINAL_20260930T200113Z/02_spikein_runs/work_v2"));RT=BASE/'SF4_REAL_WES_LENGTH_10_20261013_v1/runtime'
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(x):(OUT/'MASK_FINAL_STATUS.json').write_text(json.dumps(x,indent=2)+'\n')
lock=(OUT/'mask_final.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
rows=read(OUT/'READ_ASSOCIATION_SUMMARY.tsv');assert len(rows)==1692
dest=OUT/'mask_final_results';dest.mkdir(exist_ok=True);temp=OUT/'mask_tmp';temp.mkdir(exist_ok=True)
env=BASE/'SF4_REAL_WES_LENGTH_10_20261013_v1/runtime_rescore.env';filterpy=BASE/'Fig3_Fig4_TXT_COMPLETION_EXEC_20261010_v2/runtime_confirmed_77c4e85/workflow/build_primary_strict_qnames_v0.2.py'
try:
    for i,r in enumerate(rows):
        sid=r['run_id'];kit=r['kit'];result=dest/(kit+'__'+sid+'.json')
        if result.exists():continue
        save(dict(status='RUNNING',row=i+1,total=1692,pid=os.getpid(),loadavg=os.getloadavg()))
        p=OLD/sid;baseline=p/'07_breadth/PRIMARY/07_species_dedup_dispersion_summary.tsv';n=int(r['primary_fragments_before_product_mask']);masked=int(r['product_masked_primary_fragments']);after=int(r['strict_fragments_after_product_mask'])
        if not n:before=[]
        else:
            assert baseline.is_file(),str(baseline)
            before=read(baseline)
        if not masked:post=before;mode='REUSED_IDENTICAL_PRIMARY_QNAME_SET'
        elif not after:post=[];mode='VERIFIED_EMPTY_PRODUCT_MASK'
        else:
            qlines=subprocess.check_output(['quota','-v'],text=True).splitlines()
            vals=[line.split() for line in qlines if len(line.split())>=3 and all(v.isdigit() for v in line.split()[:3])];assert vals and (int(vals[-1][2])-int(vals[-1][0]))/1024**2>=10,'STOPPED_LOW_QUOTA'
            w=temp/(kit+'__'+sid);fdir=w/'06_filter';fdir.mkdir(parents=True,exist_ok=True);bout=w/'07_breadth'
            risk=p/'05_risk'/(sid+'.fragment_risk_annotated.tsv');bam=p/'03_joint_alignment'/(sid+'.trimmed.AB_CDE_joint.qname.bam');raw=p/'04_fragment_evidence'/(sid+'.fragment_evidence.raw.tsv')
            with (OUT/'mask_final.log').open('a') as log:
                subprocess.run(['python3',str(filterpy),str(risk),str(OUT/kit/'AB_region_candidates.bed'),str(fdir)],stdout=log,stderr=log,check=True)
                keep=fdir/'STRICT_DF2.keep_qnames.txt';assert len(keep.read_text().splitlines())==after
                subprocess.run(['nice','-n','15','ionice','-c2','-n7','bash',str(RT/'workflow/run_dedup_dispersion_portable.sh'),sid,str(bam),str(raw),str(keep),str(bout),str(env)],stdout=log,stderr=log,check=True)
            post=read(bout/'07_species_dedup_dispersion_summary.tsv');mode='SAVED_EVIDENCE_PRODUCT_MASK_DEDUP_RESCORE';shutil.rmtree(w)
        result.write_text(json.dumps(dict(**r,before_species=before,after_species=post,execution=mode,baseline='PRIMARY_no_capture_mask',mapper_reexecuted=False),indent=2)+'\n')
    flat=[]
    for p in sorted(dest.glob('*.json')):
        r=json.loads(p.read_text());before={s['species']:s for s in r['before_species']};after={s['species']:s for s in r['after_species']}
        for sp in sorted(set(before)|set(after)) or ['NO_SPECIES_EVIDENCE']:
            a=before.get(sp,{});b=after.get(sp,{})
            flat.append(dict(run_id=r['run_id'],kit=r['kit'],species=sp,before_coordinate_fragments=int(a.get('coordinate_dedup_fragments',0)),after_coordinate_fragments=int(b.get('coordinate_dedup_fragments',0)),before_tier=a.get('dispersion_tier','none'),after_tier=b.get('dispersion_tier','none'),execution=r['execution']))
    with (OUT/'MASK_BEFORE_AFTER_FINAL.tsv').open('w') as f:
        w=csv.DictWriter(f,list(flat[0]),delimiter='\t');w.writeheader();w.writerows(flat)
    save(dict(status='COMPLETE',condition_product_results=1692,metric_rows=len(flat),baseline='same_PRIMARY_before_product_only_capturemask',formal_default_unchanged=True))
except Exception as e:save(dict(status='FAILED',error=str(e)));raise
