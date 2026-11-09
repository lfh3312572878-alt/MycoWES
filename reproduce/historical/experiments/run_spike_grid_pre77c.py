#!/usr/bin/env python3
import csv, gzip, hashlib, os, random, shutil, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

OUT=Path('/data/h011/project/analysis/MycoWES_Fig3_Fig4_PRE77c4e85_FINAL_20260930T200113Z')
SOFT=OUT/'tmp/software_77c4e85'
UP=Path('/data/h011/project/analysis/PRJNA1079014_MycoWES_v0.3-beta2_BroadExclusion_FINAL_2ecc63b_ALL69')
SEEDS=[20260930,20260931,20260932,20260933,20260934]
DOSES=[0,1,2,3,5,10,20,50]
THREADS=8
SPECIES={
 'M_arginini':('Mycoplasmopsis_arginini','ERR14129335','GCF_965278885.1',Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/arginini_replacement_v2/raw_WGS/ERR14129335_1.fastq.gz'),Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/arginini_replacement_v2/raw_WGS/ERR14129335_2.fastq.gz')),
 'M_hominis':('Metamycoplasma_hominis','SRR1655810','NOT_AVAILABLE',Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR1655810_1.fastq.gz'),Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR1655810_2.fastq.gz')),
 'M_hyorhinis':('Mesomycoplasma_hyorhinis','SRR15402804','NOT_AVAILABLE',Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR15402804_1.fastq.gz'),Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR15402804_2.fastq.gz')),
 'M_salivarium':('Metamycoplasma_salivarium','SRR35637520','NOT_AVAILABLE',Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR35637520_1.fastq.gz'),Path('/data/h011/project/analysis/MycoWES_external_spikein_v1/raw_WGS/SRR35637520_2.fastq.gz')),
}

def read_tsv(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,rows,fields):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fastq_records(path):
 with gzip.open(path,'rt') as f:
  while True:
   a=f.readline()
   if not a:return
   b=f.readline();c=f.readline();d=f.readline()
   if not d:raise ValueError(f'truncated FASTQ {path}')
   yield a.rstrip(),b.rstrip(),c.rstrip(),d.rstrip()
def qbase(header):return header.split()[0].lstrip('@').removesuffix('/1').removesuffix('/2')
def count_pairs(path):
 n=0
 with gzip.open(path,'rb') as f:
  for n,_ in enumerate(f,1):pass
 if n%4:raise ValueError(path)
 return n//4
def gzwrite(path, records):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as z:
   for rec in records:z.write(('\n'.join(rec)+'\n').encode())
def combine_gz(dest,*sources):
 dest.parent.mkdir(parents=True,exist_ok=True)
 tmp=dest.with_suffix(dest.suffix+'.tmp')
 with tmp.open('wb') as w:
  for p in sources:
   with p.open('rb') as r:shutil.copyfileobj(r,w,8<<20)
 os.replace(tmp,dest)

sel=read_tsv(OUT/'01_spikein_design/sample_selection_frozen.tsv')
backgrounds=[]
for x in sel:
 if x['analysis_role'].startswith('SPIKE_BACKGROUND_'):
  sid=x['sample_id']; base=UP/'work'/sid/'02_host_depletion'
  r1=base/f'{sid}.trimmed.hg38_both_unmapped_R1.fastq.gz'; r2=base/f'{sid}.trimmed.hg38_both_unmapped_R2.fastq.gz'
  n=count_pairs(r1)
  if count_pairs(r2)!=n: raise ValueError(f'background pair mismatch {sid}')
  backgrounds.append((x['analysis_role'].rsplit('_',1)[-1],sid,r1,r2,n))
if len(backgrounds)!=3:raise SystemExit('expected three frozen spike backgrounds')

# Runtime config: identical scientific values, only absolute existing tools and lower engineering concurrency.
cfg=OUT/'01_spikein_design/runtime_77c4e85_4threads.env'
base=(SOFT/'config/default.env').read_text()
cfg.write_text(base+'\nCUTADAPT=/data/h011/project/envs/mycowes_beta3_runtime/bin/cutadapt\nMINIMAP2=/data/h011/project/envs/mycowes_beta3_runtime/bin/minimap2\nSAMTOOLS=/usr/local/bin/samtools\nPYTHON=/usr/bin/python3\nMYCOWES_HOST_MMI=/data/h011/project/ref/human/hg38.mmi\nTHREADS=4\nSAMTOOLS_THREADS=2\n')

# One scan per source creates five independent fixed-seed selections of up to 50 true pairs.
spike_dir=OUT/'01_spikein_design/spike_pairs_by_background_v2'
selected_records={}
for short,(formal,run,assembly,r1,r2) in SPECIES.items():
 n=count_pairs(r1)
 if count_pairs(r2)!=n:raise ValueError(f'pair mismatch {run}')
 selected={seed:set(random.Random(seed).sample(range(n),50)) for seed in SEEDS}
 store={seed:[] for seed in SEEDS}
 for i,(a,b) in enumerate(zip(fastq_records(r1),fastq_records(r2))):
  if qbase(a[0])!=qbase(b[0]):raise ValueError(f'qname mismatch {run} {i}')
  for seed in SEEDS:
   if i in selected[seed]:store[seed].append((i,a,b))
 for seed in SEEDS:
  store[seed].sort(key=lambda x:x[0]); selected_records[(short,seed)]=store[seed]

design=[]
truth=[]
mixed=OUT/'01_spikein_design/mixed_host_depleted_fastq_v2'
for depth,bg,br1,br2,bg_pairs in backgrounds:
 for short,(formal,run,assembly,_,__) in SPECIES.items():
  for seed in SEEDS:
   for dose in DOSES:
    sid=f'PHD_{bg}_{short}_D{dose:03d}_S{seed}'
    m1=mixed/f'{sid}_R1.fastq.gz';m2=mixed/f'{sid}_R2.fastq.gz'
    src1=[];src2=[]
    if dose:
     s1=spike_dir/f'{bg}.{short}.S{seed}.D{dose:03d}.R1.fastq.gz';s2=spike_dir/f'{bg}.{short}.S{seed}.D{dose:03d}.R2.fastq.gz';rec1=[];rec2=[]
     for j,(idx,a,b) in enumerate(selected_records[(short,seed)][:dose],1):
      q=f'SPK|{short}|{run}|BG{bg}|S{seed}|D{dose:03d}|P{j:03d}'
      rec1.append((f'@{q}',a[1],a[2],a[3]));rec2.append((f'@{q}',b[1],b[2],b[3]))
      truth.append(dict(spike_qname=q,species=formal,source_run=run,source_assembly=assembly,background_sample=bg,seed=seed,dose_pairs=dose,pair_number=j,source_pair_index_0based=idx,source_qname=qbase(a[0]),r1_length=len(a[1]),r2_length=len(b[1])))
     gzwrite(s1,rec1);gzwrite(s2,rec2);src1=[s1];src2=[s2]
    if not m1.exists():combine_gz(m1,br1,*src1)
    if not m2.exists():combine_gz(m2,br2,*src2)
    design.append(dict(run_id=sid,analysis='post_host_depletion_analytical_spikein',background_depth_group=depth,background_sample=bg,species=formal,species_short=short,source_run=run,source_assembly=assembly,dose_pairs=dose,seed=seed,r1=str(m1),r2=str(m2),background_pairs=bg_pairs,expected_input_pairs=bg_pairs+dose,scientific_commit='77c4e85fd30c59a1d482b094f253356f2676597f'))
write_tsv(OUT/'01_spikein_design/spike_pair_truth_v2.tsv',truth,list(truth[0]))
write_tsv(OUT/'01_spikein_design/spikein_design_v2.tsv',design,list(design[0]))

engine=SOFT/'workflow/run_one_sample.sh'
runs=OUT/'02_spikein_runs/work_v2';logs=OUT/'02_spikein_runs/logs_v2';runs.mkdir(parents=True,exist_ok=True);logs.mkdir(parents=True,exist_ok=True)
def run_one(x):
 root=runs/x['run_id']; manifest=root/'run_manifest.tsv'; log=logs/f"{x['run_id']}.log"
 if manifest.is_file() and 'status\tCOMPLETE' in manifest.read_text(errors='ignore'):return x['run_id'],'COMPLETE_RESUMED',0
 t=time.time()
 with log.open('a') as f:
  p=subprocess.run(['bash',str(engine),'--sample',x['run_id'],'--r1',x['r1'],'--r2',x['r2'],'--run-root',str(root),'--config',str(cfg),'--host-depleted'],stdout=f,stderr=subprocess.STDOUT)
 ok=p.returncode==0 and manifest.is_file() and 'status\tCOMPLETE' in manifest.read_text(errors='ignore')
 return x['run_id'],'COMPLETE' if ok else f'FAILED_RC_{p.returncode}',round(time.time()-t,3)

status=[]
with ThreadPoolExecutor(max_workers=THREADS) as ex:
 fut={ex.submit(run_one,x):x for x in design}
 for i,f in enumerate(as_completed(fut),1):
  rid,st,sec=f.result();status.append(dict(run_id=rid,status=st,elapsed_seconds=sec));print(i,len(design),rid,st,flush=True)
write_tsv(OUT/'02_spikein_runs/run_status_v2.tsv',sorted(status,key=lambda x:x['run_id']),['run_id','status','elapsed_seconds'])

def truth_ids(x):return {r['spike_qname'] for r in truth if r['species']==x['species'] and r['background_sample']==x['background_sample'] and int(r['seed'])==int(x['seed']) and int(r['dose_pairs'])==int(x['dose_pairs'])}
def true(v):return str(v).lower() in {'1','true','t','yes','y'}
metrics=[]
for x in design:
 root=runs/x['run_id']; ids=truth_ids(x)
 auditp=root/'04_fragment_evidence'/f"{x['run_id']}.target_candidates.audit.tsv"
 riskp=root/'05_risk'/f"{x['run_id']}.fragment_risk_annotated.tsv"
 rawp=root/'04_fragment_evidence'/f"{x['run_id']}.fragment_evidence.raw.tsv"
 primp=root/'06_filter/PRIMARY.keep_qnames.txt'; strictp=root/'06_filter/STRICT_DF2.keep_qnames.txt'
 audit=read_tsv(auditp) if auditp.is_file() else [];risk=read_tsv(riskp) if riskp.is_file() else [];raw=read_tsv(rawp) if rawp.is_file() else []
 primary=set(primp.read_text().split()) if primp.is_file() else set();strict=set(strictp.read_text().split()) if strictp.is_file() else set()
 same=[r for r in audit if r.get('target_species_label')==x['species']]
 sp_same=[r for r in same if r.get('fragment_id') in ids]
 sp_core=[r for r in sp_same if true(r.get('exploratory_core_alignment'))]
 risk_by={r.get('fragment_id'):r for r in risk}
 sp_rf=[r for r in sp_core if not true(risk_by.get(r.get('fragment_id'),{}).get('risk_any')) and not true(r.get('adapter_core_in_either'))]
 err={r.get('target_species_label') for r in audit if r.get('fragment_id') in ids and r.get('target_species_label') and r.get('target_species_label')!=x['species']}
 summ=root/'07_breadth/STRICT_DF2/07_species_dedup_dispersion_summary.tsv'; sr=read_tsv(summ) if summ.is_file() else []
 official=[r for r in sr if r.get('species')==x['species']]
 o=official[0] if official else {}
 metrics.append({**x,'run_status':next((s['status'] for s in status if s['run_id']==x['run_id']),'NOT_RUN'),'true_spiked_pairs':x['dose_pairs'],'spike_recovered_any_microbial':sum(1 for r in raw if r.get('fragment_id') in ids and r.get('initial_assignment')!='no_microbial_alignment'),'spike_same_species_candidates':len(sp_same),'spike_alignment_qc':len(sp_core),'spike_risk_free':len(sp_rf),'spike_PRIMARY':len(ids&primary),'spike_STRICT':len(ids&strict),'background_STRICT':len(strict-ids),'all_same_species_candidates':len(same),'all_PRIMARY':len(primary),'all_STRICT':len(strict),'erroneous_target_species_assignments':','.join(sorted(err)) if err else 'NONE','coordinate_unique_evidence':o.get('coordinate_dedup_fragments','0') if official else '0','occupied_250bp_bins':o.get('occupied_250bp_bins','0') if official else '0','covered_bp':o.get('covered_bp','0') if official else '0','official_final_call':'CALL' if official else 'NO_CALL','official_tier':o.get('dispersion_tier','none') if official else 'none'})
write_tsv(OUT/'03_spikein_analysis/spikein_run_metrics_v2.tsv',metrics,list(metrics[0]))
print('SPIKE_GRID_COMPLETE',sum(1 for s in status if s['status'].startswith('COMPLETE')),len(status))
