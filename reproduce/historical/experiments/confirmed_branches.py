#!/usr/bin/env python3
"""Run the four user-confirmed PRE-Broad Fig3/Fig4 branches, resumably."""
import csv, datetime, fcntl, gzip, hashlib, json, os, random, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths

P = Path(str(_paths.project_root))
ROOT = _paths.out_dir
OLD = P/'analysis/MycoWES_Fig3_Fig4_PRE77c4e85_FINAL_20260930T200113Z'
SOFT = OLD/'tmp/software_77c4e85'
RUNTIME = ROOT/'runtime_confirmed_77c4e85'
ENGINE = RUNTIME/'workflow/run_one_sample.sh'
CFG = OLD/'01_spikein_design/runtime_77c4e85_4threads.env'
COMMIT = '77c4e85fd30c59a1d482b094f253356f2676597f'
ASM_DIR = P/'analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Fig3_Fig4_TXT_COMPLETION_20261010_v1/downloads/challenge_reference_candidates/assemblies'
ASM_MAP = ROOT/'TRUTH_CHALLENGE_ASSEMBLY_MAP.tsv'
SPIKE_WORK = OLD/'02_spikein_runs/work_v2'
RECOVERED = P/'analysis/PRJNA1079014_MycoWES_PRE_BroadExclusion_77c4e85_RECOVERED'
UPSTREAM62 = P/'analysis/PRJNA1079014_MycoWES_v0.3-beta2_BroadExclusion_FINAL_2ecc63b_Validation62'
DISC42 = P/'analysis/PRJNA727917_MycoWES_v0.3-beta1_REFERENCE_QC_BETA'
BASELINE_CAPTURE = P/'ref/MycoWES_Ref_v0.1/25_reference_qc_beta_v0.28/06_capture_bias/panWES.design_family_support_ge2.reference_qc_beta_v0.28.bed'
TARGET_FASTA = P/'ref/MycoWES_Ref_v0.1/25_reference_qc_beta_v0.28/01_microbe/AB_TARGET.reference_qc_beta_v0.28.fa'
MINIMAP2 = P/'envs/mycowes_beta3_runtime/bin/minimap2'
SEEDS = [20260930, 20260931, 20260932, 20260933, 20260934]
DOSES = [1, 2, 3, 5, 10, 20, 50]
TARGET_ASMS = {'GCF_000019345.1':'Ureaplasma_parvum', 'GCF_000169535.1':'Ureaplasma_urealyticum', 'GCF_040556925.1':'Mycoplasmoides_genitalium'}
NON_TARGET_ASMS = ['GCF_000018785.1','GCF_006739385.1','GCF_016728105.1','GCF_018603625.2']
BACKGROUNDS = {
    'LOW': ('SRR28084134', 4691),
    'MID': ('SRR28084105', 13935),
    'HIGH': ('SRR28084168', 127531),
}
STATUS_FIELDS = ['branch','run_id','status','updated_utc','detail','result_path']

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read_tsv(path):
    path=Path(path)
    if not path.is_file() or not path.stat().st_size: return []
    with path.open(encoding='utf-8-sig',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def write_tsv(path, rows, fields=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or (list(rows[0]) if rows else [])
    tmp=path.with_name(path.name+f'.tmp.{os.getpid()}')
    with tmp.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    tmp.replace(path)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk:=f.read(8<<20): h.update(chunk)
    return h.hexdigest()
def run(cmd, log, env=None):
    Path(log).parent.mkdir(parents=True,exist_ok=True)
    with Path(log).open('a') as h:
        h.write('COMMAND='+json.dumps(list(map(str,cmd)))+'\n'); h.flush()
        return subprocess.run(list(map(str,cmd)),stdout=h,stderr=subprocess.STDOUT,env=env).returncode
def quota_headroom_gib():
    text=subprocess.run(['quota','-v'],capture_output=True,text=True,check=True).stdout.splitlines()
    for line in reversed(text):
        x=line.split()
        if len(x)>=3 and all(y.isdigit() for y in x[:3]): return (int(x[2])-int(x[0]))/1024**2
    raise RuntimeError('cannot parse quota')
def set_status(branch, run_id, status, detail='', result=''):
    path=ROOT/'CONFIRMED_BRANCH_STATUS.tsv'; lockp=ROOT/'CONFIRMED_BRANCH_STATUS.lock'
    with lockp.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        rows=read_tsv(path); rows=[r for r in rows if (r['branch'],r['run_id'])!=(branch,run_id)]
        rows.append(dict(branch=branch,run_id=run_id,status=status,updated_utc=now(),detail=detail,result_path=str(result)))
        write_tsv(path,sorted(rows,key=lambda r:(r['branch'],r['run_id'])),STATUS_FIELDS)
    print(branch,run_id,status,detail,flush=True)
def is_complete_run(root):
    p=Path(root)/'run_manifest.tsv'
    return p.is_file() and 'status\tCOMPLETE' in p.read_text(errors='ignore')
def qbase(s): return s.split()[0].lstrip('@').removesuffix('/1').removesuffix('/2')
def fasta(path):
    out={}; name=None; seq=[]
    with Path(path).open() as h:
        for line in h:
            if line.startswith('>'):
                if name: out[name]=''.join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line.strip())
    if name: out[name]=''.join(seq).upper()
    return out
def rc(seq): return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def gz_records(path):
    with gzip.open(path,'rt') as h:
        while True:
            x=[h.readline() for _ in range(4)]
            if not x[0]: return
            if not all(x): raise RuntimeError('truncated FASTQ '+str(path))
            yield [z.rstrip('\n') for z in x]
def concat_gz(dest, *srcs):
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True); tmp=dest.with_suffix(dest.suffix+'.tmp')
    with tmp.open('wb') as o:
        for src in srcs:
            with Path(src).open('rb') as i: shutil.copyfileobj(i,o,8<<20)
    tmp.replace(dest)

def freeze():
    ROOT.mkdir(parents=True,exist_ok=True)
    if not RUNTIME.exists(): shutil.copytree(SOFT,RUNTIME,symlinks=True)
    # Engineering-only cap: the frozen scripts hard-code 16/8 threads in two places.
    for rel in ('workflow/run_posttrim_portable.sh','workflow/run_one_sample.sh','workflow/run_dedup_dispersion_portable.sh'):
        p=RUNTIME/rel; text=p.read_text(); text=text.replace('-t 16 \\\n','-t 4 \\\n').replace('-j 8 \\\n','-j 4 \\\n').replace('-@ 8 \\\n','-@ 4 \\\n'); p.write_text(text)
    decisions=[
      dict(branch='F4A_ST4',decision='FULL plus five one-at-a-time removals; cumulative all-removed is an additional endpoint',implementation='common paired inputs within each module; actual rerun versus saved-evidence rescoring labelled'),
      dict(branch='F4D',decision='28 runs',implementation='2 biological WES backgrounds x (blank + dose1 x3 technical seeds + dose50 x3 technical seeds) x 2 host references; pre-host input'),
      dict(branch='F3E_F4B_SF3',decision='PE150 analytic truth controls',implementation='empirical insert distribution; fixed seeds; qname carries assembly/contig/coordinates; correct/incorrect/unresolved attribution'),
      dict(branch='SF4',decision='all PRE-Broad Strong/Supportive; fill Weak only to total 10 if needed',implementation='fixed evidence-descending sample list before perturbation; NM step 0.02'),
    ]
    write_tsv(ROOT/'CONFIRMED_IMPLEMENTATION_DECISIONS.tsv',decisions)
    (ROOT/'CONFIRMED_RESOURCE_BUDGET.json').write_text(json.dumps(dict(max_total_threads=16,truth_lanes=2,threads_per_truth_lane=4,sf4_threads_max=8,f4d_threads_per_lane=4,f4d_lanes_after_download=2,nice=10,ionice='class2_priority7',quota_stop_new_gib=10,quota_pause_gib=8,real_WES_negative_claim='NOT_ALLOWED'),indent=2)+'\n')
    amap={r['assembly_accession']:r for r in read_tsv(ASM_MAP)}
    for asm in list(TARGET_ASMS)+NON_TARGET_ASMS:
        if asm not in amap or not (ASM_DIR/f'{asm}.fna').is_file(): raise RuntimeError('missing truth assembly '+asm)
    truth=[]
    for asm,expected in TARGET_ASMS.items():
        for depth,(bg,bgpairs) in BACKGROUNDS.items():
            for seed in SEEDS:
                for dose in DOSES:
                    rid=f'TRUTH_{asm}_{bg}_D{dose:03d}_S{seed}_JOINT'
                    truth.append(dict(run_id=rid,panel_scope='F3E/F4B/SF3',reference_mode='JOINT_TARGET_COMPETITOR',truth_class='target',assembly=asm,species_taxid=amap[asm]['species_taxid'],expected_target_species=expected,background_depth=depth,background_sample=bg,background_pairs=bgpairs,dose_pairs=dose,seed=seed,read_length=150,insert_source='EMPIRICAL_FORMAL_SPIKE_FRAGMENT_SPANS',replicate_type='fixed_seed_technical'))
    for asm in NON_TARGET_ASMS:
        for depth,(bg,bgpairs) in BACKGROUNDS.items():
            for seed in SEEDS:
                dose=50; rid=f'TRUTH_{asm}_{bg}_D050_S{seed}_JOINT'
                truth.append(dict(run_id=rid,panel_scope='F3E/F4B/SF3',reference_mode='JOINT_TARGET_COMPETITOR',truth_class='non_target_or_near_neighbor',assembly=asm,species_taxid=amap[asm]['species_taxid'],expected_target_species='',background_depth=depth,background_sample=bg,background_pairs=bgpairs,dose_pairs=dose,seed=seed,read_length=150,insert_source='EMPIRICAL_FORMAL_SPIKE_FRAGMENT_SPANS',replicate_type='fixed_seed_technical'))
    for asm in list(TARGET_ASMS)+NON_TARGET_ASMS:
        for depth,(bg,bgpairs) in BACKGROUNDS.items():
            for seed in SEEDS:
                rid=f'TRUTH_{asm}_{bg}_D050_S{seed}_TARGET_ONLY'
                truth.append(dict(run_id=rid,panel_scope='F3E/F4B/SF3',reference_mode='TARGET_ONLY_ACTUAL_RERUN',truth_class='target' if asm in TARGET_ASMS else 'non_target_or_near_neighbor',assembly=asm,species_taxid=amap[asm]['species_taxid'],expected_target_species=TARGET_ASMS.get(asm,''),background_depth=depth,background_sample=bg,background_pairs=bgpairs,dose_pairs=50,seed=seed,read_length=150,insert_source='EMPIRICAL_FORMAL_SPIKE_FRAGMENT_SPANS',replicate_type='fixed_seed_technical'))
    assert len(truth)==480 and len({r['run_id'] for r in truth})==480
    write_tsv(ROOT/'F3E_F4B_SF3_CONDITION_MANIFEST.tsv',truth)
    f4d=[]
    bgdefs=[('NA12878','SRR1611183+SRR1611184','two technical WES runs pooled; one biological DNA source'),('NA18507','SRR2775522','one WES run; one biological DNA source')]
    for sample,accessions,repnote in bgdefs:
        states=[('blank',0,'NA')]+[(f'low_S{s}',1,str(s)) for s in SEEDS[:3]]+[(f'high_S{s}',50,str(s)) for s in SEEDS[:3]]
        for condition,dose,seed in states:
            pair=f'F4D_{sample}_{condition}'
            for ref in ('GRCh38','GRCh38_T2T_CHM13v2_maskedY_rCRS'):
                f4d.append(dict(run_id=pair+'_'+ref,pair_id=pair,biological_background=sample,source_accessions=accessions,reference=ref,condition=condition,dose_pairs=dose,seed=seed,spike_assembly='GCF_040556925.1',spike_taxid='2097',independent_biological_input='YES',replicate_nature='technical spike seed' if dose else 'single blank',replicate_note=repnote,input_scope='pre_host_depletion_public_WES; not microbiologically certified negative'))
    assert len(f4d)==28 and all(sum(x['pair_id']==r['pair_id'] for x in f4d)==2 for r in f4d)
    write_tsv(ROOT/'F4D_28_RUN_MANIFEST.tsv',f4d)
    candidates=[]
    for r in read_tsv(DISC42/'PRJNA727917_MycoWES_42sample_summary.tsv'):
        tier=r['best_dispersion_tier'].lower(); ev=int(r['top_coordinate_unique'] or 0)
        if tier in ('strong','supportive','weak'): candidates.append(dict(cohort='discovery42',sample_id=r['sample_id'],formal_tier=tier,independent_evidence=ev,summary_source=str(DISC42/'PRJNA727917_MycoWES_42sample_summary.tsv')))
    s62=P/'analysis/MycoWES_Fig6_TN_validation_77c4e85_v1/01_frozen_tables/Fig6_sample_summary.tsv'
    for r in read_tsv(s62):
        tier=r['MycoWES_tier'].lower(); ev=int(r['MycoWES_independent_evidence'] or 0)
        if tier in ('strong','supportive','weak'): candidates.append(dict(cohort='validation62',sample_id=r['sample_id'],formal_tier=tier,independent_evidence=ev,summary_source=str(s62)))
    selected=[dict(r,selection_reason='all_PRE_Broad_Strong_or_Supportive') for r in candidates if r['formal_tier'] in ('strong','supportive')]
    if len(selected)<10:
        weak=sorted((r for r in candidates if r['formal_tier']=='weak'),key=lambda r:(-r['independent_evidence'],r['sample_id']))
        selected += [dict(r,selection_reason='pre_fixed_evidence_ranked_Weak_fill_to_total_10') for r in weak[:10-len(selected)]]
    assert len(selected)>=10 and len({(r['cohort'],r['sample_id']) for r in selected})==len(selected)
    write_tsv(ROOT/'SF4_SAMPLE_SELECTION_FROZEN.tsv',selected)
    params=[dict(axis='FORMAL',value='default',default='YES',interpretation='formal final call')]
    params += [dict(axis='minimum_aligned_length_bp',value=str(v),default='YES' if v==100 else 'NO',interpretation='full final-call rescoring') for v in range(60,141,10) if v!=100]
    params += [dict(axis='minimum_MAPQ',value=str(v),default='YES' if v==20 else 'NO',interpretation='full final-call rescoring') for v in range(0,41,5) if v!=20]
    params += [dict(axis='maximum_NM_rate',value=f'{v/100:.2f}',default='YES' if v==8 else 'NO',interpretation='full final-call rescoring') for v in range(2,15,2) if v!=8]
    params += [dict(axis='dispersion_bin_bp',value=str(v),default='YES' if v==250 else 'NO',interpretation='descriptive occupied-bin sensitivity; official tier not reassigned') for v in (100,500,1000)]
    write_tsv(ROOT/'SF4_PARAMETER_MANIFEST.tsv',params)
    arch=[
      ('FULL','all five modules enabled','baseline; module-specific denominator'),
      ('REMOVE_ONLY_COMPETITOR','actual target-only reruns; others retained','actual rerun'),
      ('REMOVE_ONLY_CAPTURE_RISK','bypass only DF2 capture exclusion','saved-evidence rescoring'),
      ('REMOVE_ONLY_DEDUP','retain duplicates; other modules retained','saved-evidence counterfactual'),
      ('REMOVE_ONLY_DISPERSION','retain formal evidence before dispersion tier gate','saved-evidence counterfactual'),
      ('REMOVE_ONLY_ENHANCED_HOST','GRCh38 arm paired against GRCh38+T2T full arm','actual rerun; 28-run F4D denominator'),
      ('CUMULATIVE_ALL_REMOVED','all five removals as an added endpoint','secondary endpoint; never substitutes one-at-a-time effects')]
    write_tsv(ROOT/'F4A_ST4_ARCHITECTURE_FROZEN.tsv',[dict(architecture=a,definition=b,execution=c) for a,b,c in arch])
    for branch,total in [('F3E_F4B_SF3',480),('F4D',28),('SF4',len(selected)*len(params)),('F4A_ST4',7)]: set_status(branch,'QUEUE','FROZEN',f'planned_conditions={total}')

def empirical_spans():
    out=[]
    for p in SPIKE_WORK.glob('PHD_*_D050_*/04_fragment_evidence/*.fragment_evidence.raw.tsv'):
        for r in read_tsv(p):
            try: n=int(r.get('primary_fragment_span',''))
            except: continue
            if 300<=n<=1000 and r.get('pair_support_class')=='paired_target_consistent': out.append(n)
    if not out: raise RuntimeError('no empirical formal fragment spans')
    counts={n:out.count(n) for n in sorted(set(out))}
    write_tsv(ROOT/'EMPIRICAL_INSERT_DISTRIBUTION.tsv',[dict(fragment_span=n,count=c) for n,c in counts.items()])
    return out

def truth_pair_files(asm,seed,spans):
    d=ROOT/'analyses/TRUTH_CHALLENGE/source_pairs'; d.mkdir(parents=True,exist_ok=True)
    r1=d/f'{asm}.S{seed}.PE150.R1.fastq.gz'; r2=d/f'{asm}.S{seed}.PE150.R2.fastq.gz'; meta=d/f'{asm}.S{seed}.truth.tsv'
    if r1.is_file() and r2.is_file() and meta.is_file(): return r1,r2,read_tsv(meta)
    seqs=fasta(ASM_DIR/f'{asm}.fna'); eligible={k:v for k,v in seqs.items() if len(v)>=max(spans)}
    if not eligible: raise RuntimeError('no eligible contig '+asm)
    rng=random.Random(seed+int(hashlib.sha256(asm.encode()).hexdigest()[:8],16)); rows=[]
    with gzip.open(r1,'wt',compresslevel=6) as a,gzip.open(r2,'wt',compresslevel=6) as b:
        for i in range(1,51):
            insert=rng.choice(spans); contig=rng.choice(list(eligible)); seq=eligible[contig]; start=rng.randrange(0,len(seq)-insert+1); end=start+insert
            q=f'TRUTH|asm={asm}|taxid={next(x["species_taxid"] for x in read_tsv(ASM_MAP) if x["assembly_accession"]==asm)}|contig={contig}|start0={start}|end0={end}|seed={seed}|pair={i:03d}'
            s1=seq[start:start+150]; s2=rc(seq[end-150:end]); qual='I'*150
            a.write(f'@{q}/1\n{s1}\n+\n{qual}\n'); b.write(f'@{q}/2\n{s2}\n+\n{qual}\n')
            rows.append(dict(qname=q,assembly=asm,contig=contig,start0=start,end0=end,insert_length=insert,seed=seed,pair_number=i,read_length=150,error_model='none',quality_model='fixed_Q40',source_type='analytic_simulation'))
    write_tsv(meta,rows); return r1,r2,rows

def ensure_target_only():
    d=ROOT/'analyses/TRUTH_CHALLENGE/target_only_reference'; d.mkdir(parents=True,exist_ok=True)
    mmi=d/'AB_TARGET.reference_qc_beta_v0.28.sr.mmi'; log=d/'index.log'
    if not mmi.is_file() or mmi.stat().st_size<1_000_000:
        set_status('F3E_F4B_SF3','TARGET_ONLY_INDEX','RUNNING',result=mmi)
        if run([MINIMAP2,'-x','sr','-t','4','-d',mmi,TARGET_FASTA],log): raise RuntimeError('target-only index failed')
        set_status('F3E_F4B_SF3','TARGET_ONLY_INDEX','COMPLETE',f'sha256={sha(mmi)}',mmi)
    cfg=d/'target_only.env'
    cfg.write_text(CFG.read_text()+f'\nMYCOWES_JOINT_FASTA={TARGET_FASTA}\nMYCOWES_JOINT_MMI={mmi}\nMYCOWES_MICROBIAL_PROFILE_ID=AB_TARGET_ONLY_ANALYTIC_PERTURBATION\nTHREADS=4\nSAMTOOLS_THREADS=2\n')
    return cfg

def compact_truth(runroot, condition, truth_rows):
    rid=condition['run_id']; ids={r['qname'] for r in truth_rows[:int(condition['dose_pairs'])]}
    audit=read_tsv(runroot/'04_fragment_evidence'/f'{rid}.target_candidates.audit.tsv')
    raw=read_tsv(runroot/'04_fragment_evidence'/f'{rid}.fragment_evidence.raw.tsv')
    byq={}
    for r in audit: byq.setdefault(r.get('fragment_id',''),[]).append(r)
    rawq={}
    for r in raw: rawq.setdefault(r.get('fragment_id',''),[]).append(r)
    detail=[]; correct=incorrect=unresolved=0; expected=condition['expected_target_species']
    for q in sorted(ids):
        labels={r.get('target_species_label','') for r in byq.get(q,[]) if r.get('target_species_label','')}
        assignments={r.get('initial_assignment','') for r in rawq.get(q,[])}
        if expected:
            status='correct' if expected in labels else ('incorrect' if labels else 'unresolved')
        else:
            status='incorrect' if labels else ('correct' if any(a=='competitor' for a in assignments) else 'unresolved')
        correct+=status=='correct'; incorrect+=status=='incorrect'; unresolved+=status=='unresolved'
        detail.append(dict(run_id=rid,qname=q,truth_class=condition['truth_class'],source_assembly=condition['assembly'],expected_target_species=expected,observed_target_species=';'.join(sorted(labels)),observed_assignments=';'.join(sorted(assignments)),attribution=status))
    primary=(runroot/'06_filter/PRIMARY.keep_qnames.txt').read_text().split() if (runroot/'06_filter/PRIMARY.keep_qnames.txt').is_file() else []
    strict=(runroot/'06_filter/STRICT_DF2.keep_qnames.txt').read_text().split() if (runroot/'06_filter/STRICT_DF2.keep_qnames.txt').is_file() else []
    summ=read_tsv(runroot/'07_breadth/STRICT_DF2/07_species_dedup_dispersion_summary.tsv')
    return dict(**condition,truth_pairs=len(ids),correct_attribution=correct,incorrect_attribution=incorrect,unresolved_attribution=unresolved,PRIMARY_truth_qnames=len(ids&set(primary)),STRICT_truth_qnames=len(ids&set(strict)),all_PRIMARY=len(primary),all_STRICT=len(strict),coordinate_unique=sum(int(r.get('coordinate_dedup_fragments',0) or 0) for r in summ),final_tiers=';'.join(sorted({r.get('dispersion_tier','none') for r in summ})) or 'none',species_results_json=json.dumps(summ,separators=(',',':')),qname_detail=detail)

def truth_one(condition, spans, target_cfg):
    rid=condition['run_id']; result=ROOT/'analyses/TRUTH_CHALLENGE/results'/f'{rid}.json'
    if result.is_file(): return json.loads(result.read_text())
    if quota_headroom_gib()<10: raise RuntimeError('STOPPED_LOW_QUOTA')
    set_status('F3E_F4B_SF3',rid,'RUNNING')
    sr1,sr2,trows=truth_pair_files(condition['assembly'],int(condition['seed']),spans)
    bg=condition['background_sample']; bgdir=UPSTREAM62/'work'/bg/'02_host_depletion'
    br1=bgdir/f'{bg}.trimmed.hg38_both_unmapped_R1.fastq.gz'; br2=bgdir/f'{bg}.trimmed.hg38_both_unmapped_R2.fastq.gz'
    tmp=ROOT/'analyses/TRUTH_CHALLENGE/tmp'/rid; tmp.mkdir(parents=True,exist_ok=True)
    m1=tmp/f'{rid}.R1.fastq.gz'; m2=tmp/f'{rid}.R2.fastq.gz'; dose=int(condition['dose_pairs'])
    # Write only the requested prefix of the 50 simulated pairs, then concatenate one background.
    cut1=tmp/'spike.R1.fastq.gz'; cut2=tmp/'spike.R2.fastq.gz'
    with gzip.open(cut1,'wt') as a,gzip.open(cut2,'wt') as b:
        for i,(x,y) in enumerate(zip(gz_records(sr1),gz_records(sr2))):
            if i>=dose: break
            a.write('\n'.join(x)+'\n'); b.write('\n'.join(y)+'\n')
    concat_gz(m1,br1,cut1); concat_gz(m2,br2,cut2)
    runroot=ROOT/'analyses/TRUTH_CHALLENGE/work'/rid; log=ROOT/'analyses/TRUTH_CHALLENGE/logs'/f'{rid}.log'
    cfg=target_cfg if condition['reference_mode']=='TARGET_ONLY_ACTUAL_RERUN' else CFG
    rc0=run(['bash',ENGINE,'--sample',rid,'--r1',m1,'--r2',m2,'--run-root',runroot,'--config',cfg,'--host-depleted'],log)
    if rc0 or not is_complete_run(runroot):
        set_status('F3E_F4B_SF3',rid,'FAILED',f'rc={rc0}; retained={runroot}',log); shutil.rmtree(tmp,ignore_errors=True); raise RuntimeError('truth run failed '+rid)
    payload=compact_truth(runroot,condition,trows); result.parent.mkdir(parents=True,exist_ok=True); result.write_text(json.dumps(payload,indent=2)+'\n')
    (ROOT/'analyses/TRUTH_CHALLENGE/manifests').mkdir(parents=True,exist_ok=True); shutil.copy2(runroot/'run_manifest.tsv',ROOT/'analyses/TRUTH_CHALLENGE/manifests'/f'{rid}.tsv')
    shutil.rmtree(runroot); shutil.rmtree(tmp); set_status('F3E_F4B_SF3',rid,'COMPLETE',result=result); return payload

def truth_worker():
    lock=(ROOT/'truth_worker.lock').open('w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    set_status('F3E_F4B_SF3','QUEUE','RUNNING','2 lanes x 4 threads; nice/ionice applied by launcher')
    spans=empirical_spans(); target_cfg=ensure_target_only(); conditions=read_tsv(ROOT/'F3E_F4B_SF3_CONDITION_MANIFEST.tsv')
    results=[]; failed=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        future={pool.submit(truth_one,c,spans,target_cfg):c for c in conditions}
        for f in as_completed(future):
            try: results.append(f.result())
            except Exception as e: failed.append((future[f]['run_id'],str(e)))
    flat=[]
    for r in results:
        detail=r.pop('qname_detail',[]); flat.append(r)
        p=ROOT/'analyses/TRUTH_CHALLENGE/qname_attribution'/f'{r["run_id"]}.tsv'; write_tsv(p,detail,['run_id','qname','truth_class','source_assembly','expected_target_species','observed_target_species','observed_assignments','attribution'])
    write_tsv(ROOT/'analyses/TRUTH_CHALLENGE/TRUTH_CHALLENGE_RESULTS.tsv',sorted(flat,key=lambda x:x['run_id']))
    if failed:
        write_tsv(ROOT/'analyses/TRUTH_CHALLENGE/FAILED.tsv',[dict(run_id=a,error=b) for a,b in failed]); set_status('F3E_F4B_SF3','QUEUE','PARTIAL',f'complete={len(results)}/480 failed={len(failed)}')
    else:
        (ROOT/'analyses/TRUTH_CHALLENGE/SUCCESS.json').write_text(json.dumps(dict(status='COMPLETE',conditions=480,scientific_commit=COMMIT,completed_utc=now()),indent=2)+'\n'); set_status('F3E_F4B_SF3','QUEUE','COMPLETE','480/480',ROOT/'analyses/TRUTH_CHALLENGE/SUCCESS.json')

def rewrite_core(source,dest,length=100,mapq=20,nm=.08):
    rows=read_tsv(source); fields=list(rows[0]) if rows else []
    for r in rows:
        try: ok=(r.get('pair_support_class')=='paired_target_consistent' and int(r['R1_primary_aligned_bases'])>=length and int(r['R2_primary_aligned_bases'])>=length and int(r['R1_primary_MAPQ'])>=mapq and int(r['R2_primary_MAPQ'])>=mapq and float(r['maximum_primary_NM_rate'])<=nm)
        except: ok=False
        r['exploratory_core_alignment']='1' if ok else '0'
    write_tsv(dest,rows,fields)

def extract_breadth(path):
    rows=read_tsv(Path(path)/'07_species_dedup_dispersion_summary.tsv')
    return dict(species_n=len(rows),coordinate_unique=sum(int(r.get('coordinate_dedup_fragments',0) or 0) for r in rows),occupied_bins=sum(int(r.get('occupied_250bp_bins',0) or 0) for r in rows),tier=';'.join(sorted({r.get('dispersion_tier','none') for r in rows})) or 'none',species_results_json=json.dumps(rows,separators=(',',':')))

def sf4_worker():
    lock=(ROOT/'sf4_worker.lock').open('w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    set_status('SF4','QUEUE','RUNNING','fixed selection; sequential 4-thread final-call rescoring')
    filt=RUNTIME/'workflow/build_primary_strict_qnames_v0.2.py'; breadth=RUNTIME/'workflow/run_dedup_dispersion_portable.sh'
    params=read_tsv(ROOT/'SF4_PARAMETER_MANIFEST.tsv'); results=[]
    for s in read_tsv(ROOT/'SF4_SAMPLE_SELECTION_FROZEN.tsv'):
        sid=s['sample_id']; cohort=s['cohort']
        if cohort=='discovery42':
            base=DISC42/'work'/sid; risk=base/'05_risk'/f'{sid}.fragment_risk_annotated.tsv'; bam=base/'03_joint_alignment'/f'{sid}.trimmed.AB_CDE_joint.qname.bam'; raw=base/'04_fragment_evidence'/f'{sid}.fragment_evidence.raw.tsv'; formal=base/'07_breadth/STRICT_DF2'
        else:
            base=RECOVERED/'work'/sid; risk=base/'05_risk'/f'{sid}.fragment_risk_annotated.tsv'; upstream=UPSTREAM62/'work'/sid; bam=upstream/'03_joint_alignment'/f'{sid}.trimmed.AB_CDE_joint.qname.bam'; raw=upstream/'04_fragment_evidence'/f'{sid}.fragment_evidence.raw.tsv'; formal=base/'07_breadth/STRICT_DF2'
        if not all(p.is_file() for p in (risk,bam,raw)):
            set_status('SF4',sid,'BLOCKED_MISSING_SAVED_EVIDENCE',f'risk={risk.is_file()} bam={bam.is_file()} raw={raw.is_file()}'); continue
        for p in params:
            key=f'{sid}__{p["axis"]}_{p["value"]}'; outjson=ROOT/'analyses/SF4/results'/f'{key}.json'
            if outjson.is_file(): results.append(json.loads(outjson.read_text())); continue
            if quota_headroom_gib()<10: set_status('SF4','QUEUE','STOPPED_LOW_QUOTA'); return
            set_status('SF4',key,'RUNNING'); temp=ROOT/'analyses/SF4/tmp'/key; fdir=temp/'06_filter'; bout=temp/'07_breadth'; fdir.mkdir(parents=True,exist_ok=True)
            if p['axis']=='FORMAL':
                metrics=extract_breadth(formal) if formal.is_dir() else dict(species_n=0,coordinate_unique=0,occupied_bins=0,tier='none',species_results_json='[]')
            else:
                if p['axis']=='dispersion_bin_bp':
                    keep=base/'06_filter/STRICT_DF2.keep_qnames.txt'; cfg=temp/'bin.env'; cfg.write_text(CFG.read_text()+f'\nDISPERSION_BIN_BP={p["value"]}\nTHREADS=4\nSAMTOOLS_THREADS=2\n')
                else:
                    mod=temp/'risk.tsv'; kw=dict(length=100,mapq=20,nm=.08)
                    if p['axis']=='minimum_aligned_length_bp': kw['length']=int(p['value'])
                    elif p['axis']=='minimum_MAPQ': kw['mapq']=int(p['value'])
                    elif p['axis']=='maximum_NM_rate': kw['nm']=float(p['value'])
                    rewrite_core(risk,mod,**kw); rc0=run([sys.executable,filt,mod,BASELINE_CAPTURE,fdir],ROOT/'analyses/SF4/logs'/f'{key}.log')
                    if rc0: raise RuntimeError('SF4 filter failed '+key)
                    keep=fdir/'STRICT_DF2.keep_qnames.txt'; cfg=CFG
                n=sum(1 for _ in keep.open()) if keep.is_file() else 0
                if n:
                    rc0=run(['bash',breadth,sid,bam,raw,keep,bout,cfg],ROOT/'analyses/SF4/logs'/f'{key}.log')
                    if rc0: raise RuntimeError('SF4 breadth failed '+key)
                    metrics=extract_breadth(bout)
                else: metrics=dict(species_n=0,coordinate_unique=0,occupied_bins=0,tier='none',species_results_json='[]')
                if p['axis']=='dispersion_bin_bp': metrics['tier']=s['formal_tier']; metrics['tier_interpretation']='FORMAL_250BP_REFERENCE_ONLY; nondefault bin is descriptive'
            row=dict(**s,axis=p['axis'],value=p['value'],default=p['default'],parameter_interpretation=p['interpretation'],status='COMPLETE',**metrics,scientific_commit=COMMIT,source_risk=str(risk),source_bam=str(bam),source_raw=str(raw))
            outjson.parent.mkdir(parents=True,exist_ok=True); outjson.write_text(json.dumps(row,indent=2)+'\n'); results.append(row); shutil.rmtree(temp,ignore_errors=True); set_status('SF4',key,'COMPLETE',result=outjson)
    write_tsv(ROOT/'analyses/SF4/SF4_RESULTS.tsv',sorted(results,key=lambda r:(r['cohort'],r['sample_id'],r['axis'],r['value'])))
    expected=len(read_tsv(ROOT/'SF4_SAMPLE_SELECTION_FROZEN.tsv'))*len(params)
    status='COMPLETE' if len(results)==expected else 'PARTIAL'; set_status('SF4','QUEUE',status,f'{len(results)}/{expected}',ROOT/'analyses/SF4/SF4_RESULTS.tsv')
    if status=='COMPLETE': (ROOT/'analyses/SF4/SUCCESS.json').write_text(json.dumps(dict(status='COMPLETE',conditions=expected,scientific_commit=COMMIT,completed_utc=now()),indent=2)+'\n')

def fq_record(stream):
    x=[stream.readline() for _ in range(4)]
    if not x[0]: return None
    if not all(x) or not x[0].startswith(b'@') or not x[2].startswith(b'+'): raise RuntimeError('malformed public FASTQ')
    if len(x[1].rstrip())!=len(x[3].rstrip()): raise RuntimeError('public FASTQ sequence/quality mismatch')
    return x
def open_pigz_reader(path):
    p=subprocess.Popen(['pigz','-dc','-p','2',str(path)],stdout=subprocess.PIPE); return p,p.stdout
def open_pigz_writer(path):
    h=Path(path).open('wb'); p=subprocess.Popen(['pigz','-1','-p','2','-c'],stdin=subprocess.PIPE,stdout=h); return p,p.stdin,h
def make_wes_subset(sample, sources, target_depth, source_depth):
    out=ROOT/'analyses/F4D_28_RUNS/inputs'; out.mkdir(parents=True,exist_ok=True)
    r1=out/f'{sample}.base.R1.fastq.gz'; r2=out/f'{sample}.base.R2.fastq.gz'; qc=out/f'{sample}.subset_qc.json'
    if r1.is_file() and r2.is_file() and qc.is_file(): return r1,r2,json.loads(qc.read_text())
    prob=target_depth/source_depth; rng=random.Random(20260930); selected=observed=0
    w1,o1,h1=open_pigz_writer(r1); w2,o2,h2=open_pigz_writer(r2)
    try:
        for acc,expected in sources:
            a,f1=open_pigz_reader(ROOT/'downloads/public_wes'/f'{acc}_R1.fastq.gz'); b,f2=open_pigz_reader(ROOT/'downloads/public_wes'/f'{acc}_R2.fastq.gz'); n=0
            try:
                while True:
                    x,y=fq_record(f1),fq_record(f2)
                    if x is None or y is None:
                        if x is not None or y is not None: raise RuntimeError('mate count mismatch '+acc)
                        break
                    n+=1; observed+=1
                    if qbase(x[0].decode())!=qbase(y[0].decode()): raise RuntimeError('mate qname mismatch '+acc)
                    if rng.random()>=prob: continue
                    q=qbase(x[0].decode()).encode(); x[0]=b'@WES|'+acc.encode()+b'|'+q+b'/1\n'; y[0]=b'@WES|'+acc.encode()+b'|'+q+b'/2\n'
                    o1.writelines(x); o2.writelines(y); selected+=1
            finally:
                f1.close(); f2.close()
                if a.wait() or b.wait(): raise RuntimeError('pigz decompression failed '+acc)
            if n!=expected: raise RuntimeError(f'{acc} pair count {n} != {expected}')
    finally:
        o1.close(); o2.close()
        if w1.wait() or w2.wait(): raise RuntimeError('pigz compression failed')
        h1.close(); h2.close()
    payload=dict(sample=sample,source_accessions=[a for a,_ in sources],source_pairs=observed,selection_seed=20260930,selection_probability=prob,selected_pairs=selected,target_nominal_depth=target_depth,source_nominal_depth=source_depth,biological_independence='one biological DNA source',microbiological_negative_certification='NOT_AVAILABLE')
    qc.write_text(json.dumps(payload,indent=2)+'\n'); return r1,r2,payload

def count_retained_sources(path):
    total=truth=0
    for x in gz_records(path):
        total+=1; truth+=qbase(x[0]).startswith('TRUTH|')
    return total,truth

def f4d_run_one(row, r1, r2, cfg):
    rid=row['run_id']; result=ROOT/'analyses/F4D_28_RUNS/results'/f'{rid}.json'
    if result.is_file(): return json.loads(result.read_text())
    set_status('F4D',rid,'RUNNING'); runroot=ROOT/'analyses/F4D_28_RUNS/work'/rid; log=ROOT/'analyses/F4D_28_RUNS/logs'/f'{rid}.log'
    rc0=run(['bash',ENGINE,'--sample',rid,'--r1',r1,'--r2',r2,'--run-root',runroot,'--config',cfg],log)
    if rc0 or not is_complete_run(runroot): set_status('F4D',rid,'FAILED',f'rc={rc0}; retained={runroot}',log); raise RuntimeError('F4D run failed '+rid)
    host=list((runroot/'02_host_depletion').glob('*both_unmapped_R1.fastq.gz'))
    if len(host)!=1: raise RuntimeError('F4D host-depleted output missing '+rid)
    retained,truth_retained=count_retained_sources(host[0]); dose=int(row['dose_pairs']); base_qc=json.loads((ROOT/'analyses/F4D_28_RUNS/inputs'/f'{row["biological_background"]}.subset_qc.json').read_text())
    summ=read_tsv(runroot/'07_breadth/STRICT_DF2/07_species_dedup_dispersion_summary.tsv')
    payload=dict(**row,input_background_pairs=base_qc['selected_pairs'],input_truth_pairs=dose,retained_total_pairs=retained,retained_truth_pairs=truth_retained,retained_background_pairs=retained-truth_retained,removed_truth_pairs=dose-truth_retained,removed_background_pairs=base_qc['selected_pairs']-(retained-truth_retained),truth_retention_fraction=(truth_retained/dose if dose else ''),background_removal_fraction=(1-(retained-truth_retained)/base_qc['selected_pairs'] if base_qc['selected_pairs'] else ''),STRICT_fragments=sum(int(r.get('raw_fragments',0) or 0) for r in summ),coordinate_unique=sum(int(r.get('coordinate_dedup_fragments',0) or 0) for r in summ),tier=';'.join(sorted({r.get('dispersion_tier','none') for r in summ})) or 'none',species_results_json=json.dumps(summ,separators=(',',':')),scientific_commit=COMMIT)
    result.parent.mkdir(parents=True,exist_ok=True); result.write_text(json.dumps(payload,indent=2)+'\n'); (ROOT/'analyses/F4D_28_RUNS/manifests').mkdir(parents=True,exist_ok=True); shutil.copy2(runroot/'run_manifest.tsv',ROOT/'analyses/F4D_28_RUNS/manifests'/f'{rid}.tsv'); shutil.rmtree(runroot); set_status('F4D',rid,'COMPLETE',result=result); return payload

def f4d_worker():
    lock=(ROOT/'f4d_worker.lock').open('w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    gate=ROOT/'PUBLIC_WES_DOWNLOAD_COMPLETE.json'; set_status('F4D','QUEUE','WAITING_FOR_PUBLIC_WES_QC',result=gate)
    while not gate.exists():
        if (ROOT/'PUBLIC_WES_DOWNLOAD_FAILED.json').exists(): set_status('F4D','QUEUE','BLOCKED_DOWNLOAD_FAILED'); return
        time.sleep(300)
    manifest=read_tsv(ROOT/'F4D_28_RUN_MANIFEST.tsv'); assert len(manifest)==28 and len({r['run_id'] for r in manifest})==28
    mmi=P/'ref/Fig4_external_resources_20261008_v1/T2T_CHM13v2.0/GRCh38_plus_T2T_CHM13v2.0_maskedY_rCRS.sr.mmi'
    if not mmi.is_file() or mmi.stat().st_size<10_000_000_000: set_status('F4D','QUEUE','BLOCKED_INDEX_QC'); return
    for acc in ('SRR1611183','SRR1611184','SRR2775522'):
        q=ROOT/'downloads/public_wes'/acc/'SUCCESS.json'
        if not q.is_file() or json.loads(q.read_text()).get('gzip_fastq_pair_qc')!='PASS': set_status('F4D','QUEUE','BLOCKED_INPUT_QC',acc); return
    # Avoid concurrent full reads of the same public FASTQs; continue automatically after the F3C reader.
    pilot=ROOT/'analyses/F3C_REAL_WES_PILOT/SUCCESS.json'
    set_status('F4D','QUEUE','WAITING_FOR_F3C_SUBSET_READER',f'input+index QC PASS; 28-run manifest frozen; gate={pilot}',pilot)
    while not pilot.exists():
        if (ROOT/'analyses/F3C_REAL_WES_PILOT/FAILED.json').exists(): set_status('F4D','QUEUE','BLOCKED_F3C_INPUT_PILOT_FAILED'); return
        time.sleep(300)
    capture_gate=ROOT/'analyses/F4A_ST4/CAPTURE_ONLY_RESULTS.tsv'
    set_status('F4D','QUEUE','WAITING_FOR_RESOURCE_SLOT','public WES/index/pilot PASS; waiting SF4 and capture-only lane to release CPU',capture_gate)
    while not capture_gate.exists():
        failed=next((r for r in read_tsv(ROOT/'CONFIRMED_BRANCH_STATUS.tsv') if r['branch']=='F4A_ST4' and r['run_id']=='QUEUE' and r['status'].startswith('BLOCKED')),None)
        if failed: set_status('F4D','QUEUE','BLOCKED_RESOURCE_GATE',failed['detail']); return
        time.sleep(300)
    set_status('F4D','QUEUE','GENERATING_PAIRED_BASE_INPUTS','input/index/subset pilot gates PASS')
    na12878=make_wes_subset('NA12878',[('SRR1611183',73583334),('SRR1611184',63635601)],30,241.84)
    # Nominal WES target depth from sequenced bases / frozen Agilent-V5 footprint; technical use only.
    na18507=make_wes_subset('NA18507',[('SRR2775522',81319530)],30,16263906000/50446305)
    spans=empirical_spans(); base_cfg=CFG; t2tcfg=ROOT/'analyses/F4D_28_RUNS/GRCh38_T2T.env'; t2tcfg.parent.mkdir(parents=True,exist_ok=True)
    t2tcfg.write_text(CFG.read_text()+f'\nMYCOWES_HOST_MMI={mmi}\nHOST_REFERENCE_PROFILE=GRCh38_plus_T2T_CHM13v2_maskedY_rCRS\nTHREADS=4\nSAMTOOLS_THREADS=2\n')
    results=[]
    for pair_id in sorted({r['pair_id'] for r in manifest}):
        paired=[r for r in manifest if r['pair_id']==pair_id]; first=paired[0]; sample=first['biological_background']; base=(na12878 if sample=='NA12878' else na18507)
        temp=ROOT/'analyses/F4D_28_RUNS/tmp'/pair_id; temp.mkdir(parents=True,exist_ok=True); m1=temp/f'{pair_id}.R1.fastq.gz'; m2=temp/f'{pair_id}.R2.fastq.gz'; dose=int(first['dose_pairs'])
        if dose:
            sr1,sr2,_=truth_pair_files(first['spike_assembly'],int(first['seed']),spans); c1=temp/'spike.R1.fastq.gz'; c2=temp/'spike.R2.fastq.gz'
            with gzip.open(c1,'wt') as a,gzip.open(c2,'wt') as b:
                for i,(x,y) in enumerate(zip(gz_records(sr1),gz_records(sr2))):
                    if i>=dose: break
                    a.write('\n'.join(x)+'\n'); b.write('\n'.join(y)+'\n')
            concat_gz(m1,base[0],c1); concat_gz(m2,base[1],c2)
        else: concat_gz(m1,base[0]); concat_gz(m2,base[1])
        with ThreadPoolExecutor(max_workers=2) as pool:
            jobs=[pool.submit(f4d_run_one,r,m1,m2,t2tcfg if r['reference'].startswith('GRCh38_T2T') else base_cfg) for r in paired]
            for f in jobs: results.append(f.result())
        shutil.rmtree(temp)
    write_tsv(ROOT/'analyses/F4D_28_RUNS/F4D_RESULTS.tsv',sorted(results,key=lambda r:r['run_id']))
    for p in (na12878[0],na12878[1],na18507[0],na18507[1]): p.unlink(missing_ok=True)
    (ROOT/'analyses/F4D_28_RUNS/SUCCESS.json').write_text(json.dumps(dict(status='COMPLETE',runs=28,biological_backgrounds=2,technical_pairings=14,scientific_commit=COMMIT,completed_utc=now()),indent=2)+'\n')
    set_status('F4D','QUEUE','COMPLETE','28/28; 14 paired contrasts; 2 biological backgrounds',ROOT/'analyses/F4D_28_RUNS/F4D_RESULTS.tsv')

def ablation_worker():
    lock=(ROOT/'ablation_worker.lock').open('w'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    set_status('F4A_ST4','QUEUE','RUNNING','freeze and reuse audit; capture-only rescoring scheduled')
    src=P/'analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Figure4_ABLATION_READY_jghpus0z/module_results.tsv'
    clinical=P/'analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Figure4_ABLATION_READY_jghpus0z/existing_reference_contrasts.tsv'
    rows=[]
    for r in read_tsv(src): rows.append(dict(source='existing_PRE_Broad_ablation',architecture=r.get('mode',''),execution='formal_or_saved_evidence_counterfactual',run_id=r.get('run_id',''),status='COMPLETE',source_file=str(src)))
    for r in read_tsv(clinical): rows.append(dict(source='existing_six_clinical_target_only',architecture='REMOVE_ONLY_COMPETITOR',execution='actual_rerun_without_independent_source_truth',run_id=r.get('run_id',''),status='COMPLETE_UNRESOLVED_TRUTH',source_file=str(clinical)))
    write_tsv(ROOT/'analyses/F4A_ST4/REUSED_ENDPOINT_INVENTORY.tsv',rows)
    pending=[
      dict(architecture='REMOVE_ONLY_COMPETITOR',dependency='F3E_F4B_SF3 actual truth reruns',status='RUNNING'),
      dict(architecture='REMOVE_ONLY_CAPTURE_RISK',dependency='423 saved-evidence capture-only final-call rescoring',status='SCHEDULED'),
      dict(architecture='REMOVE_ONLY_DEDUP',dependency=str(src),status='COMPLETE_REUSED'),
      dict(architecture='REMOVE_ONLY_DISPERSION',dependency=str(src),status='COMPLETE_REUSED'),
      dict(architecture='REMOVE_ONLY_ENHANCED_HOST',dependency='F4D 28 actual paired runs',status='WAITING'),
      dict(architecture='CUMULATIVE_ALL_REMOVED',dependency='all five branches; secondary endpoint',status='WAITING')]
    write_tsv(ROOT/'analyses/F4A_ST4/ARCHITECTURE_STATUS.tsv',pending)
    # Stay within the 12-thread budget: capture-only rescoring begins after SF4 releases its 4 threads.
    set_status('F4A_ST4','QUEUE','WAITING_FOR_SF4','dedup/dispersion reused; competitor and enhanced-host workers active',ROOT/'analyses/F4A_ST4/ARCHITECTURE_STATUS.tsv')
    while not (ROOT/'analyses/SF4/SUCCESS.json').exists():
        failed=next((r for r in read_tsv(ROOT/'CONFIRMED_BRANCH_STATUS.tsv') if r['branch']=='SF4' and r['run_id']=='QUEUE' and r['status'] in ('FAILED','STOPPED_LOW_QUOTA')),None)
        if failed: set_status('F4A_ST4','QUEUE','BLOCKED_SF4_GATE',failed['detail']); return
        time.sleep(300)
    set_status('F4A_ST4','QUEUE','RUNNING_CAPTURE_ONLY_RESCORE','423 saved-evidence paired conditions')
    breadth=RUNTIME/'workflow/run_dedup_dispersion_portable.sh'; capture_results=[]
    cases=read_tsv(P/'analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1/Figure4_ABLATION_READY_jghpus0z/cases.tsv')
    cases=[r for r in cases if r.get('scope')=='spikein']
    for i,c in enumerate(cases,1):
        rid=c['run_id']; outjson=ROOT/'analyses/F4A_ST4/capture_only_results'/f'{rid}.json'
        if outjson.is_file(): capture_results.append(json.loads(outjson.read_text())); continue
        base=SPIKE_WORK/rid; keep=base/'06_filter/PRIMARY.keep_qnames.txt'; bam=base/'03_joint_alignment'/f'{rid}.trimmed.AB_CDE_joint.qname.bam'; raw=base/'04_fragment_evidence'/f'{rid}.fragment_evidence.raw.tsv'; n=sum(1 for _ in keep.open()) if keep.is_file() else 0
        if n:
            temp=ROOT/'analyses/F4A_ST4/tmp'/rid; rc0=run(['bash',breadth,rid,bam,raw,keep,temp,CFG],ROOT/'analyses/F4A_ST4/logs'/f'{rid}.log')
            if rc0: raise RuntimeError('capture-only breadth failed '+rid)
            metrics=extract_breadth(temp); shutil.rmtree(temp)
        else: metrics=dict(species_n=0,coordinate_unique=0,occupied_bins=0,tier='none',species_results_json='[]')
        row=dict(run_id=rid,architecture='REMOVE_ONLY_CAPTURE_RISK',execution='saved-evidence counterfactual; PRIMARY qnames passed to full dedup/dispersion',input_fragments=n,status='COMPLETE',**metrics,scientific_commit=COMMIT)
        outjson.parent.mkdir(parents=True,exist_ok=True); outjson.write_text(json.dumps(row,indent=2)+'\n'); capture_results.append(row)
        if i%25==0: set_status('F4A_ST4','CAPTURE_ONLY_PROGRESS','RUNNING',f'{i}/{len(cases)}')
    write_tsv(ROOT/'analyses/F4A_ST4/CAPTURE_ONLY_RESULTS.tsv',sorted(capture_results,key=lambda r:r['run_id']))
    set_status('F4A_ST4','QUEUE','WAITING_DEPENDENCIES','capture-only 423 complete; waiting truth competitor and F4D host arms',ROOT/'analyses/F4A_ST4/CAPTURE_ONLY_RESULTS.tsv')
    while not ((ROOT/'analyses/TRUTH_CHALLENGE/SUCCESS.json').exists() and (ROOT/'analyses/F4D_28_RUNS/SUCCESS.json').exists()): time.sleep(300)
    pending=[dict(r,status=('COMPLETE' if r['architecture'] in ('REMOVE_ONLY_COMPETITOR','REMOVE_ONLY_CAPTURE_RISK','REMOVE_ONLY_DEDUP','REMOVE_ONLY_DISPERSION','REMOVE_ONLY_ENHANCED_HOST') else 'READY_FOR_SECONDARY_ENDPOINT_ASSEMBLY')) for r in pending]
    write_tsv(ROOT/'analyses/F4A_ST4/ARCHITECTURE_STATUS.tsv',pending)
    set_status('F4A_ST4','QUEUE','PARTIAL','five one-at-a-time endpoints available; cumulative all-removed endpoint requires final paired assembly',ROOT/'analyses/F4A_ST4/ARCHITECTURE_STATUS.tsv')

def refresh():
    checks=[
      ('F3E_F4B_SF3','truth_worker.pid','RUNNING','truth 480-condition queue active'),
      ('SF4','sf4_worker.pid','RUNNING','fixed 10-sample perturbation queue active'),
      ('F4D','f4d_worker.pid','WAITING_FOR_PUBLIC_WES_QC' if not (ROOT/'PUBLIC_WES_DOWNLOAD_COMPLETE.json').exists() else 'RUNNING','28-run worker active'),
      ('F4A_ST4','ablation_worker.pid','WAITING_FOR_SF4' if not (ROOT/'analyses/SF4/SUCCESS.json').exists() else 'RUNNING','one-at-a-time ablation worker active')]
    for branch,pidfile,status,detail in checks:
        p=ROOT/pidfile
        alive=p.is_file() and p.read_text().strip().isdigit() and Path('/proc',p.read_text().strip()).exists()
        set_status(branch,'QUEUE',status if alive else 'EXITED_NEEDS_REVIEW',detail if alive else f'pidfile={p}')

def main():
    if len(sys.argv)!=2 or sys.argv[1] not in ('freeze','truth','sf4','f4d','ablation','refresh'): raise SystemExit('usage: confirmed_branches.py freeze|truth|sf4|f4d|ablation|refresh')
    {'freeze':freeze,'truth':truth_worker,'sf4':sf4_worker,'f4d':f4d_worker,'ablation':ablation_worker,'refresh':refresh}[sys.argv[1]]()

if __name__=='__main__':
    try: main()
    except BlockingIOError: raise SystemExit('worker already active')
    except Exception as exc:
        branch=sys.argv[1] if len(sys.argv)>1 else 'UNKNOWN'; set_status(branch,'QUEUE','FAILED',repr(exc)); raise
