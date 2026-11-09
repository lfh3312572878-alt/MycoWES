import csv,json,hashlib
from pathlib import Path
from collections import defaultdict
BASE=Path('/data/h011/project/analysis/MycoWES_MANUSCRIPT_PRE_BROAD_v1');OUT=BASE/'CAPTURE_BED_REGION_20261013_v1';OLD=Path('/data/h011/project/analysis/MycoWES_Fig3_Fig4_PRE77c4e85_FINAL_20260930T200113Z/02_spikein_runs/work_v2')
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    with p.open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
conditions=read(BASE/'Batch2_analysis_20261003T173758_296088/conditions_423.tsv');assert len(conditions)==423
profiles={}
for kit in ['Agilent_V6','Agilent_V8','IDT_xGen_v2','Twist_Exome_2']:
    db=defaultdict(list)
    for line in (OUT/kit/'AB_region_candidates.bed').open():
        f=line.rstrip().split('\t');db[f[0]].append((int(f[1]),int(f[2]),f[3]))
    profiles[kit]=db
missing=[];summary=[];detail=[]
for c in conditions:
    sid=c['run_id'];p=OLD/sid/'05_risk'/f'{sid}.fragment_risk_annotated.tsv'
    if not p.is_file():missing.append(dict(run_id=sid,missing=str(p)));continue
    rows=read(p)
    for kit,db in profiles.items():
        allhit=set();primary=set();hitprimary=set()
        for r in rows:
            q=r['fragment_id'];isprimary=r.get('initial_assignment')=='target_candidate' and r.get('exploratory_core_alignment')=='1' and r.get('risk_any','') in ['','0'] and r.get('adapter_core_in_either','') in ['','0']
            if isprimary:primary.add(q)
            hits=[]
            for mate in ['R1','R2']:
                try:a=int(r[mate+'_reference_start_0based']);b=int(r[mate+'_reference_end_0based_exclusive'])
                except (ValueError,KeyError):continue
                for x,y,label in db[r.get(mate+'_primary_reference','')]:
                    n=max(0,min(b,y)-max(a,x))
                    if n:hits.append((mate,n,label))
            if hits:
                allhit.add(q)
                if isprimary:hitprimary.add(q)
                detail.append(dict(run_id=sid,kit=kit,fragment_id=q,target_species=r.get('target_species_label',''),primary_quality_risk_adapter_pass=isprimary,region_overlap_evidence=json.dumps(hits),interpretation='observed microbial alignment overlaps region-similarity candidate; no truth or hybridization claim'))
        summary.append(dict(run_id=sid,biological_background=c['background_sample'],dose_pairs=c['dose_pairs'],seed=c['seed'],kit=kit,observed_candidate_fragments=len(rows),region_overlapping_fragments=len(allhit),primary_fragments_before_product_mask=len(primary),product_masked_primary_fragments=len(hitprimary),strict_fragments_after_product_mask=len(primary-hitprimary),metric_scope='pre-dedup qnames, not final coordinate tier; product-only mask follows existing capture comparison'))
write(OUT/'READ_ASSOCIATION_SUMMARY.tsv',summary,['run_id','biological_background','dose_pairs','seed','kit','observed_candidate_fragments','region_overlapping_fragments','primary_fragments_before_product_mask','product_masked_primary_fragments','strict_fragments_after_product_mask','metric_scope'])
write(OUT/'READ_ASSOCIATION_DETAIL.tsv',detail,['run_id','kit','fragment_id','target_species','primary_quality_risk_adapter_pass','region_overlap_evidence','interpretation'])
write(OUT/'READ_ASSOCIATION_MISSING.tsv',missing,['run_id','missing'])
(OUT/'READ_ASSOCIATION_QC.json').write_text(json.dumps(dict(status='PASS' if not missing else 'PARTIAL',planned_conditions=423,complete_condition_product_rows=len(summary),missing_conditions=len(missing),products=4,final_dedup_tier_comparison='PENDING'),indent=2)+'\n')
print('rows',len(summary),'missing',len(missing),'linkedrecords',len(detail))
