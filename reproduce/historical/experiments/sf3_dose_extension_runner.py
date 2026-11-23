#!/usr/bin/env python3
"""Add only the missing 1-20 pair dose series for four frozen SF3 challenge taxa."""
import fcntl
import importlib.util
import json
from pathlib import Path
import sys as _sys
_sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from reproduction_paths import ROOTS as _paths

ROOT = _paths.out_dir
SOURCE = Path(__file__).with_name('confirmed_branches.py')
DOSES = (1, 2, 3, 5, 10, 20)


def load():
    spec = importlib.util.spec_from_file_location('confirmed_branches', SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manifest(cb):
    taxonomy = {r['assembly_accession']: r for r in cb.read_tsv(cb.ASM_MAP)}
    rows = []
    for assembly in cb.NON_TARGET_ASMS:
        for _, (background, background_pairs) in cb.BACKGROUNDS.items():
            for seed in cb.SEEDS:
                for dose in DOSES:
                    rows.append(dict(
                        run_id=f'TRUTH_{assembly}_{background}_D{dose:03d}_S{seed}_JOINT',
                        panel_scope='SF3_ADDITIONAL_TAXA_DOSE_EXTENSION',
                        reference_mode='JOINT_TARGET_COMPETITOR',
                        truth_class='non_target_or_near_neighbor', assembly=assembly,
                        species_taxid=taxonomy[assembly]['species_taxid'], expected_target_species='',
                        background_depth=next(k for k, v in cb.BACKGROUNDS.items() if v[0] == background),
                        background_sample=background, background_pairs=background_pairs,
                        dose_pairs=dose, seed=seed, read_length=150,
                        insert_source='EMPIRICAL_FORMAL_SPIKE_FRAGMENT_SPANS',
                        replicate_type='fixed_seed_technical'))
    assert len(rows) == 360 and len({r['run_id'] for r in rows}) == 360
    original = {r['run_id'] for r in cb.read_tsv(ROOT / 'F3E_F4B_SF3_CONDITION_MANIFEST.tsv')}
    assert not original.intersection(r['run_id'] for r in rows)
    cb.write_tsv(ROOT / 'SF3_DOSE_EXTENSION_MANIFEST.tsv', rows)
    return rows


def main():
    cb = load()
    lock_handle = (ROOT / 'sf3_dose_extension.lock').open('w')
    fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    rows = manifest(cb)
    out = ROOT / 'analyses/SF3_DOSE_EXTENSION'
    out.mkdir(parents=True, exist_ok=True)
    cb.set_status('SF3_DOSE_EXTENSION', 'QUEUE', 'RUNNING', '360 missing dose conditions; one 4-thread lane')
    spans = cb.empirical_spans()
    results = []
    for index, condition in enumerate(rows, 1):
        result = cb.truth_one(condition, spans, cb.CFG)
        results.append({k: v for k, v in result.items() if k != 'qname_detail'})
        if index % 10 == 0:
            cb.set_status('SF3_DOSE_EXTENSION', 'QUEUE', 'RUNNING', f'{index}/360')
    cb.write_tsv(out / 'SF3_DOSE_EXTENSION_RESULTS.tsv', sorted(results, key=lambda r: r['run_id']))
    success = dict(status='COMPLETE', conditions=360, existing_dose50_conditions_reused=60,
                   total_additional_taxon_dose_conditions=420, scientific_commit=cb.COMMIT,
                   fragmentation_effect='NOT_RUN_MISSING_FROZEN_SHEARING_DISTRIBUTIONS', completed_utc=cb.now())
    (out / 'SUCCESS.json').write_text(json.dumps(success, indent=2) + '\n')
    cb.set_status('SF3_DOSE_EXTENSION', 'QUEUE', 'COMPLETE', '360/360; combine with existing dose50 n=60', out / 'SF3_DOSE_EXTENSION_RESULTS.tsv')


if __name__ == '__main__':
    main()
