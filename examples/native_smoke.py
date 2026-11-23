#!/usr/bin/env python3
"""Generate a tiny synthetic example and run the actual Linux analysis tools."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SEED = 20261109
TARGET = 'toy_target|role=A|species=Toy_target|asm=TOY_A'
COMPETITOR = 'toy_competitor|role=C|species=Toy_competitor|asm=TOY_C'


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def random_sequence(rng, size=12000):
    return ''.join(rng.choice('ACGT') for _ in range(size))


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join('>' + name + '\n' + seq + '\n' for name, seq in records))


def write_pairs(directory, fragments):
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / 'reads_R1.fastq.gz', directory / 'reads_R2.fastq.gz']
    for mate, path in enumerate(paths):
        data = ''.join('@' + name + '\n' + pair[mate] + '\n+\n' + 'I' * len(pair[mate]) + '\n'
                       for name, pair in fragments)
        with path.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as handle:
                handle.write(data.encode('ascii'))
    return paths


def pair_at(seq, start):
    return seq[start:start + 150], revcomp(seq[start + 150:start + 300])


def make_reference(directory, target_seq, competitor_seq, anchor_seq):
    """Use the resolver's directory contract; all sequences are synthetic."""
    single = directory / 'ONE_REP_TARGET_DEVELOPMENT'
    baseline = directory / 'baseline_v028'
    microbe = baseline / '01_microbe'
    metadata = baseline / '06_metadata'
    metadata.mkdir(parents=True, exist_ok=True)
    records = [(TARGET, target_seq), (COMPETITOR, competitor_seq)]
    write_fasta(single / 'ONE_REP_TARGET.AB_CDE_JOINT.fa', records)
    write_fasta(microbe / 'AB_CDE_JOINT.reference_qc_beta_v0.28.fa', records)
    write_fasta(microbe / 'AB_TARGET.reference_qc_beta_v0.28.fa', [(TARGET, anchor_seq)])
    taxonomy = 'reference_id\treference_species\tspecies_taxid\treference_class\trole\n'
    taxonomy += TARGET + '\tToy_target\t9000001\tTARGET\tA\n'
    taxonomy += COMPETITOR + '\tToy_competitor\t9000002\tCOMPETITOR\tC\n'
    (single / 'ONE_REP_TARGET.taxonomy.tsv').write_text(taxonomy)
    (metadata / 'AB_CDE_reference_taxonomy_lineage_map.reference_qc_beta_v0.28.tsv').write_text(taxonomy)
    identity = ('resource_origin\tsynthetic_fixture_seed_' + str(SEED) + '\n'
                'reference_release\tSYNTHETIC_EXAMPLE_NOT_V028_DATABASE\n'
                'reference_profile\tSYNTHETIC_EXAMPLE\n')
    (single / 'ONE_REP_TARGET.identity.tsv').write_text(
        'profile_id\tMycoWES_AB_CDE_v0.28_ONE_REP_TARGET_DEV_20261015\n' + identity)
    (metadata / 'identity.tsv').write_text(
        'profile_id\tMycoWES_AB_CDE_v0.28_reference_qc_beta\n' + identity)
    for path in [single / 'ONE_REP_TARGET.exclusions.tsv', metadata / 'reference_exclusion_manifest.tsv']:
        path.write_text('assembly\treason\n')
    risk = baseline / '05_annotations/AB_reference_risk_annotations.reference_qc_beta_v0.28.bed'
    risk.parent.mkdir(parents=True, exist_ok=True)
    risk.write_text(TARGET + '\t4000\t4010\tsynthetic_shared_risk\n')
    df2 = baseline / '06_capture_bias/panWES.design_family_support_ge2.reference_qc_beta_v0.28.bed'
    df2.parent.mkdir(parents=True, exist_ok=True)
    df2.write_text(TARGET + '\t5000\t5010\tsynthetic_DF2\n')
    (directory / 'SYNTHETIC_ONLY.txt').write_text(
        'Generated random DNA. No human, patient, microbial database or sequencing data.\n'
        'Internal profile IDs satisfy the resolver contract; they do not certify v0.28 content.\n')
    return single / 'ONE_REP_TARGET.AB_CDE_JOINT.fa'


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def run_command(command, log, env=None):
    with log.open('w') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise RuntimeError(f'Command exited {result.returncode}; see {log}: {command}')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect_case(output):
    filters = {r['metric']: int(r['fragments']) for r in table(output / '06_filter/filter_summary.tsv')}
    branches = {}
    for branch in ('PRIMARY', 'STRICT_DF2'):
        path = output / '07_breadth' / branch / '07_species_dedup_dispersion_summary.tsv'
        branches[branch] = table(path) if path.is_file() else []
    return {'filter_counts': filters, 'branches': branches}


def comparable_case(observed):
    return {'filter_counts':observed['filter_counts'], 'branches':{
        branch:[{key:value for key,value in row.items() if not Path(value).is_absolute()} for row in rows]
        for branch,rows in observed['branches'].items()}}


def managed_quality(output, preset='standard'):
    selection=json.loads((output/'reference_selection.json').read_text())
    manifest={r['key']:r['value'] for r in table(output/'run_manifest.tsv')}
    want={'min_aligned_bp':100,'min_mapq':20,'max_nm_rate':.08}
    assert selection['preset']==preset and selection['quality_thresholds']==want
    assert selection['effective_min_aligned_bp']==100 and set(selection['config_files_sha256'])=={'default.env','relax2.env'}
    assert (int(manifest['effective_min_aligned_bp']),int(manifest['core_min_mapq']),float(manifest['core_max_nm_rate']))==(100,20,.08)
    assert manifest['quality_preset']==preset
    return {'preset':selection['preset'],'configured':selection['quality_thresholds'],'effective_min_aligned_bp':100,
            'native_min_mapq':int(manifest['core_min_mapq']),'native_max_nm_rate':float(manifest['core_max_nm_rate']),
            'config_files_sha256':selection['config_files_sha256']}


def check_case(name, observed, expected):
    filters = observed['filter_counts']
    for key, value in expected.get('filter_counts', {}).items():
        assert filters.get(key) == value, (name, key, filters.get(key), value)
    for branch, want in expected.get('branches', {}).items():
        rows = observed['branches'][branch]
        if want.get('dispersion_tier') == 'None':
            assert filters['PRIMARY' if branch == 'PRIMARY' else 'STRICT'] == 0, (name, branch)
            assert not rows or all(r['dispersion_tier'].lower() == 'none' for r in rows), (name, branch, rows)
            continue
        assert len(rows) == 1, (name, branch, rows)
        for key, value in want.items():
            if isinstance(value, float):
                assert abs(float(rows[0][key]) - value) < 1e-6, (name, branch, key, rows[0][key], value)
            else:
                assert rows[0][key] == str(value), (name, branch, key, rows[0][key], value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, help='New directory, outside the source/release tree')
    args = parser.parse_args()
    out = Path(args.out).resolve()
    if out == ROOT or ROOT in out.parents:
        parser.error('--out must be outside the source/release tree')
    if out.exists():
        parser.error('--out must not exist; example never overwrites prior evidence')
    tools = {name: shutil.which(name) for name in ('minimap2', 'samtools', 'cutadapt', 'bash')}
    if not all(tools.values()) or sys.platform != 'linux':
        parser.error('Run in Linux/WSL with minimap2, samtools, cutadapt and bash in PATH')
    engine_files = [ROOT / 'mycowes.py', *sorted((ROOT/'config').glob('*.env')), ROOT / 'VERSION'] + sorted((ROOT / 'workflow').glob('*.py')) + sorted((ROOT / 'workflow').glob('*.sh'))
    engine_sha256 = {str(path.relative_to(ROOT)): digest(path) for path in engine_files}
    out.mkdir(parents=True)
    rng = random.Random(SEED)
    target, competitor, host, wrong_anchor = [random_sequence(rng) for _ in range(4)]
    expected = json.loads((ROOT / 'examples/expected_summary.json').read_text())
    host_fasta = out / 'toy_host.fa'
    write_fasta(host_fasta, [('toy_host_NOT_GRCh38', host)])
    host_index = out / 'toy_host.sr.mmi'
    commands = []
    command = [tools['minimap2'], '-x', 'sr', '-d', str(host_index), str(host_fasta)]
    commands.append(command)
    run_command(command, out / 'host_index.log')
    observations = {}
    effective_quality = {}
    equivalence = None
    manuscript_equivalence = None
    input_digests = {}
    for name, spec in expected['cases'].items():
        case_dir = out / name
        reference = case_dir / 'reference'
        joint = make_reference(reference, target, competitor, wrong_anchor if name == 'anchor_fail' else target)
        index = joint.with_name('ONE_REP_TARGET.AB_CDE_JOINT.sr.mmi')
        command = [tools['minimap2'], '-x', 'sr', '-d', str(index), str(joint)]
        commands.append(command)
        run_command(command, case_dir / 'joint_index.log')
        fragments = [(f'{name}_{i:02d}', pair_at(target, start))
                     for i, start in enumerate(spec.get('target_starts', []))]
        if name == 'normal':
            fragments += [('normal_risk', pair_at(target, 4000)),
                          ('normal_DF2', pair_at(target, 5000)),
                          ('normal_duplicate', pair_at(target, 100)),
                          ('normal_competitor', pair_at(competitor, 100)),
                          ('normal_host', pair_at(host, 100))]
        elif name == 'zero_candidate':
            fragments = [('zero_competitor', pair_at(competitor, 100))]
        r1, r2 = write_pairs(case_dir / 'inputs', fragments)
        for path in (r1, r2, joint, reference / 'baseline_v028/01_microbe/AB_TARGET.reference_qc_beta_v0.28.fa'):
            input_digests[str(path.relative_to(out))] = digest(path)
        command = [sys.executable, str(ROOT / 'mycowes.py'), 'run', '--sample', name,
                   '--r1', str(r1), '--r2', str(r2), '--out', str(case_dir / 'run'),
                   '--reference-root', str(reference), '--nominal-pe-length', '150']
        if 'bin_size' in spec:
            command += ['--bin-size', str(spec['bin_size'])]
        commands.append(command)
        env = dict(os.environ, MYCOWES_HOST_MMI=str(host_index), MYCOWES_HOST_FASTA=str(host_fasta))
        run_command(command, case_dir / 'run.log', env)
        observations[name] = collect_case(case_dir / 'run')
        check_case(name, observations[name], spec)
        effective_quality[name] = managed_quality(case_dir/'run')
        if name == 'normal':
            explicit = command.copy();explicit[explicit.index('--out')+1]=str(case_dir/'run_explicit_standard')
            explicit += ['--preset','standard','--min-aligned-bp','100','--min-mapq','20','--max-nm-rate','0.08']
            commands.append(explicit);run_command(explicit,case_dir/'run_explicit_standard.log',env)
            explicit_result=collect_case(case_dir/'run_explicit_standard');check_case(name,explicit_result,spec)
            assert comparable_case(observations[name])==comparable_case(explicit_result), 'Default/explicit100/20/.08 native results differ'
            equivalence={'status':'PASS','same_inputs':True,'default':comparable_case(observations[name]),
                         'explicit':comparable_case(explicit_result),'explicit_quality':managed_quality(case_dir/'run_explicit_standard'),
                         'comparison':'Exact filter counts and PRIMARY/STRICT branch records, excluding absolute file paths'}
            multi=reference/'baseline_v028/01_microbe/AB_CDE_JOINT.reference_qc_beta_v0.28.fa'
            multi_index=multi.with_name('AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi')
            build=[tools['minimap2'],'-x','sr','-d',str(multi_index),str(multi)]
            commands.append(build);run_command(build,case_dir/'manuscript_index.log')
            manuscript=command.copy();manuscript[manuscript.index('--out')+1]=str(case_dir/'run_manuscript')
            manuscript+=['--preset','manuscript','--risk-policy','legacy']
            commands.append(manuscript);run_command(manuscript,case_dir/'run_manuscript.log',env)
            observed=collect_case(case_dir/'run_manuscript');check_case(name,observed,spec)
            selected=json.loads((case_dir/'run_manuscript/reference_selection.json').read_text())
            assert selected['mode']=='MULTI_TARGET_V028'
            assert selected['mapq_policy']=='assembly-top-ties-gap10'
            assert selected['risk_policy']=='legacy'
            assert comparable_case(observations[name])==comparable_case(observed)
            manuscript_equivalence={'status':'PASS','quality':managed_quality(case_dir/'run_manuscript','manuscript'),
                                    'synthetic_risk_policy':'legacy','reference_fixture':'synthetic multi-profile',
                                    'mapq_policy':selected['mapq_policy']}
    versions = {}
    for name in ('minimap2', 'samtools', 'cutadapt'):
        text = subprocess.check_output([tools[name], '--version'], text=True).splitlines()
        versions[name] = text[0]
    assert engine_sha256 == {str(path.relative_to(ROOT)): digest(path) for path in engine_files}, 'Engine changed during validation; run once against a stable revision in a new output directory'
    summary = {'status': 'PASS', 'seed': SEED, 'synthetic_only': True, 'python': sys.version,
               'tools': tools, 'versions': versions, 'commands': commands, 'input_sha256': input_digests,
               'engine_sha256': engine_sha256,
               'expected_source_sha256': digest(ROOT / 'examples/expected_summary.json'),
               'observed': observations,
               'software_version':(ROOT/'VERSION').read_text().strip(),'manuscript_preset':manuscript_equivalence,
               'managed_quality':effective_quality,'default_standard_equivalence':equivalence,
               'scope': 'Native synthetic example using the packaged analysis engine.'}
    (out / 'native_validation_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('PASS: ' + str(out / 'native_validation_summary.json'))


if __name__ == '__main__':
    main()
