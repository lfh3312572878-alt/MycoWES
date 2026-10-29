"""Microbial resource bindings only. No host or filtering settings."""
from pathlib import Path

MODES = ('ONE_REP_TARGET_DEVELOPMENT', 'MULTI_TARGET_V028')


def bindings(root, mode=MODES[0]):
    if mode not in MODES:
        raise ValueError(f'Unknown reference mode: {mode}')
    root = Path(root).expanduser().resolve()
    baseline = root / 'baseline_v028'
    microbe = baseline / '01_microbe'
    metadata = baseline / '06_metadata'
    single = root / 'ONE_REP_TARGET_DEVELOPMENT'
    one = mode == MODES[0]
    return {
        'MYCOWES_REFERENCE_RELEASE': ('MycoWES_Ref_v0.28_ONE_REP_TARGET_DEVELOPMENT' if one else 'MycoWES_Ref_v0.28_REFERENCE_QC_BETA'),
        'MYCOWES_REFERENCE_PROFILE': ('REFERENCE_STRUCTURE_RESEARCH' if one else 'REFERENCE_QC_BETA'),
        'MYCOWES_MICROBIAL_PROFILE_ID': ('MycoWES_AB_CDE_v0.28_ONE_REP_TARGET_DEV_20261015' if one else 'MycoWES_AB_CDE_v0.28_reference_qc_beta'),
        'MYCOWES_JOINT_FASTA': str(single / 'ONE_REP_TARGET.AB_CDE_JOINT.fa' if one else microbe / 'AB_CDE_JOINT.reference_qc_beta_v0.28.fa'),
        'MYCOWES_JOINT_MMI': str(single / 'ONE_REP_TARGET.AB_CDE_JOINT.sr.mmi' if one else microbe / 'AB_CDE_JOINT.reference_qc_beta_v0.28.sr.mmi'),
        'MYCOWES_TAXONOMY': str(single / 'ONE_REP_TARGET.taxonomy.tsv' if one else metadata / 'AB_CDE_reference_taxonomy_lineage_map.reference_qc_beta_v0.28.tsv'),
        'MYCOWES_MICROBIAL_IDENTITY': str(single / 'ONE_REP_TARGET.identity.tsv' if one else metadata / 'identity.tsv'),
        'MYCOWES_REFERENCE_EXCLUSION_MANIFEST': str(single / 'ONE_REP_TARGET.exclusions.tsv' if one else metadata / 'reference_exclusion_manifest.tsv'),
        # The historical single-representative joint profile uses these shared assets.
        'MYCOWES_TARGET_FASTA': str(microbe / 'AB_TARGET.reference_qc_beta_v0.28.fa'),
        'MYCOWES_RISK_BED': str(baseline / '05_annotations/AB_reference_risk_annotations.reference_qc_beta_v0.28.bed'),
        'MYCOWES_GLOBAL_CAPTURE_DF2': str(baseline / '06_capture_bias/panWES.design_family_support_ge2.reference_qc_beta_v0.28.bed'),
        'MYCOWES_PROFILE_READY': '1',
    }


def validate(values):
    required = [key for key in values if key.endswith(('_FASTA', '_MMI', '_TAXONOMY', '_IDENTITY', '_MANIFEST', '_BED', '_DF2'))]
    for key in required:
        path = Path(values[key])
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Missing or empty {key}: {path}')
    identity = {}
    for line in Path(values['MYCOWES_MICROBIAL_IDENTITY']).read_text(encoding='utf-8-sig').splitlines():
        pair = line.split('\t', 1)
        if len(pair) == 2:
            identity[pair[0]] = pair[1]
    if identity.get('profile_id') != values['MYCOWES_MICROBIAL_PROFILE_ID']:
        raise ValueError('Microbial profile identity mismatch; refusing reference fallback')
