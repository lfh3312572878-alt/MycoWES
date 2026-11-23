#!/usr/bin/env python3
"""MycoWES paired-end analysis with explicit resource identity and run state."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import hashlib
import math
import tempfile
from datetime import datetime, timezone

from workflow.reference_profiles import MODES, bindings, validate

ROOT = Path(__file__).resolve().parent
PRESETS = {'standard': (100, 20, 0.08), 'fast': (60, 5, 0.12)}
PRESET_ALIASES = {'manuscript': 'standard', 'relax2': 'fast'}
LEGACY_PRESETS = {'normal': 'standard', 'simple': 'fast', 'hard': 'standard'}
DIFFICULTIES = {name: PRESETS[preset] for name, preset in LEGACY_PRESETS.items()}
REFERENCE_MODES = {'simple': MODES[0], 'normal': MODES[1], 'hard': MODES[1]}


def resolve_run_settings(args):
    requested = PRESET_ALIASES.get(args.preset, args.preset)
    legacy = LEGACY_PRESETS.get(args.difficulty)
    if requested is not None and legacy is not None and requested != legacy:
        raise ValueError('--preset conflicts with the legacy --difficulty quality preset')
    preset = requested or legacy or 'standard'
    mode = MODES[1] if preset == 'standard' else MODES[0]
    if args.reference_mode is not None and args.reference_mode != mode:
        raise ValueError('Reference mode conflicts with preset; standard uses 67-target competition; fast uses 18 representatives (manuscript/relax2 are aliases)')
    values = tuple(default if override is None else override for default, override in
                   zip(PRESETS[preset], (args.min_aligned_bp, args.min_mapq, args.max_nm_rate)))
    length, mapq, nm = values
    if length <= 0 or not 0 <= mapq <= 255 or not math.isfinite(nm) or nm < 0:
        raise ValueError('Require aligned length >0, MAPQ in [0,255], and a finite nonnegative NM rate')
    return preset, mode, values


def file_identity(path):
    """Detect ordinary same-path replacement without repeatedly hashing large inputs."""
    p = Path(path).resolve()
    s = p.stat()
    with p.open('rb') as f:
        first = f.read(65536)
        f.seek(max(0, s.st_size - 65536))
        last = f.read(65536)
    # shortcut: edge digest + stat is not a full content hash; use full hashing for adversarial identity.
    return {'path': str(p), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns,
            'ctime_ns': s.st_ctime_ns, 'inode': s.st_ino,
            'edge_sha256': hashlib.sha256(first + last).hexdigest()}


def atomic_json(path, value):
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(value, f, indent=2)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def execute_run(output, selection, cmd, env):
    output.mkdir(parents=True, exist_ok=True)
    lock = output / '.run_lock'
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError('Active run or stale lock; inspect it manually, do not duplicate the run')
    state_path = output / 'run_state.json'
    started = False
    try:
        recorded = output / 'reference_selection.json'
        if recorded.exists():
            if json.loads(recorded.read_text()) != selection:
                raise ValueError('Output belongs to different input/reference identities or settings; use a new --out')
            state = json.loads(state_path.read_text()) if state_path.exists() else {}
            if state.get('status') != 'COMPLETE':
                raise ValueError('Previous run is incomplete or failed; preserved evidence is not reusable; use a new --out')
            for rel, identity in state.get('outputs', {}).items():
                if file_identity(output / rel) != identity:
                    raise ValueError('Completed output changed; refusing reuse: ' + rel)
            if not state.get('outputs'):
                raise ValueError('Completion has no output identity record')
            print('COMPLETE: verified existing run, no analysis repeated')
            return 0
        if any(p != lock for p in output.iterdir()):
            raise ValueError('Nonempty output has no run identity; use a new --out')
        atomic_json(recorded, selection)
        state = {'status': 'RUNNING', 'started_utc': datetime.now(timezone.utc).isoformat()}
        atomic_json(state_path, state)
        started = True
        attempt = Path(tempfile.mkdtemp(prefix='.attempt-', dir=output))
        cmd[cmd.index('--run-root') + 1] = str(attempt)
        result = subprocess.run(cmd, env=env)
        if result.returncode:
            state.update(status='FAILED', returncode=result.returncode, preserved_attempt=attempt.name)
            atomic_json(state_path, state)
            return result.returncode
        manifest = attempt / 'run_manifest.tsv'
        if not manifest.is_file() or 'status\tCOMPLETE' not in manifest.read_text():
            raise ValueError('Engine did not produce a complete native manifest')
        for mode in ('PRIMARY', 'STRICT_DF2'):
            summary = attempt / '07_breadth' / mode / '07_species_dedup_dispersion_summary.tsv'
            if not summary.is_file() or not summary.read_text().startswith('species\t'):
                raise ValueError('Missing terminal breadth summary: ' + mode)
        for item in attempt.iterdir():
            os.replace(item, output / item.name)
        attempt.rmdir()
        identities = {str(p.relative_to(output)): file_identity(p) for p in output.rglob('*')
                      if p.is_file() and p.name not in ('run_state.json', 'reference_selection.json')}
        state.update(status='COMPLETE', finished_utc=datetime.now(timezone.utc).isoformat(),
                     published_from=str(attempt), published_output_root=str(output), outputs=identities)
        atomic_json(state_path, state)
        return 0
    except BaseException as error:
        if started:
            state.update(status='FAILED', error=str(error))
            atomic_json(state_path, state)
        raise
    finally:
        lock.rmdir()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=(ROOT / 'VERSION').read_text().strip())
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('reference', 'run'):
        p = sub.add_parser(command)
        p.add_argument('--reference-mode', choices=MODES, default=MODES[1] if command == 'reference' else None)
        p.add_argument('--reference-root', default=os.environ.get('MYCOWES_REFERENCE_ROOT', str(ROOT / 'resources/reference')))
        if command == 'reference':
            p.add_argument('--validate', action='store_true')
        else:
            p.add_argument('--sample', required=True)
            p.add_argument('--r1', required=True)
            p.add_argument('--r2', required=True)
            p.add_argument('--out', required=True)
            p.add_argument('--nominal-pe-length', type=int, default=100)
            p.add_argument('--host-depleted', action='store_true', help='Input already trimmed and both-mates host-unmapped')
            p.add_argument('--bin-size', type=int, default=250)
            p.add_argument('--threads', type=int, help='Mapper and trimming threads (default 4)')
            p.add_argument('--samtools-threads', type=int, help='Additional samtools threads (default 1)')
            p.add_argument('--anchor-min-retention', type=float, default=0.90)
            p.add_argument('--preset', choices=(*PRESETS, *PRESET_ALIASES), help='Two modes: standard (default)=67-target competition/anchor, 100/20/0.08 and gap10/LC180; fast=18-target competition/67-target anchor, 60/5/0.12 and native MAPQ/legacy risk. Compatibility aliases: manuscript=standard, relax2=fast')
            p.add_argument('--min-aligned-bp', type=int, help='Explicit configured length; effective length is min(this, nominal PE length)')
            p.add_argument('--min-mapq', type=int, help='Explicit per-mate MAPQ threshold (0..255)')
            p.add_argument('--mapq-policy', choices=('native', 'assembly-ties-only', 'assembly-top-ties-gap10'), help='Default: assembly-top-ties-gap10 for standard, native for fast')
            p.add_argument('--max-nm-rate', type=float, help='Explicit finite, nonnegative normalized NM threshold')
            p.add_argument('--risk-policy', choices=('legacy', 'low-complexity-advisory', 'paired-high-cde-only', 'low-complexity-unique-pair', 'low-complexity-nonrisk100-pair', 'low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair'),
                           help='Default: low-complexity-nonrisk180-pair for standard, legacy for fast')
            p.add_argument('--difficulty', choices=tuple(LEGACY_PRESETS),
                           help='Legacy aliases: normal/hard=standard (67-target competition, gap10/LC180); simple=fast (60/5/0.12)')
    args = parser.parse_args(argv)
    try:
        mode = args.reference_mode
        if args.command == 'run':
            preset, mode, (length, mapq, nm) = resolve_run_settings(args)
            args.mapq_policy = args.mapq_policy or ('assembly-top-ties-gap10' if preset == 'standard' else 'native')
            args.risk_policy = args.risk_policy or ('low-complexity-nonrisk180-pair' if preset == 'standard' else 'legacy')
        values = bindings(args.reference_root, mode)
        if args.command == 'reference':
            if args.validate:
                validate(values)
            print(json.dumps({'mode': args.reference_mode, 'assets': values, 'validated': args.validate}, indent=2))
            return 0
        if preset == 'fast' and args.nominal_pe_length < 60:
            raise ValueError('FAST mode requires nominal PE length >=60 to preserve its 60 bp gate')
        validate(values)
        if args.nominal_pe_length <= 0:
            raise ValueError('--nominal-pe-length must be positive')
        if args.bin_size <= 0 or not 0 < args.anchor_min_retention <= 1:
            raise ValueError('Bin size must be positive and anchor retention in (0,1]')
        if not args.sample or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-' for c in args.sample) or args.sample in ('.', '..'):
            raise ValueError('--sample must be a safe filename identifier')
        for name in ('r1', 'r2'):
            path = Path(getattr(args, name))
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f'Missing or empty input: {path}')
        host = os.environ.get('MYCOWES_HOST_MMI', '')
        if not host or not Path(host).is_file() or Path(host).stat().st_size == 0:
            raise ValueError('Set MYCOWES_HOST_MMI to the existing GRCh38 host index')
        env = dict(os.environ, **values)
        for option, variable in ((args.threads, 'THREADS'), (args.samtools_threads, 'SAMTOOLS_THREADS')):
            if option is not None:
                if option < (0 if variable == 'SAMTOOLS_THREADS' else 1):
                    raise ValueError('Threads must be positive; additional samtools threads may be zero')
                env[variable] = str(option)
        env.update(CORE_MIN_ALIGNED_BP=str(length), CORE_MIN_MAPQ=str(mapq), CORE_MAX_NM_RATE=str(nm))
        env['CORE_MAPQ_POLICY'] = args.mapq_policy
        env.update(DISPERSION_BIN_BP=str(args.bin_size), ANCHOR_MIN_RETENTION=str(args.anchor_min_retention),
                   MYCOWES_MANAGED_RUN='1', MYCOWES_PRESET=preset, MYCOWES_RISK_POLICY=args.risk_policy)
        config = ROOT / 'config' / ('relax2.env' if preset == 'fast' else 'default.env')
        output = Path(args.out).resolve()
        selection = {'mode': mode, 'assets': values, 'host_index': str(Path(host).resolve()),
                     'config_sha256': hashlib.sha256(config.read_bytes()).hexdigest(),
                     'config_files_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'config').glob('*.env'))},
                     'inputs': [str(Path(args.r1).resolve()), str(Path(args.r2).resolve())],
                     'sample': args.sample, 'nominal_pe_length': args.nominal_pe_length}
        selection['difficulty'] = args.difficulty
        selection['preset'] = preset
        selection['software_version'] = (ROOT / 'VERSION').read_text().strip()
        selection['quality_thresholds'] = {'min_aligned_bp': length, 'min_mapq': mapq, 'max_nm_rate': nm}
        selection['effective_min_aligned_bp'] = min(length, args.nominal_pe_length)
        assets = [Path(args.r1), Path(args.r2), Path(host)] + [Path(v) for k,v in values.items()
                  if k.endswith(('_FASTA', '_MMI', '_TAXONOMY', '_IDENTITY', '_MANIFEST', '_BED', '_DF2'))]
        selection['local_file_identities'] = [file_identity(p) for p in assets]
        selection['execution'] = {k: env.get(k) for k in ('DISPERSION_BIN_BP', 'ANCHOR_MIN_RETENTION',
                                                        'THREADS', 'SAMTOOLS_THREADS', 'MINIMAP2', 'SAMTOOLS', 'CUTADAPT', 'PYTHON')}
        selection['host_depleted'] = args.host_depleted
        selection['risk_policy'] = args.risk_policy
        selection['mapq_policy'] = args.mapq_policy
        engine = [ROOT/'mycowes.py', *sorted((ROOT/'config').glob('*.env')), ROOT/'VERSION', *sorted((ROOT/'workflow').glob('*.py')),
                  *sorted((ROOT/'workflow').glob('*.sh'))]
        selection['engine_sha256'] = hashlib.sha256(b''.join(p.read_bytes() for p in engine)).hexdigest()
        cmd = ['bash', str(ROOT / 'workflow/run_one_sample.sh'), '--sample', args.sample,
               '--r1', str(Path(args.r1).resolve()), '--r2', str(Path(args.r2).resolve()),
               '--run-root', str(Path(args.out).resolve()), '--config', str(config),
               '--nominal-pe-length', str(args.nominal_pe_length)]
        if args.host_depleted:
            cmd.append('--host-depleted')
        return execute_run(output, selection, cmd, env)
    except (OSError, ValueError) as error:
        parser.exit(2, f'ERROR: {error}\n')


if __name__ == '__main__':
    sys.exit(main())
