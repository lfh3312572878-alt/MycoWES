#!/usr/bin/env python3
"""Small local entry point for the unchanged single-sample analysis engine."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import hashlib

from workflow.reference_profiles import MODES, bindings, validate

ROOT = Path(__file__).resolve().parent
DIFFICULTIES = {'normal': (100, 20, 0.08), 'simple': (60, 5, 0.12), 'hard': (100, 20, 0.08)}
REFERENCE_MODES = {'simple': MODES[0], 'normal': MODES[0], 'hard': MODES[1]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=(ROOT / 'VERSION').read_text().strip())
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('reference', 'run'):
        p = sub.add_parser(command)
        p.add_argument('--reference-mode', choices=MODES, default=MODES[0] if command == 'reference' else None)
        p.add_argument('--reference-root', default=os.environ.get('MYCOWES_REFERENCE_ROOT', str(ROOT / 'resources/reference')))
        if command == 'reference':
            p.add_argument('--validate', action='store_true')
        else:
            p.add_argument('--sample', required=True)
            p.add_argument('--r1', required=True)
            p.add_argument('--r2', required=True)
            p.add_argument('--out', required=True)
            p.add_argument('--nominal-pe-length', type=int, default=100)
            p.add_argument('--difficulty', choices=('simple', 'normal', 'hard'), default='normal',
                           help='simple=RELAX2 single representative; normal=baseline single representative (default); hard=baseline multi-genome reference')
    args = parser.parse_args(argv)
    try:
        mode = args.reference_mode
        if args.command == 'run':
            selected = REFERENCE_MODES[args.difficulty]
            if mode is not None and mode != selected:
                raise ValueError('Reference mode conflicts with difficulty preset; use --difficulty hard for multi-reference')
            mode = selected
        values = bindings(args.reference_root, mode)
        if args.command == 'reference':
            if args.validate:
                validate(values)
            print(json.dumps({'mode': args.reference_mode, 'assets': values, 'validated': args.validate}, indent=2))
            return 0
        length, mapq, nm = DIFFICULTIES[args.difficulty]
        if args.difficulty == 'simple' and args.nominal_pe_length < 60:
            raise ValueError('RELAX2 simple profile requires nominal PE length >=60 to preserve its 60 bp gate')
        validate(values)
        if args.nominal_pe_length <= 0:
            raise ValueError('--nominal-pe-length must be positive')
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
        env.update(CORE_MIN_ALIGNED_BP=str(length), CORE_MIN_MAPQ=str(mapq), CORE_MAX_NM_RATE=str(nm))
        output = Path(args.out).resolve()
        selection = {'mode': mode, 'assets': values, 'host_index': str(Path(host).resolve()),
                     'config_sha256': hashlib.sha256((ROOT / 'config/default.env').read_bytes()).hexdigest(),
                     'inputs': [str(Path(args.r1).resolve()), str(Path(args.r2).resolve())],
                     'sample': args.sample, 'nominal_pe_length': args.nominal_pe_length}
        selection['difficulty'] = args.difficulty
        selection['software_version'] = (ROOT / 'VERSION').read_text().strip()
        selection['quality_thresholds'] = {'min_aligned_bp': length, 'min_mapq': mapq, 'max_nm_rate': nm}
        recorded = output / 'reference_selection.json'
        if recorded.exists():
            if json.loads(recorded.read_text()) != selection:
                raise ValueError('Output belongs to different inputs/config/reference; choose a new --out')
        elif output.exists() and any(output.iterdir()):
            raise ValueError('Nonempty output has no reference selection record; choose a new --out')
        else:
            output.mkdir(parents=True, exist_ok=True)
            recorded.write_text(json.dumps(selection, indent=2), encoding='utf-8')
        cmd = ['bash', str(ROOT / 'workflow/run_one_sample.sh'), '--sample', args.sample,
               '--r1', str(Path(args.r1).resolve()), '--r2', str(Path(args.r2).resolve()),
               '--run-root', str(Path(args.out).resolve()), '--config', str(ROOT / 'config/default.env'),
               '--nominal-pe-length', str(args.nominal_pe_length)]
        return subprocess.run(cmd, env=env).returncode
    except (OSError, ValueError) as error:
        parser.exit(2, f'ERROR: {error}\n')


if __name__ == '__main__':
    sys.exit(main())
