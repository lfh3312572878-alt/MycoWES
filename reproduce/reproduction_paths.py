"""Common root arguments for preserved reference and analysis sources."""
import argparse
from pathlib import Path
import sys


def parse():
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('--project-root', required=True, type=Path,
                        help='Analysis workspace containing analysis/, software/ and envs/')
    parser.add_argument('--ref-root', required=True, type=Path,
                        help='Reference workspace containing the recorded resource subdirectories')
    parser.add_argument('--out-dir', required=True, type=Path, help='Output directory for this source entry')
    if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
        parser.print_help()
        print('Additional source-specific positional arguments follow these common root arguments.')
        raise SystemExit(0)
    args, remaining = parser.parse_known_args()
    args.project_root = args.project_root.expanduser().resolve()
    args.ref_root = args.ref_root.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()
    sys.argv[1:] = remaining
    return args


ROOTS = parse()
