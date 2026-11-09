#!/usr/bin/env python3
"""Validate paired FASTQ once, then collapse ordered exact R1/R2 sequences."""

import gzip
import os
from pathlib import Path
import sys
import tempfile
from itertools import zip_longest


def fq(filename, mate):
    opener = gzip.open if str(filename).endswith('.gz') else open
    with opener(filename, 'rt', encoding='ascii', newline='') as handle:
        record = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record += 1
            lines = [handle.readline() for _ in range(3)]
            if not all(lines):
                raise ValueError(f'{filename}: incomplete FASTQ record {record}')
            sequence, plus, quality = [line.rstrip('\r\n') for line in lines]
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'{filename}: invalid FASTQ header/separator at record {record}')
            fields = header[1:].strip().split()
            if not fields:
                raise ValueError(f'{filename}: empty read name at record {record}')
            name = fields[0]
            if name.endswith(('/1', '/2')):
                if name[-1] != str(mate):
                    raise ValueError(f'{filename}: wrong mate suffix at record {record}')
                name = name[:-2]
            if not name:
                raise ValueError(f'{filename}: empty read name at record {record}')
            if not sequence or any(c.isspace() for c in sequence):
                raise ValueError(f'{filename}: empty/whitespace sequence at record {record}')
            if len(sequence) != len(quality):
                raise ValueError(f'{filename}: sequence/quality length mismatch at record {record}')
            if any(not 33 <= ord(c) <= 126 for c in quality):
                raise ValueError(f'{filename}: invalid quality character at record {record}')
            yield name, sequence, quality


def deduplicate(r1fn, r2fn, o1fn, o2fn, groupfn):
    inputs = {Path(r1fn).resolve(), Path(r2fn).resolve()}
    outputs = [Path(p).resolve() for p in (o1fn, o2fn, groupfn)]
    if len(set(outputs)) != 3 or inputs.intersection(outputs):
        raise ValueError('FASTQ output paths must be distinct and must not replace inputs')
    groups = {}
    seen = set()
    raw = 0
    for a, b in zip_longest(fq(r1fn, 1), fq(r2fn, 2)):
        if a is None or b is None:
            raise ValueError('FASTQ pair count mismatch: one mate ended before the other')
        q1, s1, qual1 = a
        q2, s2, qual2 = b
        if q1 != q2:
            raise ValueError(f'FASTQ mismatch: {q1} != {q2}')
        if q1 in seen:
            raise ValueError(f'FASTQ duplicate read name: {q1}')
        seen.add(q1)
        raw += 1
        key = (s1.upper(), s2.upper())
        if key not in groups:
            groups[key] = {'rep': q1, 'members': [], 'r1': (s1, qual1), 'r2': (s2, qual2)}
        groups[key]['members'].append(q1)

    # No final output is touched until both streams have passed validation.
    temporary = []
    handles = []
    try:
        for output in outputs:
            handle = tempfile.NamedTemporaryFile(mode='w', encoding='ascii', newline='\n',
                                                 dir=output.parent, prefix=output.name + '.',
                                                 suffix='.tmp', delete=False)
            handles.append(handle)
            temporary.append(Path(handle.name))
        o1, o2, table = handles
        table.write('representative_qname\texact_group_size\tmember_qnames\n')
        for item in groups.values():
            rep = item['rep']
            for handle, key in ((o1, 'r1'), (o2, 'r2')):
                sequence, quality = item[key]
                handle.write(f'@{rep}\n{sequence}\n+\n{quality}\n')
            table.write(f"{rep}\t{len(item['members'])}\t{','.join(item['members'])}\n")
        for handle in handles:
            handle.close()
        for source, output in zip(temporary, outputs):
            os.replace(source, output)
    finally:
        for handle in handles:
            handle.close()
        for path in temporary:
            path.unlink(missing_ok=True)
    print('raw_paired_fragments =', raw)
    print('exact_dedup_fragments =', len(groups))
    print('exact_duplicate_excess =', raw - len(groups))


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 5:
        print('Usage: exact_pair_dedup_v0.2.py R1 R2 OUT1 OUT2 GROUPS', file=sys.stderr)
        return 2
    try:
        deduplicate(*argv)
    except (OSError, ValueError, EOFError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
