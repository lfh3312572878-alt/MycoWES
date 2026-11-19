"""Opt-in risk exclusions. Annotation, quality gates and reference data stay unchanged."""
from collections import defaultdict
import hashlib
from pathlib import Path
import re

POLICIES = ('legacy', 'low-complexity-advisory', 'paired-high-cde-only', 'low-complexity-unique-pair',
            'low-complexity-nonrisk100-pair', 'low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair')
BED_SHA256 = 'aec94df33f252c76777300f45c0e5c3bbf121aed3f9192009ba95592cde8ff8e'
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def aligned_blocks(position, cigar):
    parts = CIGAR.findall(cigar)
    if position < 1 or ''.join(n + op for n, op in parts) != cigar:
        raise ValueError('Invalid primary SAM position/CIGAR')
    cursor, result = position - 1, []
    for length, operation in parts:
        length = int(length)
        if length < 1:
            raise ValueError('Zero-length CIGAR operation')
        if operation in 'M=X':
            result.append((cursor, cursor + length))
        if operation in 'MDN=X':
            cursor += length
    return result


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError('Invalid BED interval')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def _load_intervals(path, label):
    path = Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != BED_SHA256:
        raise ValueError('Risk policies require the recorded unchanged v0.28 risk BED')
    result = defaultdict(list)
    with path.open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 4 and fields[3] == label:
                result[fields[0]].append((int(fields[1]), int(fields[2])))
    if not result:
        raise ValueError('No recorded intervals for ' + label + '; refusing a disabled risk filter')
    return {key: merge(values) for key, values in result.items()}


def load_high_cde(path):
    return _load_intervals(path, 'PE151_CDE_highrisk')


def load_low_complexity(path):
    return _load_intervals(path, 'low_complexity_risk')


def decision(row, policy, high_cde=None, low_complexity=None):
    if policy not in POLICIES:
        raise ValueError('Unknown risk policy: ' + policy)
    flag = str(row.get('risk_any', '')).strip().lower()
    if flag not in ('0', '1', 'true', 'false', 'yes', 'no', 'y', 'n', 't', 'f', ''):
        raise ValueError('Invalid native risk_any flag')
    annotated = flag in ('1', 'true', 'yes', 'y', 't')
    core = str(row.get('exploratory_core_alignment', '')).strip().lower() in ('1', 'true')
    details = {'Risk_policy': policy, 'Legacy_risk_annotation_any': int(annotated),
               'Policy_excluded': 0, 'Policy_reason': 'NO_RISK_ANNOTATION',
               'R1_aligned_reference_bp': 'NA', 'R1_high_CDE_covered_bp': 'NA',
               'R2_aligned_reference_bp': 'NA', 'R2_high_CDE_covered_bp': 'NA',
               'R1_non_LC_aligned_bp': 'NA', 'R2_non_LC_aligned_bp': 'NA'}
    if policy == 'legacy':
        details.update(Policy_excluded=int(annotated), Policy_reason='LEGACY_ANY_RISK' if annotated else 'NO_RISK_ANNOTATION')
    elif policy in ('low-complexity-advisory', 'low-complexity-unique-pair', 'low-complexity-nonrisk100-pair', 'low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair'):
        categories = {value for value in row.get('risk_categories', '').split(';') if value}
        excluded = annotated and categories != {'low_complexity'}
        details.update(Policy_excluded=int(excluded), Policy_reason=(
            'NON_LOW_COMPLEXITY_RISK_RETAINED' if excluded else
            'LOW_COMPLEXITY_ADVISORY' if annotated else 'NO_RISK_ANNOTATION'))
        if policy != 'low-complexity-advisory' and annotated and not excluded:
            species = [{s for s in row.get(mate + '_best_target_species', '').split(';') if s}
                       for mate in ('R1', 'R2')]
            unique_pair = (all(row.get(mate + '_end_class') == 'target' for mate in ('R1', 'R2'))
                           and len(species[0]) == len(species[1]) == 1 and species[0] == species[1])
            details.update(Policy_excluded=int(not unique_pair), Policy_reason=(
                'LOW_COMPLEXITY_UNIQUE_CONCORDANT_PAIR_ADVISORY' if unique_pair else
                'LOW_COMPLEXITY_WITHOUT_UNIQUE_CONCORDANT_PAIR'))
            if policy in ('low-complexity-nonrisk100-pair', 'low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair') and unique_pair:
                if not core:
                    details['Policy_reason'] = 'QUALITY_GATE_NOT_PASSED'
                else:
                    if low_complexity is None:
                        raise ValueError('Missing fixed low-complexity interval dictionary')
                    supports = []
                    for mate in ('R1', 'R2'):
                        blocks = aligned_blocks(int(row[mate + '_primary_position']), row[mate + '_primary_CIGAR'])
                        total = sum(end - start for start, end in blocks)
                        covered = sum(max(0, min(end, b) - max(start, a)) for start, end in blocks
                                      for a, b in low_complexity.get(row[mate + '_primary_reference'], []))
                        supports.append(total - covered)
                        details[mate + '_aligned_reference_bp'] = total
                        details[mate + '_non_LC_aligned_bp'] = total - covered
                    if policy in ('low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair'):
                        minimum = 180 if policy == 'low-complexity-nonrisk180-pair' else 200
                        sufficient = all(value > 0 for value in supports) and sum(supports) >= minimum
                        reason = (f'LOW_COMPLEXITY_UNIQUE_PAIR_WITH_NONRISK{minimum}_ADVISORY' if sufficient else
                                  f'LOW_COMPLEXITY_INSUFFICIENT_NONRISK{minimum}_PAIR_SUPPORT')
                    else:
                        sufficient = all(value >= 100 for value in supports)
                        reason = ('LOW_COMPLEXITY_UNIQUE_PAIR_WITH_NONRISK100_ADVISORY' if sufficient else
                                  'LOW_COMPLEXITY_INSUFFICIENT_NONRISK100_SUPPORT')
                    details.update(Policy_excluded=int(not sufficient), Policy_reason=reason)
    elif annotated and not core:
        details['Policy_reason'] = 'QUALITY_GATE_NOT_PASSED'
    elif annotated:
        if high_cde is None:
            raise ValueError('Missing high-CDE interval dictionary')
        supported = []
        for mate in ('R1', 'R2'):
            blocks = aligned_blocks(int(row[mate + '_primary_position']), row[mate + '_primary_CIGAR'])
            total = sum(end - start for start, end in blocks)
            covered = sum(max(0, min(end, b) - max(start, a)) for start, end in blocks
                          for a, b in high_cde.get(row[mate + '_primary_reference'], []))
            details[mate + '_aligned_reference_bp'] = total
            details[mate + '_high_CDE_covered_bp'] = covered
            supported.append(row.get(mate + '_primary_class') == 'TARGET' and total >= 100 and covered == total)
        excluded = core and all(supported)
        details.update(Policy_excluded=int(excluded), Policy_reason=(
            'BOTH_MATES_FULLY_COVERED_HIGH_CDE' if excluded else 'ANNOTATION_ADVISORY_ONLY'))
    return details
