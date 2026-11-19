"""Boundary tests of opt-in policy, without asserting biological performance."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('risk_policy', ROOT/'workflow/risk_filter_policy.py')
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


def pair(**changes):
    row = dict(risk_any='1', risk_categories='low_complexity', exploratory_core_alignment='1',
        R1_primary_reference='ctg', R2_primary_reference='ctg', R1_primary_position='101',
        R2_primary_position='401', R1_primary_CIGAR='100M', R2_primary_CIGAR='100M',
        R1_primary_class='TARGET', R2_primary_class='TARGET')
    return dict(row, **changes)


class RiskPolicy(unittest.TestCase):
    def test_nonrisk180_boundaries_and_other_risks(self):
        row = pair(R1_primary_CIGAR='150M', R2_primary_CIGAR='150M',
                   R1_best_target_species='a', R2_best_target_species='a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-nonrisk180-pair'
        for intervals, excluded in (([(100,160),(400,460)],0),
                                    ([(100,161),(400,460)],1),
                                    ([(100,250)],1)):
            self.assertEqual(policy.decision(row,name,low_complexity={'ctg':intervals})['Policy_excluded'],excluded)
        intervals={'ctg':[(100,160),(400,460)]}
        self.assertEqual(policy.decision(row,'low-complexity-nonrisk200-pair',low_complexity=intervals)['Policy_excluded'],1)
        for changes in ({'risk_categories':'low_complexity;human_homology'},
                        {'risk_categories':'PE151_CDE_highrisk'},
                        {'R2_best_target_species':'b'}, {'R1_best_target_species':'a;b'}):
            self.assertEqual(policy.decision(dict(row,**changes),name,low_complexity=intervals)['Policy_excluded'],1)
        with self.assertRaises(ValueError):policy.decision(row,name)
    def test_nonrisk200_pair_compensation_boundaries(self):
        row = pair(R1_primary_CIGAR='150M', R2_primary_CIGAR='150M',
                   R1_best_target_species='a', R2_best_target_species='a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-nonrisk200-pair'
        for intervals, excluded, supports in (
            ([(100, 150), (400, 450)], 0, (100, 100)),
            ([(100, 151), (400, 449)], 0, (99, 101)),
            ([(100, 151), (400, 450)], 1, (99, 100)),
        ):
            value = policy.decision(row, name, low_complexity={'ctg': intervals})
            self.assertEqual(value['Policy_excluded'], excluded)
            self.assertEqual((value['R1_non_LC_aligned_bp'], value['R2_non_LC_aligned_bp']), supports)
        compensated = {'ctg': [(100, 151), (400, 449)]}
        self.assertEqual(policy.decision(row, 'low-complexity-nonrisk100-pair', low_complexity=compensated)['Policy_excluded'], 1)
        value = policy.decision(dict(row, R2_primary_CIGAR='250M'), name,
                                low_complexity={'ctg': [(100, 250)]})
        self.assertEqual((value['Policy_excluded'], value['R1_non_LC_aligned_bp'], value['R2_non_LC_aligned_bp']), (1, 0, 250))
        for changes in ({'risk_categories': 'low_complexity;human_homology'},
                        {'R1_best_target_species': 'a;b'}, {'R2_best_target_species': 'b'}):
            self.assertEqual(policy.decision(dict(row, **changes), name)['Policy_excluded'], 1)
        self.assertEqual(policy.decision(dict(row, risk_any='0'), name)['Policy_excluded'], 0)
        with self.assertRaises(ValueError):
            policy.decision(row, name)

    def test_nonrisk100_inclusive_boundary_and_partial_mask(self):
        row = pair(R1_primary_CIGAR='150M', R2_primary_CIGAR='150M',
                   R1_best_target_species='a', R2_best_target_species='a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-nonrisk100-pair'
        value = policy.decision(row, name, low_complexity={'ctg': [(100, 150), (400, 450)]})
        self.assertEqual((value['Policy_excluded'], value['R1_non_LC_aligned_bp'], value['R2_non_LC_aligned_bp']), (0, 100, 100))
        value = policy.decision(row, name, low_complexity={'ctg': [(100, 151), (400, 450)]})
        self.assertEqual((value['Policy_excluded'], value['R1_non_LC_aligned_bp']), (1, 99))
        self.assertEqual(value['Policy_reason'], 'LOW_COMPLEXITY_INSUFFICIENT_NONRISK100_SUPPORT')

    def test_nonrisk_support_uses_measured_blocks_only(self):
        row = pair(R1_primary_CIGAR='5S100M10D50M3I', R2_primary_CIGAR='150M',
                   R1_best_target_species='a', R2_best_target_species='a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-nonrisk100-pair'
        value = policy.decision(row, name, low_complexity={'ctg': [(200, 210), (400, 450)]})
        self.assertEqual((value['Policy_excluded'], value['R1_aligned_reference_bp'], value['R1_non_LC_aligned_bp']), (0, 150, 150))
        value = policy.decision(row, name, low_complexity={'ctg': [(300, 350)]})
        self.assertEqual((value['Policy_excluded'], value['R1_non_LC_aligned_bp'], value['R2_non_LC_aligned_bp']), (0, 150, 150))

    def test_nonrisk_policy_retains_other_gates_and_fails_closed_on_missing_dictionary(self):
        row = pair(R1_best_target_species='a', R2_best_target_species='a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-nonrisk100-pair'
        with self.assertRaises(ValueError):
            policy.decision(row, name)
        self.assertEqual(policy.decision(dict(row, risk_any='0'), name)['Policy_excluded'], 0)
        self.assertEqual(policy.decision(dict(row, risk_categories='low_complexity;human_homology'), name)['Policy_excluded'], 1)
        self.assertEqual(policy.decision(dict(row, R2_best_target_species='a;b'), name)['Policy_excluded'], 1)
        value = policy.decision(dict(row, exploratory_core_alignment='0', R2_primary_position=''), name)
        self.assertEqual(value['Policy_reason'], 'QUALITY_GATE_NOT_PASSED')

    def test_unique_pair_low_complexity_gate(self):
        row = pair(R1_best_target_species='species_a', R2_best_target_species='species_a',
                   R1_end_class='target', R2_end_class='target')
        name = 'low-complexity-unique-pair'
        self.assertEqual(policy.decision(row, name)['Policy_excluded'], 0)
        for changes in ({'R2_best_target_species': 'species_a;species_b'},
                        {'R2_best_target_species': 'species_b'}, {'R2_best_target_species': ''},
                        {'R1_end_class': 'ambiguous'}, {'R2_end_class': 'competitor'}):
            changed = dict(row, **changes)
            self.assertEqual(policy.decision(changed, name)['Policy_excluded'], 1)
            self.assertEqual(policy.decision(dict(changed, risk_any='0'), name)['Policy_excluded'], 0)
        for categories in ('human_homology', 'low_complexity;PE151_CDE_highrisk', ''):
            self.assertEqual(policy.decision(dict(row, risk_categories=categories), name)['Policy_excluded'], 1)
        self.assertEqual(policy.decision(dict(row, truth_species='different', expected_taxid='999'), name),
                         policy.decision(row, name))

    def test_legacy_any_hit_and_low_complexity_warning(self):
        self.assertEqual(policy.decision(pair(), 'legacy')['Policy_excluded'], 1)
        self.assertEqual(policy.decision(pair(), 'low-complexity-advisory')['Policy_excluded'], 0)
        for categories in ('human_homology', 'rRNA', 'PE151_CDE_highrisk',
                           'low_complexity;human_homology', ''):
            self.assertEqual(policy.decision(pair(risk_categories=categories),
                                            'low-complexity-advisory')['Policy_excluded'], 1)

    def test_both_fully_high_versus_one_end_or_partial(self):
        high = {'ctg': [(100, 200), (400, 500)]}
        self.assertEqual(policy.decision(pair(), 'paired-high-cde-only', high)['Policy_excluded'], 1)
        for changed in ({'ctg': [(100, 200)]}, {'ctg': [(100, 199), (400, 500)]},
                        {'other_reference': [(0, 1000)]}, {}):
            self.assertEqual(policy.decision(pair(), 'paired-high-cde-only', changed)['Policy_excluded'], 0)

    def test_measured_blocks_no_insert_gap_or_deletions(self):
        self.assertEqual(policy.aligned_blocks(101, '5S40M10D60M2I3S'), [(100, 140), (150, 210)])
        self.assertEqual(policy.aligned_blocks(101, '50=10N50X'), [(100, 150), (160, 210)])
        row = pair(R1_primary_CIGAR='40M10D60M')
        result = policy.decision(row, 'paired-high-cde-only', {'ctg': [(140, 150), (400, 500)]})
        self.assertEqual(result['R1_high_CDE_covered_bp'], 0)
        self.assertEqual(result['Policy_excluded'], 0)
        result = policy.decision(pair(), 'paired-high-cde-only', {'ctg': [(200, 400)]})
        self.assertEqual(result['Policy_excluded'], 0)

    def test_half_open_union_and_short_mates(self):
        self.assertEqual(policy.merge([(100, 150), (120, 200), (200, 210)]), [(100, 210)])
        result = policy.decision(pair(R1_primary_CIGAR='99M'), 'paired-high-cde-only', {'ctg': [(0, 1000)]})
        self.assertEqual(result['Policy_excluded'], 0)
        result = policy.decision(pair(), 'paired-high-cde-only', {'ctg': policy.merge([(100, 150), (140, 200), (400, 500)])})
        self.assertEqual(result['R1_high_CDE_covered_bp'], 100)

    def test_quality_failure_and_invalid_inputs(self):
        result = policy.decision(pair(exploratory_core_alignment='0', R2_primary_position='',
                                      R2_primary_CIGAR='*'), 'paired-high-cde-only', {})
        self.assertEqual(result['Policy_reason'], 'QUALITY_GATE_NOT_PASSED')
        self.assertEqual(result['Policy_excluded'], 0)
        with self.assertRaises(ValueError):
            policy.decision(pair(risk_any='UNKNOWN'), 'legacy')
        with self.assertRaises(ValueError):
            policy.decision(pair(), 'disabled')
        with self.assertRaises(ValueError):
            policy.aligned_blocks(101, '100Mgarbage')


if __name__ == '__main__':
    unittest.main()
