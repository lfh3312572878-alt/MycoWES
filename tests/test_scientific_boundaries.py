"""Boundary checks of actual production rules and strict FASTQ validation."""
import ast
import gzip
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'workflow' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dedup = load('exact_dedup', 'exact_pair_dedup_v0.2.py')
risk = load('risk_boundaries', 'annotate_fragment_risk_v0.1.py')
fragment = load('fragment_boundaries', 'parse_fragment_evidence_v0.1.py')

# The legacy breadth script executes its CLI at module scope. Load its pure
# production functions without running a mapper; native tests exercise wiring.
source = ROOT / 'workflow/evidence_anchor_breadth_v0.2.py'
tree = ast.parse(source.read_text())
names = {'dispersion_tier', 'anchor_status', 'fragment_bin_start', 'merged_intervals'}
rules = {}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)
                            and n.name in names], type_ignores=[]), str(source), 'exec'), rules)


class FastqValidation(unittest.TestCase):
    def run_case(self, r1, r2, compressed=False):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        inputs = [root / ('r1.fq.gz' if compressed else 'r1.fq'),
                  root / ('r2.fq.gz' if compressed else 'r2.fq')]
        for path, content in zip(inputs, (r1, r2)):
            if compressed:
                with gzip.open(path, 'wt') as handle:
                    handle.write(content)
            else:
                path.write_text(content)
        outputs = [root / name for name in ('out1.fq', 'out2.fq', 'groups.tsv')]
        result = subprocess.run([sys.executable, str(ROOT / 'workflow/exact_pair_dedup_v0.2.py'),
                                 *map(str, inputs), *map(str, outputs)],
                                capture_output=True, text=True)
        return result, outputs

    def test_valid_gzip_mates_and_ordered_key(self):
        result, outputs = self.run_case('@a/1\nAACG\n+\nIIII\n@b/1\nAACG\n+\nIIII\n',
                                        '@a/2\nTTCC\n+\nIIII\n@b/2\nTTCC\n+\nIIII\n', True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('a\t2\ta,b', outputs[2].read_text())
        self.assertEqual(outputs[0].read_text(), '@a\nAACG\n+\nIIII\n')

    def test_reverse_complement_pair_is_not_collapsed(self):
        result, outputs = self.run_case('@a\nAACG\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
                                        '@a\nTTCC\n+\nIIII\n@b\nCGTT\n+\nIIII\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(outputs[2].read_text().splitlines()), 3)

    def test_invalid_fastq_fails_before_output(self):
        valid = '@a\nACGT\n+\nIIII\n'
        cases = [
            (valid + '@b\nACGT\n+\nIIII\n', valid, 'pair count mismatch'),
            (valid, valid + '@b\nACGT\n+\nIIII\n', 'pair count mismatch'),
            ('@a\nACGT\n+\n', valid, 'incomplete FASTQ'),
            (valid, '@a\nACGT\n+\nIII\n', 'length mismatch'),
            (valid, '@other\nACGT\n+\nIIII\n', 'FASTQ mismatch'),
            ('a\nACGT\n+\nIIII\n', valid, 'header/separator'),
            ('@a\nACGT\n-\nIIII\n', valid, 'header/separator'),
            ('@a/2\nACGT\n+\nIIII\n', valid, 'wrong mate suffix'),
            (valid + valid, valid + valid, 'duplicate read name'),
        ]
        for r1, r2, expected in cases:
            with self.subTest(expected=expected):
                result, outputs = self.run_case(r1, r2)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected, result.stderr)
                self.assertTrue(all(not path.exists() for path in outputs))
                self.assertFalse(list(outputs[0].parent.glob('*.tmp')))

    def test_empty_pair_is_valid_zero_not_truncation(self):
        result, outputs = self.run_case('', '')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(outputs[0].read_text(), '')
        self.assertEqual(outputs[2].read_text().splitlines(),
                         ['representative_qname\texact_group_size\tmember_qnames'])

    def test_invalid_input_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'r1').write_text('@a\nACGT\n+\n')
            (root / 'r2').write_text('')
            outputs = [root / name for name in ('o1', 'o2', 'groups')]
            for path in outputs:
                path.write_text('previous-complete-result')
            with self.assertRaisesRegex(ValueError, 'incomplete FASTQ'):
                dedup.deduplicate(root/'r1', root/'r2', *outputs)
            self.assertTrue(all(path.read_text() == 'previous-complete-result' for path in outputs))

    def test_corrupt_gzip_fails_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            corrupt = root / 'r1.fq.gz'
            corrupt.write_bytes(b'not a gzip stream')
            other = root / 'r2.fq'
            other.write_text('@a\nACGT\n+\nIIII\n')
            outputs = [root / name for name in ('o1', 'o2', 'groups')]
            with self.assertRaises(OSError):
                dedup.deduplicate(corrupt, other, *outputs)
            self.assertTrue(all(not path.exists() for path in outputs))


class ScientificBoundaries(unittest.TestCase):
    def test_fmax_and_span_boundaries(self):
        grade = rules['dispersion_tier']
        self.assertEqual(grade(6, 3, 300, 0.50), 'strong')
        self.assertEqual(grade(6, 3, 300, 0.500001), 'supportive')
        self.assertEqual(grade(3, 2, 300, 2/3), 'supportive')
        self.assertEqual(grade(5, 3, 299, 0.40), 'supportive')
        self.assertEqual(grade(5, 3, 300, 0.80), 'weak')

    def test_anchor_failure_is_not_none_or_weak(self):
        self.assertEqual(rules['anchor_status'](0.90, 0.90), 'PASS')
        self.assertEqual(rules['anchor_status'](0.899999, 0.90), 'FAIL')
        self.assertEqual(rules['anchor_status'](0, 0.90), 'FAIL')

    def test_bin_size_changes_consumed_spatial_metrics(self):
        start = rules['fragment_bin_start']
        fragments = [(0, 100), (120, 220), (520, 620)]
        small = {start(a, b, 100) for a, b in fragments}
        default = {start(a, b, 250) for a, b in fragments}
        self.assertEqual(len(small), 3)
        self.assertEqual(len(default), 2)

    def test_invalid_breadth_parameters_fail_before_outputs(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'not-created'
            arguments = [sys.executable, str(source)]
            for flag in ('--bam', '--ab', '--r1', '--r2', '--groups', '--reps', '--minimap2', '--samtools'):
                arguments += [flag, 'unused']
            arguments += ['--out', str(output)]
            for option in (['--bin-size', '0'], ['--retention', 'nan']):
                result = subprocess.run(arguments + option, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())

    def test_covered_bp_is_physical_span_including_gap(self):
        self.assertEqual(sum(e-s for s,e in rules['merged_intervals']([(0, 300)])), 300)
        self.assertEqual(fragment.cigar_stats('100M', 1)[0] * 2, 200)

    def test_core_quality_boundaries_and_non_filters(self):
        row = {'pair_support_class': 'paired_target_consistent',
               'R1_primary_aligned_bases': '100', 'R2_primary_aligned_bases': '100',
               'R1_primary_MAPQ': '20', 'R2_primary_MAPQ': '20',
               'maximum_primary_NM_rate': '0.08', 'R1_primary_proper_pair': '0',
               'R2_primary_proper_pair': '0', 'R1_primary_soft_clipped_bases': '50',
               'R2_primary_soft_clipped_bases': '50'}
        self.assertEqual(risk.core_filter_reason(row, 100, 20, .08), 'PASS')
        for field, value in [('R1_primary_aligned_bases', '99'),
                             ('R2_primary_MAPQ', '19'), ('maximum_primary_NM_rate', '.080001')]:
            changed = dict(row, **{field: value})
            self.assertNotEqual(risk.core_filter_reason(changed, 100, 20, .08), 'PASS')

    def test_strict_score_advantage_not_delta_five(self):
        taxonomy = {
            'target': {'reference_class': 'TARGET', 'reference_species': 'target_species',
                       'species_taxid': '1'},
            'competitor': {'reference_class': 'COMPETITOR', 'reference_species': 'competitor_species',
                           'species_taxid': '2'}}
        def record(name, flag, score):
            return fragment.parse_record(f'q\t{flag}\t{name}\t1\t20\t100M50S\t*\t0\t0\t*\t*\tAS:i:{score}\tNM:i:0', taxonomy)
        rows = [record('target', 65, 100), record('competitor', 321, 99)]
        result = fragment.summarize_read_end(rows, 'R1')
        self.assertEqual(result['end_class'], 'target')
        self.assertEqual(result['delta_AS'], 1)
        rows[1]['AS'] = 100
        self.assertEqual(fragment.summarize_read_end(rows, 'R1')['end_class'], 'target_competitor_tie')


if __name__ == '__main__':
    unittest.main()
