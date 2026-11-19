"""CLI/identity fixtures only: the engine is mocked, not biologically validated."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
import mycowes


def complete_fixture(command, env):
    """Satisfy only the real CLI transaction's completion contract."""
    attempt = Path(command[command.index('--run-root') + 1])
    (attempt / 'run_manifest.tsv').write_text('status\tCOMPLETE\n')
    for branch in ('PRIMARY', 'STRICT_DF2'):
        summary = attempt / '07_breadth' / branch / '07_species_dedup_dispersion_summary.tsv'
        summary.parent.mkdir(parents=True)
        summary.write_text('species\tdispersion_tier\nunit_fixture\tweak\n')
    return subprocess.CompletedProcess(command, 0)


class RiskCliWiring(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.assets = mycowes.bindings(self.root)
        for key, value in self.assets.items():
            if key.endswith(('_FASTA', '_MMI', '_TAXONOMY', '_IDENTITY', '_MANIFEST', '_BED', '_DF2')):
                path = Path(value)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('unit-test-only; not a native reference\n')
        Path(self.assets['MYCOWES_MICROBIAL_IDENTITY']).write_text(
            'profile_id\t' + self.assets['MYCOWES_MICROBIAL_PROFILE_ID'] + '\n')
        self.host = self.root / 'host.mmi'
        self.host.write_text('unit-test-only\n')
        for mate in ('R1', 'R2'):
            (self.root / (mate + '.fq')).write_text('@unit\nACGT\n+\nIIII\n')

    def arguments(self, output, flags=()):
        return ['run', '--sample', 'risk_unit', '--r1', str(self.root / 'R1.fq'),
                '--r2', str(self.root / 'R2.fq'), '--out', str(output),
                '--reference-root', str(self.root), '--nominal-pe-length', '150',
                '--host-depleted', *flags]

    def test_default_and_explicit_policies_reach_environment_and_identity(self):
        cases = [((), 'legacy'), (('--risk-policy', 'legacy'), 'legacy'),
                 (('--risk-policy', 'low-complexity-advisory'), 'low-complexity-advisory'),
                 (('--risk-policy', 'paired-high-cde-only'), 'paired-high-cde-only'),
                 (('--risk-policy', 'low-complexity-unique-pair'), 'low-complexity-unique-pair'),
                 (('--risk-policy', 'low-complexity-nonrisk100-pair'), 'low-complexity-nonrisk100-pair'),
                 (('--risk-policy', 'low-complexity-nonrisk200-pair'), 'low-complexity-nonrisk200-pair'),
                 (('--risk-policy', 'low-complexity-nonrisk180-pair'), 'low-complexity-nonrisk180-pair')]
        for index, (flags, policy) in enumerate(cases):
            with self.subTest(policy=policy, flags=flags):
                output = self.root / ('out' + str(index))
                with patch.dict(os.environ, {'MYCOWES_HOST_MMI': str(self.host),
                                             'MYCOWES_RISK_POLICY': 'stale-environment-value'}), \
                     patch.object(mycowes.subprocess, 'run', side_effect=complete_fixture) as engine:
                    self.assertEqual(mycowes.main(self.arguments(output, flags)), 0)
                engine.assert_called_once()
                env = engine.call_args.kwargs['env']
                selection = json.loads((output / 'reference_selection.json').read_text())
                self.assertEqual(env['MYCOWES_RISK_POLICY'], policy)
                self.assertEqual(selection['risk_policy'], policy)
                self.assertEqual(selection['quality_thresholds'],
                                 dict(min_aligned_bp=100, min_mapq=20, max_nm_rate=0.08))
                self.assertEqual(tuple(env[k] for k in ('CORE_MIN_ALIGNED_BP', 'CORE_MIN_MAPQ', 'CORE_MAX_NM_RATE')),
                                 ('100', '20', '0.08'))
                self.assertEqual(selection['effective_min_aligned_bp'], 100)
                self.assertEqual(selection['assets'], self.assets)
                self.assertEqual(json.loads((output / 'run_state.json').read_text())['status'], 'COMPLETE')

    def test_same_policy_reuses_but_changed_policy_cannot_reuse_output(self):
        output = self.root / 'one_identity'
        with patch.dict(os.environ, {'MYCOWES_HOST_MMI': str(self.host)}), \
             patch.object(mycowes.subprocess, 'run', side_effect=complete_fixture) as engine:
            self.assertEqual(mycowes.main(self.arguments(output)), 0)
            engine.assert_called_once()
            selection = (output / 'reference_selection.json').read_bytes()
            state = (output / 'run_state.json').read_bytes()
            engine.reset_mock()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(mycowes.main(self.arguments(output, ('--risk-policy', 'legacy'))), 0)
            engine.assert_not_called()
            for policy in ('low-complexity-advisory', 'paired-high-cde-only', 'low-complexity-unique-pair', 'low-complexity-nonrisk100-pair', 'low-complexity-nonrisk200-pair', 'low-complexity-nonrisk180-pair'):
                with self.subTest(policy=policy), contextlib.redirect_stderr(io.StringIO()), \
                     self.assertRaises(SystemExit) as failure:
                    mycowes.main(self.arguments(output, ('--risk-policy', policy)))
                self.assertEqual(failure.exception.code, 2)
                engine.assert_not_called()
                self.assertEqual((output / 'reference_selection.json').read_bytes(), selection)
                self.assertEqual((output / 'run_state.json').read_bytes(), state)

    def test_unknown_policy_is_rejected_before_output_or_engine(self):
        output = self.root / 'must_not_exist'
        with patch.object(mycowes.subprocess, 'run') as engine, \
             contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as failure:
            mycowes.main(self.arguments(output, ('--risk-policy', 'disabled')))
        self.assertEqual(failure.exception.code, 2)
        engine.assert_not_called()
        self.assertFalse(output.exists())

    def test_combined_policies_reach_identity_and_environment(self):
        output=self.root/'combined'
        flags=('--mapq-policy','assembly-top-ties-gap10','--risk-policy','low-complexity-nonrisk180-pair')
        with patch.dict(os.environ,{'MYCOWES_HOST_MMI':str(self.host)}), \
             patch.object(mycowes.subprocess,'run',side_effect=complete_fixture) as engine:
            self.assertEqual(mycowes.main(self.arguments(output,flags)),0)
        self.assertEqual(engine.call_args.kwargs['env']['CORE_MAPQ_POLICY'],'assembly-top-ties-gap10')
        selection=json.loads((output/'reference_selection.json').read_text())
        self.assertEqual(selection['mapq_policy'],'assembly-top-ties-gap10')
        self.assertEqual(selection['risk_policy'],'low-complexity-nonrisk180-pair')


if __name__ == '__main__':
    unittest.main()
