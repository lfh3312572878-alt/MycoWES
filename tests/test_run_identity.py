"""CLI transaction unit fixtures, separate from the real native example."""
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
sys.path.insert(0, str(ROOT))
import mycowes


def complete_fixture(command, env):
    attempt = Path(command[command.index('--run-root')+1])
    (attempt/'run_manifest.tsv').write_text('status\tCOMPLETE\n')
    for mode in ('PRIMARY', 'STRICT_DF2'):
        summary = attempt/'07_breadth'/mode/'07_species_dedup_dispersion_summary.tsv'
        summary.parent.mkdir(parents=True)
        summary.write_text('species\tdispersion_tier\nunit_fixture\tweak\n')
    return subprocess.CompletedProcess(command, 0)


class RunIdentity(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.output = self.root/'out'
        self.selection = {'test_only': 'CLI transaction, not native evidence'}

    def command(self):
        return ['synthetic-unit-engine', '--run-root', str(self.output)]

    def state(self):
        return json.loads((self.output/'run_state.json').read_text())

    def complete(self):
        with patch.object(mycowes.subprocess, 'run', side_effect=complete_fixture) as run:
            self.assertEqual(mycowes.execute_run(self.output, self.selection, self.command(), {}), 0)
            run.assert_called_once()

    def test_complete_same_identity_reuses_without_subprocess(self):
        self.complete()
        self.assertEqual(self.state()['status'], 'COMPLETE')
        self.assertTrue(self.state()['outputs'])
        with patch.object(mycowes.subprocess, 'run') as run, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(mycowes.execute_run(self.output, self.selection, self.command(), {}), 0)
            run.assert_not_called()
        self.assertFalse((self.output/'.run_lock').exists())

    def test_failed_partial_is_preserved_and_never_reused(self):
        def fail(command, env):
            attempt = Path(command[command.index('--run-root')+1])
            (attempt/'partial.tsv').write_text('nonempty is not complete\n')
            return subprocess.CompletedProcess(command, 17)
        with patch.object(mycowes.subprocess, 'run', side_effect=fail):
            self.assertEqual(mycowes.execute_run(self.output, self.selection, self.command(), {}), 17)
        state = self.state()
        self.assertEqual(state['status'], 'FAILED')
        self.assertEqual(state['returncode'], 17)
        self.assertTrue((self.output/state['preserved_attempt']/'partial.tsv').is_file())
        self.assertFalse((self.output/'partial.tsv').exists())
        with patch.object(mycowes.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'incomplete or failed'):
                mycowes.execute_run(self.output, self.selection, self.command(), {})
            run.assert_not_called()
        self.assertEqual(self.state(), state)

    def test_same_path_replacement_changes_identity(self):
        source = self.root/'input.fastq'
        source.write_bytes(b'original')
        before = mycowes.file_identity(source)
        replacement = self.root/'replacement'
        replacement.write_bytes(b'replaced')
        os.utime(replacement, ns=(source.stat().st_atime_ns, before['mtime_ns']))
        os.replace(replacement, source)
        after = mycowes.file_identity(source)
        self.assertEqual(before['path'], after['path'])
        self.assertEqual(before['bytes'], after['bytes'])
        self.assertNotEqual(before['edge_sha256'], after['edge_sha256'])
        self.assertNotEqual(before, after)

    def test_output_tampering_rejects_reuse(self):
        self.complete()
        summary = self.output/'07_breadth/PRIMARY/07_species_dedup_dispersion_summary.tsv'
        summary.write_text(summary.read_text()+'changed\n')
        with patch.object(mycowes.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'Completed output changed'):
                mycowes.execute_run(self.output, self.selection, self.command(), {})
            run.assert_not_called()
        self.assertFalse((self.output/'.run_lock').exists())

    def test_existing_lock_is_not_deleted(self):
        lock = self.output/'.run_lock'
        lock.mkdir(parents=True)
        owner = lock/'owner.txt'
        owner.write_text('another run owns this lock')
        with patch.object(mycowes.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'Active run or stale lock'):
                mycowes.execute_run(self.output, self.selection, self.command(), {})
            run.assert_not_called()
        self.assertEqual(owner.read_text(), 'another run owns this lock')
        self.assertFalse((self.output/'run_state.json').exists())

    def test_exit_zero_without_manifest_is_failed_not_complete(self):
        with patch.object(mycowes.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0)):
            with self.assertRaisesRegex(ValueError, 'complete native manifest'):
                mycowes.execute_run(self.output, self.selection, self.command(), {})
        self.assertEqual(self.state()['status'], 'FAILED')
        self.assertNotIn('outputs', self.state())
        self.assertFalse((self.output/'.run_lock').exists())
        self.assertTrue(list(self.output.glob('.attempt-*')))


if __name__ == '__main__':
    unittest.main()
