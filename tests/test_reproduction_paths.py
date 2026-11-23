"""Common path arguments can be inspected without launching historical analyses."""
from pathlib import Path
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]


class ReproductionPaths(unittest.TestCase):
    def test_common_roots_parse_and_help(self):
        entry=ROOT/'reproduce/reproduction_paths.py'
        for flags in (['--help'],['--project-root','.','--ref-root','resources','--out-dir','output']):
            result=subprocess.run([sys.executable,str(entry),*flags],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
        result=subprocess.run([sys.executable,str(entry)],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)


if __name__=='__main__':
    unittest.main()
