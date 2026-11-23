"""Exact focal-artifact null and release preset boundaries."""
import ast
import math
from pathlib import Path
from types import SimpleNamespace
import unittest

import mycowes
from workflow.reference_profiles import MODES

ROOT = Path(__file__).resolve().parents[1]
source = ROOT/'workflow/evidence_anchor_breadth_v0.2.py'
functions = {'math': math}
exec(compile(ast.Module(body=[n for n in ast.parse(source.read_text()).body
                             if isinstance(n, ast.FunctionDef) and n.name == 'exact_spatial_dispersion_stats'],
                        type_ignores=[]), str(source), 'exec'), functions)
spatial = functions['exact_spatial_dispersion_stats']


class SpatialDispersion(unittest.TestCase):
    def test_focal_pileup(self):
        self.assertEqual(spatial(5, 5, 1, 800000)['p_focal_artifact_null'], 1.0)

    def test_minimum_strong_boundaries(self):
        for n, maximum in ((5, 2), (6, 3)):
            self.assertLess(spatial(n, maximum, 3, 800000)['p_focal_artifact_null'], 1e-4)

    def test_missing_evidence_is_na(self):
        self.assertEqual(set(spatial(0, 0, 0, 800000).values()), {'NA'})

    def test_selected_assembly_length_sets_denominator(self):
        self.assertEqual(spatial(5, 2, 3, 800001)['total_effective_bins'], 3201)

    def test_release_presets_keep_default_reference(self):
        base = dict(preset=None, difficulty=None, reference_mode=None,
                    min_aligned_bp=None, min_mapq=None, max_nm_rate=None)
        resolve = lambda **kw: mycowes.resolve_run_settings(SimpleNamespace(**dict(base, **kw)))
        self.assertEqual(resolve(), ('standard', MODES[0], (100, 20, .08)))
        self.assertEqual(resolve(preset='fast'), ('fast', MODES[0], (100, 20, .08)))
        self.assertEqual(resolve(preset='manuscript'), ('manuscript', MODES[1], (100, 20, .08)))


if __name__ == '__main__':
    unittest.main()
