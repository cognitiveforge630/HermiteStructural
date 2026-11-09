"""Minimal reproducibility checks for the new convergence workflow."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from tools.convergence_study import _parse_levels, beam_displacement, main
from core.reference.awc import awc_reference_components, DEFAULT_WEIGHT_DENSITY
from core.reference.beam_theory import L, P


class ConvergenceStudyTests(unittest.TestCase):
    def test_levels_and_reference_tip(self):
        self.assertEqual(_parse_levels('2,2,4;3,3,5'), [(2,2,4),(3,3,5)])
        ref=awc_reference_components(DEFAULT_WEIGHT_DENSITY,end_force=P)
        self.assertAlmostEqual(beam_displacement(L,DEFAULT_WEIGHT_DENSITY,P),ref['total_uy_in'],places=12)

    def test_solve_report_and_view_coarse(self):
        with tempfile.TemporaryDirectory() as tmp:
            args=['--order','2','--levels','2,2,4','--out',tmp]
            main(['--stage','all']+args)
            self.assertTrue((Path(tmp)/'report.md').exists())
            self.assertTrue((Path(tmp)/'mesh_density.png').exists())
            self.assertTrue((Path(tmp)/'hex27_2x2x4_fields.png').exists())
            with np.load(Path(tmp)/'hex27_2x2x4.npz') as data:
                self.assertEqual(data['displacement'].shape[1],3)
                self.assertEqual(data['fe_stress'].shape, (3*2**3,8,6))
