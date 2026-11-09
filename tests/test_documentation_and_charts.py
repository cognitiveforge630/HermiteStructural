"""Regression gates for the traceability of printed plots and README commands."""
from __future__ import annotations

import tempfile
import json

import numpy as np
import unittest
from pathlib import Path

from tools.audit_documentation import check_docs
from tools.plot_cubic_densification import load_measured_cases, render as plot_cubic
from tools.plot_p_order_comparison import render as plot_p_order


class DocumentationAndChartTests(unittest.TestCase):
    def test_local_links_and_literal_documented_python_entrypoints(self):
        audit=check_docs()
        self.assertGreater(audit['checked_internal_links'],70)
        self.assertGreater(audit['checked_python_commands'],100)
        self.assertEqual(audit['failures'],[])

    def test_chart_reads_source_gauss_solved_records_not_retyped_values(self):
        axial, transverse=load_measured_cases()
        self.assertEqual([c['mesh'] for c in axial],
                         [[5,5,9],[5,5,17],[5,5,33],[5,5,65],[5,5,129]])
        self.assertEqual(transverse['mesh'],[7,7,17])
        self.assertTrue(all(c['positions']==[0.0] for c in [*axial,transverse]))
        self.assertAlmostEqual(axial[-1]['maximum_midspan_shear_vs_beam_pct'],.0029081495,places=8)
        self.assertAlmostEqual(transverse['maximum_midspan_shear_vs_beam_pct'],2.07985543,places=6)
        self.assertLess(transverse['tip_vs_beam_pct'], axial[-1]['tip_vs_beam_pct'])
        self.assertGreater(transverse['maximum_midspan_moment_vs_beam_pct'],
                           axial[-1]['maximum_midspan_moment_vs_beam_pct'])

    def test_shipped_nodal_smoke_fields_agree_with_saved_solved_json(self):
        root = Path(__file__).resolve().parents[1]
        for order in (27, 64):
            with self.subTest(element=f'Hex{order}'):
                base = root / 'samples' / f'hex{order}_3x3x5'
                record = json.loads(base.with_suffix('.json').read_text(encoding='utf-8'))
                with np.load(base.with_name(base.name+'_fields.npz')) as fields:
                    self.assertEqual(int(fields['order']), 2 if order==27 else 3)
                    self.assertEqual(fields['fe_stress'].shape,(16*(2 if order==27 else 3)**3,8,6))
                    self.assertEqual(str(fields['stress_sampling_scheme']),'inset_subcell_corners_v2')
                    self.assertTrue(np.isfinite(fields['displacement']).all())
                    self.assertTrue(np.isfinite(fields['fe_von_mises']).all())
                self.assertLess(record['equilibrium']['free_DOF_relative_residual_L2'],1e-7)
                self.assertFalse(record['certified_error_bound'])

    def test_plot_scripts_write_valid_svgs_without_solver(self):
        with tempfile.TemporaryDirectory() as tmp:
            for plot, name in ((plot_cubic,'cubic.svg'),(plot_p_order,'nodal.svg')):
                out=plot(Path(tmp)/name)
                self.assertTrue(out.is_file())
                self.assertIn('<svg',out.read_text(encoding='utf-8'))


if __name__=='__main__':
    unittest.main()
