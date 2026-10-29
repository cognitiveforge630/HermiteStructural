"""Recover 5x5x129 Hermite source-corner fields from a SAVED complete solution.

This command never solves FEM and never invents nodal displacements from the
scalar JSON research summaries. It requires --save-solution output from a
fresh solve. FE D@B@(T@U) and beam shear are evaluated at IDENTICAL points.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from core.solvers.solution_archive import load_cubic_solution
from core.visualization.corner_fields import prepare_corner_fields
from core.visualization.sampling import sample_volume
from core.visualization.shear_comparison import (
    FE_SHEAR, BEAM_SHEAR, DIFFERENCE, export_shear_vtu,
)


def _point_max(positions, values):
    idx = int(np.argmax(np.abs(values)))
    return {'absolute_psi': float(abs(values[idx])), 'signed_psi': float(values[idx]),
            'xyz_in': [float(x) for x in positions[idx]], 'sample_index': idx}


def summarize_stress(data, beam, meta, *, subdivisions, inset_fraction):
    """Compare matching samples, not section-resultant versus pointwise shear."""
    xyz = data.stress_evaluation_points
    fe = np.asarray(data.fields[FE_SHEAR])
    ref = np.asarray(data.fields[BEAM_SHEAR])
    delta = fe - ref
    span = beam.Lz
    area = beam.Lx * beam.Ly
    q = meta['weight_density_lbf_in3'] * area
    theoretical_clamp_max = 1.5 * abs(meta['end_force_lbf'] - q * span) / area
    if not np.all(np.isfinite(fe)) or not np.all(np.isfinite(ref)):
        raise ValueError('Nonfinite stress comparison')
    interior = ((xyz[:, 2] >= .25 * span) & (xyz[:, 2] <= .75 * span))
    if not np.any(interior):
        raise ValueError('No samples in the middle half of the beam')
    # Raw discrete L2 error, NOT an integral norm or certified error bound.
    rel_l2 = 100.0 * np.linalg.norm(delta[interior]) / max(np.linalg.norm(ref[interior]), 1.e-30)
    near = int(np.argmax(np.abs(ref)))
    raw_zz = np.asarray(data.fields['Sigma zz (psi)'])
    raw_zz_theory = np.asarray(data.fields['Beam reference sigma zz (psi)'])
    raw_vm = np.asarray(data.fields['Von Mises (psi)'])
    return {
        'status': 'RESEARCH_ONE_SIDED_CORNER_SAMPLES_NOT_CERTIFIED',
        'node_mesh': meta['node_mesh'], 'source': 'complete solved FEM DOFs',
        'stress_evaluation': 'source-element one-sided D @ B @ (T @ U)',
        'stress_sample_location': 'inset_corners',
        'stress_inset_fraction': float(inset_fraction),
        'visualization_subdivisions_per_source_axis': int(subdivisions),
        'corner_samples': int(len(fe)),
        'gravity_weight_density_lbf_in3': float(meta['weight_density_lbf_in3']),
        'end_force_lbf': float(meta['end_force_lbf']),
        'beam_maximum_shear_reference_psi': float(theoretical_clamp_max),
        'beam_maximum_shear_reference_location': 'z=0 in, y=0 in; classical beam centerline near clamp',
        'sampled_fem_tau_yz_maximum': _point_max(xyz, fe),
        'sampled_beam_tau_yz_maximum': _point_max(xyz, ref),
        'fem_tau_yz_at_beam_maximum_sample_psi': float(fe[near]),
        'difference_at_beam_maximum_sample_psi': float(delta[near]),
        'middle_half_span_raw_discrete_tau_yz_L2_difference_pct': float(rel_l2),
        'middle_half_span_max_pointwise_abs_tau_yz_difference_psi': float(np.max(np.abs(delta[interior]))),
        'sampled_fem_sigma_zz_maximum': _point_max(xyz, raw_zz),
        'sampled_beam_sigma_zz_maximum': _point_max(xyz, raw_zz_theory),
        'sampled_fem_von_mises_maximum': _point_max(xyz, raw_vm),
        'certified_3d_stress_error_bound': False,
        'within_one_percent_claim': False,
        'limitations': ('The 3.732 psi reference is the classical centerline shear maximum at the clamp; '
                        'sampled locations are inset and not at that exact point. '
                        'Do not equate integrated section-force discrepancy with local shear-stress discrepancy. '
                        '3D FE stress has clamp/end effects absent from beam theory. '
                        'The discrete L2 samples are not a volume-weighted error bound.'),
    }


def replay(archive, prefix, *, subdivisions=1, inset=0.95):
    if not 0 < inset < 1 or subdivisions < 1 or subdivisions > 6:
        raise ValueError('Choose positive subdivisions <=6 and 0<inset<1')
    beam, u, meta = load_cubic_solution(archive)
    sampled = sample_volume(beam, u, subdivisions=subdivisions,
                            stress_sample_location='inset_corners',
                            stress_inset_fraction=inset)
    data = prepare_corner_fields(sampled, beam,
                                 end_force=meta['end_force_lbf'],
                                 density=meta['weight_density_lbf_in3'])
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    vtu = prefix.with_name(prefix.name + '_one_sided_corner_fields.vtu')
    summary_path = prefix.with_name(prefix.name + '_stress_comparison.json')
    export_shear_vtu(data, vtu)
    report = summarize_stress(data, beam, meta, subdivisions=subdivisions,
                              inset_fraction=inset)
    report['vtu_file'] = str(vtu)
    report['complete_solution_state_file'] = str(archive)
    summary_path.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f'FE and beam tau_yz compared at {len(data.stress_evaluation_points):,} matching one-sided inset samples')
    print(f'Classical max centerline shear at clamp: {report["beam_maximum_shear_reference_psi"]:.6f} psi')
    print(f'FE max sampled |tau_yz|: {report["sampled_fem_tau_yz_maximum"]["absolute_psi"]:.6f} psi')
    print(f'Interior raw sample L2 discrepancy: {report["middle_half_span_raw_discrete_tau_yz_L2_difference_pct"]:.6f}%')
    print('NOT CERTIFIED: this is NOT a demonstration of 1% local-stress accuracy')
    print(f'Written: {vtu}\nWritten: {summary_path}')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('solution_state', type=Path, help='NPZ from tutorials.cubic --save-solution')
    parser.add_argument('--out-prefix', type=Path, default=Path('outputs/cubic_replayed'))
    parser.add_argument('--subdivisions', type=int, default=1,
                        help='Visualization subcells per original FE axis; 1 for manageable full-density export')
    parser.add_argument('--inset', type=float, default=0.95,
                        help='Stress point distance from subcell midpoint toward corners (default 95%%)')
    args = parser.parse_args(argv)
    return replay(args.solution_state, args.out_prefix, subdivisions=args.subdivisions,
                  inset=args.inset)


if __name__ == '__main__':
    main()
