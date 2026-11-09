"""Source-sided, directly evaluated Hermite stress fields (NO centroid values).

VTK point coordinates follow the display mesh; stress sample coordinates may
be 95%-inset toward the same corners. Repeated VTK points prevent averaging
stress across source-element interfaces. Do not use a subcell's mean as a peak.
"""
from __future__ import annotations

import numpy as np

from .shear_comparison import prepare_shear_comparison, ShearComparison
from .cubic_stress import analytical_stress_at


def prepare_corner_fields(sample, beam, *, end_force, density) -> ShearComparison:
    """Provide point arrays for all six one-sided stresses and displacements."""
    data = prepare_shear_comparison(sample, beam, end_force=end_force, density=density)
    raw = np.asarray(sample.stress_corner_values, dtype=float).reshape(-1, 6)
    sx, sy, sz, txy, tyz, txz = raw.T
    vm = np.sqrt(0.5 * ((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2)
                 + 3.0 * (txy**2 + tyz**2 + txz**2))
    tensor = np.zeros((len(raw), 3, 3), dtype=np.float64)
    tensor[:, 0, 0], tensor[:, 1, 1], tensor[:, 2, 2] = sx, sy, sz
    tensor[:, 0, 1] = tensor[:, 1, 0] = txy
    tensor[:, 1, 2] = tensor[:, 2, 1] = tyz
    tensor[:, 0, 2] = tensor[:, 2, 0] = txz
    principal = np.linalg.eigvalsh(tensor)
    sigma_beam, _ = analytical_stress_at(data.stress_evaluation_points, beam,
                                          end_force=end_force, density=density)
    arrays = {
        'Sigma xx (psi)': sx, 'Sigma yy (psi)': sy, 'Sigma zz (psi)': sz,
        'Absolute longitudinal stress (psi)': np.abs(sz),
        'Longitudinal tension (psi)': np.maximum(sz, 0),
        'Longitudinal compression (psi)': np.maximum(-sz, 0),
        'Tau xy (psi)': txy, 'Tau yz (psi)': tyz, 'Tau xz (psi)': txz,
        'Von Mises (psi)': vm,
        'Maximum principal (psi)': principal[:, 2],
        'Minimum principal (psi)': principal[:, 0],
        'Beam reference sigma zz (psi)': sigma_beam,
        'Cubic minus beam sigma zz (psi)': sz - sigma_beam,
        'Stress evaluation X (in)': data.stress_evaluation_points[:, 0],
        'Stress evaluation Y (in)': data.stress_evaluation_points[:, 1],
        'Stress evaluation Z (in)': data.stress_evaluation_points[:, 2],
    }
    # Displacement is sampled at *visual geometry vertices*. Unlike the stress
    # arrays, it is not evaluated at inset stress sample coordinates.
    indices = np.asarray(sample.connectivity).ravel()
    vec = np.asarray(sample.displacement)[indices]
    arrays['Displacement (in)'] = vec
    arrays['Displacement magnitude (in)'] = np.linalg.norm(vec, axis=1)
    arrays['Vertical displacement magnitude (in)'] = np.abs(vec[:, 1])
    arrays['Ux (in)'], arrays['Uy (in)'], arrays['Uz (in)'] = vec.T
    for mark in ('Fixed support node', 'Loaded end node'):
        arrays[mark] = np.asarray(sample.point_fields[mark])[indices]
    data.fields.update(arrays)
    return data
