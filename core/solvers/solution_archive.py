"""Shape-checked, unpickled displacement checkpoint for field recovery.

Saved scalar reports do not contain nodal DOFs; do not invent them from a tip
value. This archive contains the complete solved DOF vector, mesh and loading
metadata. Stored DOFs are replayed as given; this loader does not establish
that the originating system met equilibrium tolerances or that its FEM field
is mathematically converged. V2 locks in the complete orthotropic matrix.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from core.solvers.beam import build_beam

FORMAT_VERSION = "hermite_cubic_solved_dofs_v2"


def save_cubic_solution(path, beam, u, *, end_force, weight_density,
                        displacement_only_clamp=False, fem_tip_uy_in=None):
    """Save complete solved DOFs and physically relevant provenance."""
    if beam.__class__.__name__ != "CubicRotationHex8":
        # Class names can be refactored, so rely on underlying cubic terms too.
        if getattr(beam, 'recovery_order', None) != 3 and not getattr(beam, 'cubic_terms', None):
            raise ValueError("Only cubic recovered-rotation Hermite solutions are supported")
    vec = np.asarray(u, dtype=np.float64)
    if vec.shape != (len(beam.grid_points) * beam.ndof_per_node,) or not np.isfinite(vec).all():
        raise ValueError("Full finite solved DOF vector is required")
    metadata = {
        'format': FORMAT_VERSION,
        'formulation': 'cubic_rotation_hermite',
        'node_mesh': [int(beam.nx), int(beam.ny), int(beam.nz)],
        'beam_dimensions_in': [float(beam.Lx), float(beam.Ly), float(beam.Lz)],
        'longitudinal_E_psi': float(beam.orthotropic_constants.Ez),
        # A longitudinal E match alone does NOT establish stress-replay compatibility.
        # Record the full engineering-Voigt orthotropic law used for the solve.
        'constitutive_matrix_psi': np.asarray(beam.elastic_matrix(), dtype=float).tolist(),
        'end_force_lbf': float(end_force),
        'weight_density_lbf_in3': float(weight_density),
        'displacement_only_clamp': bool(displacement_only_clamp),
        'fem_tip_uy_in': None if fem_tip_uy_in is None else float(fem_tip_uy_in),
        'certified_3d_stress': False,
        'stress_data_saved': False,
    }
    if not np.isfinite([metadata['end_force_lbf'], metadata['weight_density_lbf_in3'],
                       metadata['longitudinal_E_psi']]).all() or metadata['weight_density_lbf_in3'] <= 0:
        raise ValueError('Invalid physical load or material')
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, U=vec, metadata=np.array(json.dumps(metadata)))
    return destination


def load_cubic_solution(path):
    """Rebuild matching FE discretization and load a complete DOF checkpoint."""
    with np.load(path, allow_pickle=False) as saved:
        if set(saved.files) != {'U', 'metadata'}:
            raise ValueError('Expected U and metadata only; not a valid solution state')
        u = np.array(saved['U'], dtype=np.float64, copy=True)
        metadata_raw = saved['metadata']
        if metadata_raw.shape != () or metadata_raw.dtype.kind not in ('U', 'S'):
            raise ValueError('Invalid metadata encoding')
        meta = json.loads(str(metadata_raw.item()))
    if meta.get('format') != FORMAT_VERSION or meta.get('formulation') != 'cubic_rotation_hermite':
        raise ValueError('Incompatible saved-solution format/formulation')
    mesh = meta.get('node_mesh')
    if not isinstance(mesh, list) or len(mesh) != 3 or any(type(n) is not int or n < 3 for n in mesh):
        raise ValueError('Invalid saved mesh')
    beam = build_beam(nx=mesh[0], ny=mesh[1], nz=mesh[2], formulation=meta['formulation'])
    expected_dimensions = np.array([beam.Lx, beam.Ly, beam.Lz])
    if not np.allclose(meta.get('beam_dimensions_in', []), expected_dimensions, atol=0, rtol=1.e-12):
        raise ValueError('Saved geometry does not match active source')
    if not np.isclose(meta['longitudinal_E_psi'], beam.orthotropic_constants.Ez, rtol=1.e-12):
        raise ValueError('Saved modulus differs from active material')
    # A future source change can modify Ex, Ey, Gxy, Gyz, Gxz, or Poisson
    # ratios without modifying Ez. A replay must fail rather than silently
    # reporting new stress values as if they came from the saved solve.
    stored_D = np.asarray(meta.get('constitutive_matrix_psi', []), dtype=float)
    current_D = np.asarray(beam.elastic_matrix(), dtype=float)
    if (stored_D.shape != (6, 6) or not np.all(np.isfinite(stored_D))
            or not np.allclose(stored_D, current_D, rtol=1.e-12, atol=1.e-9)):
        raise ValueError('Saved orthotropic constitutive matrix differs from active material')
    if u.shape != (len(beam.grid_points) * beam.ndof_per_node,) or not np.isfinite(u).all():
        raise ValueError('Saved state is incomplete or contains nonfinite DOFs')
    if not np.isfinite([meta['end_force_lbf'], meta['weight_density_lbf_in3']]).all() or meta['weight_density_lbf_in3'] <= 0:
        raise ValueError('Invalid saved loading')
    return beam, u, meta
