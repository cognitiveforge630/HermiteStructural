"""Independent kinematic completeness audit for transverse cubic shear warping.

This is a prescribed-displacement *patch test*, not a solved shear benchmark.
The exact field u_z=y**3 (with physical spin omega_x=1.5*y**2) has the
engineering shear strain gamma_yz=3*y**2. Cubic recovered Hermite stencils
cannot reproduce its y-cubic derivative when only three y node layers exist.

Do not confuse passing this audit with satisfying natural shear tractions,
solved three-dimensional equilibrium or a certified a posteriori error bound.
"""
from __future__ import annotations

import numpy as np
from core.elements.cubic import CubicRotationHex8

# Include element corners and non-centered interior parent-coordinate samples.
_SAMPLES = ((0.0, 0.0, 0.0), (0.33, -0.42, 0.2),
            (0.85, 0.9, -0.4), (-0.88, -0.7, 0.5),
            (-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))


def audit_cubic_y_warping(beam, *, points=_SAMPLES):
    """Compare source-element strains with exact strains for u_z=y^3.

    Uses the original Hermite B and T, nodal translations and nodal spins.
    It does not use a solved cantilever, a textbook shear-stress field,
    calibration parameters, an altered material, or smoothed recovery stresses.
    """
    if not isinstance(beam, CubicRotationHex8):
        raise TypeError('Cubic shear-warping completeness requires CubicRotationHex8')
    points = tuple(tuple(float(t) for t in point) for point in points)
    if not points or any(len(p) != 3 or not np.all(np.isfinite(p))
                         or any(abs(v) > 1 for v in p) for p in points):
        raise ValueError('Sample locations must be nonempty parent points inside [-1,1]^3')

    U = np.zeros(beam.grid_points.shape[0] * beam.ndof_per_node)
    y = beam.grid_points[:, 1]
    U[2::beam.ndof_per_node] = y**3
    # 0.5 curl(u) = (1.5*y^2,0,0); important for the shared-spin B/T.
    U[3::beam.ndof_per_node] = 1.5 * y**2
    errors = []
    exact_norm = []
    max_abs = 0.
    for e in range(beam.get_num_elements()):
        coords = beam.get_element_nodes(e)
        ids, T = beam.element_transform(e)
        primitive = T @ U[ids]
        for parent in points:
            B, detJ, _ = beam.primitive_B(*parent, coords)
            if not np.isfinite(detJ) or detJ <= 0 or not np.isfinite(B).all():
                raise ArithmeticError('Invalid B or Jacobian during cubic warping audit')
            xyz = beam.hex8_shape_functions(*parent) @ coords
            reference = np.array([0., 0., 0., 0., 3.*xyz[1]**2, 0.])
            difference = B @ primitive - reference
            max_abs = max(max_abs, float(np.max(np.abs(difference))))
            errors.append(float(difference @ difference))
            exact_norm.append(float(reference @ reference))

    present = (0, 3, 0) in {tuple(p) for p in beam.cubic_terms}
    return {
        'status': 'KINEMATIC_PATCH_ONLY_NOT_SHEAR_STRESS_CERTIFICATION',
        'formulation': 'cubic_rotation_hermite',
        'node_mesh': [beam.nx, beam.ny, beam.nz],
        'polynomial': 'u_z=y^3, omega_x=1.5*y^2, gamma_yz=3*y^2',
        'actual_source_B_and_T': True,
        'y_cubic_recovery_term_present': bool(present),
        'sampled_elements': int(beam.get_num_elements()),
        'samples_per_element': len(points),
        'relative_discrete_strain_L2': float(np.sqrt(sum(errors)/max(sum(exact_norm), 1e-30))),
        'maximum_absolute_strain_difference': float(max_abs),
        'samples_include_source_corners': any(all(abs(v) == 1 for v in p) for p in points),
        'certified_error_bound': False,
        'limitations': ('Prescribed-field polynomial completeness only. Not a solved-load shear accuracy '
                        'comparison, three-dimensional traction admissibility or a posteriori certification.'),
    }
