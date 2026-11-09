"""Section resultants from *source* Hermite stress, not visualization subcells.

Research diagnostic. Raw constitutive stresses need not satisfy strong
three-dimensional equilibrium and no stress equilibration/certification is
performed here.  Cross sections are inside elements rather than on interfaces.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from core.reference.beam_theory import P


def integrate_sections(beam, U, *, order=5, axial_positions=(0.0,),
                       end_force=P, density=DEFAULT_WEIGHT_DENSITY):
    """Integrate F_y=∫σ_yz dA and M_x=∫yσ_zz dA at each axial cell slice.

    ``axial_positions`` are parent ζ values strictly inside (-1,1), to avoid
    ambiguous one-sided stress traces on a mesh interface. The FE stress is
    evaluated directly as ``D @ B @ U`` at tensor Gauss points. Each result
    also gives an independent sectional axial force and transverse shear x.

    The scalar moment reference assumes the beam's y origin is its centroid.
    """
    if isinstance(order, bool) or int(order) != order or order < 2:
        raise ValueError('Quadrature order must be an integer >= 2')
    order = int(order)
    zs = tuple(float(t) for t in axial_positions)
    if not zs or not all(np.isfinite(t) and -1. < t < 1. for t in zs):
        raise ValueError('All axial_positions must lie strictly inside (-1,1)')
    u = np.asarray(U, dtype=float)
    if u.shape != (beam.grid_points.shape[0]*beam.ndof_per_node,) or not np.isfinite(u).all():
        raise ValueError('U must contain finite solver DOFs')
    if not np.isfinite(end_force) or not np.isfinite(density):
        raise ValueError('Loads must be finite')
    gp, weights = leggauss(order)
    # Accumulate per axial layer, making each cross section independent of
    # the arbitrary ordering of the source Hex8 cells.
    result = []
    D = beam.elastic_matrix()
    area = beam.Lx*beam.Ly
    q = density*area
    for iz in range(beam.nz-1):
        for zeta in zs:
            accumulated = np.zeros(4, dtype=float)  # Fx, Fy, Fz, Mx
            area_integrated = 0.
            for ix in range(beam.nx-1):
                for iy in range(beam.ny-1):
                    e = (ix*(beam.ny-1)+iy)*(beam.nz-1)+iz
                    coords = beam.get_element_nodes(e)
                    ids, T = beam.element_transform(e)
                    primitive = T @ u[ids]
                    for ia, xi in enumerate(gp):
                        for ib, eta in enumerate(gp):
                            parent = (float(xi),float(eta),float(zeta))
                            B, _, _ = beam.primitive_B(*parent, coords)
                            stress = D @ (B @ primitive)
                            J = beam.get_hermite_jacobian(*parent, coords)
                            da = np.linalg.norm(np.cross(J[0],J[1])) * weights[ia]*weights[ib]
                            pos = beam.hex8_shape_functions(*parent) @ coords
                            accumulated += da * np.array([stress[5],stress[4],stress[2],
                                                          pos[1]*stress[2]])
                            area_integrated += da
            z = (iz + (zeta+1.)/2.) * beam.Lz / (beam.nz-1)
            expected_shear = end_force-q*(beam.Lz-z)
            expected_moment = -end_force*(beam.Lz-z)+.5*q*(beam.Lz-z)**2
            result.append({
                'z_in':float(z), 'axial_parent_zeta':zeta,
                'section_area_in2':float(area_integrated),
                'FE_shear_x_lbf':float(accumulated[0]),
                'FE_shear_y_lbf':float(accumulated[1]),
                'FE_axial_force_lbf':float(accumulated[2]),
                'FE_bending_moment_x_lbf_in':float(accumulated[3]),
                'beam_shear_y_lbf':float(expected_shear),
                'beam_bending_moment_x_lbf_in':float(expected_moment),
                'shear_y_reference_difference_pct':float(100*abs(accumulated[1]-expected_shear)/max(abs(expected_shear),1e-12)),
                'moment_x_reference_difference_pct':float(100*abs(accumulated[3]-expected_moment)/max(abs(expected_moment),1e-12)),
            })
    return {
        'status':'RESEARCH_SECTIONAL_DIAGNOSTIC_NOT_A_POSTERIORI',
        'stress_evaluated_from':'D @ source_element_B @ source_element_U',
        'visualization_subdivisions_affect_results':False,
        'cross_section_quadrature_order':order,
        'positions_parent_zeta':list(zs),
        'sections':result,
        'max_midspan_shear_y_difference_pct':float(max(
            s['shear_y_reference_difference_pct'] for s in result
            if beam.Lz/4<=s['z_in']<=3*beam.Lz/4
        )),
        'max_midspan_moment_x_difference_pct':float(max(
            s['moment_x_reference_difference_pct'] for s in result
            if beam.Lz/4<=s['z_in']<=3*beam.Lz/4
        )),
        'certified_error_bound':False,
        'limitations':('Constitutive FE stress is not statically admissible by construction; '
                       'a beam-theory cross-section check is not a 3D error estimator.'),
    }
