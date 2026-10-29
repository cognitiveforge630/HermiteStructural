"""Hex27/Hex64 visualization archive with *inset-corner* stress samples.

NON-NEGOTIABLE: Do not evaluate or report centroid-based stress. Each small
visualization subcell is sampled at eight distinct points near its corners,
strictly inside its owning source element (and its own subcell). Discontinuous
one-sided stresses are intentionally not averaged between source elements.
Neither these sampled extrema nor the classical beam comparison are certified.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from core.elements.nodal_hex import strain_matrix, tensor_shapes
from core.mesh.nodal_hex import StructuredNodalHexMesh, local_nodes
from core.reference.beam_theory import B, H, L

STRESS_INSET_FRACTION = 0.95
STRESS_SAMPLING_SCHEME = 'inset_subcell_corners_v2'


def von_mises(s):
    """Vectorized von Mises for Voigt stress order xx, yy, zz, xy, yz, xz."""
    s = np.asarray(s)
    xx, yy, zz, xy, yz, xz = np.moveaxis(s, -1, 0)
    return np.sqrt(.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)
                   +3*(xy**2+yz**2+xz**2))


@lru_cache(maxsize=8)
def inset_subcell_parent_points(order, inset=STRESS_INSET_FRACTION):
    """Return (p**3, 8, 3) one-sided inset corner coordinates.

    p cells split every source element axis into equal segments. Each corner
    is at +/- inset of the *subcell half-width*, not at a cell center.
    """
    if order not in (2, 3) or not 0 < inset < 1:
        raise ValueError('Supported order is 2 or 3 and inset must be in (0,1)')
    sample = []
    for ijk in product(range(order), repeat=3):
        cell = []
        for signs in product((-1, 1), repeat=3):
            cell.append([-1.0 + (2*i + 1 + inset*sign)/order
                         for i, sign in zip(ijk, signs)])
        sample.append(cell)
    return np.asarray(sample, dtype=float)


@lru_cache(maxsize=8)
def _strain_operators(order, extent, inset=STRESS_INSET_FRACTION):
    """Precompute genuine FE strain operators at the specified corner points."""
    parent = inset_subcell_parent_points(order, inset)
    dx = np.asarray(extent, dtype=float)
    operators = []
    for local in parent.reshape(-1, 3):
        _, derivatives = tensor_shapes(order, local)
        operators.append(strain_matrix(derivatives / (dx/2)[:, None]))
    return np.asarray(operators).reshape(order**3, 8, 6, 3*(order+1)**3)


def export_fields(path, mesh: StructuredNodalHexMesh, displacements,
                  constitutive, density, end_force):
    """Export displacement and *direct inset-corner* stress, no centroids.

    Saved stress shape: (number_of_visualization_subcells, 8, 6).
    Both the FE and beam-reference tensors are evaluated at exactly the same
    eight inset-corner locations. No stress averaging or extrapolation is
    used to produce extrema. The PyVista coloring may select one of these
    eight genuine values per subcell; no interpolated nodal stresses are made.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    coords = mesh.coordinates
    div = mesh.divisions
    dx = mesh.element_extent
    parent = inset_subcell_parent_points(mesh.order)
    bmat = _strain_operators(mesh.order, tuple(float(x) for x in dx))
    base = np.asarray(mesh.origin)
    nsub = int(np.prod(div))*mesh.order**3
    physical = np.empty((nsub, 8, 3))
    stresses = np.empty((nsub, 8, 6))
    reference = np.empty_like(stresses)
    q = density*B*H
    I = B*H**3/12
    offset = 0
    for ix, iy, iz in product(*(range(n) for n in div)):
        _, dofs = local_nodes(mesh.order, ix, iy, iz, mesh.node_shape)
        n = mesh.order**3
        uelement = displacements[dofs]
        eps = np.einsum('sacj,j->sac', bmat, uelement, optimize=True)
        stresses[offset:offset+n] = eps @ constitutive.T
        pts = base + (np.asarray((ix, iy, iz))
                      + (parent+1.)/2.)*dx
        physical[offset:offset+n] = pts
        y = pts[:, :, 1]
        z = pts[:, :, 2]
        r = L-z
        moment = -end_force*r + q*r*r/2
        sigzz = moment*y/I
        shear = (end_force-q*r)/(B*H)*1.5*(1-(2*y/H)**2)
        reference[offset:offset+n] = 0
        reference[offset:offset+n, :, 2] = sigzz
        reference[offset:offset+n, :, 4] = shear
        offset += n
    vm = von_mises(stresses)
    rvm = von_mises(reference)
    np.savez_compressed(
        path, coordinates=coords, displacement=displacements.reshape(-1,3),
        stress_sample_coordinates=physical, fe_stress=stresses,
        beam_reference_stress=reference, fe_von_mises=vm,
        beam_reference_von_mises=rvm,
        stress_sampling_scheme=np.array(STRESS_SAMPLING_SCHEME),
        stress_inset_fraction=np.array(STRESS_INSET_FRACTION),
        order=np.array(mesh.order), divisions=np.asarray(div),
        size=mesh.size, origin=mesh.origin,
        density=np.array(density), end_force=np.array(end_force))
    return {
        'field_archive':str(path),
        'stress_sampling':f'{STRESS_SAMPLING_SCHEME}: 8 direct FE samples per visualization subcell, fraction {STRESS_INSET_FRACTION}',
        'sampled_stress_points':int(stresses.shape[0]*8),
        'visualization_subcells':int(stresses.shape[0]),
        'fe_max_tension_zz_psi':float(np.maximum(stresses[...,2],0).max()),
        'fe_max_compression_zz_psi':float(np.maximum(-stresses[...,2],0).max()),
        'fe_max_von_mises_psi':float(vm.max()),
        'max_displacement_magnitude_in':float(np.linalg.norm(displacements.reshape(-1,3),axis=1).max()),
        'beam_reference_max_tension_zz_psi':float(np.maximum(reference[...,2],0).max()),
        'beam_reference_max_compression_zz_psi':float(np.maximum(-reference[...,2],0).max()),
        'beam_reference_max_von_mises_psi':float(rvm.max()),
        'sample_extrema_not_global_maxima':True,
        'smoke_test_not_validated_stress':bool(np.prod(div)<=16),
    }
