"""Structured accuracy-study mesh construction, including optional z grading."""
import numpy as np
from solver import build_beam, REFERENCE_LONGITUDINAL_E


def build_accuracy_beam(mesh, formulation, modulus=REFERENCE_LONGITUDINAL_E, z_breaks=None):
    beam = build_beam(modulus, nx=mesh[0]+1, ny=mesh[1]+1, nz=mesh[2]+1,
                      formulation=formulation)
    if z_breaks is not None:
        z = np.asarray(z_breaks, dtype=float)
        if z.shape != (mesh[2]+1,) or not np.all(np.isfinite(z)) or np.any(np.diff(z) <= 0):
            raise ValueError('z_breaks must contain EZ+1 finite increasing coordinates')
        if z[0] != 0. or z[-1] != beam.Lz:
            raise ValueError('z_breaks must preserve the physical span [0,120]')
        if formulation == 'cauchy_hex27_reference':
            z = np.sort(np.r_[z, .5*(z[:-1]+z[1:])])
        original = np.unique(beam.grid_points[:, 2])
        beam.grid_points[:, 2] = z[np.searchsorted(original, beam.grid_points[:, 2])]
        if hasattr(beam, 'gradient_weights'):
            beam.gradient_weights = beam._gradient_recovery_weights()
            beam._element_transform_cache.clear()
    return beam
