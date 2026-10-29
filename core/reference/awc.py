"""Signed end-load plus distributed-weight beam-theory displacement comparison."""
import numpy as np
from core.reference.beam_theory import E, I, L, P, cantilever_point_load_deflection
DEFAULT_WEIGHT_DENSITY = 0.0174  # lbf/in^3
CROSS_SECTION_WIDTH_IN = 4.0
CROSS_SECTION_DEPTH_IN = 8.0
REFERENCE_LONGITUDINAL_E = E

def awc_reference_components(weight_density=DEFAULT_WEIGHT_DENSITY, *, end_force=P):
    """Signed Euler-Bernoulli cantilever deflections (inches), identical load inputs.

    These are classical beam-theory references, not an independent certification.
    The transverse loading axis is global -y and the beam axis is global +z.
    """
    if not np.isfinite(weight_density) or weight_density <= 0:
        raise ValueError('Tutorial 001 requires positive finite self-weight')
    if not np.isclose(I, CROSS_SECTION_WIDTH_IN*CROSS_SECTION_DEPTH_IN**3/12):
        raise ValueError('The reference inertia does not match the FEM cross section')
    end = cantilever_point_load_deflection(end_force, L, REFERENCE_LONGITUDINAL_E, I)
    q = weight_density * CROSS_SECTION_WIDTH_IN * CROSS_SECTION_DEPTH_IN
    gravity = -q * L**4 / (8.0 * REFERENCE_LONGITUDINAL_E * I)
    return {'end_load_uy_in': float(end), 'self_weight_uy_in': float(gravity),
            'total_uy_in': float(end + gravity),
            'weight_per_length_lbf_in': float(q),
            'total_weight_lbf': float(q * L)}


def combined_reference_tip(weight_density=DEFAULT_WEIGHT_DENSITY):
    return awc_reference_components(weight_density)['total_uy_in']
