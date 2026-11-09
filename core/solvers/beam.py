"""Build the three active rotation-recovery Hermite cantilever formulations.

All three use the same geometry, material, load convention and six nodal DOFs.
"""
from core.elements import RecoveredRotationHex8, QuadraticRotationHex8, CubicRotationHex8
from core.materials.redwood import redwood_orthotropic_material
from core.reference.beam_theory import E

REFERENCE_LONGITUDINAL_E = E

def build_beam(E_longitudinal=REFERENCE_LONGITUDINAL_E, nx=5, ny=9, nz=121, *,
               formulation='quadratic_rotation_hermite'):
    element_class = {
        "linear_rotation_hermite": RecoveredRotationHex8,
        "quadratic_rotation_hermite": QuadraticRotationHex8,
        "cubic_rotation_hermite": CubicRotationHex8,
    }.get(formulation)
    if element_class is None:
        raise ValueError(f"Unknown active Hermite formulation: {formulation}")
    material = redwood_orthotropic_material(E_longitudinal)
    return element_class(
        Lx=4.0,
        Ly=8.0,
        Lz=120.0,
        nx=nx,
        ny=ny,
        nz=nz,
        E=E_longitudinal,
        nu=0.3,
        gamma=0.0,
        w=0.0,
        gauss_order=4,
        orthotropic_constants=material,
    )

