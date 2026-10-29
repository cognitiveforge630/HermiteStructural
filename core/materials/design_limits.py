"""Provisional Redwood *reference-grade* screening, NOT code design approval.

Reference: Redwood Technical Data Sheet, Structural Design Values (2022),
Redwood dimension lumber 2-4 in thick, Select Structural.
https://getredwood.com/wp-content/uploads/2024/07/2022-Structural-Design-Values.pdf

These are unadjusted published design values, not ultimate material strengths.
This idealized 4 x 8 in solid has not been assigned a real grade stamp.
Moisture, size, duration, treatment, stability, connection, and application
adjustments must be assessed independently for actual structural use.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class RedwoodDesignBasis:
    grade: str = "Redwood Select Structural (dimension lumber, nominal)"
    E_longitudinal_psi: float = 1_100_000.0
    bending_psi: float = 1100.0
    tension_parallel_psi: float = 625.0
    shear_parallel_psi: float = 160.0
    compression_perpendicular_psi: float = 425.0
    compression_parallel_psi: float = 1100.0
    deflection_ratio: float = 360.0  # chosen serviceability example, NOT Redwood failure
    source: str = "https://getredwood.com/wp-content/uploads/2024/07/2022-Structural-Design-Values.pdf"


REDWOOD_BASIS = RedwoodDesignBasis()


def reference_load_screen(*, tip_force_lbf, density_lbf_in3, length_in=120., width_in=4., depth_in=8.,
                          modulus_psi=REDWOOD_BASIS.E_longitudinal_psi, basis=REDWOOD_BASIS):
    """Conservative textbook beam checks BEFORE solving; reference screening only.

    Assumes downward end load and uniform downward gravity; checks maximum
    absolute bending, tip deflection, and rectangular-beam shear. No member
    stability, connection, perpendicular tension, or full NDS design check.
    """
    if not all(math.isfinite(x) for x in (tip_force_lbf,density_lbf_in3,length_in,width_in,depth_in,modulus_psi)):
        raise ValueError("Inputs must be finite")
    if tip_force_lbf > 0 or density_lbf_in3 < 0 or min(length_in,width_in,depth_in,modulus_psi)<=0:
        raise ValueError("This downward-load cantilever check requires P<=0, weight density>=0 and positive dimensions/E")
    q = density_lbf_in3 * width_in * depth_in
    I = width_in*depth_in**3/12.
    M = -tip_force_lbf*length_in + q*length_in**2/2.
    peak_bending = M*(depth_in/2.)/I
    tip = abs(tip_force_lbf)*length_in**3/(3.*modulus_psi*I)+q*length_in**4/(8.*modulus_psi*I)
    shear = 1.5*(-tip_force_lbf+q*length_in)/(width_in*depth_in)
    lim = length_in/basis.deflection_ratio
    checks = {
        "bending": {"demand":peak_bending,"limit":basis.bending_psi,"units":"psi"},
        "parallel_tension": {"demand":peak_bending,"limit":basis.tension_parallel_psi,"units":"psi"},
        "parallel_compression": {"demand":peak_bending,"limit":basis.compression_parallel_psi,"units":"psi"},
        "parallel_shear": {"demand":shear,"limit":basis.shear_parallel_psi,"units":"psi"},
        "example_tip_serviceability": {"demand":tip,"limit":lim,"units":"in"},
    }
    for rec in checks.values():
        rec["utilization"] = rec["demand"]/rec["limit"]
        rec["passes"] = rec["utilization"]<=1.0
    return {"status":"PASS_PRELIMINARY" if all(c["passes"] for c in checks.values()) else "FAIL_PRELIMINARY",
            "checks":checks,"weight_per_length_lbf_in":q,"assumed_nominal_grade":basis.grade,
            "limits_not_failure_strength":True,"not_a_structural_design_approval":True}
