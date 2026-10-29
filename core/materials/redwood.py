"""Reference Redwood ratios; elastic moduli are prescribed, never calibrated."""
from core.elements.base import OrthotropicElasticity

REDWOOD_CLEAR_WOOD_RATIOS = {
    # Redwood-like orthotropic ratios with L along the beam/global z-axis.
    "E_T_over_E_L": 0.089,
    "E_R_over_E_L": 0.087,
    "G_LR_over_E_L": 0.066,
    "G_LT_over_E_L": 0.077,
    "G_RT_over_E_L": 0.011,
    "nu_LR": 0.360,
    "nu_LT": 0.346,
    "nu_RT": 0.373,
}


def redwood_orthotropic_material(E_longitudinal):
    """Return orthotropic constants with wood axes R->x, T->y, and L->z."""
    ratios = REDWOOD_CLEAR_WOOD_RATIOS
    E_longitudinal = float(E_longitudinal)
    Ex = ratios["E_R_over_E_L"] * E_longitudinal
    Ey = ratios["E_T_over_E_L"] * E_longitudinal
    Ez = E_longitudinal
    return OrthotropicElasticity(
        Ex=Ex,
        Ey=Ey,
        Ez=Ez,
        nuxy=ratios["nu_RT"],
        nuxz=ratios["nu_LR"] * Ex / Ez,
        nuyz=ratios["nu_LT"] * Ey / Ez,
        Gxy=ratios["G_RT_over_E_L"] * E_longitudinal,
        Gyz=ratios["G_LT_over_E_L"] * E_longitudinal,
        Gxz=ratios["G_LR_over_E_L"] * E_longitudinal,
    )

