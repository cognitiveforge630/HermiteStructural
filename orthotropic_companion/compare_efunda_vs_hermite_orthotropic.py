"""
Compare an eFunda-style orthotropic compliance matrix against the
HermiteStructural OrthotropicElasticity compliance matrix.

Stress/strain order used here:
    xx, yy, zz, xy, yz, xz

Shear convention used here:
    engineering shear strains gamma_xy, gamma_yz, gamma_xz
so the shear compliance terms are 1/G.
"""

import numpy as np

np.set_printoptions(precision=12, suppress=False, linewidth=160)

REDWOOD_CLEAR_WOOD_RATIOS = {
    # Redwood-like orthotropic ratios with wood L along beam/global z-axis.
    "E_T_over_E_L": 0.089,
    "E_R_over_E_L": 0.087,
    "G_LR_over_E_L": 0.066,
    "G_LT_over_E_L": 0.077,
    "G_RT_over_E_L": 0.011,
    "nu_LR": 0.360,
    "nu_LT": 0.346,
    "nu_RT": 0.373,
}

# This is the calibrated longitudinal modulus printed by the current Hermite tutorial.
E_LONGITUDINAL = 1_680_668.967550  # psi
# For the seed/reference run, change this to 2_000_000.0 psi.


def redwood_constants(E_longitudinal: float):
    """Map wood axes R,T,L onto solver axes x,y,z."""
    r = REDWOOD_CLEAR_WOOD_RATIOS

    Ex = r["E_R_over_E_L"] * E_longitudinal  # R -> x
    Ey = r["E_T_over_E_L"] * E_longitudinal  # T -> y
    Ez = E_longitudinal                      # L -> z

    Gxy = r["G_RT_over_E_L"] * E_longitudinal
    Gyz = r["G_LT_over_E_L"] * E_longitudinal
    Gxz = r["G_LR_over_E_L"] * E_longitudinal

    # These are the three Poisson entries stored by the Hermite implementation.
    # nuxy is direct RT. nuxz and nuyz are the reciprocal forms RL and TL
    # derived from the tabulated LR and LT values.
    nuxy = r["nu_RT"]
    nuxz = r["nu_LR"] * Ex / Ez
    nuyz = r["nu_LT"] * Ey / Ez

    # These are the reciprocal partners used by the full/eFunda-style display.
    nuyx = nuxy * Ey / Ex
    nuzx = nuxz * Ez / Ex
    nuzy = nuyz * Ez / Ey

    return {
        "Ex": Ex, "Ey": Ey, "Ez": Ez,
        "Gxy": Gxy, "Gyz": Gyz, "Gxz": Gxz,
        "nuxy": nuxy, "nuxz": nuxz, "nuyz": nuyz,
        "nuyx": nuyx, "nuzx": nuzx, "nuzy": nuzy,
    }


def S_hermite(c):
    """Compliance matrix as written in OrthotropicElasticity."""
    Ex, Ey, Ez = c["Ex"], c["Ey"], c["Ez"]
    Gxy, Gyz, Gxz = c["Gxy"], c["Gyz"], c["Gxz"]
    nuxy, nuxz, nuyz = c["nuxy"], c["nuxz"], c["nuyz"]

    return np.array([
        [1.0 / Ex,      -nuxy / Ex,    -nuxz / Ex,    0.0,        0.0,        0.0],
        [-nuxy / Ex,    1.0 / Ey,      -nuyz / Ey,    0.0,        0.0,        0.0],
        [-nuxz / Ex,   -nuyz / Ey,     1.0 / Ez,      0.0,        0.0,        0.0],
        [0.0,           0.0,            0.0,          1.0 / Gxy,  0.0,        0.0],
        [0.0,           0.0,            0.0,          0.0,        1.0 / Gyz,  0.0],
        [0.0,           0.0,            0.0,          0.0,        0.0,        1.0 / Gxz],
    ], dtype=float)


def S_efunda_full(c):
    """
    eFunda/full display-style compliance matrix.

    This writes both directions of each Poisson pair explicitly:
        S12 = -nu_yx/Ey, S21 = -nu_xy/Ex
        S13 = -nu_zx/Ez, S31 = -nu_xz/Ex
        S23 = -nu_zy/Ez, S32 = -nu_yz/Ey

    When reciprocal symmetry is enforced, these pairs are equal.
    """
    Ex, Ey, Ez = c["Ex"], c["Ey"], c["Ez"]
    Gxy, Gyz, Gxz = c["Gxy"], c["Gyz"], c["Gxz"]
    nuxy, nuxz, nuyz = c["nuxy"], c["nuxz"], c["nuyz"]
    nuyx, nuzx, nuzy = c["nuyx"], c["nuzx"], c["nuzy"]

    return np.array([
        [1.0 / Ex,      -nuyx / Ey,    -nuzx / Ez,    0.0,        0.0,        0.0],
        [-nuxy / Ex,    1.0 / Ey,      -nuzy / Ez,    0.0,        0.0,        0.0],
        [-nuxz / Ex,   -nuyz / Ey,     1.0 / Ez,      0.0,        0.0,        0.0],
        [0.0,           0.0,            0.0,          1.0 / Gxy,  0.0,        0.0],
        [0.0,           0.0,            0.0,          0.0,        1.0 / Gyz,  0.0],
        [0.0,           0.0,            0.0,          0.0,        0.0,        1.0 / Gxz],
    ], dtype=float)


def main():
    c = redwood_constants(E_LONGITUDINAL)

    print("Orthotropic constants used for the comparison")
    print("wood axes mapped as R->x, T->y, L->z")
    for key in ["Ex", "Ey", "Ez", "Gxy", "Gyz", "Gxz", "nuxy", "nuxz", "nuyz", "nuyx", "nuzx", "nuzy"]:
        print(f"  {key:>4s} = {c[key]: .12g}")

    print("\nReciprocity checks")
    print(f"  nuxy/Ex = {c['nuxy']/c['Ex']:.16e}")
    print(f"  nuyx/Ey = {c['nuyx']/c['Ey']:.16e}")
    print(f"  nuxz/Ex = {c['nuxz']/c['Ex']:.16e}")
    print(f"  nuzx/Ez = {c['nuzx']/c['Ez']:.16e}")
    print(f"  nuyz/Ey = {c['nuyz']/c['Ey']:.16e}")
    print(f"  nuzy/Ez = {c['nuzy']/c['Ez']:.16e}")

    S_e = S_efunda_full(c)
    S_h = S_hermite(c)

    print("\nS matrix: eFunda-style full Poisson-pair form")
    print(S_e)

    print("\nS matrix: HermiteStructural custom/compact form")
    print(S_h)

    diff = S_e - S_h
    print("\nS_eFunda - S_Hermite")
    print(diff)
    print(f"\nmax(abs(S_eFunda - S_Hermite)) = {np.max(np.abs(diff)):.16e}")
    print(f"allclose? {np.allclose(S_e, S_h, rtol=0.0, atol=1e-18)}")

    D_e = np.linalg.inv(S_e)
    D_h = np.linalg.inv(S_h)
    Ddiff = D_e - D_h

    print("\nD matrix from eFunda-style S inverse")
    print(D_e)

    print("\nD matrix from HermiteStructural S inverse")
    print(D_h)

    print("\nD_eFunda - D_Hermite")
    print(Ddiff)
    print(f"\nmax(abs(D_eFunda - D_Hermite)) = {np.max(np.abs(Ddiff)):.16e}")
    print(f"allclose? {np.allclose(D_e, D_h, rtol=0.0, atol=1e-8)}")


if __name__ == "__main__":
    main()
