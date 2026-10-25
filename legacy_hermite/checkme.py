"""Theoretical cantilever check for a point load at the free end.

Reference formula:
    delta_max = P * L**3 / (3 * E * I)

Units:
    P = concentrated load at the free end, lb
    L = beam length, in
    E = modulus of elasticity, psi = lb/in^2
    I = area moment of inertia, in^4
    delta_max = free-end deflection, in
"""


E = 2_000_000.0  # psi
L = 120.0  # in

# Rectangular cross section from the model.
B = 4.0  # in, cross-section width in x
H = 8.0  # in, cross-section depth in y

# Strong-axis bending inertia for transverse y loading along a z-span beam.
I = B * H**3 / 12.0

# Total point load at the free end. The solver distributes this same total P
# equally over the nodes on the z=L face.
P = -12.8  # lb


def cantilever_point_load_deflection(load, length, youngs_modulus, inertia):
    return load * length**3 / (3.0 * youngs_modulus * inertia)


def main():
    delta_theory = cantilever_point_load_deflection(P, L, E, I)

    print("Cantilever beam point-load check")
    print(f"  P: {P:.3f} lb")
    print(f"  L: {L:.3f} in")
    print(f"  E: {E:.3f} psi")
    print(f"  I: {I:.6f} in^4")
    print(f"  Theoretical free-end deflection: {delta_theory:.6e} in")


if __name__ == "__main__":
    main()
