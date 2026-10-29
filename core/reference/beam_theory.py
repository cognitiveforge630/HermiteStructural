"""Classical beam-theory reference for this reproducible Redwood example.

Reference formula:
    delta_max = P * L**3 / (3 * E * I)

Units:
    P = concentrated load at the free end, lb
    L = beam length, in
    E = modulus of elasticity, psi = lb/in^2
    I = area moment of inertia, in^4
    delta_max = free-end deflection, in
"""


from core.materials.design_limits import REDWOOD_BASIS
E = REDWOOD_BASIS.E_longitudinal_psi  # psi; unadjusted nominal grade reference
L = 120.0  # in

# Rectangular cross section from the model.
B = 4.0  # in, cross-section width in x
H = 8.0  # in, cross-section depth in y

# Strong-axis bending inertia for transverse y loading along a z-span beam.
I = B * H**3 / 12.0

# Total resultant of the uniform end-face traction. The solver assembles its
# consistent load vector by integrating the actual Hermite displacement operator
# over the z=L face; it does NOT simply divide P equally among end nodes.
P = -12.8  # lb


def cantilever_point_load_deflection(load, length, youngs_modulus, inertia):
    return load * length**3 / (3.0 * youngs_modulus * inertia)

