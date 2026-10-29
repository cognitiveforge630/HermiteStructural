"""Readable, affine, tensor-product C0 Hex27/Hex64 finite elements.

This module is **not** the recovered-rotation Hermite element family.  Each
nodal element has three displacement DOFs per Lagrange interpolation node.
Parent coordinates (xi, eta, zeta) align with physical (X, Y, Z); Y is up,
X is right, Z increases from the back/clamp to the front/free end.

Only *axis-aligned affine* elements are supported: the constant Jacobian is
J = diag(dx/2, dy/2, dz/2). Curved or skewed cells require a different
isoparametric mapping, face Jacobian, and outward-normal calculation.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import itertools
from typing import Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss


FACET_NAMES = ('-Z', '+Z', '-Y', '+Y', '-X', '+X')
AXIS_INDEX = {'X': 0, 'Y': 1, 'Z': 2}


def lagrange_1d(p: int, s: float):
    """Cardinal Lagrange values and derivatives at one parent coordinate."""
    if p not in (2, 3):
        raise ValueError('Only Hex27 (p=2) and Hex64 (p=3) are implemented')
    nodes = np.linspace(-1.0, 1.0, p + 1)
    values = np.ones(p + 1)
    derivatives = np.zeros(p + 1)
    for i in range(p + 1):
        for j in range(p + 1):
            if i != j:
                values[i] *= (s - nodes[j]) / (nodes[i] - nodes[j])
        for j in range(p + 1):
            if i == j:
                continue
            factor = 1.0 / (nodes[i] - nodes[j])
            for k in range(p + 1):
                if k != i and k != j:
                    factor *= (s - nodes[k]) / (nodes[i] - nodes[k])
            derivatives[i] += factor
    return values, derivatives


def tensor_shapes(p: int, parent):
    """(N, dN/dparent) with (i, j, k) tensor order, k fastest."""
    n = p + 1
    one_d = [lagrange_1d(p, float(t)) for t in parent]
    N = np.empty(n**3)
    derivatives = np.empty((3, n**3))
    for m, (i, j, k) in enumerate(itertools.product(range(n), repeat=3)):
        N[m] = one_d[0][0][i] * one_d[1][0][j] * one_d[2][0][k]
        derivatives[0, m] = one_d[0][1][i] * one_d[1][0][j] * one_d[2][0][k]
        derivatives[1, m] = one_d[0][0][i] * one_d[1][1][j] * one_d[2][0][k]
        derivatives[2, m] = one_d[0][0][i] * one_d[1][0][j] * one_d[2][1][k]
    return N, derivatives


def strain_matrix(dN):
    """Engineering strain B; Voigt order xx, yy, zz, xy, yz, xz."""
    nodes = dN.shape[1]
    B = np.zeros((6, 3 * nodes))
    B[0, 0::3] = dN[0]
    B[1, 1::3] = dN[1]
    B[2, 2::3] = dN[2]
    B[3, 0::3] = dN[1]
    B[3, 1::3] = dN[0]
    B[4, 1::3] = dN[2]
    B[4, 2::3] = dN[1]
    B[5, 0::3] = dN[2]
    B[5, 2::3] = dN[0]
    return B


@lru_cache(maxsize=16)
def gauss_rule(order: int):
    """Read-only 1D Gaussian abscissas and weights; reused across cells."""
    if order < 1:
        raise ValueError('Quadrature order must be >= 1')
    points, weights = leggauss(order)
    points.flags.writeable = False
    weights.flags.writeable = False
    return points, weights


@lru_cache(maxsize=16)
def volume_shape_rule(p: int, q: int):
    """Cached reference-space (N, derivative, 3D weight) at q^3 points."""
    points, weights = gauss_rule(q)
    data = []
    for i, xi in enumerate(points):
        for j, eta in enumerate(points):
            for k, zeta in enumerate(points):
                N, derivative = tensor_shapes(p, (xi, eta, zeta))
                N.flags.writeable = False
                derivative.flags.writeable = False
                data.append((N, derivative, float(weights[i]*weights[j]*weights[k])))
    return tuple(data)


@lru_cache(maxsize=36)
def face_shape_rule(p: int, q: int, name: str):
    """Cached reference-space shape values and 2D weights for one named face."""
    if name not in FACET_NAMES:
        raise KeyError(f'Unknown facet {name!r}; expected {FACET_NAMES}')
    axis = AXIS_INDEX[name[1]]
    tangents = tuple(i for i in range(3) if i != axis)
    sign = +1.0 if name[0] == '+' else -1.0
    points, weights = gauss_rule(q)
    data = []
    for i, first in enumerate(points):
        for j, second in enumerate(points):
            parent = np.empty(3)
            parent[axis] = sign
            parent[tangents[0]] = first
            parent[tangents[1]] = second
            N, derivatives = tensor_shapes(p, parent)
            N.flags.writeable = False
            derivatives.flags.writeable = False
            data.append((N, float(weights[i]*weights[j])))
    return tuple(data)


@lru_cache(maxsize=16)
def reference_volume_quadrature(order: int):
    """Parent XYZ Gauss locations (q^3,3), reference weights (q^3,)."""
    gp, gw = gauss_rule(order)
    xyz = np.array(list(itertools.product(gp, repeat=3)), dtype=float)
    weights = np.array([gw[i]*gw[j]*gw[k]
                        for i, j, k in itertools.product(range(order), repeat=3)])
    xyz.flags.writeable = False
    weights.flags.writeable = False
    return xyz, weights


@lru_cache(maxsize=96)
def reference_face_quadrature(order: int, name: str):
    """Parent XYZ locations, reference 2D weights: NO physical area factor."""
    if name not in FACET_NAMES:
        raise KeyError(name)
    gp, gw = gauss_rule(order)
    axis = AXIS_INDEX[name[1]]
    axes = tuple(i for i in range(3) if i != axis)
    side = 1.0 if name.startswith('+') else -1.0
    xyz = np.empty((order**2, 3))
    weights = np.empty(order**2)
    for a, (i, j) in enumerate(itertools.product(range(order), repeat=2)):
        xyz[a, axis] = side
        xyz[a, axes[0]] = gp[i]
        xyz[a, axes[1]] = gp[j]
        weights[a] = gw[i]*gw[j]
    xyz.flags.writeable = False
    weights.flags.writeable = False
    return xyz, weights


@dataclass(frozen=True)
class HexFacet:
    """One *local* face attached to a particular HexElement instance."""
    element: 'HexElement'
    name: str

    @property
    def axis(self) -> int:
        return AXIS_INDEX[self.name[1]]

    @property
    def normal(self):
        normal = np.zeros(3)
        normal[self.axis] = 1.0 if self.name.startswith('+') else -1.0
        return normal

    def get_gaussian_quadrature(self):
        """Cached parent (q^2,3) coordinates and 2D reference weights.

        Multiply weights by facet.area/4 to integrate on this affine face.
        """
        return reference_face_quadrature(self.element.quadrature_order, self.name)

    @property
    def local_node_indices(self):
        """Local positions in tensor (i,j,k) order, k fastest."""
        side = self.element.order if self.name.startswith('+') else 0
        indices = []
        for index, ijk in enumerate(itertools.product(range(self.element.order+1), repeat=3)):
            if ijk[self.axis] == side:
                indices.append(index)
        return tuple(indices)

    @property
    def node_ids(self):
        if self.element.node_ids is None:
            raise ValueError('This element has no global connectivity')
        return self.element.node_ids[list(self.local_node_indices)].copy()

    @property
    def area(self):
        tangent_axes = [axis for axis in range(3) if axis != self.axis]
        return float(np.prod(self.element.extent[tangent_axes]))

    def integrate_traction(self, traction):
        """Integrate constant physical traction (tx, ty, tz), add to element Fe."""
        return self.element.integrate_traction(self.name, traction)


class HexElement:
    """One affine C0 element: reusable stiffness plus *instance-local* loads.

    Coordinates are optional (for named topology only). For the mesh solver,
    pass global node_ids and optionally a previously computed `shared_stiffness`
    from an equivalent reference element, avoiding a dense Ke per cell.
    """
    ORDER: int | None = None

    def __init__(self, extent, constitutive, *, order=None, node_ids=None,
                 origin=(0., 0., 0.), quadrature_order=None,
                 shared_stiffness=None):
        self.order = self.ORDER if self.ORDER is not None else order
        if self.order not in (2, 3):
            raise ValueError('HexElement requires order 2 or 3')
        self.quadrature_order = (self.order + 1 if quadrature_order is None
                                 else int(quadrature_order))
        if self.quadrature_order < self.order + 1:
            raise ValueError('Reduced integration is not supported')
        self.extent = np.asarray(extent, dtype=float).copy()
        if self.extent.shape != (3,) or not np.isfinite(self.extent).all() or np.any(self.extent <= 0):
            raise ValueError('Affine element extent must have three positive finite lengths')
        self.origin = np.asarray(origin, dtype=float).copy()
        if self.origin.shape != (3,) or not np.isfinite(self.origin).all():
            raise ValueError('origin must be a finite XYZ point')
        self.D = np.asarray(constitutive, dtype=float)
        if self.D.shape != (6, 6) or not np.isfinite(self.D).all():
            raise ValueError('Constitutive matrix D must be finite, 6x6')
        self.node_ids = None if node_ids is None else np.asarray(node_ids, dtype=np.int64).copy()
        if self.node_ids is not None and self.node_ids.shape != (self.node_count,):
            raise ValueError('Wrong number of global node IDs')
        self._stiffness = shared_stiffness
        if shared_stiffness is not None and shared_stiffness.shape != (self.dof_count, self.dof_count):
            raise ValueError('Shared stiffness has incorrect size')
        self._force = np.zeros(self.dof_count)
        self.facets: Mapping[str, HexFacet] = {name: HexFacet(self, name) for name in FACET_NAMES}

    def get_gaussian_quadrature(self):
        """Cached parent (q^3,3) coordinates and 3D reference weights.

        Multiply weights by element volume/8 for this affine element.
        """
        return reference_volume_quadrature(self.quadrature_order)

    @property
    def node_count(self):
        return (self.order + 1) ** 3

    @property
    def dof_count(self):
        return 3 * self.node_count

    @property
    def dof_indices(self):
        if self.node_ids is None:
            raise ValueError('Cannot assemble an element without node connectivity')
        return (self.node_ids[:, None] * 3 + np.arange(3)).ravel()

    @property
    def local_node_coordinates(self):
        """XYZ physical nodal coordinates, in stiffness matrix order."""
        locations = np.array(list(itertools.product(np.linspace(0, 1, self.order+1), repeat=3)))
        return self.origin + locations * self.extent

    @property
    def corner_nodes(self):
        """Names use X right(+), Y upper(+), Z front(+) / back(-)."""
        n = self.order
        result = {}
        for yname, j in (('lower', 0), ('upper', n)):
            for zname, k in (('back', 0), ('front', n)):
                for xname, i in (('left', 0), ('right', n)):
                    local_id = (i*(n+1) + j)*(n+1) + k
                    result[f'{yname}_{zname}_{xname}'] = (
                        int(self.node_ids[local_id]) if self.node_ids is not None else local_id)
        return result

    def shape_functions(self, parent):
        return tensor_shapes(self.order, parent)[0]

    def shape_derivatives(self, parent, *, physical=False):
        derivative = tensor_shapes(self.order, parent)[1]
        return derivative / (self.extent/2)[:, None] if physical else derivative

    def get_stiffness_matrix(self):
        """Ke = integral B.T @ D @ B dV; lazily built or shared read-only."""
        if self._stiffness is None:
            self._stiffness = self.build_stiffness_matrix()
            self._stiffness.flags.writeable = False  # safely share with many cells
        return self._stiffness

    def build_stiffness_matrix(self):
        scale = self.extent / 2
        detJ = float(np.prod(scale))
        ke = np.zeros((self.dof_count, self.dof_count))
        for _, derivatives, reference_weight in volume_shape_rule(self.order, self.quadrature_order):
            b = strain_matrix(derivatives / scale[:, None])
            ke += (b.T @ self.D @ b) * (reference_weight * detJ)
        return ke

    def clear_force_vector(self):
        self._force.fill(0.0)

    def add_force_vector(self, force):
        """Add an already integrated element load vector (e.g., shared gravity)."""
        f = np.asarray(force, dtype=float)
        if f.shape != (self.dof_count,):
            raise ValueError('Force vector has incorrect dimension')
        self._force += f

    def get_force_vector(self):
        """Get a COPY so callers cannot mutate an element's internal loads."""
        return self._force.copy()

    def _add_consistent_nodal_load(self, shape_values, vector_load, integration_weight):
        """Accumulate nodal XYZ forces, keeping the DOF stride out of lessons.

        Each interpolation node receives N_a * (three-component load) * w.
        ``reshape(-1, 3)`` is the sole low-level mapping from the 1D local
        FE force vector to its (node, XYZ-direction) view.
        """
        nodal_xyz = self._force.reshape(self.node_count, 3)
        nodal_xyz += shape_values[:, None] * vector_load[None, :] * integration_weight

    def integrate_body_force(self, force_per_volume):
        """Fe += integral N.T body_force dV; e.g. (0,-weight_density,0)."""
        body = np.asarray(force_per_volume, dtype=float)
        if body.shape != (3,) or not np.isfinite(body).all():
            raise ValueError('Body force needs three finite components')
        detJ = float(np.prod(self.extent/2))
        for N, _, w in volume_shape_rule(self.order, self.quadrature_order):
            physical_weight = w * detJ
            self._add_consistent_nodal_load(N, body, physical_weight)
        return self.get_force_vector()

    def integrate_traction(self, facet: str | HexFacet, traction):
        """Fe += integral_on_facet N.T traction dA for a constant vector.

        The two tangential Gaussian axes are integrated independently; the
        normal parent coordinate is fixed to +/-1 on the named facet.
        """
        name = facet.name if isinstance(facet, HexFacet) else facet
        if name not in self.facets or (isinstance(facet, HexFacet) and facet.element is not self):
            raise ValueError(f'Facet must be one of {FACET_NAMES} on this element')
        force = np.asarray(traction, dtype=float)
        if force.shape != (3,) or not np.isfinite(force).all():
            raise ValueError('Traction needs three finite components')
        normal_axis = AXIS_INDEX[name[1]]
        tangents = tuple(axis for axis in range(3) if axis != normal_axis)
        surface_jacobian = float(np.prod((self.extent/2)[list(tangents)]))
        for N, reference_weight in face_shape_rule(self.order, self.quadrature_order, name):
            area_weight = reference_weight * surface_jacobian
            self._add_consistent_nodal_load(N, force, area_weight)
        return self.get_force_vector()


class Hex27Element(HexElement):
    ORDER = 2


class Hex64Element(HexElement):
    ORDER = 3


def element_type(order: int):
    if order == 2:
        return Hex27Element
    if order == 3:
        return Hex64Element
    raise ValueError('Only Hex27 and Hex64 are implemented')
