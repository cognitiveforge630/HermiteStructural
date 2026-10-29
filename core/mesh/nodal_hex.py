"""Structured Hex27/Hex64 connectivity; no dense element array retained."""
from __future__ import annotations

import itertools
import numpy as np
from core.elements.nodal_hex import element_type


def local_nodes(order, ix, iy, iz, shape):
    """Global node IDs and XYZ-interleaved DOFs, local k (Z) fastest."""
    sy, sz = shape[1:]
    ids = np.fromiter((
        ((order*ix+i)*sy+order*iy+j)*sz+order*iz+k
        for i, j, k in itertools.product(range(order+1), repeat=3)
    ), dtype=np.int64, count=(order+1)**3)
    return ids, (ids[:, None]*3 + np.arange(3)).ravel()


class StructuredNodalHexMesh:
    """Generate one short-lived element at a time, sharing reference Ke."""
    def __init__(self, order, divisions, size, origin=None):
        self.order = int(order)
        self.element_class = element_type(order)
        self.divisions = tuple(int(x) for x in divisions)
        self.size = np.asarray(size, dtype=float)
        self.origin = np.asarray((0., 0., 0.) if origin is None else origin, dtype=float)
        if len(self.divisions) != 3 or min(self.divisions) < 1:
            raise ValueError('Three positive integer division counts required')
        if self.size.shape != (3,) or not np.isfinite(self.size).all() or np.any(self.size <= 0):
            raise ValueError('All physical dimensions must be finite and positive')
        if self.origin.shape != (3,) or not np.isfinite(self.origin).all():
            raise ValueError('Origin must have three finite coordinates')
        self.node_shape = tuple(self.order*n+1 for n in self.divisions)
        self.element_extent = self.size / np.asarray(self.divisions)

    @property
    def coordinates(self):
        grids = [np.linspace(self.origin[i], self.origin[i]+self.size[i], self.node_shape[i])
                 for i in range(3)]
        return np.array(np.meshgrid(*grids, indexing='ij')).reshape(3, -1).T

    @property
    def global_dof_count(self):
        return 3*int(np.prod(self.node_shape))

    def iter_elements(self, constitutive, *, quadrature_order=None, shared_stiffness=None):
        for ix, iy, iz in itertools.product(*(range(x) for x in self.divisions)):
            node_ids, _ = local_nodes(self.order, ix, iy, iz, self.node_shape)
            origin = self.origin + self.element_extent * (ix, iy, iz)
            yield (ix, iy, iz), self.element_class(
                extent=self.element_extent, constitutive=constitutive, node_ids=node_ids,
                origin=origin, quadrature_order=quadrature_order,
                shared_stiffness=shared_stiffness)
