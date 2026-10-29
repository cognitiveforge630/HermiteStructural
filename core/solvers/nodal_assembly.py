"""Bounded-buffer sparse matrix assembly for nodal high-order hexahedra."""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


class GlobalAssembler:
    """Accumulate global K,F in bounded batches of local COO contributions.

    Large full-mesh arrays of (row, column, value) are never retained. Ke
    can be shared by all equivalent affine cells. Matrix addition of batch
    CSR blocks may change the last bits from another assembly order; tests
    check the actual physics, not equality of floating-point summation order.
    """
    def __init__(self, ndof: int, *, batch_elements: int = 8):
        if ndof < 1 or batch_elements < 1:
            raise ValueError('ndof and batch_elements must be positive')
        self.ndof = ndof
        self.batch_elements = batch_elements
        self.force = np.zeros(ndof)
        self._rows = []
        self._cols = []
        self._values = []
        self._pending = 0
        self._matrix = csr_matrix((ndof, ndof))
        self._finalized = False

    def add_element(self, element):
        if self._finalized:
            raise RuntimeError('Cannot add elements after finalize()')
        dofs = element.dof_indices
        n = len(dofs)
        stiffness = element.get_stiffness_matrix()
        if stiffness.shape != (n, n) or not np.all((dofs >= 0) & (dofs < self.ndof)):
            raise ValueError('Invalid element stiffness or global DOFs')
        self.force[dofs] += element.get_force_vector()
        self._rows.append(np.repeat(dofs, n))
        self._cols.append(np.tile(dofs, n))
        self._values.append(stiffness.ravel())
        self._pending += 1
        if self._pending >= self.batch_elements:
            self.flush()

    def flush(self):
        if not self._pending:
            return
        batch = coo_matrix((np.concatenate(self._values),
                            (np.concatenate(self._rows), np.concatenate(self._cols))),
                           shape=(self.ndof, self.ndof)).tocsr()
        batch.sum_duplicates()
        self._matrix += batch
        self._rows.clear()
        self._cols.clear()
        self._values.clear()
        self._pending = 0

    def finalize(self):
        self.flush()
        self._matrix.sum_duplicates()
        self._finalized = True
        return self._matrix, self.force
