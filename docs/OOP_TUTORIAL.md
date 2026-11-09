# Tutorial: Hex27 / Hex64 finite elements as Python objects

This tutorial is **independent conventional C0 Lagrange finite-element work**.
It does **not** rewrite the separate recovered-rotation Hermite Hex8 family.
Preserved mathematical convention: 3 displacement DOFs per node and engineering
Voigt strain order `(xx, yy, zz, xy, yz, xz)`.

## 1. Start at the element, not at global K

Run these **in order** from the project root in PowerShell (standard Python;
no virtual environment):

```powershell
python -m tools.demo_nodal_oop --order 3
python -m unittest tests.test_nodal_oop tests.test_p_order_study -v
python -m tools.p_order_study --mode verify --order 3
python -m tools.p_order_study --mode plan --order 3 --nx 3 --ny 5 --nz 25
# Small 16-cell study; not the expensive 192-cell study:
python -m tools.p_order_study --mode solve --order 3 --nx 3 --ny 3 --nz 5
```

## 2. Classes and ownership

| Python object | Responsibility |
|---|---|
| `HexElement` | Shared affine Hex27/Hex64 mathematics and instance-local Fe |
| `Hex27Element` | Quadratic 3 x 3 x 3 nodal element (27 nodes, 81 DOFs) |
| `Hex64Element` | Cubic 4 x 4 x 4 nodal element (64 nodes, 192 DOFs) |
| `HexFacet` | One of the six **named local faces**, its node IDs, normal, and 2D traction integral |
| `StructuredNodalHexMesh` | Shared-node coordinate grid and lazy iteration over physical cells |
| `GlobalAssembler` | Element K/F insertion in bounded COO buffers, periodic CSR flush |

`core/elements/nodal_hex.py` holds element mathematics. `core/mesh/nodal_hex.py`
holds mesh connectivity. `core/solvers/nodal_assembly.py` holds sparse assembly.
`tools/p_order_study.py` is the solved research comparison and retains its
legacy-compatible command-line interface and JSON keys. It delegates to OOP.

## 3. Coordinate names and node ordering

Our orientation is **+X = right**, **+Y = upper**, **+Z = front**, so the
cantilever's fixed end at `z=0` is the **back (-Z)** and its loaded free end
at `z=L` is the **front (+Z)**. The eight corner lookups use names such as
`upper_back_right`, `upper_front_left`, `lower_front_right`.

Local polynomial node ordering is `(i, j, k)` with `k` fastest, and each node
owns interleaved displacement DOFs `(ux, uy, uz)`:

```python
from core.elements.nodal_hex import Hex64Element
from core.materials.redwood import redwood_orthotropic_material
from core.reference.beam_theory import E

D = redwood_orthotropic_material(E).elastic_matrix()
element = Hex64Element((2., 4., 12.), D, node_ids=range(64))
print(element.corner_nodes['upper_back_right'])
print(element.facets['+Z'].node_ids)  # 16 local face interpolation nodes
print(element.facets['+Z'].normal)    # (0, 0, +1)
```

`facets` keys are **in this order**: `-Z, +Z, -Y, +Y, -X, +X`.
`facet.local_node_indices` returns indices into local `N` and `Ke`; if
connectivity is present, `facet.node_ids` returns global shared-node IDs.

## 4. Integration points are NOT interpolation nodes

Hex64's 64 Lagrange **interpolation nodes** are the `4x4x4` tensor grid.
The *Gaussian quadrature* for the default `q=4` uses 64 **volume points**
(4 per direction), and any face uses 16 **surface points** (4 per tangent
direction). Each Gaussian point has an integration weight; four one-dimensional
weights are not four coordinate directions.

Quadrature and reference shape evaluations are cached with `lru_cache` and
read-only NumPy arrays. Different `p`, `q`, or face names have separate cached
rules. The physical affine element has constant diagonal Jacobian:

\[
J = \mathrm{diag}(h_x/2, h_y/2, h_z/2),\qquad
\mathrm{d}V=\det(J)\,\mathrm{d}\xi\mathrm{d}\eta\mathrm{d}\zeta.
\]

One facet fixes its normal parent coordinate to `+1` or `-1`, then integrates
over two **tangent** parent directions; for a `+Z` face:

\[
\mathbf f_e^{(+Z)} = \int_{-1}^{1}\int_{-1}^{1}
\mathbf N(\xi,\eta,+1)^\mathsf T\mathbf t\,
(h_x/2)(h_y/2)\,\mathrm d\xi\mathrm d\eta.
\]

A face-load vector has `3*(p+1)^3` entries even though the face has only
`(p+1)^2` interpolation nodes; all non-face shape functions are zero on that
reference face (within roundoff). **The surface Jacobian is not the 3D
Jacobian**. Each of the other five faces uses its own two tangent directions.

## 5. Readable element-local loads

```python
end_face = element.facets['+Z']

# Gravity: force per unit physical volume (lbf/in^3).
element.integrate_body_force((0., -0.02, 0.))

# Constant downward traction: force per unit face area (lbf/in^2).
end_face.integrate_traction((0., -1.25, 0.))

local_F = element.get_force_vector()
local_K = element.get_stiffness_matrix()
```

These two integrations add into the **one local force vector**. `get_force_vector`
returns a copy. `clear_force_vector()` resets loads if a new load case is
needed. Tractions supplied here are **prescribed vectors in global XYZ**, not
pressure automatically converted from a face normal; for pressure use
`traction = -pressure * facet.normal`.

The stiffness is the standard finite-element linear-elastic bilinear form:

\[
K_e = \int_{\Omega_e} B^\mathsf{T} D B\,\mathrm{d}V,
\qquad \mathbf f_e^{g}=\int_{\Omega_e}N^\mathsf T \mathbf b\,\mathrm{d}V.
\]

`get_stiffness_matrix()` lazily builds a local `Ke`, or returns a shared
reference matrix for another **equivalent** affine cell. Equivalent here means
identical geometry, orientation, constitutive tensor, interpolation order,
quadrature, and DOF convention. Never share `Ke` across unlike cells.

## 6. Assemble the global system without thousands of dense matrices

```python
mesh = StructuredNodalHexMesh(3, (2, 2, 4), (2., 4., 12.))
master = Hex64Element(mesh.element_extent, D)
shared_Ke = master.get_stiffness_matrix()
gravity_Fe = master.integrate_body_force((0., -0.02, 0.))
assembler = GlobalAssembler(mesh.global_dof_count, batch_elements=8)

for (ix, iy, iz), element in mesh.iter_elements(D, shared_stiffness=shared_Ke):
    element.add_force_vector(gravity_Fe)
    if iz == mesh.divisions[2] - 1:
        element.facets['+Z'].integrate_traction((0., -1.25, 0.))
    assembler.add_element(element)

K, F = assembler.finalize()
```

For this snippet additionally import:

```python
from core.mesh.nodal_hex import StructuredNodalHexMesh
from core.solvers.nodal_assembly import GlobalAssembler
```

The code streams *element instances*; it does not create a list of all local
192-by-192 matrices. Assembly uses batches of at most 8 cells' COO triplets,
adds each block to a CSR global matrix, and discards the temporary arrays.
Python `__del__` is **not** used to force garbage collection. Actual peak RAM
also depends on the global sparse matrix, sparse solve/factorization, and
NumPy/SciPy native allocations; smaller batches do not bound *total* memory.

## 7. How the change was verified

- `tests/reference_legacy_p_order.py` is a **frozen extraction** of the
  pre-refactor shape/strain/element kernel from the provided ZIP; it is only
  for independent regression testing, not production.
- Local matrices and volume/`+Z` loads match to floating-point tolerances.
- On **every named face**, independently test traction force and moment for
  both Hex27 and Hex64, including oblique traction components.
- Test correct shared-face connectivity, Jacobian scale, body-force moments,
  quadrature cache reuse, and batch-assembled global matrices.
- Saved baseline: `docs/verification/legacy_baseline_3x3x5.json` from before
  any code changes. Test compares both modest 16-cell solves with it.

`stiffness_nnz` may not be the same integer as the old COO path: incremental
CSR assembly can eliminate exact stored zeros and produce subtly different
floating-point cancellation. The **assembled coefficients, physical quantities,
and equilibrium residuals** are the meaningful tests. See the baseline JSON;
do not compare only `nnz` to establish equivalence.

## 8. Limits (important)

1. Only axis-aligned *affine* positive-size cells, not skewed/curved Hexes.
2. Constant global-XYZ traction on a named facet; spatially varying pressure
   is not implemented here.
3. Orthotropic `D` must be expressed in the global axes used for this model.
4. This conventional C0 solver has only nodal translation DOFs. It cannot be
   substituted for recovered-rotation Hermite Hex8 without rederiving the DOFs.
5. A passing patch test, force/moment balance, or classical beam agreement
   does **not** constitute a certified 3D error bound. JSON stays
   `certified_error_bound: false`.
6. Heavy 192-cell Hex64 factorization is opt-in, not an automatic test.

Track completed and open tasks in [`OOP_TODO.md`](OOP_TODO.md).

## 9. Optional small-mesh memory experiment (not a full-size benchmark)

```powershell
python -m tools.benchmark_nodal_assembly --order 2
python -m tools.benchmark_nodal_assembly --order 3
```

Saved Linux results are in `docs/verification/assembly_benchmark_hex27_3x3x5_linux.json`
and the analogous `hex64` file. This deliberately benchmarks **K/F assembly
only**, not a solve, and launches each method in its own Python process.
The legacy-style mode accumulates whole-mesh COO arrays before conversion;
the OOP mode uses eight-element batches. Recorded `tracemalloc` estimates
may miss native allocations; Linux max-RSS reflects the whole process but
is environment-dependent. Compare physical K/F entries and resultant loads,
not stored `nnz` alone. A 192-cell Hex64 sparse factorization and Windows
peak-RAM study remain **uncompleted**.

## Independent traction integration test

Run `python -m tests.test_facet_traction_independent --show-results` from the
project root to compare the Gaussian facet integral against independent
exact-polynomial integration. See [unit testing](UNIT_TESTING.md).
