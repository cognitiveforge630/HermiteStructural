# BUILD AND VERIFICATION REPORT — 2026-10-29

## Origin and scope

- Project: **HermiteStructural_CubicHigherOrder_OOP**.
- Created as a separate copy of user's uploaded `HermiteStructural_CubicHigherOrder(2).zip`.
- SHA-256 of **input ZIP**: `2bf3af6af23e3e70bee691f7376f37041726829b804124344e72342232656c69`.
- Original archive not edited. Historical `.git/` metadata excluded from new
  distribution, rather than creating a new Git repo automatically.
- OOP refactor covers **only conventional C0 nodal Hex27/Hex64**.
  Older Hermite recovered-rotation solvers and visualization are preserved.
- Mathematical status remains **research / NOT a certified error bound**.

## Implemented

1. `core/elements/nodal_hex.py`: `HexElement`, `Hex27Element`,
   `Hex64Element`, `HexFacet`; 8 named corners, 6 named faces, cached
   Gauss locations, 3D volume and 2D facet quadrature, instance-local Fe,
   local Ke build/getter, and readable `_add_consistent_nodal_load` helper.
2. `core/mesh/nodal_hex.py`: lazy structured mesh iterator and global
   connectivity preserving original tensor node and DOF indexing.
3. `core/solvers/nodal_assembly.py`: 8-element-batch COO-to-CSR global
   assembly; the same dense Ke is shared for equivalent affine cells;
   no explicit destructors or dense Ke array for each cell.
4. `tools/p_order_study.py`: delegates to OOP classes and preserves command
   flags, result keys, gravity/face traction and section-integrated stress.
   Mathematical/stress/report code is split into named functions.
5. New `tests/test_nodal_oop.py`, test-only frozen legacy kernel,
   before-refactor solved baseline, `tools/demo_nodal_oop.py`,
   `tools/benchmark_nodal_assembly.py`, OOP tutorial, checklist, and VS tasks.

## Validation performed

- On the original 2026-10-29 OOP refactor: `python -m unittest discover -s tests -v`: **90 tests run: 88 passed, 2 skipped**. A later audit on 2026-10-29 passed **100 tests, 2 skipped**, before adding the additional documentation/graph regression tests; the final 2026-10-29 revision passed **104 tests, 2 skipped**; see [AUDIT_20261029.md](AUDIT_20261029.md) for current validation.
  (graphics driver offscreen smoke and unavailable PyVista integration).
- `python -m unittest tests.test_nodal_oop tests.test_p_order_study`: **14 passed**.
- `python -m tools.p_order_study --mode verify --order 2`: **patch and one-cell checks pass**.
- `python -m tools.p_order_study --mode verify --order 3`: **patch and one-cell checks pass**.
- `python -m tools.p_order_study --mode plan --order 3 --nx 3 --ny 5 --nz 25`: plan only; **no large solve**.
- Local Ke, gravity and +Z traction against frozen original kernel;
  all six facet tractions and physical moment resultants; global K/F
  against independently assembled legacy COO; baseline 3x3x5 solved
  displacement/section metrics for both orders.
- Python bytecode compilation succeeded for `core`, `tools`, `tests`, `tutorials`.

## Small-mesh assembly memory observations

These are two separate subprocess runs per method on this Linux environment,
*assembly only* (not solver factorization):

| Study | Legacy-style COO tracemalloc peak | OOP-batched tracemalloc peak | Legacy Linux RSS peak | OOP Linux RSS peak | Legacy/OOP elapsed |
|---|---:|---:|---:|---:|---:|
| Hex27, 16 cells | 4.96 MiB | 3.09 MiB | 126.8 MiB | 125.4 MiB | 0.054s / 0.049s |
| Hex64, 16 cells | 27.54 MiB | 17.78 MiB | 149.5 MiB | 142.4 MiB | 0.246s / 0.253s |

This is encouraging, but **not proof of full-size memory improvement** or
cross-machine performance. RSS includes interpreter/native overhead and
`tracemalloc` may not fully reflect native memory. Sample size is small,
CPU and allocations are environment-dependent, and the sparse solve itself
can dominate memory. Stored CSR `nnz` differs from the old global COO
coalescing path because exact zero entries and FP cancellation are handled
differently. Numeric matrix and physics equivalence are checked directly.

## Known limits / open TODOs

- Full **192-cell** Hex64 benchmark and peak RSS comparison not done.
- No Windows 11 run or GPU/VTK GUI visualization verification done here.
- The new element supports **axis-aligned affine hexahedra**, constant
  vector traction in global XYZ, and global-axes orthotropic constitutive D;
  non-affine/curved cells and variable pressure require additional work.
- No Python destructor is needed; finite element instances are short lived
  and K is assembled in bounded batches.
- The research tool reports `certified_error_bound: false` and has no
  mathematical a posteriori certificate.

See [`docs/OOP_TODO.md`](OOP_TODO.md) for the full checked/unchecked plan.
