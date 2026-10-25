"""Tutorial 03: trace one element as it is assembled into global K.

This script isolates the core ideas inside ``FEMHermiteBeamRegion.build_global_K``:

1. Pick one element from the mesh.
2. Print its element index, node coordinates, and global DOF mapping.
3. Build that element's local ``Ke`` matrix.
4. Walk through the ``K_row``, ``K_col``, and ``K_data`` append operations that
   place local entries into the global sparse-matrix triplet lists.

Run:
    py tutorial_03_element_assembly_trace.py

Optional:
    py tutorial_03_element_assembly_trace.py --element 10 --max-entries 24
    py tutorial_03_element_assembly_trace.py --all-entries
    py tutorial_03_element_assembly_trace.py --elements 0 1 --skip-ke
"""

from __future__ import annotations

import argparse

import numpy as np

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


DOF_LABELS = ("ux", "uy", "uz", "phix", "phiy", "phiz")

MINI_BEAM_PARAMS = {
    "Lx": 4.0,
    "Ly": 8.0,
    "Lz": 6.0,
    "nx": 3,
    "ny": 3,
    "nz": 4,
    "E": 2_000_000.0,
    "nu": 0.3,
    "gamma": 0.0174,
    "w": -0.4,
}


def build_beam() -> FEMHermiteBeamRegion:
    """Use a small mesh so the printed lesson stays readable."""
    return FEMHermiteBeamRegion(**MINI_BEAM_PARAMS)


def print_global_context(beam: FEMHermiteBeamRegion) -> None:
    n_nodes = beam.grid_points.shape[0]
    n_dofs = n_nodes * beam.ndof_per_node

    print("Global assembly context")
    print(f"  grid nodes:       {n_nodes}")
    print(f"  elements:         {beam.get_num_elements()}")
    print(f"  dofs per node:    {beam.ndof_per_node} ({', '.join(DOF_LABELS)})")
    print(f"  global K shape:   {n_dofs} x {n_dofs}")
    print()


def print_element_identity(beam: FEMHermiteBeamRegion, elem_idx: int) -> None:
    node_indices = beam.elements[elem_idx]
    coords = beam.get_element_nodes(elem_idx)
    dof_grid = beam.get_element_global_dofs(elem_idx)
    dofs = dof_grid.ravel()

    print("=" * 79)
    print(f"Element trace: elem_idx = {elem_idx}")
    print()
    print("coords = self.get_element_nodes(elem_idx)")
    print("  These are the 8 physical node coordinates used to build J, B, and Ke.")
    print("  local_node  global_node        x        y        z")
    for local_node, (global_node, xyz) in enumerate(zip(node_indices, coords)):
        print(
            f"  {local_node:>10}  {global_node:>11}  "
            f"{xyz[0]:>7.2f}  {xyz[1]:>7.2f}  {xyz[2]:>7.2f}"
        )
    print()

    print("dofs = self.get_element_global_dofs(elem_idx).ravel()")
    print("  Before ravel(), the DOFs are arranged as 8 nodes x 6 DOFs:")
    print("  local_node  global_node        ux   uy   uz  phix phiy phiz")
    for local_node, (global_node, global_dofs) in enumerate(zip(node_indices, dof_grid)):
        dof_text = " ".join(f"{int(dof):>5}" for dof in global_dofs)
        print(f"  {local_node:>10}  {global_node:>11}  {dof_text}")
    print()

    print("  After ravel(), the same DOFs become one 48-entry local-to-global map:")
    print(f"  dofs.shape = {dofs.shape}")
    print(f"  dofs = {dofs.tolist()}")
    print()


def build_element_stiffness(beam: FEMHermiteBeamRegion, elem_idx: int) -> np.ndarray:
    coords = beam.get_element_nodes(elem_idx)
    Ke = np.zeros((48, 48))

    for i, xi in enumerate(beam.gauss_points):
        for j, eta in enumerate(beam.gauss_points):
            for k, zeta in enumerate(beam.gauss_points):
                w = beam.gauss_weights[i] * beam.gauss_weights[j] * beam.gauss_weights[k]
                B, detJ, _, _ = beam.build_B_matrix(xi, eta, zeta, coords)
                Ke += B.T @ beam.D @ B * w * detJ

    return Ke


def print_ke_summary(Ke: np.ndarray) -> None:
    print("Ke = np.zeros((48, 48)), then Gauss integration fills it")
    print("  Ke is local: row 0 means this element's first local DOF, not global DOF 0.")
    print(f"  Ke.shape:       {Ke.shape}")
    print(f"  nonzero count:  {np.count_nonzero(Ke)}")
    print(f"  min / max:      {Ke.min(): .6e} / {Ke.max(): .6e}")
    print()
    print("  Top-left 6 x 6 block of Ke:")
    for row in Ke[:6, :6]:
        print("  " + " ".join(f"{value:>12.4e}" for value in row))
    print()


def explain_triplet_lists() -> None:
    print("Why K_row, K_col, and K_data exist")
    print("  scipy.sparse.coo_matrix accepts three parallel lists:")
    print("    K_row[n]  = global row index for contribution n")
    print("    K_col[n]  = global column index for contribution n")
    print("    K_data[n] = numeric stiffness value for that row/column")
    print()
    print("  The nested loops visit every entry in the 48 x 48 local Ke block.")
    print("  For each local entry Ke[a, b], the global address is:")
    print("    global row = dofs[a]")
    print("    global col = dofs[b]")
    print("    value      = Ke[a, b]")
    print()


def trace_triplet_appends(
    beam: FEMHermiteBeamRegion,
    elem_idx: int,
    Ke: np.ndarray,
    max_entries: int,
) -> tuple[list[float], list[int], list[int]]:
    dofs = beam.get_element_global_dofs(elem_idx).ravel()
    K_data: list[float] = []
    K_row: list[int] = []
    K_col: list[int] = []
    printed = 0

    print("Step-by-step triplet append trace")
    print("  local a,b -> dofs[a], dofs[b] -> append row, col, value")
    print(
        "  step  a  b  local entry       global row  global col        "
        "K_data value"
    )

    for a in range(48):
        for b in range(48):
            K_row.append(int(dofs[a]))
            K_col.append(int(dofs[b]))
            K_data.append(float(Ke[a, b]))

            if printed < max_entries:
                print(
                    f"  {len(K_data):>4} {a:>2} {b:>2}  "
                    f"Ke[{a:>2},{b:>2}]  ->  "
                    f"{K_row[-1]:>10}  {K_col[-1]:>10}  {K_data[-1]:>16.6e}"
                )
                printed += 1

    hidden = len(K_data) - printed
    if hidden:
        print(f"  ... {hidden} additional appends omitted for readability")

    print()
    print("After this element")
    print(f"  len(K_row)  = {len(K_row)}")
    print(f"  len(K_col)  = {len(K_col)}")
    print(f"  len(K_data) = {len(K_data)}")
    print("  48 local rows x 48 local columns = 2304 triplet contributions")
    print()
    return K_data, K_row, K_col


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Print a small, step-by-step trace of global K assembly."
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        type=int,
        default=None,
        help="Element indices to trace. Defaults to --element only.",
    )
    parser.add_argument(
        "--element",
        type=int,
        default=0,
        help="Single element index to trace when --elements is not supplied.",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=18,
        help="Maximum triplet append lines to print per element.",
    )
    parser.add_argument(
        "--all-entries",
        action="store_true",
        help="Print all 48 x 48 = 2304 local-to-global triplet appends.",
    )
    parser.add_argument(
        "--skip-ke",
        action="store_true",
        help="Use a zero Ke placeholder so you can focus only on the index mapping.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    beam = build_beam()
    element_indices = args.elements if args.elements is not None else [args.element]
    max_entries = 48 * 48 if args.all_entries else args.max_entries

    print_global_context(beam)
    explain_triplet_lists()

    for elem_idx in element_indices:
        if elem_idx < 0 or elem_idx >= beam.get_num_elements():
            raise ValueError(
                f"Element {elem_idx} is outside 0..{beam.get_num_elements() - 1}"
            )

        print_element_identity(beam, elem_idx)
        Ke = np.zeros((48, 48)) if args.skip_ke else build_element_stiffness(beam, elem_idx)
        print_ke_summary(Ke)
        trace_triplet_appends(beam, elem_idx, Ke, max_entries)

    print("Key idea")
    print("  Local Ke entries are not inserted at local addresses 0..47.")
    print("  The dofs array translates each local row/column into global K addresses.")
    print("  COO triplets collect those addresses first; K.tocsr() later compresses them")
    print("  into the sparse format used by the solver.")


if __name__ == "__main__":
    main()
