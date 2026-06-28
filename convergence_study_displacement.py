"""Mesh convergence study for the Hermite cantilever free-end displacement.

This script answers the practical calibration question:

    Is the displacement at the current solver mesh already converged enough
    that the calibrated/slackened Ez is meaningful, rather than a mesh fluke?

Important: the study does *not* recalibrate E independently for each mesh and
then compare those already-forced answers.  Instead it keeps the same reference
material/load/boundary path for every mesh, solves the displacement, and reports
both:

* the free-end uy displacement at the fixed reference Ez, and
* the Ez that would be implied by that mesh if one calibrated to the same
  reference beam-theory displacement target.

If both the displacement and the implied Ez stabilize as the mesh approaches the
current ``solver.py`` mesh, the calibration is anchored to a converged mesh
rather than being an arbitrary fit at a coarse density.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.sparse.linalg import spsolve

from solver import (
    LOAD_DOF,
    REFERENCE_LONGITUDINAL_E,
    REFERENCE_TARGET_UY,
    build_beam,
    build_end_face_traction_load_vector,
    calibrated_longitudinal_E,
    free_end_uy,
)

# Current solver.py mesh is nodes 5x9x121, i.e. elements 4x8x120.
# These default levels keep approximately equal element dimensions through the
# 4 in x 8 in x 120 in beam: h = 4, 2, 4/3, and 1 inch.
DEFAULT_ELEMENT_MESHES: tuple[tuple[int, int, int], ...] = (
    (1, 2, 30),
    (2, 4, 60),
    (3, 6, 90),
    (4, 8, 120),
)


@dataclass
class ConvergenceRow:
    level: int
    ex: int
    ey: int
    ez: int
    nx: int
    ny: int
    nz: int
    nodes: int
    elements: int
    dofs: int
    hx_in: float
    hy_in: float
    hz_in: float
    seed_E_psi: float
    free_end_uy_in: float
    target_uy_in: float
    abs_error_vs_target_in: float
    pct_error_vs_target: float
    pct_change_uy_from_previous: float | None
    implied_calibrated_Ez_psi: float
    pct_change_implied_Ez_from_previous: float | None
    relative_equilibrium_residual: float


def parse_element_mesh(text: str) -> tuple[int, int, int]:
    """Parse a mesh token like ``4x8x120`` as element counts."""
    parts = text.lower().replace("×", "x").split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"mesh '{text}' must look like exXeyXez, for example 4x8x120"
        )
    try:
        ex, ey, ez = (int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"mesh '{text}' must contain integer element counts"
        ) from exc
    if ex < 1 or ey < 1 or ez < 1:
        raise argparse.ArgumentTypeError(
            f"mesh '{text}' must use positive element counts"
        )
    return ex, ey, ez


def fixed_and_free_dofs(beam) -> tuple[list[int], list[int]]:
    fixed_dof: list[int] = []
    for node in beam.get_nodes_at_min_z():
        fixed_dof.extend(beam.global_dofs[node])
    fixed = set(fixed_dof)
    free_dof = [dof for dof in range(beam.grid_points.shape[0] * beam.ndof_per_node) if dof not in fixed]
    return fixed_dof, free_dof


def solve_with_equilibrium_residual(beam, force: np.ndarray) -> tuple[np.ndarray, float]:
    """Solve KU=F and return U plus a normalized free-DOF residual."""
    stiffness = beam.build_global_K()
    fixed_dof, _ = fixed_and_free_dofs(beam)
    free_mask = np.ones(stiffness.shape[0], dtype=bool)
    free_mask[fixed_dof] = False
    K_ff = stiffness[free_mask][:, free_mask]
    F_f = force[free_mask]
    U_f = spsolve(K_ff, F_f)

    U = np.zeros(stiffness.shape[0])
    U[free_mask] = U_f

    residual = K_ff @ U_f - F_f
    denom = max(np.linalg.norm(F_f), 1.0)
    relative_residual = float(np.linalg.norm(residual) / denom)
    return U, relative_residual


def pct_change(current: float, previous: float | None) -> float | None:
    if previous is None:
        return None
    denom = max(abs(current), np.finfo(float).eps)
    return abs(current - previous) / denom * 100.0


def run_one_mesh(level: int, element_mesh: tuple[int, int, int], previous: ConvergenceRow | None) -> ConvergenceRow:
    ex, ey, ez = element_mesh
    nx, ny, nz = ex + 1, ey + 1, ez + 1
    beam = build_beam(REFERENCE_LONGITUDINAL_E, nx=nx, ny=ny, nz=nz)
    force, _, _ = build_end_face_traction_load_vector(beam)
    U, rel_resid = solve_with_equilibrium_residual(beam, force)

    uy = float(free_end_uy(beam, U))
    implied_E = float(calibrated_longitudinal_E(REFERENCE_LONGITUDINAL_E, uy, REFERENCE_TARGET_UY))
    abs_error = abs(uy - REFERENCE_TARGET_UY)
    pct_error = abs_error / abs(REFERENCE_TARGET_UY) * 100.0

    prev_uy = None if previous is None else previous.free_end_uy_in
    prev_E = None if previous is None else previous.implied_calibrated_Ez_psi

    return ConvergenceRow(
        level=level,
        ex=ex,
        ey=ey,
        ez=ez,
        nx=nx,
        ny=ny,
        nz=nz,
        nodes=int(beam.grid_points.shape[0]),
        elements=int(beam.get_num_elements()),
        dofs=int(beam.grid_points.shape[0] * beam.ndof_per_node),
        hx_in=beam.Lx / ex,
        hy_in=beam.Ly / ey,
        hz_in=beam.Lz / ez,
        seed_E_psi=float(REFERENCE_LONGITUDINAL_E),
        free_end_uy_in=uy,
        target_uy_in=float(REFERENCE_TARGET_UY),
        abs_error_vs_target_in=float(abs_error),
        pct_error_vs_target=float(pct_error),
        pct_change_uy_from_previous=pct_change(uy, prev_uy),
        implied_calibrated_Ez_psi=implied_E,
        pct_change_implied_Ez_from_previous=pct_change(implied_E, prev_E),
        relative_equilibrium_residual=rel_resid,
    )


def format_optional(value: float | None, fmt: str = ".3f") -> str:
    if value is None:
        return "first"
    return format(value, fmt)


def print_table(rows: list[ConvergenceRow]) -> None:
    print("Hermite cantilever displacement convergence study")
    print("  load model              : end-face traction integrated over z=L")
    print("  boundary condition      : fixed-free cantilever; all DOFs fixed at z=0")
    print("  seed/reference Ez       : {:.6f} psi".format(REFERENCE_LONGITUDINAL_E))
    print("  target free-end uy      : {:.12e} in".format(REFERENCE_TARGET_UY))
    print("  current solver.py mesh  : 4x8x120 elements = 5x9x121 nodes")
    print()
    header = (
        "lvl  elems         nodes    dofs    h(in)       free-end uy(in)  "
        "chg uy %   implied Ez(psi)  chg Ez %   eq rel"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        h_text = f"{row.hx_in:.3g}/{row.hy_in:.3g}/{row.hz_in:.3g}"
        mesh_text = f"{row.ex}x{row.ey}x{row.ez}"
        print(
            f"{row.level:>3d}  "
            f"{mesh_text:<10s}  "
            f"{row.nodes:>7d}  "
            f"{row.dofs:>7d}  "
            f"{h_text:<10s}  "
            f"{row.free_end_uy_in:>16.9e}  "
            f"{format_optional(row.pct_change_uy_from_previous, '.3f'):>8s}  "
            f"{row.implied_calibrated_Ez_psi:>15.6f}  "
            f"{format_optional(row.pct_change_implied_Ez_from_previous, '.3f'):>8s}  "
            f"{row.relative_equilibrium_residual:>7.1e}"
        )
    print()

    if len(rows) >= 2:
        last = rows[-1]
        prev = rows[-2]
        print("Current-mesh convergence anchor:")
        print(
            "  Last two mesh levels change in free-end uy: "
            f"{last.pct_change_uy_from_previous:.6f}% "
            f"({prev.ex}x{prev.ey}x{prev.ez} -> {last.ex}x{last.ey}x{last.ez} elements)"
        )
        print(
            "  Last two mesh levels change in implied calibrated Ez: "
            f"{last.pct_change_implied_Ez_from_previous:.6f}%"
        )
        print(
            "  Calibrated Ez implied by current solver.py mesh: "
            f"{last.implied_calibrated_Ez_psi:.6f} psi"
        )
        print()


def write_outputs(rows: list[ConvergenceRow], output_prefix: Path) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    dict_rows = [asdict(row) for row in rows]

    csv_path = output_prefix.with_suffix(".csv")
    json_path = output_prefix.with_suffix(".json")
    txt_path = output_prefix.with_suffix(".txt")

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(dict_rows[0].keys()))
        writer.writeheader()
        writer.writerows(dict_rows)

    with json_path.open("w") as f:
        json.dump(
            {
                "description": "Hermite cantilever displacement convergence study",
                "seed_reference_E_psi": REFERENCE_LONGITUDINAL_E,
                "target_uy_in": REFERENCE_TARGET_UY,
                "rows": dict_rows,
            },
            f,
            indent=2,
        )

    # Also write a human-readable table to txt.
    from io import StringIO
    import contextlib

    buffer = StringIO()
    with contextlib.redirect_stdout(buffer):
        print_table(rows)
    txt_path.write_text(buffer.getvalue())

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {txt_path}")


def mesh_list_from_args(args: argparse.Namespace) -> list[tuple[int, int, int]]:
    meshes = list(args.element_meshes or DEFAULT_ELEMENT_MESHES)
    if args.quick:
        meshes = meshes[:2]
    return meshes


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--element-meshes",
        nargs="+",
        type=parse_element_mesh,
        metavar="EXxEYxEZ",
        help=(
            "Element-count mesh levels to run. Defaults to "
            "1x2x30 2x4x60 3x6x90 4x8x120, where 4x8x120 is solver.py."
        ),
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run only the first two mesh levels as a fast smoke test.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("outputs/displacement_convergence"),
        help="Path prefix for .csv, .json, and .txt outputs.",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Print the table only; do not write output files.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    rows: list[ConvergenceRow] = []
    for level, mesh in enumerate(mesh_list_from_args(args), start=1):
        print(f"Solving mesh level {level}: {mesh[0]}x{mesh[1]}x{mesh[2]} elements", flush=True)
        previous = rows[-1] if rows else None
        rows.append(run_one_mesh(level, mesh, previous))

    print()
    print_table(rows)
    if not args.no_write:
        write_outputs(rows, args.output_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
