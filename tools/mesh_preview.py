"""Preview supported structured beam meshes without assembling or solving.

Use --summary-only on systems without PyVista. Vertex counts (--nx etc.)
are not physical element counts: divisions = counts - 1.
"""
from __future__ import annotations

import argparse
import math

import numpy as np

from core.reference.beam_theory import B, H, L

FAMILIES = {"hex8": 1, "hex27": 2, "hex64": 3, "hermite": 1}


def mesh_information(family: str, vertex_counts: tuple[int, int, int]) -> dict:
    """Pure mesh geometry and counts, with no visualization dependency."""
    if family not in FAMILIES:
        raise ValueError(f"Unknown mesh family: {family}")
    if len(vertex_counts) != 3 or min(vertex_counts) < 2:
        raise ValueError("nx, ny, nz must each be >= 2")
    order = FAMILIES[family]
    divisions = tuple(v - 1 for v in vertex_counts)
    node_shape = tuple(order * d + 1 for d in divisions)
    return {
        "family": family,
        "polynomial_order": order if family != "hermite" else None,
        "vertex_counts": vertex_counts,
        "physical_element_divisions": divisions,
        "physical_elements": math.prod(divisions),
        "geometric_nodes_per_element": (order + 1) ** 3,
        "global_geometric_nodes": math.prod(node_shape),
        "node_shape": node_shape,
        "note": (
            "Hermite: this is the shared Hex8 geometric mesh; additional enriched "
            "unknowns are formulation-dependent, not represented by these points."
            if family == "hermite" else
            "Hex8 geometry is previewable; polynomial-order solver supports only Hex27/Hex64."
            if family == "hex8" else
            "C0 tensor-product interpolation nodes; mesh preview does not perform a solve."
        ),
    }


def coordinates(shape: tuple[int, int, int]) -> np.ndarray:
    """Return xyz array of a regular beam lattice."""
    axes = (
        np.linspace(-B / 2, B / 2, shape[0]),
        np.linspace(-H / 2, H / 2, shape[1]),
        np.linspace(0.0, L, shape[2]),
    )
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)


def integration_locations(order: int, divisions: tuple[int, int, int], q: int) -> np.ndarray:
    """Physical tensor-Gauss positions for every affine cell (visualization only)."""
    from numpy.polynomial.legendre import leggauss

    xi, _ = leggauss(q)
    extent = np.array((B, H, L), dtype=float) / divisions
    low = np.array((-B / 2, -H / 2, 0.0))
    local = np.stack(np.meshgrid(xi, xi, xi, indexing="ij"), axis=-1).reshape(-1, 3)
    indices = np.stack(np.meshgrid(*(np.arange(n) for n in divisions), indexing="ij"), axis=-1).reshape(-1, 3)
    return (low + (indices[:, None, :] + (local[None, :, :] + 1) / 2) * extent).reshape(-1, 3)


def render(info: dict, *, nodes: bool, gauss_points: bool, quadrature_order: int | None) -> None:
    try:
        import pyvista as pv
    except ImportError as exc:
        raise SystemExit("Interactive previews require PyVista: python -m pip install pyvista") from exc

    vertex = info["vertex_counts"]
    coarse_points = coordinates(vertex)
    fine_points = coordinates(info["node_shape"])

    def grid(shape, positions):
        # StructuredGrid accepts coordinate arrays with the same ordering as meshgrid.
        axes = (
            np.linspace(-B / 2, B / 2, shape[0]),
            np.linspace(-H / 2, H / 2, shape[1]),
            np.linspace(0, L, shape[2]),
        )
        return pv.StructuredGrid(*np.meshgrid(*axes, indexing="ij"))

    coarse = grid(vertex, coarse_points)
    plotter = pv.Plotter()
    plotter.add_mesh(coarse, color="lightsteelblue", opacity=0.22, show_edges=False)
    plotter.add_mesh(coarse.extract_all_edges(), color="black", line_width=1.5)
    if nodes:
        plotter.add_points(fine_points, color="crimson", point_size=5,
                           render_points_as_spheres=True)
    if gauss_points:
        q = quadrature_order or (info["polynomial_order"] or 1) + 1
        points = integration_locations(info["polynomial_order"] or 1,
                                       info["physical_element_divisions"], q)
        plotter.add_points(points, color="darkorange", point_size=4,
                           render_points_as_spheres=True)
        print(f"Gaussian quadrature: {q} per axis, {len(points):,} plotted positions")
    plotter.add_axes()
    plotter.show_grid()
    plotter.add_text(f"{info['family']}: {info['physical_elements']:,} elements; "
                     f"{info['global_geometric_nodes']:,} geometric nodes", font_size=10)
    plotter.show()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=tuple(FAMILIES), default="hex64")
    parser.add_argument("--nx", type=int, default=3, help="Original vertex count in X")
    parser.add_argument("--ny", type=int, default=5, help="Original vertex count in Y")
    parser.add_argument("--nz", type=int, default=25, help="Original vertex count in Z")
    parser.add_argument("--nodes", action="store_true", help="Show geometric/interpolation nodes in red")
    parser.add_argument("--gauss-points", action="store_true", help="Show quadrature sample positions in orange")
    parser.add_argument("--quadrature-order", type=int, help="Gaussian points per axis; default p+1")
    parser.add_argument("--summary-only", action="store_true", help="Print mesh counts without opening PyVista")
    args = parser.parse_args(argv)
    if min(args.nx, args.ny, args.nz) < 2:
        parser.error("nx, ny, nz must be at least 2")
    if args.quadrature_order is not None and args.quadrature_order < 1:
        parser.error("--quadrature-order must be >= 1")
    info = mesh_information(args.family, (args.nx, args.ny, args.nz))
    for key, value in info.items():
        print(f"{key}: {value}")
    if not args.summary_only:
        render(info, nodes=args.nodes, gauss_points=args.gauss_points,
               quadrature_order=args.quadrature_order)


if __name__ == "__main__":
    main()
