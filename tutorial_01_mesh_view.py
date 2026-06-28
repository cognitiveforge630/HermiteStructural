"""Tutorial 01: visualize the beam mesh and pinned boundary nodes.

This script is intentionally focused on geometry and boundary conditions.
It does not assemble the stiffness matrix or solve the finite element system.

Run:
    py tutorial_01_mesh_view.py

Optional:
    py tutorial_01_mesh_view.py --screenshot outputs/beam_mesh_pinned_nodes.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


BEAM_PARAMS = {
    "Lx": 4.0,
    "Ly": 8.0,
    "Lz": 120.0,
    "nx": 5,
    "ny": 9,
    "nz": 121,
    "E": 2_000_000.0,
    "nu": 0.3,
    "gamma": 0.0174,
    "w": -0.4,
}


def build_beam() -> FEMHermiteBeamRegion:
    """Create the same beam used by the solver examples."""
    return FEMHermiteBeamRegion(**BEAM_PARAMS)


def build_unstructured_grid(beam: FEMHermiteBeamRegion) -> pv.UnstructuredGrid:
    """Convert the beam's hex element connectivity into a PyVista grid."""
    n_elements = beam.get_num_elements()
    cell_connectivity = np.hstack(
        [np.full((n_elements, 1), 8, dtype=int), beam.elements]
    )
    cells = cell_connectivity.ravel()
    cell_types = np.full(n_elements, pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, beam.get_all_points())


def describe_mesh(beam: FEMHermiteBeamRegion) -> None:
    """Print the counts and spacings students should verify first."""
    n_nodes = beam.get_all_points().shape[0]
    n_elements = beam.get_num_elements()
    pinned_nodes = beam.get_nodes_at_min_z()

    dx = beam.Lx / (beam.nx - 1)
    dy = beam.Ly / (beam.ny - 1)
    dz = beam.Lz / (beam.nz - 1)

    print("Beam dimensions")
    print(f"  Lx x Ly x Lz: {beam.Lx:g} x {beam.Ly:g} x {beam.Lz:g}")
    print("Mesh resolution")
    print(f"  nodes in x/y/z: {beam.nx} x {beam.ny} x {beam.nz}")
    print(f"  element size dx/dy/dz: {dx:g}, {dy:g}, {dz:g}")
    print(f"  total nodes: {n_nodes}")
    print(f"  total hex elements: {n_elements}")
    print("Boundary condition")
    print(f"  pinned face: z = 0")
    print(f"  pinned nodes: {len(pinned_nodes)}")
    print(f"  fixed dofs: {len(pinned_nodes) * beam.ndof_per_node}")
    print(f"  dofs per node: {beam.ndof_per_node} (ux, uy, uz, phix, phiy, phiz)")


def add_pinned_node_points(
    plotter: pv.Plotter,
    beam: FEMHermiteBeamRegion,
) -> None:
    """Draw the actual fixed-node coordinates on the z=0 face."""
    fixed_points = beam.get_all_points()[beam.get_nodes_at_min_z()]

    print("Fixed node coordinates:")
    for point in fixed_points:
        print(point)

    plotter.add_points(
        fixed_points,
        color="#111111",
        point_size=18,
        render_points_as_spheres=True,
    )


def add_pinned_node_cubes(
    plotter: pv.Plotter,
    beam: FEMHermiteBeamRegion,
    cube_size: float | None = None,
) -> None:
    if cube_size is None:
        dx = beam.Lx / (beam.nx - 1)
        dy = beam.Ly / (beam.ny - 1)
        dz = beam.Lz / (beam.nz - 1)
        cube_size = 0.35 * min(dx, dy, dz)

    fixed_points = beam.get_all_points()[beam.get_nodes_at_min_z()]

    offset = np.array([0, 0, 0])

    for point in fixed_points:
        cube = pv.Cube(
            center=point + offset,
            x_length=cube_size,
            y_length=cube_size,
            z_length=cube_size,
        )
        plotter.add_mesh(
            cube,
            color="#e84a5f",
            opacity=0.35,
            show_edges=True,
            edge_color="#7f1d1d",
            line_width=1.0,
        )


    plotter.add_points(
        fixed_points,
        color="#111111",
        point_size=14,
        render_points_as_spheres=True,
    )

def make_plot(screenshot: Path | None = None) -> None:
    beam = build_beam()
    describe_mesh(beam)

    grid = build_unstructured_grid(beam)
    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(1600, 1000))

    plotter.add_mesh(
        grid,
        style="wireframe",
        color="#263238",
        line_width=1,
        opacity=0.35,
    )
    add_pinned_node_cubes(plotter, beam)


    plotter.add_text(
        "Hermite beam mesh: pinned nodes at z = 0",
        position="upper_left",
        font_size=8,
        color="#111111",
    )
    plotter.show_bounds(
        bounds=grid.bounds,
        grid="front",
        location="outer",
        xtitle="x",
        ytitle="y",
        ztitle="z",
    )
    plotter.add_axes()
    plotter.camera_position = [
        (95.0, -135.0, 95.0),
        (0.0, 0.0, 60.0),
        (0.0, 0.0, 1.0),
    ]

    if screenshot:
        screenshot.parent.mkdir(parents=True, exist_ok=True)
        plotter.show(screenshot=str(screenshot))
        print(f"Saved screenshot: {screenshot}")
    else:
        plotter.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Show the Hermite beam mesh and the pinned z=0 nodes."
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Optional PNG path for an off-screen screenshot.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    make_plot(args.screenshot)
