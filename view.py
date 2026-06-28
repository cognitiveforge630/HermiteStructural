import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

from checkme import E, I, L, P, cantilever_point_load_deflection
from FEMHermiteBeamRegion import FEMHermiteBeamRegion


DEFAULT_RESULTS_PATH = Path("U.npy")
LOAD_DOF = 1  # uy: same transverse point-load direction used by solver.py


def build_beam():
    return FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=120.0,
        nx=5,
        ny=9,
        nz=121,
        E=2000000,
        nu=0.3,
        gamma=0.0174,
        w=-0.4,
    )


def load_displacement_vector(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find {path}. Run solver.py first to create the displacement file."
        )

    try:
        return np.load(path)
    except ValueError:
        # Older versions of this project wrote text data to a .npy file.
        return np.loadtxt(path, delimiter=",")


def build_unstructured_grid(beam):
    points = beam.get_all_points()
    num_elems = beam.get_num_elements()
    cell_connectivity = np.hstack(
        [np.full((num_elems, 1), 8, dtype=int), beam.elements]
    )
    cells = cell_connectivity.ravel()
    cell_types = np.full(num_elems, pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, points)


def validate_displacement_vector(beam, U):
    expected_size = beam.get_all_points().shape[0] * beam.ndof_per_node
    if U.size != expected_size:
        raise ValueError(
            f"Expected {expected_size} displacement values, but {U.size} were loaded."
        )


def add_pinned_node_cubes(plotter, beam, cube_size=None):
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


def get_displacement_summary(beam, U_reshaped):
    end_nodes = beam.get_nodes_at_max_z()
    end_uy = U_reshaped[end_nodes, LOAD_DOF]
    actual_free_end_uy = end_uy[np.argmax(np.abs(end_uy))]
    theoretical_free_end_uy = cantilever_point_load_deflection(P, L, E, I)
    percent_difference = (
        abs(theoretical_free_end_uy - actual_free_end_uy)
        / abs(theoretical_free_end_uy)
        * 100
    )

    return theoretical_free_end_uy, actual_free_end_uy, percent_difference


def show_view(results_path=DEFAULT_RESULTS_PATH, warp_factor=1.0):
    beam = build_beam()
    U = load_displacement_vector(results_path)
    validate_displacement_vector(beam, U)

    grid = build_unstructured_grid(beam)
    points = beam.get_all_points()
    n_nodes = points.shape[0]

    U_reshaped = U.reshape(n_nodes, beam.ndof_per_node)
    displacements = U_reshaped[:, :3]
    theoretical_uy, actual_uy, percent_difference = get_displacement_summary(
        beam,
        U_reshaped,
    )

    grid.point_data["displacements"] = displacements
    warped_grid = grid.warp_by_vector("displacements", factor=warp_factor)

    disp_magnitude = np.linalg.norm(displacements, axis=1)
    warped_grid.point_data["disp_magnitude"] = disp_magnitude

    plotter = pv.Plotter()
    plotter.add_mesh(
        warped_grid,
        scalars="disp_magnitude",
        cmap="coolwarm",
        scalar_bar_args={
            "title": "Displacement Magnitude",
            "title_font_size": 22,
            "label_font_size": 10,
            "n_labels": 5,
            "fmt": "%.3g",
            "height": 0.08,
            "width": 0.55,
            "position_x": 0.32,
            "position_y": 0.05,
            "vertical": False,
        },
        show_edges=True,
        edge_color="#263238",
        line_width=0.6,
    )
    add_pinned_node_cubes(plotter, beam)
    plotter.add_text(
        "Solved beam displacement with pinned nodes at z = 0",
        position="upper_left",
        font_size=3,
        color="#111111",
    )
    plotter.add_text(
        "\n".join(
            [
                "Cantilever point-load displacement",
                f"Theoretical uy: {theoretical_uy:.3e} in",
                f"Actual FEA uy:  {actual_uy:.3e} in",
                f"Difference:     {percent_difference:.2f}%",
            ]
        ),
        position="upper_right",
        font_size=3,
        color="#111111",
    )
    plotter.show_bounds(
        bounds=warped_grid.bounds,
        grid="front",
        location="outer",
        xtitle="x",
        ytitle="y",
        ztitle="z",
        font_size=8,
        bold=False,
        n_xlabels=3,
        n_ylabels=3,
        n_zlabels=5,
        fmt="%.0f",
    )
    plotter.add_axes()
    plotter.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description="View saved Hermite beam results from solver.py."
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS_PATH,
        help="Path to the displacement vector written by solver.py.",
    )
    parser.add_argument(
        "--warp-factor",
        type=float,
        default=1.0,
        help="Scale factor applied to displacement vectors in the warped view.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    show_view(args.results, args.warp_factor)
