import argparse
from pathlib import Path

import numpy as np

from checkme import E, I, L, P, cantilever_point_load_deflection
from FEMHermiteBeamRegion import FEMHermiteBeamRegion


DEFAULT_RESULTS_PATH = Path("U.npy")
LOAD_DOF = 1  # uy: same transverse point-load direction used by solver.py


def build_beam():
    # Construct with the identical mesh and orthotropic axes as the solve.
    from solver import build_beam as build_solver_beam
    return build_solver_beam()


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
    import pyvista as pv
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
    import pyvista as pv
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


def show_view(results_path=DEFAULT_RESULTS_PATH, warp_factor=1.0, *, hud=True, hud_font_size=10, title_font_size=10, scalar_title_font_size=10, scalar_label_font_size=10, bounds_font_size=10, screenshot=None, off_screen=False, window_width=1600, window_height=900):
    beam = build_beam()
    U = load_displacement_vector(results_path)
    validate_displacement_vector(beam, U)

    import pyvista as pv
    grid = build_unstructured_grid(beam)
    import json
    result_file = Path('hermite_result.json')
    if not result_file.exists():
        raise FileNotFoundError('Missing hermite_result.json: run python .\\run_legacy_hermite.py before viewing results.')
    result = json.loads(result_file.read_text())
    points = beam.get_all_points()
    n_nodes = points.shape[0]

    U_reshaped = U.reshape(n_nodes, beam.ndof_per_node)
    displacements = U_reshaped[:, :3]
    theoretical_uy, actual_uy, percent_difference = get_displacement_summary(beam, U_reshaped)
    if result:
        theoretical_uy = float(result['awc_point_plus_gravity_uy_in'])
        actual_uy = float(result['fem_tip_uy_in'])
        percent_difference = abs((actual_uy-theoretical_uy)/theoretical_uy)*100

    grid.point_data["displacements"] = displacements
    warped_grid = grid.warp_by_vector("displacements", factor=warp_factor)

    disp_magnitude = np.linalg.norm(displacements, axis=1)
    warped_grid.point_data["disp_magnitude"] = disp_magnitude

    plotter = pv.Plotter(off_screen=off_screen, window_size=(window_width, window_height))
    plotter.add_mesh(
        warped_grid,
        scalars="disp_magnitude",
        cmap="coolwarm",
        scalar_bar_args={
            "title": "Displacement Magnitude",
            "title_font_size": scalar_title_font_size,
            "label_font_size": scalar_label_font_size,
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
    if hud:
        plotter.add_text(
            "HermiteStructural | Historical Hermite Hex8 | 6 DOF/node",
            position="upper_left", font_size=title_font_size, color="#111111",
            name="hermite_title",
        )
        plotter.add_text(
            "\n".join([
                "Orthotropic Redwood | Cantilever",
                f"Gravity specific weight: {result['weight_density_lbf_in3']:g} lbf/in^3" if result else "Gravity: result metadata unavailable",
                f"End-face traction resultant: {result['end_traction_y_lbf']:+g} lbf" if result else "End-face traction: Gauss integrated",
                f"End-face force: {result['end_traction_y_lbf']:+.3f} lbf" if result else "End load: unavailable",
                f"Self-weight: {result['self_weight_y_lbf']:+.3f} lbf" if result else "Self-weight: unavailable",
                f"AWC end-load uy: {result['awc_components']['end_load_uy_in']:+.7f} in" if result else "AWC end load: unavailable",
                f"AWC gravity uy: {result['awc_components']['self_weight_uy_in']:+.7f} in" if result else "AWC self-weight: unavailable",
                f"AWC total uy: {theoretical_uy:+.7f} in",
                f"Hermite FEM tip uy:   {actual_uy:+.7f} in",
                f"Absolute difference:  {percent_difference:.3f}%",
                f"Warp factor: {warp_factor:g}x",
                "AWC reference includes gravity; FEM modulus fixed",
                "Research demonstration — not a certified error bound",
            ]),
            position="upper_right", font_size=hud_font_size, color="#111111",
            name="hermite_hud",
        )
    plotter.show_bounds(
        bounds=warped_grid.bounds,
        grid="front",
        location="outer",
        xtitle="x",
        ytitle="y",
        ztitle="z",
        font_size=bounds_font_size,
        bold=False,
        n_xlabels=3,
        n_ylabels=3,
        n_zlabels=5,
        fmt="%.0f",
    )
    plotter.add_axes()
    plotter.show(screenshot=str(screenshot) if screenshot else None)


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
    parser.add_argument("--hud", action=argparse.BooleanOptionalAction, default=True,
                        help="Show/hide the Hermite vs AWC comparison overlay.")
    parser.add_argument("--hud-font-size", type=int, default=10)
    parser.add_argument("--title-font-size", type=int, default=10)
    parser.add_argument("--scalar-title-font-size", type=int, default=10)
    parser.add_argument("--scalar-label-font-size", type=int, default=10)
    parser.add_argument("--bounds-font-size", type=int, default=10)
    parser.add_argument("--window-width", type=int, default=1600)
    parser.add_argument("--window-height", type=int, default=900)
    parser.add_argument("--screenshot", type=Path, default=None)
    parser.add_argument("--off-screen", action="store_true",
                        help="Render screenshot without opening a window (requires working VTK backend).")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.screenshot:
        args.screenshot.parent.mkdir(parents=True, exist_ok=True)
    for value in (args.hud_font_size, args.title_font_size, args.scalar_title_font_size,
                  args.scalar_label_font_size, args.bounds_font_size,
                  args.window_width, args.window_height):
        if value <= 0:
            raise SystemExit("Font sizes and window dimensions must be positive.")
    show_view(args.results, args.warp_factor, hud=args.hud,
              hud_font_size=args.hud_font_size,
              title_font_size=args.title_font_size,
              scalar_title_font_size=args.scalar_title_font_size,
              scalar_label_font_size=args.scalar_label_font_size,
              bounds_font_size=args.bounds_font_size,
              screenshot=args.screenshot, off_screen=args.off_screen,
              window_width=args.window_width, window_height=args.window_height)
