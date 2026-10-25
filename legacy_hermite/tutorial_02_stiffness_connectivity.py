"""Tutorial 02: show how element connectivity becomes global K connectivity.

The solver eventually computes numeric stiffness values, but the first idea is
simpler: an element connects its 8 nodes, and each node owns 6 degrees of
freedom. One hexahedral element therefore contributes a 48 x 48 block into the
global stiffness matrix.

Run:
    py tutorial_02_stiffness_connectivity.py

Optional:
    py tutorial_02_stiffness_connectivity.py --screenshot outputs/k_connectivity.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


DOF_LABELS = ("ux", "uy", "uz", "phix", "phiy", "phiz")

NODE_SPHERE_RADIUS = 0.1125
ELEMENT_NODE_RADIUS_SCALE = 1.08
SELECTED_NODE_RADIUS_SCALE = 1.45
CONNECTION_LINE_WIDTH = 3

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


def build_beam(use_full_mesh: bool = True) -> FEMHermiteBeamRegion:
    params = BEAM_PARAMS if use_full_mesh else MINI_BEAM_PARAMS
    return FEMHermiteBeamRegion(**params)


def build_unstructured_grid(beam: FEMHermiteBeamRegion) -> pv.UnstructuredGrid:
    num_elems = beam.get_num_elements()
    cell_connectivity = np.hstack(
        [np.full((num_elems, 1), 8, dtype=int), beam.elements]
    )
    cells = cell_connectivity.ravel()
    cell_types = np.full(num_elems, pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, beam.get_all_points())


def describe_global_counts(beam: FEMHermiteBeamRegion) -> None:
    n_nodes = beam.get_all_points().shape[0]
    n_dofs = n_nodes * beam.ndof_per_node

    print("Global mesh")
    print(f"  nodes: {n_nodes}")
    print(f"  elements: {beam.get_num_elements()}")
    print(f"  dofs per node: {beam.ndof_per_node} ({', '.join(DOF_LABELS)})")
    print(f"  global K size: {n_dofs} x {n_dofs}")
    print()


def describe_element_mapping(beam: FEMHermiteBeamRegion, elem_idx: int) -> None:
    node_indices = beam.elements[elem_idx]
    coords = beam.get_element_nodes(elem_idx)
    dofs = beam.get_element_global_dofs(elem_idx)

    print(f"Element {elem_idx} contributes one local 48 x 48 stiffness block")
    print("  8 element nodes x 6 dofs per node = 48 local dofs")
    print()
    print("  local_node  global_node        x        y        z        global_dofs")

    for local_idx, (global_node, xyz, global_dofs) in enumerate(
        zip(node_indices, coords, dofs)
    ):
        dof_text = ", ".join(str(int(dof)) for dof in global_dofs)
        print(
            f"  {local_idx:>10}  {global_node:>11}  "
            f"{xyz[0]:>7.2f}  {xyz[1]:>7.2f}  {xyz[2]:>7.2f}    {dof_text}"
        )

    print()
    print("How this enters K")
    print("  local Ke[row, col] is added to global K[global_row, global_col]")
    print("  every local row/column maps through the global_dofs listed above")
    print()


def get_nodes_connected_by_element(beam: FEMHermiteBeamRegion, elem_idx: int) -> set[int]:
    return set(int(node) for node in beam.elements[elem_idx])


def get_elements_touching_node(beam: FEMHermiteBeamRegion, node_idx: int) -> list[int]:
    return [
        elem_idx
        for elem_idx, element in enumerate(beam.elements)
        if node_idx in element
    ]


def describe_node_connections(beam: FEMHermiteBeamRegion, node_idx: int) -> None:
    touching_elements = get_elements_touching_node(beam, node_idx)
    connected_nodes: set[int] = set()

    for elem_idx in touching_elements:
        connected_nodes.update(get_nodes_connected_by_element(beam, elem_idx))

    connected_nodes.discard(node_idx)

    print(f"Node {node_idx} appears in elements: {touching_elements}")
    print(
        f"Node {node_idx} can couple to {len(connected_nodes)} neighboring nodes "
        "through those element stiffness blocks."
    )
    print(f"Connected nodes: {sorted(connected_nodes)}")
    print()


def write_local_element_pattern_svg(path: Path) -> None:
    cell = 9
    margin = 36
    size = 48 * cell
    width = size + 2 * margin
    height = size + 2 * margin
    colors = ("#2f80ed", "#56cc9d")

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin}" y="24" font-family="Arial" font-size="16">'
        "One element: local 48 x 48 K block</text>",
    ]

    for row_node in range(8):
        for col_node in range(8):
            color = colors[(row_node + col_node) % 2]
            x = margin + col_node * 6 * cell
            y = margin + row_node * 6 * cell
            parts.append(
                f'<rect x="{x}" y="{y}" width="{6 * cell}" height="{6 * cell}" '
                f'fill="{color}" opacity="0.72"/>'
            )

    for i in range(9):
        pos = margin + i * 6 * cell
        stroke = "#111111" if i in (0, 8) else "#ffffff"
        parts.append(
            f'<line x1="{margin}" y1="{pos}" x2="{margin + size}" y2="{pos}" '
            f'stroke="{stroke}" stroke-width="1"/>'
        )
        parts.append(
            f'<line x1="{pos}" y1="{margin}" x2="{pos}" y2="{margin + size}" '
            f'stroke="{stroke}" stroke-width="1"/>'
        )

    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def write_global_pattern_svg(path: Path) -> None:
    beam = build_beam(use_full_mesh=False)
    n_dofs = beam.get_all_points().shape[0] * beam.ndof_per_node
    cell = 3
    margin = 42
    size = n_dofs * cell
    width = size + 2 * margin
    height = size + 2 * margin

    occupied: set[tuple[int, int]] = set()
    for elem_idx in range(beam.get_num_elements()):
        dofs = beam.get_element_global_dofs(elem_idx).ravel()
        for row in dofs:
            for col in dofs:
                occupied.add((int(row), int(col)))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin}" y="25" font-family="Arial" font-size="16">'
        "Small mesh: global K nonzero pattern from element connectivity</text>",
        f'<rect x="{margin}" y="{margin}" width="{size}" height="{size}" '
        'fill="#f4f6f8" stroke="#111111" stroke-width="1"/>',
    ]

    for row, col in sorted(occupied):
        x = margin + col * cell
        y = margin + row * cell
        parts.append(
            f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="#263238"/>'
        )

    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def add_node_sphere(
    plotter: pv.Plotter,
    point: np.ndarray,
    color: str,
    radius: float,
    opacity: float = 0.95,
) -> None:
    sphere = pv.Sphere(
        radius=radius,
        center=point,
        theta_resolution=24,
        phi_resolution=24,
    )
    plotter.add_mesh(sphere, color=color, opacity=opacity)


def add_connection_line(
    plotter: pv.Plotter,
    start: np.ndarray,
    end: np.ndarray,
    color: str = "#f2994a",
) -> None:
    line = pv.Line(start, end)
    plotter.add_mesh(line, color=color, line_width=CONNECTION_LINE_WIDTH, opacity=0.65)


def make_plot(
    elem_idx: int,
    node_idx: int,
    screenshot: Path | None,
    sphere_radius: float,
) -> None:
    beam = build_beam(use_full_mesh=True)
    grid = build_unstructured_grid(beam)
    element_nodes = get_nodes_connected_by_element(beam, elem_idx)
    touching_elements = get_elements_touching_node(beam, node_idx)
    connected_nodes = set()

    for touching_elem in touching_elements:
        connected_nodes.update(get_nodes_connected_by_element(beam, touching_elem))
    connected_nodes.discard(node_idx)

    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(1600, 1000))
    plotter.add_mesh(
        grid,
        style="wireframe",
        color="#263238",
        line_width=1,
        opacity=0.22,
    )
    for touching_elem in touching_elements:
        plotter.add_mesh(
            grid.extract_cells(touching_elem),
            color="#f2994a",
            opacity=0.12,
            show_edges=True,
            edge_color="#f2994a",
            line_width=2,
        )

    plotter.add_mesh(
        grid.extract_cells(elem_idx),
        color="#56cc9d",
        opacity=0.18,
        show_edges=True,
        edge_color="#2f80ed",
        line_width=3,
    )

    points = beam.get_all_points()
    selected_point = points[node_idx]

    for connected_node in connected_nodes:
        add_connection_line(plotter, selected_point, points[connected_node])
        add_node_sphere(plotter, points[connected_node], "#f2994a", sphere_radius)

    for element_node in element_nodes:
        if element_node != node_idx:
            add_node_sphere(
                plotter,
                points[element_node],
                "#2f80ed",
                sphere_radius * ELEMENT_NODE_RADIUS_SCALE,
            )

    add_node_sphere(
        plotter,
        selected_point,
        "#eb5757",
        sphere_radius * SELECTED_NODE_RADIUS_SCALE,
    )

    plotter.add_text(
        "K connectivity: shared nodes create shared global rows and columns",
        position="upper_left",
        font_size=14,
        color="#111111",
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
        description="Explain and visualize how element connectivity builds global K."
    )
    parser.add_argument("--element", type=int, default=0, help="Element to highlight.")
    parser.add_argument(
        "--node",
        type=int,
        default=1211,
        help="Node to trace through K. The default belongs to several elements.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Optional PNG path for the PyVista connectivity view.",
    )
    parser.add_argument(
        "--local-svg",
        type=Path,
        default=Path("outputs/local_element_k_pattern.svg"),
        help="SVG showing the 48 x 48 local element K block.",
    )
    parser.add_argument(
        "--global-svg",
        type=Path,
        default=Path("outputs/global_k_pattern_small_mesh.svg"),
        help="SVG showing a small global K sparsity pattern.",
    )
    parser.add_argument(
        "--sphere-radius",
        type=float,
        default=NODE_SPHERE_RADIUS,
        help="Radius for node spheres in the 3D connectivity view.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    beam = build_beam(use_full_mesh=True)

    describe_global_counts(beam)
    describe_element_mapping(beam, args.element)
    describe_node_connections(beam, args.node)

    write_local_element_pattern_svg(args.local_svg)
    write_global_pattern_svg(args.global_svg)
    print(f"Saved local element K pattern: {args.local_svg}")
    print(f"Saved small global K pattern: {args.global_svg}")

    make_plot(args.element, args.node, args.screenshot, args.sphere_radius)
