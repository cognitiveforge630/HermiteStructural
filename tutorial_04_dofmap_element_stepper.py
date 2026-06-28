"""Tutorial 04: step through element DOF maps in PyVista.

Each hexahedral element has 8 nodes and each node has 6 DOFs, so every element
has a 48-entry local-to-global DOF map. This visual stepper highlights one
element at a time and prints the corresponding node and DOF table.

Run:
    py tutorial_04_dofmap_element_stepper.py

Controls:
    Slider: choose an element index
    n: next element
    p: previous element

Optional:
    py tutorial_04_dofmap_element_stepper.py --full-mesh
    py tutorial_04_dofmap_element_stepper.py --element 5
    py tutorial_04_dofmap_element_stepper.py --screenshot outputs/dofmap_element_0.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


DOF_LABELS = ("ux", "uy", "uz", "phix", "phiy", "phiz")

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


def build_beam(use_full_mesh: bool) -> FEMHermiteBeamRegion:
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


def get_previous_shared_nodes(beam: FEMHermiteBeamRegion, elem_idx: int) -> set[int]:
    if elem_idx == 0:
        return set()
    return {
        int(node)
        for node in set(beam.elements[elem_idx]).intersection(beam.elements[elem_idx - 1])
    }


def format_element_dof_table(beam: FEMHermiteBeamRegion, elem_idx: int) -> str:
    node_indices = beam.elements[elem_idx]
    coords = beam.get_element_nodes(elem_idx)
    dof_grid = beam.get_element_global_dofs(elem_idx)
    dofs = dof_grid.ravel()
    shared_previous = get_previous_shared_nodes(beam, elem_idx)

    lines = [
        f"Element {elem_idx}",
        "local_node  global_node        x        y        z        ux    uy    uz  phix  phiy  phiz  shared_prev",
    ]

    for local_node, (global_node, xyz, global_dofs) in enumerate(
        zip(node_indices, coords, dof_grid)
    ):
        dof_text = " ".join(f"{int(dof):>5}" for dof in global_dofs)
        shared_text = "yes" if int(global_node) in shared_previous else ""
        lines.append(
            f"{local_node:>10}  {int(global_node):>11}  "
            f"{xyz[0]:>7.2f}  {xyz[1]:>7.2f}  {xyz[2]:>7.2f}  "
            f"{dof_text}  {shared_text:>11}"
        )

    if shared_previous:
        shared_dofs = beam.global_dofs[sorted(shared_previous)].ravel()
        lines.extend(
            [
                "",
                f"Shared nodes with element {elem_idx - 1}: {sorted(shared_previous)}",
                f"Repeated global DOFs from those shared nodes: {shared_dofs.tolist()}",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "Shared nodes with previous element: none",
            ]
        )

    lines.extend(
        [
            "",
            f"Flattened dofs length: {len(dofs)}",
            f"Flattened dofs: {dofs.tolist()}",
            "",
            "This element contributes Ke[a,b] to global K[dofs[a], dofs[b]].",
            "Neighboring elements share node numbers, so their dof maps reuse",
            "some of the same global row/column addresses.",
        ]
    )
    return "\n".join(lines)


def print_element_dof_table(beam: FEMHermiteBeamRegion, elem_idx: int) -> None:
    print()
    print("=" * 100)
    print(format_element_dof_table(beam, elem_idx))


def make_element_label_text(beam: FEMHermiteBeamRegion, elem_idx: int) -> str:
    node_indices = beam.elements[elem_idx]
    dof_grid = beam.get_element_global_dofs(elem_idx)
    shared_previous = get_previous_shared_nodes(beam, elem_idx)
    lines = [
        f"element {elem_idx} / {beam.get_num_elements() - 1}",
        "local -> global node -> dofs",
    ]

    for local_node, (global_node, global_dofs) in enumerate(zip(node_indices, dof_grid)):
        first_dof = int(global_dofs[0])
        last_dof = int(global_dofs[-1])
        shared_marker = " shared" if int(global_node) in shared_previous else ""
        lines.append(
            f"{local_node}: node {int(global_node)} -> {first_dof}-{last_dof}{shared_marker}"
        )

    return "\n".join(lines)


def add_node_spheres(
    plotter: pv.Plotter,
    points: np.ndarray,
    colors: list[str],
    radius: float,
) -> list:
    actors = []
    for point, color in zip(points, colors):
        sphere = pv.Sphere(
            radius=radius,
            center=point,
            theta_resolution=24,
            phi_resolution=24,
        )
        actors.append(
            plotter.add_mesh(
                sphere,
                color=color,
                opacity=0.92,
                name=None,
            )
        )
    return actors


def add_element_node_labels(
    plotter: pv.Plotter,
    beam: FEMHermiteBeamRegion,
    elem_idx: int,
) -> object:
    node_indices = beam.elements[elem_idx]
    points = beam.get_all_points()[node_indices]
    labels = [
        f"L{local}: N{int(global_node)}"
        for local, global_node in enumerate(node_indices)
    ]
    return plotter.add_point_labels(
        points,
        labels,
        font_size=16,
        point_color="#111111",
        text_color="#111111",
        shape_color="white",
        shape_opacity=0.78,
        always_visible=True,
    )


def make_plot(
    elem_idx: int,
    use_full_mesh: bool,
    screenshot: Path | None,
    sphere_radius: float,
) -> None:
    beam = build_beam(use_full_mesh=use_full_mesh)
    if elem_idx < 0 or elem_idx >= beam.get_num_elements():
        raise ValueError(f"Element {elem_idx} is outside 0..{beam.get_num_elements() - 1}")

    grid = build_unstructured_grid(beam)
    plotter = pv.Plotter(off_screen=screenshot is not None, window_size=(1700, 1000))
    state = {
        "elem_idx": elem_idx,
        "selected_actor": None,
        "node_actors": [],
        "label_actor": None,
        "text_actor": None,
    }

    plotter.add_mesh(
        grid,
        style="wireframe",
        color="#263238",
        line_width=1,
        opacity=0.24,
    )

    def remove_current_overlay() -> None:
        for key in ("selected_actor", "label_actor", "text_actor"):
            actor = state[key]
            if actor is not None:
                plotter.remove_actor(actor)
                state[key] = None

        for actor in state["node_actors"]:
            plotter.remove_actor(actor)
        state["node_actors"] = []

    def show_element(value: float | int) -> None:
        selected_idx = int(round(float(value)))
        selected_idx = max(0, min(selected_idx, beam.get_num_elements() - 1))
        state["elem_idx"] = selected_idx
        remove_current_overlay()

        selected = grid.extract_cells(selected_idx)
        state["selected_actor"] = plotter.add_mesh(
            selected,
            color="#56cc9d",
            opacity=0.32,
            show_edges=True,
            edge_color="#2f80ed",
            line_width=5,
        )

        node_indices = beam.elements[selected_idx]
        node_points = beam.get_all_points()[node_indices]
        shared_previous = get_previous_shared_nodes(beam, selected_idx)
        colors = [
            "#f2994a" if int(node_idx) in shared_previous else "#eb5757"
            for node_idx in node_indices
        ]
        state["node_actors"] = add_node_spheres(
            plotter,
            node_points,
            colors,
            sphere_radius,
        )
        state["label_actor"] = add_element_node_labels(plotter, beam, selected_idx)
        state["text_actor"] = plotter.add_text(
            make_element_label_text(beam, selected_idx),
            position="upper_left",
            font_size=11,
            color="#111111",
            font="courier",
        )

        print_element_dof_table(beam, selected_idx)
        plotter.render()

    def next_element() -> None:
        show_element(state["elem_idx"] + 1)

    def previous_element() -> None:
        show_element(state["elem_idx"] - 1)

    print("Element DOF map stepper")
    print(f"  mesh nodes:       {beam.get_all_points().shape[0]}")
    print(f"  elements:         {beam.get_num_elements()}")
    print(f"  dofs per node:    {beam.ndof_per_node} ({', '.join(DOF_LABELS)})")
    print("  controls: slider, n = next, p = previous")

    plotter.add_slider_widget(
        show_element,
        rng=[0, beam.get_num_elements() - 1],
        value=elem_idx,
        title="element index",
        pointa=(0.33, 0.08),
        pointb=(0.88, 0.08),
        style="modern",
        fmt="%.0f",
    )
    plotter.add_key_event("n", next_element)
    plotter.add_key_event("p", previous_element)
    plotter.add_axes()
    plotter.camera_position = [
        (12.0, -18.0, 12.0) if not use_full_mesh else (95.0, -135.0, 95.0),
        (0.0, 0.0, 3.0) if not use_full_mesh else (0.0, 0.0, 60.0),
        (0.0, 0.0, 1.0),
    ]

    if state["selected_actor"] is None:
        show_element(elem_idx)

    if screenshot:
        screenshot.parent.mkdir(parents=True, exist_ok=True)
        plotter.show(screenshot=str(screenshot))
        print(f"Saved screenshot: {screenshot}")
    else:
        plotter.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Step through elements and inspect each local-to-global DOF map."
    )
    parser.add_argument("--element", type=int, default=0, help="Initial element index.")
    parser.add_argument(
        "--full-mesh",
        action="store_true",
        help="Use the full solver mesh instead of the small teaching mesh.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Optional PNG path for an off-screen screenshot.",
    )
    parser.add_argument(
        "--sphere-radius",
        type=float,
        default=0.11,
        help="Radius for highlighted element-node spheres.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    make_plot(
        elem_idx=args.element,
        use_full_mesh=args.full_mesh,
        screenshot=args.screenshot,
        sphere_radius=args.sphere_radius,
    )
