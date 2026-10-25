"""Interactive PyVista inspector for Hex8 elements in the structural mesh.

Run::

    py element_mesh_viewer.py

Controls:
    Right / Up     next element
    Left / Down    previous element
    R              reset the selected element to fill the view
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Protocol

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


class HexMesh(Protocol):
    """The small mesh interface consumed by :class:`HexElementViewer`."""

    elements: np.ndarray

    def get_all_points(self) -> np.ndarray: ...

    def get_num_elements(self) -> int: ...

    def get_element_global_dofs(self, elem_idx: int) -> np.ndarray: ...


def build_beam() -> FEMHermiteBeamRegion:
    """Build the same mesh used by ``solver.py`` and ``view.py``."""
    return FEMHermiteBeamRegion(**BEAM_PARAMS)


def build_grid(mesh: HexMesh) -> pv.UnstructuredGrid:
    """Convert the mesh's global Hex8 connectivity to a PyVista grid."""
    element_count = mesh.get_num_elements()
    connectivity = np.column_stack(
        (np.full(element_count, 8, dtype=np.int64), mesh.elements)
    ).ravel()
    cell_types = np.full(
        element_count, pv.CellType.HEXAHEDRON, dtype=np.uint8
    )
    return pv.UnstructuredGrid(connectivity, cell_types, mesh.get_all_points())


def _coordinate_text(value: float) -> str:
    """Keep coordinate labels compact without hiding small values."""
    return f"{value:.6g}"


class HexElementViewer:
    """Display and navigate one element from a global hexahedral mesh."""

    def __init__(
        self,
        mesh: HexMesh,
        initial_element: int = 0,
        *,
        off_screen: bool = False,
        window_size: tuple[int, int] = (1500, 900),
    ) -> None:
        self.mesh = mesh
        self.grid = build_grid(mesh)
        self.element_count = mesh.get_num_elements()
        if self.element_count == 0:
            raise ValueError("The mesh does not contain any elements.")
        if not 0 <= initial_element < self.element_count:
            raise ValueError(
                f"Element {initial_element} is outside 0..{self.element_count - 1}."
            )

        self.element_index = initial_element
        self.plotter = pv.Plotter(
            off_screen=off_screen,
            window_size=window_size,
            border=False,
        )
        self._dynamic_actors: list[object] = []

        self.plotter.set_background("#f4f7fb")
        self.plotter.add_text(
            "Hex8 element inspector",
            position="upper_edge",
            font_size=15,
            color="#172033",
            name="viewer-title",
        )
        self.plotter.add_text(
            "Arrow keys: previous / next     R: fit element",
            position="lower_edge",
            font_size=10,
            color="#34445d",
            name="viewer-controls",
        )

        # A corner orientation widget acts as the UCS icon and follows the camera.
        self.plotter.add_axes(
            interactive=True,
            line_width=3,
            x_color="#d62828",
            y_color="#2a9d3f",
            z_color="#246bfd",
            xlabel="X",
            ylabel="Y",
            zlabel="Z",
        )

        self.plotter.add_key_event("Right", self.next_element)
        self.plotter.add_key_event("Up", self.next_element)
        self.plotter.add_key_event("Left", self.previous_element)
        self.plotter.add_key_event("Down", self.previous_element)
        self.plotter.add_key_event("r", self.fit_element)

        initial_points = np.asarray(mesh.get_all_points())[mesh.elements[initial_element]]
        center = initial_points.mean(axis=0)
        extent = max(float(np.ptp(initial_points, axis=0).max()), 1.0)
        self.plotter.camera.position = center + extent * np.array([2.4, -3.2, 2.0])
        self.plotter.camera.focal_point = center
        self.plotter.camera.up = (0.0, 0.0, 1.0)
        self.show_element(initial_element, reset_view=True)

    @property
    def node_indices(self) -> np.ndarray:
        return np.asarray(self.mesh.elements[self.element_index], dtype=int)

    @property
    def node_points(self) -> np.ndarray:
        return np.asarray(self.mesh.get_all_points())[self.node_indices]

    @property
    def node_global_dofs(self) -> np.ndarray:
        return np.asarray(
            self.mesh.get_element_global_dofs(self.element_index), dtype=int
        )

    def _remove_dynamic_actors(self) -> None:
        for actor in self._dynamic_actors:
            self.plotter.remove_actor(actor, render=False)
        self._dynamic_actors.clear()

    def _add_element_geometry(self) -> None:
        selected_cell = self.grid.extract_cells([self.element_index])
        element_actor = self.plotter.add_mesh(
            selected_cell,
            color="#77bdfb",
            opacity=0.25,
            show_edges=True,
            edge_color="#123b66",
            line_width=4,
            lighting=True,
        )

        lengths = np.ptp(self.node_points, axis=0)
        nonzero_lengths = lengths[lengths > 0.0]
        scale = float(nonzero_lengths.min()) if nonzero_lengths.size else 1.0
        node_actor = self.plotter.add_points(
            self.node_points,
            color="#f05a47",
            point_size=max(12.0, min(28.0, 20.0 * scale / lengths.max())),
            render_points_as_spheres=True,
            lighting=True,
        )
        self._dynamic_actors.extend((element_actor, node_actor))

    def _add_node_labels(self) -> None:
        labels = []
        for local_node, (global_node, global_dofs, point) in enumerate(
            zip(self.node_indices, self.node_global_dofs, self.node_points)
        ):
            xyz = ", ".join(_coordinate_text(value) for value in point)
            dofs = ", ".join(str(dof) for dof in global_dofs)
            labels.append(
                f"local index: ({local_node})\n"
                f"global index: ({global_node})\n"
                f"global dofs: ({dofs})\n"
                f"global xyz: ({xyz})"
            )

        label_actor = self.plotter.add_point_labels(
            self.node_points,
            labels,
            font_size=10,
            text_color="#101820",
            point_size=0,
            shape="rounded_rect",
            shape_color="white",
            shape_opacity=0.88,
            fill_shape=True,
            always_visible=True,
            margin=4,
        )
        self._dynamic_actors.append(label_actor)

    def _add_coordinate_panel(self) -> None:
        lines = [
            "Element node map",
            "local : global : global dofs                 : global (x, y, z)",
        ]
        for local_node, (global_node, global_dofs, point) in enumerate(
            zip(self.node_indices, self.node_global_dofs, self.node_points)
        ):
            dofs = ",".join(str(dof) for dof in global_dofs)
            xyz = ", ".join(f"{value:.4g}" for value in point)
            lines.append(
                f"  {local_node}   : {global_node:6d} : "
                f"({dofs:<27}) : ({xyz})"
            )

        coordinate_actor = self.plotter.add_text(
            "\n".join(lines),
            position="upper_right",
            font_size=8,
            color="#172033",
            font="courier",
            shadow=True,
        )
        self._dynamic_actors.append(coordinate_actor)

    def _add_element_legend(self) -> None:
        legend_actor = self.plotter.add_text(
            f"ELEMENT INDEX  {self.element_index} / {self.element_count - 1}",
            position="lower_right",
            font_size=12,
            color="#246bfd",
            font="courier",
            shadow=True,
        )
        self._dynamic_actors.append(legend_actor)

    def show_element(self, element_index: int, *, reset_view: bool = True) -> None:
        """Change the selected element and refresh all dependent overlays."""
        self.element_index = element_index % self.element_count
        self._remove_dynamic_actors()
        self._add_element_geometry()
        self._add_node_labels()
        self._add_coordinate_panel()
        self._add_element_legend()
        if reset_view:
            self.fit_element(render=False)
        self.plotter.render()

    def next_element(self) -> None:
        self.show_element(self.element_index + 1)

    def previous_element(self) -> None:
        self.show_element(self.element_index - 1)

    def fit_element(self, *, render: bool = True) -> None:
        """Center the selected Hex8 and scale it to fill the client area."""
        center = self.node_points.mean(axis=0)
        minimum = self.node_points.min(axis=0)
        maximum = self.node_points.max(axis=0)
        bounds = (
            minimum[0],
            maximum[0],
            minimum[1],
            maximum[1],
            minimum[2],
            maximum[2],
        )
        self.plotter.camera.focal_point = center
        self.plotter.camera.view_angle = 30.0
        self.plotter.reset_camera(bounds=bounds)
        # Leave a margin for the eight callout labels while still making the
        # element the dominant object in the client area.
        self.plotter.camera.zoom(0.60)
        if render:
            self.plotter.render()

    def show(self, screenshot: Path | None = None) -> None:
        if screenshot is not None:
            screenshot.parent.mkdir(parents=True, exist_ok=True)
            self.plotter.show(screenshot=str(screenshot))
        else:
            self.plotter.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect individual Hex8 elements from the Hermite beam mesh."
    )
    parser.add_argument("--element", type=int, default=0, help="Initial element index.")
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Optional PNG output; enables off-screen rendering.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mesh = build_beam()
    viewer = HexElementViewer(
        mesh,
        initial_element=args.element,
        off_screen=args.screenshot is not None,
    )
    viewer.show(args.screenshot)


if __name__ == "__main__":
    main()
