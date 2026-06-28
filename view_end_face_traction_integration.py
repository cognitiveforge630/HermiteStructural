"""Interactive PyVista view of the Hermite end-face traction integration.

This viewer is intentionally diagnostic: it does not solve the structure and it
never changes the solver algorithm.  It repeats the same surface-integration
steps used by ``solver.build_end_face_traction_load_vector`` so the loaded
Hex8 end faces, their Gauss points, and the nodal generalized load quantities
can be inspected one face at a time.

Run::

    py view_end_face_traction_integration.py

Useful variants::

    py view_end_face_traction_integration.py --face 12
    py view_end_face_traction_integration.py --screenshot outputs/end_face_traction_face_12.png --face 12

Controls:
    Right / Up / N       next integrated end face
    Left / Down / P      previous integrated end face
    Home                 first integrated end face
    End                  last integrated end face
    R                    reset camera around the selected face
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyvista as pv

from solver import LOAD_DOF, P, build_beam


LOCAL_FACE_NODES = np.array([4, 5, 6, 7], dtype=int)
DOF_NAMES = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
DOF_UNITS = ("lb", "lb", "lb", "lb-in", "lb-in", "lb-in")


@dataclass(frozen=True)
class EndFaceLoad:
    """One element's zeta=+1 face and its integrated element-load vector."""

    face_index: int
    elem_idx: int
    local_face_nodes: np.ndarray
    global_face_nodes: np.ndarray
    face_points: np.ndarray
    gauss_points_xyz: np.ndarray
    gauss_weighted_areas: np.ndarray
    element_load: np.ndarray

    @property
    def local_nodal_loads(self) -> np.ndarray:
        """Return the 8 x 6 element load vector for local node/dof lookup."""
        return self.element_load.reshape((8, 6))

    @property
    def integrated_area(self) -> float:
        return float(np.sum(self.gauss_weighted_areas))

    def translational_component_sum(self, load_dof: int) -> float:
        return float(np.sum(self.local_nodal_loads[:, load_dof]))

    @property
    def center(self) -> np.ndarray:
        return self.face_points.mean(axis=0)


def build_unstructured_grid(beam) -> pv.UnstructuredGrid:
    """Convert the beam's Hex8 connectivity into a PyVista grid."""
    n_elements = beam.get_num_elements()
    cell_connectivity = np.column_stack(
        (np.full(n_elements, 8, dtype=np.int64), beam.elements)
    ).ravel()
    cell_types = np.full(n_elements, pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(cell_connectivity, cell_types, beam.get_all_points())


def build_end_face_polydata(beam, face_loads: list[EndFaceLoad]) -> pv.PolyData:
    """Return all integrated z=L quadrilateral faces as a single PolyData."""
    faces = []
    for face_load in face_loads:
        faces.extend([4, *face_load.global_face_nodes.tolist()])

    surface = pv.PolyData(beam.get_all_points(), np.asarray(faces, dtype=np.int64))
    surface.cell_data["face_index"] = np.arange(len(face_loads), dtype=int)
    surface.cell_data["element_index"] = np.array(
        [face_load.elem_idx for face_load in face_loads], dtype=int
    )
    return surface


def selected_face_polydata(face_load: EndFaceLoad) -> pv.PolyData:
    """Return the currently selected quadrilateral face as standalone geometry."""
    return pv.PolyData(
        face_load.face_points,
        np.array([4, 0, 1, 2, 3], dtype=np.int64),
    )


def integrate_one_end_face(
    beam,
    elem_idx: int,
    face_index: int,
    *,
    traction: float,
    load_dof: int,
) -> EndFaceLoad:
    """Integrate ``N_disp.T @ traction`` over one element's zeta=+1 face.

    This is the same numerical integration used in ``solver.py``:

    ``Fe += N_disp.T @ unit_traction * weight * detJs``
    """
    coords = beam.get_element_nodes(elem_idx)
    Fe = np.zeros(beam.ndof_per_node * 8)
    gauss_points_xyz: list[np.ndarray] = []
    gauss_weighted_areas: list[float] = []

    for i, xi in enumerate(beam.gauss_points):
        for j, eta in enumerate(beam.gauss_points):
            zeta = 1.0
            weight = beam.gauss_weights[i] * beam.gauss_weights[j]
            N_disp, _, _, _ = beam.get_hermite_displacement_matrices(
                xi, eta, zeta, coords
            )
            J = beam.get_hermite_jacobian(xi, eta, zeta, coords)
            dx_dxi = J[0, :]
            dx_deta = J[1, :]
            detJs = np.linalg.norm(np.cross(dx_dxi, dx_deta))

            unit_traction = np.zeros(3)
            unit_traction[load_dof] = traction
            Fe += N_disp.T @ unit_traction * weight * detJs

            N_hex = beam.hex8_shape_functions(xi, eta, zeta)
            gauss_points_xyz.append(N_hex @ coords)
            gauss_weighted_areas.append(float(weight * detJs))

    global_nodes = np.asarray(beam.elements[elem_idx], dtype=int)
    global_face_nodes = global_nodes[LOCAL_FACE_NODES]
    face_points = beam.get_all_points()[global_face_nodes]

    return EndFaceLoad(
        face_index=face_index,
        elem_idx=elem_idx,
        local_face_nodes=LOCAL_FACE_NODES.copy(),
        global_face_nodes=global_face_nodes,
        face_points=face_points,
        gauss_points_xyz=np.asarray(gauss_points_xyz, dtype=float),
        gauss_weighted_areas=np.asarray(gauss_weighted_areas, dtype=float),
        element_load=Fe,
    )


def integrate_all_end_faces(
    beam,
    *,
    total_load: float = P,
    load_dof: int = LOAD_DOF,
) -> tuple[list[EndFaceLoad], float]:
    """Return all z=L face integrations and the uniform traction value."""
    face_area = beam.Lx * beam.Ly
    traction = total_load / face_area
    face_loads = []

    for face_index, elem_idx in enumerate(beam.get_elements_at_max_z()):
        face_loads.append(
            integrate_one_end_face(
                beam,
                elem_idx,
                face_index,
                traction=traction,
                load_dof=load_dof,
            )
        )

    if not face_loads:
        raise ValueError("No elements were found on the beam's z=L end face.")

    return face_loads, traction


def assembled_load_sum(face_loads: list[EndFaceLoad], load_dof: int) -> float:
    """Return the scalar force component assembled by all displayed end faces."""
    return float(
        sum(face.translational_component_sum(load_dof) for face in face_loads)
    )


def force_label(value: float, unit: str) -> str:
    """Format tiny generalized load values without hiding their sign."""
    if abs(value) < 5e-13:
        value = 0.0
    return f"{value:+.6e} {unit}"


def node_face_load_label(
    beam,
    face_load: EndFaceLoad,
    local_node: int,
    global_node: int,
    *,
    load_dof: int,
) -> str:
    """Build one compact node label for the selected face."""
    local_load = face_load.local_nodal_loads[local_node]
    lines = [
        f"local {local_node} / global {global_node}",
        f"global dof {beam.global_dofs[global_node, load_dof]} ({DOF_NAMES[load_dof]})",
        f"{DOF_NAMES[load_dof]} face = {force_label(local_load[load_dof], DOF_UNITS[load_dof])}",
    ]

    # Surface traction through the Hermite displacement matrix can create
    # rotational generalized loads too.  Show the nonzero ones so they are not
    # mistaken for missing force when checking the 48-component Fe vector.
    for dof in range(3, 6):
        if abs(local_load[dof]) > 5e-13:
            lines.append(f"{DOF_NAMES[dof]} face = {force_label(local_load[dof], DOF_UNITS[dof])}")

    return "\n".join(lines)


class EndFaceTractionViewer:
    """Interactive keyboard viewer for end-face traction integration."""

    def __init__(
        self,
        beam,
        face_loads: list[EndFaceLoad],
        traction: float,
        *,
        total_load: float,
        load_dof: int,
        initial_face: int = 0,
        off_screen: bool = False,
        window_size: tuple[int, int] = (1600, 1000),
    ) -> None:
        if not 0 <= initial_face < len(face_loads):
            raise ValueError(
                f"Face {initial_face} is outside 0..{len(face_loads) - 1}."
            )

        self.beam = beam
        self.grid = build_unstructured_grid(beam)
        self.face_loads = face_loads
        self.traction = float(traction)
        self.total_load = float(total_load)
        self.load_dof = int(load_dof)
        self.face_index = int(initial_face)
        self.all_face_surface = build_end_face_polydata(beam, face_loads)
        self.plotter = pv.Plotter(off_screen=off_screen, window_size=window_size)
        self._dynamic_actors: list[object] = []

        self.plotter.set_background("#f6f8fb")
        self._add_static_geometry()
        self._register_keys()
        self.show_face(self.face_index, reset_view=True)

    @property
    def selected_face(self) -> EndFaceLoad:
        return self.face_loads[self.face_index]

    def _add_static_geometry(self) -> None:
        self.plotter.add_mesh(
            self.grid,
            style="wireframe",
            color="#263238",
            opacity=0.18,
            line_width=0.6,
        )
        self.plotter.add_mesh(
            self.all_face_surface,
            color="#2b6cb0",
            opacity=0.20,
            show_edges=True,
            edge_color="#0f2742",
            line_width=1.3,
        )
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
        self.plotter.add_text(
            "End-face traction integration viewer",
            position="upper_edge",
            font_size=14,
            color="#172033",
            name="viewer-title",
        )
        self.plotter.add_text(
            "Right/Up/N: next    Left/Down/P: previous    Home/End: jump    R: reset view",
            position="lower_edge",
            font_size=9,
            color="#34445d",
            name="viewer-controls",
        )

    def _register_keys(self) -> None:
        self.plotter.add_key_event("Right", self.next_face)
        self.plotter.add_key_event("Up", self.next_face)
        self.plotter.add_key_event("n", self.next_face)
        self.plotter.add_key_event("Left", self.previous_face)
        self.plotter.add_key_event("Down", self.previous_face)
        self.plotter.add_key_event("p", self.previous_face)
        self.plotter.add_key_event("Home", self.first_face)
        self.plotter.add_key_event("End", self.last_face)
        self.plotter.add_key_event("r", self.reset_selected_face_view)

    def _remove_dynamic_actors(self) -> None:
        for actor in self._dynamic_actors:
            self.plotter.remove_actor(actor, render=False)
        self._dynamic_actors.clear()

    def _add_selected_cell(self) -> None:
        selected_cell = self.grid.extract_cells([self.selected_face.elem_idx])
        actor = self.plotter.add_mesh(
            selected_cell,
            color="#90cdf4",
            opacity=0.18,
            show_edges=True,
            edge_color="#123b66",
            line_width=3.0,
        )
        self._dynamic_actors.append(actor)

    def _add_selected_face(self) -> None:
        # The selected integration face is deliberately red and transparent so
        # it can be seen both as a surface and as part of its owning element.
        face_actor = self.plotter.add_mesh(
            selected_face_polydata(self.selected_face),
            color="#e3342f",
            opacity=0.55,
            show_edges=True,
            edge_color="#7f1d1d",
            line_width=4.0,
        )
        point_actor = self.plotter.add_points(
            self.selected_face.face_points,
            color="#111111",
            point_size=16,
            render_points_as_spheres=True,
        )
        self._dynamic_actors.extend((face_actor, point_actor))

    def _add_gauss_points(self) -> None:
        actor = self.plotter.add_points(
            self.selected_face.gauss_points_xyz,
            color="#ffd166",
            point_size=13,
            render_points_as_spheres=True,
        )
        self._dynamic_actors.append(actor)

    def _add_load_arrow(self) -> None:
        direction = np.zeros(3)
        sign = 1.0 if self.traction >= 0.0 else -1.0
        element_lengths = np.ptp(self.selected_face.face_points, axis=0)
        visible_scale = max(float(element_lengths.max()), 1.0)
        direction[self.load_dof] = sign * visible_scale * 0.85
        actor = self.plotter.add_arrows(
            self.selected_face.center[None, :],
            direction[None, :],
            mag=1.0,
            color="#b91c1c",
        )
        self._dynamic_actors.append(actor)

    def _add_node_labels(self) -> None:
        labels = []
        for local_node, global_node in zip(
            self.selected_face.local_face_nodes,
            self.selected_face.global_face_nodes,
        ):
            labels.append(
                node_face_load_label(
                    self.beam,
                    self.selected_face,
                    int(local_node),
                    int(global_node),
                    load_dof=self.load_dof,
                )
            )

        actor = self.plotter.add_point_labels(
            self.selected_face.face_points,
            labels,
            font_size=9,
            text_color="#101820",
            point_size=0,
            shape="rounded_rect",
            shape_color="white",
            shape_opacity=0.90,
            fill_shape=True,
            always_visible=True,
            margin=5,
        )
        self._dynamic_actors.append(actor)

    def _add_summary_panel(self) -> None:
        face = self.selected_face
        face_component = face.translational_component_sum(self.load_dof)
        total_component = assembled_load_sum(self.face_loads, self.load_dof)
        expected_face_component = self.traction * face.integrated_area
        weighted_area_min = float(np.min(face.gauss_weighted_areas))
        weighted_area_max = float(np.max(face.gauss_weighted_areas))

        lines = [
            f"selected face index:       {face.face_index} / {len(self.face_loads) - 1}",
            f"owning element index:      {face.elem_idx}",
            f"local face nodes:          {face.local_face_nodes.tolist()}  (zeta = +1)",
            f"global face nodes:         {face.global_face_nodes.tolist()}",
            f"load component:            {DOF_NAMES[self.load_dof]}",
            f"uniform traction:          {self.traction:+.6e} lb/in^2",
            f"integrated face area:      {face.integrated_area:.6e} in^2",
            f"weight*detJs range:        {weighted_area_min:.6e} .. {weighted_area_max:.6e}",
            f"face {DOF_NAMES[self.load_dof]} sum:          {face_component:+.6e} lb",
            f"traction*face area:        {expected_face_component:+.6e} lb",
            f"all end-face {DOF_NAMES[self.load_dof]} sum:  {total_component:+.6e} lb",
            f"requested total load:      {self.total_load:+.6e} lb",
        ]
        actor = self.plotter.add_text(
            "\n".join(lines),
            position="upper_right",
            font_size=8,
            color="#172033",
            font="courier",
            shadow=True,
        )
        self._dynamic_actors.append(actor)

    def show_face(self, face_index: int, *, reset_view: bool = False) -> None:
        self.face_index = face_index % len(self.face_loads)
        self._remove_dynamic_actors()
        self._add_selected_cell()
        self._add_selected_face()
        self._add_gauss_points()
        self._add_load_arrow()
        self._add_node_labels()
        self._add_summary_panel()
        if reset_view:
            self.reset_selected_face_view(render=False)
        self.plotter.render()

    def next_face(self) -> None:
        self.show_face(self.face_index + 1)

    def previous_face(self) -> None:
        self.show_face(self.face_index - 1)

    def first_face(self) -> None:
        self.show_face(0)

    def last_face(self) -> None:
        self.show_face(len(self.face_loads) - 1)

    def reset_selected_face_view(self, *, render: bool = True) -> None:
        points = self.selected_face.face_points
        center = points.mean(axis=0)
        lengths = np.ptp(points, axis=0)
        max_extent = max(float(lengths.max()), 1.0)
        minimum = points.min(axis=0) - max_extent * 0.75
        maximum = points.max(axis=0) + max_extent * 0.75
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
        self.plotter.camera.zoom(0.55)
        if render:
            self.plotter.render()

    def show(self, screenshot: Path | None = None) -> None:
        if screenshot is not None:
            screenshot.parent.mkdir(parents=True, exist_ok=True)
            self.plotter.show(screenshot=str(screenshot))
            print(f"Saved screenshot: {screenshot}")
        else:
            self.plotter.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visualize the z=L Hermite end-face traction integration, selected "
            "faces, Gauss points, and resulting nodal generalized loads."
        )
    )
    parser.add_argument(
        "--face",
        type=int,
        default=0,
        help="Initial integrated face index, from 0 to number_of_end_faces-1.",
    )
    parser.add_argument(
        "--total-load",
        type=float,
        default=P,
        help="Total load in lb distributed over the whole z=L face.",
    )
    parser.add_argument(
        "--load-dof",
        type=int,
        choices=(0, 1, 2),
        default=LOAD_DOF,
        help="Translational component receiving the traction: 0=Fx, 1=Fy, 2=Fz.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Optional PNG output; enables off-screen rendering.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    beam = build_beam()
    face_loads, traction = integrate_all_end_faces(
        beam,
        total_load=args.total_load,
        load_dof=args.load_dof,
    )

    total_component = assembled_load_sum(face_loads, args.load_dof)
    print("End-face traction integration")
    print(f"  integrated faces:       {len(face_loads)}")
    print(f"  total load requested:   {args.total_load:+.12e} lb")
    print(f"  uniform traction:       {traction:+.12e} lb/in^2")
    print(f"  assembled component:    {total_component:+.12e} lb")
    print(f"  component error:        {total_component - args.total_load:+.12e} lb")
    print("  keyboard: Right/Up/N next, Left/Down/P previous, Home/End jump, R reset")

    viewer = EndFaceTractionViewer(
        beam,
        face_loads,
        traction,
        total_load=args.total_load,
        load_dof=args.load_dof,
        initial_face=args.face,
        off_screen=args.screenshot is not None,
    )
    viewer.show(args.screenshot)


if __name__ == "__main__":
    main()
