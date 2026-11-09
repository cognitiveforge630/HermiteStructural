"""Interactive PyVista viewer and ParaView-compatible VTU exporter.

Plot input is a *sampled* high-order Hermite field, not an independent
low-order FE recomputation. Exported geometry remains undeformed, with actual
Hermite displacement vectors and recovered Cauchy stress evaluated at subcell corners.

Importing this module does not require PyVista; the dependency is required
only for rendering and VTU output.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from .sampling import DEFAULT_FIELD, CELL_FIELDS, POINT_FIELDS, SampledVolume

VECTOR_FIELD = "Displacement (in)"


def _pyvista():
    try:
        import pyvista as pv
    except ImportError as exc:
        raise RuntimeError(
            "PyVista is required for the viewer / VTU export. Install it with "
            "'python -m pip install -r requirements.txt' from the workspace root."
        ) from exc
    return pv


def as_pyvista_grid(sample: SampledVolume):
    """Convert a validated sampled high-order field to VTK hexahedral cells."""
    pv = _pyvista()
    if len(sample.connectivity)==0:
        raise ValueError("There are no cells to render")
    nc=len(sample.connectivity)
    cells=np.column_stack((np.full(nc,8,dtype=np.int64),sample.connectivity))
    types=np.full(nc,pv.CellType.HEXAHEDRON,dtype=np.uint8)
    grid=pv.UnstructuredGrid(cells.ravel(),types,sample.points.copy())
    grid.point_data[VECTOR_FIELD]=sample.displacement.copy()
    for name,values in sample.point_fields.items():
        grid.point_data[name]=np.asarray(values)
    for name,values in sample.cell_fields.items():
        grid.cell_data[name]=np.asarray(values)
    grid.field_data["Hermite visual subdivisions per element"] = np.array([sample.subdivisions])
    grid.field_data["Original Hermite element count"] = np.array([sample.element_count])
    return grid


def export_vtu(sample: SampledVolume, path):
    """Write a true VTK unstructured hexahedron file even on headless computers.

    Use PyVista when available; fall back to the installed vtk Python package
    (same VTK format) for research installations without a graphics frontend.
    Neither path interpolates results using ordinary FE Hex8 postprocessing.
    """
    path=Path(path)
    if path.suffix.lower()!=".vtu":
        raise ValueError("VTU export requires a .vtu filename")
    path.parent.mkdir(parents=True,exist_ok=True)
    try:
        import pyvista
    except ImportError:
        try:
            import vtk
            from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        except ImportError as exc:
            raise RuntimeError("VTU export requires either pyvista or vtk") from exc
        grid=vtk.vtkUnstructuredGrid()
        coords=vtk.vtkPoints()
        coords.SetData(numpy_to_vtk(np.ascontiguousarray(sample.points),deep=True))
        grid.SetPoints(coords)
        for conn in sample.connectivity:
            cell=vtk.vtkHexahedron()
            for k,node in enumerate(conn):
                cell.GetPointIds().SetId(k,int(node))
            grid.InsertNextCell(cell.GetCellType(),cell.GetPointIds())
        arrays=[(grid.GetPointData(),VECTOR_FIELD,sample.displacement)]
        arrays += [(grid.GetPointData(),name,arr) for name,arr in sample.point_fields.items()]
        arrays += [(grid.GetCellData(),name,arr) for name,arr in sample.cell_fields.items()]
        for holder,name,values in arrays:
            arr=numpy_to_vtk(np.ascontiguousarray(values),deep=True)
            arr.SetName(name)
            holder.AddArray(arr)
        # Preserve provenance metadata on the direct VTK export path too.
        for name, values in (("Hermite visual subdivisions per element", [sample.subdivisions]),
                             ("Original Hermite element count", [sample.element_count])):
            arr=numpy_to_vtk(np.asarray(values,dtype=np.int64),deep=True)
            arr.SetName(name)
            grid.GetFieldData().AddArray(arr)
        writer=vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(str(path))
        writer.SetInputData(grid)
        if writer.Write()!=1:
            raise IOError(f"Could not export VTU: {path}")
    else:
        as_pyvista_grid(sample).save(str(path))
    return path


def _validate_field(grid,field):
    if field not in grid.point_data and field not in grid.cell_data:
        known=", ".join([*POINT_FIELDS,*CELL_FIELDS])
        raise ValueError(f"Unknown field '{field}'. Available fields: {known}")


def _auto_scale(grid):
    vec=np.asarray(grid.point_data[VECTOR_FIELD])
    largest=float(np.max(np.linalg.norm(vec,axis=1)))
    extent=float(np.ptp(grid.points,axis=0).max())
    return (0.12*extent/largest if largest>0 else 1.0)


# All numeric plots use blue = lower scalar and red = higher scalar.
# For the DEFAULT positive-magnitude displacement, that also means
# blue = least movement and red = greatest movement.
COLOR_MAP = "coolwarm"
FONT_SIZE = 10
VERTICAL_MAGNITUDE = "Vertical displacement magnitude (in)"
FIXED_FLAG = "Fixed support node"
LOADED_FLAG = "Loaded end node"


def _plot_field_class(field):
    """Human-readable legend semantics; no inferred wood strength claims."""
    if field in (VERTICAL_MAGNITUDE, "Displacement magnitude (in)"):
        return "Magnitude | blue: least movement | red: greatest movement"
    if field == "Absolute longitudinal stress (psi)":
        return "Absolute longitudinal stress | blue: least | red: greatest"
    if "utilization" in field.lower():
        return "Utilization | blue: lowest | red: highest (research screening)"
    if field == "Von Mises (psi)":
        return "Equivalent stress | blue: lowest | red: highest (NOT wood failure)"
    return "Signed or scalar value | blue: lower/negative | red: higher/positive"


def _prepare_grid_fields(grid):
    """Keep old VTUs viewable while using positive displacement magnitudes."""
    if "Uy (in)" in grid.point_data and VERTICAL_MAGNITUDE not in grid.point_data:
        grid.point_data[VERTICAL_MAGNITUDE] = np.abs(np.asarray(grid.point_data["Uy (in)"]))
    if "Sigma zz (psi)" in grid.cell_data and "Absolute longitudinal stress (psi)" not in grid.cell_data:
        grid.cell_data["Absolute longitudinal stress (psi)"] = np.abs(np.asarray(grid.cell_data["Sigma zz (psi)"]))
    if VECTOR_FIELD in grid.point_data and "Displacement magnitude (in)" not in grid.point_data:
        grid.point_data["Displacement magnitude (in)"] = np.linalg.norm(
            np.asarray(grid.point_data[VECTOR_FIELD]), axis=1)
    return grid


def _end_source_points(grid, name):
    """Prefer explicit source-node flags; fall back to the four corner nodes.

    Older standalone VTUs did not record the original node grid. The fallback
    marks four unambiguous clamp corners rather than guessing a source mesh.
    """
    points = np.asarray(grid.points)
    end_z = float(np.min(points[:, 2]) if name == FIXED_FLAG else np.max(points[:, 2]))
    face = points[np.isclose(points[:, 2], end_z, atol=1e-8, rtol=0.)]
    if name in grid.point_data:
        mask = np.asarray(grid.point_data[name]).ravel() > 0.5
        selected = points[mask]
        if len(selected) and np.all(np.isclose(selected[:, 2], end_z, atol=1e-8)):
            return selected.copy()
    xy = face[:, :2]
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    corners = np.array([[xmin, ymin],[xmin, ymax],[xmax,ymin],[xmax,ymax]])
    return np.array([face[np.argmin(np.sum((xy - corner)**2,axis=1))]
                     for corner in corners])


def _display_limits(grid, field):
    values = np.asarray(grid.point_data[field] if field in grid.point_data
                        else grid.cell_data[field], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Non-finite plot field: {field}")
    lo, hi = float(np.min(values)), float(np.max(values))
    if field in (VERTICAL_MAGNITUDE, "Displacement magnitude (in)",
                 "Absolute longitudinal stress (psi)", "Von Mises (psi)") or "utilization" in field.lower():
        lo = 0.0
    if hi <= lo:
        hi = lo + max(1e-12, abs(lo)*1e-9)
    return (lo, hi)


def show_grid(grid, *, formulation="Recovered-rotation Hermite", field=DEFAULT_FIELD,
              deform_scale=None, show_edges=False, full_screen=False,
              off_screen=False, screenshot=None):
    """Interactive PyVista contour viewer with visible supports and a clear legend.

    1: vertical displacement MAGNITUDE |Uy| (blue=least, red=greatest);
    2: displacement magnitude |U|; U: signed Uy; A: absolute sigma zz; 3: von Mises;
    4/5: principal stress; 6/7/8: normal stress; 9: Uz; 0: beam reference;
    B: FE minus beam; T: wood longitudinal tension utilization;
    D: deformation; +/-: magnification; C: section; E: edges;
    M: support marker visibility; S: screenshot; R: camera.
    """
    pv = _pyvista()
    _prepare_grid_fields(grid)
    _validate_field(grid, field)
    if VECTOR_FIELD not in grid.point_data:
        raise ValueError(f"The mesh must contain vector point array '{VECTOR_FIELD}'")
    if deform_scale is not None and (not np.isfinite(deform_scale) or deform_scale < 0):
        raise ValueError("deform_scale must be finite and nonnegative")
    default_screenshot = Path("outputs") / "viewer_screenshot.png"
    screenshot = Path(screenshot) if screenshot is not None else (
        default_screenshot if off_screen else None)
    if screenshot is not None:
        screenshot.parent.mkdir(parents=True, exist_ok=True)

    state = {"field": field, "scale": float(_auto_scale(grid) if deform_scale is None else deform_scale),
             "deformed": True, "edges": bool(show_edges), "section": False,
             "markers": True, "scalar_bar_title": None}
    points = np.asarray(grid.points).copy()
    vectors = np.asarray(grid.point_data[VECTOR_FIELD]).copy()
    mesh = grid.copy(deep=True)
    fixed_nodes = _end_source_points(grid, FIXED_FLAG)
    loaded_nodes = _end_source_points(grid, LOADED_FLAG)
    z0, zL = float(points[:, 2].min()), float(points[:, 2].max())
    length = zL - z0
    cross_depth = max(float(np.ptp(points[:,1])), 1e-6)
    cube_size = min(cross_depth * .075, max(length * .006, cross_depth * .03))
    # Distinctive physical cube markers at the fixed-end SOURCE FE nodes.
    support = pv.PolyData(fixed_nodes).glyph(
        geom=pv.Cube(x_length=cube_size, y_length=cube_size, z_length=cube_size),
        scale=False, orient=False)

    p = pv.Plotter(off_screen=off_screen, window_size=(1500, 900))
    p.set_background("#0c1726")
    p.add_axes()
    keys = {"1": VERTICAL_MAGNITUDE, "2": "Displacement magnitude (in)",
            "u": "Uy (in)", "a": "Absolute longitudinal stress (psi)",
            "3": "Von Mises (psi)",
            "4": "Maximum principal (psi)", "5": "Minimum principal (psi)",
            "6": "Sigma zz (psi)", "7": "Sigma xx (psi)",
            "8": "Sigma yy (psi)", "9": "Uz (in)",
            "0": "Beam reference sigma zz (psi)",
            "b": "Cubic minus beam sigma zz (psi)",
            "t": "Longitudinal tension utilization"}

    def reset_view():
        # Look from -X so +Z runs from left to right and +Y points up.
        center = mesh.center
        p.camera_position = [
            (center[0] - 1.65*length, center[1] + .28*length, center[2] + .1*length),
            (center[0], center[1], center[2]),
            (0., 1., 0.),
        ]
        p.reset_camera()

    # Retain returned VTK actors directly.  PyVista's actor-name mapping is
    # version-dependent, especially for point labels (vtkActor2D).
    support_actor = p.add_mesh(support, color="#ffd166", lighting=False,
                               show_edges=True, edge_color="#20252d", line_width=1,
                               name="fixed_node_cubes", render=False)
    annotations = {"labels": None}

    def draw_annotations():
        support_actor.SetVisibility(bool(state["markers"]))
        if annotations["labels"] is not None:
            p.remove_actor(annotations["labels"], reset_camera=False, render=False)
            annotations["labels"] = None
        free_offset = np.array([0., cross_depth*.9, 0.])
        support_anchor = fixed_nodes.mean(axis=0) + free_offset
        loaded_anchor = loaded_nodes.mean(axis=0) + free_offset
        if state["deformed"]:
            # Label the loaded end where it appears on the deformed image.
            full_index = int(np.argmin(np.linalg.norm(points - loaded_nodes.mean(axis=0),axis=1)))
            loaded_anchor = loaded_anchor + state["scale"] * vectors[full_index]
        label_actor = p.add_point_labels(np.array([support_anchor, loaded_anchor]),
                           [f"FIXED SUPPORT  z={z0:g} in", f"FREE / LOADED  z={zL:g} in"],
                           font_size=FONT_SIZE, text_color="white",
                           shape_color="#263a50", shape="rect",
                           always_visible=True, show_points=False,
                           reset_camera=False, name="end_labels", render=False)
        annotations["labels"] = label_actor
        if label_actor is not None:
            label_actor.SetVisibility(bool(state["markers"]))

    def redraw(*, reset=False):
        camera = None if reset else p.camera_position
        mesh.points = points + state["scale"]*vectors if state["deformed"] else points.copy()
        scalar = state["field"]
        drawn = (mesh.slice(normal=(0,0,1), origin=(0,0,float(np.median(points[:,2]))))
                 if state["section"] else mesh)
        actor = p.add_mesh(drawn, scalars=scalar,
                           clim=_display_limits(mesh, scalar),
                           show_edges=state["edges"], cmap=COLOR_MAP,
                           lighting=True, show_scalar_bar=False,
                           name="hermite_mesh", render=False)
        previous = state["scalar_bar_title"]
        if previous is not None:
            try:
                p.remove_scalar_bar(title=previous, render=False)
            except (KeyError, ValueError):
                pass
        # Explicit WHITE scalar-bar text. On the previous dark background,
        # inherited black text made the legend practically unreadable.
        bar_title = f"{scalar} | blue = low, red = high"
        p.add_scalar_bar(title=bar_title, mapper=actor.mapper, color="white",
                         title_font_size=FONT_SIZE, label_font_size=FONT_SIZE,
                         font_family="arial", fmt="%.3g", n_labels=5,
                         vertical=False, position_x=.29, position_y=.045,
                         width=.55, height=.10, render=False)
        state["scalar_bar_title"] = bar_title
        p.add_text(
            f"{formulation} | {scalar}\n"
            f"{_plot_field_class(scalar)}\n"
            f"Deformation: {'ON' if state['deformed'] else 'OFF'}  "
            f"magnification: {state['scale']:.3g}x (display only)\n"
            "1 |Uy|  2 |U|  U signed Uy  A |sigma zz|  3 VM  4/5 principal\n"
            "6-8 signed stress  9 Uz  0 beam stress  B cubic-minus-beam  T wood util\n"
            "C section  D deform  +/- scale  E edges  M fixed markers\n"
            "S screenshot  R camera\n"
            "RESEARCH: stress reference agreement does not certify the 3D field",
            position="upper_left", font_size=FONT_SIZE,
            color="white", name="hermite_hud")
        draw_annotations()
        if reset:
            reset_view()
        elif camera is not None:
            p.camera_position = camera
        p.render()

    for key, field_name in keys.items():
        if field_name in grid.point_data or field_name in grid.cell_data:
            def make_callback(which):
                def callback():
                    state["field"] = which
                    redraw()
                return callback
            p.add_key_event(key, make_callback(field_name))

    def change_scale(factor):
        state["scale"] *= factor
        redraw()
    def take_screenshot():
        target = screenshot or default_screenshot
        target.parent.mkdir(parents=True, exist_ok=True)
        p.screenshot(str(target))
        print(f"Saved screenshot: {target}", flush=True)
    p.add_key_event("c", lambda: (state.update(section=not state["section"]), redraw()))
    p.add_key_event("d", lambda: (state.update(deformed=not state["deformed"]), redraw()))
    p.add_key_event("e", lambda: (state.update(edges=not state["edges"]), redraw()))
    p.add_key_event("m", lambda: (state.update(markers=not state["markers"]), redraw()))
    p.add_key_event("plus", lambda: change_scale(1.25))
    p.add_key_event("equal", lambda: change_scale(1.25))
    p.add_key_event("minus", lambda: change_scale(1/1.25))
    p.add_key_event("s", take_screenshot)
    p.add_key_event("r", lambda: redraw(reset=True))
    redraw(reset=True)
    p.show(full_screen=(full_screen and not off_screen),
           screenshot=str(screenshot) if screenshot is not None else None,
           auto_close=True)
    if screenshot is not None:
        print(f"Saved screenshot: {screenshot}")
    return screenshot

def show_sample(sample, **kwargs):
    """Create the PyVista grid and display it, preserving true Hermite fields."""
    return show_grid(as_pyvista_grid(sample), **kwargs)


def main(argv=None):
    """View a previously exported VTU without re-running a costly FEM solve."""
    p=argparse.ArgumentParser(description="Open an exported HermiteStructural VTU in PyVista")
    p.add_argument("vtu",type=Path)
    p.add_argument("--field",default=DEFAULT_FIELD,help="Array name, including cubic beam-theory diagnostic fields")
    p.add_argument("--deform-scale",type=float,default=None)
    p.add_argument("--edges",action="store_true")
    p.add_argument("--full-screen",action="store_true")
    p.add_argument("--off-screen",action="store_true")
    p.add_argument("--screenshot",type=Path,default=None)
    a=p.parse_args(argv)
    pv=_pyvista()
    grid=pv.read(str(a.vtu))
    show_grid(grid,formulation=a.vtu.stem,field=a.field,
              deform_scale=a.deform_scale,show_edges=a.edges,
              full_screen=a.full_screen,off_screen=a.off_screen,
              screenshot=a.screenshot)

if __name__ == "__main__":
    main()
