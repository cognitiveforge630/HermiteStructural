"""Two-pane classical-vs-Hermite transverse shear visualization.

Reference: approximate Saint-Venant rectangular-beam shear tau_yz(z,y)
computed independently from the applied end force and distributed body force.
FE: unmodified, one-sided source-Hermite constitutive stress component tau_yz,
D @ (B @ (T @ U)), evaluated at visual subcell corners, NEVER fitted to the
reference. The reference is NOT a solution of exact 3D orthotropic elasticity.

Each display cell has its own eight VTK points, even at shared geometry. This
preserves discontinuous, one-sided FE stresses instead of using point averaging
across FE element interfaces. Original source mesh is sampled into rendering
subcells, not replaced by a new finite element solution. Both panels use the
same signed scale in psi, undeformed geometry and linked camera.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from .cubic_stress import analytical_stress_at
from .pyvista_viewer import _pyvista
from core.reference.beam_theory import P
from core.reference.awc import DEFAULT_WEIGHT_DENSITY

FE_SHEAR = "Hermite FE tau_yz (psi)"
BEAM_SHEAR = "Beam theory tau_yz (psi)"
DIFFERENCE = "FE minus beam tau_yz (psi)"


@dataclass
class ShearComparison:
    points: np.ndarray             # duplicated 8-per-cell rendering geometry
    connectivity: np.ndarray       # one independent point set per display hex
    fields: dict[str, np.ndarray]  # per-corner point arrays (signed psi)
    sample_location: str
    stress_evaluation_points: np.ndarray  # may be inset from display geometry


def prepare_shear_comparison(sample, beam, *, end_force=P, density=DEFAULT_WEIGHT_DENSITY):
    """Collect independent FE and classical stresses at IDENTICAL stress points.

    The FE values are sourced from sample.stress_corner_values, which are
    computed by original Hermite B matrices. No smoothing/fitting/constitutive
    substitution happens in this function.
    """
    n = len(sample.connectivity)
    if sample.stress_corner_values.shape != (n, 8, 6):
        raise ValueError("Sample must retain 8 one-sided, six-component stress values per cell")
    stress_positions = np.asarray(sample.stress_corner_points, dtype=float).reshape(-1, 3)
    xyz = np.asarray(sample.points, dtype=float)[sample.connectivity].reshape(-1, 3)
    fe = np.asarray(sample.stress_corner_values, dtype=float)[:, :, 4].ravel().copy()
    _, theory = analytical_stress_at(stress_positions, beam, end_force=end_force, density=density)
    if not np.all(np.isfinite(fe)) or not np.all(np.isfinite(theory)):
        raise ArithmeticError("Non-finite FE or theoretical tau_yz")
    # No point sharing: neighboring one-sided Hermite stresses remain separate.
    return ShearComparison(
        xyz.copy(), np.arange(n * 8, dtype=np.int64).reshape(n, 8),
        {FE_SHEAR: fe, BEAM_SHEAR: theory.copy(), DIFFERENCE: fe - theory},
        sample.stress_sample_location, stress_positions.copy(),
    )


def shared_signed_limits(data):
    """Same zero-centered limits for both panels (never independently scaled)."""
    raw = max(float(np.max(np.abs(data.fields[FE_SHEAR]))),
              float(np.max(np.abs(data.fields[BEAM_SHEAR]))))
    if not np.isfinite(raw):
        raise ValueError("Non-finite comparison data")
    return (-max(raw, 1e-12), max(raw, 1e-12))


def as_shear_pyvista_grid(data):
    """Construct an intentionally discontinuous point-stress VTK hex volume."""
    pv = _pyvista()
    n = len(data.connectivity)
    if not n or data.points.shape != (8*n, 3):
        raise ValueError("Expected eight independent corner points per display cell")
    cells = np.column_stack((np.full(n, 8, dtype=np.int64), data.connectivity)).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(n, pv.CellType.HEXAHEDRON, dtype=np.uint8),
                               data.points)
    for name, vals in data.fields.items():
        grid.point_data[name] = np.asarray(vals, dtype=float)
    grid.field_data["Shear comparison stress-source declaration"] = np.array([
        "FE: Hermite D @ B @ (T @ U); reference: rectangular beam theory;"])
    grid.field_data["Stress sample location"] = np.array([data.sample_location])
    grid.field_data["Continuous nodal stress averaging"] = np.array([0], dtype=np.int64)
    grid.field_data["Certified 3D error bound"] = np.array([0], dtype=np.int64)
    return grid


def export_shear_vtu(data, path):
    """VTU export of three separate per-corner fields, even without PyVista."""
    dest = Path(path)
    if dest.suffix.lower() != '.vtu':
        raise ValueError("Shear comparison export requires a .vtu path")
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        import pyvista  # noqa: F401
    except ImportError:
        try:
            import vtk
            from vtk.util.numpy_support import numpy_to_vtk
        except ImportError as exc:
            raise RuntimeError("Shear VTU export requires PyVista or VTK") from exc
        grid = vtk.vtkUnstructuredGrid()
        pts = vtk.vtkPoints()
        pts.SetData(numpy_to_vtk(np.ascontiguousarray(data.points), deep=True))
        grid.SetPoints(pts)
        for conn in data.connectivity:
            cell = vtk.vtkHexahedron()
            for k, ix in enumerate(conn):
                cell.GetPointIds().SetId(k, int(ix))
            grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        for name, values in data.fields.items():
            arr = numpy_to_vtk(np.ascontiguousarray(values), deep=True)
            arr.SetName(name)
            grid.GetPointData().AddArray(arr)
        for name, value in (("Continuous nodal stress averaging", 0),
                            ("Certified 3D error bound", 0)):
            arr = numpy_to_vtk(np.asarray([value], dtype=np.int64), deep=True)
            arr.SetName(name)
            grid.GetFieldData().AddArray(arr)
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(str(dest))
        writer.SetInputData(grid)
        if writer.Write() != 1:
            raise IOError(f"Could not write {dest}")
    else:
        as_shear_pyvista_grid(data).save(dest)
    return dest


def show_shear_grid(grid, *, title="Cubic recovered-rotation Hermite",
                    full_screen=False, off_screen=False, screenshot=None,
                    show_edges=True):
    """Synchronized two-pane signed tau_yz with exactly the same color scale."""
    pv = _pyvista()
    for key in (BEAM_SHEAR, FE_SHEAR):
        if key not in grid.point_data:
            raise ValueError(f"Missing comparison field {key!r}; regenerate shear comparison VTU")
    fields = {key: np.asarray(grid.point_data[key], float) for key in (BEAM_SHEAR, FE_SHEAR)}
    if any(not np.all(np.isfinite(v)) for v in fields.values()):
        raise ValueError("Nonfinite shear comparison values")
    m = max((float(np.max(np.abs(v))) for v in fields.values()), default=0.)
    clim = (-max(m, 1e-12), max(m, 1e-12))
    if screenshot is not None:
        screenshot = Path(screenshot)
        screenshot.parent.mkdir(parents=True, exist_ok=True)
    p = pv.Plotter(shape=(1, 2), window_size=(1900, 850), off_screen=off_screen)
    p.set_background('#0c1726', all_renderers=True)
    labels = ((BEAM_SHEAR, "CLASSICAL BEAM REFERENCE  /  tau_yz (psi)"),
              (FE_SHEAR, "HERMITE SOURCE-FE STRESS  /  tau_yz (psi)"))
    pts = np.asarray(grid.points)
    length = float(np.ptp(pts[:, 2]))
    target = np.mean(pts, axis=0)
    camera = [(target[0]-1.8*length,target[1]+.24*length,target[2]+.12*length),
              tuple(target), (0, 1, 0)]
    for slot, (field, label) in enumerate(labels):
        p.subplot(0, slot)
        p.add_mesh(grid, scalars=field, preference='point', cmap='coolwarm',
                   clim=clim, show_edges=show_edges, edge_color='#273544',
                   line_width=0.25, lighting=False,
                   scalar_bar_args={'title': 'Signed shear tau_yz (psi)',
                                    'title_font_size': 10, 'label_font_size': 10,
                                    'color': 'white', 'fmt': '%.2f'})
        p.add_text(label, position='upper_left', font_size=11, color='white')
        p.add_axes()
        p.camera_position = camera
    p.link_views()
    p.subplot(0, 0)
    p.add_text('Research comparison: beam theory is NOT exact 3D elasticity; FE remains unmodified.',
               position='lower_left', font_size=9, color='#e8e8e8')
    p.show(full_screen=full_screen and not off_screen,
           screenshot=str(screenshot) if screenshot else None, auto_close=True)
    return screenshot


def main(argv=None):
    """Reopen saved comparison VTU without solving the beam again."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vtu', type=Path, help='Dual shear VTU generated with --shear-comparison')
    parser.add_argument('--off-screen', action='store_true')
    parser.add_argument('--screenshot', type=Path)
    parser.add_argument('--full-screen', action='store_true')
    parser.add_argument('--no-edges', action='store_true')
    args = parser.parse_args(argv)
    pv = _pyvista()
    grid = pv.read(str(args.vtu))
    show_shear_grid(grid, title=args.vtu.stem, off_screen=args.off_screen,
                    screenshot=args.screenshot, full_screen=args.full_screen,
                    show_edges=not args.no_edges)


if __name__ == '__main__':
    main()
