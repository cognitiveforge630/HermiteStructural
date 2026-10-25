"""Self-contained PyVista viewer for Tutorial 001 Cauchy Hex8 results."""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--screenshot', default=None)
    p.add_argument('--warp', type=float, default=50.)
    p.add_argument('--hud', action='store_true')
    p.add_argument('--hud-font-size', type=int, default=10)
    p.add_argument('--title-font-size', type=int, default=10)
    p.add_argument('--scalar-title-font-size', type=int, default=10)
    p.add_argument('--scalar-label-font-size', type=int, default=10)
    p.add_argument('--window-width', type=int, default=1280)
    p.add_argument('--window-height', type=int, default=800)
    a = p.parse_args()
    import pyvista as pv
    base = Path(__file__).resolve().parent
    path = base / 'outputs' / 'tutorial001_solution.npz'
    if not path.exists():
        p.error('Missing outputs/tutorial001_solution.npz; run python .\\tutorial_001.py first.')
    with np.load(path) as d:
        xyz, cells, u = (d[k] for k in ('xyz', 'cells', 'u'))
    grid = pv.UnstructuredGrid(np.c_[np.full(len(cells), 8), cells].astype(np.int64).ravel(),
                               np.full(len(cells), pv.CellType.HEXAHEDRON, dtype=np.uint8), xyz)
    grid.point_data['displacement_in'] = u
    grid.point_data['displacement_magnitude_in'] = np.linalg.norm(u, axis=1)
    plot = pv.Plotter(off_screen=bool(a.screenshot), window_size=(a.window_width, a.window_height))
    plot.add_mesh(grid.warp_by_vector('displacement_in', factor=a.warp),
                  scalars='displacement_magnitude_in', show_edges=True, cmap='viridis',
                  scalar_bar_args={'title': 'Displacement (in)', 'title_font_size': a.scalar_title_font_size,
                                   'label_font_size': a.scalar_label_font_size})
    plot.add_axes()
    if a.hud:
        result_path = base / 'outputs' / 'tutorial001_result.json'
        if result_path.exists():
            r = json.loads(result_path.read_text())
            fem = float(r['tip_mean_uy_in'])
            ref = float(r['reference_EB_tip_uy_in'])
            error = abs((fem - ref) / ref) * 100 if ref else float('nan')
            label = (f'Cauchy Hex8 (not Hermite) | {r["elements"]} elements\n'
                     f'FEM tip uy: {fem:.8f} in\n'
                     f'Euler-Bernoulli tip uy: {ref:.8f} in\n'
                     f'Difference: {error:.2f}% | Warp: {a.warp:g}x\n'
                     'Engineering reference; not a certified bound')
        else:
            label = 'Cauchy Hex8 | Run tutorial_001.py to generate reference data'
        plot.add_text(label, position='upper_left', font_size=a.hud_font_size, name='hud')
    else:
        plot.add_text('Tutorial 001: Cauchy Hex8 | exaggerated displacement',
                      font_size=a.title_font_size)
    if a.screenshot:
        target = Path(a.screenshot)
        if not target.is_absolute():
            target = base / target
        target.parent.mkdir(parents=True, exist_ok=True)
        plot.show(screenshot=str(target))
        print('Saved', target)
    else:
        plot.show()


if __name__ == '__main__':
    main()
