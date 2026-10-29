"""View Hex27/Hex64 *inset-corner sampled* stress; no centroid stress allowed.

IMPORTANT: small 16-element examples are SMOKE TESTS, not stress validation.
Stresses are directly calculated at 8 inset corners per visualization subcell;
each plotted cell uses an actual sampled extreme (never a centroid evaluation).
Do NOT interpret displayed values as exact or certified 3D stress maxima.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from core.visualization.nodal_convergence import STRESS_SAMPLING_SCHEME

STRESS_COMPONENTS = ('Sigma xx (psi)', 'Sigma yy (psi)', 'Sigma zz (psi)',
                     'Tau xy (psi)', 'Tau yz (psi)', 'Tau xz (psi)')
DISPLACEMENT_COMPONENTS = ('Ux (in)', 'Uy (in)', 'Uz (in)')
FIELDS = ('Displacement magnitude (in)', *DISPLACEMENT_COMPONENTS,
          'Tension zz (psi)', 'Compression zz (psi)', 'Von Mises (psi)',
          'Shear magnitude (psi)', *STRESS_COMPONENTS,
          'Beam reference Sigma zz (psi)', 'Beam reference Von Mises (psi)')


def require_inset_archive(data):
    """Make old center-sample NPZ files fail closed, not silently displayable."""
    if ('stress_sampling_scheme' not in data or
            str(data['stress_sampling_scheme']) != STRESS_SAMPLING_SCHEME):
        raise ValueError('OBSOLETE CENTROID-STRESS ARCHIVE: regenerate with '
                         'python -m tools.p_order_study --mode solve ... --field-output ... '
                         'or python -m tools.convergence_study --stage solve. '
                         'Centroid-based stress viewing is prohibited.')


def sampled_extreme(s, axis=-1):
    """Choose the actual signed sample with greatest absolute magnitude."""
    idx = np.argmax(np.abs(s), axis=axis)
    return np.take_along_axis(s, np.expand_dims(idx,axis=axis),axis=axis).squeeze(axis=axis)


def build_grid(path):
    """Build Hex visualization subcells and assign directly sampled stress extrema."""
    try:
        import pyvista as pv
    except ImportError as exc:
        raise RuntimeError('Install PyVista: python -m pip install -r requirements.txt') from exc

    with np.load(path) as archive:
        require_inset_archive(archive)
        order = int(archive['order'])
        div = tuple(int(v) for v in archive['divisions'])
        shape = tuple(order*d+1 for d in div)
        points = np.asarray(archive['coordinates'], dtype=float)
        displacement = np.asarray(archive['displacement'], dtype=float)
        stress = np.asarray(archive['fe_stress'], dtype=float)
        theory = np.asarray(archive['beam_reference_stress'], dtype=float)
        vm = np.asarray(archive['fe_von_mises'], dtype=float)
        rvm = np.asarray(archive['beam_reference_von_mises'], dtype=float)
    nsub=int(np.prod(div))*order**3
    if points.shape != (int(np.prod(shape)), 3) or displacement.shape != points.shape:
        raise ValueError('Invalid nodal displacement coordinates/shape in archive')
    if stress.shape != (nsub,8,6) or theory.shape != stress.shape:
        raise ValueError('Invalid inset-corner stress sample shape in archive')
    if vm.shape != (nsub,8) or rvm.shape != vm.shape:
        raise ValueError('Invalid inset-corner von Mises shape in archive')

    fine = np.arange(int(np.prod(shape))).reshape(shape)
    conn = []
    indices = []
    for i in range(shape[0]-1):
        for j in range(shape[1]-1):
            for k in range(shape[2]-1):
                conn.append((fine[i,j,k], fine[i+1,j,k], fine[i+1,j+1,k], fine[i,j+1,k],
                             fine[i,j,k+1], fine[i+1,j,k+1], fine[i+1,j+1,k+1], fine[i,j+1,k+1]))
                indices.append(np.ravel_multi_index(
                    (i//order,j//order,k//order,i%order,j%order,k%order),
                    (*div,order,order,order)))
    conn = np.asarray(conn, dtype=np.int64)
    cells = np.column_stack((np.full(len(conn),8,dtype=np.int64),conn)).reshape(-1)
    grid = pv.UnstructuredGrid(cells, np.full(len(conn),pv.CellType.HEXAHEDRON,dtype=np.uint8),points)
    grid.point_data['Displacement (in)'] = displacement
    grid.point_data['Displacement magnitude (in)'] = np.linalg.norm(displacement,axis=1)
    for j,name in enumerate(DISPLACEMENT_COMPONENTS):
        grid.point_data[name] = displacement[:,j]
    s = stress[np.asarray(indices,dtype=np.int64)]
    r = theory[np.asarray(indices,dtype=np.int64)]
    for j,name in enumerate(STRESS_COMPONENTS):
        grid.cell_data[name] = sampled_extreme(s[:,:,j])
    grid.cell_data['Tension zz (psi)'] = np.maximum(s[:,:,2],0).max(axis=1)
    grid.cell_data['Compression zz (psi)'] = np.maximum(-s[:,:,2],0).max(axis=1)
    grid.cell_data['Von Mises (psi)'] = vm[indices].max(axis=1)
    grid.cell_data['Shear magnitude (psi)'] = np.linalg.norm(s[:,:,3:6],axis=2).max(axis=1)
    grid.cell_data['Beam reference Sigma zz (psi)'] = sampled_extreme(r[:,:,2])
    grid.cell_data['Beam reference Von Mises (psi)'] = rvm[indices].max(axis=1)
    grid.field_data['Direct inset corner stress samples per subcell'] = np.array([8])
    return grid


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path,help='Current-version NPZ archive from --field-output')
    parser.add_argument('--field',choices=FIELDS,default='Displacement magnitude (in)')
    parser.add_argument('--mesh-only',action='store_true',help='Show undeformed mesh edges without field coloring')
    parser.add_argument('--vtu',type=Path,help='Export all fields to a ParaView/PyVista VTU file')
    parser.add_argument('--off-screen',action='store_true',help='Render without an interactive display')
    parser.add_argument('--screenshot',type=Path,help='Save a PNG snapshot')
    args=parser.parse_args(argv)
    grid=build_grid(args.archive)
    with np.load(args.archive) as archive:
        div=np.asarray(archive['divisions'],dtype=int)
    smoke=int(np.prod(div))<=16
    banner = ('SMOKE TEST ONLY - STRESS NOT VALIDATED - NOT FOR DESIGN' if smoke else
              'RESEARCH ONLY - SAMPLED STRESS; NOT CERTIFIED')
    detail = ('8 DIRECT INSET-CORNER SAMPLES PER SUBCELL - NO CENTROID STRESS')
    print(f'WARNING: {banner}. {detail}. Sampled peaks are not true global maxima.',flush=True)
    if args.vtu:
        args.vtu.parent.mkdir(parents=True,exist_ok=True)
        grid.save(args.vtu)
        print('Exported:',args.vtu)
    if args.mesh_only or not args.vtu or args.screenshot:
        import pyvista as pv
        plotter=pv.Plotter(off_screen=args.off_screen)
        if args.mesh_only:
            plotter.add_mesh(grid,show_edges=True,color='lightsteelblue',opacity=0.75)
        else:
            plotter.add_mesh(grid,scalars=args.field,show_edges=False)
        plotter.add_axes()
        plotter.add_text(banner+'\n'+detail+'\nSampled extrema only; not exact 3D maxima.',font_size=10)
        if args.screenshot:
            args.screenshot.parent.mkdir(parents=True,exist_ok=True)
        plotter.show(screenshot=str(args.screenshot) if args.screenshot else None)
    return grid


if __name__=='__main__':
    main()
