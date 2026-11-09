"""Sample actual Hermite displacement and Cauchy stress for volume visualization.

The source element has eight mesh vertices but *nonlinear recovered Hermite*
interpolation. An ordinary VTK Hex8 representation of nodal values would not
reproduce its interior field. Subdivide each cell in reference coordinates and
sample the source element's nonlocal primitive reconstruction at the new points.

The resulting VTK hexes are rendering/export approximations, NOT replacement
FE elements; convergence/certification must be assessed on the original FE
solution. Stress is evaluated at the eight corners of every visualization subcell
(default: exact corners, with optional inset corners). Cell arrays retain
the arithmetic mean of those eight one-sided source-element stress samples
for compatibility with VTU cell coloring. The individual corner samples are
retained separately for stress-error and peak checks: cell means MUST NOT
be reported as point stresses or as extrema. No averaging is done between
distinct source elements at a shared node.
"""
from dataclasses import dataclass
import numpy as np

POINT_FIELDS = ("Displacement magnitude (in)", "Vertical displacement magnitude (in)",
                "Ux (in)", "Uy (in)", "Uz (in)",
                "Fixed support node", "Loaded end node")
CELL_FIELDS = (
    "Von Mises (psi)", "Sigma xx (psi)", "Sigma yy (psi)",
    "Sigma zz (psi)", "Absolute longitudinal stress (psi)", "Tau xy (psi)", "Tau yz (psi)",
    "Tau xz (psi)", "Maximum principal (psi)", "Minimum principal (psi)",
)
DEFAULT_FIELD = "Vertical displacement magnitude (in)"


@dataclass
class SampledVolume:
    points: np.ndarray
    connectivity: np.ndarray       # Ncells x 8, VTK_HEXAHEDRON ordering
    displacement: np.ndarray       # Npoints x 3, exact Hermite evaluation at samples
    point_fields: dict[str, np.ndarray]
    cell_fields: dict[str, np.ndarray]
    subdivisions: int
    element_count: int
    stress_sample_location: str
    stress_inset_fraction: float
    stress_corner_points: np.ndarray   # Ncells x 8 x 3, physical coordinates
    stress_corner_values: np.ndarray   # Ncells x 8 x 6, one-sided Cauchy stress


def sample_volume(beam, U, subdivisions=2, *, max_cells=250000,
                  stress_sample_location="exact_corners", stress_inset_fraction=0.95):
    """Build a connected structured rendering mesh from true Hermite samples.

    This is intentionally independent of PyVista, enabling numerical tests in
    non-graphics environments. Node indices are shared across element faces;
    each sample is evaluated using the source FE element it belongs to.

    Parameters
    ----------
    beam: the solved linear/quadratic/cubic recovered-rotation Hermite element
    U: complete solver DOF vector, alternating 3 translations and 3 rotations
    subdivisions: number of display hexes along each FE element axis, >=1
    max_cells: practical guard against accidentally excessive rendering meshes
    stress_sample_location: 'exact_corners' or 'inset_corners'. Stress is never
        calculated at subcell centroids in this visualization/recovery path.
    stress_inset_fraction: for inset mode, 0 < fraction < 1 scales each
        corner's distance from the subcell center (1 means exact corners).
    """
    if isinstance(subdivisions, bool) or int(subdivisions) != subdivisions or subdivisions < 1:
        raise ValueError("subdivisions must be a positive integer")
    subdivisions = int(subdivisions)
    if stress_sample_location not in ("exact_corners", "inset_corners"):
        raise ValueError('stress_sample_location must be exact_corners or inset_corners')
    if not np.isfinite(stress_inset_fraction) or not 0 < stress_inset_fraction < 1:
        raise ValueError('stress_inset_fraction must be finite and strictly between 0 and 1')
    corner_fraction = 1. if stress_sample_location == 'exact_corners' else float(stress_inset_fraction)
    nx, ny, nz = beam.nx, beam.ny, beam.nz
    n_vis_cells = (nx-1)*(ny-1)*(nz-1)*subdivisions**3
    if n_vis_cells > max_cells:
        raise ValueError(f"Visualization would contain {n_vis_cells:,} cells; lower --subdivisions")
    U = np.asarray(U, dtype=float)
    n_dof = len(beam.grid_points)*beam.ndof_per_node
    if U.shape != (n_dof,) or not np.all(np.isfinite(U)):
        raise ValueError(f"Expected a finite full solver displacement vector of shape ({n_dof},)")

    npx, npy, npz = ((n-1)*subdivisions + 1 for n in (nx,ny,nz))
    dx = np.linspace(-beam.Lx/2, beam.Lx/2, npx)
    dy = np.linspace(-beam.Ly/2, beam.Ly/2, npy)
    dz = np.linspace(0, beam.Lz, npz)
    X, Y, Z = np.meshgrid(dx, dy, dz, indexing="ij")
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    displacement = np.zeros_like(points)
    visited = np.zeros(len(points), dtype=bool)
    connectivity = np.empty((n_vis_cells,8),dtype=np.int64)
    stresses = np.empty((n_vis_cells, 6), dtype=float)
    corner_points = np.empty((n_vis_cells,8,3), dtype=float)
    corner_stresses = np.empty((n_vis_cells,8,6), dtype=float)
    D = beam.elastic_matrix()

    def vid(i,j,k):
        return (i*npy+j)*npz+k

    parents = np.linspace(-1.0,1.0,subdivisions+1)
    cell_id=0
    # Match the FE loop's i,j,k element ordering; reuse each element's
    # nonlocal support transformation for the visualization evaluations.
    for ix in range(nx-1):
        for iy in range(ny-1):
            for iz in range(nz-1):
                e = (ix*(ny-1)+iy)*(nz-1)+iz
                xyz = beam.get_element_nodes(e)
                ids,T = beam.element_transform(e)
                primitive_u = T@U[ids]
                for ai,xi in enumerate(parents):
                    for aj,eta in enumerate(parents):
                        for ak,zeta in enumerate(parents):
                            i,j,k = ix*subdivisions+ai, iy*subdivisions+aj, iz*subdivisions+ak
                            pi=vid(i,j,k)
                            if visited[pi]:
                                continue
                            op,_=beam.primitive_matrices(float(xi),float(eta),float(zeta),xyz)
                            displacement[pi]=op@primitive_u
                            visited[pi]=True
                for ai in range(subdivisions):
                    for aj in range(subdivisions):
                        for ak in range(subdivisions):
                            i,j,k=ix*subdivisions+ai,iy*subdivisions+aj,iz*subdivisions+ak
                            connectivity[cell_id] = (
                                vid(i,j,k),vid(i+1,j,k),vid(i+1,j+1,k),vid(i,j+1,k),
                                vid(i,j,k+1),vid(i+1,j,k+1),vid(i+1,j+1,k+1),vid(i,j+1,k+1))
                            # These are ONE-SIDED source-element stress samples
                            # at subcell corners, not gradient estimates from
                            # the visualization geometry and not cell centroids.
                            bounds = tuple((parents[q],parents[q+1]) for q in (ai,aj,ak))
                            midpoint = np.array([0.5*(lo+hi) for lo,hi in bounds])
                            for idx,(sx,sy,sz) in enumerate(((-1,-1,-1),(1,-1,-1),
                                    (1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),
                                    (1,1,1),(-1,1,1))):
                                parent = midpoint + corner_fraction * np.array((
                                    sx*(bounds[0][1]-bounds[0][0])/2,
                                    sy*(bounds[1][1]-bounds[1][0])/2,
                                    sz*(bounds[2][1]-bounds[2][0])/2))
                                B,detJ,_=beam.primitive_B(*parent,xyz)
                                if not np.isfinite(detJ) or detJ <= 0:
                                    raise ArithmeticError('Invalid stress-sampling Jacobian')
                                corner_stresses[cell_id,idx]=D@(B@primitive_u)
                                corner_points[cell_id,idx]=beam.hex8_shape_functions(*parent)@xyz
                            # Only for cell coloring: the individual corner
                            # values remain available to independent audits.
                            stresses[cell_id]=corner_stresses[cell_id].mean(axis=0)
                            cell_id+=1
    assert visited.all() and cell_id==n_vis_cells
    if not np.all(np.isfinite(displacement)) or not np.all(np.isfinite(stresses)) or not np.all(np.isfinite(corner_stresses)):
        raise ArithmeticError("Nonfinite recovered Hermite displacement or stress")
    sx,sy,sz,txy,tyz,txz=stresses.T
    vm=np.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(txy**2+tyz**2+txz**2))
    tensor=np.empty((n_vis_cells,3,3),dtype=float)
    tensor[:,0,0],tensor[:,1,1],tensor[:,2,2]=sx,sy,sz
    tensor[:,0,1]=tensor[:,1,0]=txy
    tensor[:,1,2]=tensor[:,2,1]=tyz
    tensor[:,0,2]=tensor[:,2,0]=txz
    eig=np.linalg.eigvalsh(tensor)
    # The visualization lattice is subdivided uniformly along each original
    # source FE edge.  Mark only original source nodes, not interpolation nodes.
    # These metadata fields are *not* force/reaction vectors.
    source_x = np.arange(npx) % subdivisions == 0
    source_y = np.arange(npy) % subdivisions == 0
    source_z = np.arange(npz) % subdivisions == 0
    source_grid = (source_x[:,None,None] & source_y[None,:,None] &
                   source_z[None,None,:])
    fixed_flag = (source_grid & (np.arange(npz)[None,None,:] == 0)).ravel()
    loaded_flag = (source_grid & (np.arange(npz)[None,None,:] == npz-1)).ravel()
    point_fields={"Displacement magnitude (in)":np.linalg.norm(displacement,axis=1),
                  "Vertical displacement magnitude (in)":np.abs(displacement[:,1]),
                  "Ux (in)":displacement[:,0],"Uy (in)":displacement[:,1],
                  "Uz (in)":displacement[:,2],
                  "Fixed support node":fixed_flag.astype(np.uint8),
                  "Loaded end node":loaded_flag.astype(np.uint8)}
    cell_fields={"Von Mises (psi)":vm,
                 "Sigma xx (psi)":sx,"Sigma yy (psi)":sy,"Sigma zz (psi)":sz,
                 "Absolute longitudinal stress (psi)":np.abs(sz),
                 "Tau xy (psi)":txy,"Tau yz (psi)":tyz,"Tau xz (psi)":txz,
                 "Maximum principal (psi)":eig[:,2],
                 "Minimum principal (psi)":eig[:,0]}
    return SampledVolume(points, connectivity, displacement,
                         point_fields,cell_fields,subdivisions,beam.get_num_elements(),
                         stress_sample_location,corner_fraction,corner_points,corner_stresses)
