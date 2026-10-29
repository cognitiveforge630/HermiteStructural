"""Strong face displacement constraints for nonlocal Hermite interpolation.

Nodal value and rotation constraints alone do *not* fully constrain the
polynomial trace when recovered derivatives depend on the next element layer.
This research KKT routine imposes face displacement at enough quadrature
points to constrain the full trace polynomial. Linearly dependent constraints
are rank filtered with a pivoted QR, and residuals are verified afterward.
It is intended for modest research meshes, NOT a scalable production solver.
"""
from __future__ import annotations
import numpy as np
import time
from scipy.linalg import qr
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import spsolve
from core.elements.linear import RecoveredRotationHex8
from .linear import solve_sparse_system


def face_trace_constraints(beam, prescribed_translation, order=4):
    if not isinstance(beam, RecoveredRotationHex8):
        raise TypeError('Face trace constraint requires the nonlocal Hermite formulation')
    gp,_=np.polynomial.legendre.leggauss(order)
    total=beam.grid_points.shape[0]*6
    rows=[]; targets=[]
    for e in beam.get_elements_at_min_z():
        xyz=beam.get_element_nodes(e)
        for xi in gp:
            for eta in gp:
                p=(float(xi),float(eta),-1.)
                ids,N=beam.displacement_operator(e,p)
                x=beam.hex8_shape_functions(*p)@xyz
                value=np.asarray(prescribed_translation(x),float)
                if value.shape!=(3,):
                    raise ValueError('Prescribed displacement must have 3 components')
                for c in range(3):
                    row=np.zeros(total)
                    row[ids]=N[c]
                    rows.append(row)
                    targets.append(value[c])
    return np.stack(rows),np.asarray(targets)


def solve_with_face_dirichlet(beam,force,prescribed_translation,
                              prescribed_rotation=None,order=4,*,
                              constrain_rotations=True,stiffness=None,linear_solver='scipy'):
    """KKT solve with a full displacement trace.

    The default also prescribes nodal spins. For classical Cauchy
    displacement Dirichlet data use constrain_rotations=False: zero face
    displacement does not imply zero tangential spin (normal derivatives are
    free). The choice of rotational constraints must be recorded explicitly.
    """
    if not constrain_rotations and prescribed_rotation is not None:
        raise ValueError('Rotation data conflicts with displacement-only constraints')
    started=time.perf_counter()
    K=(beam.build_global_K() if stiffness is None else stiffness).tocsr()
    stiffness_seconds=time.perf_counter()-started
    rhs=np.asarray(force,float)
    if rhs.shape!=(K.shape[0],):
        raise ValueError('Wrong external load shape')
    fixed=[]; known=np.zeros(K.shape[0]);
    for n in beam.get_nodes_at_min_z():
        ids=beam.global_dofs[n]
        fixed.extend((ids if constrain_rotations else ids[:3]).tolist())
        known[ids[:3]]=prescribed_translation(beam.grid_points[n])
        if prescribed_rotation is not None:
            known[ids[3:]]=prescribed_rotation(beam.grid_points[n])
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(len(known)),fixed)
    started=time.perf_counter()
    C,g=face_trace_constraints(beam,prescribed_translation,order)
    A=C[:,free]; target=g-C[:,fixed]@known[fixed]
    trace_seconds=time.perf_counter()-started
    # Select an independent row basis; keep exact target consistency explicit.
    started=time.perf_counter()
    _q,R,piv=qr(A.T,pivoting=True,mode='economic')
    qr_seconds=time.perf_counter()-started
    diag=np.abs(np.diag(R))
    rank=int(np.sum(diag > (diag[0]*1e-10 if len(diag) else 0)))
    select=piv[:rank]
    Af=A[select];bf=target[select]
    cf=csr_matrix(Af)
    system=bmat([[K[free][:,free], cf.T],[cf,csr_matrix((rank,rank))]],format='csr')
    vector=np.r_[rhs[free]-K[free][:,fixed]@known[fixed], bf]
    solution,solve_info=solve_sparse_system(system,vector,linear_solver)
    beam.last_linear_solve=solve_info
    known[free]=solution[:len(free)]
    max_trace=float(np.max(np.abs(C@known-g)))
    scale=max(1.,float(np.max(np.abs(g))))
    if not np.all(np.isfinite(known)) or max_trace>1e-8*scale:
        raise ArithmeticError(f'Face trace constraints unsatisfied: {max_trace}')
    return known,{'constrain_rotations':bool(constrain_rotations),
                  'independent_trace_constraints':rank,
                  'sampled_trace_constraints':len(g),
                  'maximum_trace_error_in':max_trace,
                  'linear_solve':solve_info,
                  'phase_timings_seconds':{'stiffness_assembly':stiffness_seconds,
                      'trace_constraints':trace_seconds,'rank_filter_qr':qr_seconds,
                      'linear_solve':solve_info['seconds']},
                  'kkt_size':system.shape[0]}
