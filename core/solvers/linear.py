"""Explicit optional sparse CPU backends, with independent residual checks."""
import time
import threading
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

_pardiso_state=threading.local()


def solve_sparse_system(matrix,rhs,backend='scipy'):
    if backend not in ('scipy','pypardiso'):
        raise ValueError(f'Unknown linear solver: {backend}')
    A=csr_matrix(matrix,dtype=np.float64,copy=True)
    A.sum_duplicates();A.eliminate_zeros();A.sort_indices()
    b=np.array(rhs,dtype=np.float64,copy=True)
    if A.shape[0]!=A.shape[1] or b.shape!=(A.shape[0],):
        raise ValueError('Expected square matrix and matching RHS vector')
    if not np.all(np.isfinite(A.data)) or not np.all(np.isfinite(b)):
        raise ValueError('Nonfinite linear system')
    if np.any(np.diff(A.indptr)==0):
        raise ValueError('Linear system has an empty row')
    started=time.perf_counter()
    if backend=='scipy':
        solution=spsolve(A,b)
    else:
        try:
            from pypardiso import PyPardisoSolver
        except (ImportError,OSError) as error:
            raise RuntimeError('PyPardiso/MKL unavailable; select --linear-solver scipy') from error
        # Full general storage supports both reduced elasticity matrices and
        # indefinite KKT systems. Never assume the clamp system is SPD.
        # Reuse the library wrapper per calling thread, but release all
        # factorization memory after every solve. No stale matrix factors or
        # mutable PARDISO handles are shared between threads.
        if not hasattr(_pardiso_state,'solver'):
            _pardiso_state.solver=PyPardisoSolver(mtype=11)
        solver=_pardiso_state.solver
        try:
            solution=solver.solve(A,b)
        finally:
            solver.free_memory(everything=True)
    elapsed=time.perf_counter()-started
    if not np.all(np.isfinite(solution)):
        raise ArithmeticError('Linear solve produced nonfinite values')
    residual=float(np.max(np.abs(A @ solution-b)))
    scale=float(np.max(np.asarray(abs(A).sum(axis=1))))*float(np.max(np.abs(solution)))+float(np.max(np.abs(b)))
    backward=residual/scale if scale else residual
    if backward>1e-10:
        raise ArithmeticError(f'{backend} linear backward error {backward:g} exceeds tolerance')
    return solution,{'backend':backend,'seconds':elapsed,'system_size':A.shape[0],
                     'nonzeros':A.nnz,'absolute_residual_inf':residual,
                     'rhs_relative_residual_inf':residual/max(float(np.max(np.abs(b))),1e-30),
                     'relative_backward_error':backward}
