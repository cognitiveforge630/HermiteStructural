"""Optional modest-mesh legacy COO vs OOP batched CSR memory/time comparison.

Command: python -m tools.benchmark_nodal_assembly --order 3
No sparse factorization or structural solve. Separate Python processes are used
so the legacy allocation doesn't contaminate the OOP process's peak reading.

Peak tracemalloc is a Python allocator indicator, NOT a guaranteed complete
accounting of SciPy native allocations or the Windows process's peak RAM.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import tracemalloc

import numpy as np
from scipy.sparse import coo_matrix

from core.elements.nodal_hex import element_type
from core.materials.redwood import redwood_orthotropic_material
from core.mesh.nodal_hex import StructuredNodalHexMesh
from core.reference.beam_theory import B, H, L, E, P
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from core.solvers.nodal_assembly import GlobalAssembler


def assemble(order, nx, ny, nz, method):
    mesh = StructuredNodalHexMesh(order, (nx-1, ny-1, nz-1), (B, H, L))
    D = redwood_orthotropic_material(E).elastic_matrix()
    template = element_type(order)(mesh.element_extent, D)
    ke = template.get_stiffness_matrix()
    fg = template.integrate_body_force((0., -DEFAULT_WEIGHT_DENSITY, 0.))
    template.clear_force_vector()
    ft = template.facets['+Z'].integrate_traction((0., P/(B*H), 0.))
    if method == 'oop':
        assembler = GlobalAssembler(mesh.global_dof_count, batch_elements=8)
        for (_, _, iz), element in mesh.iter_elements(D, shared_stiffness=ke):
            element.add_force_vector(fg)
            if iz == mesh.divisions[2] - 1:
                # True facet-owned integral, not an ad hoc y DOF stride.
                element.facets['+Z'].integrate_traction((0., P/(B*H), 0.))
            assembler.add_element(element)
        return assembler.finalize()
    if method == 'legacy':
        rows, cols, vals = [], [], []
        f = np.zeros(mesh.global_dof_count)
        for (_, _, iz), element in mesh.iter_elements(D):
            ids = element.dof_indices
            rows.append(np.repeat(ids, len(ids)))
            cols.append(np.tile(ids, len(ids)))
            vals.append(ke.ravel())
            f[ids] += fg + (ft if iz == mesh.divisions[2] - 1 else 0)
        K = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(mesh.global_dof_count,)*2).tocsr()
        K.sum_duplicates()
        return K, f
    raise ValueError(method)


def single_process(args):
    tracemalloc.start()
    began = time.perf_counter()
    K, F = assemble(args.order, args.nx, args.ny, args.nz, args.method)
    seconds = time.perf_counter()-began
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    # Native process RSS where a portable function exists; no Windows
    # process-RSS claim without a separately validated Windows measurement.
    peak_rss_kib = None
    if sys.platform.startswith('linux'):
        import resource
        peak_rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({
        'method':args.method,
        'seconds':seconds,
        'tracemalloc_peak_MiB':peak/2**20,
        'linux_process_peak_RSS_MiB':None if peak_rss_kib is None else peak_rss_kib/1024,
        'matrix_shape':K.shape,
        'stored_nnz':int(K.nnz),
        'numerical_nonzero_count':int(K.count_nonzero()),
        'sum_global_F':F.reshape(-1,3).sum(axis=0).tolist(),
    }))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', choices=(2,3), type=int, default=3)
    parser.add_argument('--nx', type=int, default=3)
    parser.add_argument('--ny', type=int, default=3)
    parser.add_argument('--nz', type=int, default=5)
    parser.add_argument('--method', choices=('legacy','oop'), default=None,
                        help='Internal single-process mode; default compares both')
    parser.add_argument('--accept-expensive-run', action='store_true')
    args=parser.parse_args(argv)
    mesh=StructuredNodalHexMesh(args.order, (args.nx-1,args.ny-1,args.nz-1), (B,H,L))
    if mesh.global_dof_count>8000 and not args.accept_expensive_run:
        parser.error('Large K may exhaust RAM: --accept-expensive-run required')
    if args.method:
        single_process(args)
        return
    rows = []
    for method in ('legacy','oop'):
        cmd = [sys.executable, '-m', 'tools.benchmark_nodal_assembly',
               '--order', str(args.order), '--nx', str(args.nx), '--ny', str(args.ny),
               '--nz', str(args.nz), '--method', method]
        if args.accept_expensive_run:
            cmd.append('--accept-expensive-run')
        output=subprocess.run(cmd, check=True, text=True, capture_output=True)
        rows.append(json.loads(output.stdout))
    print(json.dumps({'mesh':[args.nx,args.ny,args.nz], 'order':args.order,
                      'results':rows,
                      'caveat':'Do not compare nnz alone; COO/CSR may store different zero entries. '
                               'RSS is Linux-only here; tracemalloc is not a full RAM accounting. '
                               'Results depend on environment and do not include the linear solve.'}, indent=2))


if __name__=='__main__':
    main()
