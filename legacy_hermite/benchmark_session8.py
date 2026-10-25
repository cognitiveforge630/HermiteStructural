"""Short optimization benchmark; no large refinement study.

Run before/after changes with distinct labels. Timing samples are illustrative,
not a statistical scaling study. Existing solution archives are read only.
"""
import argparse
import cProfile
import io
import json
import os
from pathlib import Path
import pstats
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label',default='current')
    parser.add_argument('--backend',choices=['scipy','pypardiso'],default='scipy')
    parser.add_argument('--profile',action='store_true')
    parser.add_argument('--blas-threads',type=int)
    parser.add_argument('--mkl-threads',type=int)
    parser.add_argument('--compare',type=Path,help='Verify invariants against an earlier benchmark JSON')
    parser.add_argument('--out',type=Path,default=Path('outputs/session8_benchmarks'))
    args=parser.parse_args()
    if any(v is not None and v<1 for v in (args.blas_threads,args.mkl_threads)):
        parser.error('Thread counts must be positive')
    if args.blas_threads is not None:
        os.environ['OPENBLAS_NUM_THREADS']=str(args.blas_threads)
    if args.mkl_threads is not None:
        os.environ['MKL_NUM_THREADS']=str(args.mkl_threads)
    import numpy as np
    from run_session6_accuracy import run_case
    from compare_session6_fields import compare_fields
    from run_session7_equilibrated import evaluate_archive
    args.out.mkdir(parents=True,exist_ok=True)
    profiler=cProfile.Profile() if args.profile else None
    if profiler:
        profiler.enable()
    report={'label':args.label,'backend':args.backend,'timings_seconds':{}}
    report['thread_environment']={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}
    def timed(name,callback):
        start=time.perf_counter()
        result=callback()
        report['timings_seconds'][name]=time.perf_counter()-start
        return result
    kwargs={'linear_solver':args.backend} if args.backend!='scipy' else {}
    hermite=timed('cantilever_case',lambda:run_case('cantilever',(2,2,4),
        'cubic_rotation_hermite','displacement_only',order=6,**kwargs))
    reference=timed('hex27_case',lambda:run_case('cantilever',(2,2,4),
        'cauchy_hex27_reference','displacement_only',order=6,**kwargs))
    hrow,hbeam,hu,_=hermite
    rrow,rbeam,ru,_=reference
    difference=timed('field_comparison',lambda:compare_fields((hbeam,hu,hrow),(rbeam,ru,rrow),order=5,
        exclude_clamp_lengths=[8.,112.]))
    sine=timed('sine_case',lambda:run_case('sine',(2,2,4),'cubic_rotation_hermite',
        'displacement_only',order=6,**kwargs))
    path=Path('outputs/session7_user/quick/solutions/cantilever_2x2x4_cubic_rotation_hermite_displacement_only.npz')
    estimate=timed('enriched_estimator',lambda:evaluate_archive(path,order=6,transverse_degree=5,span_degree=5))
    report.update(cantilever_mean_uy=hrow['end_face_average_uy_in'],
        sine_true_error_percent=sine[0]['true_energy_error_percent'],
        reference_difference_percent=difference['reference_relative_difference_percent'],
        estimator_percent=estimate['estimated_energy_error_percent_computed_normalization'],
        cantilever_phase_timings=hrow.get('phase_timings_seconds'),
        assembly_statistics=hrow.get('assembly_statistics'),
        linear_solve=hrow.get('linear_solve'),estimator_phase_timings=estimate.get('phase_timings_seconds'))
    if profiler:
        profiler.disable()
        profiler.dump_stats(str(args.out/f'{args.label}.prof'))
        output=io.StringIO()
        pstats.Stats(profiler,stream=output).sort_stats('cumtime').print_stats(25)
        (args.out/f'{args.label}_profile.txt').write_text(output.getvalue(),encoding='utf-8')
    report['total_timed_seconds']=sum(report['timings_seconds'].values())
    if args.compare:
        previous=json.loads(args.compare.read_text(encoding='utf-8'))
        for key in ('cantilever_mean_uy','sine_true_error_percent','reference_difference_percent','estimator_percent'):
            np.testing.assert_allclose(report[key],previous[key],rtol=1e-8,atol=1e-9,err_msg=key)
        report['numerical_comparison']='PASS'
        report['compared_with']=str(args.compare)
    (args.out/f'{args.label}.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
