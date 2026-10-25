"""Prepare reproducible user-run studies. Execution requires explicit --run.

Each subprocess has a user-adjustable timeout and its own log. Successful case
archives are saved by the underlying driver; completed steps are checkpointed.
Large reference/refinement studies are deliberately not run by default.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def archive(folder,mesh,form,clamp='displacement_only',graded=False):
    return str(folder/f'cantilever_{mesh}_{form}_{clamp}{"_graded" if graded else ""}.npz')


def plan(stage,root):
    base=root/stage
    solutions=base/'solutions'
    steps=[]
    def add(label,*args):
        steps.append({'label':label,'argv':['python',*map(str,args)]})
    def solve(label,meshes,forms,*extra):
        add(label,'run_session6_accuracy.py','--problem','cantilever','--meshes',*meshes,
            '--formulations',*forms,'--clamps','displacement_only',*extra,'--out',solutions)
    def matrix(label,files,reference):
        add(label,'run_session7_matrix.py','--solutions',*files,'--reference',reference,
            '--out',base/'matrix')
    quad,cubic,q2,q1='quadratic_rotation_hermite','cubic_rotation_hermite','cauchy_hex27_reference','cauchy_hex8_reference'
    # Reuse the finest completed reference-stage candidate when available.
    # This remains an unqualified numerical reference, never exact truth.
    preferred_reference=Path(archive(root/'reference'/'solutions','4x4x12',q2))
    if stage=='quick':
        add('fixed_mesh_formulation_and_clamp_matrix','run_session6_accuracy.py','--problem','cantilever',
            '--meshes','2x2x4','--formulations',quad,cubic,q2,q1,'--out',solutions)
        files=[archive(solutions,'2x2x4',f,c) for f,c in
               ((quad,'displacement_only'),(quad,'historical_zero_spin'),
                (cubic,'displacement_only'),(cubic,'historical_zero_spin'),(q1,'displacement_only'))]
        matrix('common_reference_comparison',files,archive(solutions,'2x2x4',q2))
        for form in (quad,cubic):
            add(f'equilibrated_{form}','run_session7_equilibrated.py',archive(solutions,'2x2x4',form),
                '--order','6','--out',base/f'{form}_equilibrated.json')
    elif stage=='reference':
        meshes=['2x2x4','3x3x8','4x4x12']
        solve('higher_order_reference_refinement',meshes,[q2])
        for a,b in zip(meshes[:-1],meshes[1:]):
            add(f'reference_change_{a}_to_{b}','compare_session6_fields.py',archive(solutions,a,q2),
                archive(solutions,b,q2),'--out',base/f'reference_{a}_to_{b}.json')
    elif stage=='refinement':
        meshes=['2x2x4','3x3x6','4x4x8']
        solve('all_direction_refinement',meshes,[quad,cubic,q2])
        matrix('one_common_reference',
            [archive(solutions,m,f) for m in meshes for f in (quad,cubic)],
            str(preferred_reference) if preferred_reference.is_file() else archive(solutions,meshes[-1],q2))
        for mesh in meshes:
            add(f'equilibrated_cubic_{mesh}','run_session7_equilibrated.py',archive(solutions,mesh,cubic),
                '--order','6','--out',base/f'equilibrated_cubic_{mesh}.json')
    elif stage=='graded':
        solve('uniform_control',['2x2x8'],[cubic,q2])
        solve('graded_same_element_count',['2x2x8'],[cubic],
              '--z-breaks','0','4','12','28','48','72','96','112','120')
        matrix('graded_vs_uniform_with_same_reference',
               [archive(solutions,'2x2x8',cubic),archive(solutions,'2x2x8',cubic,graded=True)],
               str(preferred_reference) if preferred_reference.is_file() else archive(solutions,'2x2x8',q2))
        for name,graded in (('uniform',False),('graded',True)):
            add(f'equilibrated_{name}','run_session7_equilibrated.py',
                archive(solutions,'2x2x8',cubic,graded=graded),'--order','6',
                '--out',base/f'equilibrated_{name}.json')
    elif stage=='manufactured':
        add('nonpolynomial_exact_sine_refinement','run_session6_accuracy.py','--problem','sine',
            '--meshes','2x2x4','3x3x6','4x4x8','--formulations',quad,cubic,q2,
            '--clamps','displacement_only','--quadrature-order','6','--out',solutions)
        add('regional_true_error_cubic','report_session6_regions.py',
            solutions/'sine_4x4x8_cubic_rotation_hermite_displacement_only.npz',
            '--exclude-clamp-lengths','8','16','24','--out',base/'regional_true_error_cubic.json')
    elif stage=='estimator':
        completed_fine=Path(archive(root/'refinement'/'solutions','4x4x8',cubic))
        if completed_fine.is_file():
            estimator_solution=str(completed_fine)
        else:
            solve('estimator_control_solution',['2x2x4'],[cubic])
            estimator_solution=archive(solutions,'2x2x4',cubic)
        for order in (4,6,8):
            add(f'quadrature_{order}','run_session7_equilibrated.py',estimator_solution,
                '--order',order,'--out',base/f'gap_order_{order}.json')
        add('stress_space_enrichment','run_session7_equilibrated.py',estimator_solution,
            '--order','8','--transverse-degree','5','--span-degree','5',
            '--out',base/'gap_enriched_modes.json')
        add('known_nonzero_shear_error','run_session6_accuracy.py','--problem','shear','--meshes','2x2x2',
            '--formulations',quad,'--clamps','historical_zero_spin','--out',solutions)
        shear=str(solutions/'shear_2x2x2_quadratic_rotation_hermite_historical_zero_spin.npz')
        add('shear_effectivity','run_session7_equilibrated.py',shear,'--order','6',
            '--out',base/'shear_effectivity.json')
    return steps


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['quick','reference','refinement','graded','manufactured','estimator'],default='quick')
    parser.add_argument('--out',type=Path,default=Path('outputs/session7_user'))
    parser.add_argument('--run',action='store_true',help='Execute the prepared commands on this machine')
    parser.add_argument('--linear-solver',choices=['scipy','pypardiso'],default='scipy')
    parser.add_argument('--blas-threads',type=int,help='OpenBLAS threads per child process')
    parser.add_argument('--mkl-threads',type=int,help='MKL/PyPardiso threads per child process')
    parser.add_argument('--timeout-seconds',type=float,default=1800.)
    args=parser.parse_args()
    if args.timeout_seconds<=0:
        parser.error('Timeout must be positive')
    if any(v is not None and v<1 for v in (args.blas_threads,args.mkl_threads)):
        parser.error('Thread counts must be positive')
    steps=plan(args.stage,args.out)
    for step in steps:
        if step['argv'][1]=='run_session6_accuracy.py':
            step['argv']+=['--linear-solver',args.linear_solver]
    child_environment=os.environ.copy()
    if args.blas_threads is not None:
        child_environment['OPENBLAS_NUM_THREADS']=str(args.blas_threads)
    if args.mkl_threads is not None:
        child_environment['MKL_NUM_THREADS']=str(args.mkl_threads)
    folder=args.out/args.stage
    folder.mkdir(parents=True,exist_ok=True)
    (folder/'manifest.json').write_text(json.dumps(steps,indent=2),encoding='utf-8')
    (folder/'execution_environment.json').write_text(json.dumps({
        'linear_solver':args.linear_solver,'thread_environment':{key:child_environment.get(key)
        for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')}},indent=2),encoding='utf-8')
    for step in steps:
        print(subprocess.list2cmdline(step['argv']),flush=True)
    if not args.run:
        print('Manifest prepared. Add --run to execute; no solver was launched.')
        return
    status=[]
    for index,step in enumerate(steps):
        log=folder/f'{index+1:02d}_{step["label"]}.log'
        started=time.perf_counter()
        record={'label':step['label'],'log':str(log),'status':'RUNNING'}
        status.append(record)
        (folder/'run_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
        try:
            with log.open('w',encoding='utf-8') as fp:
                completed=subprocess.run([sys.executable,*step['argv'][1:]],
                    cwd=Path(__file__).resolve().parent,stdout=fp,stderr=subprocess.STDOUT,
                    timeout=args.timeout_seconds,check=False,env=child_environment)
            record['exit_code']=completed.returncode
            record['status']='PASS' if completed.returncode==0 else 'FAILED'
        except subprocess.TimeoutExpired:
            record['status']='TIMED_OUT'
        except KeyboardInterrupt:
            record['status']='INTERRUPTED'
        finally:
            record['seconds']=time.perf_counter()-started
            (folder/'run_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
        if record['status']!='PASS':
            raise SystemExit(f'{record["status"]}: see {log}')
        print(f'Completed {step["label"]} ({record["seconds"]:.1f}s)',flush=True)


if __name__=='__main__':
    main()
