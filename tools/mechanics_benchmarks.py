"""Fast Beer-style independent mechanics benchmarks for the cubic source formulation.

Examples (from repository root):
    python -m tools.mechanics_benchmarks
    python -m tools.mechanics_benchmarks --out outputs/mechanics_benchmarks.json
No large structural solve, GUI, PyPardiso or calibration is required.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from core.solvers.beam import build_beam
from core.reference.mechanics_of_materials import (
    cantilever_fields, rectangular_stresses, euler_bernoulli_bending_energy,
)
from core.verification.patch_benchmarks import test_modes
from core.verification.stress_tensor import (
    stress_matrix, traction, rotate_stress,
    principal_stresses, strain_matrix, elastic_energy_density, cube_force_and_moment,
)


def build_report():
    b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
    patches=test_modes(b)
    L=b.Lz; E=b.orthotropic_constants.Ez
    I=b.Lx*b.Ly**3/12.
    P=-12.8; w=-.0174*b.Lx*b.Ly
    z=np.linspace(0.,L,17)
    both=cantilever_fields(z,length=L,E=E,I=I,P=P,w=w)
    point=cantilever_fields(z,length=L,E=E,I=I,P=P,w=0.)
    weight=cantilever_fields(z,length=L,E=E,I=I,P=0.,w=w)
    superposition=max(float(np.max(np.abs(both[field]-point[field]-weight[field])))
                      for field in ('uy_in','slope_y_per_z','shear_y_lbf','moment_x_lbf_in'))
    energy=euler_bernoulli_bending_energy(length=L,E=E,I=I,P=P,w=w)
    from numpy.polynomial.legendre import leggauss
    xp,wp=leggauss(5)  # Exactly integrates quartic Euler--Bernoulli displacement.
    zg=.5*L*(xp+1.)
    ug=cantilever_fields(zg,length=L,E=E,I=I,P=P,w=w)['uy_in']
    work=float(P*both['uy_in'][-1] + .5*L*np.dot(wp,w*ug))
    # Analytically integrated 1D consistency, explicitly NOT a FE error estimate.
    integration_work_delta=abs(work-2*energy)
    sigma=np.array((12.,-3.,25.,5.,-7.,4.))
    theta=.37
    R=np.array([[np.cos(theta),-np.sin(theta),0.],
                [np.sin(theta),np.cos(theta),0.],[0.,0.,1.]])
    rotated=rotate_stress(sigma,R)
    stress_rotation_residual=float(np.max(np.abs(stress_matrix(rotated)-R@stress_matrix(sigma)@R.T)))
    principal_rotation_residual=float(np.max(np.abs(principal_stresses(sigma)-principal_stresses(rotated))))
    normal=np.array((1.,2.,3.))/np.sqrt(14.)
    traction_rotation_residual=float(np.max(np.abs(traction(rotated,R@normal)-R@traction(sigma,normal))))
    cube_force,cube_moment=cube_force_and_moment(sigma,edge=0.125)
    eps=np.array((.002,-.001,.003,.004,-.002,.001))
    energy_tensor=0.5*np.sum(stress_matrix(sigma)*strain_matrix(eps))
    energy_voigt=elastic_energy_density(sigma,eps)
    stress_sample=rectangular_stresses(np.array([-b.Ly/2,0.,b.Ly/2]),L/2,
                            width=b.Lx,depth=b.Ly,length=L,P=P,w=w)
    observations={
        'patch_modes':patches,
        'beam_theory':{
            'tip_uy_in':float(both['uy_in'][-1]),
            'clamp_displacement_in':float(both['uy_in'][0]),
            'clamp_slope':float(both['slope_y_per_z'][0]),
            'free_end_shear_lbf':float(both['shear_y_lbf'][-1]),
            'free_end_moment_lbf_in':float(both['moment_x_lbf_in'][-1]),
            'superposition_max_absolute_residual':superposition,
            'one_dimensional_bending_energy_lbf_in':float(energy),
            'one_dimensional_work_minus_2_energy_lbf_in':float(integration_work_delta),
            'sampled_bending_stress_psi':stress_sample['sigma_zz_psi'].tolist(),
            'sampled_parabolic_shear_psi':stress_sample['tau_yz_psi'].tolist(),
        },
        'stress_cube':{
            'rotation_tensor_max_residual_psi':stress_rotation_residual,
            'principal_eigenvalue_rotation_max_residual_psi':principal_rotation_residual,
            'traction_rotation_max_residual_psi':traction_rotation_residual,
            'tensor_vs_voigt_energy_density_residual_psi':float(abs(energy_tensor-energy_voigt)),
            'constant_stress_cube_force_resultant_lbf':cube_force.tolist(),
            'constant_stress_cube_moment_resultant_lbf_in':cube_moment.tolist(),
        },
    }
    passes={
        'source_polynomial_modes':all(v['stress_field_exactly_reproduced_at_test_points'] for v in patches['modes'].values()),
        'beam_superposition':superposition<1e-9,
        'clamped_boundary_reference':both['uy_in'][0]==0 and both['slope_y_per_z'][0]==0,
        'free_end_reference_moment':abs(both['moment_x_lbf_in'][-1])<1e-12,
        'stress_cube_rotation':stress_rotation_residual<1e-10 and principal_rotation_residual<1e-10 and traction_rotation_residual<1e-10,
        'constant_stress_cube_equilibrium':np.linalg.norm(cube_force)<1e-11 and np.linalg.norm(cube_moment)<1e-11,
        'engineering_shear_energy_convention':abs(energy_tensor-energy_voigt)<1e-12,
        'beam_stress_sign_and_surface_shear':stress_sample['sigma_zz_psi'][0]<0 and stress_sample['sigma_zz_psi'][-1]>0 and abs(stress_sample['tau_yz_psi'][0])<1e-12 and abs(stress_sample['tau_yz_psi'][-1])<1e-12,
        'beam_energy_quadrature_diagnostic':integration_work_delta/max(2*energy,1e-30)<1e-12,
    }
    passes={key:bool(value) for key,value in passes.items()}
    return {'status':('PASS_RESEARCH_BENCHMARKS_NOT_CERTIFIED' if all(passes.values())
                      else 'FAIL_RESEARCH_BENCHMARKS_NOT_CERTIFIED'),
            'reference':'Mechanics of Materials, 6th edition (Beer, Johnston, DeWolf, Mazurek), chapters 1, 4-6, 9, 11',
            'units':'in, lbf, psi, lbf*in',
            'checks':passes, 'observations':observations, 'certified_error_bound':False,
            'limitations':('Prescribed polynomial stress patches and elementary 1D reference checks; '
                           'not independent 3D solved-model validation, boundary traction admissibility, '
                           'or a rigorous a posteriori numerical error bound.')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('outputs/mechanics_benchmarks.json'))
    args=parser.parse_args()
    report=build_report()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    for name,passed in report['checks'].items():
        print(f'{"PASS" if passed else "FAIL"}: {name}')
    print(report['status'])
    print(f'Wrote {args.out}')
    if not all(report['checks'].values()):
        raise SystemExit(1)

if __name__=='__main__':
    main()
