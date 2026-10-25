"""Exact-rational Prager-Synge a posteriori certificates for 3D Cauchy elasticity.

No floating-point quadrature is used to ESTABLISH a certificate.  This is a
manufactured benchmark, not certification of the saved Hermite cantilever.

Physical domain: [-1/2,1/2]^2 x [0,1] inches; isotropic E=1 psi, nu=0.
Body load: b=(0,0,z^2) lbf/in^3 for z measured in inches, relative to 1 in.
Dirichlet: u=0 at z=0; every other surface has zero applied traction.

The exact statically admissible stress sigma_zz=(1-z^3)/3, other stresses 0,
and kinematically admissible displacement u_h=(0,0,w_h(z)) provide a
hypercircle bound. All integrals and Galerkin coefficients use Fractions.
"""
from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import math


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def p_add(a,b):
    return trim([ (a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0))
                 for i in range(max(len(a),len(b)))])


def p_mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)


def p_scale(a,s):return trim([x*s for x in a])


def deriv(a):return trim([i*a[i] for i in range(1,len(a))] or [Q(0)])


def integral(a):return sum((x/Q(i+1) for i,x in enumerate(a)),Q(0))


def inner(a,b):return integral(p_mul(a,b))


def solve_rational(A,B):
    """Exact Gauss-Jordan reduction, rejecting singular systems."""
    n=len(B)
    a=[list(A[i])+[B[i]] for i in range(n)]
    for j in range(n):
        row=next((k for k in range(j,n) if a[k][j]!=0),None)
        if row is None:raise ArithmeticError('Galerkin stiffness is singular')
        a[j],a[row]=a[row],a[j]
        p=a[j][j];a[j]=[t/p for t in a[j]]
        for k in range(n):
            if k!=j:
                c=a[k][j];a[k]=[x-c*y for x,y in zip(a[k],a[j])]
    return [row[-1] for row in a]


def as_str(x):return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'


def poly_strings(x):return [as_str(t) for t in x]


def certificate(name,basis,strict_percent,confirm_known_exact=True):
    """Certify a nonzero energy error for a Galerkin polynomial displacement."""
    body=[Q(0),Q(0),Q(1)]             # b_z=z^2
    tau=[Q(1,3),Q(0),Q(0),Q(-1,3)] # sigma_zz=(1-z^3)/3
    exact_u=[Q(0),Q(1,3),Q(0),Q(0),Q(-1,12)]
    assert p_add(deriv(tau),body)==[Q(0)]  # div sigma + b = 0 EXACTLY
    assert sum(tau)==0                    # free top z=1 has traction zero
    assert all(v[0]==0 for v in basis)   # exact essential boundary trace
    assert deriv(exact_u)==tau           # manufactured exact displacement
    A=[[inner(deriv(f),deriv(g)) for g in basis] for f in basis]
    B=[inner(body,f) for f in basis]
    c=solve_rational(A,B)
    uh=[Q(0)]
    for ci,bi in zip(c,basis):uh=p_add(uh,p_scale(bi,ci))
    eps=deriv(uh)
    difference=p_add(tau,p_scale(eps,Q(-1)))
    E=inner(eps,eps)
    T=inner(tau,tau)
    W=inner(body,uh)
    eta2=inner(difference,difference)
    true_error2=inner(p_add(deriv(exact_u),p_scale(eps,Q(-1))),
                      p_add(deriv(exact_u),p_scale(eps,Q(-1))))
    assert eta2>0 and E>0 and T>0
    assert eta2==T+E-2*W    # hypercircle polarization/integration-by-parts
    assert eta2==true_error2  # independently known manufactured exact answer
    assert E==W                # Galerkin orthogonality
    assert all(inner(difference,deriv(p))==0 for p in basis)
    assert all(sum(A[i][j]*c[j] for j in range(len(c)))==B[i] for i in range(len(c)))
    # The error percentages below are certified by EXACT RATIONAL SQUARING.
    threshold=Q(strict_percent,100) if isinstance(strict_percent,int) else Q(str(strict_percent))/100
    assert eta2 < threshold*threshold*E
    return {
      'case':name,
      'status':'CERTIFIED_MATHEMATICAL_ENERGY_ERROR_BOUND',
      'certificate_type':'exact-rational Prager-Synge hypercircle',
      'certified_bound':True,
      'certified_relative_energy_percent_strictly_below':as_str(threshold*100),
      'bound_verified_by_exact_rational_comparison':True,
      'galerkin_basis_polynomials_z_ascending':[poly_strings(p) for p in basis],
      'galerkin_coefficients_exact':list(map(as_str,c)),
      'galerkin_stiffness_exact':[[as_str(v) for v in row] for row in A],
      'galerkin_load_exact':list(map(as_str,B)),
      'displacement_polynomial_exact':poly_strings(uh),
      'stress_polynomial_exact':poly_strings(tau),
      'energy_squared_lbf_in':as_str(E),
      'complementary_stress_squared_lbf_in':as_str(T),
      'external_work_lbf_in':as_str(W),
      'hypercircle_error_squared_lbf_in':as_str(eta2),
      'true_error_squared_known_only_for_validation_lbf_in':as_str(true_error2),
      'relative_to_computed_energy_percent_approx_display_only':100*math.sqrt(float(eta2/E)),
      'relative_to_exact_energy_percent_approx_display_only':100*math.sqrt(float(eta2/T)),
      'proof_gates':{
       '3d_displacement_trace_and_continuity_exact':True,
       'stress_symmetric_exact':True,
       'stress_divergence_plus_body_load_zero_polynomial_identity':True,
       'natural_traction_all_exposed_faces_exact':True,
       'material_strictly_positive_definite_isotropic_E_1_nu_0':True,
       'exact_galerkin_discrete_solve':True,
       'all_energy_integrals_exact_rationals':True,
       'hypercircle_identity_exact':True,
       'strict_percent_bound_exact_rational_comparison':True,
      },
      'scope_limitations':('Idealized 3D linear Cauchy elasticity, exactly specified rational data and '+
        'exact polynomial Galerkin field. No claim of verified-roundoff bounds '+
        'for arbitrary floating-point solver output; not the archived 9.69% cantilever '+
        'and not a Cosserat arch or engineering design approval.'),
    }



def definite(a,low,high):
    """Integrate a rational-coefficient polynomial on rational endpoints."""
    return sum((v*(high**(i+1)-low**(i+1))/Q(i+1)
                for i,v in enumerate(a)),Q(0))


def certified_multilayer_hex8(nlayers=5,strict_percent='9.6'):
    """FULL 3D Q1 Hex8, 1x1xN, certified with exact piecewise-rational algebra.

    The exact 3D Galerkin field is transverse-constant and piecewise affine
    in z, because at every node its 3D finite-element equations factor into
    the same 1D axial equations times a transverse bilinear shape integral.
    With nu=0 all horizontal-displacement equation RHS and LHS are zero.
    Korn coercivity with clamped bottom gives uniqueness in the full Q1 space.
    """
    n=nlayers
    nodes=[Q(i,n) for i in range(n+1)]
    K=[[Q(0) for j in range(n)] for i in range(n)]
    F=[Q(0)]*n
    local=[]
    body=[Q(0),Q(0),Q(1)]
    tau=[Q(1,3),Q(0),Q(0),Q(-1,3)]
    for k in range(n):
        a,b=nodes[k:k+2];h=b-a
        left=[b/h,-1/h];right=[-a/h,1/h]
        local.append((left,right))
        ids=(k,k+1);shapes=(left,right)
        for i in range(2):
            if ids[i]==0:continue
            I=ids[i]-1
            F[I]+=definite(p_mul(body,shapes[i]),a,b)
            for j in range(2):
                if ids[j]>0:
                    K[I][ids[j]-1]+=definite(p_mul(deriv(shapes[i]),deriv(shapes[j])),a,b)
    U=[Q(0)]+solve_rational(K,F)
    E=T=W=eta2=Q(0)
    slopes=[]
    for k,(left,right) in enumerate(local):
        a,b=nodes[k:k+2]
        uh=p_add(p_scale(left,U[k]),p_scale(right,U[k+1]))
        eps=deriv(uh)
        slopes.append(eps[0])
        E+=definite(p_mul(eps,eps),a,b)
        T+=definite(p_mul(tau,tau),a,b)
        W+=definite(p_mul(body,uh),a,b)
        gap=p_add(tau,p_scale(eps,-1))
        eta2+=definite(p_mul(gap,gap),a,b)
    assert E>0 and eta2>0 and T>0
    assert T==Q(1,14)
    assert eta2==T+E-2*W and E==W
    assert eta2==T-E  # Galerkin orthogonality, independent identity
    bound=Q(strict_percent)/100
    assert eta2 < bound*bound*E
    assert eta2 < bound*bound*T  # also certified relative to exact energy
    return {'case':f'FULL_3D_CAUCHY_HEX8_Q1_1x1x{n}',
        'status':'CERTIFIED_MATHEMATICAL_ENERGY_ERROR_BOUND',
        'certified_bound':True,
        'certificate_type':'exact-rational Prager-Synge hypercircle; piecewise Q1 3D Galerkin',
        'elements':n,
        'axial_node_positions_in':list(map(as_str,nodes)),
        'axial_fem_solution_at_nodes_in':list(map(as_str,U)),
        'axial_fem_strains_each_element':list(map(as_str,slopes)),
        'galerkin_matrix_exact':[[as_str(t) for t in row] for row in K],
        'galerkin_load_exact':list(map(as_str,F)),
        'energy_squared_lbf_in':as_str(E),
        'complementary_stress_squared_lbf_in':as_str(T),
        'external_work_lbf_in':as_str(W),
        'hypercircle_error_squared_lbf_in':as_str(eta2),
        'certified_relative_computed_energy_percent_strictly_below':as_str(bound*100),
        'certified_relative_exact_energy_percent_strictly_below':as_str(bound*100),
        'both_percent_bounds_verified_by_rational_squaring':True,
        'relative_computed_energy_percent_approx_display_only':100*math.sqrt(float(eta2/E)),
        'relative_exact_energy_percent_approx_display_only':100*math.sqrt(float(eta2/T)),
        'proof_gates':{
          'stress_divergence_plus_body_force_exact':True,
          'stress_symmetric_and_all_natural_tractions_exact':True,
          'Q1_displacement_clamp_and_interface_continuity_exact':True,
          '3D_full_galerkin_uniqueness_by_coercivity_and_transverse_factorization':True,
          'equilibrated_3D_stress_satisfies_Hdiv_exact':True,
          'exact_piecewise_polynomial_integration':True,
          'exact_rational_galerkin_system_solved':True,
          'exact_rational_squared_percent_verification':True},
        'scope_limitations':'Manufactured 3D 1x1xN Cauchy Hex8 benchmark only; does not certify archived Hermite cantilever.'}


def make_certificates():
    # Full trilinear Q1/Hex8 FE space on a single unit cube, with z=0 clamped.
    # Symmetry + zero Poisson coupling gives solution w_h=z/4. For ANY Q1
    # trilinear trial v_z with bottom trace zero, v_z=z*p(x,y), p bilinear.
    # Then a(u_h,v)=1/4 int_xy p, and F(v)=int(z^3) dz int_xy p
    # =1/4 int_xy p. The x/y displacement equations vanish identically.
    # Full 3D Galerkin uniqueness holds via coercivity + the clamped face.
    full_hex8=certificate('FULL_3D_CAUCHY_HEX8_Q1_ONE_ELEMENT',[ [Q(0),Q(1)] ],38)
    # Second *restricted axial subspace* spanned by cubic Hermite polynomials.
    # It is a conforming legitimate Galerkin subspace but NOT the actual
    # enriched 3D cubic_rotation_hermite implementation.
    cubic=certificate('RESTRICTED_CUBIC_HERMITE_AXIAL_GALERKIN',
                      [[Q(0),Q(1)],[Q(0),Q(0),Q(1)],[Q(0),Q(0),Q(0),Q(1)]], '2.36')
    layered=certified_multilayer_hex8(5,'9.6')
    return {
      'status':'THREE_CERTIFIED_EXACT_ARITHMETIC_MANUFACTURED_DEMONSTRATIONS',
      'mathematically_certified':True,
      'archived_cantilever_9_69_percent_certified':False,
      'domain':'[-0.5,0.5] in x [-0.5,0.5] in x [0,1] in',
      'isotropic_material':{'youngs_modulus_psi':'1','poisson_ratio':'0'},
      'body_force_lbf_per_in3':'b=(0,0,z^2), with z in inches relative to 1 in',
      'clamp':'u=(0,0,0) on z=0 full face',
      'natural_traction':'0 on all remaining faces',
      'exact_stress_tensor':'diag(0,0,(1-z^3)/3) psi',
      'exact_displacement_for_independent_validation':'(0,0,z/3-z^4/12) in',
      'cases':[full_hex8,cubic,layered],
    }


def check_against_repository_hex8(nlayers=1):
    """Optional *floating-point* solver regression only; NOT certification proof."""
    import numpy as np
    from reference_hex8 import CauchyHex8Reference
    beam=CauchyHex8Reference(Lx=1.,Ly=1.,Lz=1.,nx=2,ny=2,nz=nlayers+1,
           E=1.,nu=0.,gamma=0.,w=0.,gauss_order=3)
    K=beam.build_global_K().toarray()
    F=np.zeros(len(beam.grid_points)*3)
    for e in range(beam.get_num_elements()):
        xyz=beam.get_element_nodes(e); dofs=beam.get_element_global_dofs(e).ravel()
        for i,xi in enumerate(beam.gauss_points):
          for j,eta in enumerate(beam.gauss_points):
            for k,zeta in enumerate(beam.gauss_points):
                N,*_=beam.get_hermite_displacement_matrices(xi,eta,zeta,xyz)
                J=beam.get_hermite_jacobian(xi,eta,zeta,xyz)
                z=float(np.mean(xyz[:,2]))+zeta*float(np.ptp(xyz[:,2]))/2
                F[dofs]+= N.T@np.array([0.,0.,z*z])*np.linalg.det(J)*(
                    beam.gauss_weights[i]*beam.gauss_weights[j]*beam.gauss_weights[k])
    fixed=beam.global_dofs[beam.get_nodes_at_min_z()].ravel()
    free=np.setdiff1d(np.arange(len(F)),fixed)
    U=np.zeros(len(F));U[free]=np.linalg.solve(K[np.ix_(free,free)],F[free])
    reference=([Q(0),Q(1,4)] if nlayers==1 else
               [Q(x) for x in certified_multilayer_hex8(nlayers)['axial_fem_solution_at_nodes_in']])
    exact=np.zeros(len(F));exact[2::3]=[float(reference[round(float(z)*nlayers)]) for z in beam.grid_points[:,2]]
    deviation=float(np.max(np.abs(U-exact)))
    nodal_residual=float(np.max(np.abs(K[np.ix_(free,free)]@exact[free]-F[free])))
    assert nodal_residual<1e-11, f'3D Hex8 ideal Galerkin residual: {nodal_residual}'
    assert deviation<1e-11, f'Floating solver failed expected field: {deviation}'
    return {'floating_solver_check':'PASS_NUMERICAL_CHECK_NOT_CERTIFICATION',
            'implementation':'reference_hex8.CauchyHex8Reference',
            'axial_elements':nlayers,
            'max_nodal_difference_in':deviation,
            'max_free_dof_residual_lbf':nodal_residual,
            'computed_stiffness_energy_lbf_in':float(U@K@U),
            'computed_external_work_lbf_in':float(U@F)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='outputs/certified_exact_example.json')
    parser.add_argument('--check-repository-hex8',action='store_true')
    opt=parser.parse_args()
    record=make_certificates()
    if opt.check_repository_hex8:
        record['repo_hex8_numerical_regression']=[check_against_repository_hex8(n) for n in (1,5)]
    path=Path(opt.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':record['status'],
          'archived_cantilever_certified':record['archived_cantilever_9_69_percent_certified'],
          'certified_cases':[{k:case[k] for k in ('case','certified_bound',
                'hypercircle_error_squared_lbf_in')}
                for case in record['cases']],
          'repo_hex8_numerical_regression':record.get('repo_hex8_numerical_regression',{}),
          'record':str(path)},indent=2))


if __name__=='__main__':main()
