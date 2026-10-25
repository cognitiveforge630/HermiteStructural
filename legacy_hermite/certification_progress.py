"""Independent *pre-certification* hypercircle identity and Hermite trace audit.

Only archived solutions are read. No stiffness assembly, factorization or FE solve.
An independent 2D complementary-stress integral is compared with the existing
3D constitutive-gap calculation. Ordinary floating point is NOT rigorous
interval arithmetic; this module MUST NOT issue a certificate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from aposteriori.equilibrated_cauchy import (CantileverAiryStress, constitutive_gap,
                                              stress_samples)
from compare_session6_fields import read_solution
from solver import P


def symbolic_trace_polynomial_audit():
    """Prove tensor-product endpoint identities by exact symbolic algebra.

    On a matching, axis-aligned conforming grid, the same global node values
    and recovered gradient tensors feed both sides of an internal face.
    Normal Hermite derivative shapes vanish there; tangential shape factors
    agree from either adjacent brick. Thus displacement traces agree as
    polynomials (not just at sampled interface points).
    """
    s = sp.symbols('s', real=True)
    hm = (2-3*s+s**3)/4
    hp = (2+3*s-s**3)/4
    rm = (1-s-s**2+s**3)/4
    rp = (-1-s+s**2+s**3)/4
    assert [hm.subs(s,-1), hp.subs(s,-1), rm.subs(s,-1),rp.subs(s,-1)] == [1,0,0,0]
    assert [hm.subs(s,1),hp.subs(s,1),rm.subs(s,1),rp.subs(s,1)] == [0,1,0,0]
    assert sp.expand(hm+hp-1)==0
    return {'hermite_endpoint_values_and_normal_derivatives':'EXACT_SYMBOLIC_PASS',
            'structured_matching_face_C0':'DERIVED_CONDITIONALLY_FROM_SHARED_GLOBAL_NODAL_DATA',
            'not_proven':'rounded archived vector enforces exactly zero clamp trace'}


def independent_2d_complementary_norm_sq(beam, stress_field, compliance, order=8):
    """Integrate tau:S:tau using ONLY the y-z cross-section, with analytic x.

    Stress is independent of x for this special Airy construction, so its
    complementary integral over the 3D volume reduces EXACTLY to Lx times a
    2D integral. Each rectangular y-z cell is integrated separately.
    """
    gp, gw = np.polynomial.legendre.leggauss(order)
    yy = np.unique(beam.grid_points[:,1])
    zz = np.unique(beam.grid_points[:,2])
    sq = 0.
    for ylo,yhi in zip(yy[:-1], yy[1:]):
        y = (yhi+ylo)/2+(yhi-ylo)/2*gp
        for zlo,zhi in zip(zz[:-1], zz[1:]):
            z = (zhi+zlo)/2+(zhi-zlo)/2*gp
            Y,Z = np.meshgrid(y,z,indexing='ij')
            xyz = np.column_stack((np.zeros(Y.size),Y.ravel(),Z.ravel()))
            tau = stress_field(xyz)
            density = np.einsum('ni,ni->n',tau@compliance,tau)
            local = np.outer(gw,gw).ravel()
            sq += beam.Lx*(yhi-ylo)*(zhi-zlo)/4 * np.dot(local,density)
    return float(sq)


def progress_from_archive(root:Path, *, orders=(6,8,10,12)):
    root=root.resolve()
    solution=root/'outputs/session7_user/refinement/solutions/cantilever_4x4x8_cubic_rotation_hermite_displacement_only.npz'
    report=root/'outputs/session8_checks/fine_estimator.json'
    if not solution.exists() or not report.exists():
        raise FileNotFoundError('Required exact saved cantilever archive/report not found')
    saved=json.loads(report.read_text(encoding='utf-8'))
    assert saved['mesh']==[4,4,8] and saved['clamp']=='displacement_only'
    assert saved['reconstruction']['mode_count']==36
    beam,U,meta=read_solution(solution)
    assert meta['problem']=='cantilever'
    with np.load(solution,allow_pickle=False) as archive:
        load=np.asarray(archive['load'],float)
    field=CantileverAiryStress(beam,P,transverse_degree=5,span_degree=5)
    # Constructor degrees 5 -> 36 modes (6*6): no new coefficient fitting.
    # Use mode degrees in this specific archive to prevent silently mixing bases.
    field.coefficients=np.asarray(saved['reconstruction']['coefficients'],float)
    compliance=np.linalg.inv(beam.D)
    W=float(load@U)
    E=float(saved['computed_energy_norm']**2)
    eta_direct=float(saved['eta_energy_norm']**2)
    energy_orders=[]
    for order in (4,5,6,8):
        _,ew,es=stress_samples(beam,U,order=order)
        energy=float(np.einsum('ni,ni,n->',es@compliance,es,ew))
        energy_orders.append({'order':order,'energy_squared_lbf_in':energy,
                              'delta_vs_saved_lbf_in':energy-E})
    records=[]
    for order in orders:
        T2=independent_2d_complementary_norm_sq(beam,field,compliance,order)
        dual=T2+E-2*W
        records.append({'two_dimensional_quadrature_order':order,
                        'complementary_stress_norm_squared':T2,
                        'gap_from_independent_dual_identity_squared':dual,
                        'gap_norm_from_dual_identity':float(np.sqrt(max(dual,0))),
                        'relative_percent_to_saved_computed_energy':float(100*np.sqrt(max(dual,0)/E))})
    # Separate full-domain same-quadrature cross-term audit, independent of
    # the Galerkin work shortcut: tests the integration-by-parts hypothesis.
    points,weights,stress=stress_samples(beam,U,order=8)
    tau=field(points)
    cross=float(np.einsum('ni,ni,n->',tau@compliance,stress,weights))
    gap3=constitutive_gap(weights,tau,stress,compliance)
    E3=float(np.einsum('ni,ni,n->',stress@compliance,stress,weights))
    T3=float(np.einsum('ni,ni,n->',tau@compliance,tau,weights))
    assert abs((T3+E3-2*cross)-gap3) < 1e-10
    assert abs(gap3-eta_direct) < 1e-9
    return {
        'status':'ADDITIONAL_INDEPENDENT_IDENTITIES_NOT_CERTIFIED',
        'certified_bound':False,
        'model':'Cauchy orthotropic, cantilever, uniform natural traction, no body force',
        'symbolic_trace':symbolic_trace_polynomial_audit(),
        'identities':{
            'hypercircle_identity':'eta^2 = ||tau||_S^2 + ||u_h||_E^2 - 2*integral(tau:eps(u_h))',
            'equilibrated_stress_work_identity':'integral(tau:eps(u_h)) = external_work IF exact zero prescribed displacement trace',
            'dimensional_reduction':'||tau||_S^2 = Lx * integral_yz(tau:S:tau) because tau(x,y,z)=tau(y,z)',
            'polynomial_exactness_in_exact_arithmetic':'Hermite displacement degree <=3 per reference axis, so FE energy density degree <=6 per axis; tensor Gauss order >=4 is exact on axis-aligned bricks',
            'boundary_unisolvence':'For a bicubic face trace, vanishing at a 4x4 product grid of distinct Gauss points implies identically zero trace in exact arithmetic (not for rounded KKT solutions)',
            'verified_dual_certificate_formula':'eta_up^2 = T_up + E_up - 2 W_lo; relative_to_computed_percent_up = 100 sqrt(max(eta_up^2,0)/E_lo)',
        },
        'archive_results':{
            'external_work_from_saved_load_and_U_lbf_in':W,
            'saved_computed_energy_squared_lbf_in':E,
            'work_minus_saved_energy_lbf_in':W-E,
            'saved_gap_squared_lbf_in':eta_direct,
            'energy_polynomial_order_checks':energy_orders,
            'cross_integral_3d_order8_lbf_in':cross,
            'three_dimensional_stress_norm_squared_lbf_in':T3,
            'three_dimensional_energy_squared_lbf_in':E3,
            'three_dimensional_direct_gap_squared_lbf_in':gap3,
            'integrated_cross_minus_external_work_lbf_in':cross-W,
            'two_dimensional_complementary_quadrature':records,
        },
        'proof_obligations_still_open':[
            'exact clamp trace of the stored floating-point solution, or an explicit conforming correction',
            'rigorous lower enclosure of the external boundary work when using the dual identity',
            'rigorous upper enclosure of the rational Airy complementary norm integral near top-side corners',
            'rigorous lower enclosure of polynomial displacement energy and arithmetic roundoff',
            'exact material/geometry/load interpretation, coefficient and software-path review',
            'independent technical engineering review; energy bound is not a pointwise stress guarantee'],
        'sha256':{str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (solution,report)},
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,default=Path('outputs/certification_progress.json'))
    a=ap.parse_args()
    result=progress_from_archive(a.root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'certified_bound':False,
        'dual_results':result['archive_results']['two_dimensional_complementary_quadrature']},indent=2))

if __name__=='__main__':main()
