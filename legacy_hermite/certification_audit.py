"""Independent, fast certification-readiness audit (no new FEM solves).

This audit proves symbolic interior/face identities for the *specified* Airy
family and independently checks arithmetic in saved reports. It is NOT an
interval-arithmetic enclosure, a function-space proof, or an engineering
certificate. The output deliberately fails closed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp


def symbolic_airy_audit():
    t,q,r=sp.symbols('t q r', real=True)
    T,H,L=sp.symbols('T H L', positive=True)
    F=t-sp.Rational(1,2)-t*(1-t)**2*q**2/(t**2+q**2)+(1-t)*t**2*q**2/((1-t)**2+q**2)
    # Dimensionless stresses: [yy, zz, yz]/T.
    yy=2*sp.diff(F,q)+q*sp.diff(F,q,2)
    zz=q*sp.diff(F,t,2)
    yz=sp.diff(F,t)+q*sp.diff(F,t,q)
    # Physical divergence: d/dy=(1/H) d/dt, d/dz=-(1/H)d/dq.
    div_y=sp.diff(yy,t)-sp.diff(yz,q)
    div_z=sp.diff(yz,t)-sp.diff(zz,q)
    assert sp.simplify(div_y)==0 and sp.simplify(div_z)==0
    for edge in (0,1):
        # These equalities concern open boundary faces (q>0), not corners.
        assert sp.simplify(yy.subs(t,edge))==0
        assert sp.simplify(yz.subs(t,edge))==0
    assert sp.simplify(yz.subs(q,0)-1)==0
    assert sp.simplify(zz.subs(q,0))==0
    # Independent check of mode family. One symbolic representative of
    # each nontrivial degree checks the double-zero boundary mechanism.
    a=sp.symbols('a')
    f=t**2*(1-t)**2*(1+a*(2*t-1))
    g=r**2*(1+a*(2*r-1))
    # Phi_mode = H^2 f(t)g(r), r=(L-z)/L.
    m_yy=(H/L)**2*f*sp.diff(g,r,2)
    m_zz=sp.diff(f,t,2)*g
    m_yz=(H/L)*sp.diff(f,t)*sp.diff(g,r)
    assert sp.simplify(sp.diff(m_yy,t)/H-sp.diff(m_yz,r)/L)==0
    assert sp.simplify(sp.diff(m_yz,t)/H-sp.diff(m_zz,r)/L)==0
    for edge in (0,1):
        assert sp.simplify(m_yy.subs(t,edge))==0
        assert sp.simplify(m_yz.subs(t,edge))==0
    assert sp.simplify(m_zz.subs(r,0))==0
    assert sp.simplify(m_yz.subs(r,0))==0
    return {'airy_interior_divergence_identity':True,
            'airy_natural_faces_identity':True,
            'polynomial_mode_interior_divergence_identity':True,
            'polynomial_mode_natural_faces_identity':True,
            'domain_of_claim':'interior and open natural faces; corners excluded'}


def implementation_agreement_audit():
    """Compare production evaluator with independently differentiated Phi.

    Probe points are diagnostic only; they do not prove an infinite-domain
    software identity or provide a verified integration enclosure.
    """
    import numpy as np
    from types import SimpleNamespace
    from aposteriori.equilibrated_cauchy import CantileverAiryStress
    t,q=sp.symbols('t q', real=True)
    F=t-sp.Rational(1,2)-t*(1-t)**2*q**2/(t**2+q**2)+(1-t)*t**2*q**2/((1-t)**2+q**2)
    sy=2*sp.diff(F,q)+q*sp.diff(F,q,2)
    sz=q*sp.diff(F,t,2)
    ss=sp.diff(F,t)+q*sp.diff(F,t,q)
    evaluator=sp.lambdify((t,q),[sy,sz,ss],'numpy',cse=True)
    beam=SimpleNamespace(Lx=4.,Ly=8.,Lz=120.)
    field=CantileverAiryStress(beam,total_load=-12.8)
    xyz=np.array([[0.,0.,1.],[1.,-3.0,110.],[-1.,3.2,119.],
                  [1.0,-1.8,66.],[0.,2.0,120.]])
    tt=xyz[:,1]/beam.Ly+.5
    qq=(beam.Lz-xyz[:,2])/beam.Ly
    yy,zz,yz=evaluator(tt,qq)
    independently=np.zeros((len(xyz),6))
    independently[:,1]=field.T*yy
    independently[:,2]=field.T*zz
    independently[:,4]=field.T*yz
    np.testing.assert_allclose(field.base(xyz), independently, rtol=3e-11,atol=3e-11)
    return {'independent_symbolic_vs_production_base_at_five_points':True,
            'verification_class':'sampling comparison, not exhaustive certification'}


def run_audit(root:Path):
    estimator_path=root/'outputs/session8_checks/fine_estimator.json'
    reference_path=root/'outputs/session8_checks/fine_comparison/comparison_1.json'
    e=json.loads(estimator_path.read_text())
    ref=json.loads(reference_path.read_text())
    eta=float(e['eta_energy_norm']);uh=float(e['computed_energy_norm'])
    relative=100*eta/uh
    converted=100*eta/(uh-eta) if uh>eta else None
    check=lambda a,b:abs(a-b)<1e-10*max(1,abs(a),abs(b))
    assert check(relative,e['estimated_energy_error_percent_computed_normalization'])
    assert check(converted,e['estimated_conservative_exact_normalization_percent'])
    assert e['certified_bound'] is False
    assert ref['certified_bound'] is False
    assert e['quadrature_order']==8 and e['reconstruction']['mode_count']==36
    assert e['reconstruction']['mode_rank']==36
    assert e['mesh']==[4,4,8] and ref['mesh']==[4,4,8]
    assert e['formulation']==ref['formulation'] and e['clamp']==ref['clamp']
    symbolic=symbolic_airy_audit()
    implementation=implementation_agreement_audit()
    files=[estimator_path,reference_path,root/'aposteriori/equilibrated_cauchy.py',
           root/'field_sampling.py',root/'run_session7_equilibrated.py',
           root/'recovered_rotation_hex8.py']
    digests={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    return {
      'status':'RESEARCH_ESTIMATE_AUDITED_NOT_CERTIFIED',
      'certified_bound':False,
      'audit_scope':'saved numerical reports, symbolic Airy identities, artifact hashes; no new FEM solves',
      'symbolic_checks':symbolic,
      'implementation_checks':implementation,
      'numerical_checks':{
        'stored_percentages_recomputed':True,
        'case_and_boundary_condition_match':True,
        'stress_mode_full_rank_reported':True,
        'eta_energy_norm':eta,
        'computed_energy_norm':uh,
        'research_percent_computed_energy':relative,
        'conditional_percent_exact_energy':converted,
        'reference_difference_percent':ref['reference_relative_difference_percent'],
        'sampled_traction_mismatch_psi':e['audit']['maximum_natural_traction_mismatch_psi'],
        'sampled_divergence_lbf_per_in3':e['audit']['maximum_sampled_divergence_lbf_in3'],
      },
      'certificate_gates':{
        'analytic_stress_equilibrium_and_open_face_tractions':'symbolically checked',
        'global_H1_conformity_and_exact_dirichlet_trace':'not formally proven',
        'stress_L2_Hdiv_corner_regularities':'reasoned in docs, not independently proven',
        'rigorous_quadrature_enclosure_for_constitutive_gap':'missing',
        'rigorous_numerical_roundoff_and_linear_algebra_enclosure':'missing',
        'load_geometry_material_and_boundary_model_independent_review':'missing',
        'external_independent_engineering_signoff':'missing',
      },
      'sha256':digests,
      'note':'No positive certification finding is made by this script.'
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out',type=Path,default=Path('outputs/certification_audit.json'))
    a=p.parse_args()
    result=run_audit(a.root.resolve())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'], 'checks':result['symbolic_checks'],
                      'certified_bound':result['certified_bound']},indent=2))

if __name__=='__main__':
    main()
