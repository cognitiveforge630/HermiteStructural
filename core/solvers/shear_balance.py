"""Axial-section shear/moment audit of *raw* cubic Hermite constitutive stresses.

No stress is modified, equilibrated, substituted, or used to correct the solve.
In particular, averaging is an additional diagnostic and never replaces the
worst pointwise sectional discrepancy. This is not a posteriori certification.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from core.solvers.section_equilibrium import integrate_sections


def _differentiation_matrix(x):
    """Derivative of the interpolating polynomial at distinct sample abscissae."""
    x=np.asarray(x,dtype=float)
    if x.ndim != 1 or len(x)<3 or not np.isfinite(x).all() or len(set(x))!=len(x):
        raise ValueError('Expected distinct finite axial points')
    w=np.empty(len(x),float)
    for j in range(len(x)):
        w[j]=1./np.prod(x[j]-np.delete(x,j))
    D=np.empty((len(x),len(x)),float)
    for i in range(len(x)):
        for j in range(len(x)):
            if i!=j:
                D[i,j]=w[j]/(w[i]*(x[i]-x[j]))
        D[i,i]=-sum(D[i,j] for j in range(len(x)) if j!=i)
    return D


def audit_axial_shear_balance(beam,U,*,order=4,section_order=4,end_force,density):
    """Contrast section-by-section raw FE shear, its cell-average, and dM/dz.

    The axial integration uses a separate Gauss rule; stress on each section
    uses the original source B operator and cross-section Gauss quadrature.
    dM/dz is differentiated from the axial interpolation of **FE** section
    moments (not from the reference), and is diagnostic, not a pointwise exact
    derivative. The signed moment and shear convention follows beam_theory.
    """
    if type(order)!=int or not 3<=order<=8:
        raise ValueError('Axial quadrature order must be an integer 3..8')
    gp,gw=leggauss(order)
    sections=integrate_sections(beam,U,order=section_order,axial_positions=tuple(gp),
                               end_force=end_force,density=density)
    D=_differentiation_matrix(gp)
    length=beam.Lz/(beam.nz-1)
    cell_reports=[]
    for iz in range(beam.nz-1):
        rows=sections['sections'][iz*order:(iz+1)*order]
        feV=np.array([row['FE_shear_y_lbf'] for row in rows])
        refV=np.array([row['beam_shear_y_lbf'] for row in rows])
        feM=np.array([row['FE_bending_moment_x_lbf_in'] for row in rows])
        refM=np.array([row['beam_bending_moment_x_lbf_in'] for row in rows])
        # Convert parent-coordinate derivative to physical z derivative.
        moment_slope=(2./length)*(D@feM)
        moment_slope_ref=(2./length)*(D@refM)
        scale=np.maximum(abs(refV),1e-12)
        errors=100*abs(feV-refV)/scale
        mean_fe=float(np.dot(gw,feV)/2)
        mean_ref=float(np.dot(gw,refV)/2)
        mean_error=100*abs(mean_fe-mean_ref)/max(abs(mean_ref),1e-12)
        rel_residual=100*abs(feV-moment_slope)/scale
        # All extrema reference sampled locations, NOT an assertion about
        # unsampled locations; the continuum maximum remains unknown.
        cell_reports.append({
            'axial_cell_index':iz,
            'z_interval_in':[float(iz*length),float((iz+1)*length)],
            'source_section_count':order,
            'samples':[{
                'z_in':r['z_in'],
                'parent_zeta':r['axial_parent_zeta'],
                'raw_FE_shear_y_lbf':float(feV[i]),
                'reference_shear_y_lbf':float(refV[i]),
                'raw_FE_moment_x_lbf_in':float(feM[i]),
                'reference_moment_x_lbf_in':float(refM[i]),
                'raw_shear_reference_difference_pct':float(errors[i]),
                'interpolated_dFE_moment_dz_lbf':float(moment_slope[i]),
                'FE_shear_minus_interpolated_dFE_moment_dz_lbf':float(feV[i]-moment_slope[i]),
                'interpolated_dBeam_moment_dz_lbf':float(moment_slope_ref[i]),
                'shear_moment_balance_difference_pct_of_reference_shear':float(rel_residual[i]),
            } for i,r in enumerate(rows)],
            'axially_averaged_raw_FE_shear_y_lbf':mean_fe,
            'axially_averaged_beam_shear_y_lbf':mean_ref,
            'axially_averaged_shear_reference_difference_pct':float(mean_error),
            'max_sampled_shear_reference_difference_pct':float(np.max(errors)),
            'max_sampled_shear_moment_balance_difference_pct':float(np.max(rel_residual)),
            'FE_shear_sample_min_lbf':float(np.min(feV)),
            'FE_shear_sample_max_lbf':float(np.max(feV)),
            'raw_FE_shear_sample_range_lbf':float(np.ptp(feV)),
        })
    interior=[c for c in cell_reports if c['z_interval_in'][0]>=beam.Lz/4
              and c['z_interval_in'][1]<=3*beam.Lz/4]
    if not interior:  # Retain a meaningful sample on very coarse grids.
        interior=[c for c in cell_reports if c['z_interval_in'][0]<3*beam.Lz/4
                  and c['z_interval_in'][1]>beam.Lz/4]
    return {
        'status':'RESEARCH_AXIAL_SHEAR_MOMENT_DIAGNOSTIC_NOT_CERTIFIED',
        'formulation':type(beam).__name__,
        'axial_gauss_order':order,
        'cross_section_gauss_order':int(section_order),
        'uses_actual_source_FE_stress':True,
        'does_not_modify_solver_or_stress':True,
        'FE_moment_derivative_method':'Barycentric polynomial differentiation of section-integrated FE moments at axial Gauss points',
        'max_interior_sampled_shear_reference_difference_pct':float(max(c['max_sampled_shear_reference_difference_pct'] for c in interior)),
        'max_interior_axially_averaged_shear_reference_difference_pct':float(max(c['axially_averaged_shear_reference_difference_pct'] for c in interior)),
        'max_interior_sampled_shear_moment_balance_difference_pct':float(max(c['max_sampled_shear_moment_balance_difference_pct'] for c in interior)),
        'cells':cell_reports,
        'certified_error_bound':False,
        'limitations':('Sampled sectional shear and approximate FE moment derivative are diagnostics. '
                       'Beam reference is not exact 3D elasticity. Axial averaging can conceal '
                       'oscillations and must not be interpreted as correction or certification. '
                       'No rigorous local maxima, traction admissibility, or a posteriori bound.'),
    }
