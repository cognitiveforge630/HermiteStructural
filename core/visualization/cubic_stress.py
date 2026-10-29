"""Independent classical beam reference and research-grade cubic stress audit.

The actual FE stress is always D @ B_Hermite @ U at source-element
corner (or explicitly inset-corner) samples. Analytical beam stresses are stored separately. In particular, a
smooth theoretical field must NEVER be relabeled as an FE-derived result.

Local +z is the grain direction; bending is about global x (transverse y).
The comparison is Saint-Venant beam theory, not an exact 3D solution and
not an a posteriori certificate.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from core.materials.design_limits import REDWOOD_BASIS, reference_load_screen
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from core.reference.beam_theory import P


def analytical_stress_at(points, beam, *, end_force=P, density=DEFAULT_WEIGHT_DENSITY):
    """Return beam-theory σzz and approximate rectangular-section τyz (psi)."""
    points=np.asarray(points,float)
    if points.ndim!=2 or points.shape[1]!=3:
        raise ValueError("Expected (n,3) points")
    y,z=points[:,1],points[:,2]
    remaining=beam.Lz-z
    q=density*beam.Lx*beam.Ly
    inertia=beam.Lx*beam.Ly**3/12.
    # Downward end force and gravity; clamp top fiber is in tension.
    M=-end_force*remaining + 0.5*q*remaining**2
    V=end_force-q*remaining
    sigma=M*y/inertia
    tau=1.5*V/(beam.Lx*beam.Ly)*(1.-(2.*y/beam.Ly)**2)
    return sigma,tau


def attach_cubic_reference_fields(sample, beam, *, end_force=P, density=DEFAULT_WEIGHT_DENSITY):
    """Adds DISTINCT reference and residual arrays to the existing source-FE stress field."""
    # Reference must be sampled at the IDENTICAL physical corner locations
    # as the source FE stress. Mean is for VTU cell coloring only.
    corner=sample.stress_corner_points.reshape(-1,3)
    sigma,tau=analytical_stress_at(corner,beam,end_force=end_force,density=density)
    ref_z=sigma.reshape(-1,8).mean(axis=1)
    ref_shear=tau.reshape(-1,8).mean(axis=1)
    actual_z=np.asarray(sample.cell_fields['Sigma zz (psi)'])
    actual_shear=np.asarray(sample.cell_fields['Tau yz (psi)'])
    sample.cell_fields['Beam reference sigma zz (psi)']=ref_z
    sample.cell_fields['Cubic minus beam sigma zz (psi)']=actual_z-ref_z
    sample.cell_fields['Beam reference tau yz (psi)']=ref_shear
    sample.cell_fields['Cubic minus beam tau yz (psi)']=actual_shear-ref_shear
    sample.cell_fields['Longitudinal tension utilization']=np.maximum(actual_z,0.)/REDWOOD_BASIS.tension_parallel_psi
    sample.cell_fields['Longitudinal compression utilization']=np.maximum(-actual_z,0.)/REDWOOD_BASIS.compression_parallel_psi
    sample.cell_fields['Shear yz utilization']=np.abs(actual_shear)/REDWOOD_BASIS.shear_parallel_psi
    return sample


def cubic_stress_audit(sample, beam, *, end_force=P, density=DEFAULT_WEIGHT_DENSITY):
    """Review sample stresses and section force/moment, without overclaiming accuracy."""
    attach_cubic_reference_fields(sample,beam,end_force=end_force,density=density)
    center=sample.points[sample.connectivity].mean(axis=1)
    actual=sample.cell_fields['Sigma zz (psi)']
    expected=sample.cell_fields['Beam reference sigma zz (psi)']
    # Stress comparisons are from the individual node/corner-side samples,
    # not the possibly cancellation-prone cell averages.
    stress_points=sample.stress_corner_points.reshape(-1,3)
    raw_stress=sample.stress_corner_values.reshape(-1,6)
    raw_ref,raw_tau_ref=analytical_stress_at(stress_points,beam,end_force=end_force,density=density)
    interior=(stress_points[:,2]>=beam.Lz/4.)&(stress_points[:,2]<=beam.Lz*3./4.)
    if not interior.any():
        raise ValueError('Sampling contains no interior quarter-span to three-quarter-span stress corners')
    delta=raw_stress[interior,2]-raw_ref[interior]
    relative_l2=np.linalg.norm(delta)/max(np.linalg.norm(raw_ref[interior]),1e-30)
    stress_peaks={
        'longitudinal_tension_psi':float(max(0.,np.max(raw_stress[:,2]))),
        'longitudinal_compression_psi':float(max(0.,-np.min(raw_stress[:,2]))),
        'shear_yz_psi':float(np.max(np.abs(raw_stress[:,4]))),
        'shear_xz_psi':float(np.max(np.abs(raw_stress[:,5]))),
        'transverse_compression_xx_psi':float(max(0.,-np.min(raw_stress[:,0]))),
        'transverse_compression_yy_psi':float(max(0.,-np.min(raw_stress[:,1]))),
        'transverse_tension_xx_psi':float(max(0.,np.max(raw_stress[:,0]))),
        'transverse_tension_yy_psi':float(max(0.,np.max(raw_stress[:,1]))),
    }
    limitations={
        'longitudinal_tension':REDWOOD_BASIS.tension_parallel_psi,
        'longitudinal_compression':REDWOOD_BASIS.compression_parallel_psi,
        'shear_yz':REDWOOD_BASIS.shear_parallel_psi,
        'shear_xz':REDWOOD_BASIS.shear_parallel_psi,
        'transverse_compression_xx':REDWOOD_BASIS.compression_perpendicular_psi,
        'transverse_compression_yy':REDWOOD_BASIS.compression_perpendicular_psi,
    }
    # Published source has no perpendicular-to-grain tensile strength value.
    utilization={k:stress_peaks[k+'_psi']/v for k,v in limitations.items()}
    # The spatial-accuracy indicator is a beam-theory comparison, not a
    # guarantee for transverse/shear stress or near-clamp stress concentrations.
    area=beam.Lx*beam.Ly
    cross_section_subcells=(beam.nx-1)*(beam.ny-1)*sample.subdivisions**2
    da=area/cross_section_subcells
    slices=[]
    for z in np.unique(np.round(center[:,2],decimals=9)):
        mask=np.isclose(center[:,2],z,atol=1.e-7)
        mx=float(np.sum(actual[mask]*center[mask,1])*da)
        ax=float(np.sum(actual[mask])*da)
        shear=float(np.sum(sample.cell_fields['Tau yz (psi)'][mask])*da)
        q=density*area
        mex=-end_force*(beam.Lz-z)+q*(beam.Lz-z)**2/2
        vex=end_force-q*(beam.Lz-z)
        slices.append({'z_in':float(z),'FE_bending_moment_lbf_in':mx,
                       'beam_bending_moment_lbf_in':float(mex),
                       'moment_error_pct':float(100*abs(mx-mex)/max(abs(mex),1e-12)),
                       'FE_axial_force_lbf':ax,'FE_shear_y_lbf':shear,
                       'beam_shear_y_lbf':float(vex),
                       'shear_error_pct':float(100*abs(shear-vex)/max(abs(vex),1e-12))})
    mid_slices=[v for v in slices if beam.Lz/4<=v['z_in']<=beam.Lz*3/4]
    ref=reference_load_screen(tip_force_lbf=end_force,density_lbf_in3=density,
                              length_in=beam.Lz,width_in=beam.Lx,depth_in=beam.Ly,
                              modulus_psi=beam.orthotropic_constants.Ez)
    return {
        'status':('PASS_PRELIMINARY_SAMPLE_SCREEN' if ref['status']=='PASS_PRELIMINARY'
                   and all(v<=1. for v in utilization.values()) else 'FAIL_PRELIMINARY_SCREEN'),
        'model':'research cubic Hermite; one-sided constitutive FE stress at subcell corners',
        'stress_sample_location':sample.stress_sample_location,
        'stress_inset_fraction':float(sample.stress_inset_fraction),
        'raw_corner_sample_count':int(raw_stress.shape[0]),
        'raw_corner_sigma_zz_relative_L2_pct':float(100*relative_l2),
        'stress_peaks_are_from':'individual one-sided subcell corner values (not cell mean)',
        'display_cell_stresses_are':'arithmetic mean of eight one-sided source-element corner samples',
        'stress_distribution_certified':False,
        'sample_subdivisions':int(sample.subdivisions),
        'material_reference_grade':REDWOOD_BASIS.grade,
        'material_source':REDWOOD_BASIS.source,
        'material_E_longitudinal_psi':beam.orthotropic_constants.Ez,
        'screening_limitations':'Unadjusted grade values and sampling only; no tension-perpendicular, stability, connection or full NDS verification',
        'beam_reference_screen':ref,
        'sampled_stress_peaks':stress_peaks,
        'sampled_utilization':utilization,
        'quarter_to_threequarter_span_sigma_zz_relative_L2_pct':float(100*relative_l2),
        'quarter_to_threequarter_span_max_moment_error_pct':float(max(v['moment_error_pct'] for v in mid_slices)),
        'quarter_to_threequarter_span_max_shear_error_pct':float(max(v['shear_error_pct'] for v in mid_slices)),
        'section_resultants':slices,
        'caveat':'Raw corner-sample L2 is not directly comparable with prior centroid-sample L2. Display cell values are corner means, NOT nodal extrema. Beam theory is not exact 3D elasticity; no certified stress bound.'
    }


def write_stress_audit(audit, destination):
    destination=Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(audit,indent=2)+'\n',encoding='utf8')
    csvpath=destination.with_name(destination.stem+'_sections.csv')
    with csvpath.open('w',newline='',encoding='utf8') as f:
        sections=audit['section_resultants']
        writer=csv.DictWriter(f, fieldnames=list(sections[0]))
        writer.writeheader()
        writer.writerows(sections)
    return destination,csvpath
