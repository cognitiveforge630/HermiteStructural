"""One-sided stress checks at the ORIGINAL source-Hex8 vertices.

No nodal stress averaging, extrapolation from Gauss points, or stress
smoothing is performed. Each element supplies eight independently evaluated
stress values from D @ B_Hermite @ source_U. The exact-corner mode evaluates
(+-1,+-1,+-1); the inset mode evaluates ( +-f, +-f, +-f ) with 0 < f < 1.
Shared mesh vertices may have several element-sided stress values; these
must never be silently averaged into one ostensibly exact nodal value.
"""
from __future__ import annotations

import numpy as np
from core.visualization.cubic_stress import analytical_stress_at

_SIGNS = ((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
          (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))


def audit_source_vertex_stress(beam, U, *, location='exact_corners',
                               inset_fraction=0.95, end_force, density):
    """Return original-vertex one-sided stresses and signed beam comparisons.

    The L2 comparison is a *discrete sample norm*, not a certified 3D norm.
    An exact-corner stress may be nonunique across source-element sides and
    physical corners/supports can feature high stress gradients/singularities.
    """
    if location not in ('exact_corners','inset_corners'):
        raise ValueError('location must be exact_corners or inset_corners')
    if not np.isfinite(inset_fraction) or not 0 < inset_fraction < 1:
        raise ValueError('inset_fraction must be finite and strictly between 0 and 1')
    f=1.0 if location=='exact_corners' else float(inset_fraction)
    u=np.asarray(U,dtype=float)
    if u.shape != (len(beam.grid_points)*beam.ndof_per_node,) or not np.isfinite(u).all():
        raise ValueError('Expected finite full solver DOFs')
    if not np.isfinite(end_force) or not np.isfinite(density):
        raise ValueError('Expected finite reference loads')
    D=beam.elastic_matrix()
    record=[]
    points=[]
    stresses=[]
    sides={}
    for e in range(beam.get_num_elements()):
        coords=beam.get_element_nodes(e)
        ids,T=beam.element_transform(e)
        primitive=T@u[ids]
        node_ids=beam.elements[e]
        for vertex,s in enumerate(_SIGNS):
            parent=tuple(float(f*x) for x in s)
            B,detJ,_=beam.primitive_B(*parent,coords)
            if not np.isfinite(detJ) or detJ <= 0 or not np.all(np.isfinite(B)):
                raise ArithmeticError(f'Nonfinite/singular evaluation at element {e}, vertex {vertex}')
            p=beam.hex8_shape_functions(*parent)@coords
            sig=D@(B@primitive)
            points.append(p)
            stresses.append(sig)
            node_id=int(node_ids[vertex])
            if location=='exact_corners':
                sides.setdefault(node_id,[]).append(sig)
            record.append({'source_element_index':int(e),
                           'source_corner_index':vertex,
                           'source_node_index':node_id if location=='exact_corners' else None,
                           'parent_coordinates':list(parent),
                           'physical_xyz_in':[float(q) for q in p],
                           'stress_psi':[float(v) for v in sig]})
    points=np.asarray(points)
    stresses=np.asarray(stresses)
    beam_sigma,beam_tau=analytical_stress_at(points,beam,end_force=end_force,density=density)
    for row,sig,tau in zip(record,beam_sigma,beam_tau):
        row['beam_reference_sigma_zz_psi']=float(sig)
        row['beam_reference_tau_yz_psi']=float(tau)
        row['sigma_zz_minus_beam_psi']=row['stress_psi'][2]-float(sig)
        row['tau_yz_minus_beam_psi']=row['stress_psi'][4]-float(tau)
    interior=(points[:,2]>=beam.Lz/4)&(points[:,2]<=3*beam.Lz/4)
    if not interior.any():
        raise ValueError('No interior quarter-span to three-quarter-span vertex samples')
    spread=max((np.max(np.ptp(np.stack(v),axis=0)) for v in sides.values() if len(v)>1), default=0.)
    def relative_l2(actual,ref):
        return float(100*np.linalg.norm(actual[interior]-ref[interior])/
                     max(np.linalg.norm(ref[interior]),1e-30))
    return {'status':'RESEARCH_SOURCE_VERTEX_STRESS_NOT_CERTIFIED',
            'formulation':type(beam).__name__,
            'stress_evaluated_from':'D @ source_element_B @ source_element_U',
            'sampling':'source-element-sided original vertices' if location=='exact_corners' else 'inset towards original source-element vertices',
            'sampling_location':location,
            'inset_fraction':float(f),
            'stress_component_order':['sigma_xx','sigma_yy','sigma_zz','tau_xy','tau_yz','tau_xz'],
            'sample_count':len(record),
            'unique_original_node_count':int(len(beam.grid_points)),
            'interior_sample_count':int(interior.sum()),
            'interior_sigma_zz_discrete_relative_L2_pct':relative_l2(stresses[:,2],beam_sigma),
            'interior_tau_yz_discrete_relative_L2_pct':relative_l2(stresses[:,4],beam_tau),
            'raw_max_abs_sigma_zz_psi':float(np.max(abs(stresses[:,2]))),
            'raw_max_abs_tau_yz_psi':float(np.max(abs(stresses[:,4]))),
            'max_shared_vertex_element_side_spread_psi':float(spread) if location=='exact_corners' else None,
            'no_cross_element_nodal_stress_averaging':True,
            'certified_error_bound':False,
            'limitations':('Discrete, equally weighted original-source-vertex sample L2 is NOT a volume norm, '
                           'may be sensitive to mesh density, load/clamp singularities and element-side differences; '
                           'a beam-theory shear profile is not a unique exact 3D solution; '
                           'these stresses are not equilibrated and no certification is established.'),
            'corner_samples':record}
