"""Internal-face source-stress traction jump *diagnostic*, not an estimator.

The Cauchy stress reconstructed from a conforming FE displacement is may or may not have jumps depending on the interpolation. For shared-vertex
recovered Hermite shapes in this project, it is numerically C1, so the raw
constitutive stress appears continuous across the internal faces. Face jumps
alone cannot establish equilibrium. This module is research verification,
not a certified bound or an equilibrated stress reconstruction.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss


def source_traction_jumps(beam, U, *, order=4):
    """RMS [sigma]n on internal structured Hex8 faces using source B matrices.

    Each face is evaluated with its positive-coordinate normal orientation
    (x/y/z), with the left-side stress minus right-side stress. Units are psi.
    """
    if isinstance(order,bool) or int(order)!=order or order<2:
        raise ValueError('quadrature order must be an integer >= 2')
    u=np.asarray(U,float)
    if u.shape!=(len(beam.grid_points)*6,) or not np.isfinite(u).all():
        raise ValueError('Expected finite complete six-DOF displacement field')
    gp,w=leggauss(int(order))
    D=beam.elastic_matrix()
    limits=(beam.nx-1,beam.ny-1,beam.nz-1)
    normals=((0,3,5),(3,1,4),(5,4,2))
    stats={}
    # Element IDs follow the source structured grid's ix/iy/iz ordering.
    def elem_id(ix,iy,iz):
        return (ix*limits[1]+iy)*limits[2]+iz
    for axis in range(3):
        aggregate_squared=0.
        aggregate_area=0.
        highest=(0.,None)
        face_count=0
        for ix in range(limits[0]):
            for iy in range(limits[1]):
                for iz in range(limits[2]):
                    index=(ix,iy,iz)
                    if index[axis]>=limits[axis]-1:
                        continue
                    following=list(index);following[axis]+=1
                    e0,e1=elem_id(*index),elem_id(*following)
                    coords0,coords1=beam.get_element_nodes(e0),beam.get_element_nodes(e1)
                    dofs0,T0=beam.element_transform(e0)
                    dofs1,T1=beam.element_transform(e1)
                    vec0,vec1=T0@u[dofs0],T1@u[dofs1]
                    integrated=0.
                    area=0.
                    other=[j for j in range(3) if j!=axis]
                    for a,qa in enumerate(gp):
                        for b,qb in enumerate(gp):
                            parent0=np.zeros(3)
                            parent1=np.zeros(3)
                            parent0[axis],parent1[axis]=1.,-1.
                            parent0[other]=qa,qb
                            parent1[other]=qa,qb
                            B0,_,_=beam.primitive_B(*parent0,coords0)
                            B1,_,_=beam.primitive_B(*parent1,coords1)
                            sig0=D@(B0@vec0)
                            sig1=D@(B1@vec1)
                            traction_jump=(sig0-sig1)[list(normals[axis])]
                            J=beam.get_hermite_jacobian(*parent0,coords0)
                            dA=np.linalg.norm(np.cross(J[other[0]],J[other[1]]))*w[a]*w[b]
                            integrated+=float(traction_jump@traction_jump)*dA
                            area+=dA
                    rms=float(np.sqrt(max(integrated,0.)/area))
                    if rms>highest[0]:
                        highest=(rms,[int(ix),int(iy),int(iz)])
                    aggregate_squared+=integrated
                    aggregate_area+=area
                    face_count+=1
        stats['xyz'[axis]]={
            'faces':face_count,
            'area_in2':float(aggregate_area),
            'area_weighted_traction_jump_rms_psi':float(np.sqrt(aggregate_squared/aggregate_area)) if aggregate_area else 0.,
            'maximum_face_traction_jump_rms_psi':highest[0],
            'maximum_face_lower_element_ijk':highest[1],
        }
    return {
        'status':'RESEARCH_RAW_STRESS_FACE_JUMPS_NOT_CERTIFIED',
        'source_stress_evaluated_from':'D @ B @ U on both sides of face',
        'source_face_quadrature_order':int(order),
        'stress_jump_units':'psi',
        'per_face_normal_axis':stats,
        'certified_bound':False,
        'a_posteriori_terms_not_included':[
            'volume equilibrium residual', 'boundary natural traction residual',
            'data oscillation', 'constitutive compliance scaling',
            'reliability/efficiency proof', 'numerical enclosure'],
    }
