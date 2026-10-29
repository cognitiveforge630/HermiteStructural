"""Tiny readable Hex27/Hex64 OOP demonstration; no linear solve or viewer.

From the project root:
    python -m tools.demo_nodal_oop --order 3
"""
from __future__ import annotations

import argparse
import numpy as np

from core.elements.nodal_hex import element_type, face_shape_rule, volume_shape_rule
from core.materials.redwood import redwood_orthotropic_material
from core.reference.beam_theory import E, P
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from core.mesh.nodal_hex import StructuredNodalHexMesh
from core.solvers.nodal_assembly import GlobalAssembler


def show_one_element(order):
    material = redwood_orthotropic_material(E)
    extent = (2., 4., 12.)  # inches; axis-aligned affine example
    element = element_type(order)(extent, material.elastic_matrix(),
                                  node_ids=np.arange((order+1)**3))
    print(f'\nHex{element.node_count}: {element.dof_count} displacement DOFs per cell')
    print('Coordinates: X right, Y up, Z front (free end); back = -Z.')
    print('Corner name -> local/global node ID:')
    for name, node in element.corner_nodes.items():
        print(f'  {name:20s} : {node}')
    print('\nFaces (-Z, +Z, -Y, +Y, -X, +X):')
    for name, facet in element.facets.items():
        print(f'  {name:2s}: normal={facet.normal}, area={facet.area:.6g} in^2, '
              f'{len(facet.node_ids)} interpolation nodes')

    q = element.quadrature_order
    print(f'\nGaussian rule: q={q} per parent coordinate; '
          f'{len(volume_shape_rule(order,q))} volume integration points; '
          f'{len(face_shape_rule(order,q,"+Z"))} face integration points.')
    print('Shape-function sum at parent origin:',
          element.shape_functions((0., 0., 0.)).sum())

    gravity = (0., -DEFAULT_WEIGHT_DENSITY, 0.)  # lbf per cubic inch
    element.integrate_body_force(gravity)
    end_traction = (0., P/(extent[0]*extent[1]), 0.)  # lbf per square inch
    # Facet integration is now owned by the +Z face of THIS element.
    element.facets['+Z'].integrate_traction(end_traction)
    resultant = element.get_force_vector().reshape(-1, 3).sum(axis=0)
    volume = np.prod(extent)
    print('\nConsistent nodal body + end traction resultant:', resultant, 'lbf')
    print('Analytical resultant:', volume*np.array(gravity) +
          element.facets['+Z'].area*np.array(end_traction), 'lbf')
    print('Local stiffness shape:', element.get_stiffness_matrix().shape)
    return element


def show_two_element_assembly(order):
    material = redwood_orthotropic_material(E)
    mesh = StructuredNodalHexMesh(order, divisions=(1,1,2), size=(2.,4.,12.))
    master = element_type(order)(mesh.element_extent, material.elastic_matrix())
    shared_Ke = master.get_stiffness_matrix()
    master.integrate_body_force((0., -DEFAULT_WEIGHT_DENSITY, 0.))
    gravity_Fe = master.get_force_vector()
    assembler = GlobalAssembler(mesh.global_dof_count, batch_elements=1)
    for (ix, iy, iz), element in mesh.iter_elements(material.elastic_matrix(),
                                                     shared_stiffness=shared_Ke):
        element.add_force_vector(gravity_Fe)
        if iz == mesh.divisions[2]-1:
            element.facets['+Z'].integrate_traction((0., P/8., 0.))
        assembler.add_element(element)
    K, F = assembler.finalize()
    print('\nTwo-element assembled sparse K:', K.shape, ', stored entries:', K.nnz)
    print('Global consistent load resultant:', F.reshape(-1,3).sum(axis=0), 'lbf')
    print('No solve or certification performed.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, choices=(2,3), default=3)
    args = parser.parse_args(argv)
    show_one_element(args.order)
    show_two_element_assembly(args.order)


if __name__ == '__main__':
    main()
