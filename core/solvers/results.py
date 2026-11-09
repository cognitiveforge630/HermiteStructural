"""Postprocessing of nodal transverse displacements."""
import numpy as np
from .loads import LOAD_DOF

def free_end_uy(beam, U):
    U_by_node = U.reshape((-1, beam.ndof_per_node))
    end_nodes = beam.get_nodes_at_max_z()
    end_uy = U_by_node[end_nodes, LOAD_DOF]
    return end_uy[np.argmax(np.abs(end_uy))]

