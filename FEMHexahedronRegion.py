# FEMHexahedronRegion.py
import numpy as np
from FEMPointGrid import FEMPointGrid

# DOF constants
DOF_UX = 0
DOF_UY = 1
DOF_UZ = 2
DOF_MX = 3
DOF_MY = 4
DOF_MZ = 5

class FEMHexahedronRegion(FEMPointGrid):
    # standard parametric directions and labels
    _std_signs = [
        (-1, -1, -1), ( 1, -1, -1), ( 1,  1, -1), (-1,  1, -1),
        (-1, -1,  1), ( 1, -1,  1), ( 1,  1,  1), (-1,  1,  1),
    ]
    _std_labels = [
        ["left",  "bottom", "front"],
        ["right", "bottom", "front"],
        ["right", "top",    "front"],
        ["left",  "top",    "front"],
        ["left",  "bottom", "back"],
        ["right", "bottom", "back"],
        ["right", "top",    "back"],
        ["left",  "top",    "back"],
    ]

    def __init__(self, Lx, Ly, Lz, nx, ny, nz, ndof_per_node=6):
        super().__init__(Lx, Ly, Lz, nx, ny, nz)
        self.ndof_per_node = ndof_per_node
        self._generate_elements()
        self._generate_global_dofs()

    def _node_index(self, i, j, k):
        return i * (self.ny * self.nz) + j * self.nz + k

    def _generate_elements(self):
        elems = []
        for i in range(self.nx - 1):
            for j in range(self.ny - 1):
                for k in range(self.nz - 1):
                    n0 = self._node_index(i,   j,   k)
                    n1 = self._node_index(i+1, j,   k)
                    n2 = self._node_index(i+1, j+1, k)
                    n3 = self._node_index(i,   j+1, k)
                    n4 = self._node_index(i,   j,   k+1)
                    n5 = self._node_index(i+1, j,   k+1)
                    n6 = self._node_index(i+1, j+1, k+1)
                    n7 = self._node_index(i,   j+1, k+1)
                    elems.append([n0, n1, n2, n3, n4, n5, n6, n7])
        self.elements = np.array(elems, dtype=int)

    def _generate_global_dofs(self):
        n_nodes = self.grid_points.shape[0]
        all_dofs = np.arange(n_nodes * self.ndof_per_node)
        self.global_dofs = all_dofs.reshape(n_nodes, self.ndof_per_node)

    def get_num_elements(self):
        return self.elements.shape[0]

    def get_element_node_indices(self, elem_idx):
        return self.elements[elem_idx]

    def get_element_nodes(self, elem_idx):
        idx = self.get_element_node_indices(elem_idx)
        return self.grid_points[idx]

    def get_global_dofs(self):
        return self.global_dofs

    def get_element_global_dofs(self, elem_idx):
        """Returns the global DOFs for all 8 nodes of a given element."""
        node_indices = self.get_element_node_indices(elem_idx)  # shape: (8,)
        return self.global_dofs[node_indices]  # shape: (8, ndof_per_node)

    # ———————————————————————————————————————————————————————
    # Bounding‐box extremeties for the entire region
    # ———————————————————————————————————————————————————————

    def get_extremeties(self):
        """
        Returns the global bounding box as:
        (min_x, max_x, min_y, max_y, min_z, max_z).
        
        Here we center X and Y about zero and let Z start at zero.
        """
        # segment distances
        segX = self.Lx / (self.nx - 1) if self.nx > 1 else 0
        segY = self.Ly / (self.ny - 1) if self.ny > 1 else 0
        segZ = self.Lz / (self.nz - 1) if self.nz > 1 else 0

        min_x, max_x = -self.Lx/2, self.Lx/2
        min_y, max_y = -self.Ly/2, self.Ly/2
        min_z, max_z = 0, segZ * (self.nz - 1)

        return min_x, max_x, min_y, max_y, min_z, max_z

    # ———————————————————————————————————————————————————————
    # Face/edge labeling based on min/max coords
    # ———————————————————————————————————————————————————————

    def get_element_node_labels(self, elem_idx):
        coords = self.get_element_nodes(elem_idx)
        min_x, max_x = coords[:,0].min(), coords[:,0].max()
        min_y, max_y = coords[:,1].min(), coords[:,1].max()
        min_z, max_z = coords[:,2].min(), coords[:,2].max()

        labels = []
        for x, y, z in coords:
            lbl = []
            lbl.append("left"   if x == min_x else "right")
            lbl.append("bottom" if y == min_y else "top")
            lbl.append("front"  if z == min_z else "back")
            labels.append(lbl)
        return labels

    # ———————————————————————————————————————————————————————
    # Winding‐order comparison against ξηζ standard
    # ———————————————————————————————————————————————————————

    def get_element_winding_permutation(self, elem_idx):
        coords = self.get_element_nodes(elem_idx)
        min_x, max_x = coords[:,0].min(), coords[:,0].max()
        min_y, max_y = coords[:,1].min(), coords[:,1].max()
        min_z, max_z = coords[:,2].min(), coords[:,2].max()

        perm = []
        for x, y, z in coords:
            xi   = -1 if x == min_x else 1
            eta  = -1 if y == min_y else 1
            zeta = -1 if z == min_z else 1
            perm.append(self._std_signs.index((xi, eta, zeta)))
        return perm

    def is_standard_winding(self, elem_idx):
        return self.get_element_winding_permutation(elem_idx) == list(range(8))


# ———————————————————————————————————————————————————————
# Quick tests
# ———————————————————————————————————————————————————————

def test_extremeties():
    region = FEMHexahedronRegion(4.0, 8.0, 120.0, 5, 3, 3)
    print("Extremeties:", region.get_extremeties())
    # → (-2.0, 2.0, -4.0, 4.0, 0.0, 120.0)

def test_labels_and_winding():
    region = FEMHexahedronRegion(2.0, 1.0, 1.0, 5, 3, 3)
    e0 = 0
    print("Node labels for elem 0:")
    for i, lbl in enumerate(region.get_element_node_labels(e0)):
        print(f"  Node {i}: {lbl}")
    print("Winding perm:", region.get_element_winding_permutation(e0))
    print("Standard winding?", region.is_standard_winding(e0))

def test_basic():
    region = FEMHexahedronRegion(2.0, 1.0, 1.0, 5, 3, 3)
    print("Grid shape:", region.get_grid_shape())
    print("Total nodes:", region.get_all_points().shape)
    print("Num elems:", region.get_num_elements())
    print("Elem 0 idxs:", region.get_element_node_indices(0))
    print("DOFs at node 0:", region.get_global_dofs()[0])

if __name__ == "__main__":
    test_extremeties()
    test_labels_and_winding()
    test_basic()

