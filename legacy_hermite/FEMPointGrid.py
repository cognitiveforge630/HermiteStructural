# FEMPointGrid.py
import numpy as np

class FEMPointGrid:
    def __init__(self, Lx, Ly, Lz, nx, ny, nz):
        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz
        self.nx, self.ny, self.nz = nx, ny, nz
        
        self._generate_grid()
    
    def _generate_grid(self):
        x = np.linspace(-self.Lx/2, self.Lx/2, self.nx)
        y = np.linspace(-self.Ly/2, self.Ly/2, self.ny)
        z = np.linspace(0, self.Lz, self.nz)

        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        self.grid_points = np.column_stack((self.X.ravel(), self.Y.ravel(), self.Z.ravel()))

    def get_grid_shape(self):
        return self.X.shape

    def get_point(self, i, j, k):
        return (self.X[i, j, k], self.Y[i, j, k], self.Z[i, j, k])

    def get_all_points(self):
        return self.grid_points

    def get_total_points(self):
        return self.grid_points.shape[0]

if __name__ == "__main__":
    region = FEMPointGrid(Lx=2.0, Ly=1.0, Lz=1.0, nx=5, ny=3, nz=3)
    
    print("Grid shape:", region.get_grid_shape())        # → (5, 3, 3)
    print("Corner point:", region.get_point(0, 0, 0))     # → (-1.0, -0.5, 0.0)
    print("Total nodes:", region.get_all_points().shape) # → (45, 3)
