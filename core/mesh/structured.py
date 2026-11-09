"""Regular Hex8 connectivity: shared 8 vertices, higher-order recovered displacement."""
import numpy as np

def structured_nodes(nx, ny, nz, lx, ly, lz):
    """Return original x/y/z node ordering, preserving numerical reference results."""
    x = np.linspace(-lx / 2, lx / 2, nx)
    y = np.linspace(-ly / 2, ly / 2, ny)
    z = np.linspace(0, lz, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

def structured_hex8_cells(nx, ny, nz):
    """Return standard node ordering used by HermiteStructural element kernels."""
    elems = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                n0 = i * (ny * nz) + j * nz + k
                n1 = (i + 1) * (ny * nz) + j * nz + k
                n2 = (i + 1) * (ny * nz) + (j + 1) * nz + k
                n3 = i * (ny * nz) + (j + 1) * nz + k
                n4 = i * (ny * nz) + j * nz + (k + 1)
                n5 = (i + 1) * (ny * nz) + j * nz + (k + 1)
                n6 = (i + 1) * (ny * nz) + (j + 1) * nz + (k + 1)
                n7 = i * (ny * nz) + (j + 1) * nz + (k + 1)
                elems.append([n0, n1, n2, n3, n4, n5, n6, n7])
    return np.asarray(elems, dtype=int)
