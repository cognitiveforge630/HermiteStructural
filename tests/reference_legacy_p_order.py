"""FROZEN BEFORE-REFACTOR REFERENCE from input tools/p_order_study.py.

Purpose: independent numerical regression for shape functions, B, Ke, body
loading, and +Z end traction. NOT used by production OOP solver.
"""
import itertools
import numpy as np
from numpy.polynomial.legendre import leggauss

def lagrange_1d(p: int, s: float):
    if p not in (2, 3):
        raise ValueError('Only Hex27 (p=2) and Hex64 (p=3) are implemented')
    nodes = np.linspace(-1.0, 1.0, p + 1)
    values = np.ones(p + 1)
    derivs = np.zeros(p + 1)
    for i in range(p + 1):
        for j in range(p + 1):
            if i != j:
                values[i] *= (s - nodes[j]) / (nodes[i] - nodes[j])
        for j in range(p + 1):
            if i == j:
                continue
            factor = 1.0 / (nodes[i] - nodes[j])
            for k in range(p + 1):
                if k != i and k != j:
                    factor *= (s - nodes[k]) / (nodes[i] - nodes[k])
            derivs[i] += factor
    return values, derivs


def tensor_shapes(p, parent):
    """Values and parent derivatives, tensor index order (i,j,k), k fastest."""
    n = p + 1
    a = [lagrange_1d(p, float(t)) for t in parent]
    val = np.empty(n**3)
    ds = np.empty((3, n**3))
    for m, (i, j, k) in enumerate(itertools.product(range(n), repeat=3)):
        val[m] = a[0][0][i] * a[1][0][j] * a[2][0][k]
        ds[0, m] = a[0][1][i] * a[1][0][j] * a[2][0][k]
        ds[1, m] = a[0][0][i] * a[1][1][j] * a[2][0][k]
        ds[2, m] = a[0][0][i] * a[1][0][j] * a[2][1][k]
    return val, ds


def strain_matrix(dN):
    """Engineering-Voigt order xx, yy, zz, xy, yz, xz."""
    nodes = dN.shape[1]
    b = np.zeros((6, 3 * nodes))
    b[0, 0::3] = dN[0]
    b[1, 1::3] = dN[1]
    b[2, 2::3] = dN[2]
    b[3, 0::3] = dN[1]
    b[3, 1::3] = dN[0]
    b[4, 1::3] = dN[2]
    b[4, 2::3] = dN[1]
    b[5, 0::3] = dN[2]
    b[5, 2::3] = dN[0]
    return b


def element_kernel(p, xyz_extent, D, density, traction, q):
    """One affine cell; cache Ke, gravity volume load, and end-face traction."""
    scale=np.asarray(xyz_extent,float)*0.5
    detJ=float(np.prod(scale))
    gp,gw=leggauss(q)
    n=3*(p+1)**3
    ke=np.zeros((n,n))
    fg=np.zeros(n)
    ft=np.zeros(n)
    for i,xi in enumerate(gp):
        for j,eta in enumerate(gp):
            for k,zeta in enumerate(gp):
                N,ds=tensor_shapes(p,(xi,eta,zeta))
                b=strain_matrix(ds/scale[:,None])
                w=gw[i]*gw[j]*gw[k]*detJ
                ke+=(b.T@D@b)*w
                fg[1::3]+=-density*N*w
            face_N,_=tensor_shapes(p,(xi,eta,1.0))
            ft[1::3]+=traction*face_N*gw[i]*gw[j]*scale[0]*scale[1]
    return ke,fg,ft


