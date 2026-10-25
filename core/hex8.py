"""First-degree Cauchy Hex8; same node order, Voigt convention and Gauss integration as reference_hex8.py in the October certificate research build.
Independent educational extraction: no legacy six-DOF Hermite imports.
Units: inches, lbf, psi.
"""
import numpy as np
from itertools import product
SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],dtype=float)
def shape(xi):
 a=1+SIGNS*np.asarray(xi)
 N=np.prod(a,axis=1)/8
 dN=np.empty((8,3))
 for k in range(3):
  other=[j for j in range(3) if j!=k]
  dN[:,k]=SIGNS[:,k]*a[:,other[0]]*a[:,other[1]]/8
 return N,dN
def B_matrix(xyz,xi):
 N,dN=shape(xi)
 J=dN.T@xyz # parent rows, physical columns, as in reference_hex8
 det=float(np.linalg.det(J))
 if det<=0: raise ValueError('Nonpositive Hex8 Jacobian')
 grad=dN@np.linalg.inv(J)
 B=np.zeros((6,24))
 for i,(dx,dy,dz) in enumerate(grad):
  j=3*i
  B[0,j]=dx; B[1,j+1]=dy; B[2,j+2]=dz
  B[3,j]=dy; B[3,j+1]=dx
  B[4,j+1]=dz; B[4,j+2]=dy
  B[5,j]=dz; B[5,j+2]=dx
 return B,det

def elasticity(E,nu):
 lam=E*nu/((1+nu)*(1-2*nu));mu=E/(2*(1+nu))
 D=np.zeros((6,6));D[:3,:3]=lam;np.fill_diagonal(D[:3,:3],lam+2*mu);D[3,3]=mu;D[4,4]=mu;D[5,5]=mu
 return D

def stiffness(xyz,D):
 K=np.zeros((24,24)); q=1/np.sqrt(3)
 for xi in product((-q,q),repeat=3):
  B,det=B_matrix(xyz,xi);K+=B.T@D@B*det
 return K

def mesh(nx=2,ny=2,nz=4,L=8.0):
 # 1 in by 1 in by L in cantilever, z longitudinal
 coords=np.array([(x,y,z) for z in np.linspace(0,L,nz+1) for y in np.linspace(-.5,.5,ny+1) for x in np.linspace(-.5,.5,nx+1)])
 idx=lambda i,j,k:k*(ny+1)*(nx+1)+j*(nx+1)+i
 cells=[]
 for k in range(nz):
  for j in range(ny):
   for i in range(nx):
    cells.append([idx(i,j,k),idx(i+1,j,k),idx(i+1,j+1,k),idx(i,j+1,k),idx(i,j,k+1),idx(i+1,j,k+1),idx(i+1,j+1,k+1),idx(i,j+1,k+1)])
 return coords,np.asarray(cells,dtype=int)

def solve(nx=2,ny=2,nz=4,L=8.,E=10_000_000.,nu=.3,tip_force=-100.,return_reactions=False,weight_density=0.0):
 from scipy.sparse import lil_matrix
 from scipy.sparse.linalg import spsolve
 xyz,cells=mesh(nx,ny,nz,L);D=elasticity(E,nu)
 K=lil_matrix((3*len(xyz),3*len(xyz)))
 for cell in cells:
  ke=stiffness(xyz[cell],D);dofs=np.array([3*n+a for n in cell for a in range(3)])
  for i,di in enumerate(dofs):
   for j,dj in enumerate(dofs):K[di,dj]+=ke[i,j]
 # Consistent weak-form loads: F = integral_V N^T b dV + integral_end N^T t dS.
 F=np.zeros(3*len(xyz));end=np.flatnonzero(np.isclose(xyz[:,2],L))
 q=1/np.sqrt(3)
 gauss=(-q,q)
 traction=np.array([0.,tip_force/1.0,0.]) # cross-sectional area is 1 in^2
 body=np.array([0.,-weight_density,0.])
 for cell in cells:
  coords=xyz[cell]
  fe=np.zeros((8,3))
  for xi in product(gauss,repeat=3):
   N,dN=shape(xi);jac=dN.T@coords;det=float(np.linalg.det(jac))
   if det<=0: raise ValueError('Invalid volume Jacobian')
   fe+=N[:,None]*body[None,:]*det
  if np.all(np.isclose(coords[[4,5,6,7],2],L)):
   for xi,eta in product(gauss,repeat=2):
    N,dN=shape((xi,eta,1.))
    jac=dN.T@coords
    area=float(np.linalg.norm(np.cross(jac[0],jac[1])))
    fe+=N[:,None]*traction[None,:]*area
  for n,fn in zip(cell,fe):F[3*n:3*n+3]+=fn
 fixed=np.flatnonzero(np.isclose(xyz[:,2],0));locked=np.array([3*n+a for n in fixed for a in range(3)])
 free=np.setdiff1d(np.arange(len(F)),locked);u=np.zeros_like(F);u[free]=spsolve(K.tocsr()[free,:][:,free],F[free])
 reactions=np.asarray(K.tocsr()@u-F)
 if return_reactions: return xyz,cells,u.reshape(-1,3),F,reactions.reshape(-1,3)
 return xyz,cells,u.reshape(-1,3),F
