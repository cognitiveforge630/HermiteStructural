"""Batched field evaluation on the supported axis-aligned structured bricks.

Element support is transformed once; primitive matrices are evaluated once
per geometry/parent-point batch. Bounded caches avoid retaining all matrices
in a large comparison. Each sampler belongs to one fixed beam and U snapshot.
"""
from collections import OrderedDict
import numpy as np


def tensor_quadrature(order):
    gp,gw=np.polynomial.legendre.leggauss(order)
    points=np.stack(np.meshgrid(gp,gp,gp,indexing='ij'),axis=-1).reshape(-1,3)
    weights=np.prod(np.stack(np.meshgrid(gw,gw,gw,indexing='ij'),axis=-1),axis=-1).ravel()
    return points,weights


class StructuredField:
    def __init__(self,beam,U,cache_bytes=64*1024*1024):
        self.beam=beam
        self.U=np.array(U,copy=True)
        self.coords=[beam.get_element_nodes(e) for e in range(beam.get_num_elements())]
        self.low=np.array([xyz.min(axis=0) for xyz in self.coords])
        self.high=np.array([xyz.max(axis=0) for xyz in self.coords])
        self.size=self.high-self.low
        if np.any(self.size<=0):
            raise ValueError('Degenerate brick')
        self.planes=[np.unique(np.r_[self.low[:,k],self.high[:,k]]) for k in range(3)]
        shape=tuple(len(p)-1 for p in self.planes)
        if np.prod(shape)!=len(self.coords):
            raise ValueError('Field sampling requires a complete Cartesian brick grid')
        self.lookup=np.full(shape,-1,dtype=int)
        for e,lo in enumerate(self.low):
            index=tuple(np.searchsorted(self.planes[k],lo[k]) for k in range(3))
            if self.lookup[index]!=-1:
                raise ValueError('Overlapping elements')
            self.lookup[index]=e
        self.coefficients={}
        self.cache=OrderedDict()
        self.cache_bytes=0
        self.limit=max(0,int(cache_bytes))
        self.hits=self.misses=0

    def locate(self,point):
        point=np.asarray(point)
        index=tuple(int(np.searchsorted(self.planes[k],point[k],side='right')-1) for k in range(3))
        if any(i<0 or i>=self.lookup.shape[k] for k,i in enumerate(index)):
            raise ValueError('Point outside the element grid')
        e=int(self.lookup[index])
        if e<0 or np.any(point<=self.low[e]) or np.any(point>=self.high[e]):
            raise ValueError('Lookup requires an element-interior point')
        return e,self.low[e],self.high[e]

    def primitive_coefficients(self,e):
        if e not in self.coefficients:
            if hasattr(self.beam,'element_transform'):
                ids,T=self.beam.element_transform(e)
                self.coefficients[e]=T @ self.U[ids]
            else:
                self.coefficients[e]=self.U[self.beam.get_element_global_dofs(e).ravel()]
        return self.coefficients[e]

    def tables(self,e,parents,displacement=False):
        parents=np.asarray(parents,dtype=float).reshape(-1,3)
        # Parent differences below roundoff already share keys in the legacy
        # comparator. Points themselves are evaluated as provided on a miss.
        key=(tuple(np.round(self.size[e],12)),np.round(parents,14).tobytes(),bool(displacement))
        if key in self.cache:
            self.hits+=1
            self.cache.move_to_end(key)
            return self.cache[key][0]
        self.misses+=1
        strains,values,determinants=[],[],[]
        for p in parents:
            if hasattr(self.beam,'element_transform'):
                B,det,N=self.beam.primitive_B(*p,self.coords[e])
            else:
                B,det,N,_=self.beam.build_B_matrix(*p,self.coords[e])
            strains.append(B)
            if displacement:
                values.append(N)
            determinants.append(det)
        result=(np.array(strains),np.array(determinants),np.array(values) if displacement else None)
        size=sum(v.nbytes for v in result if v is not None)
        while self.cache and self.cache_bytes+size>self.limit:
            _,(_,old_size)=self.cache.popitem(last=False)
            self.cache_bytes-=old_size
        if size<=self.limit:
            self.cache[key]=(result,size)
            self.cache_bytes+=size
        return result

    def strain(self,e,parents):
        B,_,_=self.tables(e,parents)
        coefficients=self.primitive_coefficients(e)
        return (B.reshape(-1,len(coefficients)) @ coefficients).reshape(-1,6)

    def displacement_and_strain(self,e,parents):
        B,_,N=self.tables(e,parents,displacement=True)
        coefficients=self.primitive_coefficients(e)
        return ((N.reshape(-1,len(coefficients)) @ coefficients).reshape(-1,3),
                (B.reshape(-1,len(coefficients)) @ coefficients).reshape(-1,6))

    def statistics(self):
        return {'basis_cache_hits':self.hits,'basis_cache_misses':self.misses,
                'basis_cache_bytes':self.cache_bytes,'basis_cache_limit_bytes':self.limit}
