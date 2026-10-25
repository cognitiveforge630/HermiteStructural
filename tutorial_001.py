"""Run a minimal first-degree Cauchy Hex8 cantilever and save results."""
import argparse,json
from pathlib import Path
import numpy as np
from core.hex8 import solve

def main():
 p=argparse.ArgumentParser();p.add_argument('--nx',type=int,default=2);p.add_argument('--ny',type=int,default=2);p.add_argument('--nz',type=int,default=16);p.add_argument('--length',type=float,default=8.);p.add_argument('--E',type=float,default=10_000_000.);p.add_argument('--nu',type=float,default=.3);p.add_argument('--force',type=float,default=-100.);p.add_argument('--weight-density',type=float,default=0.0174,help='specific weight in lbf/in^3, downward y; use 0 to disable');a=p.parse_args()
 xyz,cells,u,F,R=solve(nx=a.nx,ny=a.ny,nz=a.nz,L=a.length,E=a.E,nu=a.nu,tip_force=a.force,return_reactions=True,weight_density=a.weight_density)
 out=Path('outputs');out.mkdir(exist_ok=True)
 np.savez(out/'tutorial001_solution.npz',xyz=xyz,cells=cells,u=u,F=F,reactions=R)
 tip=u[np.isclose(xyz[:,2],a.length),1].mean();r={'element':'Cauchy Hex8 trilinear','nodes':len(xyz),'elements':len(cells),'tip_mean_uy_in':float(tip),'force_y_lbf':float(F[1::3].sum()),'end_traction_force_y_lbf':float(a.force),'self_weight_force_y_lbf':float(-a.weight_density*a.length), 'weight_density_lbf_in3':float(a.weight_density), 'traction_integration':'2x2 Gauss on end facets', 'body_force_integration':'2x2x2 Gauss in every cell','reaction_y_lbf':float(R[:,1].sum()),'reference_EB_tip_uy_in':float(a.force*a.length**3/(3*a.E*(1/12))-a.weight_density*1.0*a.length**4/(8*a.E*(1/12))),'mesh':{'nx':a.nx,'ny':a.ny,'nz':a.nz},'units':'in,lbf,psi'}
 (out/'tutorial001_result.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));print('Saved outputs/tutorial001_solution.npz')
if __name__=='__main__':main()
