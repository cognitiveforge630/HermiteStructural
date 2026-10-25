"""Small portable mesh study. Uses the exact Tutorial001 element and SciPy solver."""
import argparse,csv,json,time
from pathlib import Path
import numpy as np
from core.hex8 import solve

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--levels',type=int,default=3,help='1 to 4; default 3')
 p.add_argument('--E',type=float,default=10_000_000.)
 p.add_argument('--nu',type=float,default=.3)
 p.add_argument('--force',type=float,default=-100.)
 p.add_argument('--length',type=float,default=8.)
 a=p.parse_args()
 if not 1<=a.levels<=4: p.error('--levels must be 1..4')
 if a.E<=0 or not (-1<a.nu<.5) or a.length<=0 or a.force==0:p.error('Invalid physical parameters')
 out=Path('outputs');out.mkdir(exist_ok=True)
 # Cross section 1x1 in; cube-shaped elements throughout.
 levels=[(2,2,16),(3,3,24),(4,4,32),(6,6,48)][:a.levels]
 I=1/12; ref=a.force*a.length**3/(3*a.E*I)
 rows=[];prev=None
 for nx,ny,nz in levels:
  start=time.perf_counter()
  xyz,cells,u,F,R=solve(nx,ny,nz,L=a.length,E=a.E,nu=a.nu,tip_force=a.force,return_reactions=True)
  tip=float(u[np.isclose(xyz[:,2],a.length),1].mean())
  change=None if prev is None else abs(tip-prev)/abs(tip)*100
  row={'nx':nx,'ny':ny,'nz':nz,'elements':len(cells),'tip_uy_in':tip,'euler_bernoulli_uy_in':ref,'reference_difference_pct':abs(tip-ref)/abs(ref)*100,'mesh_change_pct':change,'force_balance_y_lbf':float(F[1::3].sum()+R[:,1].sum()),'seconds':round(time.perf_counter()-start,3)}
  rows.append(row);prev=tip
  print(json.dumps(row),flush=True)
 with (out/'convergence.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 (out/'convergence.json').write_text(json.dumps({'reference':'Euler-Bernoulli cantilever, isotropic rectangular cross section 1x1 in; compare with applicable AWC formulas separately','rows':rows},indent=2))
 print('Saved outputs/convergence.csv and outputs/convergence.json')
if __name__=='__main__':main()
