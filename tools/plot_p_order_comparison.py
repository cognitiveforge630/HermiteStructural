"""Replot the published two-mesh Hex27/Hex64 summary without a large solve.

This chart uses RECORDed numerical summaries, not raw 3x5x25 solve archives.
To reproduce those numbers independently, run the expensive verified solver
commands in docs/HIGHER_NODE_HEX27_HEX64_STUDY_20261027.md.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'samples/hex27_hex64_recorded_mesh_summary.json'


def render(out: Path) -> Path:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    source=json.loads(SOURCE.read_text(encoding='utf-8'))
    if source['status'] != 'RECORDED_RESEARCH_SUMMARY_NOT_RAW_SOLVE_OUTPUT':
        raise ValueError('Source is not the vetted recorded summary')
    data={(c['order'], tuple(c['mesh'])): c for c in source['cases']}
    levels=((3,3,5),(3,5,25))
    fig, axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    metrics=(('tip_pct','Tip displacement'),('section_shear_pct','Midpoint section shear'),
             ('section_moment_pct','Midpoint section moment'))
    for ax, (metric, title) in zip(axes,metrics):
        for order, label, color in ((2,'Hex27','#287c9b'),(3,'Hex64','#bd6225')):
            values=[data[order,mesh][metric] for mesh in levels]
            ax.plot([0,1],values,'o-',label=label,color=color)
            for x, val in enumerate(values):
                ax.annotate(f'{val:.4g}%',(x,val),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
        ax.set(xticks=[0,1],xticklabels=['3×3×5','3×5×25'],xlabel='Vertex mesh',
               ylabel='Absolute discrepancy vs beam (%)',title=title)
        ax.set_yscale('log')
        ax.grid(axis='y',alpha=.2)
        ax.spines[['right','top']].set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle('Recorded Hex27 / Hex64 comparison — 2 meshes only (not a convergence rate)')
    fig.supxlabel('Based on saved summary values; full 3×5×25 raw solver archives are not distributed.',fontsize=9)
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=170)
    plt.close(fig)
    return out


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'docs/figures/hex27_hex64_mesh_comparison.svg')
    args=parser.parse_args(argv)
    print('Replotted saved summary:',render(args.out))


if __name__=='__main__':
    main()
