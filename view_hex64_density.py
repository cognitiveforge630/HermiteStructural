"""Compatibility entry point for Hex27/Hex64 mesh density previews.

Prefer `python -m tools.mesh_preview --family hex64` for new work.
"""
from __future__ import annotations

import argparse

from tools.mesh_preview import main as preview_main


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, choices=(2, 3), default=3)
    parser.add_argument('--nx', type=int, default=3)
    parser.add_argument('--ny', type=int, default=5)
    parser.add_argument('--nz', type=int, default=25)
    parser.add_argument('--nodes', action='store_true')
    parser.add_argument('--gauss-points', action='store_true')
    parser.add_argument('--quadrature-order', type=int)
    parser.add_argument('--summary-only', action='store_true')
    args = parser.parse_args()
    forwarded = ['--family', f'hex{(args.order+1)**3}', '--nx', str(args.nx),
                 '--ny', str(args.ny), '--nz', str(args.nz)]
    if args.nodes:
        forwarded.append('--nodes')
    if args.gauss_points:
        forwarded.append('--gauss-points')
    if args.quadrature_order is not None:
        forwarded.extend(['--quadrature-order', str(args.quadrature_order)])
    if args.summary_only:
        forwarded.append('--summary-only')
    preview_main(forwarded)


if __name__ == '__main__':
    main()
