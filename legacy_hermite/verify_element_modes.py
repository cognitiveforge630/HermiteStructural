"""Independent affine/rigid-body strain reproduction diagnostic.

Do not replace nodal Hermite interpolation with a different element. The
reported failures represent properties of the current implementation.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from check_hermite_element import build_single_element
from solver import build_beam


def evaluate_modes(tolerance=1e-10, amplitude=1e-4, formulation="legacy_hermite"):
    if formulation == "legacy_hermite":
        beam = build_single_element()
    elif formulation in ("cauchy_hex8_reference", "affine_enriched_hex8",
                         "recovered_rotation_hermite", "quadratic_rotation_hermite"):
        n = 3 if formulation == "quadratic_rotation_hermite" else 2
        beam = build_beam(nx=n, ny=n, nz=(5 if n == 3 else 2), formulation=formulation)
    else:
        raise ValueError("Unknown formulation")
    coords = beam.get_element_nodes(0)
    nonlocal_stencil=hasattr(beam,'element_transform')
    pts = ([(-.75,.4,-.1),(.3,-.5,.7),(0.,0.,0.)] if nonlocal_stencil
           else [(x, y, z) for x in beam.gauss_points
                 for y in beam.gauss_points for z in beam.gauss_points])
    # Engineering strain uses gamma_ij = grad_i u_j + grad_j u_i.
    def check(name, gradient, offset=None, omega=None):
        G = np.asarray(gradient, dtype=float) * amplitude
        offset = np.zeros(3) if offset is None else np.asarray(offset) * amplitude
        omega = np.zeros(3) if omega is None else np.asarray(omega) * amplitude
        ndof = beam.ndof_per_node
        if nonlocal_stencil:
            U = np.zeros(beam.grid_points.shape[0] * ndof)
            for n,xyz in enumerate(beam.grid_points):
                U[n*ndof:n*ndof+3] = offset + G @ xyz
                U[n*ndof+3:n*ndof+6] = omega
        else:
            U = np.zeros(8 * ndof)
            for n, xyz in enumerate(coords):
                U[n*ndof:n*ndof+3] = offset + G @ xyz
                if ndof == 6 and formulation == "legacy_hermite":
                    U[n*ndof+3:n*ndof+6] = omega
        expected = np.array([G[0,0], G[1,1], G[2,2],
                             G[0,1]+G[1,0], G[1,2]+G[2,1], G[0,2]+G[2,0]])
        errors = []
        for e in (range(beam.get_num_elements()) if nonlocal_stencil else [0]):
            for p in pts:
                if nonlocal_stencil:
                    measured=beam.strain_at(e,U,p)
                else:
                    B, *_ = beam.build_B_matrix(*p, coords)
                    measured=B @ U
                errors.append(np.max(np.abs(measured-expected)))
        error = float(max(errors))
        return {'mode':name, 'max_absolute_strain_error':error,
                'tolerance':tolerance, 'pass':error <= tolerance}
    modes = [check('translation_x', np.zeros((3,3)), [1,0,0]),
             check('translation_y', np.zeros((3,3)), [0,1,0]),
             check('translation_z', np.zeros((3,3)), [0,0,1])]
    for a in range(3):
        w = np.zeros(3); w[a] = 1
        cross = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        modes.append(check(f'rigid_rotation_{"xyz"[a]}', cross, omega=w))
    for a in range(3):
        G=np.zeros((3,3)); G[a,a]=1
        modes.append(check(f'extension_{"xyz"[a]}', G))
    for a,b in [(0,1),(1,2),(0,2)]:
        G=np.zeros((3,3));G[a,b]=G[b,a]=0.5
        modes.append(check(f'symmetric_shear_{"xyz"[a]}{"xyz"[b]}', G))
    return {'method':'independent_single_element_gauss_point_affine_reproduction',
            'auxiliary_dofs_are_physical_rotations': formulation in ('legacy_hermite','recovered_rotation_hermite','quadratic_rotation_hermite'),
            'formulation': formulation,
            'amplitude':amplitude, 'modes':modes,
            'pass_count':sum(m['pass'] for m in modes),
            'fail_count':sum(not m['pass'] for m in modes),
            'overall_pass':all(m['pass'] for m in modes)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formulation', choices=('legacy_hermite', 'cauchy_hex8_reference', 'affine_enriched_hex8', 'recovered_rotation_hermite', 'quadratic_rotation_hermite'),
                        default='legacy_hermite')
    parser.add_argument('--report-only', action='store_true', help='Return exit status 0 even if some modes fail')
    parser.add_argument('--out', type=Path, default=Path('outputs/element_modes_baseline.json'))
    args=parser.parse_args()
    result=evaluate_modes(formulation=args.formulation)
    if args.formulation == "affine_enriched_hex8":
        print("NOTE: DOFs 3:6 are interior bubble amplitudes, NOT nodal rotations.")
    for m in result['modes']:
        status='PASS' if m['pass'] else 'FAIL'
        print(f'{status:4} {m["mode"]:28} max strain error={m["max_absolute_strain_error"]:.9e}')
    print(f'OVERALL {result["pass_count"]} PASS / {result["fail_count"]} FAIL')
    if args.formulation != 'legacy_hermite' and args.out == Path('outputs/element_modes_baseline.json'):
        args.out = Path(f'outputs/element_modes_{args.formulation}.json')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if not args.report_only and not result['overall_pass']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
