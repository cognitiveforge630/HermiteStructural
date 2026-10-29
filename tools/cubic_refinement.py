"""Compare saved cubic beam runs without re-solving or claiming certification.

Reports beam-reference discrepancies and mesh-to-mesh scalar QoI changes. This
is a *convergence diagnostic*, NOT an a posteriori error estimator or bound.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _finite_number(value, name):
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _audit_location(result_file: Path, content: dict) -> Path:
    expected = content.get('cubic_stress_audit_file')
    candidates = [result_file.with_name(result_file.stem + '_stress_audit.json')]
    if expected:
        ref = Path(expected)
        candidates += [ref, result_file.parent / ref, ROOT / ref,
                       result_file.parent / ref.name]
    for p in candidates:
        if p.is_file():
            return p
    raise FileNotFoundError(f'Missing matching cubic stress audit for {result_file}; looked in {candidates}')


def read_case(filename):
    file = Path(filename)
    data = json.loads(file.read_text(encoding='utf-8'))
    if data.get('formulation') != 'cubic_rotation_hermite':
        raise ValueError(f'{file}: only cubic_rotation_hermite results are comparable')
    dims = data.get('node_mesh')
    if not isinstance(dims, list) or len(dims) != 3 or any(type(n) is not int or n < 3 for n in dims):
        raise ValueError(f'{file}: invalid node_mesh')
    audit_file = _audit_location(file, data)
    audit = json.loads(audit_file.read_text(encoding='utf-8'))
    if audit.get('stress_distribution_certified') is not False:
        raise ValueError(f'{audit_file}: unexpected certification status')
    if not isinstance(audit.get('section_resultants'), list):
        raise ValueError(f'{audit_file}: no section resultants')
    if not math.isclose(_finite_number(audit['material_E_longitudinal_psi'], 'audit E'),
                        _finite_number(data['longitudinal_E_psi'], 'result E'),
                        rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(f'{file}: stress audit modulus does not match the solution')
    n_sections = len(audit['section_resultants'])
    num_z_elements = dims[2] - 1
    sampled_n = audit.get('sample_subdivisions')
    if sampled_n is None:
        # Shipped October 2026 audits predate explicit sample_subdivisions.
        # Each structured z-element yields one section per sample subcell.
        if n_sections % num_z_elements != 0:
            raise ValueError(f'{audit_file}: cannot infer original sampling density')
        sampled_n = n_sections // num_z_elements
    if type(sampled_n) is not int or sampled_n < 1 or n_sections != num_z_elements * sampled_n:
        raise ValueError(f'{audit_file}: inconsistent longitudinal sampling metadata')
    tips = _finite_number(data['fem_tip_uy_in'], 'tip')
    beam_tip = _finite_number(data['awc_reference']['total_uy_in'], 'beam reference tip')
    if not beam_tip:
        raise ValueError(f'{file}: zero beam reference tip')
    stress = _finite_number(audit['quarter_to_threequarter_span_sigma_zz_relative_L2_pct'], 'stress difference')
    return {
        'result_file': str(file), 'audit_file': str(audit_file),
        'mesh': dims, 'subdivisions': sampled_n,
        'tip_uy_in': tips, 'beam_tip_uy_in': beam_tip,
        'tip_reference_discrepancy_pct': 100 * abs(tips-beam_tip)/abs(beam_tip),
        'stress_reference_discrepancy_pct': stress,
        'moment_reference_discrepancy_max_pct': _finite_number(audit['quarter_to_threequarter_span_max_moment_error_pct'], 'moment discrepancy'),
        'shear_reference_discrepancy_max_pct': _finite_number(audit['quarter_to_threequarter_span_max_shear_error_pct'], 'shear discrepancy'),
        'longitudinal_tension_sampled_peak_psi': _finite_number(audit['sampled_stress_peaks']['longitudinal_tension_psi'], 'peak tension'),
        'end_force_lbf': _finite_number(data['end_force_lbf'], 'end force'),
        'total_y_force_lbf': _finite_number(data['total_y_force_lbf'], 'total force'),
        'longitudinal_E_psi': _finite_number(data['longitudinal_E_psi'], 'E'),
        'weight_per_length_lbf_in': _finite_number(data['awc_reference']['weight_per_length_lbf_in'], 'beam weight per length'),
        'boundary_condition': data['boundary_condition'],
    }


def compare(cases):
    if not cases:
        raise ValueError('Provide at least one case')
    primary = cases[0]
    for case in cases[1:]:
        for key in ('end_force_lbf', 'total_y_force_lbf', 'longitudinal_E_psi',
                    'weight_per_length_lbf_in'):
            if not math.isclose(case[key], primary[key], rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f'Incompatible cases: {key} differs between {case["result_file"]} and {primary["result_file"]}')
        if case['boundary_condition'] != primary['boundary_condition']:
            raise ValueError('Incompatible cases: boundary conditions differ')
        if case['beam_tip_uy_in'] != primary['beam_tip_uy_in']:
            raise ValueError('Incompatible cases: analytical reference differs')
    pairs = []
    for older, newer in zip(cases, cases[1:]):
        if not all(b >= a for a, b in zip(older['mesh'], newer['mesh'])) or older['mesh'] == newer['mesh']:
            raise ValueError('Provide distinct meshes in nondecreasing nx,ny,nz order')
        ratios = [(a-1)/(b-1) for a,b in zip(older['mesh'],newer['mesh'])]
        uniform = max(ratios)-min(ratios) < 1.e-12 and ratios[0] < 1
        equal_sampling = older['subdivisions'] == newer['subdivisions']
        pairs.append({
            'from_mesh': older['mesh'], 'to_mesh': newer['mesh'],
            'relative_element_size_ratios_xyz': ratios,
            'uniform_h_refinement': uniform,
            'matched_sampling_subdivisions': equal_sampling,
            'tip_change_pct_of_fine': 100*abs(newer['tip_uy_in']-older['tip_uy_in'])/max(abs(newer['tip_uy_in']),1e-30),
            'sampled_peak_tension_change_pct_of_fine': 100*abs(newer['longitudinal_tension_sampled_peak_psi']-older['longitudinal_tension_sampled_peak_psi'])/max(abs(newer['longitudinal_tension_sampled_peak_psi']),1e-30),
            'caution': ('Matched sampling; still not a stress-field norm or certified bound.' if equal_sampling
                        else 'Sampling densities DIFFER; stress peak/QoI changes are confounded. Re-run with matched --subdivisions.'),
        })
    return {
        'status': 'RESEARCH_MESH_COMPARISON_NOT_CERTIFIED',
        'certified_bound': False,
        'interpretation': 'Beam-reference discrepancies and intermesh changes of scalar QoIs only; no guaranteed true-solution error bound.',
        'cases': cases, 'pairs': pairs,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, nargs='+', help='Cubic result JSONs ordered coarse-to-fine')
    parser.add_argument('--json-out', type=Path, default=None, help='Save a machine-readable comparison report')
    a = parser.parse_args(argv)
    result = compare([read_case(p) for p in a.results])
    print('RESEARCH COMPARISON — NOT A CERTIFIED A POSTERIORI ESTIMATE')
    print('Mesh        sample   tip vs beam %   sigma_zz vs beam L2 %   shear resultant max %')
    for case in result['cases']:
        print(f"{'x'.join(map(str,case['mesh'])):11s} {case['subdivisions']:6d} "
              f"{case['tip_reference_discrepancy_pct']:15.4f} "
              f"{case['stress_reference_discrepancy_pct']:22.4f} "
              f"{case['shear_reference_discrepancy_max_pct']:23.4f}")
    for pair in result['pairs']:
        print(f"{pair['from_mesh']} -> {pair['to_mesh']}: "
              f"tip intermesh change {pair['tip_change_pct_of_fine']:.4f}%; "
              f"peak sampled tension change {pair['sampled_peak_tension_change_pct_of_fine']:.4f}%; "
              f"uniform h: {pair['uniform_h_refinement']}; {pair['caution']}")
    if a.json_out:
        a.json_out.parent.mkdir(parents=True,exist_ok=True)
        a.json_out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'Wrote {a.json_out}')
    return result


if __name__ == '__main__':
    main()
