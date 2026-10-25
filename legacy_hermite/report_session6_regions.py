"""Global and regional energy reports from saved solutions, without re-solving.

Manufactured shear/sine archives have known exact strain fields. Cantilever
archives require a reference archive; reference uncertainty remains unqualified.
Regional percentages describe strain energy, not tip-displacement errors.
"""
import argparse
import json
from pathlib import Path

from compare_session6_fields import read_solution, compare_fields
from run_session6_accuracy import exact_fields, exact_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('solution', type=Path)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--exclude-clamp-lengths', nargs='+', type=float, required=True,
                        help='Prespecified physical exclusion distances in inches')
    parser.add_argument('--order', type=int, default=6)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.order < 4:
        parser.error('Use quadrature order >=4')
    left = read_solution(args.solution)
    if args.reference is not None:
        result = compare_fields(left, read_solution(args.reference), args.order,
                                exclude_clamp_lengths=args.exclude_clamp_lengths)
        result['reference'] = str(args.reference.resolve())
    else:
        beam, U, meta = left
        if meta['problem'] not in ('shear', 'sine'):
            parser.error('Cantilever has no known exact field; supply --reference')
        strain = exact_fields(beam, meta['problem'])[1]
        result = exact_error(beam, U, strain, args.order,
                             exclude_clamp_lengths=args.exclude_clamp_lengths)
        result.update(status='TRUE_MANUFACTURED_ENERGY_ERROR', certified_bound=False,
                      a_posteriori_estimator_available=False, quadrature_order=args.order)
    result['solution'] = str(args.solution.resolve())
    result['regional_metric_definition'] = (
        '100*sqrt(integral_region(strain_difference^T D strain_difference))'
        '/sqrt(integral_region(reference_strain^T D reference_strain))')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
