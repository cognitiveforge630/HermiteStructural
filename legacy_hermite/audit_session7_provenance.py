"""Trace the approximate 43% figure to an existing reproducible result."""
import argparse
import json
from pathlib import Path
import numpy as np


def audit(root):
    path = root/'outputs/session6_regions/shear_historical.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    region = data['regional_true_errors'][0]
    retained = region['retained']
    computed = 100*retained['difference_energy_norm']/retained['reference_energy_norm']
    if not np.isclose(computed,retained['relative_energy_difference_percent'],rtol=1e-12):
        raise ArithmeticError('Stored percentage does not match its numerator and denominator')
    return {'source':str(path.resolve()),'reported_approximate_percent':43.,
            'nearest_traceable_value_percent':computed,
            'classification':'TRUE_REGIONAL_MANUFACTURED_SHEAR_ERROR',
            'clamp':'historical_zero_spin','mesh':[2,2,2],
            'region_z_range_in':region['retained_z_range_in'],
            'numerator':retained['difference_energy_norm'],
            'denominator':retained['reference_energy_norm'],
            'denominator_definition':'Exact shear solution energy norm over the SAME retained region',
            'formula':'100 * true_regional_error_norm / exact_regional_solution_norm',
            'cantilever_error_established':False,'a_posteriori_estimator_value':False,
            'qualification':'This traces the saved 42.9034% result; no separately measured 43% cantilever error was found.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('outputs/session7_audit/provenance.json'))
    args=parser.parse_args()
    result=audit(Path(__file__).resolve().parent)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
