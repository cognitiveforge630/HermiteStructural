"""Quick nodal interpolation/recovery audits; no global factorization."""
import json
import numpy as np
from core.solvers.beam import build_beam


def main():
    report = {}
    for degree in ("linear", "quadratic", "cubic"):
        beam = build_beam(nx=3, ny=3, nz=5, formulation=f"{degree}_rotation_hermite")
        transform_dofs, T = beam.element_transform(0)
        coords = beam.get_element_nodes(0)
        B, jac, N = beam.primitive_B(0., 0., 0., coords)
        patch = beam.displacement_operator(0, (0.,0.,0.))[1]
        report[degree] = {
            "nodes": int(len(beam.grid_points)), "elements": int(beam.get_num_elements()),
            "element_support_dofs": int(len(transform_dofs)), "primitive_dofs": int(beam.primitive_dofs),
            "positive_jacobian": bool(jac > 0), "finite_operator": bool(np.all(np.isfinite(B@T)) and np.all(np.isfinite(patch))),
            "quadratic_patch_min_rank": int(min(beam.quadratic_rank)) if hasattr(beam,"quadratic_rank") else None,
            "cubic_patch_min_rank": int(min(beam.cubic_rank)) if hasattr(beam,"cubic_rank") else None,
        }
        assert report[degree]["positive_jacobian"] and report[degree]["finite_operator"]
    print(json.dumps(report, indent=2))
    print("Research checks only: no convergence or error-bound certification")

if __name__ == "__main__":
    main()
