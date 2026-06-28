"""Orthotropic Redwood cantilever study with AWC-target E calibration.

The beam axis is the finite-element z-axis, so the wood longitudinal material
axis L is mapped to z. The radial/transverse directions are mapped as
R -> x and T -> y.

The calibration case solves the finite-element model once with a seed
longitudinal modulus, then uses the linear-elastic scaling law

    U(E_calibrated) = U(seed_E) * seed_E / E_calibrated

so the reported calibrated ``Ez`` lands exactly on the requested AWC/theory
free-end displacement target. A verification solve can be requested with
``--verify-calibrated``.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from FEMHermiteBeamRegion import FEMHermiteBeamRegion, OrthotropicElasticity
from checkme import E, I, L, P, cantilever_point_load_deflection
from solver import LOAD_DOF, build_end_face_traction_load_vector, solve_with_force_vector


QUICK_MESH = (3, 5, 31)
CURRENT_MESH = (5, 9, 121)
HALF_INCH_MESH = (9, 17, 241)
AWC_REDWOOD_SELECT_STRUCTURAL_E = 1_100_000.0
DEFAULT_OUTPUT = Path("outputs/redwood_awc_calibration_summary.json")
THREE_POINT_GAUSS_ORDER = 3


REDWOOD_CLEAR_WOOD_RATIOS = {
    # Approximate Redwood clear-wood orthotropic ratios. The longitudinal
    # modulus is the scalar we intentionally slacken/tune in this study.
    "E_T_over_E_L": 0.089,
    "E_R_over_E_L": 0.087,
    "G_LR_over_E_L": 0.066,
    "G_LT_over_E_L": 0.077,
    "G_RT_over_E_L": 0.011,
    # Major Poisson ratios in common wood notation. For example, nu_LR is the
    # radial strain produced by longitudinal stress. The helper below converts
    # these to the minor ratios required by the FE x/y/z compliance convention.
    "nu_LR": 0.360,
    "nu_LT": 0.346,
    "nu_RT": 0.373,
}


def awc_theory_target_uy(awc_E=AWC_REDWOOD_SELECT_STRUCTURAL_E):
    """Beam-table target deflection for the selected AWC longitudinal E."""
    return cantilever_point_load_deflection(P, L, awc_E, I)


def redwood_orthotropic_material(E_longitudinal=E):
    """Return a Redwood-like orthotropic material with L along global z.

    Mapping used by the Hermite model:

    * x = R, radial direction
    * y = T, tangential/transverse direction
    * z = L, longitudinal beam direction

    ``OrthotropicElasticity`` expects ``nuxz`` to mean z-strain from x-stress.
    Published wood constants often give the reciprocal major value ``nu_LR``
    instead: x-strain from z-stress. Reciprocity gives
    ``nu_RL / E_R = nu_LR / E_L``, so ``nu_RL = nu_LR * E_R / E_L``.
    The same conversion is applied to the T/L pair.
    """
    E_longitudinal = float(E_longitudinal)
    ratios = REDWOOD_CLEAR_WOOD_RATIOS

    Ex = ratios["E_R_over_E_L"] * E_longitudinal
    Ey = ratios["E_T_over_E_L"] * E_longitudinal
    Ez = E_longitudinal

    return OrthotropicElasticity(
        Ex=Ex,
        Ey=Ey,
        Ez=Ez,
        nuxy=ratios["nu_RT"],
        nuxz=ratios["nu_LR"] * Ex / Ez,
        nuyz=ratios["nu_LT"] * Ey / Ez,
        Gxy=ratios["G_RT_over_E_L"] * E_longitudinal,
        Gyz=ratios["G_LT_over_E_L"] * E_longitudinal,
        Gxz=ratios["G_LR_over_E_L"] * E_longitudinal,
    )


# Backward-compatible name used by earlier notes/scripts.
def redwood_orthotropic_constants(E_longitudinal=E):
    return redwood_orthotropic_material(E_longitudinal).as_dict()


def build_beam(mesh, gauss_order, material, E_longitudinal):
    return FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=120.0,
        nx=mesh[0],
        ny=mesh[1],
        nz=mesh[2],
        E=E_longitudinal,
        nu=0.3,
        gamma=0.0,
        w=0.0,
        gauss_order=gauss_order,
        orthotropic_constants=material,
    )


def mesh_stats(mesh):
    nx, ny, nz = mesh
    nodes = nx * ny * nz
    elements = (nx - 1) * (ny - 1) * (nz - 1)
    dofs = nodes * 6
    return nodes, elements, dofs


def estimate_assembly_storage_mb(mesh):
    _, elements, _ = mesh_stats(mesh)
    entries = elements * 48 * 48
    return entries * (8 + 8 + 8) / 1024**2


def material_for_name(material_name, E_longitudinal):
    if material_name == "isotropic":
        return None
    if material_name == "redwood":
        return redwood_orthotropic_material(E_longitudinal)
    raise ValueError(f"Unknown material name: {material_name}")


def run_case(
    name,
    mesh,
    gauss_order,
    material_name,
    E_longitudinal,
    target_uy=None,
    target_label=None,
):
    material = material_for_name(material_name, E_longitudinal)
    nodes, elements, dofs = mesh_stats(mesh)
    beam = build_beam(mesh, gauss_order, material, E_longitudinal)
    force, _, _ = build_end_face_traction_load_vector(beam)

    start = time.perf_counter()
    U = solve_with_force_vector(beam, force)
    elapsed = time.perf_counter() - start

    U_by_node = U.reshape((-1, beam.ndof_per_node))
    end_nodes = beam.get_nodes_at_max_z()
    end_uy = U_by_node[end_nodes, LOAD_DOF]
    fea_end_uy = end_uy[np.argmax(np.abs(end_uy))]
    max_disp = np.linalg.norm(U_by_node[:, :3], axis=1).max()

    current_theory = cantilever_point_load_deflection(P, L, E, I)
    material_theory = cantilever_point_load_deflection(P, L, E_longitudinal, I)
    if target_uy is None:
        target_uy = material_theory
        target_label = target_label or "same-E beam theory"
    else:
        target_label = target_label or "requested target"

    current_diff = percent_difference(fea_end_uy, current_theory)
    material_diff = percent_difference(fea_end_uy, material_theory)
    target_diff = percent_difference(fea_end_uy, target_uy)

    result = {
        "name": name,
        "material": material_name,
        "gauss": gauss_order,
        "mesh": "x".join(str(v) for v in mesh),
        "nodes": nodes,
        "elements": elements,
        "dofs": dofs,
        "assembly_storage_mb": estimate_assembly_storage_mb(mesh),
        "E_longitudinal": float(E_longitudinal),
        "fea_end_uy": float(fea_end_uy),
        "max_disp": float(max_disp),
        "current_theory": float(current_theory),
        "current_diff_percent": float(current_diff),
        "material_theory": float(material_theory),
        "material_diff_percent": float(material_diff),
        "target_uy": float(target_uy),
        "target_label": target_label,
        "target_diff_percent": float(target_diff),
        "elapsed_s": float(elapsed),
    }
    if isinstance(material, OrthotropicElasticity):
        result["orthotropic_constants"] = material.as_dict()
    return result


def percent_difference(value, target):
    if target == 0.0:
        return np.nan
    return abs(value - target) / abs(target) * 100.0


def calibrated_longitudinal_E(seed_E, seed_fea_uy, target_uy):
    """Return the Ez value that makes the linearly scaled FEA uy equal target."""
    seed_E = float(seed_E)
    seed_fea_uy = float(seed_fea_uy)
    target_uy = float(target_uy)
    if seed_E <= 0.0:
        raise ValueError("seed_E must be positive")
    if target_uy == 0.0:
        raise ValueError("target_uy must be nonzero")
    ratio = seed_fea_uy / target_uy
    if ratio <= 0.0:
        raise ValueError(
            "seed FEA displacement and target displacement must have the same sign "
            f"for positive E calibration; got seed={seed_fea_uy}, target={target_uy}"
        )
    return seed_E * ratio


def run_awc_calibration(
    mesh=QUICK_MESH,
    gauss_order=THREE_POINT_GAUSS_ORDER,
    seed_E=AWC_REDWOOD_SELECT_STRUCTURAL_E,
    target_uy=None,
    verify_calibrated=False,
):
    target_uy = awc_theory_target_uy() if target_uy is None else float(target_uy)
    target_label = "AWC/theory free-end uy"

    seed = run_case(
        "Redwood orthotropic calibration seed",
        mesh,
        gauss_order,
        "redwood",
        seed_E,
        target_uy=target_uy,
        target_label=target_label,
    )
    calibrated_E = calibrated_longitudinal_E(seed_E, seed["fea_end_uy"], target_uy)

    scaled_prediction_uy = seed["fea_end_uy"] * seed_E / calibrated_E
    calibrated = {
        "name": "Redwood orthotropic calibrated to AWC/theory target",
        "material": "redwood",
        "gauss": gauss_order,
        "mesh": seed["mesh"],
        "nodes": seed["nodes"],
        "elements": seed["elements"],
        "dofs": seed["dofs"],
        "assembly_storage_mb": seed["assembly_storage_mb"],
        "seed_E_longitudinal": float(seed_E),
        "E_longitudinal": float(calibrated_E),
        "current_theory": float(cantilever_point_load_deflection(P, L, E, I)),
        "material_theory": float(cantilever_point_load_deflection(P, L, calibrated_E, I)),
        "target_uy": float(target_uy),
        "target_label": target_label,
        "scaled_prediction_uy": float(scaled_prediction_uy),
        "scaled_prediction_diff_percent": float(percent_difference(scaled_prediction_uy, target_uy)),
        "orthotropic_constants": redwood_orthotropic_material(calibrated_E).as_dict(),
        "verification_solve_ran": False,
    }

    if verify_calibrated:
        verified = run_case(
            "Redwood orthotropic calibrated verification solve",
            mesh,
            gauss_order,
            "redwood",
            calibrated_E,
            target_uy=target_uy,
            target_label=target_label,
        )
        calibrated["verification_solve_ran"] = True
        calibrated["verification_fea_end_uy"] = verified["fea_end_uy"]
        calibrated["verification_diff_percent"] = verified["target_diff_percent"]
        calibrated["verification_elapsed_s"] = verified["elapsed_s"]

    return {"seed": seed, "calibrated": calibrated}


def print_result(result):
    print(f"{result['name']}")
    print(f"  material:                  {result['material']}")
    print(f"  mesh nodes:                {result['mesh']} ({result['nodes']} nodes)")
    print(f"  elements / dofs:           {result['elements']} / {result['dofs']}")
    print(f"  Gauss order:               {result['gauss']}")
    print(f"  assembly arrays estimate:  {result['assembly_storage_mb']:.1f} MB")
    print(f"  E_longitudinal / Ez:       {result['E_longitudinal']:.3f} psi")
    if "fea_end_uy" in result:
        print(f"  FEA free-end uy:           {result['fea_end_uy']:.6e} in")
    if "scaled_prediction_uy" in result:
        print(f"  scaled FEA uy prediction:  {result['scaled_prediction_uy']:.6e} in")
    print(f"  current theory uy:         {result.get('current_theory', float('nan')):.6e} in")
    if "current_diff_percent" in result:
        print(f"  diff vs current theory:    {result['current_diff_percent']:.2f}%")
    print(f"  material theory uy:        {result.get('material_theory', float('nan')):.6e} in")
    if "material_diff_percent" in result:
        print(f"  diff vs material theory:   {result['material_diff_percent']:.2f}%")
    print(f"  target ({result['target_label']}): {result['target_uy']:.6e} in")
    if "target_diff_percent" in result:
        print(f"  diff vs target:            {result['target_diff_percent']:.6f}%")
    if "scaled_prediction_diff_percent" in result:
        print(f"  scaled diff vs target:     {result['scaled_prediction_diff_percent']:.6e}%")
    if result.get("verification_solve_ran"):
        print(f"  verification FEA uy:       {result['verification_fea_end_uy']:.6e} in")
        print(f"  verification diff:         {result['verification_diff_percent']:.6e}%")
    if "elapsed_s" in result:
        print(f"  elapsed solve path:        {result['elapsed_s']:.1f} s")
    elif "verification_elapsed_s" in result:
        print(f"  verification elapsed:      {result['verification_elapsed_s']:.1f} s")
    print()


def parse_mesh(text):
    parts = text.lower().replace("×", "x").split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("mesh must look like 5x9x121")
    try:
        mesh = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("mesh entries must be integers") from exc
    if any(value < 2 for value in mesh):
        raise argparse.ArgumentTypeError("each mesh count must be at least 2")
    return mesh


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        choices=[
            "baseline",
            "isotropic-2pt-reference",
            "isotropic-3pt",
            "redwood-same-e",
            "redwood-awc-e",
            "redwood-awc-calibrated",
            "redwood-same-e-half-inch",
        ],
        help="Case to run. Repeat to run multiple cases.",
    )
    parser.add_argument(
        "--mesh",
        type=parse_mesh,
        default=QUICK_MESH,
        help="Mesh node counts as nxXnyXnz. Default: 3x5x31 quick mesh; use 5x9x121 for the current full mesh.",
    )
    parser.add_argument(
        "--awc-e",
        type=float,
        default=AWC_REDWOOD_SELECT_STRUCTURAL_E,
        help="AWC/reference longitudinal E used for the default target deflection.",
    )
    parser.add_argument(
        "--target-uy",
        type=float,
        default=None,
        help="Optional explicit AWC/theory free-end uy target in inches.",
    )
    parser.add_argument(
        "--seed-e",
        type=float,
        default=None,
        help="Seed longitudinal E for calibration. Default: --awc-e.",
    )
    parser.add_argument(
        "--verify-calibrated",
        action="store_true",
        help="Run a second solve at calibrated E to verify the exact match.",
    )
    parser.add_argument(
        "--write-json",
        nargs="?",
        const=str(DEFAULT_OUTPUT),
        default=None,
        help="Write case results to JSON. Optional path defaults to outputs/redwood_awc_calibration_summary.json.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    case_names = args.case or ["baseline", "redwood-awc-e", "redwood-awc-calibrated"]
    seed_E = args.awc_e if args.seed_e is None else args.seed_e
    target_uy = args.target_uy
    if target_uy is None:
        target_uy = awc_theory_target_uy(args.awc_e)

    results = []
    cases = {
        "baseline": ("baseline isotropic 3pt", args.mesh, THREE_POINT_GAUSS_ORDER, "isotropic", E),
        "isotropic-2pt-reference": ("isotropic 2pt reference", args.mesh, 2, "isotropic", E),
        "isotropic-3pt": ("isotropic 3pt", args.mesh, THREE_POINT_GAUSS_ORDER, "isotropic", E),
        "redwood-same-e": ("Redwood orthotropic 3pt, same E", args.mesh, THREE_POINT_GAUSS_ORDER, "redwood", E),
        "redwood-awc-e": (
            "Redwood orthotropic 3pt, AWC/reference E",
            args.mesh,
            THREE_POINT_GAUSS_ORDER,
            "redwood",
            args.awc_e,
        ),
        "redwood-same-e-half-inch": (
            "Redwood orthotropic 3pt, same E, half-inch mesh",
            HALF_INCH_MESH,
            THREE_POINT_GAUSS_ORDER,
            "redwood",
            E,
        ),
    }

    for case_name in case_names:
        if case_name == "redwood-awc-calibrated":
            calibration = run_awc_calibration(
                mesh=args.mesh,
                gauss_order=THREE_POINT_GAUSS_ORDER,
                seed_E=seed_E,
                target_uy=target_uy,
                verify_calibrated=args.verify_calibrated,
            )
            print_result(calibration["seed"])
            print_result(calibration["calibrated"])
            results.append(calibration)
            continue

        result = run_case(
            *cases[case_name],
            target_uy=target_uy,
            target_label="AWC/theory free-end uy",
        )
        print_result(result)
        results.append(result)

    if args.write_json:
        output = Path(args.write_json)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(results, indent=2))
        print(f"Saved JSON summary to {output}")


if __name__ == "__main__":
    main()
