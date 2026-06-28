"""Compare coded shape-function derivatives against SymPy derivatives.

The script differentiates the original Lagrange Hex8 and Hermite basis
functions with SymPy, evaluates the symbolic derivatives and the derivatives
implemented in ``FEMHermiteBeamRegion.py``, then saves overlay plots.

Run from the repository root:

    py plot_shape_derivative_comparison.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sympy as sp
except ModuleNotFoundError as exc:
    missing = exc.name
    raise SystemExit(
        f"Missing dependency: {missing}. Install the project requirements with "
        "`py -m pip install -r requirements.txt`, then run this script again."
    ) from exc

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


NODE_SIGNS = np.array(
    [
        [-1, -1, -1],
        [1, -1, -1],
        [1, 1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
        [1, -1, 1],
        [1, 1, 1],
        [-1, 1, 1],
    ]
)


def sympy_hermite_1d(s: sp.Symbol) -> list[sp.Expr]:
    return [
        sp.Rational(1, 4) * (2 - 3 * s + s**3),
        sp.Rational(1, 4) * (1 - s - s**2 + s**3),
        sp.Rational(1, 4) * (2 + 3 * s - s**3),
        sp.Rational(1, 4) * (-1 - s + s**2 + s**3),
    ]


def sympy_lagrange_hex8(xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> list[sp.Expr]:
    return [
        sp.Rational(1, 8) * (1 - xi) * (1 - eta) * (1 - zeta),
        sp.Rational(1, 8) * (1 + xi) * (1 - eta) * (1 - zeta),
        sp.Rational(1, 8) * (1 + xi) * (1 + eta) * (1 - zeta),
        sp.Rational(1, 8) * (1 - xi) * (1 + eta) * (1 - zeta),
        sp.Rational(1, 8) * (1 - xi) * (1 - eta) * (1 + zeta),
        sp.Rational(1, 8) * (1 + xi) * (1 - eta) * (1 + zeta),
        sp.Rational(1, 8) * (1 + xi) * (1 + eta) * (1 + zeta),
        sp.Rational(1, 8) * (1 - xi) * (1 + eta) * (1 + zeta),
    ]


def axis_shape_exprs(s: sp.Symbol, sign: int) -> tuple[sp.Expr, sp.Expr]:
    H = sympy_hermite_1d(s)
    if sign < 0:
        return H[0], H[1]
    return H[2], H[3]


def sympy_hermite_hex8(
    xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol
) -> tuple[list[sp.Expr], list[list[sp.Expr]]]:
    NH = []
    RH = [[], [], []]

    for sx, sy, sz in NODE_SIGNS:
        vx, rx = axis_shape_exprs(xi, int(sx))
        vy, ry = axis_shape_exprs(eta, int(sy))
        vz, rz = axis_shape_exprs(zeta, int(sz))

        NH.append(vx * vy * vz)
        RH[0].append(rx * vy * vz)
        RH[1].append(vx * ry * vz)
        RH[2].append(vx * vy * rz)

    return NH, RH


def make_region() -> FEMHermiteBeamRegion:
    return FEMHermiteBeamRegion(
        Lx=1.0,
        Ly=1.0,
        Lz=1.0,
        nx=2,
        ny=2,
        nz=2,
        E=1.0,
        nu=0.3,
        gamma=1.0,
        w=1.0,
    )


def eval_expr(expr: sp.Expr, symbols: tuple[sp.Symbol, ...], *values: np.ndarray | float) -> np.ndarray:
    fn = sp.lambdify(symbols, expr, "numpy")
    evaluated = fn(*values)
    array = np.asarray(evaluated, dtype=float)
    if array.shape == ():
        for value in values:
            value_array = np.asarray(value, dtype=float)
            if value_array.shape != ():
                return np.full_like(value_array, float(array), dtype=float)
    return array


def plot_grid(
    output_path: Path,
    title: str,
    rows: int,
    cols: int,
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.7, rows * 2.25), squeeze=False)
    fig.suptitle(title, fontsize=14)

    for ax, (panel_title, x, sympy_y, code_y) in zip(axes.ravel(), panels):
        ax.plot(x, sympy_y, color="#1f77b4", linewidth=2.0, label="SymPy derivative")
        ax.plot(x, code_y, color="#d62728", linestyle="--", linewidth=1.6, label="Code derivative")
        ax.set_title(panel_title, fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.axhline(0.0, color="black", linewidth=0.5, alpha=0.35)

    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def constant_like(x: np.ndarray, value: float) -> np.ndarray:
    return np.full_like(x, value, dtype=float)


def plot_hermite_1d(region: FEMHermiteBeamRegion, output_dir: Path, samples: np.ndarray) -> float:
    s = sp.Symbol("s")
    H = sympy_hermite_1d(s)
    dH = [sp.diff(expr, s) for expr in H]
    code_values = np.array([region.dhermite_H_functions(float(value)) for value in samples]).T

    panels = []
    max_err = 0.0
    for index, expr in enumerate(dH):
        sympy_values = eval_expr(expr, (s,), samples)
        err = float(np.max(np.abs(sympy_values - code_values[index])))
        max_err = max(max_err, err)
        panels.append((f"dH{index + 1}/ds", samples, sympy_values, code_values[index]))

    plot_grid(
        output_dir / "hermite_1d_derivatives.png",
        "1D Hermite Derivatives: SymPy vs Code",
        rows=2,
        cols=2,
        panels=panels,
    )
    return max_err


def plot_lagrange_hex8(region: FEMHermiteBeamRegion, output_dir: Path, samples: np.ndarray) -> float:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    symbols = (xi, eta, zeta)
    variables = (xi, eta, zeta)
    axis_names = ("xi", "eta", "zeta")
    N = sympy_lagrange_hex8(xi, eta, zeta)

    panels = []
    max_err = 0.0
    zeros = constant_like(samples, 0.0)
    for axis, variable in enumerate(variables):
        args = [zeros, zeros, zeros]
        args[axis] = samples
        for node, expr in enumerate(N):
            sympy_values = eval_expr(sp.diff(expr, variable), symbols, *args)
            code_values = np.array(
                [
                    region.hex8_shape_derivatives(float(args[0][i]), float(args[1][i]), float(args[2][i]))[
                        axis, node
                    ]
                    for i in range(samples.size)
                ]
            )
            err = float(np.max(np.abs(sympy_values - code_values)))
            max_err = max(max_err, err)
            panels.append((f"dN{node + 1}/d{axis_names[axis]}", samples, sympy_values, code_values))

    plot_grid(
        output_dir / "lagrange_hex8_derivatives.png",
        "Lagrange Hex8 Derivatives on Centerline Slices",
        rows=3,
        cols=8,
        panels=panels,
    )
    return max_err


def plot_hermite_nh(region: FEMHermiteBeamRegion, output_dir: Path, samples: np.ndarray) -> float:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    symbols = (xi, eta, zeta)
    variables = (xi, eta, zeta)
    axis_names = ("xi", "eta", "zeta")
    NH, _ = sympy_hermite_hex8(xi, eta, zeta)

    panels = []
    max_err = 0.0
    zeros = constant_like(samples, 0.0)
    for axis, variable in enumerate(variables):
        args = [zeros, zeros, zeros]
        args[axis] = samples
        for node, expr in enumerate(NH):
            sympy_values = eval_expr(sp.diff(expr, variable), symbols, *args)
            code_values = np.array(
                [
                    region.get_hermite_shape_functions_and_derivatives(
                        float(args[0][i]), float(args[1][i]), float(args[2][i])
                    )[2][axis, node]
                    for i in range(samples.size)
                ]
            )
            err = float(np.max(np.abs(sympy_values - code_values)))
            max_err = max(max_err, err)
            panels.append((f"dNH{node + 1}/d{axis_names[axis]}", samples, sympy_values, code_values))

    plot_grid(
        output_dir / "hermite_NH_derivatives.png",
        "Hermite Translation Shape Derivatives on Centerline Slices",
        rows=3,
        cols=8,
        panels=panels,
    )
    return max_err


def plot_hermite_rh(region: FEMHermiteBeamRegion, output_dir: Path, samples: np.ndarray) -> float:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    symbols = (xi, eta, zeta)
    variables = (xi, eta, zeta)
    axis_names = ("xi", "eta", "zeta")
    rotation_names = ("x", "y", "z")
    _, RH = sympy_hermite_hex8(xi, eta, zeta)

    max_err = 0.0
    zeros = constant_like(samples, 0.0)
    for rotation_axis in range(3):
        panels = []
        for derivative_axis, variable in enumerate(variables):
            args = [zeros, zeros, zeros]
            args[derivative_axis] = samples
            for node, expr in enumerate(RH[rotation_axis]):
                sympy_values = eval_expr(sp.diff(expr, variable), symbols, *args)
                code_values = np.array(
                    [
                        region.get_hermite_shape_functions_and_derivatives(
                            float(args[0][i]), float(args[1][i]), float(args[2][i])
                        )[3][derivative_axis, rotation_axis, node]
                        for i in range(samples.size)
                    ]
                )
                err = float(np.max(np.abs(sympy_values - code_values)))
                max_err = max(max_err, err)
                panels.append(
                    (
                        f"dRH{rotation_names[rotation_axis]}{node + 1}/d{axis_names[derivative_axis]}",
                        samples,
                        sympy_values,
                        code_values,
                    )
                )

        plot_grid(
            output_dir / f"hermite_RH_{rotation_names[rotation_axis]}_derivatives.png",
            f"Hermite Rotation-{rotation_names[rotation_axis].upper()} Shape Derivatives",
            rows=3,
            cols=8,
            panels=panels,
        )

    return max_err


def main() -> None:
    repo_root = Path(__file__).resolve().parent
    output_dir = repo_root / "outputs" / "shape_derivative_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    region = make_region()
    samples = np.linspace(-1.0, 1.0, 201)

    errors = {
        "1D Hermite dH/ds": plot_hermite_1d(region, output_dir, samples),
        "Lagrange Hex8 dN/d(xi,eta,zeta)": plot_lagrange_hex8(region, output_dir, samples),
        "Hermite NH dNH/d(xi,eta,zeta)": plot_hermite_nh(region, output_dir, samples),
        "Hermite RH dRH/d(xi,eta,zeta)": plot_hermite_rh(region, output_dir, samples),
    }

    print("Derivative comparison complete.")
    for label, error in errors.items():
        print(f"  {label}: max abs error = {error:.3e}")
    print(f"Plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
