# HermiteStructural

Finite element tutorials and source code for experimenting with a Hermite
hexahedral structural model in Python.

The current assembly follows a Hermitian Hexa8-style formulation with six
degrees of freedom per node. Translation DOFs use `NH` shape functions, rotation
DOFs use `RH` shape functions, and stiffness is assembled from the standard
six-component small-strain matrix. Rotational DOFs are coupled into
displacement through the small-rotation relation `u = theta x r`.

The element now uses three-point Gauss integration by default. The
orthotropic Redwood/AWC calibration path intentionally keeps this full
integration order on so the orthotropic stiffness matrix `D` is integrated
consistently. Use the explicit two-point comparison case only when you are
studying integration sensitivity.

This repository starts with a cantilever beam example, a theoretical point-load
check, a solver, and PyVista visualization tools. The goal is to keep the math,
code, and tutorial writeups close together so each result can be inspected
locally.

## Setup

If Python or Visual Studio Code are not set up yet, start here:

[Installing Python on Windows with Python Manager and Visual Studio Code](https://cognitiveforge630.github.io/cognitiveforge.github.io/2026/05/10/install-python-windows-vscode.html)

Then create and activate a virtual environment from this repository folder:

```powershell
py -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
.\.venv\Scripts\Activate.ps1
```

Install dependencies:

```powershell
py -m pip install -r requirements.txt
```

If `py` is not available, use the Python executable from your virtual
environment instead.

## First Checks

Run the theoretical cantilever point-load check:

```powershell
py checkme.py
```

Run the finite element solver:

```powershell
py solver.py
```

View the saved displacement result:

```powershell
py view.py
```

The viewer shows the solved beam, Dirichlet fixed nodes at `z = 0`, displacement
magnitude, and a small overlay comparing theoretical and actual free-end
displacement.


## Orthotropic Redwood / AWC Calibration

The project now includes an orthotropic 3D elasticity matrix for material axes
`x, y, z` with stress/strain order `xx, yy, zz, xy, yz, xz`. The helper class is
`OrthotropicElasticity` in `FEMHermiteBeamRegion.py`; it builds the compliance
matrix first, checks that it is positive definite, and then inverts it to get
the stiffness matrix `D`. The AWC/Redwood path uses three-point Gauss
integration by default and the check script verifies that this stays turned on.

For the Redwood study, `redwood_orthotropic_study.py` maps the wood axes as
`R -> x`, `T -> y`, and `L -> z`, so the beam-length direction uses the
longitudinal modulus `Ez`. The script can then slacken/tune that longitudinal
`E` so the finite-element free-end displacement lands exactly on the requested
AWC/theory displacement target.

Fast matrix check:

```powershell
py check_orthotropic_elasticity.py
```

Quick calibration smoke test on a tiny mesh; this uses 3-point Gauss integration:

```powershell
py redwood_orthotropic_study.py --mesh 2x2x3 --case redwood-awc-calibrated --verify-calibrated
```

Full current cantilever mesh calibration:

```powershell
py redwood_orthotropic_study.py --mesh 5x9x121 --case redwood-awc-calibrated --write-json
```

To match a specific AWC displacement number directly, pass it in inches:

```powershell
py redwood_orthotropic_study.py --mesh 5x9x121 --case redwood-awc-calibrated --target-uy -3.927272727e-2 --write-json
```

The calibration uses linear-elastic scaling. If a seed solve gives
`u_seed` at `E_seed`, the exact matching modulus is
`E_calibrated = E_seed * u_seed / u_target`. Because all orthotropic moduli and
shear moduli scale together, the stiffness matrix scales linearly and the
displacement scales as `1/E`. Use `--verify-calibrated` when you want the script
to run a second solve at the calibrated value and print the residual difference.

The main `solver.py` path now does this automatically too. It first runs the
orthotropic 3-point seed solve at the AWC/reference modulus, computes
`calibrated_longitudinal_E(seed_E, seed_fea_uy, target_uy)`, then re-solves at
the calibrated orthotropic `Ez`. The printed target still uses the
AWC/reference `E`; the calibrated `Ez` is the solver-side stiffness adjustment.

## Tutorial Files

- `checkme.py` computes the beam-table theoretical deflection for a free-end
  point load.
- `solver.py` integrates the same total load over the free-end face,
  solves the model, and prints displacement and strain checks.
- `view.py` opens the PyVista displacement view.
- `tutorial_01_mesh_and_boundary_conditions.md` introduces the mesh and fixed
  boundary.
- `tutorial_02_how_k_connects_nodes.md` explains how element connectivity maps
  into the global stiffness matrix.

## Notes

The current model is intentionally exploratory. The theoretical check and finite
element result are kept side by side so stiffness, loading, boundary condition,
and Hermite rotation-scaling assumptions can be reviewed as the formulation
evolves.

## Current State

[Hermite beam theoretical check](https://cognitiveforge630.github.io/cognitiveforge.github.io/2026/05/10/hermite-beam-theoretical-check.html?v=2697cb6)

![Current HermiteStructural state](Screenshot%202026-05-11%20174755.png)
