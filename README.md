# HermiteStructural

Finite element tutorials and source code for experimenting with a Hermite
hexahedral structural model in Python.

The current assembly follows a Hermitian Hexa8-style formulation with six
degrees of freedom per node. Translation DOFs use `NH` shape functions, rotation
DOFs use `RH` shape functions, and stiffness is assembled from the standard
six-component small-strain matrix. Rotational DOFs are coupled into
displacement through the small-rotation relation `u = theta x r`.

The element uses two-point Gauss integration by default. Three-point integration
over-stiffens this Hermitian bending formulation on the cantilever benchmark,
while one-point integration is unstable.

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
