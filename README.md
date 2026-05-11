# HermiteStructural

Finite element tutorials and source code for experimenting with a Hermite-style
structural beam model in Python.

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

The viewer shows the solved beam, pinned nodes at `z = 0`, displacement
magnitude, and a small overlay comparing theoretical and actual free-end
displacement.

## Tutorial Files

- `checkme.py` computes the beam-table theoretical deflection for a free-end
  point load.
- `solver.py` distributes the same total point load across the free-end nodes,
  solves the model, and prints displacement and strain checks.
- `view.py` opens the PyVista displacement view.
- `tutorial_01_mesh_and_boundary_conditions.md` introduces the mesh and fixed
  boundary.
- `tutorial_02_how_k_connects_nodes.md` explains how element connectivity maps
  into the global stiffness matrix.

## Notes

The current model is intentionally exploratory. The theoretical check and finite
element result are kept side by side so stiffness, loading, and boundary
condition assumptions can be reviewed as the formulation evolves.
