# Tutorial 01: Beam Mesh and Boundary Conditions

Before solving a finite element model, verify the geometry, mesh, and boundary
conditions. This project starts with a rectangular beam:

```text
Lx = 4
Ly = 8
Lz = 120
```

The mesh uses:

```text
nx = 5
ny = 9
nz = 121
```

That means the model has 5 x 9 x 121 nodes and 4 x 8 x 120 hexahedral
elements. Each node has six degrees of freedom:

```text
ux, uy, uz, phix, phiy, phiz
```

The first boundary condition is a fully pinned face at `z = 0`. Since that face
contains `nx * ny = 45` nodes, the model fixes `45 * 6 = 270` degrees of
freedom before solving.

## Run the Mesh View

From this folder, run:

```powershell
py tutorial_01_mesh_view.py
```

To save a screenshot for the blog tutorial, run:

```powershell
py tutorial_01_mesh_view.py --screenshot outputs/beam_mesh_pinned_nodes.png
```

The mesh view draws:

- the hexahedral beam mesh as a wireframe,
- red cubes at the pinned nodes on `z = 0`.

This gives students a concrete first picture of what the solver is about to
turn into matrices.
