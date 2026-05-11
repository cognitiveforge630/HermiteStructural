# Tutorial 02: How the Stiffness Matrix Connects Nodes

The global stiffness matrix `K` is not just a table of numbers. It is a record
of which degrees of freedom can influence each other through the finite element
mesh.

For this beam, each node has six degrees of freedom:

```text
ux, uy, uz, phix, phiy, phiz
```

Each hexahedral element has eight nodes. So one element works with:

```text
8 nodes * 6 dofs per node = 48 local dofs
```

That is why the element stiffness matrix `Ke` is `48 x 48`.

## The Local-to-Global Idea

For one element, the code gets the global degrees of freedom with:

```python
dofs = self.get_element_global_dofs(elem_idx).ravel()
```

Then every local stiffness entry is placed into the global matrix:

```python
K_row.append(dofs[a])
K_col.append(dofs[b])
K_data.append(Ke[a, b])
```

Read that as:

```text
local Ke[a, b] becomes global K[dofs[a], dofs[b]]
```

The element does not connect every node in the whole beam. It only connects the
eight nodes that belong to that element. Since neighboring elements share nodes,
those local blocks overlap and assemble into one global matrix.

## Run the Connectivity Tutorial

From this folder, run:

```powershell
py tutorial_02_stiffness_connectivity.py
```

By default, the script traces node `1211`, which belongs to several elements
near the start of the beam. That makes the shared-node coupling easier to see
than tracing a corner node.

To save the 3D connectivity view:

```powershell
py tutorial_02_stiffness_connectivity.py --screenshot outputs/k_connectivity.png
```

The script also writes:

```text
outputs/local_element_k_pattern.svg
outputs/global_k_pattern_small_mesh.svg
```

The local SVG shows one dense `48 x 48` element block. The global SVG shows how
many small element blocks become a sparse global matrix on a tiny demonstration
mesh.

## What the Colors Mean

In the 3D view:

- faint green shows the selected element,
- faint orange shows elements that touch the selected node,
- red sphere shows the selected node being traced,
- blue spheres show the other nodes on the selected element,
- orange spheres show other nodes that can couple to the selected node through
  shared element stiffness blocks,
- orange lines connect the selected node to those coupled nodes.

This is the picture to keep in mind before looking at the numeric values in
`K`: elements create small local connections, and assembly places those
connections into the correct global rows and columns.

The spheres are the important part of the view. The transparent element
surfaces are only context. A neighboring element matters to the selected node
because it shares one or more node numbers with it, not because every node in
the visible volume is pinned to every other node in a physical way.
