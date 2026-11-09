# Linear Hermite: planned stiffness investigation

**Status: not accepted for cantilever bending.** With a 3×3×5-node mesh, the linear recovered-rotation element gives `−0.005385648 in` versus `−0.116149091 in` from the beam reference (approximately **95.36% smaller deflection**). This is severe apparent over-stiffness, not evidence of accuracy. The number is a single-mesh diagnostic and cannot by itself identify which mechanism causes the behavior.

## Next investigations, in a sensible order

1. **Patch and algebra audit.** Check rigid-body modes and every affine displacement/rotation state, identify nullspace and positive modes, confirm that the primitive B and `T` mapping have the intended signs and scale factors.
2. **Boundary conditions.** Compare full nonlocal face displacement constraints with rotations fixed versus rotations released (`constrain_rotations=False`); measure constraint rank, exact trace residual and rigid-motion suppression. Do **not** soften the clamp merely to match a theory target.
3. **Bending and shear modes.** Apply a manufactured constant-moment bending field with analytically prescribed boundary trace. Inspect recovered vertex gradients, symmetric/skew decomposition, energy split and transverse/shear strain locking.
4. **Mesh sequence.** Use identical geometry, prescribed material and loads for linear/quadratic/cubic grids. Track tip uy, energy, equilibrium and stress QoIs versus h rather than using a single difference.
5. **Independent controls.** Compare to high-order Cauchy or analytical benchmark under identical mathematical boundary conditions; document the model assumptions before attributing errors to elements.
6. **Only then adjust formulation.** Fix interpolation/recovery/constraint issues and add tests. Do not alter `E`, gravity or the load to make a graph look better.

Target acceptance should include stable patch tests, reliable equilibrated loads, a documented convergence study, justified QoI tolerances and eventually a rigorous a posteriori bound. The quadratic success does not automatically make the current linear variant acceptable.
