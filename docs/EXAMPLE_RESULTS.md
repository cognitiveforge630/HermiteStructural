# Current Redwood reference-grade example results (2026-10-26)

The **active** model uses `E_L=1,100,000 psi` from Redwood Select Structural
dimension-lumber unadjusted published design values. The loads are `P=-12.8 lbf`
and `γ=0.0174 lbf/in³` integrated over the volume; total `Fy=-79.616 lbf`.
The table below is the original recovered-rotation Hermite family on Hex8
connectivity; these rows are not the newer conventional Hex27/Hex64 elements.

| Method | Nodes | Tip uy (in) | Euler-Bernoulli beam uy (in) | Difference |
|---|---|---:|---:|---:|
| linear | 3×3×5 | −0.005385648 | −0.116149091 | 95.363160% |
| quadratic | 3×3×5 | −0.094852324 | −0.116149091 | 18.335716% |
| cubic | 3×3×5 | −0.114796439 | −0.116149091 | 1.164582% |
| cubic | 3×3×9 | −0.115608854 | −0.116149091 | 0.465124% |

Cubic raw **σzz** interior relative L² disagreement with beam theory:
3×3×5 subdivision 2 ≈ 9.12%; 3×3×9 subdivision 3 ≈ 2.16%.
These are *sampled stress discrepancies*, not rigorous 3D FEM error bounds.
Shear resultant comparison remains worse (12.23% on 3×3×9 subdiv 2).
All preliminary grade/demand checks passed on those cubic samples; see their
`_stress_audit.json` and `_sections.csv` for actual units and field samples.

**Historical, incompatible elastic modulus:** a quadratic 4×5×17 run at
`E_L=2,000,000 psi` agreed with its original theoretical reference by ≈0.483%.
That original JSON is saved under `docs/history/`; the datum is *not* a
new grade-aligned 4×5×17 result, and no large cubic run is claimed.

Neither a displacement agreement nor visual stress smoothness establishes a
correct, mesh-converged 3D stress distribution. For physical limitations,
source references and output interpretation read
[Cubic stress and Redwood](CUBIC_STRESS_AND_REDWOOD.md).
For the separate conventional-element results, see the
[Hex27/Hex64 higher-node study](HIGHER_NODE_HEX27_HEX64_STUDY_20261027.md).
