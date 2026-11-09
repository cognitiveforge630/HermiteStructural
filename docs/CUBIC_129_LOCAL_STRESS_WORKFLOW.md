# Reproduce the 5×5×129 cubic Hermite local-shear study

**Research comparison only, not a certified 3D stress bound.** Never mistake the
0.002908% *integrated section-shear force* discrepancy for pointwise shear
stress accuracy. The **3.732 psi** classical quantity is the approximate
rectangular-section **maximum** `|tau_yz|` at `z=0`, `y=0`, due to the 12.8 lbf
end force **plus** 66.816 lbf self-weight. The local distribution varies with
`y` and `z`; compare FE and theoretical stresses at **identical** coordinates.

## Why a new solve is required

`samples/cubic_5x5x129.json` and
`samples/cubic_5x5x129_source_sections.json` contain scalar measurements and
section integrals, **not the full FEM displacement vector**. There is no
recoverable solution checkpoint for that saved run. **Do not reconstruct or
invent** a displacement/stress field from the saved tip displacement alone.

## Commands (Windows PowerShell, project root)

First, verify the new saved-solution / direct-corner export workflow without
any expensive densification solve:

```powershell
python -m unittest tests.test_cubic_replay_fields -v
```

Then run the dense cubic recovered-rotation Hermite **Hex8 source-element**
model. This is NOT a 5×5×129 Hex64 mesh. The following produces the exact
same *midpoint* section definitions used for the recorded 5×5 axial series
and saves a reusable solution state before postprocessing:

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 129 --section-only --section-axial-positions 0 --save-solution outputs/cubic_5x5x129_solved_dofs.npz --out outputs/cubic_5x5x129.json
```

**Resource warning:** this large nonlocal cubic Hermite solve requires
substantial system RAM. It was **not rerun** in the supplied limited-memory
environment. Do not treat the historical JSON as proof the above solve has
just succeeded. The solved DOF state can be reused after one successful run. Since the
2026-10-29 publication audit, its checkpoint format is
`hermite_cubic_solved_dofs_v2`. It retains the complete **6×6 orthotropic
constitutive matrix** and refuses replay if the active material is changed,
even if longitudinal $E_z$ stays the same. A legacy v1 checkpoint must be
re-solved or carefully reviewed/migrated; never silently reinterpret it.

After the solve, recover the actual original FE element's constitutive stress
from `D @ B @ (T @ U)` at **eight 95%-inset corners per source cell**:

```powershell
python -m tools.replay_cubic_fields outputs/cubic_5x5x129_solved_dofs.npz --subdivisions 1 --inset 0.95 --out-prefix outputs/cubic_5x5x129
```

This writes:

- `outputs/cubic_5x5x129_one_sided_corner_fields.vtu` — all directly
  evaluated pointwise stress components, tension, compression, von Mises,
  displacement, and separate beam reference stresses. Repeated display
  vertices preserve **one-sided source-element** stresses without
  cross-element averaging. Inset stress-evaluation coordinates are stored as
  three point arrays; visual geometry corners differ slightly.
- `outputs/cubic_5x5x129_stress_comparison.json` — theoretical 3.732 psi
  clamp-center maximum, FE peak sample and its coordinates, **beam reference
  at the same sample positions**, signed differences, and a raw discrete
  interior (30–90 in) sample-L2 discrepancy. This is not a certified L2 norm.

View the existing full-field VTU without repeating the solve:

```powershell
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "Tau yz (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "Beam theory tau_yz (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "FE minus beam tau_yz (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "Sigma zz (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "Von Mises (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic_5x5x129_one_sided_corner_fields.vtu --field "Vertical displacement magnitude (in)"
```

**Interpretation:** the 3.732 psi theory peak is precisely at the original
clamped face/section centerline; 95%-inset sampling does **not** include this
exact location. Compare to the theoretical shear evaluated at each inset
stress location, rather than dividing an arbitrary peak by 3.732 psi. The
3D clamp region can have genuine local effects not captured by the 1D
reference. To show convergence, repeat the saved-solution and replay process
for at least 5×5×33, 5×5×65 and 5×5×129 with the same sampling convention,
then compare a fixed physical location or matching interior-region norms.

No element or subcell centroid FE stress is evaluated. These saved fields are
research diagnostics, **not final validated local stress, not a certified
bound, and not for structural design**.

## Verification performed on the submitted project

The new checkpoint/replay process was verified on a **3×3×5 / 16-source-element
SMOKE TEST**, not the requested 5×5×129. The smoke solve recorded a net
vertical load of `−79.616 lbf` (end load plus gravity), and its NPZ state was
successfully replayed into a one-sided corner-field VTU and JSON report. This
small result is only an interface/integrity check; its local stresses are
**not** used to make density-accuracy claims. The 5×5×17 solve was also
attempted in a container restricted to 4 GiB of RAM, but the operating
system terminated it for exceeding the available memory. Therefore no
5×5×17 or 5×5×129 new numerical stress result is claimed here. These
memory-limited attempts do not imply a solver failure on an adequately sized
workstation.
