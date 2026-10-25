# Commands in order — Windows PowerShell

Open `HermiteStructural-Tutorial001.code-workspace` in VS Code. Run commands from this repository root, without a virtual environment.

## 1. Historical orthotropic Redwood Hermite cantilever (main tutorial)

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python .\tools\verify_source_provenance.py
python .\solver.py
python .\view.py --hud
```

The reference includes both the integrated end-face force and the Redwood beam self-weight. The historical Hermite solver holds its longitudinal modulus fixed; it does not calibrate its stiffness to force agreement. The reference and FEM both include gravity. Agreement does not certify the element. Solver outputs are in `legacy_hermite/`.

For a blog image:

```powershell
python .\view.py --hud --hud-font-size 10 --title-font-size 10 --scalar-title-font-size 10 --scalar-label-font-size 10 --bounds-font-size 10 --screenshot outputs/hermite_redwood_blog.png
```

The Hermite viewer resolves screenshot paths relative to `legacy_hermite/`.

## 2. Orthotropic material matrix lesson

```powershell
cd .\orthotropic_companion
python -m unittest discover -s tests -v
python .\compare_efunda_vs_hermite_orthotropic.py
cd ..
```

You can also open `orthotropic_companion/HermiteStructural-Orthotropic-Companion.code-workspace` independently.

## 3. Optional Cauchy Hex8 reference

```powershell
python .\tutorial_001.py
python .\view_cauchy.py --hud
python .\convergence.py
```

Cauchy is optional; it is not the featured historical Hermite demonstration.
