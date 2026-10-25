# Windows smoke commands

Run in workspace root, using the installed system Python:

```powershell
python -m unittest discover -s tests -v
python .\tutorial_001.py
python .\view.py
python .\view.py --screenshot .\outputs\tutorial001.png
python .\convergence.py
```

Report any failures and attach `outputs/convergence.csv` plus a PyVista screenshot. No PyPardiso dependency. The historical blog's `solver.py` and `U.npy` contract is **not** implemented in this separate educational workspace; do not overwrite the public repository yet.

For current Redwood Hermite primary demonstration use `docs/RUN_COMMANDS_IN_ORDER.md`; gravity is mandatory.
