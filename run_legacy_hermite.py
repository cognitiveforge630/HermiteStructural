"""Run the historical Hermite implementation in its own directory."""
import subprocess
import sys
from pathlib import Path

legacy = Path(__file__).resolve().parent / "legacy_hermite"
raise SystemExit(subprocess.call([sys.executable, str(legacy / "solver.py"), *sys.argv[1:]], cwd=legacy))
