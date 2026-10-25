"""Published tutorial compatibility entry point."""
from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable,str(ROOT/"legacy_hermite"/"tutorial_04_dofmap_element_stepper.py"),*sys.argv[1:]],cwd=ROOT/"legacy_hermite"))
