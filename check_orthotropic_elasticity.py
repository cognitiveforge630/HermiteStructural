"""Compatibility entry point for published HermiteStructural blog commands."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable,str(ROOT/'legacy_hermite'/'check_orthotropic_elasticity.py'),*sys.argv[1:]],cwd=ROOT/'legacy_hermite'))
