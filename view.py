"""Published Hermite cantilever viewer; Cauchy viewer is view_cauchy.py."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable,str(ROOT/'legacy_hermite'/'view.py'),*sys.argv[1:]],cwd=ROOT/'legacy_hermite'))
