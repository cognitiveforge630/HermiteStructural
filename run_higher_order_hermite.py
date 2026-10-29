"""Published command kept active; defaults to the quadratic Hermite example.

Use --formulation linear_rotation_hermite|quadratic_rotation_hermite|
cubic_rotation_hermite. For independent saved outputs prefer python -m tutorials.<order>.
"""
from pathlib import Path
from tutorials.cantilever import main

if __name__ == "__main__":
    main(default_out=Path("hermite_result.json"))
