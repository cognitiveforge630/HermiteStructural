"""Display linear/quadratic/cubic JSON results side by side."""
import json
from pathlib import Path


def main():
    print(f"{'Order':12} {'Grid':12} {'Tip uy [in]':>16} {'vs AWC [%]':>14} {'OK <1%?':>9}")
    print("-"*70)
    root=Path("outputs")
    for order in ("linear","quadratic","cubic"):
        p=root/f"{order}.json"
        if not p.is_file():
            print(f"{order:12} no result yet: python -m tutorials.{order}")
            continue
        r=json.loads(p.read_text(encoding="utf-8"))
        grid="x".join(map(str,r["node_mesh"]))
        if r["formulation"] != f"{order}_rotation_hermite":
            raise ValueError(f"Wrong formulation in {p}: {r['formulation']}")
        print(f"{order:12} {grid:12} {r['fem_tip_uy_in']:16.9f} {r['difference_vs_awc_pct']:14.6f} {str(r['under_one_percent']):>9}")
    print("Comparison with beam theory is not a rigorous a posteriori error bound.")

if __name__ == "__main__":
    main()
