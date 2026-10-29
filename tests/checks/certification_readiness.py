"""Fail-closed status of numerical a posteriori certification obligations."""
import argparse
import json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--require-certified",action="store_true",help="exit nonzero unless all gates are complete")
    args=ap.parse_args()
    data=json.loads((Path(__file__).resolve().parents[2]/"docs"/"CERTIFICATION_STATUS.json").read_text())
    print(json.dumps(data,indent=2))
    complete=bool(data["certified_bound"]) and all(data["gates"].values())
    if args.require_certified and not complete:
        raise SystemExit("NOT CERTIFIED: at least one required gate remains unsatisfied")
    return complete

if __name__ == "__main__":
    main()
