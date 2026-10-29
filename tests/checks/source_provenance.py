"""Verify shipped active modules against the pinned reorganization snapshot.

The pin is a source-integrity *snapshot*, not a mathematical certificate.
If changing intentional source, review then explicitly repin.
"""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[2]
expected = json.loads((root / "docs" / "ACTIVE_SOURCE_SHA256.json").read_text())
failures = []
for name, digest in expected.items():
    path = root / name
    if not path.is_file():
        failures.append(name)
        print(f"FAIL missing: {name}")
        continue
    normalized = path.read_bytes().replace(b"\r\n", b"\n")
    actual = hashlib.sha256(normalized).hexdigest()
    if actual != digest:
        failures.append(name)
        print(f"FAIL altered: {name}: {actual}")
    else:
        print(f"PASS source: {name}")
if failures:
    raise SystemExit(f"FAIL: {len(failures)} active sources differ from snapshot")
print(f"PASS: {len(expected)} active sources identical to distribution snapshot")
print("NOTE: Identical bytes do not establish a certified a posteriori error bound.")
