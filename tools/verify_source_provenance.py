"""Verify the unchanged Hermite core files against pinned certified-build hashes."""
from pathlib import Path
import hashlib, json
root=Path(__file__).resolve().parents[1]
expected=json.loads((root/'docs/HERMITE_SOURCE_SHA256.json').read_text())
for name,sha in expected.items():
    path=root/'legacy_hermite'/name
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=sha:
        raise SystemExit(f'FAIL {name}: core differs from certified-build source')
    print(f'PASS core provenance: {name}')
print('NOTE: identical historical Hermite source does not constitute a certified Hermite a posteriori bound.')
