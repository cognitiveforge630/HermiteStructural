import json
import unittest
from pathlib import Path

class CertificationLabelTests(unittest.TestCase):
    def test_unproven_gates_are_not_reported_certified(self):
        p=Path(__file__).resolve().parents[1]/"docs"/"CERTIFICATION_STATUS.json"
        r=json.loads(p.read_text())
        self.assertFalse(r["certified_bound"])
        self.assertIn(False,r["gates"].values())
