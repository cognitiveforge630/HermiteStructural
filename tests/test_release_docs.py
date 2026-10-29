"""Check that the public archive keeps one README and consistent labels."""
import unittest
from pathlib import Path
from core.materials.design_limits import REDWOOD_BASIS
from core.visualization.sampling import DEFAULT_FIELD


class ReleasePresentationTests(unittest.TestCase):
    def test_one_readme_and_current_baseline(self):
        root=Path(__file__).resolve().parents[1]
        readmes=sorted(root.rglob('README.md'))
        self.assertEqual(readmes,[root/'README.md'])
        content=(root/'README.md').read_text(encoding='utf-8')
        self.assertIn('1,100,000 psi',content)
        self.assertIn('95.36%',content)
        self.assertIn('2.16%',content)
        self.assertIn('12.23%',content)
        self.assertIn('19.82%',content)
        self.assertIn('not a certified',content.lower())
        self.assertEqual(REDWOOD_BASIS.E_longitudinal_psi,1_100_000.)
        self.assertEqual(DEFAULT_FIELD,'Vertical displacement magnitude (in)')
        self.assertTrue((root/'docs'/'INDEX.md').is_file())


if __name__=='__main__':
    unittest.main()
