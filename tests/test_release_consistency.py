"""Distribution-level checks independent of a PyVista graphical session."""
import hashlib
import json
from pathlib import Path
import unittest
import numpy as np
from core.visualization.sampling import sample_volume, DEFAULT_FIELD
from core.solvers.beam import build_beam

ROOT = Path(__file__).resolve().parents[1]

class ReleaseConsistencyTests(unittest.TestCase):
    def test_source_sha_snapshot_is_current(self):
        hashes=json.loads((ROOT / 'docs' / 'ACTIVE_SOURCE_SHA256.json').read_text())
        self.assertGreaterEqual(len(hashes), 19)
        for path,digest in hashes.items():
            with self.subTest(file=path):
                payload=(ROOT / path).read_bytes().replace(b'\r\n',b'\n')
                self.assertEqual(hashlib.sha256(payload).hexdigest(),digest)

    def test_source_support_flags_are_original_nodes(self):
        beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        data=sample_volume(beam,np.zeros(len(beam.grid_points)*6),subdivisions=2)
        self.assertEqual(data.point_fields['Fixed support node'].sum(),9)
        self.assertEqual(data.point_fields['Loaded end node'].sum(),9)
        np.testing.assert_allclose(data.points[data.point_fields['Fixed support node']>0,2],0)
        np.testing.assert_allclose(data.points[data.point_fields['Loaded end node']>0,2],120)
        np.testing.assert_allclose(data.point_fields[DEFAULT_FIELD],np.abs(data.point_fields['Uy (in)']))
        np.testing.assert_allclose(data.cell_fields['Absolute longitudinal stress (psi)'],
                                   np.abs(data.cell_fields['Sigma zz (psi)']))

    def test_matched_face_displacement_continuity(self):
        """Random DOFs must not open gaps across neighboring source Hex8 faces.

        A local numerical test, not a general conforming-method proof.
        """
        rng=np.random.default_rng(173)
        for form in ('linear','quadratic','cubic'):
            with self.subTest(formulation=form):
                beam=build_beam(nx=3,ny=3,nz=5,
                                formulation=f'{form}_rotation_hermite')
                U=rng.normal(size=beam.get_all_points().shape[0]*6)
                # Match physical x/y/z face locations from neighboring
                # structured source cells, with non-vertex face samples.
                pairs=((0,1,(.3,-.4,1.),(.3,-.4,-1.)),
                       (0,8,(1.,-.2,.4),(-1.,-.2,.4)),
                       (0,4,(.2,1.,.4),(.2,-1.,.4)))
                for e0,e1,p0,p1 in pairs:
                    ix0,N0=beam.displacement_operator(e0,p0)
                    ix1,N1=beam.displacement_operator(e1,p1)
                    np.testing.assert_allclose(N0@U[ix0],N1@U[ix1],
                                               atol=1e-10,rtol=1e-10)

    def test_saved_vtu_metadata_and_cantilever_trend(self):
        try:
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
        except ImportError:
            self.skipTest('VTK reader unavailable; test the VTUs on the desktop')
        for name in ('linear_3x3x5','quadratic_3x3x5','cubic_3x3x5','cubic_3x3x9'):
            with self.subTest(sample=name):
                reader=vtk.vtkXMLUnstructuredGridReader()
                reader.SetFileName(str(ROOT / 'samples' / (name+'.vtu')))
                reader.Update()
                g=reader.GetOutput(); pts=vtk_to_numpy(g.GetPoints().GetData())
                pd=g.GetPointData();cd=g.GetCellData()
                u=vtk_to_numpy(pd.GetArray('Uy (in)'))
                mag=vtk_to_numpy(pd.GetArray('Vertical displacement magnitude (in)'))
                np.testing.assert_allclose(mag,np.abs(u),atol=1e-12)
                fixed=vtk_to_numpy(pd.GetArray('Fixed support node'))>0
                loaded=vtk_to_numpy(pd.GetArray('Loaded end node'))>0
                self.assertEqual(int(fixed.sum()),9)
                self.assertEqual(int(loaded.sum()),9)
                self.assertTrue(np.allclose(pts[fixed,2],0))
                self.assertTrue(np.allclose(pts[loaded,2],120))
                self.assertLess(np.max(mag[fixed]),1e-9)
                self.assertGreater(np.mean(mag[loaded]),.001)
                sig=vtk_to_numpy(cd.GetArray('Sigma zz (psi)'))
                sig_abs=vtk_to_numpy(cd.GetArray('Absolute longitudinal stress (psi)'))
                np.testing.assert_allclose(sig_abs,np.abs(sig),atol=1e-9)
                self.assertTrue(np.all(np.isfinite(sig)))

    def test_documents_have_one_readme_and_no_unsafe_claim(self):
        readmes=sorted(ROOT.rglob('README.md'))
        self.assertEqual(readmes,[ROOT/'README.md'])
        main=(ROOT/'README.md').read_text(encoding='utf8')
        self.assertIn('2.16%',main)
        self.assertIn('19.82%',main)
        self.assertIn('no',main.lower())
        self.assertIn('not a certified',main.lower())

if __name__ == '__main__': unittest.main()
