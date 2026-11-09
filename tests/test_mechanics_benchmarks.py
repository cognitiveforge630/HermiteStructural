"""Book-aligned analytical mechanics and actual cubic-element polynomial modes.

These are deterministic verification regressions, never certification gates.
"""
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from core.reference.mechanics_of_materials import (
    cantilever_fields, rectangular_stresses, euler_bernoulli_bending_energy,
)
from core.reference.awc import awc_reference_components
from core.solvers.beam import build_beam
from core.verification.stress_tensor import (
    stress_matrix, stress_components, strain_matrix, rotate_stress, traction,
    principal_stresses, elastic_energy_density, cube_force_and_moment,
)
from core.verification.patch_benchmarks import test_modes as evaluate_patch_modes
from core.verification.solved_audit import source_elastic_strain_energy, solved_mechanics_report
from tools.mechanics_benchmarks import build_report


class MechanicsReferences(unittest.TestCase):
    def test_point_and_uniform_load_combine_in_displacement_and_force(self):
        L=120.; E=1.1e6; I=4*8**3/12; P=-12.8; w=-.0174*32
        z=np.linspace(0,L,13)
        full=cantilever_fields(z,length=L,E=E,I=I,P=P,w=w)
        a=cantilever_fields(z,length=L,E=E,I=I,P=P)
        b=cantilever_fields(z,length=L,E=E,I=I,w=w)
        for name in full:
            np.testing.assert_allclose(full[name],a[name]+b[name],atol=1e-12,rtol=1e-12)
        self.assertAlmostEqual(float(full['uy_in'][-1]),awc_reference_components()['total_uy_in'])
        self.assertAlmostEqual(float(full['shear_y_lbf'][-1]),P)
        self.assertAlmostEqual(float(full['shear_y_lbf'][0]),P+w*L)
        self.assertAlmostEqual(float(full['moment_x_lbf_in'][-1]),0)
        self.assertEqual(float(full['uy_in'][0]),0)
        self.assertEqual(float(full['slope_y_per_z'][0]),0)
        for zz in (20.,57.,90.):
            h=1.e-3
            m=lambda x:cantilever_fields(x,length=L,E=E,I=I,P=P,w=w)['moment_x_lbf_in']
            v=float(cantilever_fields(zz,length=L,E=E,I=I,P=P,w=w)['shear_y_lbf'])
            self.assertAlmostEqual(float((m(zz+h)-m(zz-h))/(2*h)),v,places=6)

    def test_bending_tension_compression_and_parabolic_shear(self):
        f=rectangular_stresses(np.array((-4.,0.,4.)),40,width=4,depth=8,length=120,P=-12.8,w=-.5568)
        self.assertLess(f['sigma_zz_psi'][0],0)
        self.assertGreater(f['sigma_zz_psi'][2],0)
        self.assertAlmostEqual(f['sigma_zz_psi'][0],-f['sigma_zz_psi'][2])
        self.assertEqual(f['tau_yz_psi'][0],0.)
        self.assertEqual(f['tau_yz_psi'][2],0.)
        self.assertAlmostEqual(f['tau_yz_psi'][1],1.5*(-12.8-.5568*80)/32)
        with self.assertRaises(ValueError):
            rectangular_stresses(5,40,width=4,depth=8,length=120)

    def test_work_energy_by_independent_gauss_integral(self):
        L=120.; E=1.1e6; I=4*8**3/12; P=-12.8; w=-.5568
        gp,gw=leggauss(5)
        z=L*(gp+1)/2
        displacement=cantilever_fields(z,length=L,E=E,I=I,P=P,w=w)['uy_in']
        tip=float(cantilever_fields(L,length=L,E=E,I=I,P=P,w=w)['uy_in'])
        work=P*tip+L/2*float(np.dot(gw,w*displacement))
        energy=euler_bernoulli_bending_energy(length=L,E=E,I=I,P=P,w=w)
        self.assertAlmostEqual(work,2*energy,places=12)
        self.assertGreater(energy,0)
        with self.assertRaises(ValueError):
            cantilever_fields(121,length=L,E=E,I=I)


class StressCubeInvariants(unittest.TestCase):
    def test_cauchy_symmetry_and_rotation_of_traction(self):
        s=np.array((2.,5.,11.,-3.,7.,-4.))
        theta=.612
        R=np.array(((np.cos(theta),-np.sin(theta),0),
                    (np.sin(theta),np.cos(theta),0),(0,0,1)))
        n=np.array((1.,2.,-2.))/3
        tr=rotate_stress(s,R)
        np.testing.assert_allclose(stress_matrix(tr),R@stress_matrix(s)@R.T,atol=1e-13)
        np.testing.assert_allclose(traction(tr,R@n),R@traction(s,n),atol=1e-13)
        np.testing.assert_allclose(principal_stresses(tr),principal_stresses(s),atol=1e-13)
        np.testing.assert_allclose(stress_components(stress_matrix(s)),s)
        self.assertEqual(len(principal_stresses(s)),3)

    def test_beer_uniaxial_stress_on_forty_five_degree_plane(self):
        # Independent closed form: an initially uniaxial 100 psi stress at 45 degrees
        # gives 50 psi normal, 50 psi complementary shear on rotated faces.
        angle=np.pi/4
        R=np.array(((np.cos(angle),-np.sin(angle),0),
                    (np.sin(angle),np.cos(angle),0),(0,0,1)))
        transformed=rotate_stress([100,0,0,0,0,0],R)
        np.testing.assert_allclose(transformed,[50,50,0,50,0,0],atol=2e-14)
        np.testing.assert_allclose(principal_stresses([0,0,0,0,30,0]),[30,0,-30])
        f,m=cube_force_and_moment([10,20,30,3,5,7],edge=.15)
        np.testing.assert_allclose(f,[0,0,0],atol=1e-14)
        np.testing.assert_allclose(m,[0,0,0],atol=1e-14)

    def test_engineering_shear_energy_double_counting_is_avoided(self):
        s=np.array((14.,18.,22.,4.,-3.,7.))
        eps=np.array((.001,.002,.003,.004,.005,.006))
        self.assertAlmostEqual(elastic_energy_density(s,eps),
                               0.5*np.sum(stress_matrix(s)*strain_matrix(eps)),places=13)
        self.assertAlmostEqual(elastic_energy_density(np.array((0,0,0,0,6.,0.)),
                                                       np.array((0,0,0,0,.01,0))),.03)
        with self.assertRaises(ValueError):
            stress_matrix([1,2,3])
        with self.assertRaises(ValueError):
            traction(s,(0,0,2))
        with self.assertRaises(ValueError):
            rotate_stress(s,np.diag((1.,1.,-1.)))


class SourceElementBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')

    def test_physical_source_stress_modes_and_pure_bending(self):
        report=evaluate_patch_modes(self.beam)
        self.assertFalse(report['certified_error_bound'])
        for case,entry in report['modes'].items():
            with self.subTest(case=case):
                self.assertLess(entry['maximum_absolute_stress_difference_psi'],1e-8)
                self.assertTrue(entry['stress_field_exactly_reproduced_at_test_points'])

    def test_constitutive_energy_of_exact_uniaxial_patch(self):
        b=self.beam
        strain=.0001
        props=b.orthotropic_constants
        n_zx=props.nuxz*props.Ez/props.Ex
        n_zy=props.nuyz*props.Ez/props.Ey
        U=np.zeros((len(b.grid_points),6))
        U[:,0]=-n_zx*strain*b.grid_points[:,0]
        U[:,1]=-n_zy*strain*b.grid_points[:,1]
        U[:,2]=strain*b.grid_points[:,2]
        expected=.5*props.Ez*strain**2*b.Lx*b.Ly*b.Lz
        actual=source_elastic_strain_energy(b,U.ravel())
        self.assertAlmostEqual(actual/expected,1.,places=10)

    def test_signed_source_samples_are_not_replaced_by_beam_analytical_stress(self):
        b=self.beam
        zeros=np.zeros(len(b.grid_points)*6)
        report=solved_mechanics_report(b,zeros,zeros,end_force=-12.8,weight_density=.0174)
        self.assertFalse(report['certified_error_bound'])
        self.assertEqual(report['source_FE_strain_energy_lbf_in'],0)
        samples=report['section_samples'][1]['local_stress_samples']
        self.assertEqual(samples[0]['FE_stress_voigt_psi'],[0.]*6)
        self.assertNotEqual(samples[0]['beam_sigma_zz_psi'],0.)
        self.assertLess(samples[0]['beam_sigma_zz_psi'],0.)

    def test_cli_benchmark_report_and_status_are_fail_closed(self):
        report=build_report()
        self.assertFalse(report['certified_error_bound'])
        self.assertTrue(all(report['checks'].values()))
        self.assertEqual(report['status'],'PASS_RESEARCH_BENCHMARKS_NOT_CERTIFIED')

if __name__=='__main__': unittest.main()
