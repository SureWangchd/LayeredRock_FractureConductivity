"""Self-tests of analytical coefficients/algebra, NOT COMSOL verification."""
from pathlib import Path
import json
import math
import unittest
import numpy as np
from periphery_modal_response import (
    Layer, input_admittance, modal_admittances, equivalent_tensor,
    equivalent_modal_admittances, hybrid_admittances, band_layers,
    original_transfer_check, cosine_norms, modal_update_weights,
    apply_modal_update, solve_supplied_modal_core,
)

class TestReference(unittest.TestCase):
    def test_dc_series(self):
        L=(Layer(.01,.1),Layer(.02,.01))
        self.assertAlmostEqual(input_admittance(L,0),1/(.01/.1+.02/.01),14)
    def test_uniform_analytic(self):
        for n in (1,2,8,32,256):
            k=n*math.pi/.12
            self.assertAlmostEqual(input_admittance([Layer(.03,.1)],k)/(.1*k/math.tanh(k*.03)),1,14)
    def test_split_identical_layer(self):
        one=[Layer(.03,.1)]; split=[Layer(.008,.1),Layer(.009,.1),Layer(.013,.1)]
        np.testing.assert_allclose(modal_admittances(one,.12,64),modal_admittances(split,.12,64),rtol=1e-13)
    def test_independent_transfer(self):
        for L in (band_layers(.003,.01,.1),band_layers(.005,.1,.01),
                  (Layer(.002,.07),Layer(.007,.3),Layer(.021,.011))):
            for n in (0,1,3,8):
                k=n*math.pi/.12
                self.assertAlmostEqual(input_admittance(L,k)/original_transfer_check(L,k),1,12)
    def test_small_k_limit(self):
        L=band_layers(.005,.1,.01)
        self.assertAlmostEqual(input_admittance(L,1e-7)/input_admittance(L,0),1,12)
    def test_large_k_stable(self):
        L=band_layers(.003,.01,.1)
        y=modal_admittances(L,.12,2048)
        self.assertTrue(np.isfinite(y).all() and (y>0).all())
        self.assertAlmostEqual(y[-1]/(.01*2048*math.pi/.12),1,12)
    def test_order_matters_only_nonzero(self):
        a=modal_admittances(band_layers(.003,.01,.1),.12,8)
        b=modal_admittances(band_layers(.005,.01,.1),.12,8)
        self.assertAlmostEqual(a[0],b[0],14)
        self.assertGreater(np.max(abs(a[1:]-b[1:])),1e-3)
    def test_band_monotonicity_in_this_family(self):
        for host,band,direction in ((.01,.1,1),(.1,.01,-1)):
            a=modal_admittances(band_layers(.001,host,band),.12,16)
            b=modal_admittances(band_layers(.007,host,band),.12,16)
            self.assertTrue(np.all(direction*(a[1:]-b[1:])>0))
    def test_equivalent_dc(self):
        for host,band in ((.01,.1),(.1,.01)):
            L=band_layers(.004,host,band)
            self.assertAlmostEqual(equivalent_modal_admittances(L,.12,8)[0],input_admittance(L,0),14)
    def test_equivalent_is_position_blind(self):
        a=band_layers(.001,.01,.1);b=band_layers(.007,.01,.1)
        np.testing.assert_allclose(equivalent_tensor(a),equivalent_tensor(b),rtol=1e-14)
        np.testing.assert_allclose(equivalent_modal_admittances(a,.12,32),equivalent_modal_admittances(b,.12,32),rtol=1e-14)
    def test_hybrid_reference_tail(self):
        y=modal_admittances(band_layers(.003,.01,.1),.12,8)
        r=equivalent_modal_admittances(band_layers(.003,.01,.1),.12,8)
        h=hybrid_admittances(y,r,2)
        np.testing.assert_array_equal(h[:3],y[:3]);np.testing.assert_array_equal(h[3:],r[3:])
        self.assertTrue(np.all(h>0))
    def test_lowrank_quadratic_energy(self):
        rng=np.random.default_rng(4);C=rng.normal(size=(9,4))*.01;v=rng.normal(size=9)
        dy=np.array([0.,.2,-.1,.4]);N=cosine_norms(.12,3)
        w=modal_update_weights(dy,.12,.02)
        left=v@apply_modal_update(C,w,v)
        right=.02*np.sum(N*dy*(C.T@v/N)**2)
        self.assertAlmostEqual(left,right,14)
    def test_affine_port_series(self):
        W=.12;t=1.;h=.06;s=.03;k=t*W*s/h
        core=k*np.array([[1.,-1.],[-1.,1.]])
        yl=np.array([2.]);yu=np.array([.4])
        r=solve_supplied_modal_core(core,yl,yu,W,t)
        expected=t*W/(1/yl[0]+h/s+1/yu[0])
        self.assertAlmostEqual(r['G_S'],expected,14)
        self.assertAlmostEqual(r['upper_inward_current_A'],expected,14)
        self.assertAlmostEqual(r['lower_outward_current_A'],expected,14)
    def test_two_field_identity(self):
        # Artificial SPD/PSD core, just to check signs and affine energy algebra.
        rng=np.random.default_rng(5);B=rng.normal(size=(8,6));c=np.array([1,0,0,1,0,0.])
        B-=np.outer(B@c,c)/(c@c);S=.04*B.T@B
        yl=np.array([2.,3.,4.]);yu=np.array([.4,1.,2.])
        zl=yl*np.array([1.,.95,.90]);zu=yu*np.array([1.,1.1,1.2])
        r1=solve_supplied_modal_core(S,yl,yu,.12)
        r2=solve_supplied_modal_core(S,zl,zu,.12)
        N=cosine_norms(.12,2);dw=np.r_[N*(zl-yl),N*(zu-yu)]
        identity=float(np.dot(dw*r1['trace_coefficients_V'],r2['trace_coefficients_V']))
        self.assertAlmostEqual(r2['G_S']-r1['G_S'],identity,13)
        self.assertAlmostEqual(r2['G_S'],r2['upper_inward_current_A'],13)
    def test_directional_derivative(self):
        rng=np.random.default_rng(6);B=rng.normal(size=(8,6));c=np.array([1,0,0,1,0,0.])
        B-=np.outer(B@c,c)/(c@c);S=.04*B.T@B
        yl=np.array([2.,3.,4.]);yu=np.array([.4,1.,2.]);d=np.array([0.,.1,-.05]);h=1e-5
        r=solve_supplied_modal_core(S,yl,yu,.12)
        p=solve_supplied_modal_core(S,yl,yu+h*d,.12)
        m=solve_supplied_modal_core(S,yl,yu-h*d,.12)
        exact=np.dot(cosine_norms(.12,2)*d,r['trace_coefficients_V'][3:]**2)
        self.assertAlmostEqual((p['G_S']-m['G_S'])/(2*h),exact,10)
    def test_bad_inputs(self):
        with self.assertRaises(ValueError):Layer(0,.1)
        with self.assertRaises(ValueError):input_admittance([],0)
        with self.assertRaises(ValueError):input_admittance([Layer(.1,.1)],-1)
        with self.assertRaises(ValueError):band_layers(.029,.1,.01)
        with self.assertRaises(ValueError):cosine_norms(.12,-1)

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(TestReference)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    record={'scope':'Analytical layered-exterior coefficients and algebra only; NOT COMSOL, not new rock experiments.',
            'tests_run':result.testsRun,'passed':result.testsRun-len(result.failures)-len(result.errors),
            'failures':len(result.failures),'errors':len(result.errors),
            'numpy_version':np.__version__}
    Path(__file__).with_name('reference_selftest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
