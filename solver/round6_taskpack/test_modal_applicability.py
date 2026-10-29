"""Algebraic tests only: no COMSOL or empirical-physics validation."""
import unittest
import numpy as np
from modal_applicability import modal_interval, contrast_interval

class TestInterval(unittest.TestCase):
    def test_zero_change(self):
        v=np.ones((2,5));y=np.ones((2,5))
        x=modal_interval(.01,v,y,y)
        self.assertEqual(x.first_variation_S,0)
        self.assertEqual(x.remainder_bound_S,0)
        self.assertEqual(x.target_gain_lower_S,.01)
    def test_zero_drive_at_modes(self):
        v=np.zeros((2,5));y=np.ones((2,5));z=y.copy();z[:,1:]=3
        x=modal_interval(.01,v,y,z)
        self.assertEqual(x.first_variation_S,0)
    def test_forbid_changed_zero_mode(self):
        y=np.ones((2,5));z=y.copy();z[0,0]=2
        with self.assertRaises(ValueError):modal_interval(.01,y,y,z)
    def test_forbid_nonpositive(self):
        y=np.ones((2,5));z=y.copy();z[0,2]=0
        with self.assertRaises(ValueError):modal_interval(.01,y,y,z)
    def test_signed_contrast(self):
        self.assertEqual(contrast_interval(np.array([1,3]),np.array([2,5]),np.array([-1,1])),(1.,4.))
    def test_random_positive_energy_models(self):
        rng=np.random.default_rng(20261026)
        for _ in range(120):
            nfree=18; nmodes=8;w=.06
            S=rng.normal(size=(nmodes,nfree))
            T=rng.normal(size=(nfree,nfree)); core=T.T@T+np.eye(nfree)*.3
            yr=np.exp(rng.normal(size=nmodes)); yt=yr*np.exp(rng.normal(size=nmodes))
            A0=core+S.T@np.diag(w*yr)@S
            A1=core+S.T@np.diag(w*yt)@S
            b=rng.normal(size=nfree)
            u0=np.linalg.solve(A0,b);u1=np.linalg.solve(A1,b)
            # Same affine drive and constant; their stationary-energy difference.
            exact=float(b@(u0-u1))
            vv=np.zeros((2,5));a=np.ones((2,5));c=np.ones((2,5))
            vv[:,1:]=(S@u0).reshape(2,4);a[:,1:]=yr.reshape(2,4);c[:,1:]=yt.reshape(2,4)
            x=modal_interval(1.,vv,a,c)
            self.assertGreaterEqual(exact+1e-10,x.correction_lower_S)
            self.assertLessEqual(exact-1e-10,x.correction_upper_S)

if __name__=='__main__':unittest.main()
