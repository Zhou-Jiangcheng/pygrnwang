"""Regression coverage for stable double-couple projection and EDCMP readers."""
import json
from pathlib import Path
import tempfile
import unittest
import warnings

import numpy as np

from pygrnwang.focal_mechanism import mt2plane, plane2mt, tensor2full_tensor_matrix
from pygrnwang.geo import create_rotate_z_mat
from pygrnwang.read_edcmp import seek_edcmp2, seek_edcmp2_bulk


def projected(mt, plane=0):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        return plane2mt(1, *mt2plane(mt)[plane])


class DoubleCoupleProjectionTests(unittest.TestCase):
    def test_roundoff_in_rotated_strike_slip(self):
        mt = plane2mt(1, 90, 90, 180)
        a = create_rotate_z_mat(np.deg2rad(30.0))
        matrix = a.T @ tensor2full_tensor_matrix(mt) @ a
        scalar = matrix[(0, 0, 0, 1, 1, 2), (0, 1, 2, 1, 2, 2)]
        c, s = np.cos(np.deg2rad(30.0)), np.sin(np.deg2rad(30.0))
        x, y, z, u, v, w = mt
        batch = np.array([c*c*x+2*c*s*y+s*s*u, -c*s*x+(c*c-s*s)*y+c*s*u,
                          c*z+s*v, s*s*x-2*c*s*y+c*c*u, -s*z+c*v, w])
        self.assertLess(np.max(np.abs(scalar-batch)), 1e-15)
        for tensor in (scalar, batch, np.nextafter(batch, np.inf), np.nextafter(batch, -np.inf)):
            for plane in (0, 1):
                with self.subTest(tensor=tensor, plane=plane):
                    np.testing.assert_allclose(projected(tensor, plane), tensor, rtol=0, atol=2e-14)

    def test_roundtrip_axis_aligned_and_near_horizontal_planes(self):
        for strike in (0, 30, 90, 180, 270, 359.9):
            for dip in (0, 1e-8, 45, 90-1e-8, 90):
                for rake in (-180, -90, -45, 0, 45, 90, 180):
                    mt = plane2mt(1, strike, dip, rake)
                    for plane in (0, 1):
                        with self.subTest(strike=strike, dip=dip, rake=rake, plane=plane):
                            np.testing.assert_allclose(projected(mt, plane), mt, rtol=0, atol=2e-12)

    def test_roundtrip_random_double_couples_preserves_unit_moment(self):
        rng = np.random.default_rng(20261029)
        for _ in range(200):
            strike, dip, rake = rng.uniform(0,360), rng.uniform(0,90), rng.uniform(-180,180)
            expected = plane2mt(1, strike, dip, rake)
            for scale in (1e-12, 1, 1e20):
                for plane in (0, 1):
                    np.testing.assert_allclose(projected(expected*scale, plane), expected, rtol=0, atol=3e-12)

    def test_non_double_couple_projection_retains_principal_axes(self):
        rng = np.random.default_rng(19)
        for _ in range(20):
            q, _ = np.linalg.qr(rng.standard_normal((3,3)))
            matrix = q @ np.diag([-3.0,-1.0,4.0]) @ q.T
            expected = np.outer(q[:,2],q[:,2])-np.outer(q[:,0],q[:,0])
            mt = matrix[(0,0,0,1,1,2),(0,1,2,1,2,2)]
            for plane in (0,1):
                np.testing.assert_allclose(tensor2full_tensor_matrix(projected(mt,plane)),expected,rtol=0,atol=2e-14)


class StaticReaderProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.library = Path(cls.temp.name)
        info = dict(grn_source_depth_range=[1,60], grn_source_delta_depth=59,
                    grn_dist_range=[1,500], grn_delta_dist=249.5, obs_depth_list=[0.5])
        (cls.library/'green_lib_info.json').write_text(json.dumps(info))
        (cls.library/'noQ.nd').write_text('0 6 3.5 2.7\n80 6 3.5 2.7\n')
        rng = np.random.default_rng(31)
        for output,width in [('disp',3),('strain',6),('stress',6),('tilt',2)]:
            values = rng.normal(size=(2,1,5,3,width)).astype(np.float32)
            values.tofile(cls.library/f'edcmp2_{output}.bin')
            for i,depth in enumerate([1,60]):
                for basis in range(5):
                    sub=cls.library/'edcmp2'/f'{depth:.2f}'/'0.50'/str(basis)
                    sub.mkdir(parents=True,exist_ok=True)
                    values[i,0,basis].tofile(sub/f'{output}.bin')

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_defaults_match_batch_and_direct_double_couple_tensor(self):
        angles=np.array([[90,90,180],[0,90,-90],[25,0,45],[128,38,-61]]*3,dtype=float)
        az=np.array([30,0,90,180,270,45,60,120,30,90,270,359.9])
        depth=np.array([1,60]*6,dtype=float)
        distance=np.array([1,250.5,500]*4)
        tensors=np.array([plane2mt(1,*fm) for fm in angles])
        forms=[angles,np.column_stack([np.full(12,1e18),angles]),tensors,
               np.column_stack([np.full(12,1e18),tensors])]
        for fm in forms:
            for output in ('disp','strain','stress','tilt'):
                for times_mu in (False,True):
                    kwargs=dict(path_green=str(self.library),model_name=str(self.library/'noQ.nd'),
                                output_type=output,times_mu=times_mu)
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore',UserWarning)
                        batch=seek_edcmp2_bulk(event_depth_km_arr=depth,receiver_depth_km_arr=np.full(12,.5),
                            az_deg_arr=az,dist_km_arr=distance,focal_mechanism_arr=fm,**kwargs)
                        direct=seek_edcmp2_bulk(event_depth_km_arr=depth,receiver_depth_km_arr=np.full(12,.5),
                            az_deg_arr=az,dist_km_arr=distance,focal_mechanism_arr=fm,
                            check_convert_pure_dp=False,**kwargs)
                        scalar=np.array([seek_edcmp2(event_depth_km=depth[i],receiver_depth_km=.5,
                            az_deg=az[i],dist_km=distance[i],focal_mechanism=fm[i],**kwargs) for i in range(12)])
                    with self.subTest(format=fm.shape[1],output=output,times_mu=times_mu):
                        scale=max(float(np.max(np.abs(scalar))),float(np.max(np.abs(direct))))
                        np.testing.assert_allclose(batch,scalar,rtol=5e-6,atol=5e-7*scale)
                        np.testing.assert_allclose(batch,direct,rtol=5e-6,atol=5e-7*scale)


if __name__ == '__main__':
    unittest.main(verbosity=2)
