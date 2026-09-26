"""Regression tests for the EDCMP homogeneous half-space Green's functions.

A 1 m x 1 m rectangle with strike 180 once gave stress up to 1e5 Pa at
hundreds of km when source and receiver depths were equal (Okada DC3D
singular-branch test on the fault-edge extension). Far-field stress of a
unit-potency source must stay below a few mu / r^3.

Run with: python -m unittest discover -s tests -p test_edcmp_halfspace.py -v
"""
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

from pygrnwang.create_edcmp import (
    call_edcmp2,
    convert_edcmp2,
    create_inp_edcmp2,
    fm_base_list,
)
from pygrnwang.focal_mechanism import plane2mt

EXE = Path(sys.exec_prefix) / (
    "Scripts/edcmp2.exe" if sys.platform == "win32" else "bin/edcmp2.bin"
)
LAM, MU = 25.21168e9, 31.12616e9


class MechanismTests(unittest.TestCase):
    def test_base_mechanisms_are_unit_tensors(self):
        expected = np.array(
            [
                [1, 0, 0, -1, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, -1],
                [0, 0, 0, 0, 1, 0],
            ]
        )
        mts = np.array([plane2mt(1, *fm) for fm in fm_base_list])
        np.testing.assert_allclose(mts, expected, atol=1e-12)

    def test_no_strike_180_rectangle_degeneracy(self):
        for fm in fm_base_list:
            self.assertNotAlmostEqual(abs(fm[0]), 180.0)

    def test_source_size_by_medium(self):
        with tempfile.TemporaryDirectory() as tmp:
            for layered, size in ((True, "1"), (False, "0")):
                create_inp_edcmp2(tmp, 15, 15, (0, 10), 1, 3, (0, 0, 1, 0),
                                  layered=layered, lam=LAM, mu=MU)
                inp = Path(tmp, "edcmp2", "15.00", "15.00", "3", "grn.inp")
                row = inp.read_text().splitlines()[96].split()
                self.assertEqual(row[5:7], [size, size])


@unittest.skipUnless(EXE.is_file(), "edcmp2 executable is not installed")
class HalfSpaceFarFieldTests(unittest.TestCase):
    def test_far_field_stress_decays(self):
        dist_km = np.arange(0, 1001, 1.0)
        far = dist_km >= 10
        cwd = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                for depth in (14.0, 15.0, 16.0):
                    for mt_ind in range(5):
                        with self.subTest(source_depth=depth, mt_ind=mt_ind):
                            create_inp_edcmp2(tmp, depth, 15.0, (0, 1000), 1,
                                              mt_ind, (0, 0, 1, 0),
                                              layered=False, lam=LAM, mu=MU)
                            call_edcmp2(depth, 15.0, mt_ind, tmp)
                            sub = os.path.join(tmp, "edcmp2", "%.2f" % depth,
                                               "15.00", "%d" % mt_ind, "")
                            stress = convert_edcmp2(sub, 2, remove=True)
                            os.chdir(cwd)
                            self.assertEqual(stress.shape, (dist_km.size, 6))
                            self.assertTrue(np.isfinite(stress).all())
                            r = dist_km[far] * 1e3
                            scaled = np.abs(stress[far]).max(axis=1) * r**3 / MU
                            self.assertLess(scaled.max(), 10.0)
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    unittest.main(verbosity=2)
