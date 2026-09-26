"""Regression tests for TauP subprocess queries and wheel resource installation.

Run with: python -m unittest discover -s tests -p test_pytaup_subprocess.py -v
"""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pygrnwang import pytaup as taup


class BackendTests(unittest.TestCase):
    def test_import_does_not_start_java_or_import_obspy_or_jpype(self):
        spec = importlib.util.spec_from_file_location("taup_import_test", taup.__file__)
        module = importlib.util.module_from_spec(spec)
        original_import = __import__

        def guarded_import(name, *args, **kwargs):
            if name.split(".")[0] in {"jpype", "obspy"}:
                raise AssertionError("Backend imported eagerly: " + name)
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=guarded_import), patch(
            "subprocess.run", side_effect=AssertionError("Java started on import")
        ):
            spec.loader.exec_module(module)

    def test_missing_javac_selects_obspy(self):
        with patch.object(taup.shutil, "which", side_effect=lambda n: None if n == "javac" else n):
            self.assertEqual(taup._detect_java_backend(), (False, None))

    def test_environment_jar_fallback(self):
        expected = str(Path(sys.exec_prefix) / (
            "Scripts" if sys.platform == "win32" else "bin"
        ) / "TauP.jar")
        with patch.object(taup.shutil, "which", return_value="java"), patch.object(
            taup.os.path, "isfile", side_effect=lambda path: path == expected
        ):
            self.assertEqual(taup._detect_java_backend(), (True, expected))

    def test_obspy_fallback(self):
        original_run = subprocess.run

        def guard_java(command, *args, **kwargs):
            if command[0] in {"java", "javac"}:
                raise AssertionError("Java was invoked")
            return original_run(command, *args, **kwargs)

        with patch.object(taup, "_USE_JAVA", False), patch.object(
            taup.subprocess, "run", side_effect=guard_java
        ):
            p, s = taup.cal_first_p_s(10, 1000)
            self.assertAlmostEqual(p, 129.8924, delta=0.02)
            self.assertAlmostEqual(s, 231.2147, delta=0.02)

    def test_installed_wheel_contains_both_jar_copies(self):
        package = Path(taup.__file__).parent
        if (package.parent / "pyproject.toml").is_file():
            self.skipTest("Run against an installed wheel to check its scripts directory")
        installed = Path(sysconfig.get_path("scripts")) / "TauP.jar"
        self.assertEqual(installed.read_bytes(), (package / "exec" / "TauP.jar").read_bytes())
        self.assertTrue((package / "java" / "PygrnwangTauPBatch.java").is_file())


@unittest.skipUnless(taup._USE_JAVA, "A JDK is needed for Java integration tests")
class JavaTests(unittest.TestCase):
    def test_full_arrivals_and_ray_parameter_units(self):
        from obspy.taup import TauPyModel

        expected = TauPyModel("ak135").get_travel_times(
            source_depth_in_km=10, distance_in_degree=1000 * taup._DEG_PER_KM,
            phase_list=["P", "S"],
        )
        actual = taup.taup_time_java(10, 1000, ["P", "S"])
        self.assertEqual(actual["phase"], [a.name for a in expected])
        self.assertEqual(actual["puristphase"], [a.purist_name for a in expected])
        np.testing.assert_allclose(actual["time"], [a.time for a in expected], atol=0.02)
        np.testing.assert_allclose(
            actual["rayparameter"], [a.ray_param for a in expected], atol=0.1,
        )

    def test_empty_arrivals(self):
        actual = taup.taup_time_java(10, 0, ["PKIKP"])
        self.assertTrue(all(value == [] for value in actual.values()))
        self.assertTrue(np.isnan(taup._first_arrival_java(10, 0, 0, "ak135", ["PKIKP"])))

    def test_source_receiver_reciprocity(self):
        np.testing.assert_allclose(taup.cal_first_p_s(5, 1000, 20), taup.cal_first_p_s(20, 1000, 5))

    def test_nd_path_with_spaces(self):
        import obspy.taup

        source = Path(obspy.taup.__file__).parent / "data" / "ak135f_no_mud.nd"
        with tempfile.TemporaryDirectory(prefix="taup model ") as folder:
            model = Path(folder) / "custom model.nd"
            shutil.copyfile(source, model)
            actual = taup.cal_first_p_s(10, 1000, model_name=str(model))
            expected = taup._cal_first_p_s_obspy(10, 1000, model_name=str(model))
            np.testing.assert_allclose(actual, expected, atol=0.03)

    def test_table_uses_one_java_process_and_preserves_order_and_dtype(self):
        distances = [1000, 500, 1000]
        expected = np.array([taup._cal_first_p_s_obspy(10, d) for d in distances])
        with tempfile.TemporaryDirectory() as folder, patch.object(
            taup.subprocess, "run", wraps=subprocess.run
        ) as run:
            taup.create_tpts_table(folder, 10, 0, distances)
            java_calls = [c for c in run.call_args_list if c.args[0][0] == "java"]
            self.assertEqual(len(java_calls), 1)
            directory = Path(folder) / "10.00" / "0.00"
            for index, name in enumerate(["tp_table.bin", "ts_table.bin"]):
                actual = np.fromfile(directory / name, dtype=np.float32)
                np.testing.assert_allclose(actual, expected[:, index], atol=0.02)
            before = run.call_count
            taup.create_tpts_table(folder, 10, 0, distances, check_finished=True)
            self.assertEqual(run.call_count, before)

    def test_empty_table_does_not_start_java(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(
            taup.subprocess, "run", side_effect=AssertionError("Empty query started Java")
        ):
            taup.create_tpts_table(folder, 10, 0, [])
            directory = Path(folder) / "10.00" / "0.00"
            self.assertEqual((directory / "tp_table.bin").stat().st_size, 0)
            self.assertEqual((directory / "ts_table.bin").stat().st_size, 0)

    def test_java_failure_reports_stderr(self):
        with patch.object(taup, "_compile_java_bridge") as compile_bridge, patch.object(
            taup.subprocess, "run", return_value=subprocess.CompletedProcess([], 1, "", "bad model")
        ):
            compile_bridge.return_value.name = "unused"
            with self.assertRaisesRegex(RuntimeError, "bad model"):
                taup.taup_time_java(10, 1000, ["P"])

    def test_incomplete_output_is_rejected(self):
        with patch.object(taup, "_compile_java_bridge") as compile_bridge, patch.object(
            taup.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")
        ):
            compile_bridge.return_value.name = "unused"
            with self.assertRaisesRegex(RuntimeError, "incomplete"):
                taup.taup_time_java(10, 1000, ["P"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
