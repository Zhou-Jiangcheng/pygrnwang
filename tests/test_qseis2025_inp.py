"""Regression tests for reading the QSEIS2025 output_observables flags.

A custom wavelet (wavelet_type=0) inserts its samples into grn.inp above the
Green's-function output section. The flags were once read from a fixed line,
so a finished custom-wavelet library was reported as unreadable and rerun.

Run with: python -m unittest discover -s tests -p test_qseis2025_inp.py -v
"""
import tempfile
import unittest

from pygrnwang.create_qseis2025 import (
    create_dir_qseis2025,
    create_inp_qseis2025,
    read_output_observables_qseis2025,
)

FLAGS = [1, 0, 1, 1, 0]


def write_inp(folder, wavelet_type):
    N_dist, N_dist_group = create_dir_qseis2025(
        folder, 10.0, 0.0, [300.0, 900.0], 300.0, N_each_group=3
    )
    return create_inp_qseis2025(
        folder, 10.0, 0.0, [300.0, 900.0], 300.0, N_dist, N_dist_group, 3,
        time_window=4092.0, sampling_interval=4.0, output_observables=FLAGS,
        wavelet_duration=16, wavelet_type=wavelet_type,
    )


class ReadOutputObservablesTests(unittest.TestCase):
    def test_generated_input(self):
        with tempfile.TemporaryDirectory() as folder:
            self.assertEqual(read_output_observables_qseis2025(write_inp(folder, 1)), FLAGS)

    def test_custom_wavelet_samples_above_the_flags(self):
        with tempfile.TemporaryDirectory() as folder:
            path_inp = write_inp(folder, 0)
            with open(path_inp) as fr:
                lines = fr.readlines()
            # 1024 samples, eight per line, after the "16 0" selector
            selector = lines.index("16 0\n")
            block = ["1024\n"] + [" ".join(["0.5"] * 8) + "\n"] * 128
            with open(path_inp, "w") as fw:
                fw.writelines(lines[:selector + 1] + block + lines[selector + 1:])
            self.assertEqual(read_output_observables_qseis2025(path_inp), FLAGS)

    def test_missing_flags_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path_inp = write_inp(folder, 1)
            with open(path_inp) as fr:
                lines = fr.readlines()
            with open(path_inp, "w") as fw:
                fw.writelines(line for line in lines
                              if line.strip() != " ".join(map(str, FLAGS)))
            with self.assertRaises(ValueError):
                read_output_observables_qseis2025(path_inp)


if __name__ == "__main__":
    unittest.main(verbosity=2)
