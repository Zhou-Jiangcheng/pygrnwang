"""Prepare the regional QSEIS example's physical 64 s moment-rate pulse.

QSEIS transforms a user wavelet at real frequency, although its Green functions
use f + i*fi. Its final inverse transform removes exp(2*pi*fi*t). Consequently,
we supply r(t)*exp(2*pi*fi*t), where r(t)=2/64*sin(pi*t/64)**2 on [0, 64] s.
The input's integral is about 0.9647: renormalizing that input to unit area would
undo the physical compensation. The effective moment-rate integral is one.

Type 0 returns rate kernels. Read velocity, strain_rate, or stress_rate and
integrate once in the calling example; the reader does not integrate type 0.
This helper is deliberately restricted to the regional tutorial's one group.
"""
import hashlib
import json
from pathlib import Path
import re

import numpy as np


SCHEME_ID = "normalized_sin_squared_64s_damping_compensated_v1"
EXPECTED_SETTINGS = {
    "event_depth_list": [10.0], "receiver_depth_list": [0.0],
    "grn_dist_range": [300.0, 900.0], "grn_delta_dist": 300.0,
    "N_dist": 3, "N_dist_group": 1, "N_each_group": 3,
    "sampling_interval": 4.0, "sampling_num": 1024, "time_window": 4092.0,
    "wavelet_type": 0, "wavelet_duration": 16, "anti_alias": 0.01,
    "time_reduction_velo": 0, "flat_earth_transform": True,
    "earth_model_layer_num": 24,
}


def _guarded_path(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("STF target escapes its intended directory: %s" % path)
    return path


def _paths(library):
    library = Path(library).expanduser().resolve(strict=True)
    if not library.is_dir():
        raise ValueError("Expected an existing library directory")
    output = library.parent
    paths = {
        "input": _guarded_path(library, "10.00/0.00/0_0/grn.inp"),
        "finished": _guarded_path(library, "10.00/0.00/0_0/.finished"),
        "library_info": _guarded_path(library, "green_lib_info.json"),
        "metadata": _guarded_path(library, "stf.json"),
        "samples": _guarded_path(output, "source_time_function.npz"),
        "description": _guarded_path(output, "source_time_function.json"),
    }
    info = json.loads(paths["library_info"].read_text(encoding="utf-8"))
    mismatches = [name for name, expected in EXPECTED_SETTINGS.items()
                  if info.get(name) != expected]
    if mismatches:
        raise ValueError("STF helper requires the regional type-0 library settings: "
                         + ", ".join(mismatches))
    return paths


def _read_input(path):
    data = path.read_bytes()
    encoding = "utf-8-sig" if data.startswith(b"\xef\xbb\xbf") else "utf-8"
    text = data.decode(encoding)
    newline = "\r\n" if "\r\n" in text else "\n"
    lines = text.splitlines(keepends=True)
    headings = [index for index, line in enumerate(lines)
                if line.lstrip("# \t").startswith("SOURCE TIME FUNCTION (WAVELET) PARAMETERS")]
    selectors = [index for index, line in enumerate(lines) if line.strip() == "16 0"]
    if len(headings) != 1 or len(selectors) != 1 or selectors[0] <= headings[0]:
        raise ValueError("Expected exactly one SOURCE TIME FUNCTION heading and one '16 0' selector")
    first_data = next((index for index in range(headings[0] + 1, len(lines))
                       if lines[index].strip() and not lines[index].lstrip().startswith("#")), None)
    if first_data != selectors[0]:
        raise ValueError("Unexpected data before the custom-wavelet selector")
    return data, lines, selectors[0], encoding, newline


def _integral(values, times):
    return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(times)))


def _source_definition():
    duration, count = 64.0, 1024
    fi = float(np.log(0.01) / (2 * np.pi * 4092.0))
    time = np.linspace(0.0, duration, count)
    target = 2.0 / duration * np.sin(np.pi * time / duration) ** 2
    target[[0, -1]] = 0.0
    supplied = target * np.exp(2 * np.pi * fi * time)
    fine_time = np.linspace(0.0, duration, 8193)
    effective = np.interp(fine_time, time, supplied) * np.exp(-2 * np.pi * fi * fine_time)
    fine_target = 2.0 / duration * np.sin(np.pi * fine_time / duration) ** 2
    fine_target[[0, -1]] = 0.0
    # Integrate the actual piecewise-linear input after removing damping, using
    # eight Gauss nodes per segment instead of assuming equality at the knots.
    nodes, weights = np.polynomial.legendre.leggauss(8)
    fraction = 0.5 * (nodes + 1)
    step = time[1] - time[0]
    quadrature_times = time[:-1, None] + step * fraction
    quadrature_input = supplied[:-1, None] + np.diff(supplied)[:, None] * fraction
    quadrature_rates = quadrature_input * np.exp(-2 * np.pi * fi * quadrature_times)
    effective_area = float(np.sum(quadrature_rates * weights) * step * 0.5)
    centroid = float(np.sum(quadrature_times * quadrature_rates * weights) * step * 0.5 / effective_area)
    # Reproduce qswavelet's analytic transform of linear segments and compare
    # it with the desired pulse at complex frequency, over the spherical band.
    frequency = np.arange(257, dtype=float) / 4096.0
    spectrum = np.zeros(frequency.size, dtype=complex)
    spectrum[0] = _integral(supplied, time)
    for index, value in enumerate(frequency[1:], 1):
        omega = 2 * np.pi * value
        alpha = np.exp(-1j * omega * step)
        beta = (alpha - 1) * 1j / omega
        gamma = alpha * 1j / omega - beta * 1j / omega / step
        phase = np.exp(-1j * omega * time[:-1])
        spectrum[index] = np.sum(phase * (supplied[:-1] * (beta - gamma) + supplied[1:] * gamma))
    z = (frequency + 1j * fi) * duration
    expected_spectrum = 1j * (np.exp(-2j * np.pi * z) - 1) / (2 * np.pi * z * (1 + z) * (1 - z))
    spectral_error = float(np.linalg.norm(spectrum - expected_spectrum) / np.linalg.norm(expected_spectrum))
    time_error = float(np.linalg.norm(effective - fine_target) / np.linalg.norm(fine_target))
    if (abs(effective_area - 1.0) >= 1e-8 or abs(centroid - 32.0) >= 1e-5
            or spectral_error >= 1e-5 or time_error >= 1e-5):
        raise ValueError("Custom source failed physical area, centroid, or pulse-shape validation")
    arrays = {
        "time_s": time, "target_rate": target, "input_rate": supplied,
        "effective_time_s": fine_time, "effective_rate_timegrid": effective,
        "target_rate_timegrid": fine_target,
    }
    metadata = {
        "scheme_id": SCHEME_ID, "duration_s": duration, "samples": count,
        "sample_spacing_s": float(step), "wavelet_type": 0, "wavelet_duration_samples": 16,
        "sampling_interval_s": 4.0, "native_time_window_s": 4092.0,
        "anti_alias": 0.01, "fi_hz": fi, "rate_unit": "1/s", "time_unit": "s",
        "physical_definition": "r(t) = 2/64 * sin(pi*t/64)^2 for 0 <= t <= 64 s; zero otherwise",
        "input_definition": "piecewise-linear samples of r(t) * exp(2*pi*fi*t); no input-area renormalization",
        "effective_definition": "piecewise-linear input(t) * exp(-2*pi*fi*t)",
        "target_rate_integral": _integral(target, time),
        "input_rate_integral": _integral(supplied, time),
        "effective_rate_integral": effective_area,
        "target_centroid_s": 32.0, "effective_centroid_s": centroid,
        "effective_rate_relative_l2": time_error,
        "spectral_relative_l2_0_to_0_0625_hz": spectral_error,
        "spectral_max_absolute_error_0_to_0_0625_hz": float(np.max(np.abs(spectrum - expected_spectrum))),
        "physical_validation_limits": {"area_absolute_error": 1e-8, "centroid_absolute_error_s": 1e-5,
                                       "rate_relative_l2": 1e-5, "spectral_relative_l2": 1e-5},
        "reader_instruction": "Read velo/strain_rate/stress_rate; integrate once with cumsum(rate)*4 s in the example",
    }
    return arrays, metadata


def prepare_qseis_stf(library, duration_s=64.0, samples=1024):
    """Insert the compensated rate into a fresh regional type-0 input.

    Call after preprocessing and before running QSEIS. This tutorial helper
    accepts only 64 s and 1024 custom-wavelet nodes. It refuses completed runs,
    existing custom blocks, and existing STF archives. The input's encoding and
    line endings are preserved. Return the metadata also saved in library/stf.json.
    """
    if duration_s != 64.0 or samples != 1024:
        raise ValueError("This helper is restricted to the 64 s, 1024-node regional example")
    paths = _paths(library)
    if paths["finished"].exists():
        raise ValueError("Refusing to change a library containing .finished; use a fresh output directory")
    for name in ("metadata", "samples", "description"):
        if paths[name].exists():
            raise ValueError("STF archive already exists; use validation for reuse: %s" % paths[name])
    original, lines, selector, encoding, newline = _read_input(paths["input"])
    following = next((line.strip() for line in lines[selector + 1:] if line.strip()), "")
    if re.fullmatch(r"#-+", following) is None:
        raise ValueError("Expected the untouched separator after '16 0'; a custom block may already exist")
    arrays, metadata = _source_definition()
    block = [str(samples) + newline]
    for start in range(0, samples, 8):
        block.append(" ".join(format(float(value), ".17g")
                              for value in arrays["input_rate"][start:start + 8]) + newline)
    prepared = "".join(lines[:selector + 1] + block + lines[selector + 1:]).encode(encoding)
    metadata.update(input_relative_path="10.00/0.00/0_0/grn.inp",
                    input_original_sha256=hashlib.sha256(original).hexdigest(),
                    input_prepared_sha256=hashlib.sha256(prepared).hexdigest(),
                    samples_relative_path="../source_time_function.npz",
                    description_relative_path="../source_time_function.json")
    # All targets and old structures were checked above. Recheck immediately
    # before changing the sole existing file; archives are created exclusively.
    if paths["finished"].exists() or paths["input"].read_bytes() != original:
        raise ValueError("The library input or completion status changed during preparation")
    with paths["samples"].open("xb") as stream:
        np.savez_compressed(stream, **arrays, physical_definition=metadata["physical_definition"],
                            scheme_id=SCHEME_ID, rate_unit="1/s", time_unit="s")
    metadata["samples_sha256"] = hashlib.sha256(paths["samples"].read_bytes()).hexdigest()
    description = json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    with paths["description"].open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(description)
    paths["input"].write_bytes(prepared)
    with paths["metadata"].open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(description)
    return validate_qseis_stf(library)


def validate_qseis_stf(library):
    """Validate an archived compensated source and unchanged native input.

    Completed runs are allowed here. Physical parameters, pulse arrays, file
    hashes and the inserted native block must agree before reusing the library.
    """
    paths = _paths(library)
    metadata_bytes = paths["metadata"].read_bytes()
    if paths["description"].read_bytes() != metadata_bytes:
        raise ValueError("The source description differs from library/stf.json")
    metadata = json.loads(metadata_bytes.decode("utf-8"))
    arrays, expected = _source_definition()
    for key, value in expected.items():
        actual = metadata.get(key)
        if isinstance(value, float):
            valid = isinstance(actual, (float, int)) and np.isclose(actual, value, rtol=1e-12, atol=1e-14)
        else:
            valid = actual == value
        if not valid:
            raise ValueError("Archived STF physical parameter changed: %s" % key)
    current, lines, selector, encoding, _ = _read_input(paths["input"])
    if hashlib.sha256(current).hexdigest() != metadata.get("input_prepared_sha256"):
        raise ValueError("The prepared native input has changed; the saved calculation cannot be reused")
    if hashlib.sha256(paths["samples"].read_bytes()).hexdigest() != metadata.get("samples_sha256"):
        raise ValueError("The archived source samples have changed")
    block_end = selector + 1 + 1 + 1024 // 8
    if lines[selector + 1].strip() != "1024":
        raise ValueError("Native custom-wavelet sample count changed")
    block_lines = lines[selector + 2:block_end]
    if len(block_lines) != 128 or any(len(line.split()) != 8 for line in block_lines):
        raise ValueError("Expected 128 native sample lines containing eight values each")
    actual_input = np.array([float(value) for line in block_lines for value in line.split()])
    if not np.allclose(actual_input, arrays["input_rate"], rtol=1e-12, atol=1e-16):
        raise ValueError("Native custom-wavelet samples do not match the physical definition")
    restored = "".join(lines[:selector + 1] + lines[block_end:]).encode(encoding)
    if hashlib.sha256(restored).hexdigest() != metadata.get("input_original_sha256"):
        raise ValueError("Original native input hash does not match the archived preparation")
    with np.load(paths["samples"], allow_pickle=False) as archive:
        for key, value in arrays.items():
            if (key not in archive.files or archive[key].shape != value.shape
                    or not np.allclose(archive[key], value, rtol=1e-12, atol=1e-16)):
                raise ValueError("Archived source array changed: %s" % key)
        for key, value in (("scheme_id", SCHEME_ID), ("physical_definition", expected["physical_definition"]),
                           ("rate_unit", "1/s"), ("time_unit", "s")):
            if key not in archive.files or archive[key].shape != () or archive[key].item() != value:
                raise ValueError("Archived source description changed: %s" % key)
    return metadata
