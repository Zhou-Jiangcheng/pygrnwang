"""Apply a physical 64 s moment-rate pulse to a full SPGRN2012 impulse record.

The native source_duration=0 branch has unit spectrum. We restore its damping,
transform the complete 1024 samples, multiply by the analytic source spectrum
at f+i*fi, invert, and remove damping. This is forward convolution, with no
spectral division, fitted amplitude, time shift, or sampled-pulse approximation.
The local time axis starts at each record's first sample; its original native
origin-time offset stays unchanged. Crop and integrate only after this step.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import FortranEOFError, FortranFile

from common import (REGIONAL_MAX_FREQUENCY_HZ, REGIONAL_SAMPLING_INTERVAL_S,
                    REGIONAL_SOURCE_DURATION_S, REGIONAL_STF)


SCHEME_ID = "normalized_sin_squared_64s_complex_frequency_forward_v1"


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inspect_spgrn2012_native(library):
    """Check the actual spectrum headers and all three complete velocity blocks."""
    library = Path(library).resolve()
    info = json.loads((library / "green_lib_info.json").read_text(encoding="utf-8"))
    expected = {"source_duration": 0.0, "spec_time_window": 4092.0, "time_window": 4092.0,
                "sampling_interval": REGIONAL_SAMPLING_INTERVAL_S, "samples_num": 1024,
                "max_frequency": REGIONAL_MAX_FREQUENCY_HZ, "max_slowness": 0.0,
                "anti_alias": 0.01, "dist_list": [300.0, 600.0, 900.0]}
    if any(info.get(key) != value for key, value in expected.items()):
        raise ValueError("SPGRN2012 requires the complete 1024-sample, Nyquist-band impulse library")
    dtype = np.dtype([("nt", "<i4"), ("ntcut", "<i4"), ("dt", "<f8"),
                      ("nf", "<i4"), ("nfcut", "<i4"), ("df", "<f8"), ("ldegup", "<i4")])
    headers = []
    for component in "RTP":
        path = library / "GreenSpec" / "10.00" / "0.00" / (component + "_grn_d10.00")
        with FortranFile(path, "r") as stream:
            record = stream.read_record(dtype)
        if record.size != 1:
            raise ValueError("Unexpected native spectrum header size")
        headers.append({name: record[0][name].item() for name in dtype.names})
    header = headers[0]
    expected_header = {"nt": 1024, "ntcut": 1024, "dt": 4.0, "nf": 512,
                       "nfcut": 512, "df": 1.0 / 4096.0}
    if (any(item != header for item in headers[1:])
            or any(header[key] != value for key, value in expected_header.items())):
        raise ValueError("Native spectral sampling does not match the requested Nyquist band")
    native_path = library / "GreenFunc" / "10.00" / "0.00" / "grn_d10.00"
    input_text = (native_path.parent / "grn.inp").read_text(encoding="utf-8")
    marker = "SPACE-TIME DOMAIN GREEN'S FUNCTIONS"
    if input_text.count(marker) != 1:
        raise ValueError("Expected one native space-time input section")
    records = [line.strip() for line in input_text.split(marker, 1)[1].splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    if (len(records) < 7 or list(map(float, records[2].split())) != [4092.0, 4.0]
            or list(map(float, records[3].split())) != [-40.0, 10.0]
            or float(records[4]) != 0.0 or float(records[5].split()[0]) != 0.0
            or list(map(float, records[6].split())) != [300.0, 900.0, 300.0, 300.0]):
        raise ValueError("Native input must request the unfiltered, complete impulse velocity")
    marker = "TIME (FREQUENCY) SAMPLING"
    if input_text.count(marker) != 1:
        raise ValueError("Expected one native time-frequency input section")
    records = [line.strip() for line in input_text.split(marker, 1)[1].splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    if (len(records) < 4 or list(map(float, records[0].split())) != [4092.0, 4.0]
            or float(records[1]) != REGIONAL_MAX_FREQUENCY_HZ
            or float(records[2]) != 0.0 or float(records[3]) != 0.01):
        raise ValueError("Native input must request the full-wavefield Nyquist-band spectrum")
    starts = []
    with FortranFile(native_path, "r") as stream:
        for _ in range(3):
            start = stream.read_reals(np.float32)
            if start.size != 1 or not np.isfinite(start).all():
                raise ValueError("Invalid native time-origin record")
            starts.append(float(start[0]))
            for _ in range(10):
                values = stream.read_reals(np.float32)
                if values.size != header["nt"] or not np.isfinite(values).all():
                    raise ValueError("Native impulse velocity must contain 1024 finite samples")
        try:
            stream.read_record(np.float32)
        except FortranEOFError:
            pass
        else:
            raise ValueError("Unexpected extra native distance records")
    if starts != [-10.0, 20.0, 50.0]:
        raise ValueError("Unexpected native starts for t0=-40 s and v0=10 km/s")
    return {
        "fft_header": header, "fft_period_s": header["nt"] * header["dt"],
        "fi_hz": float(np.log(info["anti_alias"]) * header["df"] / (2 * np.pi)),
        "requested_max_frequency_hz": REGIONAL_MAX_FREQUENCY_HZ,
        "last_retained_frequency_hz": (header["nfcut"] - 1) * header["df"],
        "nyquist_bin_zero": True, "trace_start_times_s": starts,
        "input_sha256": file_sha256(native_path.parent / "grn.inp"),
        "native_velocity_sha256": file_sha256(native_path),
    }


def sin_squared_spectrum(frequency_hz, fi_hz):
    """Return the unit-area target pulse's exact transform at complex frequency."""
    z = (np.asarray(frequency_hz) + 1j * fi_hz) * REGIONAL_SOURCE_DURATION_S
    # Here fi is strictly negative, so z never equals the removable poles 0, +/-1.
    if fi_hz >= 0:
        raise ValueError("This tutorial transform requires the native negative damping frequency")
    return 1j * np.expm1(-2j * np.pi * z) / (2 * np.pi * z * (1 + z) * (1 - z))


def apply_spgrn2012_stf(velocity, native):
    """Forward-convolve the complete physical impulse velocity, retaining its time axis."""
    values = np.asarray(velocity, dtype=float)
    if values.shape != (3, 3, 1024) or not np.isfinite(values).all():
        raise ValueError("Expected three distances by ENU by 1024 full-period velocity samples")
    dt = native["fft_header"]["dt"]
    frequency = np.fft.rfftfreq(1024, dt)
    transfer = sin_squared_spectrum(frequency, native["fi_hz"])
    damping = np.exp(2 * np.pi * native["fi_hz"] * np.arange(1024) * dt)
    spectrum = np.fft.rfft(values * damping, axis=-1)
    # The native inverse FFT explicitly sets Nyquist to zero. Remove only any
    # roundoff introduced there by storing the impulse velocity as float32.
    spectrum[..., -1] = 0.0
    matched = np.fft.irfft(spectrum * transfer, n=1024, axis=-1) / damping
    if not np.isfinite(matched).all():
        raise ValueError("Non-finite velocity after forward source convolution")
    # Independent quadrature checks the analytic transform over the whole band.
    nodes, weights = np.polynomial.legendre.leggauss(64)
    times = (nodes + 1) * REGIONAL_SOURCE_DURATION_S / 2
    rate = 2 / REGIONAL_SOURCE_DURATION_S * np.sin(np.pi * times / REGIONAL_SOURCE_DURATION_S) ** 2
    quadrature = (np.exp(-2j * np.pi * (frequency[:, None] + 1j * native["fi_hz"]) * times)
                  @ (weights * rate)) * REGIONAL_SOURCE_DURATION_S / 2
    error = float(np.linalg.norm(quadrature - transfer) / np.linalg.norm(transfer))
    area = float(np.sum(weights * rate) * REGIONAL_SOURCE_DURATION_S / 2)
    if error >= 1e-12 or abs(area - 1.0) >= 1e-12:
        raise ValueError("Analytic source transform failed independent full-band quadrature validation")
    metadata = {
        "scheme_id": SCHEME_ID, "physical_source_time_function": dict(REGIONAL_STF),
        "native_source_duration_s": 0.0, "effective_source_duration_s": REGIONAL_SOURCE_DURATION_S,
        "method": "full-period damped-domain analytic forward convolution of native impulse velocity",
        "spectral_division": False, "fitted_amplitude": False, "fitted_time_shift": False,
        "native": native, "physical_rate_integral": area,
        "analytic_transform_relative_l2_vs_quadrature": error,
        "damped_dc_coefficient": float(transfer[0].real),
        "integration": "cumsum of matched velocity times 4 s, then retain the first 256 samples",
    }
    return matched, metadata


def save_spgrn2012_stf(output, native, raw_velocity, matched_velocity, metadata):
    """Save full velocity records and source evidence; calculations remain untouched."""
    output = Path(output).resolve()
    time = (np.asarray(native["trace_start_times_s"])[:, None]
            + np.arange(1024)[None, :] * REGIONAL_SAMPLING_INTERVAL_S)
    paths = {}
    for filename, values in (("velocity-impulse.npz", raw_velocity),
                             ("velocity-matched.npz", matched_velocity)):
        path = (output / filename).resolve()
        if not path.is_relative_to(output):
            raise ValueError("Velocity archive target escaped the example output directory")
        np.savez_compressed(path, values=values, time_s=time, distance_km=[300., 600., 900.],
                            components=["E", "N", "U"], unit="m/s")
        paths[filename] = file_sha256(path)
    source_path = (output / "source_time_function.npz").resolve()
    description_path = (output / "source_time_function.json").resolve()
    if not source_path.is_relative_to(output) or not description_path.is_relative_to(output):
        raise ValueError("Source archive target escaped the example output directory")
    source_time = np.linspace(0, REGIONAL_SOURCE_DURATION_S, 1024)
    rate = 2 / REGIONAL_SOURCE_DURATION_S * np.sin(np.pi * source_time / REGIONAL_SOURCE_DURATION_S) ** 2
    rate[[0, -1]] = 0.0
    frequency = np.fft.rfftfreq(1024, REGIONAL_SAMPLING_INTERVAL_S)
    np.savez_compressed(source_path, time_s=source_time, target_rate=rate,
                        frequency_hz=frequency, transfer_function=sin_squared_spectrum(frequency, native["fi_hz"]),
                        scheme_id=SCHEME_ID, rate_unit="1/s")
    paths["source_time_function.npz"] = file_sha256(source_path)
    paths["disp.npz"] = file_sha256(output / "disp.npz")
    metadata["archive_sha256"] = paths
    description_path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return metadata


def validate_spgrn2012_stf(library):
    """Reject old native inputs, frequency limits, pulse definitions, or altered archives."""
    library = Path(library).resolve()
    metadata = json.loads((library.parent / "source_time_function.json").read_text(encoding="utf-8"))
    native = inspect_spgrn2012_native(library)
    if (metadata.get("scheme_id") != SCHEME_ID
            or metadata.get("physical_source_time_function") != REGIONAL_STF
            or metadata.get("native_source_duration_s") != 0.0
            or metadata.get("effective_source_duration_s") != REGIONAL_SOURCE_DURATION_S
            or metadata.get("native") != native):
        raise ValueError("SPGRN2012 source evidence no longer matches the impulse library and target pulse")
    for filename in ("velocity-impulse.npz", "velocity-matched.npz", "source_time_function.npz", "disp.npz"):
        if file_sha256(library.parent / filename) != metadata.get("archive_sha256", {}).get(filename):
            raise ValueError("Changed SPGRN2012 source/output archive: %s" % filename)
    with np.load(library.parent / "source_time_function.npz", allow_pickle=False) as archive:
        source_time = np.linspace(0, REGIONAL_SOURCE_DURATION_S, 1024)
        rate = 2 / REGIONAL_SOURCE_DURATION_S * np.sin(np.pi * source_time / REGIONAL_SOURCE_DURATION_S) ** 2
        rate[[0, -1]] = 0.0
        if (not np.array_equal(archive["time_s"], source_time)
                or not np.allclose(archive["target_rate"], rate, rtol=1e-14, atol=1e-16)
                or archive["scheme_id"].item() != SCHEME_ID
                or archive["rate_unit"].item() != "1/s"):
            raise ValueError("Archived physical source pulse changed")
        frequency = np.fft.rfftfreq(1024, REGIONAL_SAMPLING_INTERVAL_S)
        if (not np.array_equal(archive["frequency_hz"], frequency)
                or not np.allclose(archive["transfer_function"], sin_squared_spectrum(frequency, native["fi_hz"]),
                                   rtol=1e-12, atol=1e-14)):
            raise ValueError("Archived forward source transform changed")
    return metadata
