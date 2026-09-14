"""Verify the frequency grid of the five regional comparison libraries."""
import json
import math
from pathlib import Path
import struct

from common import REGIONAL_MAX_FREQUENCY_HZ, REGIONAL_SAMPLING_INTERVAL_S


_EXPECTED_NT = 1024
_EXPECTED_NF = 512
_EXPECTED_WINDOW = 4092.0
_EXPECTED_DF = 1.0 / 4096.0
_EXPECTED_ANTI_ALIAS = 0.01
_HEADER_NAMES = ("nt", "ntcut", "dt", "nf", "nfcut", "df", "ldegup")


def _metadata(library):
    return json.loads((library / "green_lib_info.json").read_text(encoding="utf-8"))


def _check_values(actual, expected, description):
    mismatches = []
    for key, wanted in expected.items():
        value = actual.get(key)
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            matches = False
        else:
            matches = math.isclose(value, wanted, rel_tol=1e-12, abs_tol=1e-12)
        if not matches:
            mismatches.append("%s=%r (expected %r)" % (key, value, wanted))
    if mismatches:
        raise ValueError("%s: %s. Rebuild in a fresh --output-dir without --reuse."
                         % (description, ", ".join(mismatches)))


def _report(header, fi_window, path, source):
    dt = header["dt"]
    df = header["df"]
    return {
        "requested_max_frequency_hz": REGIONAL_MAX_FREQUENCY_HZ,
        "nyquist_frequency_hz": 0.5 / dt,
        "highest_computed_frequency_hz": (header["nfcut"] - 1) * df,
        "nyquist_bin_zeroed": True,
        "nyquist_bin_zeroed_source": "Fortran inverse FFT sets bin nf + 1 to zero",
        "sampling_interval_s": dt,
        "fft_samples": header["nt"],
        "frequency_spacing_hz": df,
        "computed_frequency_samples": header["nfcut"],
        "spec_time_window_s": _EXPECTED_WINDOW,
        "fft_period_s": header["nt"] * dt,
        "anti_alias": _EXPECTED_ANTI_ALIAS,
        "fi_hz": math.log(_EXPECTED_ANTI_ALIAS) / (2.0 * math.pi * fi_window),
        "fi_window_s": fi_window,
        "fi_source": "Validated library anti_alias and the solver window formula",
        "ldegup": header["ldegup"],
        "native_header": header,
        "native_header_path": str(path),
        "validation_source": source,
    }


def verify_spherical_spectrum(library, backend):
    """Validate the native spectrum header for a 4 s, 0.125 Hz tutorial run.

    SPGRN2012/2020 and QSSP2020 write seven values in the first unformatted
    Fortran record. The frequency cutoff is checked before the metadata so
    an old 0.0625 Hz spectrum cannot pass merely by editing its JSON file.
    """
    backend = backend.lower()
    names = {"spgrn2012": "T_grn_d10.00", "spgrn2020": "T_grn_d10.00",
             "qssp2020": "U_Green_10.00km"}
    if backend not in names:
        raise ValueError("Unsupported spherical backend: %s" % backend)
    library = Path(library).expanduser().resolve()
    path = library / "GreenSpec" / "10.00" / "0.00" / names[backend]
    with path.open("rb") as native:
        marker = native.read(4)
        if marker == struct.pack("<I", 36):
            endian = "<"
        elif marker == struct.pack(">I", 36):
            endian = ">"
        else:
            raise ValueError("Expected a 36-byte Fortran spectrum header in %s" % path)
        payload = native.read(36)
        if len(payload) != 36 or native.read(4) != marker:
            raise ValueError("Truncated or invalid Fortran spectrum header in %s" % path)
    header = dict(zip(_HEADER_NAMES, struct.unpack(endian + "iidiidi", payload)))
    _check_values(header, {
        "nt": _EXPECTED_NT, "ntcut": _EXPECTED_NT,
        "dt": REGIONAL_SAMPLING_INTERVAL_S, "nf": _EXPECTED_NF,
        "nfcut": _EXPECTED_NF, "df": _EXPECTED_DF,
    }, "Native %s spectrum differs" % backend)
    if header["ldegup"] < 0:
        raise ValueError("Invalid negative harmonic cutoff in %s" % path)
    _check_values(_metadata(library), {
        "spec_time_window": _EXPECTED_WINDOW,
        "sampling_interval": REGIONAL_SAMPLING_INTERVAL_S,
        "max_frequency": REGIONAL_MAX_FREQUENCY_HZ,
        "anti_alias": _EXPECTED_ANTI_ALIAS,
    }, "%s library metadata differs" % backend)
    return _report(header, header["nt"] * header["dt"], path,
                   "Native Fortran spectrum header and library metadata")


def qseis_spectral_settings(library):
    """Verify QSEIS06/2025 input and all native output time labels.

    The text header names T_sec but contains no dt or sample count. Those
    values are verified from the output rows and native input, then nf/df
    are derived using qsgetinp.f. This helper is for the regional examples.
    """
    library = Path(library).expanduser().resolve()
    info = _metadata(library)
    _check_values(info, {
        "time_window": _EXPECTED_WINDOW,
        "sampling_interval": REGIONAL_SAMPLING_INTERVAL_S,
        "sampling_num": _EXPECTED_NT,
        "anti_alias": _EXPECTED_ANTI_ALIAS,
    }, "QSEIS library metadata differs")
    group = library / "10.00" / "0.00" / "0_0"
    input_path = group / "grn.inp"
    with input_path.open(encoding="utf-8") as native:
        records = []
        for line in native:
            line = line.strip()
            if line and not line.startswith("#"):
                records.append(line)
            if len(records) == 6:
                break
    if len(records) != 6 or len(records[5].split()) != 3:
        raise ValueError("Missing QSEIS time-sampling input record in %s" % input_path)
    start, window, nt = map(float, records[5].split())
    _check_values({"nt": nt, "time_window": window},
                  {"nt": _EXPECTED_NT, "time_window": _EXPECTED_WINDOW},
                  "Native QSEIS input differs")
    if not math.isfinite(start):
        raise ValueError("Nonfinite QSEIS start time in %s" % input_path)
    path = group / "ss.tz"
    with path.open(encoding="utf-8") as native:
        columns = native.readline().split()
        if not columns or columns[0] != "T_sec":
            raise ValueError("Missing QSEIS T_sec output header in %s" % path)
        times = [float(line.split()[0]) for line in native if line.strip()]
    _check_values({"nt": len(times)}, {"nt": _EXPECTED_NT},
                  "Native QSEIS output differs")
    dt = window / (int(nt) - 1)
    for index, timestamp in enumerate(times):
        if not math.isclose(timestamp, start + index * dt,
                            rel_tol=0.0, abs_tol=1e-5):
            raise ValueError("Native QSEIS time grid differs at row %d in %s"
                             % (index + 1, path))
    fft_nt = 1 << (int(nt) - 1).bit_length()
    header = {"nt": fft_nt, "ntcut": int(nt), "dt": dt, "nf": fft_nt // 2,
              "nfcut": fft_nt // 2, "df": 1.0 / (fft_nt * dt), "ldegup": None}
    _check_values(header, {
        "nt": _EXPECTED_NT, "dt": REGIONAL_SAMPLING_INTERVAL_S,
        "nf": _EXPECTED_NF, "nfcut": _EXPECTED_NF, "df": _EXPECTED_DF,
    }, "Derived QSEIS spectrum differs")
    report = _report(
        header, window, path,
        "Native input, all native T_sec output rows, and library metadata; "
        "nf, nfcut and df derived from qsgetinp.f (not stored in the text header)")
    report["native_input_path"] = str(input_path)
    report["native_output_columns"] = columns
    return report
