import os
import platform
import shutil
import subprocess
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools import setup, Distribution

project_root = os.path.dirname(os.path.abspath(__file__))
platform_exec = "exe" if platform.system() == "Windows" else "bin"


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


def _compile_dir(src_dir: str, out_bin: str, extra_flags: list[str]) -> None:
    cmd = "gfortran ./*.f -O3 %s -o %s" % (
        " ".join("%s" % extra_flags[_] for _ in range(len(extra_flags))),
        out_bin,
    )
    print(src_dir)
    print(cmd)
    proc = subprocess.run(cmd, cwd=src_dir, shell=True, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"[pygrnwang] Fortran compile failed for {src_dir}\n"
            f"Error:\n{proc.stderr or proc.stdout}"
        )


class CustomBuildPy(_build_py):
    def run(self):
        _build_py.run(self)

        if not shutil.which("gfortran"):
            raise ValueError(
                r"Please install gfortran and ensure "
                r"that command \'gfortran\' can be directly called"
            )
        if not shutil.which("java"):
            raise ValueError(
                r"Please install java and ensure that "
                r"command \'java\' can be directly called"
            )

        exec_dir = os.path.join(project_root, "pygrnwang", "exec")
        os.makedirs(exec_dir, exist_ok=True)

        fortran_src_root = os.path.join(project_root, "fortran_src_codes")
        fortran_subdirs = {
            "edgrn2.0_src": f"edgrn2.{platform_exec}",
            "edcmp2.0_src": f"edcmp2.{platform_exec}",
            "qssp2020_src": f"qssp2020.{platform_exec}",
            "spgrn2020_src": f"spgrn2020.{platform_exec}",
            "qseis06_src": f"qseis06.{platform_exec}",
            "qseis2025_src": f"qseis2025.{platform_exec}",
        }

        for src_folder, bin_name in fortran_subdirs.items():
            fortran_src_dir = os.path.join(fortran_src_root, src_folder)
            output_binary = os.path.join(exec_dir, bin_name)

            extra = []
            env_fflags = os.environ.get("PYGRNWANG_FFLAGS", "")
            if env_fflags:
                extra += env_fflags.split()
            if src_folder == "qseis2025_src":
                extra += ["-ffixed-line-length-none"]

            print(f"[pygrnwang] Compiling {src_folder} -> {output_binary}")
            _compile_dir(fortran_src_dir, output_binary, extra)


setup(cmdclass={"build_py": CustomBuildPy}, distclass=BinaryDistribution)  # type:ignore
