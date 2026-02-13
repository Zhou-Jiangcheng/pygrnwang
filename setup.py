import os
import platform
import shutil
import subprocess
import sys
from setuptools import setup, Distribution
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop

# Try to import editable_wheel (for `pip install -e .`)
try:
    from setuptools.command.editable_wheel import editable_wheel as _editable_wheel
except ImportError:
    # If the setuptools version is too old and does not support PEP 660, set it to None
    _editable_wheel = None

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
    print(cmd)
    proc = subprocess.run(cmd, cwd=src_dir, shell=True, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"[pygrnwang] Fortran compile failed for {src_dir}\n"
            f"Error:\n{proc.stderr or proc.stdout}"
        )


def install_binaries(target_exec_dir):
    """
    Core logic: compile Fortran and copy all binaries files into conda's bin directory
    target_exec_dir: directory to store compilation outputs (build directory or source directory)
    """
    print(f"[pygrnwang] Starting custom installation logic...")
    
    if not shutil.which("gfortran"):
        raise ValueError(r"Please install gfortran.")
    
    # Ensure the output directory exists
    os.makedirs(target_exec_dir, exist_ok=True)

    # Get the bin directory of the current conda environment
    if platform.system() == "Windows":
        env_bin_dir = os.path.join(sys.exec_prefix, 'Scripts')
    else:
        env_bin_dir = os.path.join(sys.exec_prefix, 'bin')
        
    print(f"[pygrnwang] Target environment bin: {env_bin_dir}")

    # 1. Compile Fortran
    fortran_src_root = os.path.join(project_root, "fortran_src_codes")
    fortran_subdirs = {
        "edgrn2.0_src": f"edgrn2.{platform_exec}",
        "edcmp2.0_src": f"edcmp2.{platform_exec}",
        "qssp2020_src": f"qssp2020.{platform_exec}",
        "spgrn2012_src": f"spgrn2012.{platform_exec}",
        "spgrn2020_src": f"spgrn2020.{platform_exec}",
        "qseis06_src": f"qseis06.{platform_exec}",
        "qseis2025_src": f"qseis2025.{platform_exec}",
    }

    for src_folder, bin_name in fortran_subdirs.items():
        fortran_src_dir = os.path.join(fortran_src_root, src_folder)
        output_binary = os.path.join(target_exec_dir, bin_name)

        extra = []
        env_fflags = os.environ.get("PYGRNWANG_FFLAGS", "")
        if env_fflags:
            extra += env_fflags.split()
        if src_folder == "qseis2025_src":
            extra += ["-ffixed-line-length-none"]

        print(f"[pygrnwang] Compiling {src_folder} -> {output_binary}")
        _compile_dir(fortran_src_dir, output_binary, extra)

        # Copy the executable to conda/bin
        if os.path.exists(output_binary):
            dest_link = os.path.join(env_bin_dir, bin_name)
            print(f"[pygrnwang] Installing binary to {dest_link}")
            try:
                shutil.copy2(output_binary, dest_link)
                st = os.stat(dest_link)
                os.chmod(dest_link, st.st_mode | 0o111)
            except Exception as e:
                print(f"[Warning] Could not copy binary to {env_bin_dir}: {e}")

# --- Custom command classes ---

class CustomBuildPy(_build_py):
    """Controls `pip install .`"""
    def run(self):
        _build_py.run(self)
        # Standard install: compile into build/lib/pygrnwang/exec
        exec_dir = os.path.join(self.build_lib, "pygrnwang", "exec")
        install_binaries(exec_dir)


class CustomDevelop(_develop):
    """Controls `python setup.py develop`"""
    def run(self):
        _develop.run(self)
        # Development install: compile directly into the source directory pygrnwang/exec
        exec_dir = os.path.join(project_root, "pygrnwang", "exec")
        install_binaries(exec_dir)


# Collect cmdclass
cmd_classes = {
    "build_py": CustomBuildPy,
    "develop": CustomDevelop,
}

# Extra handling for `pip install -e .` (PEP 660)
if _editable_wheel:
    class CustomEditableWheel(_editable_wheel):
        """Controls `pip install -e .`"""
        def run(self):
            _editable_wheel.run(self)
            # Editable install is also treated as development mode; compile into the source directory
            exec_dir = os.path.join(project_root, "pygrnwang", "exec")
            install_binaries(exec_dir)
    
    cmd_classes["editable_wheel"] = CustomEditableWheel


setup(
    cmdclass=cmd_classes,
    distclass=BinaryDistribution
)
