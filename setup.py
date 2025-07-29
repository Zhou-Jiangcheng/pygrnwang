import platform
import shutil
import os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py

# Read dependencies from requirements.txt
with open("requirements.txt", encoding="utf-8") as f:
    # Filter out comments and empty lines
    requirements = [
        line.strip() for line in f if line.strip() and not line.startswith("#")
    ]

if platform.system() == "Windows":
    platform_exec = "exe"
else:
    platform_exec = "bin"


class CustomBuildPy(_build_py):
    def run(self):
        # Check the environment
        if not shutil.which("gfortran"):
            raise ValueError(
                "Please install gcc and ensure that command 'gfortran' can be directly called"
            )
        if not shutil.which("jar"):
            raise ValueError(
                "Please install java and ensure that command 'jar' can be directly called"
            )

        # Run the standard build process first.
        _build_py.run(self)

        exec_dir = os.path.join(os.getcwd(), "pycfs", "exec")
        fortran_subdirs = {
            "edgrn2.0_src": "edgrn2.%s" % platform_exec,
            "edcmp2.0_src": "edcmp2.%s" % platform_exec,
            "qssp2020_src": "qssp2020.%s" % platform_exec,
            "qseis_stress_src": "qseis_stress.%s" % platform_exec,
            "spgrn2020_src": "spgrn2020.%s" % platform_exec,
        }
        for src_folder, bin_name in fortran_subdirs.items():
            fortran_src_dir = os.path.join(os.getcwd(), "fortran_src_codes", src_folder)
            output_binary = os.path.join(exec_dir, bin_name)
            if src_folder == 'qseis_stress_src':
                compile_command = (f"gfortran {fortran_src_dir}/*.f -O3 "
                                   f"-ffixed-line-length-none -std=legacy "
                                   f"-fPIC -Wl,--no-relax -o {output_binary}")
            else:
                compile_command = f"gfortran {fortran_src_dir}/*.f -O3 -o {output_binary}"
            print(
                f"Compiling Fortran sources in {src_folder} with command: {compile_command}"
            )
            result = subprocess.run(compile_command, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"Fortran compilation failed for {src_folder}")


setup(
    name="pygrnwang",
    version="1.0.0",
    author="Zhou Jiangcheng",
    author_email="zhoujcpku@outlook.com",
    description="This Python package serves as the frontend for calculating and building "
    "Green's function libraries for synthetic seismograms.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Zhou-Jiangcheng/pygrnwang",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "pygrnwang": ["exec/*"],
    },
    # Pass the dependencies read from requirements.txt to install_requires
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "pygrnwang=pygrnwang.main:main",
        ],
    },
    cmdclass={"build_py": CustomBuildPy},  # type:ignore
)
