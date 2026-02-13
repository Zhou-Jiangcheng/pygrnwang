import os
import sys
import subprocess
import platform

def _get_binary_path(name):
    """
    Locate the binary executable inside the package installation directory.
    """
    # Determine extension based on current OS (matches setup.py logic)
    ext = "exe" if platform.system() == "Windows" else "bin"
    filename = f"{name}.{ext}"
    
    # Path relative to this file: ./exec/filename
    # This works for both editable install and wheel install
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_dir, "exec", filename)

def _run(name):
    """
    Generic runner function.
    """
    bin_path = _get_binary_path(name)
    if not os.path.exists(bin_path):
        print(f"Error: Binary not found at {bin_path}")
        print("Please ensure the package was installed correctly with binary extensions.")
        sys.exit(1)
    
    # Pass arguments from sys.argv to the binary
    # sys.argv[1:] contains command line arguments passed by the user
    args = sys.argv[1:]
    
    # On Windows, ensure the environment's Library/bin is in PATH for DLLs (MinGW/Conda)
    if platform.system() == "Windows":
        conda_lib_bin = os.path.join(sys.prefix, 'Library', 'bin')
        if os.path.exists(conda_lib_bin):
            os.environ["PATH"] = conda_lib_bin + os.pathsep + os.environ.get("PATH", "")

    try:
        # Ensure executable permission on POSIX systems (Linux/Mac)
        if platform.system() != "Windows":
            import stat
            st = os.stat(bin_path)
            os.chmod(bin_path, st.st_mode | stat.S_IEXEC)

        # Run the binary and wait for it to finish
        cmd = [bin_path] + args
        ret = subprocess.call(cmd)
        sys.exit(ret)
        
    except KeyboardInterrupt:
        sys.exit(1)
    except Exception as e:
        print(f"Execution failed: {e}")
        sys.exit(1)

def run_edgrn2(): _run("edgrn2")
def run_edcmp2(): _run("edcmp2")
def run_qssp2020(): _run("qssp2020")
def run_spgrn2012(): _run("spgrn2012")
def run_spgrn2020(): _run("spgrn2020")
def run_qseis06(): _run("qseis06")
def run_qseis2025(): _run("qseis2025")
