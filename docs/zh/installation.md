# 安装

## 安装发布包

在已激活的 Python 环境中执行：

```bash
python -m pip install --upgrade pip
python -m pip install pygrnwang
python -c "import pygrnwang; print(pygrnwang.__version__)"
```

Python 最低版本为 **3.9**。发布配置覆盖 Linux x86-64、Windows x86-64 和 macOS Apple Silicon；具体版本能否直接安装 wheel，应以该次发布的文件为准。没有匹配 wheel 时，pip 可能转为源码编译，需要 `gfortran`。可用 `--only-binary=pygrnwang` 要求仅安装本库的 wheel。

```{important}
当前 Unix 系统上的 Python 批量计算函数查找环境 `bin/<solver>.bin`，
而标准 wheel 将原生程序保存在包内 `pygrnwang/exec`，因此可能报告找不到程序。
命令行包装器能找到包内副本；本教程的完整建库流程请使用下方源码可编辑安装。
这是已有打包路径限制，本轮文档工作不修改计算行为。
```

## 准备教程环境

教程脚本在 Git 仓库中。建议用源码对应的环境运行当前文档：

```bash
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
conda create -n pygrnwang -c conda-forge python=3.12 numpy scipy pandas obspy tqdm matplotlib gfortran
conda activate pygrnwang
python -m pip install -e .
```

源码构建要求 `setuptools>=77`，pip 会处理隔离构建依赖。安装过程中将编译七个 Fortran 程序；修改 Fortran 后需要重新构建。

Windows 下必须在已激活的 Conda 环境中运行，以便找到数值库和编译器的 DLL。非交互调用可写成：

```powershell
conda run -n pygrnwang python examples/qseis2025.py
```

不要通过未激活环境中的 `python.exe` 绝对路径启动计算。

## 可选的 Java 与 MPI

- **Java TauP：**需要 JDK 提供的 `java` 和 `javac` 同时位于 PATH。首次查询才编译桥接程序；导入 pygrnwang 不启动 JVM。通用走时接口在缺少 Java 后端时使用 ObsPy，显式 `taup_time_java` 调用则要求 Java 可用。
- **TauP.jar：**wheel 包内保留副本，并将另一份安装到环境的 `Scripts`（Windows）或 `bin`（Unix）；程序优先读取包内副本。
- **MPI：**仅多节点计算需要 `mpi4py` 和匹配的 MPI 运行时。默认串行教程无需 MPI。

运行 `pygrnwang` 命令只检查入口能否启动。验证真正的 Fortran 计算，请继续[快速开始](quickstart.md)。更多平台说明见[英文安装指南](../installation.md)。
