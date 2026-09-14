# 中文入门

pygrnwang 为 Wang 系列 Fortran 程序提供 Python 前端，用于建立格林函数库、合成地震波形及计算静态形变。

本节提供中文安装和快速开始。后端细节、科学约定和 API 以英文维护，示例脚本与英文文档共用。

```{toctree}
:maxdepth: 1

installation
quickstart
```

## 从哪里开始

1. 按[安装说明](installation.md)准备环境。
2. 完成 [QSEIS2025 小算例](quickstart.md)，得到位移波形和图。
3. 阅读[科学约定](../conventions.md)，核对坐标、分量、单位、震源归一化和时间零点。
4. 按[后端选择表](../backends/index.md)选择静态、分层或球对称模型计算流程。

程序支持 Python 3.9 及以上版本。TauP 已通过 Java 子进程调用，不需要 JPype；缺少 JDK 时通用走时接口使用 ObsPy。
