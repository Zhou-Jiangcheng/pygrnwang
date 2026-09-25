# QSEIS2025：从模型到位移波形

此教程使用仓库提供的 AK135 弹性参数、一个震源深度、一个接收深度和三个距离。默认串行计算，不需要 MPI，也不依赖本机 `test/` 目录。

## 1. 运行位移示例

先完成[安装](installation.md)，进入仓库根目录并激活环境：

```bash
python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional
```

Windows 非交互调用可使用：

```powershell
conda run -n pygrnwang python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional
```

脚本自动写出模型、预处理输入、写入 64 秒震源样本、运行 QSEIS2025、转换结果，并用 `output_type="disp"` 读取合成位移。示例使用 10 km 深度震源、地表接收、300/600/900 km 距离、4 秒采样间隔（Nyquist 频率 0.125 Hz），并启用平地球变换。底层格林函数库计算 4092 秒、1024 个采样点；合成后保存发震后 0–1020 秒（包含两端，共 256 点）的数组和绘图。求解器运行需要几分钟。模型取前 24 个数值行，常数 Qp=600、Qs=300 是教程设置，并不代表完整 AK135-F 衰减模型。

## 2. 查看结果

输出目录中包含 `disp.npz`、`disp.png`、`source_time_function.json`、`summary.json`、模型文件与 `library/`。NPZ 保存 0–1020 秒的波形、时间、距离、分量标签及单位，位移数组形状为 `(3, 3, 256)`；摘要记录环境、运行时间、输出形状和文件大小。`library/` 保留完整的 1024 点原生计算结果。

```{figure} ../_static/examples/qseis2025-regional.png
:alt: QSEIS2025 小算例在 300、600、900 公里处的三分量位移波形。
:width: 100%

10 km 深度震源在 300、600、900 km 处的合成位移，发震后 0–1020 秒；精确参数以共用脚本为准。
```

此示例的 `rotate=True` 对应 **东、北、上** 三分量。图中位移单位为米，震源矩为脚本中指定的 $10^{15}$ N m。图的横轴使用示例设定的时间零点；改变约化时间或平移选项后，应重新核对横轴含义。

震源采用 `wavelet_type=0` 和 1024 个自定义矩率样本：64 秒归一化 sin² 脉冲，预先补偿 QSEIS 的数值阻尼，使实际 STF 面积为 1、质心为 32 秒，与 SPGRN2020 一致。格林函数库保存速率核，读取程序在 `output_type="disp"` 时积分一次。

## 3. 计算应变与应力

在另一个目录中启用完整的教程输出：

```bash
python examples/qseis2025.py --regional --observables all --output-dir examples/output/qseis2025-regional-tensors
```

新增 `strain` 和 `stress` 的数组与图，脚本直接向读取程序请求位移、应变和应力，均保存 0–1020 秒；应变与应力数组形状为 `(3, 6, 256)`。启用地理旋转时，对称张量排列为 `[EE, EN, EU, NN, NU, UU]`，U 表示向上，也对应代码中的 Z。应变无量纲，应力单位为 Pa。不同后端未旋转张量的排列并不统一，详见[科学约定](../conventions.md)。

## 4. 复用已完成的计算

```bash
python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional --reuse
```

复用运行会另写 `summary-reuse.json`，保留首次计算的摘要。仅对同一模型、网格和输出设置复用。改变计算参数时使用新目录，以免将旧库误认为新参数的结果。

模型取前 24 个数值行形成分层半空间，与完整球形地球是不同近似，不能只凭相同距离判断结果应完全相等。参见[远距离教程](../backends/qseis2025.md#regional-waveforms-at-300-600-and-900-km)和[后端对比](backend-comparison.md)。

## 更短的入门计算

不加 `--regional` 时，同一脚本计算 30/60/90 km、0.5 秒采样的轻量算例：原生库为 127.5 秒、256 点，保存 0–100 秒（201 点），几秒内完成：

```bash
python examples/qseis2025.py --output-dir examples/output/qseis2025
```

## 完整脚本

中英文页面直接引用同一脚本：

```{literalinclude} ../../examples/qseis2025.py
:language: python
:linenos:
```

下一步可阅读 [QSEIS2025 后端指南](../backends/qseis2025.md)、[读取与后处理](../guides/reading.md)及[验证记录](../validation.md)。小算例用于学习完整流程；正式研究还需检查空间采样、频率范围和数值收敛。
