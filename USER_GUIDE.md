# AbacusCopilot 使用指南

> **版本**: v0.1.33 (2026-08-31)  
> **开发者**: Xu Yang (xuyangmark@foxmail.com)、Rong-yu Zhang  
> **简介**: AbacusCopilot 是 ABACUS DFT 软件的前后处理 CLI 工具包，灵感来源于 VASPKIT。

---

## 目录

1. [安装与配置](#1-安装与配置)
2. [快速开始](#2-快速开始)
3. [交互式菜单](#3-交互式菜单)
4. [命令行模式](#4-命令行模式)
5. [INPUT 文件生成 (任务 101–109)](#5-input-文件生成)
6. [STRU 结构文件 (任务 201–210)](#6-stru-结构文件)
7. [K 点生成 (任务 301–304)](#7-k-点生成)
8. [结构编辑 (任务 401–408)](#8-结构编辑)
9. [对称性分析 (任务 501–502)](#9-对称性分析)
10. [批处理 (任务 601–604)](#10-批处理)
11. [SCF 分析 (任务 701–703)](#11-scf-分析)
12. [能带结构 (任务 801–802)](#12-能带结构)
13. [态密度 DOS/PDOS (任务 901–903)](#13-态密度-dospdos)
14. [电荷密度 (任务 1001–1003)](#14-电荷密度)
15. [功函数 (任务 1101–1102)](#15-功函数)
16. [力学性质 (任务 1201–1202)](#16-力学性质)
17. [布居分析 (任务 1301)](#17-布居分析)
18. [MD 轨迹分析 (任务 3101–3108)](#18-md-轨迹分析)
19. [反应动力学 (任务 3301–3302)](#19-反应动力学)
20. [配置文件](#20-配置文件)
21. [常见问题](#21-常见问题)
22. [系统与配置 (任务 9901–9904, 9906)](#22-系统与配置)

---

## 1. 安装与配置

### 环境要求

- Python ≥ 3.9
- 可选依赖: ASE, spglib, seekpath, PyTorch (Apple Silicon GPU 加速)

### 安装

```bash
pip install -e .
# 或
bash setup.sh
```

### 初始化配置

```bash
abacuscopilot
# 首次运行自动进入 System Setup (任务 3201)
# 或手动: abacuscopilot -task 9901
```

配置项保存在 `~/.abacuscopilot/config.yaml`：

```yaml
libraries:
  # 支持单个目录或目录列表(按序查找,自动检测 PP-Orb/ 下的库,失效路径自动剔除)
  pseudo_library:
    - /path/to/pseudopotentials
    - /path/to/lanthanides
  orbital_library:
    - /path/to/orbitals
    - /path/to/lanthanides
paths:
  abacus_binary: abacus
  mpirun: mpirun
  sub_script: /path/to/submit_script.sh
defaults:
  kspacing: 0.14
  ecutwfc: 100.0
  calculation: scf
  basis_type: pw
```

赝势/轨道库放在项目根目录 `PP-Orb/` 下(如 `SG15-Version1p0_Pseudopotential`、
`SG15-Version1p0__StandardOrbitals-Version2p0`、`lanthanides-f--core.icmod1`),
首次运行自动检测。文件名大小写不敏感,容忍 `Sm3+_…` 这种带价态前缀的命名;
库支持递归查找(APNS 镧系包的 `{元素}/{基组}/` 嵌套布局可直接用)。当一个元素
有多个轨道时默认优先 DZP 基组 `4s2p2d1f` @ 7 au。

---

## 2. 快速开始

### 典型工作流

```bash
# 1. 从 CIF 生成 STRU
abacuscopilot -task 201                      # 交互式: 输入 CIF 路径

# 2. 生成 K 点
abacuscopilot -task 301                      # 自动 MP 网格

# 3. 生成 SCF INPUT
abacuscopilot -task 101 --basis pw          # PW basis SCF

# 4. 检查计算结果
abacuscopilot -task 703                      # 离子/电子步汇总

# 5. 分析 SCF 收敛
abacuscopilot -task 701                      # 能量收敛曲线

# 6. 画能带
abacuscopilot -task 801                      # 能带图

# 7. 画 DOS
abacuscopilot -task 901                      # 态密度
```

---

## 3. 交互式菜单

直接运行 `abacuscopilot` 进入交互式 TUI：

```
 =================== Structural Utilities ===================
   1)  INPUT File Generation
   2)  STRU File Generation
   3)  K-Point Generation
   ...

 =================== Electronic Utilities ===================
   7)  SCF Convergence Analysis
   8)  Band Structure Visualization
   ...
```

- 输入菜单编号进入子菜单
- 输入任务编号直接执行
- `0` 退出，`9` 返回上级

---

## 4. 命令行模式

```bash
abacuscopilot -task <任务号> [参数...]
abacuscopilot --list-tasks    # 列出所有任务
abacuscopilot --version       # 版本信息
abacuscopilot --clean         # 清理目录 (保留 STRU/INPUT/KPT)
```

**非交互模式**: 带 `--` 参数调用时跳过交互提示，使用模板默认值或命令行参数。

---

## 5. INPUT 文件生成

任务号 101–109，放在 `INPUT` 菜单下。

### 101 — SCF INPUT

生成自洽计算 INPUT。

```bash
abacuscopilot -task 101              # 交互式
abacuscopilot -task 101 --basis pw   # 非交互: PW basis
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw` | 基组类型 |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

### 102 — Relax INPUT

生成结构弛豫 INPUT。

```bash
abacuscopilot -task 102 --basis pw --relax-type cell-relax
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw` | 基组类型 |
| `--relax-type cell-relax\|relax` | 弛豫类型 |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

### 103 — MD INPUT

生成分子动力学 INPUT。

```bash
# DFT-PW NVT 系综
abacuscopilot -task 103 --basis pw --ensemble nvt

# DP NPT 100 kbar 各向同性
abacuscopilot -task 103 --basis dp --ensemble npt --press 100 --pmode iso

# LCAO Langevin
abacuscopilot -task 103 --basis lcao --ensemble langevin
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw\|dp` | 基组类型 (dp = Deep Potential) |
| `--ensemble nvt\|npt\|nve\|langevin\|fire\|msst` | MD 系综 |
| `--nstep N` | 总步数 (默认 10000) |
| `--dt N` | 时间步长 fs (默认 1.0) |
| `--tfirst N` | 初始温度 K (默认 300) |
| `--tlast N` | 最终温度 K (默认 300) |
| `--press N` | 目标压强 kbar (NPT 专用) |
| `--pmode iso\|aniso\|tri` | 控压模式 (NPT 专用) |
| `--dumpfreq N` | MD_dump 输出频率 (DP 专用, 默认 100) |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

**DP (Deep Potential) 模式**:
- 自动生成 STRU-dp (无赝势/轨道)
- 自动清理废弃的 upf/orb/KPT 文件
- 模型文件默认 `graph-compress.pb`
- 初始速度从 Maxwell-Boltzmann 分布随机生成

### 104 — Band Structure INPUT

```bash
abacuscopilot -task 104 --basis pw --proj    # PW + 投影能带
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw` | 基组类型 |
| `--proj` | 输出投影能带 (out_proj_band=1) |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

### 105 — DOS INPUT

```bash
abacuscopilot -task 105 --basis pw --pdos    # PW PDOS
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw` | 基组类型 |
| `--pdos` | 投影 DOS 而非总 DOS |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

### 106 — Work Function INPUT

```bash
abacuscopilot -task 106 --basis pw
```

### 110 — NEB INPUT

生成 NEB 过渡态计算 INPUT。

```bash
abacuscopilot -task 110 --basis pw --images 7
```

| CLI 参数 | 说明 |
|----------|------|
| `--basis lcao\|pw` | 基组类型 |
| `--images N` | NEB 中间像数量 (默认 5) |
| `--solver genelpa\|cusolver` | LCAO 求解器 |

**NEB 工作流**：
1. 准备初态 (`STRU_ini` 或 `init/STRU`) 和末态 (`STRU_fin` 或 `final/STRU`)
2. `abacuscopilot -task 3301` 或 `3302` 生成 NEB 路径
3. `abacuscopilot -task 110` 生成 INPUT
4. 将 INPUT/KPT 复制到每个 `00/` → `NN/` 目录及初末态目录

---

### 107 — Ecutwfc 收敛性测试

生成不同截断能的 INPUT 用于收敛性测试。

```bash
abacuscopilot -task 107    # 交互式
```

交互式提示输入起始/终止/步长 (默认 40/100/10 Ry)。自动生成 `ecutwfc_N/` 目录。

### 108 — Kspacing 收敛性测试

```bash
abacuscopilot -task 108    # 交互式
```

自动检测重复网格 (相同 kspacing 产生相同 k 点)，询问是否跳过。

### 109 — 收敛性分析

```bash
abacuscopilot -task 109    # 自动扫描 ecutwfc_*/ kspacing_*/
```

输出收敛曲线图 + 数据文件，能量自动归一化为 eV/atom (原子数从 log 读取)。

---

## 6. STRU 结构文件

任务号 201–210。

### 201 — CIF to STRU

```bash
abacuscopilot -task 201    # 交互式输入 CIF 路径
```

支持 ASE 可读的所有 CIF 格式。坐标类型可选 Direct (fractional) 或 Cartesian (Å)。

### 202 — POSCAR to STRU

```bash
abacuscopilot -task 202    # 输入 VASP POSCAR/CONTCAR
```

### 203 — Coord Convert

Direct ↔ Cartesian (Å) 坐标转换。输出 `STRU_Direct` 或 `STRU_Cartesian_angstrom`。

```bash
abacuscopilot -task 203
```

### 204 — STRU to CIF

```bash
abacuscopilot -task 204    # STRU → CIF
```

### 205 — STRU to POSCAR

```bash
abacuscopilot -task 205    # STRU → POSCAR
```

### 206 — STRU to PDB

```bash
abacuscopilot -task 206    # STRU → PDB (VMD/PyMOL)
```

### 207 — PDB to STRU(分子)

```bash
abacuscopilot -task 207 water.pdb   # PDB → STRU
```

孤立分子用(如 H2O)。PDB 自带 CRYST1 盒子则沿用;否则输入立方盒子边长(默认 15 Å),分子自动居中到盒子中心。

### 208 — STRU to LAMMPS

```bash
abacuscopilot -task 208    # STRU → LAMMPS data 文件
```

### 209 — LAMMPS to STRU

```bash
abacuscopilot -task 209    # LAMMPS data → STRU
```

### 210 — View Structure

在 ASE 3D 交互式窗口中可视化结构。

```bash
abacuscopilot -task 210          # 自动检测 STRU/POSCAR/CONTCAR/*.cif
abacuscopilot -task 210 STRU     # 指定 STRU 文件
abacuscopilot -task 210 POSCAR   # 指定 POSCAR 文件
```

**依赖**: `pip install ase pyqt5`

**功能**:
- 旋转/缩放/平移结构
- 周期性超胞显示 (`View → Repeat` 或 `Ctrl+R`)
- 切换配色方案 (`View → Colors → jmol`)
- 显示/隐藏键、晶胞边框、坐标轴
- 调整原子大小、键粗细
- 截图保存 (`File → Save Image`)

STRU 文件会自动转换为 POSCAR（临时文件，ASE 不支持 ABACUS 原生格式）。

---

## 7. K 点生成

任务号 301–304。

### 301 — Auto KPT (MP mesh)

自动 Monkhorst-Pack 网格。

```bash
abacuscopilot -task 301
```

默认 kspacing = 0.14 2π/Å，可通过配置修改。

### 302 — KPT (band path)

使用 seekpath 自动生成高对称 k 路径。

```bash
abacuscopilot -task 302
```

### 303 — Custom KPT

手动输入 k 点列表。

```bash
abacuscopilot -task 303
```

### 304 — KPT (phonon)

生成 phonopy 兼容的 `band.conf`。

```bash
abacuscopilot -task 304
```

---

## 8. 结构编辑

任务号 401–408。

### 401 — Build Supercell

```bash
abacuscopilot -task 401    # Nx×Ny×Nz 超胞
```

### 402 — Redefine Lattice

3×3 变换矩阵重定义晶胞。输出 `Redefined.STRU`。

```bash
abacuscopilot -task 402
```

### 403 — Coord Conversion

```bash
abacuscopilot -task 403    # Direct ↔ Cartesian (Å)
```

### 404 — Fix/Unfix Atoms

选择性固定/释放原子。

```bash
abacuscopilot -task 404
```

选项: Fix all / Unfix all / Fix by species / Fix by coordinate range。

### 405 — Slab Builder

Miller 指数切面。

```bash
abacuscopilot -task 405
```

参数: Miller (hkl), 原子层数, 真空层厚度 (Å)。

### 406 — Vacuum Layer

```bash
abacuscopilot -task 406
```

沿 a/b/c 轴方向添加真空。输出 `Vacuum_X_NN.N.STRU`。原子 Cartesian 坐标保持不变。

### 407 — Shift Atoms

沿指定方向平移所有原子。

```bash
abacuscopilot -task 407
```

支持 +a/-a/+b/-b/+c/-c 方向及自定义向量。输出 `Shifted.STRU`。

### 408 — Sort Atoms

按 x/y/z 坐标升序或降序排列原子 (保持物种分组)。

```bash
abacuscopilot -task 408
```

---

## 9. 对称性分析

任务号 501–502。

### 501 — Symmetry Analysis

使用 spglib 分析空间群和 Wyckoff 位置。

```bash
abacuscopilot -task 501
```

### 502 — Primitive Cell

将结构约化到原胞。

```bash
abacuscopilot -task 502
```

---

## 10. 批处理

任务号 601–604。

### 601 — PBS Script

```bash
abacuscopilot -task 601    # 生成 PBS/Torque 提交脚本
```

### 602 — SLURM Script

```bash
abacuscopilot -task 602    # 生成 SLURM 提交脚本
```

### 603 — Validate Inputs

```bash
abacuscopilot -task 603    # 检查 INPUT/STRU/KPT 完整性
```

### 604 — Batch Submit

```bash
abacuscopilot -task 604    # 批量提交子目录任务
```

---

## 11. SCF 分析

任务号 701–704。

### 701 — SCF Convergence

```bash
abacuscopilot -task 701
```

解析 `running_*.log`，显示步数、能量、收敛状态。可选画收敛图。

- 支持 v3.10 CU 格式 (`CU1`, `CU2`, ...)
- 支持 `E_KohnSham` 格式
- 能量自动识别 Ry/eV 单位

### 702 — SCF Compare

```bash
abacuscopilot -task 702    # 比较多个 SCF 结果
```

自动发现 `kspacing_*/` 或 `ecutwfc_*/` 目录，生成对比表格。

### 703 — Ion Steps

```bash
abacuscopilot -task 703    # 离子步汇总表 (relax/MD/SCF)
```

显示完成状态、收敛、能量、原子数。

### 704 — MD Monitor

```bash
abacuscopilot -task 704    # 实时监控 MD 模拟
```

实时显示 MD 步数、能量、温度。

---

## 12. 能带结构

任务号 801–802。

### 801 — Plot Band Structure

```bash
abacuscopilot -task 801 --bands BANDS_1.dat --kpt KPT --erange -5,5
```

| CLI 参数 | 说明 |
|----------|------|
| `--bands PATH` | BANDS 文件路径 |
| `--kpt PATH` | KPT 文件路径 |
| `--erange MIN,MAX` | 能量范围 (eV) |

### 802 — Fat-band Plot

```bash
abacuscopilot -task 802    # 投影能带图
```

---

## 13. 态密度 DOS/PDOS

任务号 901–903。

### 901 — Plot DOS

```bash
abacuscopilot -task 901
```

### 902 — Plot PDOS

```bash
abacuscopilot -task 902
```

### 903 — DOS+Bands Combined

```bash
abacuscopilot -task 903    # DOS + 能带双拼图
```

---

## 14. 电荷密度

任务号 1001–1003。

### 1001 — 1D Planar Avg Charge

```bash
abacuscopilot -task 1001   # z 方向一维平面平均电荷
```

### 1002 — Export Cube/XSF

```bash
abacuscopilot -task 1002   # CHGCAR → Cube/XSF
```

### 1003 — Diff. Charge Density

```bash
abacuscopilot -task 1003   # 差分电荷密度
```

Δρ = ρ(AB) − ρ(A) − ρ(B)。依次输入三个体系的电荷密度文件路径。

**输入格式支持两种**（自动识别）：

1. **Gaussian Cube 文本**（ABACUS `out_chg 1` 输出的 `OUT.*/SPIN*_CHG.cube`）。
2. **ABACUS 二进制 rhog 重启文件**（`OUT.*/ABACUS-CHARGE-DENSITY.restart`，或拷贝改名后的 `CHG`/`chg1` 等）。无需重新计算，直接读 SCF 生成的二进制电荷密度。

三个体系的 FFT 网格需一致（否则自动重采样到 ρ(AB) 的网格）。结果写入 `diff_charge_density.cube`，用 VESTA/Jmol 打开。

---

## 15. 功函数

> ⏳ **待开发 (coming soon)**：任务 1101–1102 功能可用性尚未验证，当前已隐藏（交互菜单进入 `11) Work Function Analysis` 会显示"此模块待开发"）。验证通过后重新开放。

任务号 1101–1102（预留）。

### 1101 — Work Function

```bash
abacuscopilot -task 1101   # 功函数计算（待开发）
```

### 1102 — Macroscopic Avg

```bash
abacuscopilot -task 1102   # 宏观平均法功函数（待开发）
```

---

## 16. 力学性质

任务号 1201–1202。

### 1201 — Elastic Constants

```bash
abacuscopilot -task 1201   # 弹性常数 (需 ELASTIC_CONSTANTS 输出)
```

### 1202 — EOS Fitting

```bash
abacuscopilot -task 1202   # 状态方程拟合 (Birch-Murnaghan)
```

---

## 17. 布居分析

任务号 1301。

### 1301 — Mulliken Analysis

```bash
abacuscopilot -task 1301
```

> 注 1:键级分析(任务 **1401 Mulliken Bond Order**,即 Mulliken 重叠布居类型,来自 ABACUS `out_mul 1`)已独立为主菜单 **`14) Bond Order`** 模块。后续 Mayer/Wiberg 等其他键级类型可在 1402 之后追加。
>
> 注 2:Lowdin(原 1302)已移除——ABACUS v3.10 不输出 Lowdin 数据,无法分析。

### Bader(1304)/ Hirshfeld(1305)/ RESP(1306)均已隐藏

> **Bader (1304)**:ABACUS `SPIN1_CHG.cube` 数据是 z-fastest 循环,读取器曾按标准 x-fastest 解析导致密度错位、Bader 电荷错误(水 O=8/H=0);转置修复后网格仍难分辨 H 小盆地。下架待修。
>
> **Hirshfeld (1305)**:ABACUS v3.10 不支持 `out_hirshfeld`,永远找不到数据。
>
> **RESP (1306)**:多原子分子拟合病态。三者待后续版本修复后再开放。

---

## 18. MD 轨迹分析

任务号 3101–3108。支持 ABACUS MD_dump、LAMMPS dump、VASP XDATCAR 三种格式。

### 3101 — MD Trajectory → PDB

```bash
abacuscopilot -task 3101
```

将 MD_dump 转为 VMD 可读的 PDB 格式。支持单帧/全部帧/范围导出。

### 3102 — Extract Frames

```bash
abacuscopilot -task 3102   # 按步长采样
```

### 3103 — MSD

```bash
abacuscopilot -task 3103
```

均方根位移 + 扩散系数。支持 x/y/z 方向分离。计算使用 GPU 加速 (macOS M-chip PyTorch MPS，自动检测)。

### 3104 — RDF

```bash
abacuscopilot -task 3104
```

径向分布函数。可选计算配位数 CN(r)。支持双 Y 轴合并图或分别作图。

### 3105 — Probability Density

```bash
abacuscopilot -task 3105
```

3-D 概率密度 → CHGCAR 格式。VESTA 可直接打开看等值面。

### 3106 — van Hove

```bash
abacuscopilot -task 3106
```

van Hove 关联函数 (Gs + 方向分量 + Gd + NGP)。GPU 加速。

输出:
- `vanHove_X_Gs_heatmap.png` + `vanHove_X_Gd_heatmap.png` — 热力图，用户自定义颜色范围
- `vanHove_X_Gs_X/Y/Z_heatmap.png` — 方向自关联 P(|Δx|,t)、P(|Δy|,t)、P(|Δz|,t)，用于分析 x/y/z 方向离子 hopping
- `vanHove_X_NGP.png` — 非高斯参数
- `vanHove_X_Gs_slices.png` — 分时切片 (可选平滑)

`vanHove_X.npz` 同时保存全部结果 (r, time_lags, gs, gsx/gsy/gsz, gd, ngp) 及计算参数 (r_max, dr, stride, 帧范围)。再次运行时若当前目录检测到有效缓存，会交互式询问是否复用——选择复用则**跳过漫长的轨迹读取**，直接从缓存重建所有图 (可继续调整颜色范围)。旧版本 (无方向分量的) npz 会自动忽略并重算。

### 3107 — LAMMPS → MD_dump

```bash
abacuscopilot -task 3107 new.dump
```

将 LAMMPS 轨迹 dump 转为 ABACUS MD_dump 格式。自动检测原子类型，交互式映射元素符号。

- LAMMPS "real" units: Å, fs
- 原子按 LAMMPS ID 排序确保跨帧一致性

### 3108 — XDATCAR → MD_dump

```bash
abacuscopilot -task 3108 XDATCAR
```

VASP XDATCAR 转 ABACUS MD_dump。

- 自动识别 NVT (固定晶胞) / NPT (每帧晶胞变化)
- Direct 坐标自动转 Cartesian (Å)

---

## 22. 系统与配置

任务号 9901–9904、9906。在交互式菜单中显示为 `99) System & Configuration`。

| 任务 | 说明 |
|------|------|
| 9901 | System Setup — 配置路径和默认值 |
| 9902 | Show Config — 显示当前配置 |
| 9903 | Check Environment — 检查环境和依赖 |
| 9904 | Clean Directory — 清理 (保留 INPUT/STRU/KPT/upf/orb) |
| 9906 | Set Submit Script — 配置提交脚本路径 |

---

## 19. 反应动力学

任务号 3301–3302。NEB 路径生成。

### 3301 — NEB Path (Linear)

线性插值生成 NEB 中间像。

```bash
abacuscopilot -task 3301
```

自动检测 `init/STRU` 和 `final/STRU`（或 `STRU_ini`/`STRU_fin`、`POSCAR_ini`/`POSCAR_fin`）。计算最大原子位移并建议图像数（d_max / 0.8 Å）。

输出：
- `00/` → `NN/` — 每个中间像一个子目录
- `trj.STRU` + `trj.vasp` — 合并视图（移动原子多帧展示，静态原子只出现一次）
- `path_Nframes.traj` — ASE 动画文件（用 `abacuscopilot -task 206` 打开）

### 3302 — NEB Path (IDPP)

IDPP 原子对距离优化插值。先用线性插值创建初始路径，再用 `ase.mep.idpp_interpolate` 优化键长畸变。

```bash
abacuscopilot -task 3302
```

参数和输出同 3301。IDPP 避免原子碰撞，产生物理合理的中间结构。

### 输出格式选择

两个任务均支持选择输出格式：
- `STRU (ABACUS)` — 默认，ABACUS 原生格式
- `POSCAR (VASP)` — VASP 格式

---

## 20. 配置文件

`~/.abacuscopilot/config.yaml` 完整示例:

```yaml
defaults:
  pseudo_dir: ./
  orbital_dir: ./
  kspacing: 0.14
  ecutwfc: 100.0
  scf_thr: 1.0e-07
  force_thr_ev: 0.01
  calculation: scf
  basis_type: pw
  dft_functional: pbe
  mixing_beta: 0.8
  smearing_sigma: 0.015
  smearing_method: gauss
  relax_method: cg
  relax_nmax: 60
  ks_solver: genelpa

plotting:
  style: abacuscopilot
  dpi: 300
  figure_format: png
  figure_size: [8, 6]
  color_cycle: tab10
  font_size: 12
  show_fermi: true

paths:
  abacus_binary: abacus
  mpirun: mpirun
  sub_script: /home/user/bin/sub.abacus

libraries:
  pseudo_library:
    - /home/user/PP/SG15_ONCV_v1.0_upf
    - /home/user/PP/lanthanides-f--core.icmod1
  orbital_library:
    - /home/user/Orb/SG15_StandardOrbitals_v2.0
    - /home/user/PP/lanthanides-f--core.icmod1

user_presets: {}
```

---

## 21. 常见问题

### Q: "No log file found" 但计算已完成

OUT.ABACUS 在当前目录的子目录中。检查是否有 `ABACUS-CHARGE-DENSITY.restart` 文件导致误判 (已修复)。NSCF 计算会将该文件拷到工作目录，程序已改为优先使用 `OUT.*/` 子目录中的日志。

### Q: MSD/RDF/van Hove 很慢

- 使用 frame interval 和 time stride 减少抽样
- macOS M-chip 安装 PyTorch: `pip install torch` (自动 MPS GPU 加速)
- van Hove 的 Gd 计算最有 GPU 加速效果

### Q: 热力图中看不到颜色变化

颜色范围太大导致细节被冲淡。程序会显示数据范围，输入 0 使用最大值 1/40 作为自动上限，或手动输入合适的 vmax 值。

### Q: STRU-dp 生成出来原子全是空的

STRU 文件的 `ATOMIC_SPECIES` 段至少要有两列 (元素名 质量)。如果只有元素名没有质量和赝势路径，旧版 `read_stru` 会跳过去 (已修复)。

### Q: 命令行 `--pw` 没生效

使用 `--basis pw` 而非 `--pw`。程序会检测这种错误并提供正确语法提示。

### Q: 不认识的参数被静默忽略

已修复。现在 `parse_known_args` 的 unknown args 会输出 warning 并提示使用 `--help`。
