# AbacusCopilot 使用指南

> **版本**: v0.1.37 (2026-09-12)  
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
  # 支持单个目录或目录列表,按序查找。自动检测 PP-Orb/ 下声明的系列;任一路径
  # 失效(如升级后目录搬家)会重新探测并与幸存的路径合并,不会整表丢失。
  # 这两个列表是“当前生效”的目录:当你在下面选定某个系列(family)时会自动重指。
  pseudo_library:
    - .../PP-Orb/SG15-Version1p0/SG15-Version1p0_Pseudopotential
    - .../PP-Orb/lanthanides-f--core.icmod1/PD04.3+f--core.icmod1
  orbital_library:
    - .../PP-Orb/SG15-Version1p0/SG15-Version1p0__StandardOrbitals-Version2p0
    - .../PP-Orb/lanthanides-f--core.icmod1/PD04.3+f--core.icmod1

  # 当前生效的赝势/轨道“系列”: sg15(默认) | apns | apns/<variant>
  # | dojoncfr | dojoncfr/<sz|dzp|tzdp> | custom。
  # 在 INPUT / 全计算配置流程中交互选择,持久化为全局默认(可随时切回 SG15)。
  family: sg15

  # 当某系列把“同一套基组”按多个截断半径发行时(如 Dojo-NC-FR 每个 tier
  # 提供 6-12 au),选哪个 rcut 副本。交互选择时会被问一次,值持久化在这里;
  # 仅 dojoncfr 系列生效。
  #   "7"   (默认,推荐) 最接近内置 SG15 轨道通行的 7 au
  #   "min"              最小 rcut(最快/最紧凑)
  #   "max"              最大 rcut(最完备/最接近平面波极限)
  rcut_policy: "7"

  # 用户自配的外部系列(不打进软件包)。把你自己下载的目录填到下面
  # (或让 INPUT 的“库家族”向导帮你填写),生成文件时选该系列即可。
  # abacuscopilot 不会跨系列混用文件。
  families:
    # apns:  # ABACUS-APNS-PPORBs-v1
    #   pseudo_dir: /path/to/ABACUS-APNS-PPORBs-v1/apns-pseudopotentials-v1
    #   orbital_dirs:
    #     efficiency: /path/to/.../apns-orbitals-efficiency-v1
    #     precision:  /path/to/.../apns-orbitals-precision-v1
    #
    # dojoncfr:  # Dojo-NC-FR(全相对论 NC 赝势+轨道, 用于 SOC)
    #   pseudo_dir: /path/to/Dojo-NC-FR/Pseudopotential
    #   orbital_dirs:  # 三个 tier(SZ/DZP/TZDP)共用一个 Orbitals_v2.0 根目录
    #     sz:   /path/to/Dojo-NC-FR/Orbitals_v2.0
    #     dzp:  /path/to/Dojo-NC-FR/Orbitals_v2.0
    #     tzdp: /path/to/Dojo-NC-FR/Orbitals_v2.0
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

赝势/轨道库放在项目根目录 `PP-Orb/` 下,**一个顶层目录 = 一个系列**,系列内部再分
赝势与轨道:

```
PP-Orb/
├── SG15-Version1p0/                                  # 默认系列
│   ├── SG15-Version1p0_Pseudopotential/              #   .upf ×69
│   └── SG15-Version1p0__StandardOrbitals-Version2p0/ #   .orb ×68
├── lanthanides-f--core.icmod1/                       # 4f 补充系列
│   └── PD04.3+f--core.icmod1/<元素>/                 #   .upf + .orb ×14 (Ce–Lu)
├── ABACUS-APNS-PPORBs-v1/                            # 外部系列(不随包发布)
├── Dojo-NC-FR/                                       # 外部系列(不随包发布)
└── README.md
```

首次运行自动检测。`PP-Orb/` 内附 `README.md`:把你的库(含 `*.upf` / `*.orb` 的
文件夹)整体放进去即可被自动检测——发布版 tarball 会随包内置上面的 SG15 与镧系两个
系列,git 仓库只保留这份说明文件。**两个系列的化学元素零重叠**(SG15 覆盖不到 Ce–Lu,
镧系只有那 14 个),所以镧系常驻检测范围也不会抢 SG15 的任何元素。配置里登记的库路径
**任何一条失效**(例如升级换了布局)都会自动重新探测并与幸存的路径合并——不会让整个
SG15 列表静默消失;反过来,全部路径仍然有效时你的配置一个字都不动,此时新丢进来的顶层
库需要重选一次家族才会纳入。文件名大小写不敏感,容忍 `Sm3+_…` 这种带价态前缀的命名,
也支持 APNS 里 Hf/Os 的半芯命名 `Hf-sp.PD04.PBE.UPF` / `Os-sp.PD04.PBE.UPF`
(`-sp` 分隔符)。库支持递归查找(APNS 镧系包的 `{元素}/{基组}/` 嵌套布局可直接用)。
当一个元素有多个轨道时,按所在系列取**确定性默认档**: SG15 优先 DZP `4s2p2d1f` @
7 au;ABACUS-APNS 的 `efficiency` 取更小 rcut(Cs→`10au`),`precision` 取最完备基组
(B→`4s4p3d2f`,K/Cs/Na/Rb→`5s4p3d2f`,Sb→`4s4p4d3f2g`);Dojo-NC-FR 尊重所选
tier(`family: dojoncfr/<tier>`),只在对应 `{元素}_{SZ|DZP|TZDP}` 子目录里查找,
某 tier 缺失时**绝不**静默拿其它 tier 顶替;同一 tier 内默认取最接近 7 au 的 rcut
(DZP Hf→`Hf_DZP/Hf_gga_7au_100Ry_4s2p2d1f.orb`)。由于 Dojo 把同一基组按多个
rcut 发行,选择 Dojo 家族时会被问一句“用哪个 rcut 副本”——选择结果存为
`libraries.rcut_policy`(`"7"`|`"min"`|`"max"`,见上),7 au 那档标为推荐。

某元素在**当前生效的系列里**没有文件时(最典型的例子:`La` 在 SG15 下有赝势、
没有轨道),写 STRU 的警告除了列出缺失项,还会**点名哪个系列有**这个文件并给出
文件名(例如 `Tip: ABACUS-APNS-PPORBs-v1 has La_gga_9au_100Ry_4s2p2d1f.orb`)。
这是纯提示,不会自动切换:要换就回到流程开头的“赝势/轨道库家族”提示里**整体**
换一个系列。**不要**把同一个元素的赝势与轨道拆到两个系列——数值原子轨道是为特定
赝势生成的赝波函数,两者必须同源,混用得到的结果没有意义;所以本工具也不提供
“只给某个元素换系列”的入口。

**外部库系列(如 ABACUS-APNS-PPORBs-v1、Dojo-NC-FR)。** 其它系列**不随软件包分发**。
把下载好的目录整个丢进 `PP-Orb/`,然后在 INPUT / 全计算配置流程开头的“赝势/轨道库家族”
提示里选择它(这一步会自动登记到 `libraries.families.<id>`);也可以先手工登记再选。
默认回车保持 SG15(内置),选中后该选择持久化为全局默认,可随时切回。
**Dojo-NC-FR** 提供全相对论(NC,`relativistic="full"`,
`has_so=1`)赝势及配套轨道——做自旋轨道耦合(SOC,`lspinorb 1`)就该用它;普通非 SOC
计算同样可用(ABACUS 会自动退化为标量相对论)。当生效的是已登记的外部系列而某元素在其
库中缺失时,会**报错终止并提示换用其它系列**(绝不跨系列混用文件);SG15 下仍沿用原有
的黄字软提醒。凡是登记在 `libraries.families` 里的目录都会被**按路径**从 SG15 检测中
排除;而上面两个已知外部系列**按名字、在登记之前**就已排除——它们的文件与 SG15 完全
同名(Dojo 甚至有个裸 `Li.upf`,字典序还排在 SG15 的 `Li_ONCV_PBE-1.0.upf` 前面),混进来
会让 SG15 作业**静默**选错赝势。所以把这两个目录放在 `PP-Orb/` 下也无害,切回 SG15
绝不会混入其它系列。

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

> **库系列选择。** 交互式生成 INPUT (101–109) 时，开头会询问"赝势/轨道库家族"
> (当前 SG15 内置)。回车保持当前系列;输入序号可切换到其它系列 (如
> ABACUS-APNS-PPORBs-v1——选中后再选 efficiency/precision 轨道档;Dojo-NC-FR——SOC 用,
> 选中后再选 SZ/DZP/TZDP)。**只要该系列的目录已经放在 `PP-Orb/` 下,菜单里就会直接
> 列出并自动登记**(如 `Dojo-NC-FR  (found in PP-Orb — will auto-register)`),
> 不必事先手工写 config。选择持久化为全局默认，可随时切回 SG15。
> 若当前系列缺少某元素，自动补文件会**报错终止**并建议换系列 (绝不跨库混用);
> SG15 下仍为黄字软提醒。详见 §1。

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

默认 kspacing = 0.14（单位 1/bohr，ABACUS 约定，值越小 k 点网格越密），可通过配置修改。

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

任务号 1101–1102 已开放。两项任务直接读取 ABACUS 原始
`ElecStaticPot.cube`：从 cube 中的原子坐标识别最大的周期性无原子真空间隙，
在排除表面邻近区域后得到 `V_vacuum`；再从 `running_scf.log` 读取最终
`E_Fermi`，按 `Φ = V_vacuum − E_Fermi` 计算功函数。整个流程不依赖 zstar。

### 1101 — Work Function

```bash
abacuscopilot -task 1101   # 功函数计算
# 也可显式指定：--file ElecStaticPot.cube --log OUT.ABACUS/running_scf.log
```

### 1102 — Macroscopic Avg

```bash
abacuscopilot -task 1102   # 周期双窗口宏观平均法功函数
# --period 可指定平滑周期（Å）；默认使用 c 轴长度的四分之一
```

### 真实 ABACUS 输出验证

仓库提供了一个不携带大体积计算结果的可选端到端测试。准备一个已经完成的
ABACUS SCF 目录（需要 `OUT.ABACUS/ElecStaticPot.cube` 和
`OUT.ABACUS/running_scf.log`），然后设置目录并运行：

```bash
export ABACUS_WORKFUNC_CASE=/path/to/completed_scf
python -m pytest -q tests/test_workfunc_real_case.py
```

该测试会同时检查原始 cube 真空能级路径和 1102 周期宏观平均路径；未设置环境变量时自动跳过。

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

### 3101 — Trajectory → PDB

```bash
abacuscopilot -task 3101
```

将 ABACUS MD_dump 或 ASE `.traj`（MD/NEB 链，如 `converged.traj`/`path_*.traj`）转为 VMD 可读的 PDB 格式。支持单帧/全部帧/范围导出。

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
| 9901 | System Setup — 配置 ABACUS 二进制路径与计算默认值(不再询问赝势/轨道目录) |
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

配置文件是 `~/.abacuscopilot/config.yaml`,首次运行(或 `-task 9901` 配置向导)
时自动生成。

每个值右边带一句短语说明(样式仿 VASPKIT 的 `~/.vaspkit`):**注释不占行、
不遮值**,一行一个参数。那句话说删就删,删了不影响任何东西;本节表格才是完整
说明,缩写放不下的都写在这里。

生成出来的完整样子(所有键都在;路径按你机器上的实际探测结果):

```yaml
# abacuscopilot configuration
#
# Plain YAML -- edit any value in place. The notes on the right are for
# reference only; deleting them changes nothing. Longer descriptions and
# the meaning of every key: USER_GUIDE.md, section 20.
#
# Generated files are only written when the file is missing or you run the
# wizard (task 9901); reading the config never rewrites it.

# Defaults for generated ABACUS INPUT files -- plain ABACUS keywords.
defaults:
  pseudo_dir: ./          # pseudopotential dir; './' = job dir (matching files are copied in)
  orbital_dir: ./         # orbital dir; './' = job dir (matching files are copied in)
  kspacing: 0.14          # k-spacing in 1/Bohr; ignored when a KPT file is present
  ecutwfc: 100.0          # plane-wave cutoff, Ry
  scf_thr: 1.0e-07        # SCF convergence threshold
  force_thr_ev: 0.01      # force convergence, eV/Angstrom
  calculation: scf        # default task: scf | relax | cell-relax | md | ...
  basis_type: lcao        # lcao (numeric orbitals) | pw (plane waves)
  dft_functional: pbe     # exchange-correlation functional
  mixing_beta: 0.8        # charge mixing; lower it (0.2) when SCF will not converge
  smearing_sigma: 0.015   # smearing width, Ry
  smearing_method: gauss  # smearing scheme
  relax_method: cg        # relaxation algorithm
  relax_nmax: 60          # maximum ionic steps
  ks_solver: genelpa      # eigenvalue solver

# Figure defaults.
plotting:
  style: abacuscopilot    # figure style; names are defined in plotting/style.py
  dpi: 300                # figure resolution
  figure_format: png      # png | pdf | eps | svg | ...
  figure_size: [8, 6]     # figure size in inches, [width, height]
  color_cycle: tab10      # NOT USED YET -- nothing in the code reads this
  font_size: 12           # base font size
  show_fermi: true        # NOT USED YET -- nothing in the code reads this

# External programs abacuscopilot runs.
paths:
  abacus_binary: abacus   # ABACUS executable; a bare name is looked up on $PATH
  mpirun: mpirun          # MPI launcher, used for parallel runs
  abacus_source: ''       # ABACUS source tree, for the abacuslite Python interface
  slurm_env_file: ''      # sourced at the top of every SLURM job
  sub_script: ''          # sbatch template, copied into the job dir
  sub_script_dp: ''       # as above for ABACUS-DP jobs; empty or missing = use sub_script
  deepmd_python: ''       # NOT USED YET -- nothing in the code reads this
  abacus_dp_binary: ''    # ABACUS-DP binary; empty = probe abacus_binary
  abacus_dp_version: ''   # DeepMD-kit version it was built against -- NOT ABACUS's

# Pseudopotentials and orbitals.
libraries:
  pseudo_library:         # searched in order; dirs that no longer exist are dropped
    - /home/user/PP/SG15_ONCV_v1.0_upf
    - /home/user/PP/lanthanides-f--core.icmod1
  orbital_library:        # as above, for *.orb
    - /home/user/Orb/SG15_StandardOrbitals_v2.0
    - /home/user/Orb/lanthanides-f--core.icmod1
  family: sg15            # sg15 | apns | apns/<variant> | dojoncfr[/<tier>] | custom
  rcut_policy: '7'        # '7' | 'min' | 'max'; dojoncfr only -- quote it
  families: {}            # external series you downloaded yourself, keyed by id
```

**这份文件什么时候被写。** 只有两种情况:文件不存在时自动生成,以及你跑配置
向导(`-task 9901`)时。**读取配置从不回写** —— 你自己加的注释、你排的版式、
你写的值,加载过程一概不动。发布版本新加的键会先在内存里生效,要让它出现在
文件里,重跑一次 9901 即可(键名和含义本节表格里都有)。

### 20.1 `defaults` —— 生成 INPUT 时的默认参数

这些是 ABACUS 关键字本身,只是预先存了一份,生成 INPUT 时用得上 —— 交互流程里
你随时可以当场改。

| 键 | 默认 | 含义 |
|----|------|------|
| `pseudo_dir` | `./` | 赝势文件目录。`./` 表示作业目录 —— 匹配到的赝势会被拷进作业目录。9901 向导**不再询问**此项;若旧值指向一个已不存在的路径,读取配置时会自动重置回 `./` |
| `orbital_dir` | `./` | 轨道文件目录,同上 |
| `kspacing` | `0.14` | k 点间距,单位 1/Bohr。作业目录里已有 KPT 文件时不生效 |
| `ecutwfc` | `100.0` | 平面波截断能(Ry) |
| `scf_thr` | `1e-7` | SCF 收敛阈值 |
| `force_thr_ev` | `0.01` | 受力收敛判据(eV/Å),弛豫类任务用 |
| `calculation` | `scf` | 默认计算类型:`scf` / `relax` / `cell-relax` / `md` … |
| `basis_type` | `lcao` | 基组:`lcao`(数值原子轨道)或 `pw`(平面波) |
| `dft_functional` | `pbe` | 交换关联泛函 |
| `mixing_beta` | `0.8` | 电荷密度混合系数。难收敛时调小(如 0.2) |
| `smearing_sigma` | `0.015` | 展宽宽度(Ry) |
| `smearing_method` | `gauss` | 展宽方式 |
| `relax_method` | `cg` | 结构弛豫算法 |
| `relax_nmax` | `60` | 最大离子步数 |
| `ks_solver` | `genelpa` | 本征值求解器 |

### 20.2 `plotting` —— 出图默认值

| 键 | 默认 | 含义 |
|----|------|------|
| `style` | `abacuscopilot` | 样式名。除 `abacuscopilot` 外可写内置样式名(见 `plotting/style.py` 的 `COLOR_SCHEMES`) |
| `dpi` | `300` | 出图与保存的 DPI |
| `figure_format` | `png` | 保存格式(png / pdf / eps / svg …) |
| `figure_size` | `[8, 6]` | 图尺寸,单位英寸 |
| `color_cycle` | `tab10` | ⚠️ **目前不生效** —— 没有任何代码读它 |
| `font_size` | `12` | 字号 |
| `show_fermi` | `true` | ⚠️ **目前不生效** —— 没有任何代码读它 |

### 20.3 `paths` —— 外部程序

| 键 | 默认 | 含义 |
|----|------|------|
| `abacus_binary` | `abacus` | ABACUS 可执行文件。写裸名字走 `$PATH`,写绝对路径按原样用 |
| `mpirun` | `mpirun` | MPI 启动器,并行运行时用 |
| `abacus_source` | `''` | ABACUS 源码树路径。用 abacuslite 的 Python 接口(`interfaces/ASE_interface/`)时需要 |
| `slurm_env_file` | `''` | SLURM 作业开头 `source` 的环境脚本(CUDA / 编译器 / `module load`) |
| `sub_script` | `''` | 普通 ABACUS 作业的 sbatch 模板,会随 INPUT 一起拷进作业目录 |
| `sub_script_dp` | `''` | ABACUS-DP(`esolver_type dp`)作业的模板;**留空或文件不存在则回退到 `sub_script`** |
| `deepmd_python` | `''` | ⚠️ **目前不生效** —— 没有任何代码读它 |
| `abacus_dp_binary` | `''` | ABACUS-DP 二进制。留空 = 去探测 `abacus_binary` |
| `abacus_dp_version` | `''` | 上面这个二进制编译时链接的 **DeepMD-kit 版本**(如 `'2.2.11'`),**不是 ABACUS 的版本** |

**关于 `abacus_dp_binary`,值得单独说清楚。** DP 支持是 ABACUS **编译期**决定的:
一个二进制要么带 DP、要么不带,运行时换不了。机器上 DP 版和非 DP 版是两个不同
编译产物并存时(例如 `abacus-LTSv3.10.0-dp` 与 `abacus-LTSv3.10.0`),必须在这里
指定 DP 那个 —— 这个键决定了能读哪些模型格式(`.pb` / `.pth` / `.ckpt`)。

`abacus_dp_version` 覆盖自动探测,**只在该探测不可靠时手填**(例如静态链接、
二进制里没有留版本串)。写的时候注意加引号:`abacus_dp_version: '2.10'` ——
不加引号会被 YAML 读成浮点数 `2.1`,版本号就悄悄变了。

### 20.4 `libraries` —— 赝势与轨道

| 键 | 默认 | 含义 |
|----|------|------|
| `pseudo_library` | 自动探测 | 赝势搜索目录列表,按序查找。首次运行时按 `PP-Orb/` 下声明的系列探测;**任一路径失效即重新探测并与幸存的路径合并**——全部路径仍有效时你的配置一个字都不动 |
| `orbital_library` | 自动探测 | 轨道文件(`.orb`)搜索目录,同上 |
| `family` | `sg15` | 当前生效的库“系列” |
| `rcut_policy` | `"7"` | 同一套基组有多个 rcut 副本时选哪个,**仅 `dojoncfr` 系列生效** |
| `families` | `{}` | 用户自备的外部系列(不随包分发) |

`family` 的取值:

| 取值 | 说明 |
|------|------|
| `sg15` | 默认,随包内置 |
| `apns` / `apns/<variant>` | ABACUS-APNS-PPORBs 系列 |
| `dojoncfr` / `dojoncfr/<sz\|dzp\|tzdp>` | Dojo-NC-FR 全相对论系列,做 SOC 用 |
| `custom` | 由 `families` 里注册的内容决定 |

`rcut_policy` 的取值:`"7"`(默认,取最接近 7 au 的)/ `"min"` / `"max"`。
**注意要写引号** —— `7` 会被 YAML 读成整数。

`families` 是用户自己下载的系列,在 INPUT / STRU 流程里选择时会用到:

```yaml
libraries:
  families:
    apns:                       # ABACUS-APNS-PPORBs-v1
      pseudo_dir: /path/to/ABACUS-APNS-PPORBs-v1/apns-pseudopotentials-v1
      orbital_dirs:
        efficiency: /path/to/.../apns-orbitals-efficiency-v1
        precision:  /path/to/.../apns-orbitals-precision-v1
    dojoncfr:                   # Dojo-NC-FR(全相对论 NC 赝势+轨道,用于 SOC)
      pseudo_dir: /path/to/Dojo-NC-FR/Pseudopotential
      orbital_dirs:             # SZ/DZP/TZDP 共用一个 Orbitals_v2.0 根目录
        sz:   /path/to/Dojo-NC-FR/Orbitals_v2.0
        dzp:  /path/to/Dojo-NC-FR/Orbitals_v2.0
        tzdp: /path/to/Dojo-NC-FR/Orbitals_v2.0
```

细节(自动探测规则、文件名匹配、各系列选哪个基组)见 `README.md` 的
Configuration 一节和 `PP-Orb/README.md`。

> **本文件不认识的键不会被删掉。** 你另外加的顶层键(例如自己的脚本用的配置)
> 会原样保留 —— abacuscopilot 只在读的时候做合并,不会把文件重写掉。

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
