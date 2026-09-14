# PP-Orb/ — 赝势 / 数值轨道库目录
## Pseudopotential & numerical-orbital library directory

这里存放 abacuscopilot 在生成 INPUT / 自动同步 STRU 时按元素查找并拷贝的
**赝势（PP，`*.upf`）** 与 **数值原子轨道（NAO，`*.orb`）** 库树。

**一个顶层目录 = 一个系列（series）**，系列内部再分赝势与轨道。默认只用 SG15
系列；遇到 4f 元素时镧系补充系列会被自动用上（见下）。

- 发布版的 release tarball 随包内置两个系列（SG15 与镧系 f-core），开箱即用。
- git 仓库本身**不**携带这些大体积库，只保留本说明文件。

### 目录形状

```
PP-Orb/
├── SG15-Version1p0/                                  # 默认系列
│   ├── SG15-Version1p0_Pseudopotential/              #   .upf ×69
│   └── SG15-Version1p0__StandardOrbitals-Version2p0/ #   .orb ×68
├── lanthanides-f--core.icmod1/                       # 4f 补充系列
│   └── PD04.3+f--core.icmod1/<元素>/                 #   .upf + .orb ×14 (Ce–Lu)
├── ABACUS-APNS-PPORBs-v1/                            # 外部系列（不随包发布）
├── Dojo-NC-FR/                                       # 外部系列（不随包发布）
└── README.md
```

### 规则

- **默认查询范围**：SG15 系列的两个子目录 + 镧系补充系列，外加 `PP-Orb/` 下
  任何**未知**的顶层文件夹。把你手头含 `*.upf` / `*.orb` 的目录整个复制进来就会
  被自动探测到，无需改配置。
- SG15 与镧系的**元素零重叠**（SG15 覆盖不到 Ce–Lu 的 4f；镧系只有那 14 个元素），
  所以镧系常驻列表也不会抢 SG15 的任何元素。
- **已知外部系列例外**：`ABACUS-APNS-PPORBs-v1`、`Dojo-NC-FR` 这两个顶层目录
  **无论是否登记进配置，都绝不会并入上面两个默认列表**——它们含有与 SG15
  同名的文件（例如 Dojo 的裸 `Li.upf`，字典序还排在 `Li_ONCV_PBE-1.0.upf` 前面），
  混进来会让 SG15 作业**静默**选错赝势与轨道。要用它们请从家族提示里选择。
- 元素→文件匹配**大小写不敏感**，容忍价态前缀
  （`Sm3+_f--core-icmod1.PD04.PBE.UPF`）与半芯 `-sp` 命名
  （`Hf-sp.PD04.PBE.UPF` / `Os-sp.PD04.PBE.UPF`）。
- 一个元素对应多个轨道时，按当前**库家族（family）**取确定性默认档
  （SG15：DZP `4s2p2d1f`@7au；APNS：`efficiency` 更小 rcut、`precision` 最完备基组）。
- 配置里登记的库路径**任何一条失效**（例如升级后目录搬家）都会自动重新探测，
  并与仍然有效的路径合并——不会让整个 SG15 列表静默消失。全部路径仍然有效时，
  你的配置一个字都不动；此时新丢进来的顶层库需要重选一次 **SG15** 家族
  （或手工在配置里登记新目录）才会纳入。
- `La` 在 LCAO 下**只有赝势、没有轨道**（SG15 有 La 赝势，两个系列都不含 La
  轨道；只有 `ABACUS-APNS-PPORBs-v1` 有）——会如实报缺，而不是拿别的系列的
  轨道顶上。
- 某元素在生效的系列里缺文件时，写 STRU 的警告会**点名哪个系列有**该文件
  （`Tip: ABACUS-APNS-PPORBs-v1 has La_gga_9au_100Ry_4s2p2d1f.orb`）。这是
  纯提示：要换就在流程开头的**库家族提示**里整体换系列，不会也不该把同一个
  元素的赝势与轨道拆到两个系列去——数值原子轨道是为特定赝势生成的赝波函数，
  两者必须同源。

### 外部系列

`ABACUS-APNS-PPORBs-v1`、`Dojo-NC-FR` 这类**不随包分发**的系列，两种用法：

1. 把下载好的目录**整个放进 `PP-Orb/`**，然后在 INPUT / STRU 开头的库家族提示里
   选中它——会自动登记到 `~/.abacuscopilot/config.yaml` 的 `libraries.families.<id>`；
2. 或者在配置里手工登记 `libraries.families.<id>`（指向它的目录树）。

注意：登记在 `libraries.families` 里的目录（哪怕物理上就放在 `PP-Orb/` 下）会被
按**路径**排除出 SG15 内置检测；上面两个已知外部系列则是按**名字**、在登记之前就
已经排除——两道机制合起来，切回 SG15 时两个系列绝不混用。

---
More: 项目根 README 的 §Configuration，或 `abacuscopilot -task 9902` 查看当前生效的库列表。
