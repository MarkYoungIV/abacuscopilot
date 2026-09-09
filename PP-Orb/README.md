# PP-Orb/ — 赝势 / 数值轨道库目录
## Pseudopotential & numerical-orbital library directory

这里存放 abacuscopilot 在生成 INPUT / 自动同步 STRU 时按元素查找并拷贝的
**赝势（PP，`*.upf`）** 与 **数值原子轨道（NAO，`*.orb`）** 库树。

- 发布版的 release tarball 会随包内置一组默认库（SG15 赝势、SG15 标准轨道、
  镧系 f-core 包），开箱即用。
- git 仓库本身**不**携带这些大体积库，只保留本说明文件。

### 把「你自己的」库放进这里（推荐做法）

把你手头的库目录（只要里面含 `*.upf` / `*.orb` 文件）**整个复制为 `PP-Orb/`
下的一个顶层文件夹**即可，abacuscopilot 下次运行会自动检测它——无需改配置：

```
PP-Orb/
├── SG15-Version1p0_Pseudopotential/             # 自动检测（赝势）
├── SG15-Version1p0__StandardOrbitals-Version2p0 # 自动检测（轨道）
├── lanthanides-f--core.icmod1/                  # 同时含 .upf 与 .orb
└── <你自己的系列目录>/                          # 例如另一套 SG15 / PBE-dojo …
```

规则：

- 顶层目录被递归搜索，元素→文件匹配**大小写不敏感**，容忍价态前缀
  （`Sm3+_f--core-icmod1.PD04.PBE.UPF`）与半芯 `-sp` 命名
  （`Hf-sp.PD04.PBE.UPF` / `Os-sp.PD04.PBE.UPF`）。
- 一个元素对应多个轨道时，按当前**库家族（family）**取确定性默认档
  （SG15：DZP `4s2p2d1f`@7au；APNS：`efficiency` 更小 rcut、`precision` 最完备基组）。
- 若 `~/.abacuscopilot/config.yaml` 的 `libraries.pseudo_library` /
  `orbital_library` 里已登记了旧目录，新增的顶层库不会自动并入——在 INPUT
  流程里重新选一次 **SG15** 家族（或手工在配置里登记新目录）即可触发重检。

### 外部系列（不放进 PP-Orb）

`ABACUS-APNS-PPORBs-v1` 这类**不随包分发**的系列，建议把下载好的目录登记到
`~/.abacuscopilot/config.yaml` 的 `libraries.families.<id>`（或在 `abacuscopilot
-task 9901` 向导 / INPUT 开头的家族提示里现场填写路径），再在生成时选择该系列。

注意：凡是登记在 `libraries.families` 里的目录（哪怕物理上放在 `PP-Orb/` 下），
都会被自动从 SG15 内置检测中**排除**，确保切回 SG15 时两个系列绝不混用。

---
More: 项目根 README 的 §Configuration，或 `abacuscopilot -task 9901`。
