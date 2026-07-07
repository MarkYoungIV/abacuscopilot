# ABACUS DFT Software Knowledge for abacuskit Development

---

## ABACUS User Manual V0.2 — Formal Reference

*Extracted from: ABACUS User Manual V0.2 (2026-02-06, Peking University)*
*Source: `/Users/young/young/claude-code_file/young_pc_file/Markdown_knoledge/abacus/ABACUS_用户手册_V0p2.md`*

---

### A. INPUT File Parameters (Complete Specification)

#### A.1 General Control Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `suffix` | str | `"autotest"` | Output directory suffix: `OUT.<suffix>/` |
| `calculation` | str | `"scf"` | `scf`, `relax`, `cell-relax`, `md`, `nscf`, `get_pchg`, `get_wf` |
| `stru_file` | str | `"STRU"` | Path to structure file |
| `kpoint_file` | str | `"KPT"` | Path to k-point file |
| `pseudo_dir` | str | `"./"` | Pseudopotential file directory |
| `orbital_dir` | str | `"./"` | Numerical atomic orbital file directory |
| `read_file_dir` | str | `OUT.<suffix>` | Directory to read restart files from |
| `ntype` | int | — | Number of element types |
| `symmetry` | int | `0` | `-1`=off, `0`=time-reversal only, `1`=full symmetry |
| `init_chg` | str | `"atomic"` | Initial charge density: `"atomic"` or `"file"` |
| `init_wfc` | str | `"atomic"` | Initial wavefunction: `"atomic"`, `"nao"`, `"file"` |
| `kspacing` | float/array | — | Auto KPT mesh; 1D or 3D value. Overrides KPT file. |
| `gamma_only` | bool | `0` | Gamma-point only (LCAO fast mode, KPT ignored) |
| `device` | str | `"cpu"` | `"cpu"` or `"gpu"` |
| `precision` | str | `"double"` | `"single"`, `"double"`, `"mixing"` |

#### A.2 Basis Set Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `basis_type` | str | — | `"pw"` (plane wave), `"lcao"` (NAO), `"lcao_in_pw"` (testing) |
| `ecutwfc` | float | — | Plane wave cutoff in **Ry**. For LCAO: grid integration precision only. |
| `ks_solver` | str | PW: `cg`; LCAO: `genelpa` | PW: `cg`, `dav`, `bpcg`; LCAO CPU: `genelpa`, `scalapack_gvx`, `lapack`; LCAO GPU: `cusolver`, `cusolvermp`, `elpa` |
| `nbands` | int | auto | Number of KS orbitals |
| `nspin` | int | `1` | `1`=non-magnetic, `2`=collinear spin, `4`=non-collinear |
| `lspinorb` | bool | `0` | Spin-orbit coupling (needs fully relativistic PPs) |
| `noncolin` | bool | `0` | Non-collinear magnetism |
| `lmaxmax` | int | `2` | Max angular momentum for pseudopotential |

#### A.3 SCF Iteration & Mixing

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scf_nmax` | int | `100` | Max SCF iterations |
| `scf_thr` | float | `1e-7` | Charge density convergence. PW: `1e-8`–`1e-9`; LCAO: `1e-7` |
| `mixing_type` | str | `"broyden"` | `"broyden"`, `"pulay"` (DIIS), `"plain"` (linear) |
| `mixing_beta` | float | `0.8` (nspin=1); `0.4` (nspin=2/4) | New charge fraction |
| `mixing_ndim` | int | `8` | DIIS history length; higher=better convergence |
| `mixing_gg0` | float | `1.0` | Kerker preconditioning. `0.0`=off. |
| `mixing_gg0_min` | float | `0.1` | Min |q| for Kerker |
| `mixing_beta_mag` | float | `4*mixing_beta` (max 1.6) | Magnetic density mixing (nspin=2/4) |
| `mixing_gg0_mag` | float | `0.0` | Kerker for magnetic density (default off) |
| `mixing_angle` | float | `0.0` | Angle mixing for non-collinear. Set `1.0` to enable. |
| `mixing_tau` | bool | `false` | Mix kinetic energy density (meta-GGA only) |
| `mixing_dmr` | bool | `false` | Mix density matrix (DFT+U/EXX/DeePKS) |
| `mixing_restart` | int | `0` | Restart mixing at this SCF step |
| `chg_extrap` | str | — | Charge extrapolation: `"second-order"` for relax/MD |

#### A.4 Smearing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `smearing_method` | str | `"gauss"` | `"gauss"`, `"mp"`/`"mp2"` (Methfessel-Paxton), `"fd"` (Fermi-Dirac), `"fixed"`, `"mv"` (Marzari-Vanderbilt) |
| `smearing_sigma` | float | `0.015` Ry | Smearing width in **Ry** |

#### A.5 Relaxation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `relax_method` | str | `"bfgs"` | `"cg"`, `"bfgs"`, `"bfgs_trad"`, `"cg_bfgs"`, `"sd"`, `"fire"` |
| `relax_nmax` | int | `50` | Max ionic steps |
| `force_thr_ev` | float | `0.01` | Force threshold in **eV/Angstrom** |
| `stress_thr` | float | `1.0` | Stress threshold in **kbar** (cell-relax) |
| `cal_force` | bool | `0` | Compute forces |
| `cal_stress` | bool | `0` | Compute stress |
| `cal_syns` | bool | `0` | Compute synced forces |
| `fixed_axes` | str | `"None"` | Fix axes: `"xy"`, `"z"`, `"volume"`, `"shape"`, `"a"`, `"b"`, `"c"` |
| `fixed_atoms` | bool | — | Fix atoms (also per-atom in STRU) |
| `out_stru` | bool | `0` | Output intermediate structures |

#### A.6 DFT Functional

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dft_functional` | str | `"default"` (PBE) | `"pbe"`, `"lda"`, `"hse"`, `"pbe0"`, `"hf"`, `"scan0"`, or Libxc: `"XC_GGA_X_PBE+XC_GGA_C_PBE"` |
| `esolver_type` | str | `"ksdft"` | `"ksdft"`, `"sdft"` (stochastic), `"tddft"` (real-time TDDFT), `"ofdft"` (orbital-free), `"lj"` (Lennard-Jones), `"dp"` (deep potential) |

#### A.7 Hybrid Functional (EXX) — LCAO + LibRI

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `exx_hybrid_alpha` | float | `0.25` (1.0 for HF) | Fock exchange fraction |
| `exx_separate_loop` | bool | `true` | Double-loop (inner GGA + outer hybrid) |
| `exx_hybrid_step` | int | `100` | Max outer loop iterations |
| `exx_mixing_beta` | float | `1.0` | DM mixing beta for inner loop |
| `exx_pca_threshold` | float | `1e-4` | PCA threshold |
| `exx_cauchy_threshold` | float | `1e-7` | Cauchy-Schwarz truncation |
| `exx_ccp_rmesh_times` | float | `1.5` (HSE); `5` (others) | Coulomb cutoff radius multiplier |
| `exx_erfc_alpha` | float | varies | Screening mixing fraction |
| `exx_erfc_omega` | float | varies | Screening parameter (HSE03: 0.15, HSE06: 0.11) |

#### A.8 Molecular Dynamics Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md_type` | str | `"nvt"` | `"nve"`, `"nvt"`, `"npt"`, `"f1"`, `"msst"` |
| `md_nstep` | int | `10` | Total MD steps |
| `md_dt` | float | `1.0` | Time step in **fs** |
| `md_tfirst` | float | `300.0` | Initial temperature in **K** |
| `md_tlast` | float | — | Final temperature for ramping |
| `md_dumpfreq` | int | `1` | MD_dump output frequency |
| `md_restartfreq` | int | `1` | Restart/STRU output frequency |
| `init_vel` | bool | `0` | Read velocities from STRU |
| `dump_force` | bool | — | Include forces in MD_dump |
| `dump_vel` | bool | — | Include velocities in MD_dump |
| `dump_virial` | bool | — | Include virial in MD_dump |

#### A.9 DOS and Band Structure Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `out_dos` | int | `0` | `0`=off, `1`=total DOS, `2`=DOS+PDOS (LCAO XML) |
| `out_band` | bool | `0` | Output band structure (nscf + line-mode KPT) |
| `out_bandgap` | bool | `0` | Output band gap |
| `out_proj_band` | bool | `0` | Output projected band structure |
| `dos_emin_ev` | float | — | DOS min energy in **eV** |
| `dos_emax_ev` | float | — | DOS max energy in **eV** |
| `dos_edelta_ev` | float | `0.1` | DOS energy spacing in **eV** |
| `dos_sigma` | float | `0.07` | DOS Gaussian broadening in **eV** |
| `dos_nche` | int | `100` | Chebyshev expansion order (SDFT) |

#### A.10 Output Control Flags

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `out_chg` | bool | `0` | Total charge density (cube: `SPIN*_CHG.cube`) |
| `out_pot` | int | `0` | `1`=total local pot; `2`=electrostatic pot |
| `out_wfc_norm` | bool | `0` | Wavefunction modulus (PW, cube) |
| `out_wfc_lcao` | bool | `0` | LCAO wavefunction coefficients |
| `out_elf` | int | `0` | Electron Localization Function |
| `out_mul` | bool | `0` | Mulliken population analysis |
| `out_stru` | bool | `0` | Structure files during relax/MD |
| `out_dm` | bool | `0` | Density matrix |
| `out_dm1` | bool | `0` | First-order density matrix |
| `out_mat_hs` | bool | `0` | Hamiltonian and overlap matrices |
| `out_mat_r` | bool | `0` | r(R) matrices |
| `out_mat_t` | bool | `0` | Kinetic energy matrix |
| `out_mat_dh` | bool | `0` | dH matrix |
| `out_bandgap` | bool | `0` | Band gap info |
| `out_band` | bool | `0` | Band structure |
| `out_proj_band` | bool | `0` | Projected bands |
| `out_app_flag` | bool | `0` | Application-specific output |
| `out_interval` | int | `1` | Output interval |
| `out_alllog` | bool | `0` | All log files |
| `out_level` | str | `"ie"` | Verbosity level |
| `restart_save` | bool | `0` | Save restart files |
| `restart_load` | bool | `0` | Load restart files |
| `mem_saver` | int | `0` | Memory saving mode |

#### A.11 Band-Decomposed Charge Density

- `bands_to_print`: format `2*1 2*0 1` = first 2 bands included, next 2 skipped, next 1 included
- `if_separate_k`: separate k-point contributions
- Output: `BAND*_K*_SPIN*_CHG.cube`

#### A.12 Parallelization

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `kpar` | int | `1` | K-point parallelization groups |
| `bndpar` | int | `1` | Band parallelization groups |

---

### B. STRU File Format (Complete)

```
ATOMIC_SPECIES
<Element> <Mass> <PP_File> [pp_type]   # pp_type: upf/upf201/vwr/blps/auto

NUMERICAL_ORBITAL              # LCAO only; omit for PW
<Orbital_File>                  # One per species

LATTICE_CONSTANT
<value>                         # Scale factor (1.889726... = 1 Angstrom)

LATTICE_VECTORS
<v1x> <v1y> <v1z>
<v2x> <v2y> <v2z>
<v3x> <v3y> <v3z>

ATOMIC_POSITIONS
<coord_type>                    # Direct / Cartesian / Cartesian_au / Cartesian_angstrom
<Element>
<default_mag>                   # 0.0; per-atom mag overrides
<n_atoms>
<x> <y> <z> [m] <fx> <fy> <fz> [mag <val>] [v <vx> <vy> <vz>] [angle1 <deg>] [angle2 <deg>] [lambda <val>] [sc <val>]
...
```

**Coordinate conventions:**
- `Direct` = fractional (0-1), scaled by lattice vectors
- `Cartesian` = Bohr (multiplied by LATTICE_CONSTANT)
- `Cartesian_au` = Bohr directly
- `Cartesian_angstrom` = Angstrom directly
- Cell vectors: `cell[i] = LATTICE_CONSTANT * LATTICE_VECTORS[i]`

**Per-atom tokens:**
- `m` (optional): precedes movement flags `fx fy fz` (1=free, 0=fixed)
- `mag <val>`: collinear magnetic moment; `magmom <mx> <my> <mz>` for non-collinear (nspin=4)
- `v <vx> <vy> <vz>`: initial velocities for MD restart
- `angle1/angle2`: for constrained magnetization
- `lambda`: constraint lambda
- `sc`: self-consistent magnetization value

**Magnetism rules:**
- nspin=1: mag values ignored
- nspin=2: single scalar; if all atoms mag=0, ABACUS auto-sets 1.0
- nspin=4: three components; `mag <mx> <my> <mz>` per atom

**Empty/ghost atoms:** if element name contains `empty` suffix (e.g., `H_empty`), atom has no nuclear charge.
**MD extension:** STRU includes velocity info per atom.

---

### C. KPT File Format (Complete)

**Mode 1 -- Auto Monkhorst-Pack mesh (line 2 = `0`):**
```
K_POINTS
0
Gamma|mp
<nx> <ny> <nz> <sx> <sy> <sz>
```
- `Gamma` = Gamma-centered; `mp` = standard MP; `sx sy sz` typically `0 0 0`
- Total k-points = nx*ny*nz before symmetry

**Mode 2 -- Explicit list (line 2 = `N > 0`):**
```
K_POINTS
<N>
Direct|Cartesian
<kx> <ky> <kz> <weight>
... (N lines)
```

**Mode 3 -- Line mode for band structure (line 3 = `Line`):**
```
K_POINTS
<N_segments>
Line|Line_Cartesian
<kx> <ky> <kz> <n_pts>  #<Label>
... (N_segments endpoints, last one n_pts=1)
```
ABACUS interpolates `n_pts` k-points between consecutive endpoints.

**KPT override precedence:**
1. `gamma_only=1` → KPT ignored, Gamma only
2. `kspacing > 0` → KPT ignored, auto-generated
3. Otherwise → KPT file read as-is

---

### D. Output Files Conventions

| File | Content | Trigger |
|------|---------|---------|
| `OUT.<suffix>/running_scf.log` | SCF convergence, energy, forces, Fermi energy | Always |
| `OUT.<suffix>/istate.info` | Band/k-point occupations | Always |
| `OUT.<suffix>/kpoints` | Actual k-points (after symmetry) | Always |
| `OUT.<suffix>/STRU_ION_D` | Final structure after relax | relax/cell-relax |
| `OUT.<suffix>/SPIN*_CHG.cube` | Total charge density | `out_chg=1` |
| `OUT.<suffix>/BAND*_K*_SPIN*_CHG.cube` | Band-decomposed charge density | `bands_to_print` |
| `OUT.<suffix>/DOS1_smearing.dat` | Total DOS (eV, DOS, integrated DOS) | `out_dos=1` |
| `OUT.<suffix>/DOS2_smearing.dat` | DOS spin-down | `out_dos=1` + `nspin=2` |
| `OUT.<suffix>/PDOS` | Partial DOS XML | `out_dos=2` (LCAO) |
| `OUT.<suffix>/BANDS_1.dat` | Band structure (k-index, k-dist, energies[eV]) | `out_band=1` |
| `OUT.<suffix>/BANDS_2.dat` | Bands spin-down | `out_band=1` + `nspin=2` |
| `OUT.<suffix>/ELF.cube` | Electron Localization Function | `out_elf=1` |
| `OUT.<suffix>/STRU/STRU_MD_$istep` | MD trajectory structures | MD |
| `OUT.<suffix>/Restart_md.dat` | MD restart file | MD |
| `OUT.<suffix>/MD_dump` | MD trajectory (positions, forces, velocities, virial) | MD |
| `OUT.<suffix>/orb_matrix.0.dat` | Overlap matrix | orbital generation |
| `OUT.<suffix>/ElecStaticPot.cube` | Electrostatic potential | `out_pot=2` |

**Cube file header format:**
```
Line 1: Comment (free text)
Line 2: nspin  fermi_energy(Ry)
Line 3: natoms  origin_x  origin_y  origin_z  (Bohr)
Line 4: nx  v1_x  v1_y  v1_z
Line 5: ny  v2_x  v2_y  v2_z
Line 6: nz  v3_x  v3_y  v3_z
Line 7+: atomic_number  valence  x  y  z
Grid data: 6 per line, z-fastest (column-major)
```

---

### E. Unit Systems and Conversions

| Quantity | ABACUS Unit | VASP Unit | Conversion |
|----------|------------|-----------|------------|
| Length (default) | **Bohr** | Angstrom | 1 Bohr = 0.529177 A; 1 A = 1.889726 Bohr |
| Energy cutoff | **Ry** | eV | 1 Ry = 13.6057 eV |
| Smearing sigma | **Ry** | eV | Same as above |
| Force | **eV/A** | eV/A | Same |
| Stress | **kbar** | kbar | Same |
| MD time step | **fs** | fs | Same |
| MD temperature | **K** | K | Same |
| DOS energy range | **eV** | eV | Same |
| Charge density (cube) | **1/Bohr^3** | e/A^3 | Different |

**Key VASP migration traps:**
1. `ecutwfc=50` means 50 Ry (~680 eV), NOT 50 eV
2. ABACUS uses NCPP which requires higher ecutwfc than VASP's PAW
3. Smearing `smearing_sigma=0.015` means 0.015 Ry (~0.2 eV), NOT 0.015 eV
4. Default length unit is Bohr (not Angstrom)

---

### F. Pseudopotential & Orbital File Naming

**Pseudopotential files (UPF format):**
```
<Element>_<PP_Type>_<XC>_<Version>.upf
Example: Si_ONCV_PBE-1.0.upf
```
- PBE PPs can be used for meta-GGA (SCAN) and hybrid (PBE0, HSE06)
- LDA PPs cannot be used for SCAN/PBE0/HSE06
- SOC calculations (`lspinorb=1`) require fully relativistic PPs

**Orbital files (LCAO only):**
```
<Element>_<XC>_<RCut>au_<Ecut>Ry_<Basis>.orb
Example: Si_gga_8au_60Ry_2s2p1d.orb
- gga = functional type
- 8au = orbital cutoff radius (in Bohr)
- 60Ry = recommended ecutwfc
- 2s2p1d = DZP basis (2s + 2p + 1d)
```
- Orbitals MUST match pseudopotential (generated from it)
- ecutwfc in LCAO only improves grid integration, NOT basis completeness
- Upgrade SZ→DZP→TZDP for better accuracy (not ecutwfc)
- BSSE affects molecules more than solids in LCAO

---

### G. SCF Convergence Troubleshooting Guide

**Non-magnetic (nspin=1):**
- Decrease `mixing_beta` (0.8→0.3)
- Increase `mixing_ndim` (8→20)
- Molecules/insulators: `mixing_gg0=0.0` (Kerker off)
- Metals: adjust `mixing_gg0` and `mixing_gg0_min`

**Collinear magnetic (nspin=2):**
- Decrease `mixing_beta` and `mixing_beta_mag` proportionally
- Extreme case: `mixing_beta=0.4`, `mixing_beta_mag=0.4`, `mixing_gg0=0.0`

**Non-collinear (nspin=4):**
- Try `mixing_angle=1.0` for stubborn cases

**DFT+U:** `mixing_dmr=true`, `mixing_restart=10`
**meta-GGA:** `mixing_tau=true`

---

### H. ABACUS vs VASP Quick Reference

| Aspect | ABACUS | VASP |
|--------|--------|------|
| PP format | UPF (NCPP/USPP) | PAW (POTCAR) |
| Energy cutoff | **Ry** | eV |
| Smearing sigma | **Ry** | eV |
| Default length | **Bohr** | Angstrom |
| Basis types | PW + LCAO (NAO) | PW only |
| Orbital files | .orb (for LCAO) | N/A |
| KPT file | KPT | KPOINTS |
| Structure file | STRU | POSCAR |
| Parameter file | INPUT | INCAR |
| Output dir | OUT.suffix/ | current dir |
| PW eigensolver | cg, dav, bpcg | Davidson, RMM-DIIS |
| LCAO eigensolver | genelpa, scalapack_gvx, cusolver | N/A |
| Hybrid func | LibRI (LCAO) or native (PW) | Native |
| GPU | PW + LCAO (CUDA, DCU) | Limited |

---

### I. Tools Ecosystem

| Tool | Purpose |
|------|---------|
| **ATOMKIT** | Structure conversion, symmetry, k-path, format export (direct VASPKIT analog) |
| **ASE-ABACUS** | Python ABACUS I/O |
| **abacustest** | High-throughput INPUT generation and job management |
| **PYATB** | Berry curvature, band unfolding |
| **abacus-plot** | PDOS plotting (in ABACUS repo `tools/plot-tools/`) |

---

## Tutorial Insights

Practical knowledge extracted from three ABACUS tutorials (v3.11.0-beta.3 / LTSv3.10.0), covering Li6PS5Cl (electrolyte) and SiC (semiconductor) as example materials.

---

### 1. Structure Optimization Workflow

#### 1.1 Standard cell-relax INPUT parameters (PW basis, Li6PS5Cl)

```
calculation         cell-relax
symmetry            1
kspacing            0.14          # unit: 1/bohr; 0.14 = medium density
precision           double        # LTS version: always double; latest PW can use single
ecutwfc             80            # unit: Ry (~1088 eV); moderate accuracy
basis_type          pw
ks_solver           dav_subspace  # Davidson subspace for PW
smearing_method     gauss
smearing_sigma      0.015         # unit: Ry
mixing_type         broyden
mixing_beta         0.8           # reduce to 0.3-0.5 if SCF diverges
scf_nmax            100
scf_thr             1e-08         # PW: 1e-8 to 1e-9; LCAO: 1e-7
relax_method        cg            # or bfgs, bfgs_trad, cg_bfgs, sd, fire
relax_nmax          60
cal_force           1
force_thr_ev        0.01          # unit: eV/A
cal_stress          1
stress_thr          0.5           # unit: kbar
fixed_axes          None          # or volume, shape, a, b, c, ab, ac, bc
```

#### 1.2 Key differences for LCAO basis

| Parameter | PW value | LCAO value |
|---|---|---|
| `basis_type` | pw | lcao |
| `ks_solver` | dav_subspace | genelpa (CPU multi-core) / cusolver (single GPU) |
| `scf_thr` | 1e-08 to 1e-09 | 1e-07 |
| Orbital files | Not needed | Required (*.orb) |

Even with LCAO, `ecutwfc` is still set (e.g., 80 Ry) because the code uses plane waves for grid-related operations.

#### 1.3 Two strategies for structure optimization

1. **EOS fitting + relax (ISIF=2)**: Preserves space group symmetry but overestimates lattice parameters vs experiment.
2. **cell-relax directly (ISIF=3)**: Closer to experimental lattice, especially combined with `dft_functional pbesol`. May break symmetry for doped systems or cause ionic convergence issues.

---

### 2. MD Simulation Workflow

#### 2.1 PW-based AIMD INPUT parameters

Key changes from relaxation to MD:

| Parameter | Relaxation | MD (PW) | Reason |
|---|---|---|---|
| `calculation` | cell-relax | md | — |
| `symmetry` | 1 | 0 | Must be off for MD |
| `ecutwfc` | 80 | 60 | Reduced for speed |
| `gamma_only` | — | 1 | Only gamma point, big speedup |
| `mixing_type` | broyden | pulay | More stable for MD |
| `mixing_beta` | 0.8 | 0.3 | Conservative for SCF stability |
| `chg_extrap` | — | second-order | Charge extrapolation between MD steps |
| `scf_thr` | 1e-08 | 1e-07 | Relaxed for MD |
| `kspacing` | 0.14 | commented out | Redundant when gamma_only=1 |

Standard MD parameters:
```
md_type             nvt          # nvt, npt, nve, langevin, fire, msst
md_nstep            1000
md_dt               2.0          # unit: fs
md_tfirst           300          # unit: K
md_tlast            300          # optional; omit for constant-T, include for temperature ramp
```

#### 2.2 LCAO-based AIMD (much faster)

Same as PW MD but with `basis_type lcao`, `ks_solver genelpa`, and `scf_thr 1e-05` (further relaxed). Performance on 52-atom Li6PS5Cl: PW took 2h26min for 1000 steps, LCAO took 33min — roughly 4x faster.

#### 2.3 Thermostat options

| Keyword | Method | Notes |
|---|---|---|
| nhc | Nose-Hoover chain | **Default, recommended** |
| anderson | Anderson thermostat | — |
| berendsen | Berendsen thermostat | — |
| rescaling | Velocity rescaling method 1 | — |
| rescale_v | Velocity rescaling method 2 | Recommended by some practitioners; note this is not the CSVR thermostat from cp2k |

#### 2.4 MD-specific monitoring

- Log file: `OUT.ABACUS/running_md.log`
- `grep STEP running_md.log` — check MD steps
- `tail running_md.log` — check total time

---

### 3. KPT File Generation

#### 3.1 kspacing (automatic k-point generation)

- `kspacing 0.14` (unit: 1/bohr) — medium density, suitable for insulators/semiconductors
- For metals: use smaller values (higher k-point density)
- Mutually exclusive with explicit KPT file — use one or the other
- Used directly in INPUT as a replacement for a separate KPT file

#### 3.2 Explicit KPT file format

```
K_POINTS
0
Gamma
3 3 3 0 0 0
```

- First line: number of k-points (0 = automatic gamma-centered mesh)
- Second line: mesh type (Gamma or Monkhorst-Pack)
- Third line: nx ny nz sx sy sz (mesh subdivisions and shifts)

#### 3.3 K-point settings by calculation type

| Calculation type | Typical k-point setting |
|---|---|
| scf / relax / cell-relax | kspacing 0.14 or 3 3 3 mesh |
| md (gamma_only=1) | No KPT needed (only gamma point used) |
| DP-GEN init bulk (MD) | 3 3 3 mesh (even for MD during training) |
| DP-GEN fp (scf) | kspacing 0.14 |

#### 3.4 abacustest KPT flag

```
abacustest model inputs ... --kpt 3 3 3
```
Accepts one or three integers. A single integer like `2` is expanded to `2 2 2 0 0 0`.

---

### 4. STRU File Preparation

#### 4.1 Primary method: abacustest model inputs

```bash
# For PW (no orbital files needed):
abacustest model inputs -f POSCAR --ftype poscar --jtype cell-relax \
  --pp ~/bin/abacus/SG15-Version1p0_Pseudopotential/SG15_ONCV_v1.0_upf/

# For LCAO (needs orbital files):
abacustest model inputs -f Li6PS5Cl.cif --ftype cif --jtype cell-relax \
  --pp ~/bin/abacus/SG15-Version1p0_Pseudopotential/SG15_ONCV_v1.0_upf/ \
  --orb ~/bin/abacus/SG15-Version1p0_StandardOrbitals-Version2p0/

# Batch mode (wildcards):
abacustest model inputs -f *.cif --ftype cif ...
```

Supported `--jtype` values: `scf`, `relax`, `cell-relax`, `md`, `band`, `uscf`.

Additional flags: `--dftu`, `--dftu_param`, `--init_mag`, `--afm`, `--nspin`, `--soc`, `--lcao`.

#### 4.2 Alternative: atomkit

Interactive CLI tool (similar to vaspkit) by the same author. Menu-driven conversion from CIF to STRU.

#### 4.3 STRU file structure

Four sections:
1. **ATOMIC_SPECIES**: element, atomic mass, pseudopotential filename
2. **NUMERICAL_ORBITAL**: orbital filenames (LCAO only)
3. **LATTICE_CONSTANT**: 1.889726 (Bohr-to-Angstrom conversion factor)
4. **LATTICE_VECTORS**: 3x3 matrix
5. **ATOMIC_POSITIONS**: Direct or Cartesian; each element block: label, magnetism, number of atoms, coordinates with `1 1 1 mag 0.0` triplet

#### 4.4 Important caveats

- CIF files with full symmetry: conversion works correctly.
- CIF files reduced to P1 symmetry (e.g., after optimization): may cause errors in abacustest. **Use POSCAR format instead** for post-optimization structures.
- PW needs only *.upf files; LCAO needs both *.upf and *.orb.
- Recommended pseudopotentials: SG15 ONCV v1.0 (norm-conserving).
- Recommended orbitals: SG15 StandardOrbitals v2.0.

---

### 5. Typical Workflow Patterns

#### 5.1 Standard DFT workflow

```
1. cell-relax (structure optimization, ISIF=3)
   ↓  STRU_ION_D from OUT.ABACUS/
2. scf (self-consistent charge density, out_chg=1)
   ↓  chg.cube from OUT.ABACUS/
3. band / DOS (non-self-consistent, uscf, reads charge density)
```

For step 1, if scf convergence is difficult: set `out_chg = 1` during relaxation to save charge density for restart. For band structure: the preceding SCF must output charge density with `out_chg = 1`.

#### 5.2 Parameter changes between workflow steps

| Step | Key changes |
|---|---|
| Relax → SCF | `calculation scf`, remove `cal_force`/`cal_stress` (unless training data), keep `out_chg=1` |
| SCF → band/DOS | `calculation uscf`, read charge density from SCF |
| Relax → MD | `calculation md`, `symmetry 0`, `gamma_only 1`, reduce `ecutwfc`, change `mixing_type` to `pulay`, reduce `mixing_beta` to `0.3`, add `chg_extrap second-order` |

#### 5.3 DP-GEN workflow (machine learning potential training)

```
1. cell-relax (dft_functional pbesol recommended for experimental lattice)
   ↓  STRU_ION_D
2. dpgen init_bulk
   - stages 1,2 with skip_relax=true (if already optimized)
   - then stages 3,4 for AIMD on perturbed structures
   ↓  initial training data + perturbed configurations
3. dpgen run (iterative: 00.train → 01.model_devi → 02.fp)
   ↓  iter.* directories with training data
4. dpgen collect → consolidated training set
5. auto-test → EOS, elastic constants
```

Critical for DP-GEN INPUT files: `cal_force 1` and `cal_stress 1` must be explicitly set — ABACUS does not default these to on unlike VASP. Training needs both energies and forces. Also set `symmetry 0` for MD steps.

---

### 6. Common Pitfalls and Best Practices

#### 6.1 SCF convergence

- **First resort**: reduce `mixing_beta` from 0.8 to 0.3–0.5.
- **Second resort**: check `smearing_method` and `smearing_sigma`:
  - Metals: use `mp`, `mp2`, or `gauss`
  - Insulators/semiconductors: use `fixed` or `fd` (Fermi-Dirac)
  - Wrong smearing is a common cause of SCF non-convergence.
- For difficult cases: `out_chg = 1` to save charge density for restart.

#### 6.2 Basis set choice

- **PW**: higher accuracy, slower. Good for benchmarks and training data.
- **LCAO**: faster (4x+ for MD), good for large systems and production MD.
- ABACUS uses norm-conserving pseudopotentials, which require higher ecutwfc than VASP's PAW. Do not directly compare ecutwfc values — run convergence tests for each code/pseudopotential combination.
- GPU acceleration: use `cusolver` for single-GPU LCAO calculations.

#### 6.3 Structure preparation gotchas

- P1 symmetry CIF files from previous optimizations may cause abacustest errors. Convert to POSCAR first.
- `atomkit` and `abacustest model inputs` are the two recommended conversion tools.
- Always verify the ATOMIC_SPECIES section matches your pseudopotential filenames.

#### 6.4 MD-specific pitfalls

- **Symmetry must be off** (`symmetry 0`) for MD — otherwise results are wrong.
- **cal_force=1 and cal_stress=1** must be explicitly set for DP training data generation.
- When using `gamma_only=1`, comment out or remove `kspacing`.
- `chg_extrap second-order` is important for MD stability — extrapolates charge density between steps.

#### 6.5 Precision settings

| Context | Recommendation |
|---|---|
| LTS version (any basis) | `precision double` |
| Latest version, PW | `precision single` possible |
| Latest version, LCAO | `precision double` (but use `gint_*` precision control for speed) |

#### 6.6 Material-specific notes

- **Metals**: need denser k-points (smaller kspacing), appropriate smearing (mp/mp2/gauss).
- **Insulators/semiconductors**: medium kspacing (0.14) is typically sufficient.
- **Doped systems**: cell-relax may break symmetry or fail to converge ionic steps — consider EOS fitting + relax as alternative.

---

### 7. Post-processing Patterns

#### 7.1 Key output files

| File | Location | Purpose |
|---|---|---|
| `running_cell-relax.log` | OUT.ABACUS/ | Ionic step monitoring — `grep STEP` / `grep converged` |
| `running_md.log` | OUT.ABACUS/ | MD step monitoring — `grep STEP` |
| `STRU_ION_D` | OUT.ABACUS/ | Optimized structure (input for next calculation) |
| `STRU.cif` | OUT.ABACUS/ | Original structure in CIF format |
| `ABACUS-CHARGE-DENSITY.restart` | OUT.ABACUS/ | Restart file for interrupted calculations |
| `istate.info` | OUT.ABACUS/ | Electronic state information |
| `kpoints` | OUT.ABACUS/ | K-point information |
| `chg.cube` | OUT.ABACUS/ | Charge density (only if `out_chg=1`) |

#### 7.2 Process monitoring commands

```bash
# Check ionic steps
grep STEP running_cell-relax.log

# Check convergence status
grep converged running_cell-relax.log

# Check total runtime
tail running_cell-relax.log
```

#### 7.3 abacustest data collection

```bash
# Extract key results from multiple calculations
abacustest collectdata -p collect.json

# Collect config example (collect.json):
{
  "PARAM": [
    "natom",
    "energy",
    "scf_steps",
    {"energy_per_atom": "{energy}/{natom}"}
  ]
}

# Generate HTML report
abacustest report

# Compute vacancy formation energy
abacustest model vacancy post -j 000000
```

#### 7.4 DP-GEN post-processing

```bash
# Gather training data after sampling
dpgen collect -p param.json . ./collection

# Auto-test for EOS and elastic constants (separate auto-test input files)
```

#### 7.5 Convergence testing workflow

Use `abacustest prepare` with a param_test.json for batch convergence tests:

```json
{
  "prepare": {
    "example_template": ["example_base"],
    "mix_input": {
      "ecutwfc": [50, 60, 70, 80]
    },
    "mix_kpt": [2, [3, 3, 3], [4, 4, 4, 1, 1, 1]]
  }
}
```

This generates all parameter combinations automatically in separate subdirectories.

---

## Abacustest Insights

*Analysis of the official abacustest source code at `abacustest/`*
*Version analyzed: 2026-06-30*

---

### 1. Architecture Overview

Abacustest is structured in layers:

| Layer | Location | Purpose |
|---|---|---|
| `lib_model/` | `abacustest/lib_model/` | 21 high-level "models" (subcommands) for specific workflows: converge, relax, EOS, phonon, FD force, band, etc. |
| `lib_prepare/` | `abacustest/lib_prepare/` | Core library: INPUT generator, STRU parser/writer, KPT writer, structure conversion. Contains `abacus.py` (1800 lines, the core), `stru.py` (structure reader), `comm.py` (utilities), `input-params.json` (413 ABACUS parameters with metadata). |
| `launching/` | `abacustest/launching/` | DFlow/Bohrium integration for cloud job submission. Defines workflow models: Normal, Advanced, Expert, FDForce, FDStress, FDMagForce, Phonon, AutoABACUS. |
| `myflow/` | `abacustest/myflow/` | DFlow OP definitions for Predft, Rundft, Postdft. |

---

### 2. How abacustest Generates INPUT Files

The central INPUT generation logic lives in **`lib_model/model_013_inputs.py`** (the `inputs` subcommand) and **`lib_prepare/abacus.py`** (`WriteInput` function).

#### 2.1 JOB_TYPES Dictionary — Parameter Presets by Job Type

The `JOB_TYPES` dict in `model_013_inputs.py` defines default INPUT parameters for 6 job types:

```python
JOB_TYPES = {
    "scf", "relax", "cell-relax", "md", "band"
}
```

**Common defaults across all job types** (lines 14-59 of `model_013_inputs.py`):

| Parameter | Value | Notes |
|---|---|---|
| `symmetry` | `1` | Symmetry on by default |
| `ecutwfc` | `80` (Ry) | Can be overridden by `recommand_ecutwfc` from PP library |
| `scf_thr` | `1e-8` | SCF convergence threshold (Ry) |
| `scf_nmax` | `100` | Max SCF iterations |
| `smearing_method` | `gauss` | Gaussian smearing |
| `smearing_sigma` | `0.015` | Smearing width (Ry) |
| `mixing_type` | `broyden` | Charge mixing algorithm |
| `mixing_beta` | `0.8` | Mixing parameter (reduced to 0.4 for spin-polarized/non-collinear) |
| `basis_type` | `pw` | Plane-wave basis (default); LCAO if `--lcao` flag set |
| `ks_solver` | `dav_subspace` | PW solver; `genelpa` if LCAO |
| `pw_diag_ndim` | `2` | Davidson diagonalization subspace dimension |
| `pw_diag_nmax` | `20` | Max Davidson iterations |
| `precision` | `double` | Double precision |
| `kspacing` | `0.14` (1/bohr) | K-point spacing (not explicitly written as default — used as comment) |

**Job-type-specific defaults:**

| Job Type | `calculation` | Extra Parameters | Notes |
|---|---|---|---|
| `scf` | `scf` | `#cal_force=1`, `#cal_stress=1` (commented out) | Default calculation type |
| `relax` | `relax` | `cal_force=1`, `relax_method=cg`, `relax_nmax=60`, `force_thr_ev=0.01` | Ionic relaxation only |
| `cell-relax` | `cell-relax` | `cal_force=1`, `cal_stress=1`, `relax_method=cg`, `relax_nmax=60`, `force_thr_ev=0.01`, `stress_thr=0.5`, `fixed_axes=None` | Full cell + ion relaxation |
| `md` | `md` | `md_type=nvt`, `md_nstep=10`, `md_dt=1.0`, `md_tfirst=100`, `md_tlast=100` | Molecular dynamics (NVT default) |
| `band` | `scf` | Same as `scf` defaults | Band structure — calculation=scf for charge density, then non-SCF band path |

#### 2.2 LCAO Parameter Override

When `--lcao` is specified, the following are overridden (from `LCAO_PARAM` dict):

```python
LCAO_PARAM = {
    "basis_type": "lcao",
    "ks_solver": "genelpa",
    "ecutwfc": 100,       # Higher cutoff for LCAO (was 80)
    "scf_thr": 1e-7,      # Tighter SCF convergence (was 1e-8)
}
```

Note: `pw_diag_ndim` and `pw_diag_nmax` are removed for LCAO basis.

#### 2.3 Spin/Magnetic Configuration

| Scenario | Parameters Set | Notes |
|---|---|---|
| `nspin=1` (default) | No nspin parameter | Default no-spin |
| `nspin=2` (spin-polarized) | `nspin=2`, `mixing_beta=0.4`, `symmetry=0`, `onsite_radius=3`, `out_mul=1` | `out_mul` only valid for LCAO |
| `nspin=4` / `--soc` | `nspin=4`, `noncolin=1`, `mixing_beta=0.4`, `symmetry=-1`, `onsite_radius=3`, `out_mul=1` | Non-collinear; SOC adds `lspinorb=1` |

#### 2.4 DFT+U Configuration

When `--dftu` is specified:
- `dft_plus_u=1`, `orbital_corr` and `hubbard_u` arrays set per element
- Default U values: **4 eV for d-orbital elements**, **6 eV for f-orbital elements**
- d-orbital elements: Sc-Zn, Y-Cd, Hf-Hg, La, Ac, Th
- f-orbital elements: Ce-Lu, Pa-Lr
- `uramping=max(hubbard_u)` for PW, `mixing_dmr=1` for LCAO
- `mixing_restart=0.001`

#### 2.5 INPUT Generation Flow (`PrepInput.run()`)

```
1. gen_stru() → converts structure files to ABACUS STRU, finds PP/ORB, creates job dirs
2. For each job:
   a. generate_input(element) → copies JOB_TYPES[jobtype] defaults
   b. Apply LCAO overrides if --lcao
   c. Apply nspin/soc overrides
   d. Apply DFT+U overrides if --dftu
   e. For PW basis: remove LCAO-only params (out_mul, onsite_radius when nspin=2)
   f. update_input() → merge with user template, set ecutwfc from PP library recommendation
3. WriteInput() → write INPUT file
4. set_init_mag() → write magnetic moments to STRU
5. write_kpt() → write KPT file if --kpt specified
6. Generate setting.json for abacustest submit
```

#### 2.6 ecutwfc Auto-detection

The `gen_stru()` function reads `ecutwfc.json` from the PP library directory to get recommended cutoff values per element. The highest recommended ecutwfc across all elements in the structure is used as the effective ecutwfc (unless the user explicitly sets it in the INPUT template). This is stored in `recommand_ecutwfc` and applied in `update_input()`.

---

### 3. How abacustest Handles Different Job Types

#### 3.1 All Available Models (Subcommands)

From `lib_model/model_args.py`:

| Subcommand | Class | Description |
|---|---|---|
| `committest` | `CommitTest` | Test different ABACUS commits |
| `conv` | `ConvEcutwfc` | Convergence test |
| `aserelax` | `AseRelax` | ASE+ABACUS relax |
| `eos` | `Eos` | Equation of state |
| `phonon` | `Phonon` | Phonon calculation |
| `fdforce` | `fdforce` | Finite difference of force |
| `fdstress` | `fdstress` | Finite difference of stress |
| `comparem` | `CompareMetricsModel` | Compare metrics.json files |
| `fdmagforce` | `fdmagforce` | Finite diff of magnetic force |
| `magj` | `magj` | Magnetic exchange interactions |
| `sptest` | `SPTestModel` | FP test results |
| `band` | `BandModel` | Band structure |
| `inputs` | `InputsModel` | **INPUT file generator** (central) |
| `vasp2abacus` | `Vasp2AbacusModel` | VASP to ABACUS conversion |
| `elastic` | `ElasticModel` | Elastic constants |
| `bec` | `BECModel` | Born effective charge |
| `vacancy` | `VacancyModel` | Vacancy formation energy |
| `supercell` | `SuperCellModel` | Supercell extension |
| `vibration` | `VibrationModel` | Vibration frequency |
| `workfunc` | `WorkFuncModel` | Work function |
| `dos-pdos` | `DOSPDOSModel` | DOS/PDOS post-processing |

#### 3.2 Job Type to calculation Mapping

The `inputs` model maps job types directly to ABACUS `calculation` parameter:

- `scf` → `calculation = scf`
- `relax` → `calculation = relax` (ion movement only)
- `cell-relax` → `calculation = cell-relax` (ion + cell)
- `md` → `calculation = md` (NVT by default)
- `band` → `calculation = scf` (then band structure follows)

#### 3.3 Launching Workflow (Cloud/Bohrium)

The `launching/` module defines DFlow-based workflows for cloud execution:

- **Normal**: Simple run DFT + post DFT (`01_Normal`)
- **Expert**: Complete pre DFT + run DFT + post DFT
- **Phonon** (`03-Phonon`): Phonon calculation workflow
- **FDForce** (`04-FDForce`): Finite difference force workflow
- **FDStress** (`05-FDStress`): Finite difference stress workflow
- **FDMagForce** (`06-FDMagForce`): Finite difference magnetic force workflow
- **AutoRun** (`07-AutoRun`): Auto-run ABACUS jobs
- **Vasp2Abacus** (`08-Vasp2Abacus`): Convert VASP jobs to ABACUS
- **Report** (`02-Report`): Report generation
- **Reuse** (`00-Reuse`): Reuse model for other datasets

---

### 4. How abacustest Handles STRU Files

#### 4.1 `AbacusStru` Class (`lib_prepare/abacus.py`, line 105)

The `AbacusStru` class is a comprehensive STRU representation with attributes:

| Attribute | Type | Purpose |
|---|---|---|
| `label` | `List[str]` | Element labels (per type or per atom) |
| `cell` | `List[List[float]]` | Lattice vectors (in Bohr * lattice_constant) |
| `coord` | `List[List[float]]` | Atomic coordinates |
| `lattice_constant` | `float` | Scaling factor (default 1) |
| `pp` | `List[str]` | Pseudopotential file names per type |
| `orb` | `List[str]` | Orbital file names per type |
| `paw` | `List[str]` | PAW file names per type |
| `mass` | `List[float]` | Atomic masses (auto from periodic table) |
| `element` | `List[str]` | Element symbols |
| `move` | `List[List[int]]` | Atom movement flags (1=move, 0=fixed) |
| `magmom` | `List[float]` | Magnetic moment per type |
| `magmom_atom` | `List[Union[float,List[float]]]` | Magnetic moment per atom (3-elem for non-collinear) |
| `velocity` | `List[List[float]]` | Atomic velocities (for MD restart) |
| `angle1/angle2` | `List[float]` | Magnetic moment angles |
| `constrain` | `List[Union[bool,List[bool]]]` | Magnetization constraints (sc) |
| `dpks` | `str` | DeepKS descriptor filename |
| `cartesian` | `bool` | Whether coords are Cartesian |

#### 4.2 Key STRU Operations

- **Read**: `AbacusStru.ReadStru(stru_file)` — parses STRU into object
- **Write**: `stru.write("STRU")` — writes STRU file with all metadata (mag, sc, v, angle1/2, lambda, dpks)
- **From dpdata**: `AbacusStru.FromDpdata(input_file, fmt)` — converts from POSCAR, CIF, etc.
- **Supercell**: `stru.supercell([na, nb, nc])` — generates supercell
- **Delete atom**: `stru.delete_atom(label, idx)` — removes specific atoms (for vacancy calculations)
- **Set empty atom**: `stru.set_empty_atom(label, idx)` — marks atoms as empty (for BSSE)
- **Perturb**: `stru.perturb_stru(n, cell_pert_frac, atom_pert_dist, mag_*)` — random perturbations
- **K-line generation**: `stru.get_kline()` / `stru.get_kline_ase()` — seekpath-based band path
- **Convert**: `stru.to_ase()`, `stru.to_pymatgen()`, `stru.to_phonopy()` — export to other frameworks
- **Coordinate transforms**: `set_coord()`, `set_cell()` with direct/Cartesian/bohr/Angstrom options
- **Magnetic setup**: `set_atommag()`, `set_atommag_from_mulliken()` — set per-atom magnetic moments

#### 4.3 Structure Conversion Pipeline

The `translate_strus()` function in `lib_prepare/comm.py` supports:
- **dpdata-supported formats**: POSCAR, ABACUS/STRU, VASP, QE, CP2K, LAMMPS, etc.
- **CIF format**: via ASE (`ase.io.read`) or pymatgen (`Structure.from_file`)
- Multi-frame structures are split into separate job directories

#### 4.4 Pseudopotential Collection

`collect_pp(pp_path)` scans a directory for files matching `Element_*.*` pattern, building a `{element: full_path}` dict. The PP/ORB files are symlinked (or copied if `--copy-pp-orb`) to each job directory. An `ecutwfc.json` file in the PP library provides per-element recommended cutoffs.

---

### 5. How abacustest Handles KPT Files

#### 5.1 KPT Generation (`WriteKpt` in `lib_prepare/abacus.py`, line 1567)

Supports **5 KPT modes**:

| Mode | ABACUS KPT Type | Input Format | Example |
|---|---|---|---|
| `gamma` | Gamma-centered Monkhorst-Pack | `[nx, ny, nz, shift_x, shift_y, shift_z]` | `[2,2,2,0,0,0]` |
| `mp` | Standard Monkhorst-Pack | Same as gamma | `[4,4,4,0,0,0]` |
| `direct` | Explicit direct-coordinate k-points | `[[kx,ky,kz,weight], ...]` | `[[0,0,0,1.0]]` |
| `cartesian` | Explicit Cartesian k-points | Same as direct | `[[0,0,0,1.0]]` |
| `line` | Line-mode (band structure) | `[[kx,ky,kz,npoints,"#comment"], ...]` | `[[0,0,0,10,"#Gamma"],[0.5,0,0,1]]` |

For gamma/mp:
- 3-element lists auto-expand to 6 with `[0,0,0]` shifts
- Single integer `n` expands to `[n,n,n,0,0,0]`

For direct/cartesian:
- 3-element sublists get equal weight `1/n`
- 4-element sublists use specified weights (normalized)

#### 5.2 KPT Auto-detection

`ReadKpt()` can read KPT from:
1. An ABACUS input directory (checks INPUT for `kpoint_file`, `kspacing`, `gamma_only`)
2. A KPT file directly
3. If `kspacing` is set and no KPT file exists, it computes k-points from the cell using `kspacing2kpt()`

The `inputs` model passes `--kpt` as `[int, int, int]` to `WriteKpt` in gamma mode.

---

### 6. Key ABACUS Parameters Exposed to Users

The `inputs` model CLI exposes these ABACUS parameters directly:

| CLI Flag | ABACUS Parameter(s) | Values |
|---|---|---|
| `--jtype` | `calculation` | `scf`, `relax`, `cell-relax`, `md`, `band` |
| `--lcao` | `basis_type=lcao`, `ks_solver=genelpa` | Flag |
| `--nspin` | `nspin` | `1`, `2`, `4` |
| `--soc` | `lspinorb=1`, `nspin=4` | Flag |
| `--dftu` | `dft_plus_u=1`, `orbital_corr`, `hubbard_u` | Flag (+ optional `--dftu_param`) |
| `--init_mag` | `mag` per atom in STRU | Element-moment pairs |
| `--afm` | Negates half of `init_mag` values | Flag |
| `--kpt` | KPT file (gamma mode) | 1 or 3 integers |
| `--input` | Template INPUT file (overrides defaults) | Path |
| `--pp` | `pseudo_dir` + STRU PP entries | Path |
| `--orb` | `orbital_dir` + STRU ORB entries | Path |
| `--copy-pp-orb` | Copy vs symlink PP/ORB files | Flag |

The `prepare.PrepareAbacus` class additionally exposes for programmatic use:
- `mix_input`: Parameter sweeping dict (`{"ecutwfc":[50,60,70], "mixing_beta":[0.3,0.4]}` → 9 combinations)
- `mix_kpt`: Multiple k-point settings
- `mix_stru`: Multiple structure files
- `pert_stru`: Structure perturbation (cell, atom positions, magnetic moments)

---

### 7. Parameter Presets Summary

#### 7.1 Default Preset (PW, no spin)

```
ecutwfc=80, scf_thr=1e-8, scf_nmax=100, mixing_beta=0.8, mixing_type=broyden
smearing_method=gauss, smearing_sigma=0.015, basis_type=pw, ks_solver=dav_subspace
pw_diag_ndim=2, pw_diag_nmax=20, precision=double, symmetry=1
calculation=scf, kspacing≈0.14
```

#### 7.2 LCAO Preset

```
Same as default PW but:
  basis_type=lcao, ks_solver=genelpa, ecutwfc=100, scf_thr=1e-7
  (pw_diag_ndim, pw_diag_nmax removed)
```

#### 7.3 Spin-Polarized Preset (nspin=2)

```
nspin=2, mixing_beta=0.4, symmetry=0
onsite_radius=3, out_mul=1 (LCAO only)
```

#### 7.4 Non-Collinear Preset (nspin=4 / SOC)

```
nspin=4, noncolin=1, mixing_beta=0.4, symmetry=-1
onsite_radius=3, out_mul=1
(+ lspinorb=1 if SOC)
```

#### 7.5 DFT+U Preset

```
dft_plus_u=1, orbital_corr=[per element], hubbard_u=[per element]
Default U: 4 eV (d-orbital), 6 eV (f-orbital)
uramping=max(U) (PW) or mixing_dmr=1 (LCAO)
mixing_restart=0.001
```

#### 7.6 Relax/Cell-Relax Preset

```
Same as default PW +:
  cal_force=1, relax_method=cg, relax_nmax=60, force_thr_ev=0.01
  (+ cal_stress=1, stress_thr=0.5 for cell-relax)
```

#### 7.7 MD Preset

```
Same as default PW +:
  md_type=nvt, md_nstep=10, md_dt=1.0, md_tfirst=100, md_tlast=100
```

---

### 8. INPUT Parameter Categorization

The `input-params.json` file (413 parameters in 27 categories) provides metadata for organized INPUT writing. When `WriteInput(categorized=True)`, parameters are grouped with comment headers:

| Category | Count | Example Parameters |
|---|---|---|
| System variables | 19 | `suffix`, `calculation`, `nspin`, `symmetry` |
| Plane wave related | 14 | `ecutwfc`, `ecutrho`, `pw_diag_thr` |
| Electronic structure | 41 | `basis_type`, `ks_solver`, `nbands`, `scf_thr`, `mixing_beta` |
| Geometry relaxation | 21 | `relax_method`, `relax_nmax`, `cal_force`, `force_thr_ev` |
| Molecular dynamics | 39 | `md_type`, `md_nstep`, `md_dt`, `md_thermostat` |
| DFT+U | 8 | `dft_plus_u`, `orbital_corr`, `hubbard_u`, `uramping` |
| vdW correction | 17 | `vdw_method`, `vdw_s6`, `vdw_C6_file` |
| Exact Exchange | 26 | `exx_hybrid_type`, various EXX parameters |
| TDDFT | 37 | Various time-dependent DFT parameters |
| Output | 36 | `out_chg`, `out_mul`, various output controls |
| PEXSI | 23 | PEXSI solver parameters |
| Input files | 7 | `stru_file`, `kpoint_file`, `pseudo_dir`, `orbital_dir` |
| DeePKS | 13 | `deepks_scf`, `deepks_model`, etc. |
| Others | OFDFT(16), E-field(6), Gate(6), Berry(11), Debug(8), Conductivity(9), Solvation(5), QO(5), LR-TDDFT(13), RDMFT(2) |

---

### 9. Key Design Patterns

#### 9.1 Parameter Sweeping (`mix_input`)

The `PrepareAbacus` class supports combinatorial parameter sweeping:

```python
mix_input = {
    "ecutwfc": [50, 60, 70],      # 3 values
    "mixing_beta": [0.3, 0.4]      # 2 values
}
# Produces 3 x 2 = 6 INPUT files
```

Combined parameters can also use the `|` syntax: `"ecutwfc|kspacing": "50|0.2"` sets both simultaneously (no combinatorial explosion).

#### 9.2 Structure Perturbation

`pert_stru` dictionary supports:
- `pert_number`: Number of perturbed copies
- `cell_pert_frac`: Random cell perturbation (diagonal up to frac, off-diagonal up to frac/2)
- `atom_pert_dist`: Random atom displacement
- `mag_rotate_angle`: Random rotation of constrained magnetization vectors
- `mag_tilt_angle`: Random tilt of constrained magnetization
- `mag_norm_dist`: Random norm change of constrained magnetization

#### 9.3 Example Template System

The `PrepareAbacus` class can start from an "example template" directory containing INPUT/STRU/KPT files. It then:
1. Reads these as base templates
2. Applies parameter overrides from `mix_input`
3. Creates parameter combinations in numbered subdirectories
4. Links/copies PP, ORB, and extra files

#### 9.4 Cross-Format Support

Abacustest can convert ABACUS inputs to:
- **Quantum ESPRESSO**: via `Abacus2Qe` class
- **VASP**: via `Abacus2Vasp` class (needs POTCAR)
- **CP2K**: via `Abacus2Cp2k` class

#### 9.5 Recommended Defaults (from `constant.py`)

```python
RECOMMAND_IMAGE = "registry.dp.tech/dptech/abacus:LTSv3.10.1"
RECOMMAND_COMMAND = "OMP_NUM_THREADS=1 mpirun -np 16 abacus | tee out.log"
RECOMMAND_MACHINE = "c32_m64_cpu"  # Bohrium machine type
```

The `inputs` model auto-generates a `setting.json` and `run.sh` for direct cloud submission via `abacustest submit`.

---

## User Guide Insights

*Extracted from the ABACUS user guide directory (~128 files, 28 example subdirectories, sampled via 5 parallel agents).*
*Source: `/Users/young/young/claude-code_file/young_pc_file/Markdown_knoledge/abacus/abacus-user-guide/`*

---

### 1. STRU File Format (Complete Specification)

Five sections in order:

```
ATOMIC_SPECIES
<Element> <mass> <pp_file> [pp_type]    # pp_type: upf/upf201/vwr/blps/auto

NUMERICAL_ORBITAL                       # LCAO only; omit entire block for PW
<orb_file>                              # one per species, same order as ATOMIC_SPECIES

LATTICE_CONSTANT
<value>                                 # scaling factor in Bohr (1.889726125 = 1 Angstrom)

LATTICE_VECTORS                         # omit if latname used in INPUT
<v1x> <v1y> <v1z>
<v2x> <v2y> <v2z>
<v3x> <v3y> <v3z>

ATOMIC_POSITIONS
<coord_type>                            # Direct / Cartesian / Cartesian_au / Cartesian_angstrom
<Element>
<default_mag>                           # 0.0 for non-magnetic; per-atom mag overrides
<n_atoms>
<x> <y> <z> [m fx fy fz] [mag val] [v vx vy vz] [angle1 deg] [angle2 deg] [lambda val] [sc val]
```

**Coordinate conventions:**
- `Direct` = fractional (0-1), scaled by lattice vectors
- `Cartesian` = multiples of LATTICE_CONSTANT (Bohr)
- `Cartesian_au` = Bohr directly
- `Cartesian_angstrom` = Angstrom directly
- Cell vectors: `cell[i] = LATTICE_CONSTANT * LATTICE_VECTORS[i]`

**Per-atom movement flags** (`m fx fy fz` or plain `fx fy fz`): 1=free, 0=fixed in x/y/z directions. Used during relaxation to freeze selected atoms (e.g., bottom layers of slab).

**Magnetic moment** (`mag val`): collinear; `magmom mx my mz` for non-collinear. Per-atom values override per-element default. If `nspin=2` and no atom has finite magnetization, ABACUS auto-sets `1.0` for every atom.

**Initial velocities** (`v vx vy vz`): for MD restarts.

**Example STRU (magnetic Fe, with selective dynamics):**
```
ATOMIC_SPECIES
Fe 55.85 Fe_ONCV_PBE-1.0.upf upf201

LATTICE_CONSTANT
1.88972612584

LATTICE_VECTORS
    2.8997376650     0.0000000000     0.0000000000
    0.0000000000     2.8997376650     0.0000000000
    0.0000000000     0.0000000000     2.8997376650

ATOMIC_POSITIONS
Direct

Fe
4                              # default magnetism
2                              # number of atoms
    0.0000000000     0.0000000000     0.0000000000 m  1  1  1 mag 4.0
    0.5000000000     0.5000000000     0.5000000000 m  1  1  1 mag 4.0
```
Note: `m` keyword is optional -- the three movement flags can appear directly after coordinates without `m`.

---

### 2. KPT File Format (Full Specification)

**Mode 1 -- Auto Monkhorst-Pack mesh (line 2 = `0`):**
```
K_POINTS
0
Gamma                   # or "MP"
nx ny nz sx sy sz
```
- `Gamma` = Gamma-centered (insulators/semiconductors); `MP` = standard (metals)
- `sx sy sz` = shifts, typically `0 0 0`; `0.5 0.5 0.5` for metallic shifted grids
- Total k-points = nx * ny * nz before symmetry reduction

**Mode 2 -- Explicit k-point list (line 2 = `N > 0`):**
```
K_POINTS
N
Direct                  # or "Cartesian" (1/Bohr)
kx ky kz weight          # N lines
```

**Mode 3 -- Line mode for band structure (line 2 = `N`, line 3 = `Line`):**
```
K_POINTS
N                       # number of high-symmetry endpoints
Line                    # or "Line_Cartesian"
kx ky kz n_pts   #Label # N endpoints, last one with n_pts=1
```
ABACUS interpolates `n_pts` k-points between consecutive endpoints.

**KPT override precedence (INPUT parameters):**
| Condition | Effect |
|---|---|
| `gamma_only = 1` | KPT file IGNORED; Gamma point only (LCAO only) |
| `kspacing > 0.0` | KPT file IGNORED; auto-generate from spacing (1/Bohr) |
| Otherwise | KPT file READ as-is |

`kspacing` can be a single value (isotropic) or three values (anisotropic, for 2D materials: `kspacing 0.1 0.1 1.0`).

**Standard high-symmetry paths** (fractional reciprocal coordinates):
| Lattice | Path |
|---|---|
| SC | G(0,0,0) - X(0.5,0,0) - M(0.5,0.5,0) - R(0.5,0.5,0.5) - G |
| FCC | G(0,0,0) - X(0.5,0,0.5) - W(0.5,0.25,0.75) - K(0.375,0.375,0.75) - G - L(0.5,0.5,0.5) - U(0.625,0.25,0.625) - W |
| BCC | G(0,0,0) - H(0.5,-0.5,0.5) - N(0,0,0.5) - P(0.25,0.25,0.25) - G - N |
| HEX | G - M(0.5,0,0) - K(1/3,1/3,0) - G - A(0,0,0.5) - L(0.5,0,0.5) - H(1/3,1/3,0.5) - A |

---

### 3. CUBE File Format (Charge Density, ELF, Wavefunction)

Common output format for volumetric data.

```
Line 1:    Comment (free text)
Line 2:    nspin  fermi_energy(Ry)
Line 3:    natoms  origin_x  origin_y  origin_z  (Bohr)
Line 4:    nx  v1_x  v1_y  v1_z            (grid step vector 1, Bohr)
Line 5:    ny  v2_x  v2_y  v2_z            (grid step vector 2, Bohr)
Line 6:    nz  v3_x  v3_y  v3_z            (grid step vector 3, Bohr)
Lines 7..(7+natoms-1): atomic_number  valence  x  y  z  (Bohr)
Remaining:  nx*ny*nz float values, 6 per row, z-fastest order (Fortran/column-major)
```

**Grid index mapping:** `data[ix][iy][iz] = values[ix*ny*nz + iy*nz + iz]` (0-indexed).  
**Position:** `R = origin + ix*v1 + iy*v2 + iz*v3`

**Unit conventions per file type:**
- `SPIN*_CHG.cube` (charge density): electrons/Bohr^3
- `elf.cube` / `ELF.cube`: dimensionless [0, 1] (ELF=0.5 = uniform electron gas)
- `ElecStaticPot.cube`: electrostatic potential (Hartree + external + dipole corrections)
- `wfc_realspace/*.cube`: wavefunction modulus

**Tools:** `cube_manipulator.py` (ABACUS source `tools/plot-tools/`) for scaling, adding, subtracting, slicing, and 1D-integration of cube files. Used for Bader charge analysis: sum `SPIN1_CHG.cube` + `SPIN2_CHG.cube` before running `bader.x`.

---

### 4. DOS and PDOS Output Format

**DOS (out_dos = 1):**
- `DOS1` -- unsmeared total DOS for spin 1 (no header, raw data)
- `DOS1_smearing.dat` -- columns: `Energy(eV)  DOS(states/eV)  Integrated_DOS`
- `DOS2` / `DOS2_smearing.dat` for nspin=2

Fermi energy extracted via: `grep EFERMI running_scf.log` (eV).

**PDOS (out_dos = 2, note: `2`, not `1`):**
- Single XML file `PDOS` (no extension)

```xml
<pdos>
  <nspin>1</nspin>
  <norbitals>720</norbitals>
  <energy_values units="eV">...</energy_values>
  <orbital index="1" atom_index="1" species="Si" l="0" m="0" z="1">
    <data>val_spin1  val_spin2  ...</data>
  </orbital>
</pdos>
```

Per-orbital attributes: `atom_index` (1-based), `species`, `l` (0=s,1=p,2=d), `m`, `z` (radial zeta index). Data rows: 1 column for nspin=1, 2 columns for nspin=2.

Post-processing tool: `tools/plot-tools/` reads PDOS + config.json to produce grouped `.dat` files by species, atom_index, l-channel, or individual orbital index.

---

### 5. Standard Output Log (`running_scf.log`)

**Key grep patterns:**
- `grep EFERMI` -- Fermi energy (eV)
- `grep "ETOT"` or `grep "CG"` -- Total energy per SCF iteration (eV) with `EDIFF` and `DRHO`
- `grep "TOTAL-FORCE"` -- Forces (eV/Angstrom)
- `grep "converged"` -- SCF convergence status
- `grep STEP` -- Relax/MD ionic step marker
- `tail running_scf.log` -- final timing summary

**Energy decomposition** (final SCF step):
`E_KohnSham`, `E_Harris`, `E_band`, `E_one_elec`, `E_Hartree`, `E_xc`, `E_Ewald`, `E_entropy(-TS)`, `E_descf`, `E_localpp`, `E_exx`, `E_Fermi`, `E_gap(k)`. Reported in both Ry and eV.

---

### 6. Unit Reference and VASP Migration Traps

| Quantity | ABACUS | VASP |
|---|---|---|
| Length (structure) | Bohr (default) or Angstrom | Angstrom |
| Energy cutoff | **Ry** | **eV** (ENCUT) |
| Smearing | **Ry** | **eV** (SIGMA) |
| K-point spacing | 1/Bohr | 1/Angstrom |
| Pseudopotential | **NCPP** (norm-conserving) | **PAW** (+ USPP) |
| Mixing | `mixing_beta` (0-1) | `AMIX`, `BMIX`, `AMIN` |

**Conversions:** 1 Ry = 13.605703976 eV; 1 Angstrom = 1.889726125 Bohr.

**Key traps:**
1. `ecutwfc=50` means 50 Ry (~680 eV in VASP) -- NOT 50 eV
2. ABACUS uses NCPP which requires higher ecutwfc than VASP's PAW. Do not compare ecutwfc values directly.
3. Smearing parameters are also in Ry (not eV)

---

### 7. LCAO vs PW: Critical Distinctions

| Aspect | PW | LCAO |
|---|---|---|
| Basis completeness | Increase `ecutwfc` | Upgrade orbital set: SZ -> DZP -> TZDP |
| `ecutwfc` function | Controls basis completeness | Controls grid integration precision only |
| k-parallelization (`kpar`) | Yes | No |
| Best for | <50 atoms, high precision | >50 atoms, large-scale MD |
| GPU support | Yes | Yes (single-GPU via cusolver) |
| Required files | INPUT, STRU, KPT, PP | + `.orb` per element |
| `gamma_only` | No | Yes |

**Orbital file naming convention:** `<Element>_<functional>_<R>au_<E>Ry_<shells>.orb`
- `Si_gga_8au_60Ry_2s2p1d.orb` = DZP basis (2s+2p+1d), 8 Bohr cutoff, 60 Ry recommended ecutwfc
- `Si_gga_8au_100Ry_2s2p1d.orb` = DZP, 8 Bohr, 100 Ry recommended ecutwfc
- Larger cutoff = more accurate but slower; 7-10 au typical

**For LCAO, the recommended ecutwfc printed on the orbital filename should be used.** No ecutwfc convergence test is needed for LCAO -- only PW needs ecut convergence testing.

---

### 8. 10-Test-Case Benchmark (ABACUS v3.9.0.19)

Concrete parameter values from production tests (all PW basis, `mixing_type=broyden`):

| # | System | Atoms | ecutwfc(Ry) | mix_beta | KPT | nbands | SCF | Energy (eV/atom) |
|---|---|---|---|---|---|---|---|---|
| 1 | GaAs | 8 | 60 | 0.7 | 7x7x7 | 46 | 8 | -979.5896 |
| 2 | C2H6O | 9 | 90 | 0.7 | 1x1x1 | 20 | 17 | -74.5619 |
| 3 | MoS2 (2D) | 12 | 60 | 0.7 | 1x7x7 | 62 | 18 | -808.7695 |
| 4 | Pt(111) | 12 | 50 | 0.3 | 1x7x7 | 129 | 19 | -3300.5524 |
| 5 | BaTiO3 | 15 | 90 | 0.3 | 7x7x6 | 72 | 12 | -716.7728 |
| 6 | Na (metal) | 16 | 60 | 0.3 | 7x7x7 | 86 | 11 | -1157.9277 |
| 7 | Fe (mag) | 27 | 60 | 0.3 | 6x6x6 | 286 | 35 | -3220.5415 |
| 8 | 32 H2O | 96 | 90 | 0.7 | Gamma | 128 | 14 | -155.7804 |
| 9 | Battery | 108 | 80 | 0.3 | Gamma | 490 | 23 | -1148.9175 |
| 10 | Si (large) | 216 | 60 | 0.7 | Gamma | 518 | 9 | -107.2307 |

**Key takeaways:**
- `mixing_beta`: 0.7 for insulators/semiconductors; 0.3 for metals/magnetic
- `ecutwfc`: 50-90 Ry range; 60 Ry sufficient for most materials
- Anisotropic KPT for 2D (1x7x7) and surfaces
- Gamma-only sufficient for >50-atom supercells
- `nbands` roughly 1.5-2x occupied for metals; tighter for insulators

---

### 9. Cross-Workflow Parameter Matrix

| Parameter | SCF(PW) | AIMD | Phonopy | Elastic | Slab |
|---|---|---|---|---|---|
| ecutwfc (Ry) | 60-100 | 30 | 100 | 100 | 60 |
| scf_thr | 1e-8 | 1e-5 | 1e-7 | 1e-7 | 1e-8 |
| scf_nmax | 50-100 | 100 | 50 | 50 | 100 |
| smearing | gauss | gaussian | mp | gaussian | gauss |
| sigma (Ry) | 0.002 | 0.001 | 0.015 | 0.002 | 0.007 |
| mixing_type | broyden | pulay | pulay | broyden | pulay |
| mixing_beta | 0.7 | 0.3 | 0.7 | 0.7 | 0.7 |
| mixing_gg0 | 0 | -- | 1.5 | -- | -- |
| ks_solver | cg | genelpa | genelpa | genelpa | cg |
| gamma_only | 0 | 1 | 0 | 0 | 0 |
| symmetry | 1 | 0 | 1 | 1 | 1 |
| basis_type | pw | lcao | lcao | lcao | pw |

---

### 10. Workflow Chain Patterns

**Standard SCF->DOS/Band:**
```
SCF (out_chg=1) -> NSCF for DOS (init_chg=file, out_dos=1) -> NSCF for Band (init_chg=file, out_band=1)
```

**Geometry Relaxation:**
```
cell-relax (cal_force=1, cal_stress=1, force_thr_ev=0.01, stress_thr=0.5)
  -> STRU_ION_D (final structure)
  -> Final SCF with tighter convergence (out_chg=1 for properties)
```

**AIMD:**
```
calculation=md, symmetry=0, basis_type=lcao, gamma_only=1
ecutwfc=30, scf_thr=1e-5, mixing_type=pulay, mixing_beta=0.3
chg_extrap=second-order, md_type=nvt, md_dt=1.0-2.0 fs
```

**Phonon (Phonopy):**
```
relax (tight force_thr_ev ~ 1e-4)
  -> phonopy -d --dim="2 2 2" -> STRU-001, STRU-002, ...
  -> ABACUS SCF per displacement (cal_force=1)
  -> phonopy -f -> FORCE_SET
  -> phonopy -p band.conf -> phonon dispersion
```

**Elastic Constants:**
```
cell-relax -> gene_dfm.py (24 strain configs) -> ABACUS relax per config -> compute_dfm.py
Strain magnitudes: -0.010, -0.005, 0.005, 0.010
```

**Surface/Work Function:**
```
Bulk relax -> Slab calculation (out_pot=2, KPT: 20x20x1)
  -> aveElecStatPot.py -> ElecStaticPot_AVE
  -> Work function = V_vacuum - E_Fermi
```

**Wannier90:**
```
wannier90.x -pp seedname -> SCF(out_chg=1) -> NSCF(towannier90=1, init_chg=file) -> wannier90.x seedname
```

---

### 11. Abacustest `inputs` Model Comparison

The abacustest `inputs` model's parameter presets (from `JOB_TYPES` / `LCAO_PARAM`) differ slightly from the user guide examples:

| Parameter | abacustest PW default | User guide PW example (GaAs) | User guide LCAO example |
|---|---|---|---|
| ecutwfc | 80 (or from PP lib) | 60 | 60 (LCAO: 100 via `--lcao`) |
| ks_solver | dav_subspace | cg | scalapack_gvx |
| mixing_beta | 0.8 | 0.7 | 0.7 |
| smearing_sigma | 0.015 | 0.002 | 0.002 |
| scf_thr | 1e-8 | 1e-7 | 1e-7 |
| precision | double | (not specified) | (not specified) |

The user guide examples use somewhat tighter smearing (0.002 Ry vs 0.015 Ry) and lower mixing_beta (0.7 vs 0.8) than abacustest defaults, reflecting that the examples are tuned for specific materials while abacustest uses safer/broader defaults.

---

### 12. Source File Paths

| Resource | Path |
|---|---|
| User guide root | `/Users/young/young/claude-code_file/young_pc_file/Markdown_knoledge/abacus/abacus-user-guide/` |
| Online docs | `https://abacus.deepmodeling.com/en/latest/` |
| INPUT reference | `https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html` |
| KPT reference | `https://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html` |
| STRU reference | `https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html` |
| PP downloads | `https://abacus.ustc.edu.cn/pseudo/list.htm` |
| Main GitHub | `https://github.com/deepmodeling/abacus-develop` |
| User guide GitHub | `https://github.com/MCresearch/abacus-user-guide` |
