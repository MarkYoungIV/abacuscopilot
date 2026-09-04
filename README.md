# AbacusCopilot

> ୧(๑•̀⌄•́๑)૭   A pre- & post-processing copilot for the ABACUS DFT software

**AbacusCopilot** is an interactive command-line copilot for the ABACUS DFT software, streamlining the whole workflow — generating input files, analyzing calculation outputs, and producing publication-quality figures. It provides both an interactive menu-driven interface and a command-line task mode.

## Quick Start

```bash
# Install
pip install -e .

# Launch interactive mode
abacuscopilot

# Run a task non-interactively
abacuscopilot -task 101                    # Generate SCF INPUT
abacuscopilot -task 201                    # Convert CIF → STRU
abacuscopilot -task 301                    # Generate KPT from k-spacing
abacuscopilot -task 801 --bands OUT.Si/BANDS_1.dat  # Plot band structure
abacuscopilot -task 701 --log OUT.ABACUS/running_scf.log  # SCF convergence
abacuscopilot --list-tasks                 # Show all available tasks
```

### First-time Setup

```bash
abacuscopilot -task 9901    # Configuration wizard
```

Configures pseudopotential paths, orbital directories, ABACUS binary location, and default calculation parameters. Settings are saved to `~/.abacuscopilot/config.yaml`.

> **Pseudopotential / orbital libraries (PP-Orb/)** are *not* bundled in this git repo (large + third-party redistribution licensing). For out-of-the-box library support, download the full release tarball (`abacuscopilot_v*.tar.gz`, which ships `PP-Orb/`) from the [Releases](../../releases) page and re-run `./setup.sh` from it, or place your own SG15 / lanthanide UPF + orbital files under `PP-Orb/`. The tool warns (but does not fail) when it cannot resolve a species' library files.

## Task Reference

### Pre-processing

| ID | Name | Description |
|----|------|-------------|
| 9901 | System Setup | Configure pseudopotential paths, defaults |
| 9902 | Show Config | Display current configuration |
| 9903 | Check Environment | Verify tools and dependencies |
| 9904 | Clean Directory | Clean up calculation directory |
| 9906 | Set Submit Script | Configure submission-script path |
| 101 | SCF INPUT | Generate INPUT for SCF calculation |
| 102 | Relax INPUT | Generate INPUT for geometry relaxation |
| 103 | MD INPUT | Generate INPUT for molecular dynamics |
| 104 | Band Structure INPUT | Generate INPUT for band (NSCF) |
| 105 | DOS INPUT | Generate INPUT for DOS/PDOS |
| 106 | Work Function INPUT | Generate INPUT for work function |
| 107 | Ecutwfc Test | ecutwfc convergence sweep (PW & LCAO) |
| 108 | Kspacing Test | k-spacing convergence sweep |
| 109 | Convergence Analysis | Analyze ecutwfc/kspacing sweep results |
| 110 | EOS Setup | Equation-of-state setup |
| 111 | Elastic Setup | Elastic-constants setup |
| 112 | Phonon Setup | Phonon calculation setup |
| 113 | PhonoABACUS | Phonon + ABACUS workflow |
| 201 | CIF to STRU | Convert CIF to ABACUS STRU |
| 202 | POSCAR to STRU | Convert VASP POSCAR/CONTCAR to STRU |
| 203 | Coord Convert | Convert STRU Direct ↔ Cartesian |
| 204 | STRU to CIF | Convert STRU to CIF |
| 205 | STRU to POSCAR | Convert STRU to VASP POSCAR |
| 206 | STRU to PDB | Convert STRU to PDB |
| 207 | PDB to STRU | Convert PDB (molecule) to STRU; cubic box if none |
| 208 | STRU to LAMMPS | Convert STRU to LAMMPS data file |
| 209 | LAMMPS to STRU | Convert LAMMPS data file to STRU |
| 210 | View Structure | Open structure in ASE viewer |
| 301 | Auto KPT (MP mesh) | Generate MP k-point mesh from k-spacing |
| 302 | KPT (band path) | Generate line-mode KPT for bands |
| 303 | Custom KPT | Custom k-point list |
| 304 | KPT (phonon) | K-points for phonon calculation |
| 401 | Build Supercell | Nx×Ny×Nz supercell |
| 402 | Redefine Lattice | Apply 3×3 transformation matrix |
| 403 | Coord Conversion | Direct ↔ Cartesian conversion |
| 404 | Fix/Unfix Atoms | Selective dynamics constraints |
| 405 | Slab Builder | Build slab from Miller indices |
| 406 | Vacuum Layer | Add vacuum along cell axis |
| 407 | Shift Atoms | Translate atoms |
| 408 | Sort Atoms | Reorder atoms by species |
| 409 | Reorder Species | Change species order in STRU |
| 410 | XYZ Transform | Rotate lattice so vacuum lies on X/Y |
| 501 | Symmetry Analysis | Space group, Wyckoff positions |
| 502 | Primitive Cell | Reduce to primitive cell |
| 601 | PBS Script | Generate PBS submission script |
| 602 | SLURM Script | Generate SLURM submission script |
| 603 | Validate Inputs | Validate INPUT/STRU/KPT consistency |
| 604 | Batch Submit | Submit batch of calculations |

### Post-processing

| ID | Name | Description |
|----|------|-------------|
| 701 | SCF Convergence | Check SCF convergence history |
| 702 | SCF Compare | Compare multiple SCF runs |
| 703 | Ion Steps | Per-ionic-step summary table |
| 704 | MD Monitor | Live monitor of running MD |
| 801 | Plot Band Structure | Band structure from BANDS_*.dat |
| 802 | Fat-band Plot | Projected (fat) band structure |
| 901 | Plot DOS | Total density of states |
| 902 | Plot PDOS | Projected density of states |
| 903 | DOS+Bands Combined | Side-by-side figure |
| 1001 | 1D Planar Avg Charge | Planar-averaged charge density |
| 1002 | Export Cube/XSF | Export for VESTA visualization |
| 1003 | Diff. Charge Density | Δρ = ρ_AB − ρ_A − ρ_B |
| 1101 | Work Function | ⏳ 待开发 (coming soon) — Φ = V_vacuum − E_Fermi |
| 1102 | Macroscopic Avg | ⏳ 待开发 (coming soon) — Double-averaged potential |
| 1201 | Elastic Constants | Parse elastic tensor, VRH moduli |
| 1202 | EOS Fitting | Birch-Murnaghan equation of state |
| 1301 | Mulliken Analysis | Mulliken population |
| 1401 | Mulliken Bond Order | Mulliken bond order / overlap population (ABACUS out_mul 1) |
| 3101 | Trajectory → PDB | MD_dump / ASE .traj → PDB for VMD/PyMOL |
| 3102 | Extract Frames | Subsample MD trajectory |
| 3103 | MSD | Mean square displacement + diffusion |
| 3104 | RDF | Radial distribution function |
| 3105 | Probability Density | Atomic probability density → CHGCAR |
| 3106 | van Hove | Gs(r,t), Gd(r,t), NGP |
| 3107 | LAMMPS→MD_dump | LAMMPS dump → ABACUS MD_dump |
| 3108 | XDATCAR→MD_dump | VASP XDATCAR → ABACUS MD_dump |
| 3109 | E-T vs Time | Energy/temperature vs time |
| 3201 | Phonon Analysis | Phonon dispersion |
| 3202 | Phonon DOS | Phonon DOS |
| 3203 | Phonon PDOS | Phonon projected DOS |
| 3204 | Phonon Combined | Dispersion + DOS figure |
| 3205 | MLP-SSCHA Setup | MLP-SSCHA workflow setup |
| 3206 | MLP Training | Machine-learned potential training |
| 3207 | SSCHA Run | SSCHA calculation |
| 3208 | SSCHA Plot | SSCHA results plot |
| 3301 | NEB Path (Linear) | Linear NEB path generation |
| 3302 | NEB Path (IDPP) | IDPP NEB path generation |
| 3303 | atst-tools NEB Config | atst-tools NEB configuration |
| 3304 | ASE NEB Script | ASE NEB script |
| 3305 | atst-tools NEB Analysis | atst-tools NEB analysis |
| 3306 | ASE NEB Analysis | ASE NEB analysis |

## CLI Examples

```bash
# Generate all input files for a relaxation
abacuscopilot -task 201                    # STRU from CIF
abacuscopilot -task 102                    # Relax INPUT
abacuscopilot -task 301                    # Auto KPT
abacuscopilot -task 101                    # Generate INPUT (auto-copies PP/orbital files)

# After calculation: analyze results
abacuscopilot -task 703                    # Ionic/electronic step summary
abacuscopilot -task 701                    # SCF convergence
abacuscopilot -task 801 --bands OUT.Si/BANDS_1.dat  # Band plot
abacuscopilot -task 901 --dos OUT.Si/DOS_1.dat      # DOS plot

# Output to a specific directory
abacuscopilot -task 101 -o ./my_calculation/

# Note: for flags with negative values, use = syntax
abacuscopilot -task 801 --erange=-5,5 --bands BANDS.dat
```

## Configuration

Edit `~/.abacuscopilot/config.yaml` or run `abacuscopilot -task 9901`:

```yaml
libraries:
  # Lists of directories, searched in order.  Auto-detected from the
  # bundled PP-Orb/ folder; stale paths are dropped automatically.
  pseudo_library:
    - /path/to/pseudopotentials/
    - /path/to/lanthanides/
  orbital_library:
    - /path/to/orbitals/
    - /path/to/lanthanides/

defaults:
  kspacing: 0.14
  ecutwfc: 100.0
  basis_type: lcao
  dft_functional: pbe
  scf_thr: 1e-7

plotting:
  dpi: 300
  figure_format: png
```

**PP/orbital libraries.** Bundled libraries live under `PP-Orb/` in the project
root (e.g. `PP-Orb/SG15-Version1p0_Pseudopotential`,
`PP-Orb/SG15-Version1p0__StandardOrbitals-Version2p0`,
`PP-Orb/lanthanides-f--core.icmod1`). They are auto-detected on first run and
may be added to/removed freely. `pseudo_library` / `orbital_library` accept a
single directory or a list; each is searched recursively, so the nested APNS
lanthanide layout (`{Element}/{basis}/{Element}_gga_*.orb`) works as-is.
Filenames are matched case-insensitively and tolerate charge-state prefixes
(e.g. `Sm3+_f--core-icmod1.PD04.PBE.UPF`). When a library ships several
orbitals per element, the DZP basis at 7 au (`4s2p2d1f`, `7au`) is preferred.

**Large-core lanthanide PPs.** The APNS `lanthanides-f--core.icmod1` bundle
contains f-electron-pseudized, +3-valent large-core PPs (frozen f-electrons).
When a structure uses one, abacuscopilot prints a warning stating what these
PPs are suited for (structure relaxation / MD / ionic transport; Li3MCl6-type
halides) and NOT suited for (f-electron magnetism/spectroscopy, unaries, which
usually fail to converge). During full calculation setup it also detects the
300 Ry orbital cutoff and offers to raise `ecutwfc` to 300 Ry if you confirm.

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_io.py -v
```

### Adding a New Task

```python
from abacuscopilot.tasks import task

@task(9001, category="My Category", name="My Task",
      description="Description for CLI and menu")
def my_task(args=None, interactive=True):
    console = ...
    # task logic here
```

Tasks are auto-discovered from `abacuscopilot/preprocessing/` and `abacuscopilot/postprocessing/` packages.

## License

GPL-3.0

## Acknowledgments

AbacusCopilot's interaction design (menu-driven, task-numbered workflows) follows the convention popularized by **VASPKIT** (https://vaspkit.com, by Vei Wang & Nan Xu). We adopted a similar task-menu style on purpose: many ABACUS users come from a VASP background and are already familiar with that workflow, so keeping a comparable interaction lowers the barrier to getting started, respects users' existing habits, and reduces the software learning cost (降低上手门槛、尊重用户使用习惯、降低学习成本). The code itself is an original, independent implementation.
