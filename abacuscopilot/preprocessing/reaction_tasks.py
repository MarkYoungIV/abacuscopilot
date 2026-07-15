"""Reaction pathway tasks for abacuscopilot.

Task IDs 1601-1699

NEB (Nudged Elastic Band) path generation between initial and final structures.
Supports linear interpolation and IDPP (Image-Dependent Pair Potential) methods.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.tasks import task


def _read_structure(filepath: str):
    """Read a structure file (STRU, POSCAR, CONTCAR, CIF)."""
    p = Path(filepath)

    # Detect STRU format: no suffix, or filename contains "STRU"
    is_stru = (p.suffix == "" and "STRU" in p.name.upper()) or \
              p.suffix.upper() == ".STRU"

    if is_stru:
        from abacuscopilot.io.stru_file import read_stru
        return read_stru(filepath).to_ase()
    else:
        from ase.io import read as ase_read
        return ase_read(filepath)


def _write_images(images, fmt: str = "STRU"):
    """Write images to 00/, 01/, ... subdirectories in the given format."""
    from abacuscopilot.core.models import Structure

    for i, img in enumerate(images):
        d = Path(f"{i:02d}")
        d.mkdir(exist_ok=True)
        if fmt == "STRU":
            s = Structure.from_ase(img)
            from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
            _write_stru_bare(s, is_lcao=False, filepath=str(d / "STRU"))
        else:
            from ase.io import write as ase_write
            ase_write(d / "POSCAR", img, format="vasp")


def _write_chain_structure(images, cell):
    """Merge all frames into one structure.

    Static atoms (displacement < 0.05 Å across all frames) appear once.
    Moving atoms appear once per frame.
    Cell size unchanged.
    """
    from ase import Atoms as ASEAtoms

    from abacuscopilot.core.models import Structure

    n_frames = len(images)
    n_atoms_per_frame = len(images[0])

    # Identify moving atoms: |r_final - r_init| > 0.05 Å
    moving = np.zeros(n_atoms_per_frame, dtype=bool)
    pos_init = images[0].get_positions()
    pos_final = images[-1].get_positions()
    for i in range(n_atoms_per_frame):
        if np.linalg.norm(pos_final[i] - pos_init[i]) > 0.05:
            moving[i] = True

    all_pos = []
    all_symbols = []
    for i in range(n_atoms_per_frame):
        if moving[i]:
            for img in images:
                all_pos.append(img.get_positions()[i])
                all_symbols.append(img.get_chemical_symbols()[i])
        else:
            all_pos.append(images[0].get_positions()[i])
            all_symbols.append(images[0].get_chemical_symbols()[i])

    big_atoms = ASEAtoms(symbols=all_symbols, positions=all_pos,
                          cell=cell, pbc=True)

    s = Structure.from_ase(big_atoms)
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    _write_stru_bare(s, is_lcao=False, filepath="trj.STRU", is_dp=True)

    from ase.io import write as ase_write
    ase_write("trj.vasp", big_atoms, format="vasp")


def _compute_max_displacement(atoms_init, atoms_final) -> float:
    """Compute max atomic displacement (Å) between two structures.

    For single-atom movement this is the actual displacement.
    For multi-atom/rotation this is the maximum atom-wise displacement.
    Used to estimate NEB image count: n_images ≈ d_max / 0.8.
    """
    pos_i = atoms_init.get_positions()
    pos_f = atoms_final.get_positions()
    diffs = np.linalg.norm(pos_f - pos_i, axis=1)
    return float(np.max(diffs))


def _linear_interpolate(atoms_init, atoms_final, n_images: int):
    """Linear interpolation of Cartesian coordinates between two structures."""
    from ase import Atoms as ASEAtoms

    images = []
    for i in range(n_images + 2):  # including endpoints
        alpha = i / (n_images + 1)
        pos = (1 - alpha) * atoms_init.get_positions() + alpha * atoms_final.get_positions()
        img = ASEAtoms(
            symbols=atoms_init.get_chemical_symbols(),
            positions=pos,
            cell=atoms_init.cell,
            pbc=True,
        )
        images.append(img)
    return images


def _idpp_interpolate(atoms_init, atoms_final, n_images: int):
    """IDPP (Image-Dependent Pair Potential) path generation.

    Uses ASE's idpp_interpolate to minimize bond-length distortion
    between adjacent NEB images.
    """
    from ase.mep import idpp_interpolate

    console = _get_console()
    console.print("  [dim]Optimizing path via IDPP (may take a moment)...[/dim]")
    images = _linear_interpolate(atoms_init, atoms_final, n_images)
    idpp_interpolate(images, fmax=0.1, steps=100, mic=False)
    return images


# =============================================================================
# Task 1601: NEB Path (Linear)
# =============================================================================

@task(1601, category="Reaction Dynamics", name="NEB Path (Linear)",
      description="Generate NEB path images by linear interpolation of coordinates")
def task_neb_linear(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate intermediate NEB images by linear interpolation."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== NEB Path — Linear Interpolation ===[/bold cyan]")
    console.print()

    # Find initial and final structures
    # Auto-detect: STRU_ini/STRU_fin, POSCAR_ini/fin, or relaxed STRU_ION_D
    init_path = "STRU_ini" if Path("STRU_ini").exists() else \
                ("POSCAR_ini" if Path("POSCAR_ini").exists() else "init/OUT.ABACUS/STRU_ION_D")
    final_path = "STRU_fin" if Path("STRU_fin").exists() else \
                 ("POSCAR_fin" if Path("POSCAR_fin").exists() else "final/OUT.ABACUS/STRU_ION_D")

    if interactive:
        init_path = _prompt(console, "Initial structure (POSCAR/STRU/CIF)", init_path)
        final_path = _prompt(console, "Final structure (POSCAR/STRU/CIF)", final_path)

    if not Path(init_path).exists():
        console.print(f"[red]Initial structure not found: {init_path}[/red]")
        return
    if not Path(final_path).exists():
        console.print(f"[red]Final structure not found: {final_path}[/red]")
        return

    try:
        atoms_init = _read_structure(init_path)
        atoms_final = _read_structure(final_path)
    except Exception as e:
        console.print(f"[red]Failed: {e}[/red]")
        return

    console.print(f"  Initial: {len(atoms_init)} atoms")
    console.print(f"  Final:   {len(atoms_final)} atoms")

    d_max = _compute_max_displacement(atoms_init, atoms_final)
    suggested = max(1, int(np.ceil(d_max / 0.8)))
    if suggested % 2 == 0:
        suggested += 1  # odd → middle image can sit on saddle for symmetric reactions
    console.print(f"  Max atomic displacement: {d_max:.4f} Å → suggested {suggested} images (odd)")

    if interactive:
        n_images = int(_prompt(console, "Number of intermediate images (odd recommended)", str(suggested)))
    else:
        n_images = suggested

    console.print(f"  Images: {n_images} intermediate + 2 endpoints")
    console.print("  Method: linear Cartesian interpolation")

    images = _linear_interpolate(atoms_init, atoms_final, n_images)

    # 00/01/... directories are not used by downstream NEB tasks (1603/1604
    # read path_*frames.traj instead).  They exist as a visual aid for
    # developers who want to inspect individual images manually.
    write_dirs = False
    fmt = "STRU"
    if interactive:
        want = _prompt_choice(
            console,
            "Generate per-image directories 00/ 01/ ...? (developer aid, not needed for NEB)",
            ["Yes, generate them", "No (skip)"],
            "No (skip)",
        )
        write_dirs = "Yes" in want
        if write_dirs:
            fmt = _prompt_choice(console, "Output format",
                                 ["STRU (ABACUS)", "POSCAR (VASP)"],
                                 "STRU (ABACUS)")
            fmt = "STRU" if "STRU" in fmt else "POSCAR"

    if write_dirs:
        _write_images(images, fmt)

    # Also write combined multi-frame view files (all atoms stacked in one cell)
    _write_chain_structure(images, atoms_init.cell)

    from ase.io import write as ase_write
    traj_path = f"path_{n_images + 2}frames.traj"
    ase_write(traj_path, images)
    console.print()
    if write_dirs:
        console.print(f"[green]✓ {n_images + 2} images ({fmt}) written to 00/ → {n_images + 1:02d}/[/green]")
    console.print("[green]✓ Chain view: trj.STRU + trj.vasp (all frames in one structure)[/green]")
    console.print(f"[green]✓ Trajectory: {traj_path} (open with task 206)[/green]")
    console.print("  [dim]Copy INPUT and KPT to each subdirectory before running NEB[/dim]")
    console.print("  [dim]Or use task 1603 to generate atst-tools config[/dim]")
    console.print()


# =============================================================================
# Task 1602: NEB Path (IDPP)
# =============================================================================

@task(1602, category="Reaction Dynamics", name="NEB Path (IDPP)",
      description="Generate NEB path images with IDPP pairwise-distance optimization")
def task_neb_idpp(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate NEB images using IDPP (Image-Dependent Pair Potential).

    IDPP minimizes bond-length distortion between adjacent images,
    producing physically reasonable intermediate structures that avoid
    atom clashes.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== NEB Path — IDPP ===[/bold cyan]")
    console.print()

    # Auto-detect: STRU_ini/STRU_fin, POSCAR_ini/fin, or relaxed STRU_ION_D
    init_path = "STRU_ini" if Path("STRU_ini").exists() else \
                ("POSCAR_ini" if Path("POSCAR_ini").exists() else "init/OUT.ABACUS/STRU_ION_D")
    final_path = "STRU_fin" if Path("STRU_fin").exists() else \
                 ("POSCAR_fin" if Path("POSCAR_fin").exists() else "final/OUT.ABACUS/STRU_ION_D")

    if interactive:
        init_path = _prompt(console, "Initial structure (POSCAR/STRU/CIF)", init_path)
        final_path = _prompt(console, "Final structure (POSCAR/STRU/CIF)", final_path)

    if not Path(init_path).exists():
        console.print(f"[red]Initial structure not found: {init_path}[/red]")
        return
    if not Path(final_path).exists():
        console.print(f"[red]Final structure not found: {final_path}[/red]")
        return

    try:
        atoms_init = _read_structure(init_path)
        atoms_final = _read_structure(final_path)
    except Exception as e:
        console.print(f"[red]Failed: {e}[/red]")
        return

    console.print(f"  Initial: {len(atoms_init)} atoms")
    console.print(f"  Final:   {len(atoms_final)} atoms")

    d_max = _compute_max_displacement(atoms_init, atoms_final)
    suggested = max(1, int(np.ceil(d_max / 0.8)))
    if suggested % 2 == 0:
        suggested += 1  # odd → middle image can sit on saddle for symmetric reactions
    console.print(f"  Max atomic displacement: {d_max:.4f} Å → suggested {suggested} images (odd)")

    if interactive:
        n_images = int(_prompt(console, "Number of intermediate images (odd recommended)", str(suggested)))
    else:
        n_images = suggested

    console.print(f"  Images: {n_images} intermediate + 2 endpoints")
    console.print("  Method: IDPP (Image-Dependent Pair Potential)")

    images = _idpp_interpolate(atoms_init, atoms_final, n_images)

    # 00/01/... directories are not used by downstream NEB tasks (1603/1604
    # read path_*frames.traj instead).  They exist as a visual aid for
    # developers who want to inspect individual images manually.
    write_dirs = False
    fmt = "STRU"
    if interactive:
        want = _prompt_choice(
            console,
            "Generate per-image directories 00/ 01/ ...? (developer aid, not needed for NEB)",
            ["Yes, generate them", "No (skip)"],
            "No (skip)",
        )
        write_dirs = "Yes" in want
        if write_dirs:
            fmt = _prompt_choice(console, "Output format",
                                 ["STRU (ABACUS)", "POSCAR (VASP)"],
                                 "STRU (ABACUS)")
            fmt = "STRU" if "STRU" in fmt else "POSCAR"

    if write_dirs:
        _write_images(images, fmt)

    _write_chain_structure(images, atoms_init.cell)

    from ase.io import write as ase_write
    traj_path = f"path_{n_images + 2}frames.traj"
    ase_write(traj_path, images)
    console.print()
    if write_dirs:
        console.print(f"[green]✓ {n_images + 2} images ({fmt}) written to 00/ → {n_images + 1:02d}/[/green]")
    console.print("[green]✓ Chain view: trj.STRU + trj.vasp (all frames in one structure)[/green]")
    console.print(f"[green]✓ Trajectory: {traj_path} (open with task 206)[/green]")
    console.print("  [dim]Copy INPUT and KPT to each subdirectory before running NEB[/dim]")
    console.print("  [dim]Or use task 1603 to generate atst-tools config[/dim]")
    console.print()


# =============================================================================
# Task 1603: atst-tools NEB YAML generator
# =============================================================================

@task(1603, category="Reaction Dynamics", name="atst-tools NEB Config",
      description="Generate atst-tools neb.yaml from NEB image directories")
def task_atst_neb_config(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate an atst-tools compatible YAML configuration for NEB."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate atst-tools NEB Config ===[/bold cyan]")
    console.print()

    image_dirs = sorted([d for d in Path(".").iterdir()
                         if d.is_dir() and d.name.isdigit() and len(d.name) == 2])
    n_images = len(image_dirs)
    first_dir = image_dirs[0] if image_dirs else None

    # If no 00/01/... directories, fall back to the path_*frames.traj that
    # 1601/1602 always produce (the image dirs are optional since v0.1.8).
    from_traj = False
    if n_images < 3:
        traj_files = sorted(Path(".").glob("path_*frames.traj"))
        if not traj_files:
            console.print("[red]No NEB image directories or path_*frames.traj found.[/red]")
            console.print("[dim]Run task 1601 or 1602 first to generate NEB paths.[/dim]")
            return
        from ase.io import read as ase_read
        images = ase_read(str(traj_files[0]), index=":")
        n_images = len(images)
        atoms0 = images[0]
        symbols = atoms0.get_chemical_symbols()
        seen = set()
        species = [s for s in symbols if not (s in seen or seen.add(s))]
        cell_a = atoms0.get_cell().array
        from_traj = True
        console.print(f"  Reading {n_images} images from {traj_files[0].name}")
    else:
        console.print(f"  Found {n_images} image directories: "
                      f"{first_dir.name}/ -> {image_dirs[-1].name}/")
        stru_file = first_dir / "STRU"
        poscar_file = first_dir / "POSCAR"
        if stru_file.exists():
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru(str(stru_file))
            species = structure.species_order
        elif poscar_file.exists():
            from ase.io import read as ase_read
            atoms = ase_read(str(poscar_file))
            symbols = atoms.get_chemical_symbols()
            seen = set()
            species = [s for s in symbols if not (s in seen or seen.add(s))]
        else:
            console.print(f"[red]No STRU or POSCAR found in {first_dir}/[/red]")
            return

    console.print(f"  Species: {', '.join(species)}")

    # Basis type decides the standard规范: lcao writes+copies orb, pw does not.
    from abacuscopilot.core.standards import is_lcao_basis, solver_for
    basis_type = "lcao"
    if interactive:
        basis_type = _prompt_choice(console, "Basis type", ["lcao", "pw"], "lcao")
    is_lcao = is_lcao_basis(basis_type)

    # Read cell for KPT auto-calculation
    if not from_traj:
        cell_a = np.eye(3)
        if stru_file.exists():
            cell_a = structure.lattice.cell_angstrom
        elif poscar_file.exists():
            from ase.io import read as _ase_read
            cell_a = _ase_read(str(poscar_file)).cell.array

    # Resolve pseudo/orbital filenames
    from abacuscopilot.config import load_config
    from abacuscopilot.preprocessing.system_tasks import _find_file_for_element
    config = load_config()
    libs = config.get("libraries", {})
    pseudo_lib = libs.get("pseudo_library", "")
    orbital_lib = libs.get("orbital_library", "")

    pp_map = {}
    orb_map = {}
    if not from_traj and stru_file.exists():
        pp_map = dict(structure.pseudo_files)
        if is_lcao:
            orb_map = dict(structure.orbital_files)
    for sp in species:
        if sp not in pp_map or not pp_map[sp]:
            upfs = sorted(Path(".").glob(f"{sp}_*.upf"))
            if upfs:
                pp_map[sp] = upfs[0].name
            elif pseudo_lib:
                found = _find_file_for_element(pseudo_lib, sp, ".upf")
                pp_map[sp] = found if found else f"{sp}.upf"
            else:
                pp_map[sp] = f"{sp}.upf"
        if not is_lcao:
            continue  # pw: no orbital files (标准规范)
        if sp not in orb_map or not orb_map[sp]:
            orbs = sorted(Path(".").glob(f"{sp}_*.orb"))
            if orbs:
                orb_map[sp] = orbs[0].name
            elif orbital_lib:
                found = _find_file_for_element(orbital_lib, sp, ".orb")
                orb_map[sp] = found if found else f"{sp}.orb"
            else:
                orb_map[sp] = f"{sp}.orb"

    # Auto-copy missing files
    import shutil
    paths_cfg = config.get("paths", {})
    copied = []
    for sp in species:
        pp_file = pp_map.get(sp, f"{sp}.upf")
        if not Path(pp_file).exists() and pseudo_lib:
            src = Path(pseudo_lib) / pp_file
            if src.exists():
                shutil.copy2(src, ".")
                copied.append(pp_file)
        orb_file = orb_map.get(sp) if is_lcao else None
        if orb_file and not Path(orb_file).exists() and orbital_lib:
            src = Path(orbital_lib) / orb_file
            if src.exists():
                shutil.copy2(src, ".")
                copied.append(orb_file)
    if copied:
        console.print(f"  [green]Copied: {', '.join(copied)}[/green]")

    # Distribute to image directories
    for d in image_dirs:
        for sp in species:
            pp = pp_map.get(sp, f"{sp}.upf")
            if Path(pp).exists() and not (d / pp).exists():
                shutil.copy2(pp, d / pp)
            orb = orb_map.get(sp) if is_lcao else None
            if orb and Path(orb).exists() and not (d / orb).exists():
                shutil.copy2(orb, d / orb)
    console.print(f"  [green]Files distributed to {len(image_dirs)} image dirs[/green]")

    abacus_bin = paths_cfg.get("abacus_binary", "abacus")

    # Interactive prompts
    use_kpt = False
    kspacing = "0.14"
    kpt_grid = [1, 1, 1]
    climb = True
    do_two_stage = True
    # Basis + hardware defaults. Standard规范 (ABACUS manual):
    #   pw   -> dav_subspace (cpu/gpu);  lcao -> genelpa (cpu) / cusolver (gpu single card)
    is_gpu = False
    device = "cpu"
    n_mpi = 8
    n_omp = 1

    if interactive:
        hw = _prompt_choice(console, "Target hardware",
                            ["CPU (single node)", "GPU (single card)"],
                            "CPU (single node)")
        is_gpu = "GPU" in hw
        if is_gpu:
            device = "gpu"
            n_mpi = 1             # cusolver uses one GPU card
            n_omp = int(_prompt(console, "OMP threads per image", "12"))
        else:
            device = "cpu"
            n_mpi = int(_prompt(console, "MPI cores per image", "8"))
            n_omp = 1
        kpt_mode = _prompt_choice(console, "K-point mode",
                                   ["kspacing (auto mesh)", "KPT (explicit grid)"],
                                   "kspacing (auto mesh)")
        use_kpt = "KPT" in kpt_mode
        if use_kpt:
            a = np.linalg.norm(cell_a[0])
            b = np.linalg.norm(cell_a[1])
            c = np.linalg.norm(cell_a[2])
            nkx = max(1, int(round(30.0 / a)))
            nky = max(1, int(round(30.0 / b)))
            nkz = max(1, int(round(30.0 / c)))
            kpt_default = f"{nkx} {nky} {nkz}"
            inp = _prompt(console, f"KPT grid (nkx nky nkz) [auto: {kpt_default}]", kpt_default)
            try:
                parts = [int(x) for x in inp.split()]
                kpt_grid = parts[:3] if len(parts) >= 3 else [nkx, nky, nkz]
            except ValueError:
                kpt_grid = [nkx, nky, nkz]
        else:
            kspacing = _prompt(console, "K-spacing (2pi/A)", "0.14")
        climbing = _prompt_choice(console, "CI-NEB (climbing image)?",
                                   ["Yes", "No"], "Yes")
        climb = "Yes" in climbing
        two_stage = _prompt_choice(console, "Two-stage optimization? (coarse -> fine)",
                                    ["Yes", "No"], "Yes")
        do_two_stage = "Yes" in two_stage

    solver = solver_for(basis_type, device)

    # Use existing .traj if available
    traj_files = sorted(Path(".").glob("path_*frames.traj"))
    init_chain = str(traj_files[0]) if traj_files else [str(d / "POSCAR") for d in image_dirs]
    if traj_files:
        console.print(f"  [dim]Using: {init_chain}[/dim]")

    # Build YAML config
    if not from_traj:
        n_images = len(image_dirs)
    calc_section = {
        "type": "neb",
        "init_chain": init_chain,
        "climb": climb,
    }
    if do_two_stage:
        calc_section["two_stage"] = True
        calc_section["stage1_steps"] = 20
        calc_section["stage1_fmax"] = 0.20
    calc_section["fmax"] = 0.05
    calc_section["k"] = 0.1
    calc_section["algorism"] = "improvedtangent"
    calc_section["parallel"] = not is_gpu
    calc_section["optimizer"] = "FIRE"
    calc_section["max_steps"] = 200

    abacus_params = {
        "calculation": "scf",
        "ecutwfc": 100,
        "basis_type": basis_type,
        "device": device,
        "ks_solver": solver,
        "dft_functional": "pbe",
        "scf_thr": 1e-7,
        "scf_nmax": 100,
        "smearing_method": "gaussian",
        "smearing_sigma": 0.001,
        "mixing_type": "broyden",
    }
    if not use_kpt:
        abacus_params["kspacing"] = float(kspacing)
    abacus_params.update({
        "cal_force": 1,
        "cal_stress": 1,
        "init_wfc": "atomic",
        "init_chg": "atomic",
        "out_stru": 1,
        "out_chg": 0,
        "out_mul": 0,
        "pseudo_dir": "./",
        "pseudopotentials": pp_map,
    })
    if is_lcao:
        abacus_params["out_wfc_lcao"] = 0
        abacus_params["orbital_dir"] = "./"
        abacus_params["basissets"] = orb_map

    abacus_section = {
        "command": abacus_bin,
        "mpi": n_mpi,
        "omp": n_omp,
        "directory": "run_neb",
    }
    if use_kpt:
        abacus_section["kpts"] = kpt_grid
    abacus_section["parameters"] = abacus_params

    import yaml
    neb_config = {
        "calculation": calc_section,
        "calculator": {
            "name": "abacus",
            "abacus": abacus_section,
        },
    }

    yaml_out = yaml.dump(neb_config, default_flow_style=False, sort_keys=False,
                         allow_unicode=True)
    # Post-process: blank line before calculator, inline kpts
    yaml_out = yaml_out.replace("\ncalculator:", "\n\ncalculator:")
    # Make kpts inline: "[2, 2, 2]" instead of "- 2\n- 2\n- 2"
    import re
    yaml_out = re.sub(r'kpts:\n(\s+- \d+\n)+', lambda m: 'kpts: [' +
                      ', '.join(re.findall(r'- (\d+)', m.group())) + ']\n', yaml_out)
    with open("neb.yaml", "w") as f:
        f.write(yaml_out)

    console.print()
    console.print("[green]NEB config written: neb.yaml[/green]")
    if is_gpu:
        console.print(f"  Images: {len(image_dirs)} | Basis: {basis_type} | Hardware: GPU "
                      f"(device=gpu, ks_solver={solver}, mpi={n_mpi}, omp={n_omp}, parallel=False)")
    else:
        console.print(f"  Images: {len(image_dirs)} | Basis: {basis_type} | Hardware: CPU "
                      f"(ks_solver={solver}, mpi={n_mpi}, omp={n_omp}, parallel=True)")
    if use_kpt:
        console.print(f"  K-points: {kpt_grid[0]}x{kpt_grid[1]}x{kpt_grid[2]} (explicit)")
    else:
        console.print(f"  K-spacing: {kspacing} (auto mesh)")
    console.print(f"  CI-NEB: {'Yes' if climb else 'No'}")
    console.print(f"  Two-stage: {'Yes (20 coarse + 200 fine)' if do_two_stage else 'No (200 fine only)'}")
    console.print("  Optimizer: FIRE (fmax=0.05, max_steps=200)")
    console.print("  [dim]pseudo_dir/orbital_dir = ./ (files auto-copied)[/dim]")
    console.print()
    console.print("  [bold]Next steps:[/bold]")
    console.print("    1. Review neb.yaml and edit as needed")
    console.print("    2. atst run neb.yaml")
    console.print("  [dim]atst-tools is installed automatically as a core dependency[/dim]")
    console.print("  [dim]Or use task 1604 for a self-contained ASE script (no atst-tools needed)[/dim]")
    console.print()

    # Copy the Slurm sbatch template and adapt it for atst-tools NEB
    sub_template = config.get("paths", {}).get("sub_script", "")
    if sub_template:
        sub_src = Path(sub_template)
        if sub_src.exists():
            content = sub_src.read_text()
            # Use the absolute path to atst so the SLURM job doesn't need
            # conda activate (atst lives in the abacuscopilot env).
            import shutil as _shutil
            atst_bin = _shutil.which("atst") or "atst"
            import re
            # Common patterns: "abacus", "mpirun -np N abacus", "srun abacus"
            content = re.sub(
                r"^(mpirun\s+.*\s+)?abacus\b.*$",
                f"{atst_bin} run neb.yaml",
                content,
                flags=re.MULTILINE,
            )
            # Also handle "srun abacus"
            content = re.sub(
                r"^srun\s+abacus\b.*$",
                f"{atst_bin} run neb.yaml",
                content,
                flags=re.MULTILINE,
            )
            out_name = "sub.abacus_neb_atst.sh"
            with open(out_name, "w") as f:
                f.write(content)
            console.print(f"  [green]✓ {out_name} copied from template[/green] (→ {atst_bin} run neb.yaml)")
        else:
            console.print(f"  [yellow]! sub_script not found: {sub_template}[/yellow]")
    console.print()


# =============================================================================
# Task 1604: ASE NEB standalone script
# =============================================================================

@task(1604, category="Reaction Dynamics", name="ASE NEB Script",
      description="Generate a self-contained Python script for ASE+abacuslite NEB")
def task_ase_neb_script(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a standalone Python script that runs NEB via ASE + abacuslite.

    No external tools needed — only ASE (pip install ase) and abacuslite
    (shipped with ABACUS source under interfaces/ASE_interface/).
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate ASE NEB Script ===[/bold cyan]")
    console.print()

    image_dirs = sorted([d for d in Path(".").iterdir()
                         if d.is_dir() and d.name.isdigit() and len(d.name) == 2])
    from_traj = False

    if len(image_dirs) < 3:
        traj_files = sorted(Path(".").glob("path_*frames.traj"))
        if not traj_files:
            console.print("[red]No NEB image directories or path_*frames.traj found.[/red]")
            return
        from ase.io import read as ase_read
        images_traj = ase_read(str(traj_files[0]), index=":")
        atoms0 = images_traj[0]
        symbols = atoms0.get_chemical_symbols()
        seen = set()
        species = [s for s in symbols if not (s in seen or seen.add(s))]
        # Create a dummy structure from the first frame for cell info
        from abacuscopilot.core.models import Structure
        structure = Structure.from_ase(atoms0)
        stru_file = Path()
        poscar_file = Path()
        from_traj = True
        console.print(f"  Reading {len(images_traj)} images from {traj_files[0].name}")
    else:
        console.print(f"  Found {len(image_dirs)} images: "
                      f"{image_dirs[0].name}/ -> {image_dirs[-1].name}/")
        first_dir = image_dirs[0]
        stru_file = first_dir / "STRU"
        poscar_file = first_dir / "POSCAR"
        if stru_file.exists():
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru(str(stru_file))
            species = structure.species_order
        elif poscar_file.exists():
            from ase.io import read as ase_read
            atoms = ase_read(str(poscar_file))
            symbols = atoms.get_chemical_symbols()
            seen = set()
            species = [s for s in symbols if not (s in seen or seen.add(s))]
        else:
            console.print(f"[red]No STRU or POSCAR in {first_dir}/[/red]")
            return

    from abacuscopilot.config import load_config
    config = load_config()
    libs = config.get("libraries", {})
    paths_cfg = config.get("paths", {})
    pseudo_lib = libs.get("pseudo_library", "")
    orbital_lib = libs.get("orbital_library", "")
    abacus_bin = paths_cfg.get("abacus_binary", "abacus")
    mpirun = paths_cfg.get("mpirun", "mpirun")
    abacus_src = paths_cfg.get("abacus_source", "")
    slurm_env = paths_cfg.get("slurm_env_file", "")
    n_total = len(images_traj) if from_traj else len(image_dirs)

    if interactive:
        console.print(f"  [dim]SLURM env script (CUDA/compiler/conda): {slurm_env or 'not set'}[/dim]")
        env_in = console.input("  Change? [Enter=keep]: ").strip()
        if env_in:
            slurm_env = env_in
            paths_cfg["slurm_env_file"] = env_in
            from abacuscopilot.config import save_config
            save_config(config)

    # Basis type decides the standard规范 (lcao writes orb, pw does not).
    from abacuscopilot.core.standards import is_lcao_basis, solver_for
    basis_type = "lcao"
    if interactive:
        basis_type = _prompt_choice(console, "Basis type", ["lcao", "pw"], "lcao")
    is_lcao = is_lcao_basis(basis_type)

    # Build PP/orb maps — resolve real filenames from library (标准规范).
    # Same logic as 1603: current dir first, then library, then bare fallback.
    from abacuscopilot.preprocessing.system_tasks import _find_file_for_element
    pp_map = {}
    orb_map = {}
    if stru_file.exists():
        pp_map = dict(structure.pseudo_files)
        if is_lcao:
            orb_map = dict(structure.orbital_files)
    for sp in species:
        if sp not in pp_map or not pp_map[sp]:
            upfs = sorted(Path(".").glob(f"{sp}_*.upf"))
            if upfs:
                pp_map[sp] = upfs[0].name
            elif pseudo_lib:
                found = _find_file_for_element(pseudo_lib, sp, ".upf")
                pp_map[sp] = found if found else f"{sp}.upf"
            else:
                pp_map[sp] = f"{sp}.upf"
        if is_lcao and (sp not in orb_map or not orb_map[sp]):
            orbs = sorted(Path(".").glob(f"{sp}_*.orb"))
            if orbs:
                orb_map[sp] = orbs[0].name
            elif orbital_lib:
                found = _find_file_for_element(orbital_lib, sp, ".orb")
                orb_map[sp] = found if found else f"{sp}.orb"
            else:
                orb_map[sp] = f"{sp}.orb"

    # Interactive config
    # Hardware defaults (CPU single node). device=cpu -> genelpa;
    # device=gpu (single card) -> cusolver (ABACUS manual).
    is_gpu = False
    device = "cpu"
    n_mpi = 8
    n_omp = 1
    if interactive:
        hw = _prompt_choice(console, "Target hardware",
                            ["CPU (single node)", "GPU (single card)"],
                            "CPU (single node)")
        is_gpu = "GPU" in hw
        if is_gpu:
            device = "gpu"
            n_mpi = 1             # cusolver uses one GPU card
            n_omp = int(_prompt(console, "OMP threads per image", "12"))
        else:
            device = "cpu"
            n_mpi = int(_prompt(console, "MPI cores per image", "8"))
            n_omp = 1
        kpt_mode = _prompt_choice(console, "K-point mode",
                                   ["kspacing (auto mesh)", "KPT (explicit grid)"],
                                   "kspacing (auto mesh)")
        use_kpt = "KPT" in kpt_mode
        if use_kpt:
            if from_traj:
                cell_a = structure.lattice.cell_angstrom
            elif stru_file.exists():
                cell_a = structure.lattice.cell_angstrom
            elif poscar_file.exists():
                cell_a = _ase_read(str(poscar_file)).cell.array
            else:
                cell_a = np.eye(3)
            a, b, c = np.linalg.norm(cell_a[0]), np.linalg.norm(cell_a[1]), np.linalg.norm(cell_a[2])
            nkx = max(1, int(round(30.0 / a)))
            nky = max(1, int(round(30.0 / b)))
            nkz = max(1, int(round(30.0 / c)))
            kpt_str = f"{nkx} {nky} {nkz}"
            inp = _prompt(console, f"KPT grid (nkx nky nkz) [auto: {kpt_str}]", kpt_str)
            try:
                parts = [int(x) for x in inp.split()]
                kpt_grid = parts[:3] if len(parts) >= 3 else [nkx, nky, nkz]
            except ValueError:
                kpt_grid = [nkx, nky, nkz]
            kspacing_val = 0.14
        else:
            kspacing_val = float(_prompt(console, "K-spacing (2pi/A)", "0.14"))
            kpt_grid = None
        climbing = _prompt_choice(console, "CI-NEB?", ["Yes", "No"], "Yes")
        climb = "Yes" in climbing
        two_stage = _prompt_choice(console, "Two-stage optimization?", ["Yes", "No"], "Yes")
        do_two_stage = "Yes" in two_stage
    else:
        use_kpt = False
        kspacing_val = 0.14
        kpt_grid = None
        climb = True
        do_two_stage = True

    solver = solver_for(basis_type, device)

    # Use .traj if available
    traj_files = sorted(Path(".").glob("path_*frames.traj"))
    traj_path = str(traj_files[0]) if traj_files else None
    if traj_files:
        console.print(f"  [dim]Using: {traj_path}[/dim]")

    # Build the Python script
    pp_block = "    " + "\n    ".join(f'"{sp}": "{pp_map[sp]}",' for sp in species)

    # LCAO-only blocks (标准规范: pw writes no orbital info)
    if is_lcao:
        orb_block = "    " + "\n    ".join(f'"{sp}": "{orb_map[sp]}",' for sp in species)
        orbital_dir_line = f'ORBITAL_DIR = "{orbital_lib or "./"}"\n'
        profile_orb_line = "    orbital_dir=ORBITAL_DIR,\n"
        basissets_block = f"basissets = {{\n{orb_block}\n}}\n"
        inp_orbital_lines = '    "out_wfc_lcao": 0,\n'
        calc_basissets_arg = "        basissets=basissets,\n"
    else:
        orbital_dir_line = ""
        profile_orb_line = ""
        basissets_block = ""
        inp_orbital_lines = ""
        calc_basissets_arg = ""

    kpt_lines = ""
    if use_kpt:
        kpt_lines = f'        "kpts": {kpt_grid},\n'
    else:
        kpt_lines = f'        "kspacing": {kspacing_val},\n'

    two_stage_block = ""
    par = "True" if not is_gpu else "False"
    if do_two_stage:
        two_stage_block = f"""neb = NEB(images, climb=True, k=0.1, parallel={par},
            method="improvedtangent")
# Stage 1: coarse optimization
opt1 = FIRE(neb, trajectory="neb_stage1.traj")
opt1.run(fmax=0.20, steps=20)
print("Stage 1 done — starting stage 2 (fine)")
# Stage 2: fine optimization
opt2 = FIRE(neb, trajectory="neb_stage2.traj")
opt2.run(fmax=0.05, steps=200)
"""
    else:
        two_stage_block = f"""neb = NEB(images, climb=True, k=0.1, parallel={par},
            method="improvedtangent")
opt = FIRE(neb, trajectory="neb.traj")
opt.run(fmax=0.05, steps=200)
"""

    script = f'''#!/usr/bin/env python3
"""ASE NEB script for ABACUS — generated by AbacusCopilot task 1604.

Dependencies: ase, abacuslite (from ABACUS source: interfaces/ASE_interface/)
Usage:
  python neb_run.py               # run NEB
  python neb_run.py --check       # validate structure only (no ABACUS needed)
"""
import os, sys
import numpy as np
from ase.io import read

# === Check mode: validate trajectory without ABACUS ===
CHECK_ONLY = "--check" in sys.argv
if CHECK_ONLY:
    TRAJ_FILE = "{traj_path}"
    print(f"Checking trajectory: {{TRAJ_FILE}}")
    images = read(TRAJ_FILE, index=":")
    n = len(images)
    print(f"  Images: {{n}}")
    for i, img in enumerate(images):
        print(f"  Image {{i}}: {{len(img)}} atoms, "
              f"cell={{img.cell.cellpar()[:3].round(3)}}, "
              f"formula={{img.get_chemical_formula()}}")
    # Energy span estimate (final - initial displacement)
    d = np.linalg.norm(images[-1].get_positions() - images[0].get_positions(), axis=1)
    print(f"  Max displacement: {{np.max(d):.4f}} A")
    print(f"  Hardware: {device} | MPI={n_mpi}/image, OMP={n_omp}, parallel={par}")
    print("Validation passed — ready to run on server with ABACUS.")
    sys.exit(0)

from ase.mep import NEB
from ase.optimize import FIRE

# abacuslite is bundled with AbacusCopilot (interfaces/ASE_interface/)
from abacuscopilot.interfaces.ASE_interface.abacuslite import Abacus, AbacusProfile

# === Configuration ===
N_MPI = {n_mpi}
N_OMP = {n_omp}
TRAJ_FILE = "{traj_path}"
PSEUDO_DIR = "{pseudo_lib or './'}"
{orbital_dir_line}
profile = AbacusProfile(
    command="{mpirun} -np {n_mpi} {abacus_bin}",
    omp_num_threads={n_omp},
    pseudo_dir=PSEUDO_DIR,
{profile_orb_line})

pseudopotentials = {{
{pp_block}
}}

{basissets_block}
inp_params = {{
    "calculation": "scf",
    "ecutwfc": 100,
    "basis_type": "{basis_type}",
    "device": "{device}",
    "ks_solver": "{solver}",
    "dft_functional": "pbe",
    "scf_thr": 1e-7,
    "scf_nmax": 100,
    "smearing_method": "gaussian",
    "smearing_sigma": 0.001,
    "mixing_type": "broyden",
    "cal_force": 1,
    "cal_stress": 1,
    "init_wfc": "atomic",
    "init_chg": "atomic",
    "out_stru": 1,
    "out_chg": 0,
    "out_mul": 0,
{inp_orbital_lines}{kpt_lines}}}

# === Read NEB chain ===
images = read(TRAJ_FILE, index=":")
n_images = len(images)
print(f"Loaded {{n_images}} images from {{TRAJ_FILE}}")

# === Attach ABACUS calculators ===
for i, img in enumerate(images):
    img.calc = Abacus(
        profile=profile,
        directory=f"neb-{{i}}",
        pseudopotentials=pseudopotentials,
{calc_basissets_arg}        inp=inp_params,
    )
    print(f"  Image {{i}}: directory=neb-{{i}}")

# === Run NEB ===
print(f"Running NEB ({{n_images}} images, parallel={par})...")
{two_stage_block}

# === Collect energies ===
energies = []
for i, img in enumerate(images):
    e = img.get_potential_energy()
    energies.append(e)
    print(f"  Image {{i}}: E = {{e:.6f}} eV")

energies = np.array(energies)
barrier = np.max(energies) - energies[0]
print(f"\\nEnergy barrier: {{barrier:.4f}} eV")
print("Done! Trajectory saved to neb.traj")
'''

    out_path = "neb_run.py"
    with open(out_path, "w") as f:
        f.write(script)
    import os
    os.chmod(out_path, 0o755)

    console.print()
    # Also generate SLURM script
    n_active = n_total - 2  # exclude endpoints
    # SLURM environment setup — driven by config so each server can use its
    # own CUDA/compiler/ABACUS env script.  When slurm_env_file is unset,
    # the script has a commented placeholder for the user to fill in.
    if slurm_env:
        slurm_env_block = f"source {slurm_env}\n"
    else:
        slurm_env_block = (
            "# === Environment: source your CUDA + ABACUS setup (configure in\n"
            "#     abacuscopilot config → System Setup → 配置生成) ===\n"
        )
    # abacuslite is now bundled with AbacusCopilot — no extra PYTHONPATH needed.
    slurm_pp_block = ""

    if is_gpu:
        # GPU single card: one MPI task, one GPU. NEB images run sequentially
        # (parallel=False) since a single card can't split across images.
        slurm_script = f'''#!/bin/bash
#SBATCH --job-name=neb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={n_omp}
#SBATCH --gres=gpu:1
#SBATCH --output=neb_%j.out
#SBATCH --error=neb_%j.err

{slurm_env_block}export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CUDA_VISIBLE_DEVICES=0
{slurm_pp_block}
python neb_run.py
'''
    else:
        # CPU: one MPI task per active image, n_mpi cores each (parallel=True).
        slurm_script = f'''#!/bin/bash
#SBATCH --job-name=neb
#SBATCH --nodes=1
#SBATCH --ntasks={n_active}
#SBATCH --cpus-per-task={n_mpi}
#SBATCH --output=neb_%j.out
#SBATCH --error=neb_%j.err

{slurm_env_block}export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
{slurm_pp_block}
python neb_run.py
'''
    slurm_path = "neb_slurm.sh"
    with open(slurm_path, "w") as f:
        f.write(slurm_script)
    import os
    os.chmod(slurm_path, 0o755)

    console.print(f"[green]ASE NEB script written: {out_path}[/green]")
    console.print(f"[green]SLURM script written: {slurm_path}[/green]")
    if is_gpu:
        console.print(f"  Images: {n_total} ({n_active} active) | Hardware: GPU "
                      f"(device=gpu, ks_solver={solver}, omp={n_omp}, parallel=False)")
    else:
        total_cores = n_mpi * n_active
        console.print(f"  Images: {n_total} ({n_active} active) | Hardware: CPU "
                      f"(ks_solver={solver}, parallel=True)")
        console.print(f"  Cores: {n_mpi}/image × {n_active} = {total_cores} total")
    if use_kpt and kpt_grid:
        console.print(f"  K-points: {kpt_grid[0]}x{kpt_grid[1]}x{kpt_grid[2]}")
    else:
        console.print(f"  K-spacing: {kspacing_val}")
    console.print(f"  CI-NEB: {'Yes' if climb else 'No'}")
    console.print(f"  Two-stage: {'Yes' if do_two_stage else 'No'}")
    console.print()
    console.print("  [bold]Usage (local):[/bold]")
    console.print(f"    python {out_path}")
    console.print("  [bold]Usage (SLURM):[/bold]")
    console.print(f"    sbatch {slurm_path}")
    console.print("  [dim]Requires: ase + abacuslite + ABACUS toolchain[/dim]")
    console.print()


# =============================================================================
# Task 1605: NEB result analysis
# =============================================================================

@task(1605, category="Reaction Dynamics", name="atst-tools NEB Analysis",
      description="Analyze atst-tools NEB results (from task 1603): barrier, saddle, convergence")
def task_atst_neb_analysis(args: list[str] | None = None, interactive: bool = True) -> None:
    """Analyze atst-tools NEB trajectory and plot the energy barrier profile."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== NEB Result Analysis ===[/bold cyan]")
    console.print()

    # Find NEB trajectory
    traj_path = None
    for candidate in ("neb.traj", "neb_stage2.traj"):
        p = Path(candidate)
        if p.exists():
            traj_path = str(p)
            break
    if traj_path is None:
        for p in sorted(Path(".").glob("neb*.traj")):
            if "idpp" not in p.name:  # skip idpp.traj (optimization history)
                traj_path = str(p)
                break

    if interactive and traj_path:
        console.print(f"  [dim]Default: {traj_path}[/dim]")
        inp = console.input("  NEB trajectory file: ").strip()
        if inp:
            traj_path = inp
    elif not traj_path:
        console.print("  [dim]Looking for neb.traj...[/dim]")
        inp = console.input("  NEB trajectory file: ").strip()
        traj_path = inp if inp else "neb.traj"

    if not traj_path or not Path(traj_path).exists():
        console.print("[red]NEB trajectory not found.[/red]")
        console.print("[dim]Look for neb.traj or path_*frames.traj in the working directory.[/dim]")
        return

    try:
        from ase.io import read as ase_read
        all_frames = ase_read(traj_path, index=":")
    except Exception as e:
        console.print(f"[red]Failed to read trajectory: {e}[/red]")
        return

    n_total = len(all_frames)
    console.print(f"  Total frames in trajectory: {n_total}")

    # Detect NEB chain size: read init_chain from neb.yaml if present
    n_chain = None
    yaml_path = Path("neb.yaml")
    if yaml_path.exists():
        import yaml
        try:
            cfg = yaml.safe_load(yaml_path.read_text())
            init = cfg.get("calculation", {}).get("init_chain", None)
            if isinstance(init, list):
                n_chain = len(init)
            elif isinstance(init, str) and Path(init).exists():
                chain = ase_read(init, index=":")
                n_chain = len(chain)
        except Exception:
            pass

    if n_chain is None and interactive:
        n_chain = int(_prompt(console, "Number of images in NEB chain", "5"))

    if n_chain is None or n_chain < 2 or n_chain > n_total:
        # Fallback: auto-detect by finding repeating atom count pattern
        n_atoms = len(all_frames[0])
        for period in range(2, min(n_total // 2 + 1, 50)):
            if n_total % period == 0 and all(len(f) == n_atoms for f in all_frames[-period:]):
                n_chain = period
                break

    if n_chain is None or n_chain > n_total:
        console.print("[red]Cannot determine NEB chain size.[/red]")
        console.print("[dim]neb.traj contains optimization history — specify number of images.[/dim]")
        return

    # Take only the final converged chain (last n_chain frames)
    images = all_frames[-n_chain:]
    console.print(f"  NEB chain: {n_chain} images (using last {n_chain} of {n_total} frames)")
    if n_chain < 2:
        console.print("[red]Need at least 2 images.[/red]")
        return

    # Extract energies
    energies = np.zeros(len(images))
    for i, img in enumerate(images):
        try:
            energies[i] = img.get_potential_energy()
        except Exception:
            energies[i] = np.nan

    # Relative to initial
    if np.isnan(energies[0]):
        console.print("[red]Energies not available — run with ABACUS calculator first.[/red]")
        return

    e_rel = energies - energies[0]
    barrier = np.max(e_rel)
    saddle_idx = int(np.argmax(e_rel))

    console.print(f"  Initial energy: {energies[0]:.6f} eV")
    console.print(f"  Final energy:   {energies[-1]:.6f} eV")
    console.print(f"  ΔE (final-init): {energies[-1] - energies[0]:.4f} eV")
    console.print(f"  Energy barrier: {barrier:.4f} eV")
    console.print(f"  Saddle point:   image {saddle_idx} (E = {energies[saddle_idx]:.6f} eV)")

    # Cubic spline interpolation (VTST neb.results.pl style)
    x_raw = np.arange(len(images))
    n_img = len(images)
    try:
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(x_raw, e_rel, bc_type="natural")
        x_fine = np.linspace(0, n_img - 1, (n_img - 1) * 20 + 1)
        e_fine = cs(x_fine)
        # Find precise saddle point from spline
        spline_max_idx = np.argmax(e_fine)
        spline_saddle_x = x_fine[spline_max_idx]
        spline_saddle_e = e_fine[spline_max_idx]
        has_spline = True
    except Exception:
        has_spline = False
        spline_saddle_x = float(saddle_idx)
        spline_saddle_e = e_rel[saddle_idx]

    if has_spline:
        console.print(f"  Spline-fitted saddle: image {spline_saddle_x:.2f}, barrier = {spline_saddle_e:.4f} eV")
    console.print()

    # Plot options
    show_saddle = True
    if interactive:
        show_saddle = "Yes" in _prompt_choice(console, "Mark saddle point on plot?",
                                               ["Yes", "No"], "Yes")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    # Raw data points
    ax.plot(x_raw, e_rel, "o", color="#1f77b4", markersize=8, zorder=5,
            label="NEB images")
    # Spline curve
    if has_spline:
        ax.plot(x_fine, e_fine, "-", color="#1f77b4", linewidth=1.2, alpha=0.7,
                label="Cubic spline")
        if show_saddle:
            ax.axvline(x=spline_saddle_x, color="#d62728", linestyle="--", linewidth=0.8,
                       label=f"Saddle: {spline_saddle_e:.3f} eV @ image {spline_saddle_x:.2f}")
    else:
        ax.plot(x_raw, e_rel, "-", color="#1f77b4", linewidth=1.2)
        if show_saddle:
            ax.axvline(x=saddle_idx, color="#d62728", linestyle="--", linewidth=0.8,
                       label=f"Saddle: {e_rel[saddle_idx]:.3f} eV @ image {saddle_idx}")
    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_xlabel("NEB image index")
    ax.set_ylabel("Relative energy (eV)")
    title = f"NEB Energy Barrier — {barrier:.4f} eV"
    if has_spline and show_saddle:
        title += f" (spline: {spline_saddle_e:.4f})"
    ax.set_title(title)
    ax.legend()
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)

    out_png = "neb_barrier.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save data (raw + spline)
    out_dat = "neb_barrier.dat"
    with open(out_dat, "w") as f:
        f.write("# image  E_rel(eV)  E_spline(eV)\n")
        for i in range(len(x_raw)):
            e_s = float(cs(x_raw[i])) if has_spline else e_rel[i]
            f.write(f"{i}  {e_rel[i]:.8f}  {e_s:.8f}\n")

    console.print(f"[green]✓ Energy profile: {out_png}[/green]")
    console.print(f"[green]✓ Data file: {out_dat}[/green]")
    console.print(f"  Barrier: {barrier:.4f} eV, Saddle: image {saddle_idx}/{n_img - 1}")
    if has_spline:
        console.print(f"  Spline-fit: barrier = {spline_saddle_e:.4f} eV at image {spline_saddle_x:.2f}")
    console.print()


# =============================================================================
# Task 1606: ASE NEB Analysis
# =============================================================================

@task(1606, category="Reaction Dynamics", name="ASE NEB Analysis",
      description="Analyze ASE NEB results (from task 1604): convergence, barrier, forces")
def task_ase_neb_analysis(args: list[str] | None = None, interactive: bool = True) -> None:
    """Analyze ASE NEB results generated by task 1604 (ASE NEB Script).

    Auto-detects neb_run.py + neb_stage2.traj, reads the final converged NEB
    chain, reports the energy barrier with cubic-spline interpolation, and
    produces a publication-quality barrier plot.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== ASE NEB Result Analysis ===[/bold cyan]")
    console.print()

    # --- Auto-detect 1604 output ---
    if not Path("neb_run.py").exists():
        console.print("[red]No neb_run.py found — not a task 1604 output directory.[/red]")
        console.print("[dim]Run task 1604 first, then run this analysis in the same directory.[/dim]")
        return
    console.print("  [dim]Detected task 1604 output[/dim]")

    # Find trajectory: prefer two-stage, fall back to single
    traj_path = None
    stage_label = ""
    for candidate, label in [
        ("neb_stage2.traj", "Stage 2 (fine)"),
        ("neb_stage1.traj", "Stage 1 (coarse)"),
        ("neb.traj", "single-stage"),
    ]:
        if Path(candidate).exists():
            traj_path = candidate
            stage_label = label
            break
    if traj_path is None:
        console.print("[red]No NEB trajectory found (neb*.traj).[/red]")
        return

    try:
        from ase.io import read as ase_read
        all_frames = ase_read(traj_path, index=":")
    except Exception as e:
        console.print(f"[red]Failed to read trajectory: {e}[/red]")
        return

    n_total = len(all_frames)
    console.print(f"  Trajectory: {traj_path} ({n_total} frames, {stage_label})")

    # --- Determine NEB chain size from the initial path file ---
    n_chain = None
    for pattern in ("path_*frames.traj", "path_*.traj"):
        matches = sorted(Path(".").glob(pattern))
        if matches:
            try:
                chain = ase_read(str(matches[0]), index=":")
                n_chain = len(chain)
                console.print(f"  Initial chain: {matches[0].name} ({n_chain} images)")
            except Exception:
                pass
            break

    if n_chain is None and interactive:
        n_chain = int(_prompt(console, "Number of images in NEB chain", "7"))
    if n_chain is None or n_chain < 2 or n_chain > n_total:
        console.print("[red]Cannot determine NEB chain size.[/red]")
        return

    # --- Extract the final converged chain (last N frames) ---
    images = all_frames[-n_chain:]
    n_atoms = len(images[0])

    # Stage statistics
    if n_total % n_chain == 0:
        n_steps = n_total // n_chain
        console.print(f"  NEB steps completed: {n_steps} ({stage_label})")

    # --- Energies & forces ---
    energies = np.zeros(n_chain)
    fmaxes = np.zeros(n_chain)
    for i, img in enumerate(images):
        try:
            energies[i] = img.get_potential_energy()
            forces = img.get_forces()
            fmaxes[i] = np.max(np.abs(forces))
        except Exception:
            energies[i] = np.nan
            fmaxes[i] = np.nan

    if np.isnan(energies[0]):
        console.print("[red]Energies not available — images have no calculator results.[/red]")
        return

    e_rel = energies - energies[0]
    barrier = np.max(e_rel)
    saddle_idx = int(np.argmax(e_rel))
    saddle_fmax = fmaxes[saddle_idx]

    console.print()
    console.print("  [bold]Converged NEB chain:[/bold]")
    console.print(f"  {'Image':>6s}  {'E (eV)':>14s}  {'E rel (eV)':>12s}  {'fmax (eV/A)':>13s}")
    console.print(f"  {'-'*6}  {'-'*14}  {'-'*12}  {'-'*13}")
    for i in range(n_chain):
        marker = " ← saddle" if i == saddle_idx else ""
        console.print(
            f"  {i:>6d}  {energies[i]:>14.6f}  {e_rel[i]:>12.6f}  {fmaxes[i]:>13.6f}{marker}"
        )
    console.print()
    console.print(f"  Initial energy: {energies[0]:.6f} eV")
    console.print(f"  Final energy:   {energies[-1]:.6f} eV")
    console.print(f"  ΔE (final-init): {energies[-1] - energies[0]:.4f} eV")
    console.print(f"  Energy barrier: {barrier:.4f} eV  |  atoms/image: {n_atoms}")
    console.print(f"  Saddle point:   image {saddle_idx}  |  fmax = {saddle_fmax:.4f} eV/A")

    # --- Cubic spline ---
    x_raw = np.arange(n_chain)
    n_img = n_chain
    try:
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(x_raw, e_rel, bc_type="natural")
        x_fine = np.linspace(0, n_img - 1, (n_img - 1) * 50 + 1)
        e_fine = cs(x_fine)
        spline_max_idx = np.argmax(e_fine)
        spline_barrier = e_fine[spline_max_idx]
        spline_saddle_x = x_fine[spline_max_idx]
        has_spline = True
    except Exception:
        has_spline = False
        spline_barrier = barrier
        spline_saddle_x = float(saddle_idx)

    if has_spline:
        console.print(f"  Spline barrier: {spline_barrier:.4f} eV at image {spline_saddle_x:.2f}")
    console.print()

    # --- Plot (energy profile + force bar chart) ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 8), gridspec_kw={"height_ratios": [3, 1]}
    )

    # Upper panel: energy barrier
    ax1.plot(x_raw, e_rel, "o", color="#1f77b4", markersize=8, zorder=5,
             label="ASE NEB images")
    if has_spline:
        ax1.plot(x_fine, e_fine, "-", color="#1f77b4", linewidth=1.2, alpha=0.6,
                 label="Cubic spline")
        ax1.axvline(x=spline_saddle_x, color="#d62728", linestyle="--",
                    linewidth=0.8, label=f"Barrier: {spline_barrier:.3f} eV")
    else:
        ax1.plot(x_raw, e_rel, "-", color="#1f77b4", linewidth=1.2)
        ax1.axvline(x=saddle_idx, color="#d62728", linestyle="--", linewidth=0.8)
    ax1.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
    ax1.set_ylabel("Relative energy (eV)")
    ax1.set_title(
        f"ASE NEB Barrier — {spline_barrier:.3f} eV" if has_spline
        else f"ASE NEB Barrier — {barrier:.3f} eV"
    )
    ax1.legend(loc="upper left", fontsize=9)

    # Lower panel: max-force per image
    ax2.bar(x_raw, fmaxes, color="#ff7f0e", alpha=0.7, label="|F| max")
    ax2.axhline(y=0.05, color="#2ca02c", linestyle="--", linewidth=0.8,
                label="target fmax (0.05 eV/A)")
    ax2.set_xlabel("NEB image index")
    ax2.set_ylabel("Max force (eV/A)")
    ax2.legend(loc="upper right", fontsize=9)

    for ax in (ax1, ax2):
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)

    out_png = "ase_neb_barrier.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- Save data ---
    out_dat = "ase_neb_barrier.dat"
    with open(out_dat, "w") as f:
        f.write("# image  E(eV)  E_rel(eV)  E_spline(eV)  fmax(eV/A)\n")
        for i in range(n_chain):
            e_s = float(cs(x_raw[i])) if has_spline else e_rel[i]
            f.write(f"{i}  {energies[i]:.8f}  {e_rel[i]:.8f}  {e_s:.8f}  {fmaxes[i]:.8f}\n")

    console.print(f"[green]✓ Energy profile: {out_png}[/green]")
    console.print(f"[green]✓ Data file: {out_dat}[/green]")
    console.print(f"  Barrier: {spline_barrier:.4f} eV (spline)" if has_spline
                  else f"  Barrier: {barrier:.4f} eV")
    console.print(f"  Saddle:  image {saddle_idx}/{n_img - 1}" +
                  (f" (spline: {spline_saddle_x:.2f})" if has_spline else ""))
    console.print()
