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
    # Auto-detect: STRU_ini/STRU_fin or POSCAR_ini/POSCAR_fin
    init_path = "STRU_ini" if Path("STRU_ini").exists() else \
                ("POSCAR_ini" if Path("POSCAR_ini").exists() else "POSCAR_init")
    final_path = "STRU_fin" if Path("STRU_fin").exists() else \
                 ("POSCAR_fin" if Path("POSCAR_fin").exists() else "POSCAR_final")

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
    console.print(f"  Max atomic displacement: {d_max:.4f} Å → suggested {suggested} images (d_max / 0.8)")

    if interactive:
        n_images = int(_prompt(console, "Number of intermediate images", str(suggested)))
    else:
        n_images = suggested

    console.print(f"  Images: {n_images} intermediate + 2 endpoints")
    console.print("  Method: linear Cartesian interpolation")

    images = _linear_interpolate(atoms_init, atoms_final, n_images)

    fmt = "STRU"
    if interactive:
        fmt = _prompt_choice(console, "Output format",
                             ["STRU (ABACUS)", "POSCAR (VASP)"],
                             "STRU (ABACUS)")
        fmt = "STRU" if "STRU" in fmt else "POSCAR"

    _write_images(images, fmt)

    # Also write combined multi-frame view files (all atoms stacked in one cell)
    _write_chain_structure(images, atoms_init.cell)

    from ase.io import write as ase_write
    traj_path = f"path_{n_images + 2}frames.traj"
    ase_write(traj_path, images)
    console.print()
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

    # Auto-detect: STRU_ini/STRU_fin or POSCAR_ini/POSCAR_fin
    init_path = "STRU_ini" if Path("STRU_ini").exists() else \
                ("POSCAR_ini" if Path("POSCAR_ini").exists() else "POSCAR_init")
    final_path = "STRU_fin" if Path("STRU_fin").exists() else \
                 ("POSCAR_fin" if Path("POSCAR_fin").exists() else "POSCAR_final")

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
    console.print(f"  Max atomic displacement: {d_max:.4f} Å → suggested {suggested} images (d_max / 0.8)")

    if interactive:
        n_images = int(_prompt(console, "Number of intermediate images", str(suggested)))
    else:
        n_images = suggested

    console.print(f"  Images: {n_images} intermediate + 2 endpoints")
    console.print("  Method: IDPP (Image-Dependent Pair Potential)")

    images = _idpp_interpolate(atoms_init, atoms_final, n_images)

    fmt = "STRU"
    if interactive:
        fmt = _prompt_choice(console, "Output format",
                             ["STRU (ABACUS)", "POSCAR (VASP)"],
                             "STRU (ABACUS)")
        fmt = "STRU" if "STRU" in fmt else "POSCAR"

    _write_images(images, fmt)

    _write_chain_structure(images, atoms_init.cell)

    from ase.io import write as ase_write
    traj_path = f"path_{n_images + 2}frames.traj"
    ase_write(traj_path, images)
    console.print()
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
    if len(image_dirs) < 3:
        console.print("[red]No NEB image directories found (00/, 01/, ...).[/red]")
        console.print("[dim]Run task 1601 or 1602 first to generate NEB paths.[/dim]")
        return

    console.print(f"  Found {len(image_dirs)} image directories: "
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
        console.print(f"[red]No STRU or POSCAR found in {first_dir}/[/red]")
        return

    console.print(f"  Species: {', '.join(species)}")

    # Read cell for KPT auto-calculation
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
    if stru_file.exists():
        pp_map = dict(structure.pseudo_files)
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
        orb_file = orb_map.get(sp, f"{sp}.orb")
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
            orb = orb_map.get(sp, f"{sp}.orb")
            if orb and Path(orb).exists() and not (d / orb).exists():
                shutil.copy2(orb, d / orb)
    console.print(f"  [green]Files distributed to {len(image_dirs)} image dirs[/green]")

    abacus_bin = paths_cfg.get("abacus_binary", "abacus")

    # Interactive prompts
    use_kpt = False
    kspacing = "0.14"
    kpt_grid = [1, 1, 1]
    n_cores = "8"
    climb = True
    do_two_stage = True

    if interactive:
        n_cores = _prompt(console, "MPI cores per image", "8")
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

    # Use existing .traj if available
    traj_files = sorted(Path(".").glob("path_*frames.traj"))
    init_chain = str(traj_files[0]) if traj_files else [str(d / "POSCAR") for d in image_dirs]
    if traj_files:
        console.print(f"  [dim]Using: {init_chain}[/dim]")

    # Build YAML config
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
    calc_section["parallel"] = True
    calc_section["optimizer"] = "FIRE"
    calc_section["max_steps"] = 200

    abacus_params = {
        "calculation": "scf",
        "ecutwfc": 100,
        "basis_type": "lcao",
        "ks_solver": "genelpa",
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
        "out_wfc_lcao": 0,
        "out_bandgap": 0,
        "pseudo_dir": "./",
        "orbital_dir": "./",
        "pseudopotentials": pp_map,
        "basissets": orb_map,
    })

    abacus_section = {
        "command": abacus_bin,
        "mpi": int(n_cores),
        "omp": 1,
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
    console.print(f"  Images: {len(image_dirs)} (parallel=True, MPI={n_cores})")
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
    console.print("    2. Create conda env and install atst-tools:")
    console.print("       conda create -n atst-dev python=3.12 -y")
    console.print("       conda activate atst-dev")
    console.print('       pip install "atst-tools[parallel]==2.1.2"')
    console.print("    3. atst run neb.yaml")
    console.print("  [dim]Docs: https://github.com/QuantumMisaka/atst-tools/blob/main/docs/user/USER_GUIDE_CN.md[/dim]")
    console.print("  [dim]Or use task 1604 for a self-contained ASE script (no atst-tools needed)[/dim]")
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
    if len(image_dirs) < 3:
        console.print("[red]No NEB image directories found.[/red]")
        return

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

    if interactive:
        console.print(f"  [dim]ABACUS source (for PYTHONPATH): {abacus_src or 'not set'}[/dim]")
        src_in = console.input("  ABACUS source path [Enter=skip]: ").strip()
        if src_in:
            abacus_src = src_in
            paths_cfg["abacus_source"] = src_in
            from abacuscopilot.config import save_config
            save_config(config)

    # Build PP/orb maps (same logic as 1603)
    pp_map = {}
    orb_map = {}
    if stru_file.exists():
        pp_map = dict(structure.pseudo_files)
        orb_map = dict(structure.orbital_files)
    for sp in species:
        if sp not in pp_map or not pp_map[sp]:
            upfs = sorted(Path(".").glob(f"{sp}_*.upf"))
            pp_map[sp] = upfs[0].name if upfs else f"{sp}.upf"
        if sp not in orb_map or not orb_map[sp]:
            orbs = sorted(Path(".").glob(f"{sp}_*.orb"))
            orb_map[sp] = orbs[0].name if orbs else f"{sp}.orb"

    # Interactive config
    if interactive:
        n_cores = _prompt(console, "MPI cores per image", "8")
        kpt_mode = _prompt_choice(console, "K-point mode",
                                   ["kspacing (auto mesh)", "KPT (explicit grid)"],
                                   "kspacing (auto mesh)")
        use_kpt = "KPT" in kpt_mode
        if use_kpt:
            if stru_file.exists():
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
        n_cores = "8"
        use_kpt = False
        kspacing_val = 0.14
        kpt_grid = None
        climb = True
        do_two_stage = True

    # Use .traj if available
    traj_files = sorted(Path(".").glob("path_*frames.traj"))
    traj_path = str(traj_files[0]) if traj_files else None
    if traj_files:
        console.print(f"  [dim]Using: {traj_path}[/dim]")

    # Build the Python script
    pp_block = "    " + "\n    ".join(f'"{sp}": "{pp_map[sp]}",' for sp in species)
    orb_block = "    " + "\n    ".join(f'"{sp}": "{orb_map[sp]}",' for sp in species)

    kpt_lines = ""
    if use_kpt:
        kpt_lines = f'        "kpts": {kpt_grid},\n'
    else:
        kpt_lines = f'        "kspacing": {kspacing_val},\n'

    two_stage_block = ""
    if do_two_stage:
        two_stage_block = """neb = NEB(images, climb=True, k=0.1, parallel=True,
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
        two_stage_block = """neb = NEB(images, climb=True, k=0.1, parallel=True,
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
    print(f"  Suggested MPI cores: {n_cores} per image x {{n-2}} active = {{{n_cores}*(n-2)}} total")
    print("Validation passed — ready to run on server with ABACUS.")
    sys.exit(0)

from ase.mep import NEB
from ase.optimize import FIRE

# === Path setup ===
# Add ABACUS interface to Python path (auto-filled by AbacusCopilot)
_ABACUS_SRC = "{abacus_src}"
if _ABACUS_SRC and _ABACUS_SRC not in sys.path:
    _ase_iface = os.path.join(_ABACUS_SRC, "interfaces", "ASE_interface")
    if os.path.isdir(_ase_iface):
        sys.path.insert(0, _ase_iface)
    else:
        print(f"Warning: {{_ase_iface}} not found — adjust ABACUS source path")

from abacuslite import Abacus, AbacusProfile

# === Configuration ===
N_CORES = {n_cores}
TRAJ_FILE = "{traj_path}"
PSEUDO_DIR = "{pseudo_lib or './'}"
ORBITAL_DIR = "{orbital_lib or './'}"

profile = AbacusProfile(
    command="{mpirun} -np {n_cores} {abacus_bin}",
    omp_num_threads=1,
    pseudo_dir=PSEUDO_DIR,
    orbital_dir=ORBITAL_DIR,
)

pseudopotentials = {{
{pp_block}
}}

basissets = {{
{orb_block}
}}

inp_params = {{
    "calculation": "scf",
    "ecutwfc": 100,
    "basis_type": "lcao",
    "ks_solver": "genelpa",
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
    "out_wfc_lcao": 0,
    "out_bandgap": 0,
{kpt_lines}}}

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
        basissets=basissets,
        inp=inp_params,
    )
    print(f"  Image {{i}}: directory=neb-{{i}}")

# === Run NEB ===
print(f"Running NEB ({{n_images}} images, parallel=True)...")
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
    n_active = len(image_dirs) - 2  # exclude endpoints
    total_cores = int(n_cores) * n_active
    slurm_script = f'''#!/bin/bash
#SBATCH --job-name=neb
#SBATCH --nodes=1
#SBATCH --ntasks={n_active}
#SBATCH --cpus-per-task={n_cores}
#SBATCH --output=neb_%j.out
#SBATCH --error=neb_%j.err

source ~/softwares/abacus-develop-LTSv3.10.0/toolchain/abacus_env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Add ABACUS ASE interface to Python path
export PYTHONPATH="{abacus_src}/interfaces/ASE_interface:$PYTHONPATH"

python neb_run.py
'''
    slurm_path = "neb_slurm.sh"
    with open(slurm_path, "w") as f:
        f.write(slurm_script)
    import os
    os.chmod(slurm_path, 0o755)

    console.print(f"[green]ASE NEB script written: {out_path}[/green]")
    console.print(f"[green]SLURM script written: {slurm_path}[/green]")
    console.print(f"  Images: {len(image_dirs)} ({n_active} active, parallel=True)")
    console.print(f"  Cores: {n_cores}/image × {n_active} = {total_cores} total")
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
