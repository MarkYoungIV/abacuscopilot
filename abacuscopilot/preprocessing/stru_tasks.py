"""Structure file generation tasks for ABACUS.

Task IDs 201-299

STRU file generation from various input formats, structure conversion,
and basic structure editing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.models import Structure
from abacuscopilot.tasks import task


def _write_stru_bare(structure: Structure, is_lcao: bool = False,
                     filepath: str = "STRU", is_dp: bool = False) -> None:
    """Write STRU file, resolving pseudopotential/orbital filenames from library.

    Looks up the configured pseudo_library and orbital_library directories
    to find actual filenames (e.g. 'Li_ONCV_PBE-1.0.upf' instead of 'Li.upf').
    Falls back to placeholder filenames if libraries are not configured.

    Args:
        structure: Structure to write.
        is_lcao: Whether LCAO orbital section is needed (ignored if is_dp).
        filepath: Output path (default: "STRU").
        is_dp: If True, skip pseudo/orbital filenames (Deep Potential mode).
    """
    from abacuscopilot.config import load_config
    from abacuscopilot.io.stru_file import write_stru
    from abacuscopilot.preprocessing.system_tasks import _find_file_for_element

    if not is_dp:
        config = load_config()
        libraries = config.get("libraries", {})
        pseudo_lib = libraries.get("pseudo_library", "")
        orbital_lib = libraries.get("orbital_library", "")

        for sp in structure.species_order:
            # Always resolve from library; overwrites bare filenames like "S.upf"
            actual = _find_file_for_element(pseudo_lib, sp, ".upf") if pseudo_lib else None
            structure.pseudo_files[sp] = actual if actual else structure.pseudo_files.get(sp, f"{sp}.upf")

            if is_lcao:
                actual = _find_file_for_element(orbital_lib, sp, ".orb") if orbital_lib else None
                structure.orbital_files[sp] = actual if actual else structure.orbital_files.get(sp, f"{sp}.orb")

    write_stru(structure, filepath=filepath, is_lcao=is_lcao, is_dp=is_dp)


def _offer_3d_view(console, source_path: str = "STRU") -> None:
    """Offer to open the structure in ASE's 3D GUI viewer.

    If a POSCAR/CIF exists, use it directly.  Otherwise convert STRU to
    a temp file (cleaned up when the viewer closes).
    """
    try:
        import ase
        del ase
    except ImportError:
        return

    want = _prompt_choice(console, "Open in 3D viewer? (ase gui)",
                          ["Yes", "No"], "Yes")
    if "No" in want:
        return

    # Prefer original POSCAR/CIF if present; otherwise convert STRU
    import shutil
    import tempfile
    tmpdir = None
    viewer_path = None

    for f in ("POSCAR", "CONTCAR"):
        if Path(f).exists():
            viewer_path = f
            break
    if viewer_path is None:
        for f in sorted(Path(".").glob("*.cif")):
            viewer_path = str(f)
            break

    if viewer_path is None:
        try:
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru(source_path)
            atoms = structure.to_ase()
            tmpdir = tempfile.mkdtemp(prefix="abacuscopilot_")
            viewer_path = f"{tmpdir}/POSCAR"
            from ase.io import write as ase_write
            ase_write(viewer_path, atoms, format="vasp")
        except Exception as e:
            console.print(f"  [yellow]! Cannot preview: {e}[/yellow]")
            return

    console.print(f"  [dim]Launching ase gui {viewer_path} ...[/dim]")
    console.print("  [dim](close the 3D window to return to AbacusCopilot)[/dim]")
    import subprocess
    proc = subprocess.Popen(["ase", "gui", viewer_path])
    if tmpdir:
        import threading
        threading.Thread(
            target=lambda: (proc.wait(), shutil.rmtree(tmpdir, ignore_errors=True)),
            daemon=True).start()


def _choose_coordinate_type(console, structure: Structure) -> Structure:
    """Ask the user which coordinate system to use and convert if needed.

    ASE always returns fractional (Direct) coordinates.  Some users prefer
    Cartesian (Å) for readability.  This helper prompts for the choice and
    converts positions in-place when Cartesian is selected.
    """
    console.print()
    choice = _prompt_choice(
        console,
        "Coordinate type for STRU",
        ["Direct (fractional)", "Cartesian (Å)"],
        "Direct (fractional)",
    )
    if "Cartesian" not in choice:
        return structure  # already Direct

    # Direct → Cartesian_angstrom
    cell_ang = structure.lattice.cell_angstrom
    for atom in structure.atoms:
        atom.position = atom.position @ cell_ang
    structure.coordinate_type = "Cartesian_angstrom"
    return structure


def _run_full_calculation_setup(
    console, structure, cif_path: str | None = None
) -> None:
    """After STRU is written, proceed to INPUT generation + file preparation.

    This is the "configure full calculation" branch of STRU tasks.
    It delegates to input_tasks for the INPUT generation logic.
    """
    from abacuscopilot.config import load_config
    from abacuscopilot.core.models import InputParams
    from abacuscopilot.preprocessing.input_tasks import (
        _apply_template,
        _ask_basis_and_calc,
        _ask_lcao_solver,
        _get_template,
    )
    from abacuscopilot.preprocessing.system_tasks import (
        prepare_calculation_files,
        read_species_from_stru,
    )

    console.print()
    console.print("[bold cyan]--- Configure Calculation ---[/bold cyan]")
    console.print()

    # Step 1: pick basis + calculation, apply INPUT template
    basis, calc = _ask_basis_and_calc(console, [
        ("cell-relax (atoms + cell)", "cell-relax"),
        ("relax (atoms only)", "relax"),
        ("SCF", "scf"),
        ("Band (NSCF)", "nscf"),
        ("DOS (NSCF)", "dos"),
        ("MD", "md"),
    ], default_calc_index=0)

    params = InputParams()
    params.suffix = "ABACUS"
    template = _get_template(basis, calc)
    if template:
        _apply_template(params, template)

    # Ask CPU/GPU so LCAO gets the correct ks_solver (genelpa/cusolver).
    _ask_lcao_solver(console, params)

    # Override specific params for MD
    if calc == "md":
        params.md_type = _prompt_choice(console, "MD ensemble",
                                        ["nvt", "npt", "nve", "langevin", "fire", "msst"], "nvt")
        params.md_nstep = int(_prompt(console, "Number of MD steps", 10000))
        params.md_dt = float(_prompt(console, "Time step (fs)", 1.0))
        params.md_tfirst = float(_prompt(console, "Initial temperature (K)", 300.0))
        params.md_tlast = float(_prompt(console, "Final temperature (K)", 300.0))
    elif calc == "dos":
        dos_type = _prompt_choice(console, "DOS type",
                                  ["Total DOS", "Projected DOS (PDOS)"], "Total DOS")
        params.out_dos = 2 if "Projected" in dos_type else 1

    # Write INPUT
    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ INPUT file written.[/green]")
    console.print(f"  Calculation: {params.calculation}, Basis: {params.basis_type}")

    # Step 2: fix STRU with resolved upf/orb filenames
    from abacuscopilot.io.stru_file import read_stru as _read_stru
    _stru = _read_stru("STRU")
    _write_stru_bare(_stru, is_lcao=params.basis_type == "lcao", filepath="STRU")

    # Step 3: copy pseudopotential & orbital files
    config = load_config()
    libraries = config.get("libraries", {})
    pseudo_lib = libraries.get("pseudo_library", "")
    orbital_lib = libraries.get("orbital_library", "")

    if pseudo_lib:
        species = read_species_from_stru("STRU")
        if species:
            console.print()
            console.print("[bold]Copying pseudopotential / orbital files...[/bold]")
            result = prepare_calculation_files(
                species, params.basis_type, pseudo_lib, orbital_lib, ".", dry_run=False
            )
            if result["pseudo_files"]:
                console.print("  [green]✓ PP:[/green] " + ", ".join(result["pseudo_files"]))
            if result["orbital_files"]:
                console.print("  [green]✓ Orb:[/green] " + ", ".join(result["orbital_files"]))
            if result["errors"]:
                for e in result["errors"]:
                    console.print(f"  [yellow]![/yellow] {e}")

    console.print()
    console.print("[bold green]✓ Full calculation setup complete.[/bold green]")
    console.print(f"  Files: STRU + INPUT + {params.basis_type} support files")
    console.print()


# =============================================================================
# Task 201: CIF to STRU
# =============================================================================

@task(201, category="STRU", name="CIF to STRU",
      description="Convert CIF (Crystallographic Information File) to ABACUS STRU format")
def task_stru_from_cif(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a STRU file from a CIF file.

    Offers two modes:
      1. Just convert structure — STRU only, no pseudopotential/orbital setup
      2. Configure full calculation — STRU + INPUT + auto-copy PP/orb files
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== CIF to STRU ===[/bold cyan]")
    console.print()

    # --- find CIF file ---
    cif_path = None
    if args:
        for arg in args:
            if arg.endswith(".cif") and Path(arg).exists():
                cif_path = arg
                break

    if cif_path is None and interactive:
        cif_path = _prompt(console, "Path to CIF file")
        if not cif_path or not Path(cif_path).exists():
            console.print(f"[red]File not found: {cif_path}[/red]")
            return

    # --- read CIF ---
    try:
        from ase.io import read as ase_read
    except ImportError:
        console.print("[red]This task requires ASE. Install with: pip install ase[/red]")
        return

    try:
        atoms = ase_read(cif_path, format="cif")
        console.print(f"  [dim]Loaded {len(atoms)} atoms from CIF[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read CIF file: {e}[/red]")
        return

    structure = Structure.from_ase(atoms)

    # --- coordinate type ---
    if interactive:
        _choose_coordinate_type(console, structure)

    # --- choose mode ---
    if interactive:
        console.print()
        mode = _prompt_choice(
            console,
            "What would you like to do?",
            [
                "Just convert structure",
                "Configure full calculation (STRU + INPUT + files)",
            ],
            "Just convert structure",
        )
    else:
        mode = "Just convert structure"

    # --- resolve upf/orb info (标准规范: basis_type decides orbital section) ---
    if "full" in mode.lower() or "configure" in mode.lower():
        write_orb = True  # full setup re-writes STRU with the chosen basis anyway
    else:
        from abacuscopilot.core.standards import is_lcao_basis
        from abacuscopilot.preprocessing.system_tasks import resolve_basis_type
        write_orb = is_lcao_basis(resolve_basis_type(structure, interactive))

    _write_stru_bare(structure, is_lcao=write_orb)

    console.print()
    console.print("[green]✓ STRU file written successfully.[/green]")
    console.print(f"  {structure.num_atoms} atoms, {structure.num_species} species")
    console.print(f"  Lattice constant: {structure.lattice.constant:.6f} Bohr")
    if write_orb:
        console.print("  [dim]upf/orb filenames resolved from library[/dim]")

    # 3D viewer moved to task 206

    # --- full calculation setup ---
    if "full" in mode.lower() or "configure" in mode.lower():
        _run_full_calculation_setup(console, structure, str(cif_path))


# =============================================================================
# Task 202: POSCAR to STRU (VASP format)
# =============================================================================

@task(202, category="STRU", name="POSCAR to STRU",
      description="Convert VASP POSCAR/CONTCAR to ABACUS STRU format")
def task_stru_from_poscar(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert a VASP POSCAR file to ABACUS STRU format.

    Same two-mode workflow as task 201.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Convert POSCAR to STRU ===[/bold cyan]")
    console.print()

    poscar_path = "POSCAR"
    if args:
        for arg in args:
            if Path(arg).exists():
                poscar_path = arg
                break
    if interactive:
        poscar_path = _prompt(console, "Path to POSCAR/VASP structure file", poscar_path)

    if not Path(poscar_path).exists():
        console.print(f"[red]File not found: {poscar_path}[/red]")
        return

    try:
        from ase.io import read as ase_read
    except ImportError:
        console.print("[red]This task requires ASE. Install with: pip install ase[/red]")
        return

    try:
        atoms = ase_read(poscar_path, format="vasp")
        console.print(f"  [dim]Loaded {len(atoms)} atoms from POSCAR[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read POSCAR: {e}[/red]")
        if Path(poscar_path).exists():
            console.print(f"[dim]'{poscar_path}' doesn't look like a VASP POSCAR/CONTCAR.[/dim]")
            console.print("[dim]For CIF files, use task 201 instead.[/dim]")
        return

    structure = Structure.from_ase(atoms)

    # --- coordinate type ---
    if interactive:
        _choose_coordinate_type(console, structure)

    if interactive:
        console.print()
        mode = _prompt_choice(
            console,
            "What would you like to do?",
            [
                "Just convert structure",
                "Configure full calculation (STRU + INPUT + files)",
            ],
            "Just convert structure",
        )
    else:
        mode = "Just convert structure"

    # --- resolve upf/orb info (标准规范: basis_type decides orbital section) ---
    if "full" in mode.lower() or "configure" in mode.lower():
        write_orb = True  # full setup re-writes STRU with the chosen basis anyway
    else:
        from abacuscopilot.core.standards import is_lcao_basis
        from abacuscopilot.preprocessing.system_tasks import resolve_basis_type
        write_orb = is_lcao_basis(resolve_basis_type(structure, interactive))

    _write_stru_bare(structure, is_lcao=write_orb)

    console.print()
    console.print("[green]✓ STRU file written successfully.[/green]")
    console.print(f"  {structure.num_atoms} atoms, {structure.num_species} species")
    if write_orb:
        console.print("  [dim]upf/orb filenames resolved from library[/dim]")

    # 3D viewer moved to task 206

    if "full" in mode.lower() or "configure" in mode.lower():
        _run_full_calculation_setup(console, structure)


# =============================================================================
# Task 203: Coordinate conversion (Direct ↔ Cartesian)
# =============================================================================


def _is_direct(coord_type: str) -> bool:
    """Return True if *coord_type* represents fractional (Direct) coordinates."""
    return coord_type.lower().startswith("direct")


@task(203, category="STRU", name="Coord Convert",
      description="Convert STRU atomic positions between Direct (fractional) and Cartesian (Å)")
def task_coord_convert(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert STRU coordinates between Direct and Cartesian_angstrom.

    Reads the STRU, detects the current coordinate type, and converts to
    the other type.  Output is written to a separate file — the original
    STRU is never overwritten.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Convert STRU Coordinates ===[/bold cyan]")
    console.print()

    # Determine input STRU path
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break
    if interactive:
        inp = _prompt(console, "Input STRU file", stru_path)
        if inp:
            stru_path = inp

    if not Path(stru_path).exists():
        console.print(f"[red]File not found: {stru_path}[/red]")
        return

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                      f"{structure.num_species} species from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    current = structure.coordinate_type
    console.print(f"  Current coordinate type: [bold]{current}[/bold]")
    console.print()

    # Determine target
    if _is_direct(current):
        target = "Cartesian_angstrom"
        target_label = "Cartesian (Å)"
    else:
        target = "Direct"
        target_label = "Direct (fractional)"

    if interactive:
        choice = _prompt_choice(
            console,
            "Convert to",
            [target_label, current + " (keep as-is)"],
            target_label,
        )
        if "keep" in choice.lower():
            console.print("[yellow]No conversion needed.[/yellow]")
            return

    # Convert
    cell = structure.lattice.cell_angstrom  # Å

    if target == "Direct":
        # Cartesian → Direct
        for atom in structure.atoms:
            if current == "Cartesian_angstrom":
                cart = atom.position
            else:
                # Cartesian_bohr / Cartesian_au → convert to Angstrom first
                from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
                cart = atom.position * BOHR_TO_ANGSTROM
            atom.position = np.linalg.solve(cell.T, cart)
        structure.coordinate_type = "Direct"

    else:
        # Direct → Cartesian_angstrom
        for atom in structure.atoms:
            atom.position = atom.position @ cell
        structure.coordinate_type = "Cartesian_angstrom"

    # Determine output filename (never overwrite original)
    stem = Path(stru_path).stem
    if target == "Direct":
        out_path = f"{stem}_Direct" if stem != "STRU" else "STRU_Direct"
    else:
        out_path = f"{stem}_Cartesian" if stem != "STRU" else "STRU_Cartesian"

    from abacuscopilot.core.standards import is_lcao_basis
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    from abacuscopilot.preprocessing.system_tasks import resolve_basis_type
    _write_stru_bare(structure,
                     is_lcao=is_lcao_basis(resolve_basis_type(structure, interactive)),
                     filepath=out_path)

    console.print()
    console.print(f"[green]✓ Converted to {target} → {out_path}[/green]")
    console.print(f"   Original {stru_path} is unchanged.")
    console.print()


# =============================================================================
# Task 204: STRU to CIF
# =============================================================================

@task(204, category="STRU", name="STRU to CIF",
      description="Convert ABACUS STRU file to CIF (Crystallographic Information File)")
def task_stru_to_cif(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert a STRU file to CIF format."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Convert STRU to CIF ===[/bold cyan]")
    console.print()

    # Determine STRU path
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break
    if interactive:
        inp = _prompt(console, "Path to STRU file", stru_path)
        if inp:
            stru_path = inp

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                      f"{structure.num_species} species from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        if Path(stru_path).suffix.lower() in (".cif", ".vasp", ".poscar"):
            console.print(f"[dim]'{stru_path}' looks like a CIF/POSCAR file, not a STRU file.[/dim]")
            console.print("[dim]To convert CIF → STRU, use task 201. To convert POSCAR → STRU, use task 202.[/dim]")
        return

    # Convert via ASE
    try:
        atoms = structure.to_ase()
    except Exception as e:
        console.print(f"[red]ASE conversion failed: {e}[/red]")
        console.print("[dim]This task requires ASE. Install with: pip install ase[/dim]")
        return

    # Determine output path
    cif_path = Path(stru_path).with_suffix(".cif")
    if cif_path.name == stru_path or cif_path.suffix != ".cif":
        cif_path = Path(str(stru_path) + ".cif")
    out = _prompt(console, "Output CIF file", str(cif_path)) if interactive else str(cif_path)
    if not out:
        out = str(cif_path)

    try:
        from ase.io import write as ase_write
        ase_write(out, atoms, format="cif")
    except Exception as e:
        console.print(f"[red]Failed to write CIF: {e}[/red]")
        return

    console.print()
    console.print(f"[green]✓ CIF file written: {out}[/green]")
    console.print(f"  {len(atoms)} atoms, chemical formula: {atoms.get_chemical_formula()}")
    console.print()


# =============================================================================
# Task 205: STRU to POSCAR
# =============================================================================

@task(205, category="STRU", name="STRU to POSCAR",
      description="Convert ABACUS STRU file to VASP POSCAR format")
def task_stru_to_poscar(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert a STRU file to VASP POSCAR format."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Convert STRU to POSCAR ===[/bold cyan]")
    console.print()

    # Determine STRU path
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break
    if interactive:
        inp = _prompt(console, "Path to STRU file", stru_path)
        if inp:
            stru_path = inp

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                      f"{structure.num_species} species from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    # Convert via ASE
    try:
        atoms = structure.to_ase()
    except Exception as e:
        console.print(f"[red]ASE conversion failed: {e}[/red]")
        console.print("[dim]This task requires ASE. Install with: pip install ase[/dim]")
        return

    # Determine output path
    poscar_path = "POSCAR"
    out = _prompt(console, "Output POSCAR file", poscar_path) if interactive else poscar_path
    if not out:
        out = poscar_path

    try:
        from ase.io import write as ase_write
        ase_write(out, atoms, format="vasp", direct=True)
    except Exception as e:
        console.print(f"[red]Failed to write POSCAR: {e}[/red]")
        return

    console.print()
    console.print(f"[green]✓ POSCAR file written: {out}[/green]")
    console.print(f"  {len(atoms)} atoms, chemical formula: {atoms.get_chemical_formula()}")
    console.print("  [dim]Tip: use 'abacuscopilot -task 302' to generate KPT for the band path[/dim]")
    console.print()


# =============================================================================
# Task 206: STRU to PDB
# =============================================================================

@task(206, category="STRU", name="STRU to PDB",
      description="Convert ABACUS STRU to PDB format for VMD/PyMOL visualization")
def task_stru_to_pdb(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert STRU to PDB format."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== STRU to PDB ===[/bold cyan]")
    console.print()

    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break
    if interactive:
        inp = _prompt(console, "STRU file path", stru_path)
        if inp and Path(inp).exists():
            stru_path = inp

    if not Path(stru_path).exists():
        console.print(f"[red]File not found: {stru_path}[/red]")
        return

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        atoms = structure.to_ase()
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    out = Path(stru_path).with_suffix(".pdb")
    if interactive:
        console.print(f"  [dim]Default: {out}[/dim]")
        inp = console.input("  Output PDB file: ").strip()
        out = inp if inp else str(out)
    if out == str(stru_path):
        out = str(Path(stru_path).with_suffix(".pdb"))

    try:
        from ase.io import write as ase_write
        ase_write(str(out), atoms, format="proteindatabank")
    except Exception as e:
        console.print(f"[red]Failed to write PDB: {e}[/red]")
        return

    console.print()
    console.print(f"[green]✓ PDB written: {out}[/green]")
    console.print(f"  {len(atoms)} atoms, open with VMD/PyMOL")
    console.print()


# =============================================================================
# Task 207: STRU 3D visualization
# =============================================================================

@task(207, category="STRU", name="View Structure",
      description="Open STRU/CIF/POSCAR in ASE 3D viewer for visual inspection")
def task_view_structure(args: list[str] | None = None, interactive: bool = True) -> None:
    """Open a structure file in ASE's interactive 3D GUI viewer."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== 3D Structure Viewer ===[/bold cyan]")
    console.print()

    try:
        import ase
    except ImportError:
        console.print("[red]ASE is required. Install with: pip install ase pyqt5[/red]")
        return

    viewer_path = None
    if args:
        for arg in args:
            if Path(arg).exists():
                viewer_path = arg
                break

    if viewer_path is None and interactive:
        inp = _prompt(console, "Structure file path", "STRU")
        if inp and Path(inp).exists():
            viewer_path = inp

    if viewer_path is None:
        for f in ("STRU", "POSCAR", "CONTCAR"):
            if Path(f).exists():
                viewer_path = f
                break
        if viewer_path is None:
            cif = sorted(Path(".").glob("*.cif"))
            if cif:
                viewer_path = str(cif[0])
        if viewer_path is None:
            traj = sorted(Path(".").glob("*.traj"))
            if traj:
                viewer_path = str(traj[0])

    if viewer_path is None or not Path(viewer_path).exists():
        console.print("[red]No structure file found (STRU/POSCAR/CONTCAR/*.cif/*.traj).[/red]")
        return

    console.print(f"  [dim]Opening: {viewer_path}[/dim]")

    # STRU needs conversion (ASE can't read ABACUS format natively)
    import shutil
    import tempfile
    tmpdir = None
    open_path = viewer_path
    if (Path(viewer_path).suffix == "" or "STRU" in str(viewer_path)) and \
       Path(viewer_path).suffix != ".traj":
        try:
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru(viewer_path)
            atoms = structure.to_ase()
            tmpdir = tempfile.mkdtemp(prefix="abacuscopilot_")
            open_path = f"{tmpdir}/POSCAR"
            from ase.io import write as ase_write
            ase_write(open_path, atoms, format="vasp")
            console.print("  [dim](converted STRU -> POSCAR for ASE)[/dim]")
        except Exception as e:
            console.print(f"[red]Failed to read STRU: {e}[/red]")
            return

    console.print("  [dim](close the 3D window to return)[/dim]")
    import subprocess
    import threading
    proc = subprocess.Popen(["ase", "gui", open_path])
    if tmpdir:
        threading.Thread(
            target=lambda: (proc.wait(), shutil.rmtree(tmpdir, ignore_errors=True)),
            daemon=True).start()
    console.print()
