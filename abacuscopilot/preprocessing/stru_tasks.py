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


# =============================================================================
# Task 208: STRU to LAMMPS data file
# =============================================================================

@task(208, category="STRU", name="STRU to LAMMPS",
      description="Convert STRU to a LAMMPS data file (graph.lmp)")
def task_stru_to_lammps(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert STRU to LAMMPS data format.

    Outputs a LAMMPS `data` file (atom_style full) with masses, box
    dimensions, and Cartesian coordinates in Angstrom.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== STRU → LAMMPS Data File ===[/bold cyan]")
    console.print()

    # Load STRU
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                  f"{structure.num_species} species from {stru_path}[/dim]")

    # Convert to ASE for Cartesian Angstrom coordinates + cell
    atoms_ase = structure.to_ase()
    pos_cart = atoms_ase.get_positions()       # Cartesian Angstrom
    cell_ang = atoms_ase.get_cell()[:]          # 3×3 Angstrom

    # Output filename
    if interactive:
        out_name = _prompt(console, "Output filename", "graph.lmp")
    else:
        out_name = "graph.lmp"
    out_path = Path(out_name)

    # Species → type mapping (LAMMPS type IDs start at 1)
    species = structure.species_order
    type_map = {sp: i + 1 for i, sp in enumerate(species)}

    # Box dimensions (LAMMPS upper-triangular decomposition)
    a_len, b_len, c_len, alpha, beta, gamma = atoms_ase.get_cell_lengths_and_angles()
    import numpy as np
    alpha_r, beta_r, gamma_r = np.radians(alpha), np.radians(beta), np.radians(gamma)
    xhi = a_len
    xy = b_len * np.cos(gamma_r)
    xz = c_len * np.cos(beta_r)
    yhi = np.sqrt(b_len**2 - xy**2)
    yz = (b_len * c_len * np.cos(alpha_r) - xy * xz) / yhi if yhi > 0 else 0.0
    zhi = np.sqrt(c_len**2 - xz**2 - yz**2)
    is_triclinic = abs(xy) > 1e-6 or abs(xz) > 1e-6 or abs(yz) > 1e-6

    # Atomic masses (from ABACUS STRU standard values)
    from abacuscopilot.io.stru_file import _ATOMIC_MASSES

    with open(out_path, "w") as f:
        f.write(f"# LAMMPS data file — converted from {stru_path} by AbacusCopilot\n\n")
        f.write(f"{structure.num_atoms} atoms\n")
        f.write(f"{structure.num_species} atom types\n\n")

        if is_triclinic:
            f.write(f"  0.000000  {xhi:.6f}  xlo xhi\n")
            f.write(f"  0.000000  {yhi:.6f}  ylo yhi\n")
            f.write(f"  0.000000  {zhi:.6f}  zlo zhi\n")
            f.write(f"  {xy:.6f}  {xz:.6f}  {yz:.6f}  xy xz yz\n")
        else:
            f.write(f"  0.000000  {xhi:.6f}  xlo xhi\n")
            f.write(f"  0.000000  {yhi:.6f}  ylo yhi\n")
            f.write(f"  0.000000  {zhi:.6f}  zlo zhi\n")

        f.write(f"\nMasses\n\n")
        for sp in species:
            mass = _ATOMIC_MASSES.get(sp, 0.0)
            f.write(f"  {type_map[sp]}  {mass:.6f}  # {sp}\n")

        f.write(f"\nAtoms\n\n")
        for i, atom in enumerate(structure.atoms, 1):
            t = type_map[atom.species]
            x, y, z = pos_cart[i - 1]
            f.write(f"  {i}  {t}  {x:.6f}  {y:.6f}  {z:.6f}\n")

    console.print()
    console.print(f"[green]✓ LAMMPS data file written: {out_path.name}[/green]")
    console.print(f"  {structure.num_atoms} atoms, {structure.num_species} atom types")
    console.print(f"  Box: {xhi:.4f} × {yhi:.4f} × {zhi:.4f} Å"
                  + (" (triclinic)" if is_triclinic else ""))
    console.print()


# =============================================================================
# Task 209: LAMMPS data file → STRU
# =============================================================================

@task(209, category="STRU", name="LAMMPS to STRU",
      description="Convert a LAMMPS data file to ABACUS STRU format")
def task_lammps_to_stru(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert a LAMMPS data file to STRU.

    Detects atom types from masses (IUPAC standard).  If a mass is
    ambiguous (multiple elements within tolerance), the user is prompted.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== LAMMPS → STRU ===[/bold cyan]")
    console.print()

    # Find input file
    lmp_path = "graph.lmp"
    if args:
        for arg in args:
            if Path(arg).exists() and arg != "STRU":
                lmp_path = arg
                break
    if interactive:
        inp = console.input(f"  LAMMPS data file [{lmp_path}]: ").strip()
        if inp and Path(inp).exists():
            lmp_path = inp

    if not Path(lmp_path).exists():
        console.print(f"[red]File not found: {lmp_path}[/red]")
        return

    content = Path(lmp_path).read_text(errors="ignore")
    lines = content.split("\n")

    # --- Parse LAMMPS data file ---
    import re
    from abacuscopilot.core.models import Atom, Lattice, Structure

    n_atoms = n_types = 0
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0
    xy = xz = yz = 0.0
    masses: dict[int, float] = {}
    positions: list[tuple[int, float, float, float]] = []  # (type, x, y, z)
    section = None

    for line in content.split("\n"):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            # LAMMPS comment lines may contain the atom/type counts
            m = re.search(r"(\d+)\s+atoms", stripped)
            if m and n_atoms == 0:
                n_atoms = int(m.group(1))
            m = re.search(r"(\d+)\s+atom types", stripped)
            if m and n_types == 0:
                n_types = int(m.group(1))
            continue

        # Section headers (case-insensitive)
        low = stripped.lower()

        if low in ("masses",):
            section = "masses"
            continue
        elif low in ("atoms",):
            section = "atoms"
            continue
        elif low in ("velocities", "bonds", "angles", "dihedrals", "impropers",
                     "pair coeffs", "bond coeffs", "angle coeffs",
                     "pair_coeffs", "bond_coeffs", "angle_coeffs",
                     "atoms #", "atoms  #"):
            section = None  # skip to end of file
            continue

        # Box: xlo xhi [xy xz yz]
        if section is None and n_atoms == 0:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    vals = [float(x) for x in parts[:2]]
                    if xlo == xhi == 0:
                        xlo, xhi = vals[0], vals[1]
                        continue
                except ValueError:
                    pass

        if section is None and xlo != 0 and xhi != 0:
            # Try to read ylo yhi
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    vals = [float(x) for x in parts[:2]]
                    if ylo == yhi == 0:
                        ylo, yhi = vals[0], vals[1]
                        continue
                except ValueError:
                    pass

        if section is None and ylo != 0 and yhi != 0 and zlo == zhi == 0:
            # Try to read zlo zhi
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    vals = [float(x) for x in parts[:2]]
                    if zlo == zhi == 0:
                        zlo, zhi = vals[0], vals[1]
                        continue
                except ValueError:
                    pass

        if section is None and zhi != 0:
            # xy xz yz (triclinic)
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    xy, xz, yz = float(parts[0]), float(parts[1]), float(parts[2])
                    continue
                except ValueError:
                    pass

        # Masses section: type mass
        if section == "masses":
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    masses[int(parts[0])] = float(parts[1])
                except (ValueError, IndexError):
                    pass
            continue

        # Atoms section: id type x y z [more...]
        if section == "atoms":
            parts = stripped.split()
            if len(parts) >= 5:
                try:
                    a_type = int(parts[1])
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    positions.append((a_type, x, y, z))
                except (ValueError, IndexError):
                    pass
            continue

    # Deduce atom counts from data
    if not positions:
        console.print("[red]No atom positions found in the file.[/red]")
        return
    n_atoms = max(n_atoms, len(positions))
    console.print(f"  [dim]Parsed {len(positions)} atoms, {len(masses)} types[/dim]")

    # --- Match masses to elements ---
    from abacuscopilot.io.stru_file import _ATOMIC_MASSES
    # Build reverse lookup: element → mass (keep only most common isotope)
    elem_mass: dict[str, float] = {}
    for elem, mass in sorted(_ATOMIC_MASSES.items()):
        elem_mass[elem] = mass

    TOLERANCE = 0.1  # mass tolerance for element identification
    type_to_elem: dict[int, str] = {}

    for tid in sorted(masses):
        m = masses[tid]
        # Find closest element by mass
        candidates = []
        for elem, ref_mass in elem_mass.items():
            if abs(m - ref_mass) <= TOLERANCE:
                candidates.append((abs(m - ref_mass), elem))

        if not candidates:
            # Broader search
            candidates = sorted(
                [(abs(m - ref_mass), elem) for elem, ref_mass in elem_mass.items()]
            )[:3]

        candidates.sort()
        if len(candidates) == 1 and candidates[0][0] <= 0.01:
            # Unique exact match
            type_to_elem[tid] = candidates[0][1]
            console.print(f"  Type {tid} (mass {m:.4f}) → [green]{candidates[0][1]}[/green]")
        elif interactive:
            console.print(f"\n  Type {tid} (mass {m:.4f}):")
            for i, (diff, elem) in enumerate(candidates[:5], 1):
                marker = " ← best match" if i == 1 else ""
                console.print(f"    {i}. {elem} (mass {elem_mass[elem]:.4f}, Δ={diff:.4f}){marker}")
            choice = console.input(
                f"  Enter element symbol (or 1-{min(len(candidates), 5)}): "
            ).strip()
            if choice.isdigit() and 1 <= int(choice) <= len(candidates):
                type_to_elem[tid] = candidates[int(choice) - 1][1]
            elif choice.upper() in elem_mass:
                type_to_elem[tid] = choice.upper()
            else:
                type_to_elem[tid] = candidates[0][1]  # best guess
        else:
            type_to_elem[tid] = candidates[0][1]  # non-interactive: best guess
            console.print(f"  Type {tid} (mass {m:.4f}) → [yellow]{candidates[0][1]}[/yellow] (best guess)")

    # --- Build Structure ---
    species = [type_to_elem[tid] for tid in sorted(type_to_elem)]
    species_order = list(dict.fromkeys(species))  # dedupe, preserve order

    # Lattice from box dimensions (Cartesian Angstrom)
    a = [xhi - xlo, 0.0, 0.0]
    b = [xy, yhi - ylo, 0.0]
    c = [xz, yz, zhi - zlo]
    lattice = Lattice()
    cell_ang = np.array([a, b, c], dtype=float)
    # STRU stores lattice vectors in Bohr; constant=1 keeps them as-is
    from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
    lattice.constant = 1.0
    lattice._vectors = cell_ang * ANGSTROM_TO_BOHR  # Angstrom → Bohr

    # Atoms in Cartesian Angstrom
    struct = Structure()
    struct.coordinate_type = "Cartesian_angstrom"
    struct.species_order = species_order
    struct.lattice = lattice

    # Assign species to each atom
    type_to_sp = {tid: type_to_elem[tid] for tid in sorted(type_to_elem)}
    for a_type, x, y, z in positions:
        sp = type_to_sp.get(a_type, "X")
        struct.atoms.append(Atom(species=sp, position=[x, y, z]))

    # Placeholder pseudo files
    for sp in species_order:
        struct.pseudo_files[sp] = f"{sp}.upf"

    # Write STRU
    out_path = "STRU_LAMMPS" if Path("STRU").exists() else "STRU"
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    _write_stru_bare(struct, is_lcao=False, filepath=out_path)

    console.print()
    console.print(f"[green]✓ STRU written: {out_path}[/green]")
    console.print(f"  {struct.num_atoms} atoms, {struct.num_species} species: {' '.join(species_order)}")
    if Path("STRU").exists() and out_path != "STRU":
        console.print("  [dim]Original STRU is unchanged.[/dim]")
    console.print()
