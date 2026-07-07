"""Structure editing tasks for abacuscopilot.

Task IDs 401-499

Supercell construction, slab/surface building, coordinate conversion,
atom fixing/unfixing, and other structure manipulation operations.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.models import Atom, Lattice, Structure
from abacuscopilot.tasks import task

# =============================================================================
# Task 401: Build supercell
# =============================================================================


def build_supercell(structure: Structure, nx: int, ny: int, nz: int) -> Structure:
    """Build an Nx×Ny×Nz supercell from a Structure.

    Replicates the lattice and atomic positions to construct a supercell.
    Fractional coordinates are scaled to the new supercell: new_pos = pos / (Nx,Ny,Nz).
    """
    nx, ny, nz = max(1, nx), max(1, ny), max(1, nz)

    new_structure = Structure()
    new_structure.coordinate_type = structure.coordinate_type
    new_structure.species_order = list(structure.species_order)
    new_structure.magnetism = dict(structure.magnetism)

    # Scale lattice vectors
    scale = np.diag([nx, ny, nz])
    new_structure.lattice = Lattice(
        constant=structure.lattice.constant,
        vectors=structure.lattice.vectors @ scale,
    )

    # Copy pseudopotential / orbital info
    new_structure.pseudo_files = dict(structure.pseudo_files)
    new_structure.orbital_files = dict(structure.orbital_files)

    # Replicate atoms
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                offset = np.array([ix, iy, iz], dtype=float)
                for atom in structure.atoms:
                    if structure.coordinate_type.startswith("Direct"):
                        new_pos = atom.position / np.array([nx, ny, nz]) + offset / np.array([nx, ny, nz])
                    else:
                        # Cartesian: translate by the original cell vector
                        cell = structure.lattice.cell
                        new_pos = atom.position + ix * cell[0] + iy * cell[1] + iz * cell[2]

                    new_structure.atoms.append(Atom(
                        species=atom.species,
                        position=new_pos,
                        fix=atom.fix,
                        magmom=atom.magmom,
                        velocity=atom.velocity,
                        angle1=atom.angle1,
                        angle2=atom.angle2,
                    ))

    return new_structure


@task(401, category="Structure Editing", name="Build Supercell",
      description="Build Nx×Ny×Nz supercell from current STRU")
def task_supercell(args: list[str] | None = None, interactive: bool = True) -> None:
    """Build a supercell from the current STRU file."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Build Supercell ===[/bold cyan]")
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
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                      f"{structure.num_species} species from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    console.print(f"  Original cell: {structure.lattice.cell_angstrom[0,0]:.4f} "
                  f"× {structure.lattice.cell_angstrom[1,1]:.4f} "
                  f"× {structure.lattice.cell_angstrom[2,2]:.4f} Å")

    if interactive:
        console.print()
        nx = int(_prompt(console, "Supercell Nx", "1"))
        ny = int(_prompt(console, "Supercell Ny", "1"))
        nz = int(_prompt(console, "Supercell Nz", "1"))
    else:
        nx, ny, nz = 1, 1, 1

    if nx == 1 and ny == 1 and nz == 1:
        console.print("[yellow]Nx=Ny=Nz=1 results in the same cell. Use N,N,N > 1 for a supercell.[/yellow]")
        return

    new_structure = build_supercell(structure, nx, ny, nz)

    # Build output filename: Supercell{Nx}{Ny}{Nz}.STRU (e.g. Supercell234.STRU)
    out_stru = f"Supercell{nx}{ny}{nz}.STRU"
    out_transform = f"Supercell{nx}{ny}{nz}.TRANSFORM"

    # Write supercell STRU
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    is_lcao = bool(structure.orbital_files)
    _write_stru_bare(new_structure, is_lcao=is_lcao, filepath=out_stru)

    # Write TRANSFORM matrix
    a, b, c = np.linalg.norm(structure.lattice.cell_angstrom, axis=1)
    transform_lines = [
        f"# Supercell transformation matrix: {nx}×{ny}×{nz}",
        f"# Original cell: a={a:.6f}  b={b:.6f}  c={c:.6f} Å",
        f"{nx:>4d} {0:>4d} {0:>4d}",
        f"{0:>4d} {ny:>4d} {0:>4d}",
        f"{0:>4d} {0:>4d} {nz:>4d}",
        "",
    ]
    with open(out_transform, "w") as f:
        f.write("\n".join(transform_lines))

    console.print()
    console.print(f"[green]✓ Supercell written to {out_stru} ({nx}×{ny}×{nz}).[/green]")
    console.print(f"   TRANSFORM matrix saved to {out_transform}")
    console.print("   Original STRU is unchanged.")
    console.print(f"  Atoms: {structure.num_atoms} → {new_structure.num_atoms}")
    cell_ang = new_structure.lattice.cell_angstrom
    console.print(f"  Cell: {cell_ang[0,0]:.4f} × {cell_ang[1,1]:.4f} × {cell_ang[2,2]:.4f} Å")
    console.print()


# =============================================================================
# Task 402: Redefine Lattice
# =============================================================================


def apply_lattice_transform(structure: Structure, M: np.ndarray) -> Structure:
    """Apply a 3×3 transformation matrix to redefine the lattice.

    New lattice vectors: L' = L @ M
    Fractional positions transform as: f' = M⁻¹ @ f

    When |det M| > 1 and M has integer entries, the new cell is larger than
    the old one — atoms are replicated (|det M| images per original atom) to
    fill the expanded cell.  This is the same principle as a supercell:
    the old cell maps to a sub-volume of the new cell; the remaining volume
    is filled by translating old atoms by integer multiples of the old
    lattice vectors and then mapping them into the new basis.
    """
    M_inv = np.linalg.inv(M)
    det_M = abs(np.linalg.det(M))
    det_int = int(round(det_M))

    # Check whether M is an integer matrix (within tolerance)
    is_integer = np.allclose(M, np.round(M), atol=1e-10)
    do_replicate = is_integer and det_int > 1

    new_structure = Structure()
    new_structure.coordinate_type = "Direct"  # always output Direct
    new_structure.species_order = list(structure.species_order)
    new_structure.magnetism = dict(structure.magnetism)
    new_structure.pseudo_files = dict(structure.pseudo_files)
    new_structure.orbital_files = dict(structure.orbital_files)

    # Transform lattice
    new_structure.lattice = Lattice(
        constant=structure.lattice.constant,
        vectors=structure.lattice.vectors @ M,
    )

    # --- Build tiling offsets (only needed when replicating atoms) ---
    g_offsets: list[np.ndarray] = []
    if do_replicate:
        # Find the |det M| distinct integer g vectors whose images under
        # M⁻¹ cover the new cell without overlap.
        col_bounds = np.sum(np.abs(np.round(M).astype(int)), axis=0)
        seen: set[tuple[float, float, float]] = set()
        for gx in range(col_bounds[0]):
            for gy in range(col_bounds[1]):
                for gz in range(col_bounds[2]):
                    g = np.array([gx, gy, gz], dtype=float)
                    offset_key = tuple((M_inv @ g).round(10) % 1.0)
                    if offset_key not in seen:
                        seen.add(offset_key)
                        g_offsets.append(g)
                        if len(g_offsets) >= det_int:
                            break
                if len(g_offsets) >= det_int:
                    break
            if len(g_offsets) >= det_int:
                break
        # Safety: if we didn't get exactly det_int offsets, fall back to
        # no replication (the user may have supplied a non-unimodular float
        # matrix for which the integer lattice interpretation breaks).
        if len(g_offsets) != det_int:
            do_replicate = False
            g_offsets = []

    # --- Transform atomic positions ---
    cell = structure.lattice.cell  # Bohr

    def _get_fractional(atom) -> np.ndarray:
        """Return the atom's fractional position in the OLD basis."""
        if structure.coordinate_type.startswith("Direct"):
            return atom.position.copy()
        # Cartesian → fractional via old cell
        if structure.coordinate_type == "Cartesian_angstrom":
            from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
            pos_bohr = atom.position * ANGSTROM_TO_BOHR
        else:
            pos_bohr = atom.position.copy()
        return np.linalg.solve(cell.T, pos_bohr)

    for atom in structure.atoms:
        f_old = _get_fractional(atom)

        if do_replicate:
            # Generate |det M| images
            for g in g_offsets:
                f_new = (M_inv @ (f_old + g)) % 1.0
                new_structure.atoms.append(Atom(
                    species=atom.species,
                    position=f_new,
                    fix=atom.fix,
                    magmom=atom.magmom,
                    velocity=atom.velocity,
                    angle1=atom.angle1,
                    angle2=atom.angle2,
                ))
        else:
            # Pure redefinition — just transform and wrap
            f_new = (M_inv @ f_old) % 1.0
            new_structure.atoms.append(Atom(
                species=atom.species,
                position=f_new,
                fix=atom.fix,
                magmom=atom.magmom,
                velocity=atom.velocity,
                angle1=atom.angle1,
                angle2=atom.angle2,
            ))

    return new_structure


@task(402, category="Structure Editing", name="Redefine Lattice",
      description="Apply a 3×3 transformation matrix to redefine the unit cell")
def task_redefine_lattice(args: list[str] | None = None, interactive: bool = True) -> None:
    """Redefine the lattice by applying a 3×3 transformation matrix.

    This is used to change unit-cell conventions (e.g. conventional ↔ primitive
    cell, or re-orientation) while preserving the physical crystal structure.

    The matrix M is applied as  new_lattice = old_lattice @ M.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Redefine Lattice ===[/bold cyan]")
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
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms, "
                      f"{structure.num_species} species from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    # Show current lattice
    cell = structure.lattice.cell_angstrom
    console.print("  Current lattice (Å):")
    console.print(f"    a = [{cell[0,0]:.6f}  {cell[0,1]:.6f}  {cell[0,2]:.6f}]")
    console.print(f"    b = [{cell[1,0]:.6f}  {cell[1,1]:.6f}  {cell[1,2]:.6f}]")
    console.print(f"    c = [{cell[2,0]:.6f}  {cell[2,1]:.6f}  {cell[2,2]:.6f}]")
    console.print(f"  Coordinate type: {structure.coordinate_type}")
    console.print()

    if interactive:
        console.print("  Enter the 3×3 transformation matrix M")
        console.print("  (new_lattice = old_lattice @ M):")
        console.print()
        r1 = [float(x) for x in _prompt(console, "  Row 1  [1 0 0]", "1 0 0").split()]
        r2 = [float(x) for x in _prompt(console, "  Row 2  [0 1 0]", "0 1 0").split()]
        r3 = [float(x) for x in _prompt(console, "  Row 3  [0 0 1]", "0 0 1").split()]

        if len(r1) != 3 or len(r2) != 3 or len(r3) != 3:
            console.print("[red]Each row must have exactly 3 numbers.[/red]")
            return
        M = np.array([r1, r2, r3])
    else:
        console.print("[yellow]Non-interactive mode: identity matrix (no change).[/yellow]")
        M = np.eye(3)

    det = np.linalg.det(M)
    if abs(det) < 1e-10:
        console.print("[red]Transformation matrix is singular — cannot invert.[/red]")
        return

    # Apply
    try:
        new_structure = apply_lattice_transform(structure, M)
    except Exception as e:
        console.print(f"[red]Transformation failed: {e}[/red]")
        return

    # Write
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    is_lcao = bool(structure.orbital_files)
    _write_stru_bare(new_structure, is_lcao=is_lcao, filepath="Redefined.STRU")

    # Write TRANSFORM
    a_old = np.linalg.norm(structure.lattice.cell_angstrom, axis=1)
    a_new = np.linalg.norm(new_structure.lattice.cell_angstrom, axis=1)
    transform_lines = [
        "# Lattice redefinition matrix",
        f"# Old cell: a={a_old[0]:.6f}  b={a_old[1]:.6f}  c={a_old[2]:.6f} Å",
        f"# New cell: a={a_new[0]:.6f}  b={a_new[1]:.6f}  c={a_new[2]:.6f} Å",
        f"# Volume ratio (det M): {det:.6f}",
        f"{M[0,0]:>10.6f} {M[0,1]:>10.6f} {M[0,2]:>10.6f}",
        f"{M[1,0]:>10.6f} {M[1,1]:>10.6f} {M[1,2]:>10.6f}",
        f"{M[2,0]:>10.6f} {M[2,1]:>10.6f} {M[2,2]:>10.6f}",
        "",
    ]
    with open("Redefined.TRANSFORM", "w") as f:
        f.write("\n".join(transform_lines))

    console.print()
    console.print("[green]✓ Redefined lattice written to Redefined.STRU[/green]")
    console.print("   TRANSFORM matrix saved to Redefined.TRANSFORM")
    console.print("   Original STRU is unchanged.")
    console.print(f"   Volume ratio (|det M|): {abs(det):.6f}")
    console.print(f"  Atoms: {new_structure.num_atoms}")
    console.print(f"  New cell: a={a_new[0]:.4f}  b={a_new[1]:.4f}  c={a_new[2]:.4f} Å")
    console.print()


# =============================================================================
# Task 403: Coordinate conversion (Direct ↔ Cartesian)
# =============================================================================


def _is_fractional(coord_type: str) -> bool:
    return coord_type.lower().startswith("direct")


@task(403, category="Structure Editing", name="Coord Conversion",
      description="Convert atomic positions between Direct (fractional) and Cartesian coordinates")
def task_coord_convert(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert between fractional and Cartesian coordinates in STRU."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Coordinate Conversion ===[/bold cyan]")
    console.print()

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

    current = structure.coordinate_type
    console.print(f"  Current coordinate type: [bold]{current}[/bold]")
    console.print(f"  Atoms: {structure.num_atoms}")

    if interactive:
        if _is_fractional(current):
            console.print()
            target = _prompt_choice(console, "Convert to",
                                    ["Cartesian (Å)", "Direct (keep as-is)"],
                                    "Cartesian (Å)")
            target_type = "Direct" if "Direct" in target else "Cartesian_angstrom"
        else:
            target_type = "Direct"
            console.print("  → Converting to Direct (fractional)")
    else:
        target_type = "Cartesian_angstrom" if not _is_fractional(current) else "Direct"

    if target_type == current:
        console.print("[yellow]Already in the requested coordinate type.[/yellow]")
        return

    # Perform conversion
    cell = structure.lattice.cell  # Bohr

    if _is_fractional(target_type):
        # Cartesian → Direct
        from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
        cart_positions = structure.positions * ANGSTROM_TO_BOHR
        frac = np.linalg.solve(cell.T, cart_positions.T).T
        structure.coordinate_type = "Direct"
        for i, atom in enumerate(structure.atoms):
            atom.position = frac[i]

    else:
        # Direct → Cartesian (Å)
        frac = structure.positions
        from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
        cart_final = (frac @ cell) * BOHR_TO_ANGSTROM
        structure.coordinate_type = "Cartesian_angstrom"
        for i, atom in enumerate(structure.atoms):
            atom.position = cart_final[i]

    # Write to new file (preserve original STRU)
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    if structure.coordinate_type.startswith("Direct"):
        out_path = "STRU_Direct"
    elif "bohr" in structure.coordinate_type.lower():
        out_path = "STRU_Cartesian_bohr"
    else:
        out_path = "STRU_Cartesian"
    _write_stru_bare(structure, is_lcao=bool(structure.orbital_files), filepath=out_path)

    console.print()
    console.print(f"[green]✓ Converted to {structure.coordinate_type} → {out_path}[/green]")
    console.print("   Original STRU is unchanged.")
    console.print()


# =============================================================================
# Task 404: Fix / Unfix atoms
# =============================================================================


def _choose_directions(console) -> list[int]:
    """Ask which Cartesian directions to apply the fix/unfix to.

    Returns a list of axis indices (0=x, 1=y, 2=z) that should be toggled.
    """
    choice = _prompt_choice(
        console,
        "Directions to constrain",
        ["xyz (all)", "xy", "xz", "yz", "x only", "y only", "z only"],
        "xyz (all)",
    )
    if "all" in choice:
        return [0, 1, 2]
    if choice.startswith("xy"):
        return [0, 1]
    if choice.startswith("xz"):
        return [0, 2]
    if choice.startswith("yz"):
        return [1, 2]
    if "x" in choice:
        return [0]
    if "y" in choice:
        return [1]
    if "z" in choice:
        return [2]
    return [0, 1, 2]


def _get_range_positions(structure) -> tuple[np.ndarray, str]:
    """Return a copy of atomic positions in the structure's display coordinate system.

    Returns (positions, unit_label) where positions use the same coordinate type
    as shown in the atom list (Direct → fractional, Cartesian_* → Å or Bohr).
    """
    if structure.coordinate_type.startswith("Direct"):
        return structure.positions.copy(), "fractional"
    if structure.coordinate_type == "Cartesian_angstrom":
        return structure.positions.copy(), "Å"
    # Cartesian_bohr / Cartesian_au
    return structure.positions.copy(), "Bohr"


@task(404, category="Structure Editing", name="Fix/Unfix Atoms",
      description="Selectively fix or unfix atoms for relaxation constraints")
def task_fix_atoms(args: list[str] | None = None, interactive: bool = True) -> None:
    """Fix or unfix atoms by species, coordinate range, or all.

    Supports per-direction constraints (x, y, z, or xyz).

    In ABACUS STRU, the movement flags (1/0) control whether an atom moves
    during relaxation in each Cartesian direction.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Fix / Unfix Atoms ===[/bold cyan]")
    console.print()

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

    console.print(f"  Atoms: {structure.num_atoms}, Species: {', '.join(structure.species_order)}")
    # Show coordinate type + cell info so the user knows what units the range uses
    ctype = structure.coordinate_type
    cell_a = np.linalg.norm(structure.lattice.cell_angstrom, axis=1)
    if ctype.startswith("Direct"):
        console.print("  Coord: [bold]Direct[/bold] (fractional 0–1)")
        console.print(f"  Cell:  a={cell_a[0]:.3f}  b={cell_a[1]:.3f}  c={cell_a[2]:.3f} Å")
    else:
        console.print(f"  Coord: [bold]{ctype}[/bold]")
    console.print()

    # Show current fix status
    for i, atom in enumerate(structure.atoms):
        flags = "".join("1" if f else "0" for f in atom.fix)
        console.print(f"  [{i+1:3d}] {atom.species:>3s}  {flags}  "
                      f"{atom.position[0]:.3f} {atom.position[1]:.3f} {atom.position[2]:.3f}")

    console.print()
    action = _prompt_choice(
        console, "Action",
        [
            "Fix all atoms",
            "Unfix all atoms",
            "Fix by species",
            "Fix by coordinate range",
        ],
        "Fix by species",
    )

    is_fix = action.startswith("Fix")      # True → set to 0 (fixed)
    is_range = "range" in action
    is_all = "all" in action

    # --- Determine which atoms to affect ---
    indices: list[int] = []

    if is_all:
        indices = list(range(structure.num_atoms))

    elif is_range:
        positions, unit = _get_range_positions(structure)
        axis_label = _prompt_choice(console, "Axis for range",
                                    ["z (c-axis)", "x (a-axis)", "y (b-axis)"],
                                    "z (c-axis)")
        ax = {"x (a-axis)": 0, "y (b-axis)": 1, "z (c-axis)": 2}[axis_label]
        axis_name = "xyz"[ax]

        # Build a helpful hint based on actual cell size
        cell_len = np.linalg.norm(structure.lattice.cell_angstrom[ax])
        if unit == "fractional":
            hint = f"(cell {axis_name}={cell_len:.3f} Å, range 0–1)"
        else:
            hint = f"(cell {axis_name}={cell_len:.3f} Å)"

        rng = _prompt(console,
                      f"  {axis_name} range in {unit} {hint} (min,max)", "")
        try:
            parts = [float(x.strip()) for x in rng.split(",")]
            vmin, vmax = min(parts), max(parts)
        except (ValueError, AttributeError):
            console.print("[red]Invalid range. Use format: min,max[/red]")
            return

        for i, atom in enumerate(structure.atoms):
            if vmin <= positions[i, ax] <= vmax:
                indices.append(i)

        if not indices:
            console.print(f"[yellow]No atoms found in {axis_name} range [{vmin:.4f}, {vmax:.4f}] ({unit}).[/yellow]")
            return
        console.print(f"  {len(indices)} atoms in {axis_name} ∈ [{vmin:.4f}, {vmax:.4f}] ({unit})")

    else:  # by species
        raw = _prompt(console,
                      "Species (comma or space separated, e.g., Li,P,S or Li P S)")
        targets = [s.strip() for s in raw.replace(",", " ").split() if s.strip()]
        if not targets:
            console.print("[yellow]No species entered.[/yellow]")
            return
        for target in targets:
            found = 0
            for i, atom in enumerate(structure.atoms):
                if atom.species == target:
                    indices.append(i)
                    found += 1
            if found == 0:
                console.print(f"[yellow]No atoms of species '{target}' found — skipped.[/yellow]")
        if not indices:
            console.print("[yellow]No matching atoms found for any of the specified species.[/yellow]")
            return

    # --- Choose Cartesian directions ---
    axes = _choose_directions(console)

    # --- Apply ---
    for i in indices:
        current = list(structure.atoms[i].fix)
        for ax in axes:
            current[ax] = not is_fix  # Fix → False=0, Unfix → True=1
        structure.atoms[i].fix = tuple(current)

    # --- Summary ---
    op = "fixed" if is_fix else "unfixed"
    dirs = "".join("xyz"[a] for a in axes)
    console.print()
    console.print(f"[dim]{len(indices)} atoms: {op} in {dirs} direction(s)[/dim]")

    # Write to Fixed.STRU (never overwrite original STRU)
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    _write_stru_bare(structure, is_lcao=bool(structure.orbital_files),
                     filepath="Fixed.STRU")

    console.print()
    console.print("[green]✓ Movement constraints written to Fixed.STRU[/green]")
    console.print("   Original STRU is unchanged.")
    console.print()


# =============================================================================
# Task 405: Build slab from Miller indices
# =============================================================================


@task(405, category="Structure Editing", name="Slab Builder",
      description="Build a slab/surface from Miller indices")
def task_slab_builder(args: list[str] | None = None, interactive: bool = True) -> None:
    """Build a slab structure from Miller indices using ASE."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Slab Builder (Miller Indices) ===[/bold cyan]")
    console.print()

    # Check ASE
    try:
        from ase.build import surface as ase_surface
        from ase.io import read as ase_read
    except ImportError:
        console.print("[red]This task requires ASE. Install with: pip install ase[/red]")
        return

    # Load structure
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        atoms = structure.to_ase()
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms from {stru_path}[/dim]")
    except Exception as e:
        # Try loading from CIF
        cif_files = list(Path(".").glob("*.cif"))
        if cif_files:
            try:
                atoms = ase_read(str(cif_files[0]))
                console.print(f"  [dim]Loaded {len(atoms)} atoms from {cif_files[0].name}[/dim]")
            except Exception:
                console.print(f"[red]Failed to load structure: {e}[/red]")
                return
        else:
            console.print(f"[red]Failed to load structure: {e}[/red]")
            return

    if interactive:
        h = int(_prompt(console, "Miller h", "1"))
        k = int(_prompt(console, "Miller k", "1"))
        l = int(_prompt(console, "Miller l", "1"))
        layers = int(_prompt(console, "Number of atomic layers", "4"))
        vacuum = float(_prompt(console, "Vacuum thickness (Å)", "15.0"))
    else:
        h, k, l = 1, 1, 1
        layers = 4
        vacuum = 15.0

    try:
        slab_atoms = ase_surface(atoms, (h, k, l), layers=layers, vacuum=vacuum)
    except Exception as e:
        console.print(f"[red]Failed to build slab: {e}[/red]")
        return

    console.print(f"  [dim]Slab: {len(slab_atoms)} atoms, "
                  f"cell = {slab_atoms.cell.cellpar()[:3].round(3)} Å[/dim]")

    # Convert to Structure and write
    new_structure = Structure.from_ase(slab_atoms)

    # Set pseudopotential / orbital info from original
    new_structure.pseudo_files = dict(structure.pseudo_files) if 'structure' in dir() else {}
    new_structure.orbital_files = dict(structure.orbital_files) if 'structure' in dir() else {}
    new_structure.species_order = list(structure.species_order) if 'structure' in dir() else sorted(set(slab_atoms.get_chemical_symbols()))
    new_structure.coordinate_type = "Cartesian_angstrom"

    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    is_lcao = bool(new_structure.orbital_files)
    out_path = f"Slab_{h}{k}{l}.STRU"
    _write_stru_bare(new_structure, is_lcao=is_lcao, filepath=out_path)

    console.print()
    console.print(f"[green]✓ Slab written to {out_path} ({h}{k}{l} surface, {layers} layers, {vacuum:.1f} Å vacuum).[/green]")
    console.print("   Original STRU is unchanged.")
    console.print(f"  Atoms: {len(slab_atoms)}, formula: {slab_atoms.get_chemical_formula()}")
    console.print()


# =============================================================================
# Task 406: Add vacuum layer
# =============================================================================


@task(406, category="Structure Editing", name="Vacuum Layer",
      description="Add vacuum layer along a specified cell direction")
def task_vacuum_layer(args: list[str] | None = None, interactive: bool = True) -> None:
    """Add vacuum space along a specified direction (a, b, or c axis)."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Add Vacuum Layer ===[/bold cyan]")
    console.print()

    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    cell = structure.lattice.cell_angstrom
    console.print(f"  Cell: {cell[0,0]:.4f} × {cell[1,1]:.4f} × {cell[2,2]:.4f} Å")

    if interactive:
        direction = _prompt_choice(console, "Vacuum direction",
                                    ["c (z-axis)", "a (x-axis)", "b (y-axis)"],
                                    "c (z-axis)")
        axis_map = {"a (x-axis)": 0, "b (y-axis)": 1, "c (z-axis)": 2}
        ax = axis_map.get(direction, 2)

        vacuum = float(_prompt(console, "Vacuum thickness to add (Å)", "15.0"))
    else:
        ax = 2  # c-axis
        vacuum = 15.0

    # Extend lattice vector along chosen axis, keeping atoms fixed in Cartesian
    old_len = np.linalg.norm(cell[ax])
    new_len = old_len + vacuum
    scale = new_len / old_len

    structure.lattice.vectors[ax] *= scale

    # For Direct (fractional) coords: rescale to keep Cartesian positions fixed.
    # For Cartesian coords: positions are already absolute — no change needed.
    if structure.coordinate_type.startswith("Direct"):
        for atom in structure.atoms:
            atom.position[ax] /= scale

    # Write to new file (preserve original STRU)
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    axis_label = "abc"[ax]
    out_path = f"Vacuum_{axis_label}_{vacuum:.1f}.STRU"
    _write_stru_bare(structure, is_lcao=bool(structure.orbital_files), filepath=out_path)

    new_cell = structure.lattice.cell_angstrom
    console.print()
    console.print(f"[green]✓ Vacuum layer added → {out_path}[/green]")
    console.print("   Original STRU is unchanged.")
    console.print(f"  Cell: {new_cell[0,0]:.4f} × {new_cell[1,1]:.4f} × {new_cell[2,2]:.4f} Å")
    console.print()


# =============================================================================
# Task 407: Shift all atoms
# =============================================================================

@task(407, category="Structure Editing", name="Shift Atoms",
      description="Translate all atoms by a given distance along a chosen direction")
def task_shift_atoms(args: list[str] | None = None, interactive: bool = True) -> None:
    """Shift all atomic positions by a specified displacement vector (in Å)."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Shift Atoms ===[/bold cyan]")
    console.print()

    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    is_direct = structure.coordinate_type.startswith("Direct")
    ctype = structure.coordinate_type
    console.print(f"  Coord type: {ctype}")
    console.print()

    # Direction
    if interactive:
        direction = _prompt_choice(console, "Shift direction",
                                    ["+a (x-axis)", "-a (x-axis)",
                                     "+b (y-axis)", "-b (y-axis)",
                                     "+c (z-axis)", "-c (z-axis)",
                                     "Custom vector"],
                                    "+c (z-axis)")
    else:
        direction = "+c (z-axis)"

    # Map direction to unit vector in Cartesian (Å)
    dir_map = {
        "+a (x-axis)": np.array([1, 0, 0]),
        "-a (x-axis)": np.array([-1, 0, 0]),
        "+b (y-axis)": np.array([0, 1, 0]),
        "-b (y-axis)": np.array([0, -1, 0]),
        "+c (z-axis)": np.array([0, 0, 1]),
        "-c (z-axis)": np.array([0, 0, -1]),
    }

    if "Custom" in direction:
        raw = _prompt(console, "Vector components (x,y,z in Å)")
        try:
            parts = [float(x.strip()) for x in raw.split(",")]
            if len(parts) != 3:
                console.print("[red]Need exactly 3 components (x,y,z).[/red]")
                return
            cart_shift = np.array(parts)
        except ValueError:
            console.print("[red]Invalid vector. Use format: x,y,z (e.g., 1.5,0,-2.0)[/red]")
            return
    else:
        cart_shift = dir_map[direction]

    if interactive:
        distance = float(_prompt(console, "Displacement (Å)", "1.0"))
    else:
        distance = 1.0

    cart_shift = cart_shift * distance

    # Apply shift
    if is_direct:
        # Convert Cartesian shift to fractional shift
        cell = structure.lattice.cell  # Bohr
        # Shift in Å → Bohr → fractional
        from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
        cart_shift_bohr = cart_shift * ANGSTROM_TO_BOHR
        frac_shift = np.linalg.solve(cell.T, cart_shift_bohr)
        for atom in structure.atoms:
            atom.position = atom.position + frac_shift
    else:
        for atom in structure.atoms:
            atom.position = atom.position + cart_shift

    # Write
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    _write_stru_bare(structure, is_lcao=bool(structure.orbital_files),
                     filepath="Shifted.STRU")

    console.print()
    console.print(f"[green]✓ All atoms shifted by {distance:.2f} Å → Shifted.STRU[/green]")
    console.print("   Original STRU is unchanged.")
    console.print(f"  Shift vector (Cartesian Å): ({cart_shift[0]:.4f}, {cart_shift[1]:.4f}, {cart_shift[2]:.4f})")
    console.print()


# =============================================================================
# Task 408: Sort atoms by coordinate
# =============================================================================

@task(408, category="Structure Editing", name="Sort Atoms",
      description="Sort atoms within each species by x, y, or z coordinate")
def task_sort_atoms(args: list[str] | None = None, interactive: bool = True) -> None:
    """Sort atoms within each species group by position along a chosen axis."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Sort Atoms ===[/bold cyan]")
    console.print()

    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists():
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
        console.print(f"  [dim]Loaded {structure.num_atoms} atoms from {stru_path}[/dim]")
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    console.print(f"  Species order: {', '.join(structure.species_order)}")
    console.print()

    if interactive:
        ax = _prompt_choice(console, "Sort by axis",
                            ["z (c-axis)", "x (a-axis)", "y (b-axis)"],
                            "z (c-axis)")
        axis_idx = {"x (a-axis)": 0, "y (b-axis)": 1, "z (c-axis)": 2}[ax]

        order = _prompt_choice(console, "Sort order",
                                ["Ascending (small → large)", "Descending (large → small)"],
                                "Ascending (small → large)")
        ascending = "Ascending" in order
    else:
        axis_idx = 2
        ascending = True

    # Sort within each species group (preserves species grouping for STRU format)
    new_atoms = []
    for sp in structure.species_order:
        group = [(i, a) for i, a in enumerate(structure.atoms) if a.species == sp]
        group.sort(key=lambda x: x[1].position[axis_idx], reverse=not ascending)
        new_atoms.extend([a for _, a in group])

    structure.atoms = new_atoms

    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
    _write_stru_bare(structure, is_lcao=bool(structure.orbital_files),
                     filepath="Sorted.STRU")

    axis_name = "xyz"[axis_idx]
    order_name = "ascending" if ascending else "descending"
    console.print()
    console.print(f"[green]✓ Atoms sorted by {axis_name} ({order_name}) → Sorted.STRU[/green]")
    console.print("   Original STRU is unchanged.")
    console.print()
