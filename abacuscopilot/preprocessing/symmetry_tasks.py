"""Symmetry analysis tasks for crystal structures.

Task IDs 501-599

Uses spglib for space group determination, symmetry operations,
primitive cell reduction, and conventional cell transformation.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.core.models import Atom, Lattice, Structure
from abacuscopilot.tasks import task

# spglib v2+ uses attribute access but still triggers dict-access warnings internally
warnings.filterwarnings("ignore", message=".*dict interface.*", category=DeprecationWarning)


def _get_spglib_attr(obj, key: str, default=None):
    """Get a value from a spglib result, supporting both dict and attribute access."""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


# --- Shared helpers -----------------------------------------------------------

# Element symbol → atomic number mapping (used by both symmetry tasks)
_SPECIES_TO_Z: dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
    "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77,
    "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
    "Po": 84, "At": 85, "Rn": 86,
}

# Reverse mapping for primitive cell task
_Z_TO_SPECIES: dict[int, str] = {v: k for k, v in _SPECIES_TO_Z.items()}


def _get_atomic_numbers(structure) -> list[int]:
    """Extract atomic numbers from a Structure's atoms."""
    numbers = []
    for atom in structure.atoms:
        z = _SPECIES_TO_Z.get(atom.species, 0)
        if z == 0:
            z = 1  # Unknown element → default to H
        numbers.append(z)
    return numbers


def _to_fractional(structure) -> np.ndarray:
    """Return atomic positions in fractional coordinates.

    spglib requires fractional (reduced) coordinates regardless of the
    STRU coordinate type.  This helper converts as needed.
    """
    if structure.coordinate_type.startswith("Direct"):
        # Already fractional — spglib wants these as-is
        return structure.positions.copy()

    # Cartesian → fractional
    cart = structure.positions.copy()
    cell = structure.lattice.cell_angstrom  # same unit as Cartesian_angstrom

    if structure.coordinate_type == "Cartesian_angstrom":
        pass  # cart already in Angstrom
    else:
        # Cartesian_bohr or Cartesian_au: convert to Angstrom first
        from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
        cart = cart * BOHR_TO_ANGSTROM

    # Solve: cart = frac @ cell  ⇒  frac = cart @ cell⁻¹
    return np.linalg.solve(cell.T, cart.T).T


def _get_spglib_info(structure) -> dict[str, Any] | None:
    """Get symmetry information from a Structure using spglib.

    Returns None if spglib is not available or analysis fails.
    """
    try:
        import spglib
    except ImportError:
        return None

    # spglib expects: (cell_angstrom, fractional_positions, atomic_numbers)
    cell = structure.lattice.cell_angstrom
    positions = _to_fractional(structure)       # ← was incorrectly cartesian
    numbers = _get_atomic_numbers(structure)

    cell_tuple = (cell, positions, numbers)

    try:
        dataset = spglib.get_symmetry_dataset(cell_tuple, symprec=1e-3)
        if dataset is None:
            return None

        return {
            "number": _get_spglib_attr(dataset, "number", 0),
            "international": _get_spglib_attr(dataset, "international", "unknown"),
            "hall": _get_spglib_attr(dataset, "hall", "unknown"),
            "wyckoffs": _get_spglib_attr(dataset, "wyckoffs", []),
            "equivalent_atoms": _get_spglib_attr(dataset, "equivalent_atoms", []),
            "primitive_lattice": _get_spglib_attr(dataset, "primitive_lattice", None),
            "rotations": _get_spglib_attr(dataset, "rotations", None),
            "translations": _get_spglib_attr(dataset, "translations", None),
        }

    except Exception:
        return None


# Space group names (1–230) for display
_SPACEGROUP_NAMES: dict[int, str] = {
    1: "P1", 2: "P-1", 3: "P2", 4: "P2_1", 5: "C2",
    14: "P2_1/c", 15: "C2/c",
    19: "P2_12_12_1",
    33: "Pna2_1", 36: "Cmc2_1",
    47: "Pmmm", 51: "Pmma", 57: "Pbcm", 62: "Pnma", 63: "Cmcm",
    74: "Imma",
    99: "P4mm", 123: "P4/mmm", 129: "P4/nmm", 136: "P4_2/mnm",
    139: "I4/mmm",
    143: "P3", 147: "P-3",
    149: "P312", 150: "P321",
    160: "R3m",
    166: "R-3m",
    176: "P6_3/m", 186: "P6_3mc",
    191: "P6/mmm", 193: "P6_3/mcm", 194: "P6_3/mmc",
    221: "Pm-3m", 225: "Fm-3m", 227: "Fd-3m", 229: "Im-3m",
}


def _get_spacegroup_name(number: int) -> str:
    """Get space group name from number, with common ones hard-coded."""
    return _SPACEGROUP_NAMES.get(number, f"#{number}")


# =============================================================================
# Task 501: Symmetry analysis
# =============================================================================

@task(501, category="Symmetry", name="Symmetry Analysis",
      description="Analyze crystal symmetry (space group, Wyckoff positions) using spglib")
def task_symmetry_analysis(args: list[str] | None = None, interactive: bool = True) -> None:
    """Analyze the symmetry of a crystal structure from STRU file."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Symmetry Analysis ===[/bold cyan]")
    console.print()

    # Load STRU
    stru_path = "STRU"
    if args:
        for arg in args:
            if Path(arg).exists() and arg != "STRU":
                stru_path = arg
                break

    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
    except Exception as e:
        console.print(f"[red]Failed to read STRU file '{stru_path}': {e}[/red]")
        return

    console.print(f"  [dim]Structure: {structure.num_atoms} atoms, "
                  f"{structure.num_species} species[/dim]")

    # Symmetry analysis via spglib
    info = _get_spglib_info(structure)

    if info is None:
        console.print("[yellow]spglib is not installed. Install with: pip install spglib[/yellow]")
        console.print("[dim]Only basic structure information is available.[/dim]")
        console.print()
        console.print(f"  Lattice constant: {structure.lattice.constant:.6f} Bohr")
        console.print(f"  Coordinate type: {structure.coordinate_type}")
        console.print(f"  Species: {structure.species_order}")
        return

    console.print()
    console.print(f"  [bold green]Space Group:[/bold green] {info['international']} "
                  f"({_get_spacegroup_name(info['number'])})")
    console.print(f"  [bold]Number:[/bold] {info['number']}")
    console.print(f"  [bold]Hall symbol:[/bold] {info['hall']}")
    console.print(f"  [bold]Number of symmetry operations:[/bold] {len(info.get('rotations', []))}")
    console.print()

    # Wyckoff positions
    if info.get("wyckoffs"):
        console.print("  [bold]Wyckoff positions:[/bold]")
        equivalent = info.get("equivalent_atoms", [])
        for i, atom in enumerate(structure.atoms):
            wyck = info["wyckoffs"][i] if i < len(info["wyckoffs"]) else "?"
            eq = equivalent[i] if i < len(equivalent) else i
            console.print(f"    {atom.species} atom {i+1}: {wyck} (equiv. group {eq})")
        console.print()

    # Primitive cell
    if info.get("primitive_lattice") is not None:
        prim_cell = np.array(info["primitive_lattice"])
        console.print("  [bold]Primitive cell vectors (Angstrom):[/bold]")
        for row in prim_cell:
            console.print(f"    {row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}")
        console.print()


# =============================================================================
# Task 502: Primitive cell reduction
# =============================================================================

@task(502, category="Symmetry", name="Primitive Cell",
      description="Reduce structure to its primitive cell using spglib")
def task_primitive_cell(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a primitive-cell STRU file."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Primitive Cell Reduction ===[/bold cyan]")
    console.print()

    stru_path = "STRU"
    try:
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru(stru_path)
    except Exception as e:
        console.print(f"[red]Failed to read STRU: {e}[/red]")
        return

    info = _get_spglib_info(structure)

    if info is None:
        console.print("[red]spglib is required for primitive cell reduction.[/red]")
        return

    import spglib

    # spglib.find_primitive expects: (cell_angstrom, fractional_positions, atomic_numbers)
    cell = structure.lattice.cell_angstrom
    positions = _to_fractional(structure)       # ← was incorrectly cartesian
    numbers = _get_atomic_numbers(structure)

    try:
        result = spglib.find_primitive(
            (cell, positions, numbers), symprec=1e-3
        )
        if result is None:
            console.print("[yellow]Could not find primitive cell.[/yellow]")
            return

        prim_cell, prim_scaled_pos, prim_numbers = result
    except Exception as e:
        console.print(f"[red]spglib error: {e}[/red]")
        return

    # Build new Structure
    prim_structure = Structure()
    prim_structure.lattice = Lattice.from_cell_angstrom(np.array(prim_cell))
    prim_structure.coordinate_type = "Direct"

    species_count = {}

    # Map species labels consistently
    atom_labels = []
    for z in prim_numbers:
        label = _Z_TO_SPECIES.get(z, f"Z{z}")
        atom_labels.append(label)
        if label not in species_count:
            species_count[label] = 0
            prim_structure.species_order.append(label)
        species_count[label] += 1

    for i, (label, pos) in enumerate(zip(atom_labels, prim_scaled_pos)):
        prim_structure.atoms.append(Atom(species=label, position=pos))

    # Copy pseudopotential info
    for species in prim_structure.species_order:
        if species in structure.pseudo_files:
            prim_structure.pseudo_files[species] = structure.pseudo_files[species]

    # Write to Primitive.STRU (never overwrite original STRU)
    from abacuscopilot.io.stru_file import write_stru
    write_stru(prim_structure, filepath="Primitive.STRU",
               is_lcao=bool(structure.orbital_files))

    console.print()
    console.print("[green]✓ Primitive cell written to Primitive.STRU[/green]")
    console.print("   Original STRU is unchanged.")
    console.print(f"  Original: {structure.num_atoms} atoms")
    console.print(f"  Primitive: {len(prim_structure.atoms)} atoms")
    console.print()
