"""K-point generation tasks for ABACUS calculations.

Task IDs 301-399

Automatic Monkhorst-Pack mesh generation from k-spacing,
line-mode k-path for band structure from high-symmetry points,
and custom k-point list generation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.config import load_config
from abacuscopilot.console_utils import _get_console, _prompt
from abacuscopilot.core.constants import DEFAULT_KSPACING
from abacuscopilot.core.models import KPoints, Lattice
from abacuscopilot.tasks import task

# Common high-symmetry k-path definitions for standard Bravais lattices
# Format: list of (label, coordinates) where coordinates are in fractional units

HIGH_SYMMETRY_PATHS: dict[str, list[tuple[str, tuple[float, float, float]]]] = {
    "SC": [
        ("GAMMA", (0, 0, 0)), ("X", (0.5, 0, 0)), ("M", (0.5, 0.5, 0)),
        ("R", (0.5, 0.5, 0.5)), ("GAMMA", (0, 0, 0)),
    ],
    "FCC": [
        ("GAMMA", (0, 0, 0)), ("X", (0.5, 0, 0.5)), ("U", (0.625, 0.25, 0.625)),
        ("K", (0.375, 0.375, 0.75)), ("GAMMA", (0, 0, 0)), ("L", (0.5, 0.5, 0.5)),
        ("W", (0.5, 0.25, 0.75)), ("X", (0.5, 0, 0.5)),
    ],
    "BCC": [
        ("GAMMA", (0, 0, 0)), ("H", (0.5, -0.5, 0.5)), ("N", (0, 0, 0.5)),
        ("P", (0.25, 0.25, 0.25)), ("GAMMA", (0, 0, 0)), ("N", (0, 0, 0.5)),
    ],
    "HEX": [
        ("GAMMA", (0, 0, 0)), ("M", (0.5, 0, 0)), ("K", (1/3, 1/3, 0)),
        ("GAMMA", (0, 0, 0)), ("A", (0, 0, 0.5)), ("L", (0.5, 0, 0.5)),
        ("H", (1/3, 1/3, 0.5)), ("A", (0, 0, 0.5)),
    ],
    "TET": [
        ("GAMMA", (0, 0, 0)), ("X", (0.5, 0, 0)), ("M", (0.5, 0.5, 0)),
        ("GAMMA", (0, 0, 0)), ("Z", (0, 0, 0.5)), ("R", (0.5, 0, 0.5)),
        ("A", (0.5, 0.5, 0.5)), ("Z", (0, 0, 0.5)),
    ],
    "ORC": [
        ("GAMMA", (0, 0, 0)), ("X", (0.5, 0, 0)), ("S", (0.5, 0.5, 0)),
        ("Y", (0, 0.5, 0)), ("GAMMA", (0, 0, 0)), ("Z", (0, 0, 0.5)),
        ("U", (0.5, 0, 0.5)), ("R", (0.5, 0.5, 0.5)), ("T", (0, 0.5, 0.5)),
        ("Z", (0, 0, 0.5)),
    ],
}

# Aliases
HIGH_SYMMETRY_PATHS["CUB"] = HIGH_SYMMETRY_PATHS["SC"]
HIGH_SYMMETRY_PATHS["HCP"] = HIGH_SYMMETRY_PATHS["HEX"]


def _guess_lattice_type(lattice: Lattice, structure=None) -> str:
    """Guess the Bravais lattice type from cell vectors and (optionally) spglib.

    Uses spglib for reliable SC/BCC/FCC distinction when a structure is
    available; falls back to cell-vector heuristics otherwise.
    """
    from abacuscopilot.core.constants import BOHR_TO_ANGSTROM

    cell = lattice.cell * BOHR_TO_ANGSTROM  # Angstrom
    a_vec, b_vec, c_vec = cell[0], cell[1], cell[2]

    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)

    alpha = np.degrees(np.arccos(np.dot(b_vec, c_vec) / (b * c)))
    beta = np.degrees(np.arccos(np.dot(a_vec, c_vec) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(a_vec, b_vec) / (a * b)))

    # Cubic: a == b == c, all angles 90
    if (abs(a - b) / a < 0.01 and abs(b - c) / a < 0.01 and
        abs(alpha - 90) < 1 and abs(beta - 90) < 1 and abs(gamma - 90) < 1):
        # Try spglib to distinguish SC / BCC / FCC
        if structure is not None:
            try:
                from abacuscopilot.preprocessing.symmetry_tasks import _get_spglib_info
                info = _get_spglib_info(structure)
                if info and info.get("number"):
                    sg = info["number"]
                    # FCC: 195-230 except BCC (229, 230...)
                    if sg in (216, 225, 227):
                        return "FCC"
                    # BCC: 229, 230
                    if sg in (229, 230, 217, 220):
                        return "BCC"
            except Exception:
                pass
        return "CUB"

    # Hexagonal: a == b, gamma == 120
    if abs(a - b) / a < 0.01 and abs(gamma - 120) < 1:
        return "HEX"

    # Tetragonal: a == b != c, all angles 90
    if (abs(a - b) / a < 0.01 and abs(a - c) / a > 0.01 and
        abs(alpha - 90) < 1 and abs(beta - 90) < 1 and abs(gamma - 90) < 1):
        return "TET"

    # Orthorhombic: all angles 90
    if abs(alpha - 90) < 1 and abs(beta - 90) < 1 and abs(gamma - 90) < 1:
        return "ORC"

    return "CUB"  # Default


# =============================================================================
# Task 301: Auto MP K-point mesh
# =============================================================================

@task(301, category="KPT", name="Auto KPT (MP mesh)",
      description="Generate automatic Monkhorst-Pack KPT file from k-spacing")
def task_auto_kpt(args: list[str] | None = None, interactive: bool = True) -> None:
    """Auto-generate a KPT file with Monkhorst-Pack mesh from k-spacing."""
    console = _get_console()
    config = load_config()
    defaults = config.get("defaults", {})

    console.print()
    console.print("[bold cyan]=== Generate Auto KPT (MP mesh) ===[/bold cyan]")
    console.print()

    kspacing = DEFAULT_KSPACING
    gamma_centered = True

    # Always try to load STRU for lattice information (both interactive and CLI modes)
    lattice = None
    for stru_path in ("STRU", "stru"):
        if Path(stru_path).exists():
            try:
                from abacuscopilot.io.stru_file import read_stru
                structure = read_stru(stru_path)
                lattice = structure.lattice
                console.print(f"  [dim]Loaded structure from {stru_path}[/dim]")
                break
            except Exception:
                pass

    if interactive:
        kspacing = float(_prompt(console, "K-spacing (2π/Å)",
                                  defaults.get("kspacing", DEFAULT_KSPACING)))
        gamma_centered_input = _prompt(console, "Gamma-centered? (y/n)", "y")
        gamma_centered = gamma_centered_input.lower() in ("y", "yes", "true", "1")

    # Generate
    from abacuscopilot.io.kpt_file import auto_mp_kpts

    if lattice is None:
        lattice = Lattice()

    kpts = auto_mp_kpts(lattice, kspacing=kspacing, gamma_centered=gamma_centered)

    # Write
    from abacuscopilot.io.kpt_file import write_kpt
    write_kpt(kpts)

    console.print()
    console.print("[green]✓ KPT file written successfully.[/green]")
    console.print(f"  Mode: {'Gamma' if gamma_centered else 'MP'}-centered")
    console.print(f"  K-spacing: {kspacing} 2π/Å")
    console.print(f"  Grid: {kpts.grid[0]}×{kpts.grid[1]}×{kpts.grid[2]}")
    console.print(f"  Total k-points: {kpts.grid[0] * kpts.grid[1] * kpts.grid[2]}")
    console.print()


# =============================================================================
# Task 302: Band path KPT — uses seekpath for automatic path generation
# =============================================================================


def _get_stru_for_seekpath(structure):
    """Convert a Structure into the (cell, positions, numbers) tuple seekpath expects.

    Returns None if the conversion fails.
    """
    try:
        from abacuscopilot.preprocessing.symmetry_tasks import (
            _get_atomic_numbers,
            _to_fractional,
        )
    except ImportError:
        return None

    cell = structure.lattice.cell_angstrom          # 3×3, Angstrom
    positions = _to_fractional(structure)           # N×3, fractional
    numbers = _get_atomic_numbers(structure)        # N atomic numbers
    return (cell, positions, numbers)


def _generate_kpt_via_seekpath(structure, npts: int):
    """Use seekpath to auto-generate a high-symmetry k-path.

    Returns (KPoints, path_labels_str) or raises an exception on failure.
    """
    import seekpath

    cell_tuple = _get_stru_for_seekpath(structure)
    if cell_tuple is None:
        raise RuntimeError("Cannot convert structure for seekpath")

    result = seekpath.get_path(
        cell_tuple,
        with_time_reversal=True,
        recipe="hpkot",
        threshold=1e-07,
        symprec=1e-05,
        angle_tolerance=-1.0,
    )

    point_coords = result["point_coords"]   # {"GAMMA": [0,0,0], "X": [0.5,0,0.5], ...}
    path = result["path"]                   # [("GAMMA","X"), ("X","U"), ...]

    # Build line-mode segments
    segments = []
    labels = []
    for i, (start_label, end_label) in enumerate(path):
        start = list(point_coords[start_label])
        end = list(point_coords[end_label])
        segments.append((start, end, npts))
        labels.append(start_label)
    labels.append(path[-1][1])  # final endpoint label

    from abacuscopilot.io.kpt_file import line_mode_kpts_from_path
    kpts = line_mode_kpts_from_path(segments, labels=labels)

    path_str = " — ".join(labels)
    return kpts, path_str


@task(302, category="KPT", name="KPT (band path)",
      description="Generate KPT file with high-symmetry path for band structure")
def task_band_kpt(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a line-mode KPT file for band structure using seekpath.

    Automatically detects the space group and high-symmetry k-path from STRU.
    Falls back to built-in paths if seekpath is not installed.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate Band Path KPT ===[/bold cyan]")
    console.print()

    # --- Load STRU ---
    structure = None
    for stru_path in ("STRU", "stru"):
        if Path(stru_path).exists():
            try:
                from abacuscopilot.io.stru_file import read_stru
                structure = read_stru(stru_path)
                console.print(f"  [dim]Loaded structure from {stru_path}[/dim]")
                break
            except Exception:
                pass

    if structure is None:
        console.print("[red]No STRU file found in current directory.[/red]")
        return

    # --- Determine npts ---
    if interactive:
        npts = int(_prompt(console, "Points per segment", "20"))
    else:
        npts = 20

    # --- Generate k-path ---
    use_seekpath = False
    try:
        kpts, path_str = _generate_kpt_via_seekpath(structure, npts)
        use_seekpath = True
        console.print(f"\n  [bold]High-symmetry path (seekpath):[/bold] {path_str}")
    except Exception:
        # Fallback to built-in paths
        lattice_type = _guess_lattice_type(structure.lattice, structure)
        if lattice_type not in HIGH_SYMMETRY_PATHS:
            lattice_type = "CUB"
        path = HIGH_SYMMETRY_PATHS[lattice_type]

        segments = []
        labels = []
        for i in range(len(path) - 1):
            start_label, start_coords = path[i]
            end_label, end_coords = path[i + 1]
            segments.append((list(start_coords), list(end_coords), npts))
            labels.append(start_label)
        labels.append(path[-1][0])

        from abacuscopilot.io.kpt_file import line_mode_kpts_from_path
        kpts = line_mode_kpts_from_path(segments, labels=labels)
        console.print(f"\n  [bold]High-symmetry path ({lattice_type}):[/bold] {' — '.join(labels)}")
        console.print("  [dim](install seekpath for automatic detection: pip install seekpath)[/dim]")

    console.print()

    # Write
    from abacuscopilot.io.kpt_file import write_kpt
    write_kpt(kpts)

    # Build labels string from kpts.line_path for display
    disp_labels = []
    for seg in kpts.line_path:
        lbl = seg.get("label", "")
        if lbl:
            disp_labels.append(lbl)
    if kpts.line_path:
        end_lbl = kpts.line_path[-1].get("end_label", "")
        if end_lbl:
            disp_labels.append(end_lbl)

    console.print()
    console.print("[green]✓ Band path KPT file written successfully.[/green]")
    console.print(f"  Path: {' — '.join(disp_labels)}")
    console.print(f"  Total segments: {len(kpts.line_path)}")
    console.print()


# =============================================================================
# Task 304: Phonopy band.conf generator
# =============================================================================

@task(304, category="KPT", name="KPT (phonon)",
      description="Generate phonopy band.conf for phonon dispersion from STRU")
def task_phonon_kpt(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a phonopy band.conf file for phonon dispersion calculation.

    Reads STRU, uses seekpath for the q-point path, and writes band.conf
    ready for ``phonopy band.conf --abacus``.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate Phonopy band.conf ===[/bold cyan]")
    console.print()

    # --- Load STRU ---
    structure = None
    for stru_path in ("STRU", "stru"):
        if Path(stru_path).exists():
            try:
                from abacuscopilot.io.stru_file import read_stru
                structure = read_stru(stru_path)
                console.print(f"  [dim]Loaded structure from {stru_path}[/dim]")
                break
            except Exception:
                pass

    if structure is None:
        console.print("[red]No STRU file found in current directory.[/red]")
        return

    # Species list for ATOM_NAME
    species = structure.species_order
    atom_name = " ".join(species)

    if interactive:
        dim_in = _prompt(console, "Supercell DIM (e.g. 2 2 2)", "2 2 2")
        dim = dim_in.strip()
        mesh_in = _prompt(console, "MP mesh for force constants (e.g. 8 8 8)", "8 8 8")
        mesh = mesh_in.strip()
        npts = int(_prompt(console, "Points per q-segment", "101"))
    else:
        dim, mesh, npts = "2 2 2", "8 8 8", 101

    # Generate q-path via seekpath
    try:
        kpts, path_str = _generate_kpt_via_seekpath(structure, npts)
        console.print(f"\n  [bold]High-symmetry path (seekpath):[/bold] {path_str}")
    except Exception:
        lattice_type = _guess_lattice_type(structure.lattice, structure)
        if lattice_type not in HIGH_SYMMETRY_PATHS:
            lattice_type = "CUB"
        path = HIGH_SYMMETRY_PATHS[lattice_type]
        segments = []
        labels = []
        for i in range(len(path) - 1):
            start_label, start_coords = path[i]
            end_label, end_coords = path[i + 1]
            segments.append((list(start_coords), list(end_coords), npts))
            labels.append(start_label)
        labels.append(path[-1][0])
        from abacuscopilot.io.kpt_file import line_mode_kpts_from_path
        kpts = line_mode_kpts_from_path(segments, labels=labels)
        console.print(f"\n  [bold]High-symmetry path ({lattice_type}):[/bold] {' — '.join(labels)}")

    # Build BAND line (coordinates only) and BAND_LABELS
    band_coords = []
    band_labels_list = []
    for seg in kpts.line_path:
        s = seg["start"]
        band_coords.append(f"{s[0]:.6f} {s[1]:.6f} {s[2]:.6f}")
        lbl = seg.get("label", "")
        if lbl:
            band_labels_list.append(lbl)
    if kpts.line_path:
        last = kpts.line_path[-1]
        e = last["end"]
        band_coords.append(f"{e[0]:.6f} {e[1]:.6f} {e[2]:.6f}")
        elbl = last.get("end_label", "")
        if elbl:
            band_labels_list.append(elbl)

    band_line = "  ".join(band_coords)
    labels_line = " ".join(
        lab.replace("GAMMA", "Γ") for lab in band_labels_list
    )

    # Assemble band.conf
    content = f"""ATOM_NAME = {atom_name}
DIM = {dim}
MESH = {mesh}
BAND = {band_line}
BAND_LABELS = {labels_line}
BAND_POINTS = {npts}
BAND_CONNECTION = .TRUE.
"""

    out_path = "band.conf"
    with open(out_path, "w") as f:
        f.write(content)

    console.print()
    console.print("[green]✓ band.conf written.[/green]")
    console.print(f"  ATOM_NAME = {atom_name}")
    console.print(f"  DIM = {dim}, MESH = {mesh}, BAND_POINTS = {npts}")
    console.print(f"  Run: phonopy -d --dim=\"{dim}\" --abacus && phonopy band.conf --abacus")
    console.print()


# =============================================================================
# Task 303: Custom KPT
# =============================================================================

@task(303, category="KPT", name="Custom KPT",
      description="Generate KPT file with custom k-point list")
def task_custom_kpt(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a KPT file from a user-specified custom list."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate Custom KPT ===[/bold cyan]")
    console.print()
    console.print("[dim]Enter k-point coordinates and weights (one per line).[/dim]")
    console.print("[dim]Format: kx ky kz weight[/dim]")
    console.print("[dim]Enter blank line when done.[/dim]")
    console.print()

    kpoints_list = []
    while True:
        line = console.input("  k-point: ").strip()
        if not line:
            break
        parts = line.split()
        if len(parts) == 3:
            kx, ky, kz = float(parts[0]), float(parts[1]), float(parts[2])
            kpoints_list.append((kx, ky, kz, 1.0))
        elif len(parts) >= 4:
            kx, ky, kz, w = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
            kpoints_list.append((kx, ky, kz, w))

    if not kpoints_list:
        console.print("[yellow]No k-points entered. Using Gamma-only.[/yellow]")
        kpoints_list = [(0.0, 0.0, 0.0, 1.0)]

    kpts = KPoints(mode="direct", explicit_kpoints=kpoints_list)

    from abacuscopilot.io.kpt_file import write_kpt
    write_kpt(kpts)

    console.print()
    console.print(f"[green]✓ Custom KPT file written with {len(kpoints_list)} k-point(s).[/green]")
    console.print()
