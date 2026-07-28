"""K-point generation tasks for ABACUS calculations.

Task IDs 301-399

Automatic Monkhorst-Pack mesh generation from k-spacing,
line-mode k-path for band structure from high-symmetry points,
and custom k-point list generation.
"""

from __future__ import annotations

from pathlib import Path

from abacuscopilot.config import load_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.constants import DEFAULT_KSPACING
from abacuscopilot.core.models import KPoints, Lattice
from abacuscopilot.tasks import task

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

    Returns (KPoints, path_labels_str, bz_data) or raises on failure.
    bz_data carries the reciprocal lattice / point coords / path for an
    optional Brillouin-zone plot.
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

    # Build line-mode segments directly, keeping each segment's true
    # start/end labels. seekpath breaks the path at equivalent points
    # (e.g. FCC: ...->X, X->U, K->GAMMA...), where segment i's end label
    # differs from segment i+1's start label. We must NOT flatten labels to
    # "one label per point" or the coords and labels desync (U gets mislabelled
    # as K). At a break, join the two labels as "END|START" for display.
    from abacuscopilot.core.models import KPoints

    def _disp(lbl: str) -> str:
        return "Γ" if lbl.upper() == "GAMMA" else lbl

    kpts = KPoints(mode="line")
    disp_path_labels = []
    for i, (start_label, end_label) in enumerate(path):
        start = list(point_coords[start_label])
        end = list(point_coords[end_label])
        kpts.line_path.append({
            "start": tuple(start),
            "end": tuple(end),
            "npoints": npts,
            "label": start_label,
            "end_label": end_label,
        })
        # Build a human-readable path string that shows breaks
        if i == 0:
            disp_path_labels.append(_disp(start_label))
        elif start_label != path[i - 1][1]:
            # discontinuity: previous end and this start differ
            disp_path_labels[-1] = f"{disp_path_labels[-1]}|{_disp(start_label)}"
        disp_path_labels.append(_disp(end_label))

    path_str = " — ".join(disp_path_labels)
    # Data needed to draw the path inside the Brillouin zone (optional plot).
    # Segments come straight from kpts.line_path so the plot matches the KPT
    # file exactly (incl. FCC-style U/K breaks).
    bz_data = {
        "recip_lattice": result["reciprocal_primitive_lattice"],
        "segments": kpts.line_path,
    }
    return kpts, path_str, bz_data


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

    # --- Generate k-path via seekpath (required) ---
    try:
        kpts, path_str, bz_data = _generate_kpt_via_seekpath(structure, npts)
        console.print(f"\n  [bold]High-symmetry path (seekpath):[/bold] {path_str}")
    except ImportError:
        console.print("[red]seekpath is required for automatic k-path detection.[/red]")
        console.print("[dim]Install it with: pip install seekpath[/dim]")
        return
    except Exception as e:
        console.print(f"[red]seekpath failed to generate a k-path: {e}[/red]")
        console.print("[dim]Check that STRU is a valid crystal structure.[/dim]")
        return

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

    # --- Optional: plot the path inside the Brillouin zone ---
    if interactive:
        want_bz = "Yes" in _prompt_choice(
            console, "Plot the path in the Brillouin zone?", ["Yes", "No"], "Yes"
        )
    else:
        want_bz = True
    if want_bz:
        from abacuscopilot.plotting.brillouin import plot_brillouin_zone
        out = plot_brillouin_zone(
            bz_data["recip_lattice"], bz_data["segments"],
            out_png="brillouin_zone.png", title="Brillouin zone — band path",
        )
        if out:
            console.print(f"  [green]✓ Brillouin-zone plot: {out}[/green]")
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

    # Generate q-path via seekpath (required)
    try:
        kpts, path_str, bz_data = _generate_kpt_via_seekpath(structure, npts)
        console.print(f"\n  [bold]High-symmetry path (seekpath):[/bold] {path_str}")
    except ImportError:
        console.print("[red]seekpath is required for automatic q-path detection.[/red]")
        console.print("[dim]Install it with: pip install seekpath[/dim]")
        return
    except Exception as e:
        console.print(f"[red]seekpath failed to generate a q-path: {e}[/red]")
        console.print("[dim]Check that STRU is a valid crystal structure.[/dim]")
        return

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

    # --- Optional: plot the q-path inside the Brillouin zone ---
    if interactive:
        want_bz = "Yes" in _prompt_choice(
            console, "Plot the path in the Brillouin zone?", ["Yes", "No"], "Yes"
        )
    else:
        want_bz = True
    if want_bz:
        from abacuscopilot.plotting.brillouin import plot_brillouin_zone
        out = plot_brillouin_zone(
            bz_data["recip_lattice"], bz_data["segments"],
            out_png="brillouin_zone_phonon.png",
            title="Brillouin zone — phonon q-path",
        )
        if out:
            console.print(f"  [green]✓ Brillouin-zone plot: {out}[/green]")
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
