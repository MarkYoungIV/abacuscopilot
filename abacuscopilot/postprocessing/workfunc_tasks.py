"""Work function analysis tasks for ABACUS output.

Task IDs 1101-1102.

Reads electrostatic potential from OUT.ABACUS, computes 1D planar
average along the vacuum direction, and extracts the work function.

Method:
1. Read LOCPOT / electrostatic potential cube file
2. Compute 1D planar average along vacuum direction (typically z)
3. Identify vacuum level from the plateau region
4. Compute work function:  Φ = V_vacuum − E_Fermi
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.core.constants import BOHR_TO_ANGSTROM, RY_TO_EV
from abacuscopilot.tasks import task

# =============================================================================
# Potential file finder
# =============================================================================


def _find_potential_file(root: str | Path = ".") -> Path | None:
    """Find the newest electrostatic-potential file below *root*.

    ABACUS writes ``ElecStaticPot.cube`` below ``OUT.*``.  A few older
    workflows use ``POT.cube`` or VASP's ``LOCPOT`` name, so those aliases are
    kept for compatibility with the original task.
    """
    root = Path(root)
    candidates = (
        list(root.glob("ElecStaticPot*.cube")) +
        list(root.glob("LOCPOT*")) +
        list(root.glob("POT.cube")) +
        list(root.glob("OUT.*/ElecStaticPot*.cube")) +
        list(root.glob("OUT.*/LOCPOT*")) +
        list(root.glob("OUT.*/POT.cube"))
    )
    candidates = [path for path in candidates if path.is_file()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def read_vacuum_level(filepath: str | Path) -> float | None:
    """Read a vacuum level in eV from a precomputed text override.

    The parser accepts ``E_VACUUM (eV) = ...`` and two-sided reports.  This is
    only an explicit compatibility override; normal tasks derive the level
    directly from the raw cube.
    """
    path = Path(filepath)
    if not path.is_file():
        return None

    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    direct = re.compile(
        rf"(?:E|V)[_ ]?VACUUM(?:\s*\([^)]*\))?\s*=\s*({number})", re.I
    )
    side = re.compile(
        rf"(?:LOWER|UPPER)[_ ]?VACUUM(?:\s*\([^)]*\))?\s*=\s*({number})", re.I
    )

    side_values = []
    for line in path.read_text(errors="ignore").splitlines():
        match = direct.search(line)
        if match:
            return float(match.group(1))
        match = side.search(line)
        if match:
            side_values.append(float(match.group(1)))
    if side_values:
        return float(np.mean(side_values))
    return None


def _find_scf_log(
    root: str | Path = ".",
    potential_path: str | Path | None = None,
    explicit: str | Path | None = None,
) -> Path | None:
    """Find the newest SCF log associated with a potential output."""
    if explicit:
        path = Path(explicit)
        return path if path.is_file() else None

    root = Path(root)
    search_roots = []
    if potential_path is not None:
        potential = Path(potential_path)
        search_roots.extend([potential.parent, potential.parent.parent])
    search_roots.append(root)

    candidates = []
    for directory in search_roots:
        candidates.extend(directory.glob("running_*.log"))
        candidates.extend(directory.glob("OUT.*/running_*.log"))
    candidates = [path for path in candidates if path.is_file()]
    return max(candidates, key=lambda p: p.stat().st_mtime) if candidates else None


def read_fermi_energy(filepath: str | Path) -> float | None:
    """Extract the final ABACUS Fermi energy in eV from an SCF log.

    ABACUS commonly prints ``E_Fermi <Ry> <eV>``; newer builds may spell this
    as ``Fermi energy (eV) = ...``.  The last matching line is used so relax
    and MD logs resolve to the final electronic step.
    """
    path = Path(filepath)
    if not path.is_file():
        return None

    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    ev_value = re.compile(rf"({number})\s*eV\b", re.I)
    fermi_line = re.compile(r"(?:E_Fermi|E-fermi|EFERMI|Fermi\s+(?:energy|level))", re.I)

    for line in reversed(path.read_text(errors="ignore").splitlines()):
        if not fermi_line.search(line):
            continue
        values = ev_value.findall(line)
        if values:
            return float(values[-1])

        # Standard ABACUS text has two bare numbers: Ry first, eV second.
        tail = line[fermi_line.search(line).end():]
        values = re.findall(number, tail)
        if len(values) >= 2:
            return float(values[-1])
        if values:
            return float(values[-1])
    return None


# =============================================================================
# Cube potential reader
# =============================================================================


def _read_potential_cube_details(
    filepath: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Read cube data plus atom coordinates needed for vacuum detection."""
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    lines = filepath.read_text(errors="ignore").splitlines()
    if len(lines) < 6:
        return None

    # Cube files have two comments followed by a geometry block.  Searching
    # for that block also handles ABACUS files with an optional extra comment.
    geometry = None
    for index in range(min(40, len(lines) - 3)):
        parts = lines[index].split()
        if len(parts) != 4:
            continue
        try:
            natom = abs(int(parts[0]))
            origin = np.array([float(x) for x in parts[1:4]])
            dimensions = []
            cell = np.zeros((3, 3))
            for axis in range(3):
                row = lines[index + axis + 1].split()
                if len(row) != 4:
                    raise ValueError
                n = abs(int(row[0]))
                dimensions.append(n)
                cell[axis] = np.array([float(x) for x in row[1:4]]) * n
            geometry = (index, natom, origin, dimensions, cell)
            break
        except (TypeError, ValueError, IndexError):
            continue

    if geometry is None:
        return None

    geometry_index, natom, origin, dimensions, cell = geometry
    data_start = geometry_index + 4 + natom
    if data_start >= len(lines):
        return None

    atom_positions = []
    for line in lines[geometry_index + 4:data_start]:
        parts = line.split()
        if len(parts) >= 5:
            try:
                atom_positions.append([float(parts[2]), float(parts[3]), float(parts[4])])
            except ValueError:
                return None

    raw = []
    for line in lines[data_start:]:
        raw.extend(float(x) for x in line.split())

    nx, ny, nz = dimensions
    expected = nx * ny * nz
    if len(raw) < expected:
        return None

    # Gaussian/ABACUS cube data are ordered with the last grid axis varying
    # fastest.  Keeping (x, y, z) here makes the z planar average unambiguous.
    data = np.asarray(raw[:expected], dtype=float).reshape(nx, ny, nz)
    atoms = np.asarray(atom_positions, dtype=float).reshape((-1, 3))
    return data, cell, origin, atoms


def read_potential_cube(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read electrostatic potential from a Cube-format file.

    Returns ``(data_3d, cell_bohr, origin)`` for compatibility with the
    original helper.  Task 1101/1102 additionally use the internal reader's
    atom coordinates to locate the vacuum gap directly from the raw cube.
    """
    result = _read_potential_cube_details(filepath)
    if result is None:
        return None
    data, cell, origin, _atoms = result
    return data, cell, origin


# =============================================================================
# Work function extraction
# =============================================================================


def extract_work_function(
    potential_1d: np.ndarray,
    z_axis: np.ndarray,
    e_fermi: float = 0.0,
    vacuum_fraction: float = 0.3,
) -> dict[str, Any]:
    """Extract the work function from a 1D planar-averaged potential.

    The vacuum level is estimated as the average of the potential in the
    top `vacuum_fraction` of the cell (where the vacuum layer is).

    Args:
        potential_1d: 1D potential along z, in eV.
        z_axis: z-coordinate grid.
        e_fermi: Fermi energy in eV.
        vacuum_fraction: Fraction of cell to consider as vacuum region.

    Returns:
        Dict with 'vacuum_level', 'work_function', 'e_fermi', 'z_vacuum'.
    """
    potential_1d = np.asarray(potential_1d, dtype=float)
    z_axis = np.asarray(z_axis, dtype=float)
    if potential_1d.ndim != 1 or z_axis.ndim != 1 or len(potential_1d) != len(z_axis):
        raise ValueError("potential_1d and z_axis must be one-dimensional arrays of equal length")
    if len(potential_1d) == 0:
        raise ValueError("potential_1d cannot be empty")
    if not 0.0 < vacuum_fraction <= 1.0:
        raise ValueError("vacuum_fraction must be between 0 and 1")

    npts = len(potential_1d)
    vac_start = int(npts * (1 - vacuum_fraction))
    vac_start = min(max(vac_start, 0), npts - 1)

    # Find vacuum level: average over the vacuum region
    vac_region = potential_1d[vac_start:]
    vacuum_level = float(vac_region.mean())

    # Work function
    work_function = vacuum_level - e_fermi

    return {
        "vacuum_level": vacuum_level,
        "work_function": work_function,
        "e_fermi": e_fermi,
        "z_vacuum_start": float(z_axis[vac_start]),
    }


def _largest_periodic_gap(atom_coordinates: np.ndarray, cell_length: float) -> tuple[float, float, float] | None:
    """Return the largest atom-free interval on a periodic one-dimensional cell."""
    if cell_length <= 0.0 or atom_coordinates.size == 0:
        return None
    atoms = np.sort(np.mod(np.asarray(atom_coordinates, dtype=float), cell_length))
    gaps = np.diff(np.append(atoms, atoms[0] + cell_length))
    index = int(np.argmax(gaps))
    start = float(atoms[index])
    end = float(atoms[(index + 1) % atoms.size])
    if index == atoms.size - 1:
        end += cell_length
    return start, end, float(gaps[index])


def estimate_vacuum_level_from_cube(
    potential_1d: np.ndarray,
    z_axis: np.ndarray,
    atom_z: np.ndarray,
    cell_length: float,
    *,
    exclude_distance: float = 3.0,
) -> dict[str, Any] | None:
    """Estimate ``V_vacuum`` directly from a planar-averaged ABACUS cube.

    The largest periodic atom-free gap is treated as vacuum.  Surface-adjacent
    regions are excluded on both sides before averaging the remaining plateau,
    which avoids assuming that the vacuum is always the last 30% of the cell.
    """
    potential_1d = np.asarray(potential_1d, dtype=float)
    z_axis = np.asarray(z_axis, dtype=float)
    gap = _largest_periodic_gap(atom_z, cell_length)
    if gap is None:
        return None
    gap_start, gap_end, gap_length = gap
    trim = max(0.0, float(exclude_distance))
    if gap_length <= 2.0 * trim:
        trim = 0.25 * gap_length
    start, end = gap_start + trim, gap_end - trim
    if end <= start:
        return None

    wrapped_z = np.mod(z_axis, cell_length)
    if end <= cell_length:
        mask = (wrapped_z >= start) & (wrapped_z <= end)
    else:
        mask = (wrapped_z >= start) | (wrapped_z <= end - cell_length)
    indices = np.where(mask)[0]
    if indices.size == 0:
        return None
    values = potential_1d[indices]
    return {
        "vacuum_level": float(np.mean(values)),
        "vacuum_std": float(np.std(values)),
        "vacuum_points": int(indices.size),
        "vacuum_coordinate": float(np.mod(0.5 * (start + end), cell_length)),
        "vacuum_gap_start": float(np.mod(gap_start, cell_length)),
        "vacuum_gap_end": float(np.mod(gap_end, cell_length)),
        "vacuum_gap_length": gap_length,
        "vacuum_exclude": trim,
    }


def macroscopic_average(potential_1d: np.ndarray, window_points: int) -> np.ndarray:
    """Apply the periodic double-window macroscopic average used by 1102.

    Periodic wrapping avoids the artificial zero-padding edge distortion of
    ``np.convolve(..., mode='same')`` at the cell boundary.
    """
    values = np.asarray(potential_1d, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("potential_1d must be a non-empty one-dimensional array")
    window = max(1, min(int(window_points), values.size))
    if window == 1:
        return values.copy()

    def running_average(arr: np.ndarray) -> np.ndarray:
        extended = np.concatenate([arr, arr[: window - 1]])
        cumulative = np.concatenate([[0.0], np.cumsum(extended)])
        averaged = (cumulative[window:] - cumulative[:-window]) / window
        # The forward window is shifted so it is centered at each grid point.
        return np.roll(averaged[: arr.size], window // 2)

    return running_average(running_average(values))


def _potential_to_ev(values: np.ndarray, path: Path) -> np.ndarray:
    """Convert ABACUS cube values to eV, preserving LOCPOT compatibility."""
    # ABACUS ElecStaticPot.cube is in Ry.  LOCPOT is conventionally already eV.
    return np.asarray(values, dtype=float) if path.name.upper().startswith("LOCPOT") else values * RY_TO_EV


# =============================================================================
# Task 1101: 1D planar average potential + work function
# =============================================================================


def _pick_path(parsed_args, args: list[str] | None, name: str, keywords: tuple[str, ...]) -> Path | None:
    """Resolve an explicit CLI path, then a positional path, if present."""
    explicit = getattr(parsed_args, name, None) if parsed_args is not None else None
    if explicit:
        return Path(explicit)
    for arg in args or []:
        path = Path(arg)
        if path.is_file() and any(keyword in path.name.upper() for keyword in keywords):
            return path
    return None


def _read_profile(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, float, tuple[int, int, int], np.ndarray, float] | None:
    """Read a cube and return profile, atom z-coordinates, and cell length."""
    result = _read_potential_cube_details(path)
    if result is None:
        return None
    potential, cell, origin, atoms = result
    nx, ny, nz = potential.shape
    cell_ang = cell * BOHR_TO_ANGSTROM
    c_length = float(np.linalg.norm(cell_ang[2]))
    z = np.arange(nz, dtype=float) * c_length / nz
    dz = c_length / nz
    profile = _potential_to_ev(np.mean(potential, axis=(0, 1)), path)
    if atoms.size:
        c_vector = cell[2]
        c_norm = float(np.linalg.norm(c_vector))
        atom_z = np.dot(atoms - origin, c_vector / c_norm) * BOHR_TO_ANGSTROM
    else:
        atom_z = np.empty(0, dtype=float)
    return z, profile, dz, (nx, ny, nz), atom_z, c_length


def _resolve_fermi(
    root: Path,
    potential_path: Path | None,
    parsed_args,
    args: list[str] | None,
) -> tuple[float | None, Path | None]:
    explicit = getattr(parsed_args, "fermi", None) if parsed_args is not None else None
    if explicit is not None:
        return float(explicit), None
    log_path = _pick_path(parsed_args, args, "log", ("LOG", "RUNNING"))
    log_path = _find_scf_log(root, potential_path, log_path)
    return (read_fermi_energy(log_path), log_path) if log_path else (None, None)


def _write_work_function_result(output_dir: Path, result: dict[str, Any]) -> None:
    import json

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "work_function.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _plot_work_function(
    output_path: Path,
    z_ang: np.ndarray,
    planar: np.ndarray,
    e_fermi: float,
    vacuum_level: float,
    title: str,
    macro: np.ndarray | None = None,
) -> None:
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_ang, planar, color="#1f77b4", linewidth=0.9, alpha=0.65, label="Planar average")
    if macro is not None:
        ax.plot(z_ang, macro, color="#d62728", linewidth=1.5, label="Macroscopic average")
    ax.axhline(y=vacuum_level, color="red", linestyle="--", linewidth=0.8,
               label=f"V_vac = {vacuum_level:.3f} eV")
    ax.axhline(y=e_fermi, color="blue", linestyle="--", linewidth=0.8,
               label=f"E_F = {e_fermi:.3f} eV")
    ax.set_xlabel("z (Å)")
    ax.set_ylabel("V (eV)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.set_xlim(z_ang.min(), z_ang.max())
    fig.savefig(output_path)
    plt.close(fig)


@task(
    1101,
    category="Work Function",
    name="Work Function",
    description="Compute work function from V_vacuum and the final SCF E_Fermi",
    cli_args=[
        {"name": "--file", "type": str, "default": None, "help": "Path to electrostatic potential cube file"},
        {"name": "--vacuum-file", "type": str, "default": None, "help": "Optional precomputed vacuum-level override"},
        {"name": "--vacuum-exclude", "type": float, "default": 3.0, "help": "Distance (Å) excluded next to slab atoms"},
        {"name": "--log", "type": str, "default": None, "help": "Path to SCF running log"},
        {"name": "--fermi", "type": float, "default": None, "help": "Override E_Fermi (eV)"},
        {"name": "--no-plot", "action": "store_true", "help": "Write data without a PNG plot"},
    ],
)
def task_work_function(
    args: list[str] | None = None,
    interactive: bool = True,
    parsed_args=None,
    output_dir: str = ".",
) -> dict[str, Any] | None:
    """Compute ``Phi = V_vacuum - E_Fermi`` for an ABACUS slab.

    The vacuum level is derived from the largest atom-free periodic gap in the
    raw cube.  ``--vacuum-file`` is retained only as an explicit override for
    users who already have an independently validated scalar reference.
    """
    console = _get_console()
    root = Path.cwd()
    output = Path(output_dir)
    pot_path = _pick_path(parsed_args, args, "file", ("POT", "ELEC", "LOCPOT")) or _find_potential_file(root)
    vacuum_path = Path(parsed_args.vacuum_file) if parsed_args and parsed_args.vacuum_file else None
    fermi, log_path = _resolve_fermi(root, pot_path, parsed_args, args)

    console.print()
    console.print("[bold cyan]=== Work Function Analysis ===[/bold cyan]")
    console.print()
    if fermi is None:
        console.print("[red]Could not find E_Fermi in an SCF log.[/red]")
        console.print("[dim]Pass --log running_scf.log or --fermi <eV>.[/dim]")
        return None
    if log_path:
        console.print(f"  [dim]SCF log: {log_path}[/dim]")
    console.print(f"  [bold]E_Fermi:[/bold] {fermi:.6f} eV")

    profile = _read_profile(pot_path) if pot_path else None
    if pot_path:
        console.print(f"  [dim]Potential file: {pot_path}[/dim]")
    if profile:
        z_ang, planar, _dz, shape, atom_z, cell_length = profile
        console.print(f"  [bold]Grid:[/bold] {shape[0]} × {shape[1]} × {shape[2]}")
    else:
        z_ang = planar = atom_z = cell_length = None

    vacuum_info = (
        estimate_vacuum_level_from_cube(
            planar, z_ang, atom_z, cell_length,
            exclude_distance=getattr(parsed_args, "vacuum_exclude", 3.0),
        )
        if profile and atom_z.size else None
    )
    vacuum_level = read_vacuum_level(vacuum_path) if vacuum_path else None
    if vacuum_info is not None:
        wf = {**vacuum_info, "work_function": vacuum_info["vacuum_level"] - fermi,
              "e_fermi": fermi, "vacuum_source": "largest atom-free gap in cube"}
    elif vacuum_level is not None:
        source = str(vacuum_path)
        wf = {"vacuum_level": vacuum_level, "work_function": vacuum_level - fermi,
              "e_fermi": fermi, "vacuum_source": source}
    elif profile:
        fallback = extract_work_function(planar, z_ang, e_fermi=fermi, vacuum_fraction=0.3)
        wf = {**fallback, "vacuum_source": "top 30% planar-average fallback"}
        console.print("  [yellow]Atom coordinates were unavailable; using top 30% planar average.[/yellow]")
    else:
        console.print("[red]No readable electrostatic-potential cube found.[/red]")
        return None

    console.print(f"  [bold green]Vacuum Level:[/bold green] {wf['vacuum_level']:.4f} eV")
    console.print(f"  [bold green]Work Function (Φ):[/bold green] {wf['work_function']:.4f} eV")
    if vacuum_path:
        console.print(f"  [dim]Vacuum source: {vacuum_path}[/dim]")

    output.mkdir(parents=True, exist_ok=True)
    if profile:
        np.savetxt(output / "planar_avg_potential.dat", np.column_stack([z_ang, planar]),
                   fmt="%.8f", header="z_Angstrom  V_eV")
        if not getattr(parsed_args, "no_plot", False):
            _plot_work_function(
                output / "work_function.png", z_ang, planar, fermi, wf["vacuum_level"],
                f"Planar-Averaged Potential — Work Function = {wf['work_function']:.3f} eV",
            )
    _write_work_function_result(output, wf)
    console.print(f"  [dim]Results saved to {output / 'work_function.json'}[/dim]")
    return wf


@task(
    1102,
    category="Work Function",
    name="Macroscopic Avg",
    description="Compute a periodic double-averaged potential and work function",
    cli_args=[
        {"name": "--file", "type": str, "default": None, "help": "Path to electrostatic potential cube file"},
        {"name": "--vacuum-file", "type": str, "default": None, "help": "Optional precomputed vacuum-level override"},
        {"name": "--vacuum-exclude", "type": float, "default": 3.0, "help": "Distance (Å) excluded next to slab atoms"},
        {"name": "--log", "type": str, "default": None, "help": "Path to SCF running log"},
        {"name": "--fermi", "type": float, "default": None, "help": "Override E_Fermi (eV)"},
        {"name": "--period", "type": float, "default": None, "help": "Averaging period in Å"},
        {"name": "--no-plot", "action": "store_true", "help": "Write data without a PNG plot"},
    ],
)
def task_macro_avg_potential(
    args: list[str] | None = None,
    interactive: bool = True,
    parsed_args=None,
    output_dir: str = ".",
) -> dict[str, Any] | None:
    """Compute the periodic double-window macroscopic potential (task 1102)."""
    console = _get_console()
    root = Path.cwd()
    pot_path = _pick_path(parsed_args, args, "file", ("POT", "ELEC", "LOCPOT")) or _find_potential_file(root)
    if pot_path is None:
        console.print("[red]No electrostatic potential cube file found.[/red]")
        return None
    profile = _read_profile(pot_path)
    if profile is None:
        console.print("[red]Could not read potential cube file.[/red]")
        return None
    z_ang, planar, dz, shape, atom_z, cell_length = profile
    fermi, log_path = _resolve_fermi(root, pot_path, parsed_args, args)
    if fermi is None:
        console.print("[red]Could not find E_Fermi in an SCF log.[/red]")
        console.print("[dim]Pass --log running_scf.log or --fermi <eV>.[/dim]")
        return None

    period = getattr(parsed_args, "period", None) if parsed_args is not None else None
    if period is None and interactive:
        from rich.prompt import Prompt
        period = float(Prompt.ask("  Averaging period (Å) — typically inter-layer spacing",
                                  default=f"{cell_length / 4:.2f}"))
    period = float(period if period is not None else cell_length / 4)
    if period <= 0.0:
        raise ValueError("Averaging period must be positive")
    window_points = max(1, int(round(period / dz)))
    macro = macroscopic_average(planar, window_points)

    vacuum_path = Path(parsed_args.vacuum_file) if parsed_args and parsed_args.vacuum_file else None
    vacuum_info = (
        estimate_vacuum_level_from_cube(
            macro, z_ang, atom_z, cell_length,
            exclude_distance=getattr(parsed_args, "vacuum_exclude", 3.0),
        )
        if atom_z.size else None
    )
    vacuum_level = read_vacuum_level(vacuum_path) if vacuum_path else None
    if vacuum_info is not None:
        wf = {**vacuum_info, "work_function": vacuum_info["vacuum_level"] - fermi,
              "e_fermi": fermi, "vacuum_source": "largest atom-free gap in macro profile"}
    elif vacuum_level is None:
        wf = extract_work_function(macro, z_ang, e_fermi=fermi, vacuum_fraction=0.3)
        wf["vacuum_source"] = "top 30% macroscopic-average fallback (no atom coordinates)"
    else:
        wf = {"vacuum_level": vacuum_level, "work_function": vacuum_level - fermi,
              "e_fermi": fermi, "vacuum_source": str(vacuum_path)}

    console.print()
    console.print("[bold cyan]=== Macroscopic-Averaged Potential ===[/bold cyan]")
    console.print(f"  [bold]Grid:[/bold] {shape[0]} × {shape[1]} × {shape[2]}")
    console.print(f"  [bold]Averaging period:[/bold] {period:.4f} Å ({window_points} points)")
    console.print(f"  [bold]E_Fermi:[/bold] {fermi:.6f} eV")
    console.print(f"  [bold green]Vacuum Level:[/bold green] {wf['vacuum_level']:.4f} eV")
    console.print(f"  [bold green]Work Function (Φ):[/bold green] {wf['work_function']:.4f} eV")
    if log_path:
        console.print(f"  [dim]SCF log: {log_path}[/dim]")

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "macro_avg_potential.dat", np.column_stack([z_ang, macro]),
               fmt="%.8f", header="z_Angstrom  V_eV")
    wf.update({"period_angstrom": period, "window_points": window_points})
    if not getattr(parsed_args, "no_plot", False):
        _plot_work_function(
            output / "macro_avg_potential.png", z_ang, planar, fermi, wf["vacuum_level"],
            f"Macroscopic-Averaged Potential — Φ = {wf['work_function']:.3f} eV", macro=macro,
        )
    _write_work_function_result(output, wf)
    console.print(f"  [dim]Results saved to {output / 'work_function.json'}[/dim]")
    return wf
