"""Work function analysis tasks for ABACUS output.

Task IDs 1101-1102 (currently HIDDEN — functionality not yet verified;
tasks 1101/1102 are disabled until validated, see the @task comments below).

Reads electrostatic potential from OUT.ABACUS, computes 1D planar
average along the vacuum direction, and extracts the work function.

Method:
1. Read LOCPOT / electrostatic potential cube file
2. Compute 1D planar average along vacuum direction (typically z)
3. Identify vacuum level from the plateau region
4. Compute work function:  Φ = V_vacuum − E_Fermi
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console

# =============================================================================
# Potential file finder
# =============================================================================


def _find_potential_file() -> Path | None:
    """Find electrostatic potential file in the working directory."""
    candidates = (
        list(Path().glob("ElecStaticPot*.cube")) +
        list(Path().glob("LOCPOT*")) +
        list(Path().glob("POT.cube")) +
        list(Path().glob("OUT.*/ElecStaticPot*.cube")) +
        list(Path().glob("OUT.*/LOCPOT*")) +
        list(Path().glob("OUT.*/POT.cube"))
    )
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


# =============================================================================
# Cube potential reader
# =============================================================================


def read_potential_cube(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read electrostatic potential from a Cube-format file.

    Returns (data_3d, cell_bohr, origin) or None on failure.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip() if not l.startswith("#")]

    if len(lines) < 6:
        return None

    # Skip 2 comment lines
    data_start = 2

    # Origin line: natom x0 y0 z0
    parts = lines[data_start].split()
    natom = int(parts[0])
    origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
    data_start += 1

    # Grid dimensions: nx dx dy dz
    parts = lines[data_start].split()
    nx = int(parts[0])
    data_start += 1

    cell = np.zeros((3, 3))
    for i in range(3):
        parts = lines[data_start + i].split()
        n = int(parts[0])
        cell[i] = np.array([float(p) * n for p in parts[1:4]])
    data_start += 3

    # Skip atom positions
    data_start += natom

    if data_start >= len(lines):
        return None

    # Read ny, nz
    ny = nx
    nz = nx
    if data_start < len(lines):
        parts = lines[data_start].split()
        if len(parts) == 1:
            ny = int(parts[0])
            data_start += 1
    if data_start < len(lines):
        parts = lines[data_start].split()
        if len(parts) == 1:
            nz = int(parts[0])
            data_start += 1

    # Read volumetric data
    raw = []
    for line in lines[data_start:]:
        raw.extend(float(x) for x in line.split())

    expected = nx * ny * nz
    data = np.array(raw[:expected]).reshape(nz, ny, nx)

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
    npts = len(potential_1d)
    vac_start = int(npts * (1 - vacuum_fraction))

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


# =============================================================================
# Task 1101: 1D planar average potential + work function
#
# HIDDEN (2026-08-31): functionality not yet verified — availability pending.
# Re-enable by un-commenting the @task decorator below.
# =============================================================================


# @task(1101, category="Work Function", name="Work Function",
#       description="Compute 1D planar-averaged electrostatic potential and extract work function",
#       cli_args=[
#           {"name": "--file", "type": str, "default": None, "help": "Path to electrostatic potential cube file"},
#       ])
def task_work_function(args: list[str] | None = None, interactive: bool = True,
                       parsed_args=None) -> None:
    """Compute work function from electrostatic potential.

    Reads the electrostatic potential cube file, computes a 1D planar
    average along z, and extracts the vacuum level to determine the
    work function: Φ = V_vacuum − E_Fermi.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Work Function Analysis ===[/bold cyan]")
    console.print()

    # Find potential file
    pot_path = None
    if parsed_args and parsed_args.file:
        pot_path = Path(parsed_args.file)
    if pot_path is None and args:
        for arg in args:
            if "POT" in arg.upper() or "ELEC" in arg.upper():
                p = Path(arg)
                if p.exists():
                    pot_path = p
                    break
    if pot_path is None:
        pot_path = _find_potential_file()

    if pot_path is None:
        console.print("[red]No electrostatic potential file found.[/red]")
        console.print("[dim]Run an ABACUS SCF calculation with out_pot=2.[/dim]")
        console.print("[dim]Look for files like: ElecStaticPot.cube or POT.cube in OUT.ABACUS/[/dim]")
        return

    console.print(f"  [dim]Potential file: {pot_path}[/dim]")

    # Read potential
    result = read_potential_cube(pot_path)
    if result is None:
        console.print("[red]Could not read potential cube file.[/red]")
        return

    potential_3d, cell, origin = result
    nz, ny, nx = potential_3d.shape
    console.print(f"  [bold]Grid:[/bold] {nx} × {ny} × {nz}")

    # Convert cell from Bohr to Angstrom for display
    from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
    cell_ang = cell * BOHR_TO_ANGSTROM
    c_length = np.linalg.norm(cell_ang[2])
    console.print(f"  [bold]Cell c-axis:[/bold] {c_length:.4f} Å")

    # 1D planar average along z (axis=0 in z,y,x ordering)
    avg_pot = np.mean(potential_3d, axis=(1, 2))  # average over x and y
    z_coords = origin[2] + np.linspace(0, cell[2, 2], nz, endpoint=False)
    z_ang = z_coords * BOHR_TO_ANGSTROM

    # Get E_Fermi
    e_fermi = 0.0
    # Try reading from INPUT or output
    try:
        from abacuscopilot.io.input_file import read_input
        params = read_input("INPUT")
        # Look for scf_thr in extras (some versions store E_Fermi)
    except Exception:
        pass

    # Try reading from running_scf.log for e_fermi
    log_candidates = (list(Path().glob("running_*.log")) +
                     list(Path().glob("OUT.*/running_scf.log")))
    for log_path in log_candidates:
        try:
            content = log_path.read_text()
            import re
            m = re.search(r"(?:EFERMI|E-fermi|Fermi)\s*[=:]\s*([\d.-]+)", content, re.IGNORECASE)
            if m:
                e_fermi = float(m.group(1))
                break
        except Exception:
            pass

    console.print(f"  [bold]E_Fermi:[/bold] {e_fermi:.6f} eV" if e_fermi != 0.0 else "  [yellow]E_Fermi: 0.0 (not found, using 0)[/yellow]")

    # Extract work function
    wf = extract_work_function(avg_pot, z_ang, e_fermi=e_fermi, vacuum_fraction=0.3)

    console.print()
    console.print(f"  [bold green]Vacuum Level:[/bold green] {wf['vacuum_level']:.4f} eV")
    console.print(f"  [bold green]Work Function (Φ):[/bold green] {wf['work_function']:.4f} eV")
    console.print(f"  Vacuum region: z > {wf['z_vacuum_start']:.2f} Å")

    # Plot
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_ang, avg_pot, color="#1f77b4", linewidth=1.5)
    ax.axhline(y=wf["vacuum_level"], color="red", linestyle="--", linewidth=0.8,
               label=f"V_vac = {wf['vacuum_level']:.3f} eV")
    ax.axhline(y=e_fermi, color="blue", linestyle="--", linewidth=0.8,
               label=f"E_F = {e_fermi:.3f} eV")
    ax.axvline(x=wf["z_vacuum_start"], color="gray", linestyle=":", linewidth=0.5,
               label=f"Vacuum start ({wf['z_vacuum_start']:.1f} Å)")
    ax.set_xlabel("z (Å)")
    ax.set_ylabel("V (eV)")
    ax.set_title(f"Planar-Averaged Potential — Work Function = {wf['work_function']:.3f} eV")
    ax.legend(loc="upper right")
    ax.set_xlim(z_ang.min(), z_ang.max())

    # Save data
    data_file = "planar_avg_potential.dat"
    np.savetxt(data_file, np.column_stack([z_ang, avg_pot]),
               fmt="%.8f", header="z_Angstrom  V_eV")
    console.print(f"  [dim]Data saved to {data_file}[/dim]")

    save_name = "work_function.png"
    fig.savefig(save_name)
    console.print(f"  [green]✓ Plot saved to {save_name}[/green]")
    plt.close(fig)
    console.print()


# =============================================================================
# Task 1102: Macroscopic average (double average for semiconductor slabs)
#
# HIDDEN (2026-08-31): functionality not yet verified — availability pending.
# Re-enable by un-commenting the @task decorator below.
# =============================================================================


# @task(1102, category="Work Function", name="Macroscopic Avg",
#       description="Compute macroscopic-averaged potential (double filter) for better vacuum level")
def task_macro_avg_potential(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compute macroscopic-averaged electrostatic potential.

    Uses a double-window averaging technique to remove atomic oscillations,
    yielding a smoother potential profile for more accurate work function
    extraction.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Macroscopic-Averaged Potential ===[/bold cyan]")
    console.print()

    pot_path = None
    if args:
        for arg in args:
            if "POT" in arg.upper() or "ELEC" in arg.upper():
                p = Path(arg)
                if p.exists():
                    pot_path = p
                    break
    if pot_path is None:
        pot_path = _find_potential_file()

    if pot_path is None:
        console.print("[red]No electrostatic potential file found.[/red]")
        return

    result = read_potential_cube(pot_path)
    if result is None:
        console.print("[red]Could not read potential file.[/red]")
        return

    potential_3d, cell, origin = result
    nz = potential_3d.shape[0]

    from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
    cell_ang = cell * BOHR_TO_ANGSTROM
    c_len = np.linalg.norm(cell_ang[2])

    # 1D planar average
    avg_pot = np.mean(potential_3d, axis=(1, 2))
    z_coords = origin[2] + np.linspace(0, cell[2, 2], nz, endpoint=False)
    z_ang = z_coords * BOHR_TO_ANGSTROM
    dz = c_len / nz

    # Ask for period
    if interactive:
        from rich.prompt import Prompt
        period_str = Prompt.ask("  Averaging period (Å) — typically inter-layer spacing", default=f"{c_len/4:.2f}")
        period = float(period_str)
    else:
        period = c_len / 4  # rough estimate

    n_period = max(1, int(period / dz))

    # Double running average (macroscopic average)
    def running_avg(arr: np.ndarray, n: int) -> np.ndarray:
        return np.convolve(arr, np.ones(n) / n, mode="same")

    macro_pot = running_avg(running_avg(avg_pot, n_period), n_period)

    # Extract work function
    e_fermi = 0.0
    wf = extract_work_function(macro_pot, z_ang, e_fermi=e_fermi, vacuum_fraction=0.3)

    console.print()
    console.print(f"  [bold green]Vacuum Level:[/bold green] {wf['vacuum_level']:.4f} eV")
    console.print(f"  [bold green]Work Function (Φ):[/bold green] {wf['work_function']:.4f} eV")

    # Plot
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_ang, avg_pot, color="#1f77b4", linewidth=0.8, alpha=0.5, label="Planar average")
    ax.plot(z_ang, macro_pot, color="#d62728", linewidth=1.8, label="Macroscopic average")
    ax.axhline(y=wf["vacuum_level"], color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("z (Å)")
    ax.set_ylabel("V (eV)")
    ax.set_title(f"Macroscopic-Averaged Potential — Φ = {wf['work_function']:.3f} eV")
    ax.legend(loc="upper right")

    save_name = "macro_avg_potential.png"
    fig.savefig(save_name)
    console.print(f"  [green]✓ Plot saved to {save_name}[/green]")
    plt.close(fig)
    console.print()
