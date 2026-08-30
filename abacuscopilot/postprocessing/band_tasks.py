"""Band structure postprocessing tasks for ABACUS.

Task IDs 721-729

Plots band structures, fat-bands (projected bands), and extracts
band gap information from ABACUS band structure output.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import task


def _find_bands_file() -> Path | None:
    """Find the BANDS file in the working directory.

    Handles all ABACUS versions:
    - PW / older LCAO: ``BANDS_*.dat``
    - LCAO GPU v3.x:   ``band.txt``
    """
    candidates = (list(Path().glob("band.txt")) +
                  list(Path().glob("band_*.dat")) +
                  list(Path().glob("BANDS*.dat")) +
                  list(Path().glob("BANDS*")) +
                  list(Path().glob("OUT.*/band.txt")) +
                  list(Path().glob("OUT.*/band_*.dat")) +
                  list(Path().glob("OUT.*/BANDS*.dat")) +
                  list(Path().glob("OUT.*/BANDS*")))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _find_kpt_file() -> Path | None:
    """Find KPT file for band path labels."""
    # Current dir first, then OUT.* subdirs
    for name in ("KPT", "KPT_BAND", "KPT_LINE"):
        p = Path(name)
        if p.exists():
            return p
    for d in Path().glob("OUT.*"):
        for name in ("KPT", "KPT_BAND", "KPT_LINE"):
            p = d / name
            if p.exists():
                return p
    # Also check parent dir (KPT is usually in the working dir)
    for name in ("KPT", "KPT_BAND", "KPT_LINE"):
        p = Path("..") / name
        if p.exists():
            return p
    return None


def _extract_band_gap(energies: np.ndarray, e_fermi: float,
                      occupation: np.ndarray | None = None) -> dict[str, Any]:
    """Determine the band gap from band energies.

    For insulators/semiconductors: find the VBM and CBM.
    For metals: report zero gap.

    Args:
        energies: (nbands, nkpoints) band energies in eV.
        e_fermi: Fermi energy in eV.
        occupation: Optional occupation numbers (nbands, nkpoints).

    Returns:
        Dict with 'vbm', 'cbm', 'gap', 'is_metal', 'vbm_kpt', 'cbm_kpt'.
    """
    nbands, nkpts = energies.shape
    energies_shifted = energies - e_fermi

    # Find valence and conduction band edges
    # VBM: highest occupied band (energy <= 0)
    # CBM: lowest unoccupied band (energy >= 0)

    # For each k-point, find the highest band below E_F and lowest band above E_F
    vbm = -np.inf
    cbm = np.inf
    vbm_kpt = 0
    cbm_kpt = 0
    vbm_band = 0
    cbm_band = 0

    for ik in range(nkpts):
        for ib in range(nbands):
            e = energies_shifted[ib, ik]
            if e <= 0 and e > vbm:
                vbm = e
                vbm_kpt = ik
                vbm_band = ib
            if e >= 0 and e < cbm:
                cbm = e
                cbm_kpt = ik
                cbm_band = ib

    # Determine if metal (bands crossing E_F)
    is_metal = False
    for ik in range(nkpts):
        below = any(energies_shifted[ib, ik] <= 0 for ib in range(nbands))
        above = any(energies_shifted[ib, ik] > 0 for ib in range(nbands))
        # Check if any band crosses E_F (one band with both negative and positive)
        for ib in range(nbands):
            # Simple check: is any band very close to E_F?
            if abs(energies_shifted[ib, ik]) < 0.01:
                is_metal = True
                break

    gap = max(0.0, cbm - vbm) if not is_metal else 0.0

    return {
        "vbm": vbm + e_fermi,
        "cbm": cbm + e_fermi,
        "gap": gap,
        "vbm_shifted": vbm,
        "cbm_shifted": cbm,
        "is_metal": is_metal or gap < 0.01,
        "vbm_kpt": vbm_kpt,
        "cbm_kpt": cbm_kpt,
        "vbm_band": vbm_band,
        "cbm_band": cbm_band,
    }


# =============================================================================
# Task 721: Plot band structure
# =============================================================================

@task(801, category="Band Structure", name="Plot Band Structure",
      description="Plot electronic band structure from BANDS_*.dat",
      cli_args=[
          {"name": "--bands", "type": str, "default": None, "help": "Path to BANDS file"},
          {"name": "--kpt", "type": str, "default": None, "help": "Path to KPT file"},
          {"name": "--erange", "type": str, "default": None, "help": "Energy range (min,max)"},
      ])
def task_plot_bands(args: list[str] | None = None, interactive: bool = True,
                    parsed_args=None) -> None:
    """Plot band structure from ABACUS output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Plot Band Structure ===[/bold cyan]")
    console.print()

    # Find bands file
    bands_path = None
    if parsed_args and parsed_args.bands:
        bands_path = Path(parsed_args.bands)
    elif args:
        for arg in args:
            if "BANDS" in arg and Path(arg).exists():
                bands_path = Path(arg)
                break

    if bands_path is None:
        bands_path = _find_bands_file()

    if bands_path is None:
        console.print("[red]No BANDS_*.dat file found.[/red]")
        console.print("[dim]Run an ABACUS NSCF calculation with out_band=True first.[/dim]")
        return

    console.print(f"  [dim]Bands file: {bands_path}[/dim]")

    # Find KPT
    kpt_path = None
    if parsed_args and parsed_args.kpt:
        kpt_path = Path(parsed_args.kpt)
    if kpt_path is None:
        kpt_path = _find_kpt_file()
    if kpt_path:
        console.print(f"  [dim]KPT file: {kpt_path}[/dim]")

    # Read data
    from abacuscopilot.plotting.bands import read_bands_dat
    k_dists, energies, e_fermi = read_bands_dat(bands_path)

    if energies.size == 0:
        console.print("[red]No band data could be read.[/red]")
        return

    nbands, nkpts = energies.shape
    console.print(f"  [bold]Bands:[/bold] {nbands}")
    console.print(f"  [bold]K-points:[/bold] {nkpts}")
    console.print(f"  [bold]E_Fermi:[/bold] {e_fermi:.6f} eV")

    # Band gap analysis
    gap_info = _extract_band_gap(energies, e_fermi)
    console.print()
    if gap_info["is_metal"]:
        console.print("  [bold]Character:[/bold] [yellow]Metallic[/yellow] (no gap)")
    else:
        console.print(f"  [bold]VBM:[/bold] {gap_info['vbm']:.6f} eV")
        console.print(f"  [bold]CBM:[/bold] {gap_info['cbm']:.6f} eV")
        console.print(f"  [bold]Band Gap:[/bold] [green]{gap_info['gap']:.4f} eV[/green]")
        if gap_info["vbm_shifted"] != gap_info["cbm_shifted"]:
            kpt_str = "Γ" if gap_info["vbm_kpt"] == 0 else str(gap_info["vbm_kpt"])
            console.print(f"  [bold]Gap type:[/bold] {'Direct' if gap_info['vbm_kpt'] == gap_info['cbm_kpt'] else 'Indirect'}")

    # Parse energy range from CLI or interactive prompt
    e_range = None
    if parsed_args and parsed_args.erange:
        parts = parsed_args.erange.split(",")
        if len(parts) == 2:
            e_range = (float(parts[0]), float(parts[1]))

    if interactive:
        from rich.prompt import Prompt

        custom_range = Prompt.ask(
            "\n  Custom energy range? (e.g., '-5,5' or Enter for auto)",
            default=""
        )
        if custom_range:
            parts = custom_range.split(",")
            if len(parts) == 2:
                e_range = (float(parts[0]), float(parts[1]))

    # Plot
    from abacuscopilot.plotting.bands import plot_bands
    save_name = "BAND.png"
    plot_bands(
        bands_path,
        kpt=kpt_path,
        e_range=e_range,
        save=save_name,
    )
    console.print(f"  [green]✓ Plot saved to {save_name}[/green]")

    console.print()


# =============================================================================
# Task 722: Fat-band plot
# =============================================================================

@task(802, category="Band Structure", name="Fat-band Plot",
      description="Plot fat-band (projected band) structure")
def task_fatbands(args: list[str] | None = None, interactive: bool = True) -> None:
    """Plot projected (fat) band structure."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Fat-band Plot ===[/bold cyan]")
    console.print()

    bands_path = _find_bands_file()
    if bands_path is None:
        console.print("[red]No BANDS file found.[/red]")
        return

    # Find projected bands file
    proj_path = None
    proj_candidates = (list(Path().glob("pbands*.xml")) +
                      list(Path().glob("PBANDS_1")) +
                      list(Path().glob("PBANDS_2")) +
                      list(Path().glob("ProjBands*.dat")) +
                      list(Path().glob("OUT.*/pbands*.xml")) +
                      list(Path().glob("OUT.*/PBANDS_1")) +
                      list(Path().glob("OUT.*/PBANDS_2")) +
                      list(Path().glob("OUT.*/ProjBands*.dat")) +
                      list(Path().glob("OUT.*/PBANDS*.dat")))

    if proj_candidates:
        proj_path = proj_candidates[0]
        console.print(f"  [dim]Projected bands: {proj_path}[/dim]")
    else:
        console.print("[red]No projected band data found.[/red]")
        console.print("[yellow]Fat-band plot requires projected band data.[/yellow]")
        console.print("[yellow]Add 'out_proj_band 1' to INPUT and re-run the NSCF band calculation,[/yellow]")
        console.print("[yellow]then re-run this task in the output directory.[/yellow]")
        return

    kpt_path = _find_kpt_file()

    if interactive:
        from rich.prompt import Prompt
        custom_range = Prompt.ask(
            "  Custom energy range? (e.g., '-5,5' or Enter for auto)",
            default=""
        )
        e_range = None
        if custom_range:
            parts = custom_range.split(",")
            if len(parts) == 2:
                e_range = (float(parts[0]), float(parts[1]))
    else:
        e_range = None

    from abacuscopilot.plotting.bands import plot_fatbands
    save_name = "PBAND.png"
    plot_fatbands(
        bands_path,
        proj_path or bands_path,
        kpt=kpt_path,
        e_range=e_range,
        save=save_name,
    )
    console.print(f"  [green]✓ Fat-band plot saved to {save_name}[/green]")
    console.print()
