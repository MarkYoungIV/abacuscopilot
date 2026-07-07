"""Density of states postprocessing tasks for ABACUS.

Task IDs 731-739

Plots total DOS, projected DOS (PDOS), and combined DOS+band plots
from ABACUS DOS output files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import task


def _find_dos_file() -> Path | None:
    """Find DOS file in working directory."""
    # v3.7+: DOS1_smearing.dat, DOS1; older: DOS*.dat
    candidates = (list(Path().glob("DOS1_smearing.dat")) +
                  list(Path().glob("DOS1")) +
                  list(Path().glob("DOS*.dat")) +
                  list(Path().glob("OUT.*/DOS1_smearing.dat")) +
                  list(Path().glob("OUT.*/DOS1")) +
                  list(Path().glob("OUT.*/DOS*.dat")))
    candidates = [c for c in candidates if "PDOS" not in c.name.upper()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _find_pdos_file() -> Path | None:
    """Find PDOS file in working directory."""
    # New ABACUS v3.7+: PDOS (no extension); older: PDOS_*.dat
    candidates = (list(Path().glob("PDOS")) +
                  list(Path().glob("PDOS_*.dat")) +
                  list(Path().glob("OUT.*/PDOS")) +
                  list(Path().glob("OUT.*/PDOS_*.dat")))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


# =============================================================================
# Task 731: Plot DOS
# =============================================================================

@task(901, category="DOS/PDOS", name="Plot DOS",
      description="Plot total density of states from DOS output",
      cli_args=[
          {"name": "--dos", "type": str, "default": None, "help": "Path to DOS file"},
          {"name": "--erange", "type": str, "default": None, "help": "Energy range (min,max)"},
          {"name": "--fill", "type": str, "default": "y", "help": "Fill under curve (y/n)"},
      ])
def task_plot_dos(args: list[str] | None = None, interactive: bool = True,
                  parsed_args=None) -> None:
    """Plot total density of states from ABACUS output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Plot DOS ===[/bold cyan]")
    console.print()

    dos_path = None
    if parsed_args and parsed_args.dos:
        dos_path = Path(parsed_args.dos)
    if dos_path is None and args:
        for arg in args:
            if "DOS" in arg.upper() and Path(arg).exists():
                dos_path = Path(arg)
                break

    if dos_path is None:
        dos_path = _find_dos_file()

    if dos_path is None:
        console.print("[red]No DOS_*.dat file found.[/red]")
        console.print("[dim]Run an ABACUS NSCF calculation with out_dos=1 or out_dos=2.[/dim]")
        return

    console.print(f"  [dim]DOS file: {dos_path}[/dim]")

    from abacuscopilot.plotting.dos import read_dos_dat
    energies, dos, e_fermi = read_dos_dat(dos_path)

    if energies.size == 0:
        console.print("[red]No DOS data could be read.[/red]")
        return

    console.print(f"  [bold]Energy points:[/bold] {len(energies)}")
    console.print(f"  [bold]E_Fermi:[/bold] {e_fermi:.6f} eV")
    console.print(f"  [bold]DOS at E_F:[/bold] {dos[np.argmin(np.abs(energies - e_fermi))]:.4f} states/eV")

    # Parse energy range and fill from CLI or interactive prompt
    e_range = None
    if parsed_args and parsed_args.erange:
        parts = parsed_args.erange.split(",")
        if len(parts) == 2:
            e_range = (float(parts[0]), float(parts[1]))

    fill = True
    if parsed_args:
        fill = parsed_args.fill.lower() in ("y", "yes", "true", "1")

    if interactive:
        custom_range = console.input(
            "  Custom energy range? (e.g., '-10,10' or Enter for auto): "
        ).strip()
        if custom_range:
            parts = custom_range.split(",")
            if len(parts) == 2:
                try:
                    e_range = (float(parts[0]), float(parts[1]))
                except ValueError:
                    pass

        fill = console.input("  Fill under curve? (y/n) [y]: ").strip().lower() != "n"

    # Plot (both modes)
    from abacuscopilot.plotting.dos import plot_dos
    save_name = f"{dos_path.stem}.png"
    plot_dos(
        dos_path,
        e_range=e_range,
        fill=fill,
        save=save_name,
    )
    console.print(f"  [green]✓ DOS plot saved to {save_name}[/green]")

    console.print()


# =============================================================================
# Task 732: Plot PDOS
# =============================================================================

@task(902, category="DOS/PDOS", name="Plot PDOS",
      description="Plot projected density of states (orbital-resolved)")
def task_plot_pdos(args: list[str] | None = None, interactive: bool = True) -> None:
    """Plot projected density of states."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Plot Projected DOS ===[/bold cyan]")
    console.print()

    pdos_path = None
    if args:
        for arg in args:
            if "PDOS" in arg.upper() and Path(arg).exists():
                pdos_path = Path(arg)
                break

    if pdos_path is None:
        pdos_path = _find_pdos_file()

    if pdos_path is None:
        console.print("[red]No PDOS_*.dat file found.[/red]")
        console.print("[dim]Run an ABACUS NSCF calculation with out_dos=2.[/dim]")
        return

    console.print(f"  [dim]PDOS file: {pdos_path}[/dim]")

    from abacuscopilot.plotting.dos import read_pdos_dat
    data = read_pdos_dat(pdos_path)

    if "energies" not in data or data["energies"].size == 0:
        console.print("[red]No PDOS data could be read.[/red]")
        return

    e_fermi = data.get("e_fermi", 0.0)
    energies = data["energies"]

    console.print(f"  [bold]Energy points:[/bold] {len(energies)}")
    console.print(f"  [bold]E_Fermi:[/bold] {e_fermi:.6f} eV")

    # Show available orbitals
    orbitals = [k for k in data if k not in ("energies", "e_fermi")
                and isinstance(data[k], np.ndarray)
                and data[k].shape == energies.shape]
    console.print(f"  [bold]Available orbitals:[/bold] {', '.join(orbitals)}")

    if interactive:
        # Numbered multi-select for orbitals
        console.print()
        console.print("  [bold]Select orbitals to plot:[/bold]")
        console.print(f"    0. All ({len(orbitals)} orbitals)")
        for i, orb in enumerate(orbitals, 1):
            console.print(f"    {i:>2}. {orb}")
        console.print()
        sel = console.input("  Enter numbers (comma/space separated, 0=all) [0]: ").strip()
        if not sel or "0" in sel:
            selected_orbitals = orbitals
        else:
            # Parse numbers
            import re
            nums = [int(x) for x in re.split(r"[,\s]+", sel) if x.isdigit()]
            selected_orbitals = [orbitals[i - 1] for i in nums if 1 <= i <= len(orbitals)]
            if not selected_orbitals:
                selected_orbitals = orbitals

        e_range = None
        custom_range = console.input(
            "  Custom energy range? (e.g., '-10,10' or Enter for auto): "
        ).strip()
        if custom_range:
            parts = custom_range.split(",")
            if len(parts) == 2:
                try:
                    e_range = (float(parts[0]), float(parts[1]))
                except ValueError:
                    pass

        from abacuscopilot.plotting.dos import plot_pdos
        save_name = "PDOS.png"
        plot_pdos(
            pdos_path,
            orbitals=selected_orbitals,
            e_range=e_range,
            save=save_name,
        )
        console.print(f"  [green]✓ PDOS plot saved to {save_name}[/green]")

    console.print()


# =============================================================================
# Task 733: Combined DOS + Bands
# =============================================================================

@task(903, category="DOS/PDOS", name="DOS+Bands Combined",
      description="Combined band structure and DOS plot side by side")
def task_dos_bands_combined(args: list[str] | None = None, interactive: bool = True) -> None:
    """Create a combined band structure + DOS figure."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Combined DOS + Band Structure ===[/bold cyan]")
    console.print()

    # Find files
    from abacuscopilot.postprocessing.band_tasks import _find_bands_file, _find_kpt_file
    bands_path = _find_bands_file()
    dos_path = _find_dos_file()
    kpt_path = _find_kpt_file()

    if bands_path is None:
        console.print("[red]No BANDS file found.[/red]")
        return
    if dos_path is None:
        console.print("[red]No DOS file found.[/red]")
        return

    console.print(f"  [dim]Bands: {bands_path}[/dim]")
    console.print(f"  [dim]DOS: {dos_path}[/dim]")
    if kpt_path:
        console.print(f"  [dim]KPT: {kpt_path}[/dim]")

    from abacuscopilot.plotting.dos import plot_dos_bands_combined
    save_name = "DOS_bands_combined.png"
    plot_dos_bands_combined(
        bands_path,
        dos_path,
        kpt=kpt_path,
        save=save_name,
    )
    console.print(f"  [green]✓ Combined plot saved to {save_name}[/green]")
    console.print()
