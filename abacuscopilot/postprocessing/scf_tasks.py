"""SCF convergence analysis tasks for ABACUS output.

Task IDs 711-719

Analyzes SCF convergence from OUT.ABACUS/running_*.log files,
extracts energy convergence history, and diagnoses convergence issues.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import task


def parse_scf_log(filepath: str | Path) -> dict[str, Any]:
    """Parse SCF convergence data from an ABACUS output file.

    Reads the running log or OUT.ABACUS/output file to extract:
    - Total energy per SCF step
    - Energy difference between steps
    - Charge density difference
    - Final convergence status

    Args:
        filepath: Path to OUT.ABACUS/running_scf.log or similar.

    Returns:
        Dict with keys: 'energies', 'ediffs', 'drho', 'converged', 'nsteps',
        'final_energy'.
    """
    filepath = Path(filepath)
    result = {
        "energies": [],
        "ediffs": [],
        "drho": [],
        "converged": False,
        "nsteps": 0,
        "final_energy": 0.0,
        "natom": 0,
    }

    if not filepath.exists():
        return result

    # Energy patterns — each format is handled separately to avoid unit mixing.
    # E_KohnSham: "E_KohnSham     Ry_value     eV_value" — capture first number (Ry)
    _eks_re = re.compile(r"E_KohnSham\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)", re.IGNORECASE)
    # CU format: "CU1     eV_value     dE     drho     time" — eV
    _cu_re = re.compile(r"CU\d+\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)", re.IGNORECASE)
    # Generic: "final etot is ... eV", "total energy = ... Ry", etc.
    _generic_re = re.compile(
        r"(?:final\s+etot\s+is|total energy|ETOT|ENERGY)\s*[=:\s]+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    # Convergence: "charge density convergence" or "SCF convergence"
    drho_re = re.compile(
        r"(?:drho|charge density)\s*[=:]\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    # Ediff pattern
    ediff_re = re.compile(
        r"(?:delta E|dE|energy diff)\s*[=:]\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    # Convergence markers differ between ABACUS builds:
    #   CPU build: "charge density convergence is achieved"
    #   GPU build: "#SCF IS CONVERGED#"
    # (NOTE: the old pattern "converge[nc]ed" was a typo — it matched
    #  "convergned"/"convergced", never the real word "converged", so GPU
    #  logs were wrongly reported as not converged.)
    converged_re = re.compile(
        r"converged|convergence\s+is\s+achieved|SCF\s+is\s+converged|SCF\s+done|reach",
        re.IGNORECASE,
    )

    with open(filepath, errors="ignore") as f:
        content = f.read()

    # Extract energies — try formats in priority order, never mix units
    raw_matches = _eks_re.findall(content)       # E_KohnSham → Ry
    if raw_matches:
        is_ev = False
    else:
        raw_matches = _cu_re.findall(content)     # CU\d+ → eV
        if raw_matches:
            is_ev = True
        else:
            raw_matches = _generic_re.findall(content)  # fallback
            is_ev = bool(re.search(r"final\s+etot\s+is", content, re.IGNORECASE))

    if raw_matches:
        result["energies"] = [float(m) for m in raw_matches]
        result["final_energy"] = result["energies"][-1]
        result["final_energy_is_ev"] = is_ev  # flag for downstream display
        result["nsteps"] = len(result["energies"])

    # Extract drho
    drhos = [float(m) for m in drho_re.findall(content)]
    if drhos:
        result["drho"] = drhos

    # Compute ediffs from energies
    energ = result["energies"]
    if len(energ) > 1:
        result["ediffs"] = [
            abs(energ[i] - energ[i - 1]) for i in range(1, len(energ))
        ]

    # Check convergence
    if converged_re.search(content):
        result["converged"] = True

    # Atom count — authoritative, straight from ABACUS
    natom_re = re.compile(r"TOTAL\s+ATOM\s+NUMBER\s*[=:]\s*(\d+)", re.IGNORECASE)
    natom_match = natom_re.search(content)
    if natom_match:
        result["natom"] = int(natom_match.group(1))

    return result


def _find_latest_scf_log() -> Path | None:
    """Find the most recent SCF log file in the working directory."""
    candidates = []

    for pattern in [
        "running_*.log",             # current dir (already inside OUT.*)
        "OUT.ABACUS/running_scf.log",
        "OUT.ABACUS/running*.log",
        "OUT.*/running_scf.log",
        "OUT.*/running*.log",
        "output.log",
    ]:
        matches = list(Path().glob(pattern))
        candidates.extend(matches)

    if not candidates:
        return None

    # Return the most recently modified
    return max(candidates, key=lambda p: p.stat().st_mtime)


# =============================================================================
# Task 711: SCF convergence check
# =============================================================================

@task(711, category="SCF Analysis", name="SCF Convergence",
      description="Check and visualize SCF convergence from OUT.ABACUS",
      cli_args=[
          {"name": "--log", "type": str, "default": None, "help": "Path to SCF log file"},
      ])
def task_scf_convergence(args: list[str] | None = None, interactive: bool = True,
                         parsed_args=None) -> None:
    """Analyze SCF convergence from ABACUS output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== SCF Convergence Analysis ===[/bold cyan]")
    console.print()

    # Find SCF log
    log_path = None
    if parsed_args and parsed_args.log:
        log_path = Path(parsed_args.log)
    if log_path is None and args:
        for arg in args:
            if Path(arg).exists():
                log_path = Path(arg)
                break

    if log_path is None:
        log_path = _find_latest_scf_log()

    if log_path is None or not log_path.exists():
        console.print("[yellow]No SCF log file found. Looking in OUT.ABACUS/...[/yellow]")
        # Try to be helpful
        out_dirs = list(Path().glob("OUT.*"))
        if out_dirs:
            for d in out_dirs:
                console.print(f"  Found: {d}/")
                for f in sorted(d.glob("running*")):
                    console.print(f"    - {f.name}")
        else:
            console.print("[red]No OUT.* directories found. Run an ABACUS SCF calculation first.[/red]")
        return

    console.print(f"  [dim]Reading: {log_path}[/dim]")

    data = parse_scf_log(log_path)

    if not data["energies"]:
        console.print("[red]No energy data found in log file.[/red]")
        return

    nsteps = data["nsteps"]
    final_energy = data["final_energy"]
    converged = data["converged"]
    is_ev = data.get("final_energy_is_ev", False)

    from abacuscopilot.core.constants import RY_TO_EV

    if is_ev:
        final_energy_ev = final_energy
        final_energy_ry = final_energy / RY_TO_EV
    else:
        final_energy_ry = final_energy
        final_energy_ev = final_energy * RY_TO_EV
    final_energy_ha = final_energy_ry / 2.0

    console.print()
    console.print(f"  [bold]SCF steps:[/bold] {nsteps}")
    console.print(f"  [bold]Final energy:[/bold] {final_energy_ry:.8f} Ry  =  {final_energy_ha:.8f} Ha  =  {final_energy_ev:.6f} eV")
    console.print(f"  [bold]Converged:[/bold] {'[green]Yes[/green]' if converged else '[red]No[/red]'}")

    if data["ediffs"]:
        if is_ev:
            # ediffs are in eV — convert to Ry for display
            de_ev = abs(data["ediffs"][-1])
            de_ry = de_ev / RY_TO_EV
        else:
            de_ry = abs(data["ediffs"][-1])
            de_ev = de_ry * RY_TO_EV
        console.print(f"  [bold]Final |dE|:[/bold] {de_ev:.6e} eV")
        if len(data["ediffs"]) > 1:
            max_de_ev = max(abs(d) for d in data["ediffs"])
            if is_ev:
                max_de_ev = max_de_ev
            else:
                max_de_ev = max_de_ev * RY_TO_EV
            console.print(f"  [bold]Max |dE|:[/bold] {max_de_ev:.6e} eV")

    if data["drho"]:
        console.print(f"  [bold]Final |drho|:[/bold] {data['drho'][-1]:.6e}")

    # Optionally plot (only for multi-step SCF, not NSCF)
    if nsteps > 1 and interactive:
        from rich.prompt import Prompt
        do_plot = Prompt.ask(
            "\n  Plot convergence?", choices=["y", "n"], default="y"
        )
        if do_plot == "y":
            _plot_convergence(data, console)

    console.print()


def _plot_convergence(data: dict, console) -> None:
    """Plot SCF convergence history."""
    import matplotlib.pyplot as plt

    from abacuscopilot.core.constants import RY_TO_EV
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    is_ev = data.get("final_energy_is_ev", False)
    factor = 1.0 / RY_TO_EV if is_ev else 1.0  # eV → Ry

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Energy vs step
    energies = np.array(data["energies"]) * factor
    steps = np.arange(1, len(energies) + 1)
    ax1.plot(steps, energies, "o-", color="#1f77b4", markersize=4)
    ax1.set_xlabel("SCF Step")
    ax1.set_ylabel("Total Energy (Ry)")
    ax1.set_title("Energy Convergence")

    # |dE| vs step (semilogy — skip zeros which can't be displayed on log scale)
    if data["ediffs"]:
        ediffs = np.array(data["ediffs"]) * factor
        ediff_steps = np.arange(2, len(ediffs) + 2)
        # Filter zeros: log(0) = -∞, can't plot
        mask = ediffs > 0
        if mask.any():
            ax2.semilogy(ediff_steps[mask], ediffs[mask], "s-", color="#d62728", markersize=4)
        ax2.set_xlabel("SCF Step")
        ax2.set_ylabel("|ΔE| (Ry)")
        ax2.set_title("Energy Difference")
        ax2.axhline(y=1e-7, color="gray", linestyle="--", alpha=0.5, label="1e-7 Ry")
        ax2.legend()

    for ax in (ax1, ax2):
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        ax.tick_params(axis="both", direction="out")

    fig.tight_layout(pad=1.2)
    save_path = "scf_convergence.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    console.print(f"  [dim]Plot saved to {save_path}[/dim]")
    plt.close(fig)


# =============================================================================
# Task 712: Compare SCF convergence
# =============================================================================

@task(712, category="SCF Analysis", name="SCF Compare",
      description="Compare SCF convergence between multiple calculations")
def task_scf_compare(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compare SCF convergence across multiple runs."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Compare SCF Convergence ===[/bold cyan]")
    console.print()

    log_files = []
    if args:
        for arg in args:
            p = Path(arg)
            if p.is_dir():
                # Search directory for OUT.*/running* logs
                found = list(p.glob("OUT.*/running*"))
                if found:
                    log_files.extend(found)
                else:
                    console.print(f"  [yellow]No log files in {arg}[/yellow]")
            elif p.is_file():
                log_files.append(p)
    else:
        # Auto-find: OUT.*/ + any_subdir/OUT.*/ running* logs
        seen: set[str] = set()
        for pattern in ("OUT.*/running*", "*/OUT.*/running*"):
            for f in sorted(Path().glob(pattern)):
                if str(f) not in seen:
                    log_files.append(f)
                    seen.add(str(f))

    if interactive:
        if log_files:
            console.print(f"  Auto-found {len(log_files)} log(s):")
            for f in log_files:
                console.print(f"    - {f}")
            choice = console.input("  Use these? (y/n) [y]: ").strip().lower()
            if choice and choice not in ("y", "yes"):
                log_files = []
            if choice.lower() not in ("y", "yes"):
                log_files = []
        if not log_files:
            console.print("[dim]Enter paths (space-separated, blank to finish, q to quit):[/dim]")
            while True:
                raw = console.input("  > ").strip()
                if not raw:
                    break
                if raw.lower() in ("q", "quit"):
                    break
                for item in raw.split():
                    p = Path(item)
                    if p.is_dir():
                        found = list(p.glob("OUT.*/running*"))
                        if found:
                            log_files.extend(found)
                            console.print(f"  [dim]Found {len(found)} log(s) in {item}[/dim]")
                        else:
                            console.print(f"  [red]No log files in {item}[/red]")
                    elif p.is_file():
                        log_files.append(p)
                    else:
                        console.print(f"  [red]Not found: {item}[/red]")

    if len(log_files) < 2:
        console.print("[yellow]Need at least 2 log files to compare.[/yellow]")
        return

    console.print(f"  Comparing {len(log_files)} runs:")
    for f in log_files:
        console.print(f"    - {f}")

    # Parse all
    all_data = {}
    for f in log_files:
        all_data[str(f)] = parse_scf_log(f)

    # Print comparison table
    from rich.table import Table

    table = Table(title="SCF Convergence Comparison")
    table.add_column("Run", style="cyan", min_width=18, no_wrap=True)
    table.add_column("Steps", justify="right", min_width=7)
    table.add_column("Final E (Ry)", justify="right", min_width=16, no_wrap=True)
    table.add_column("Final E (eV)", justify="right", min_width=16, no_wrap=True)
    table.add_column("Conv.", justify="center", min_width=5)
    table.add_column("Final |dE|", justify="right", min_width=8, no_wrap=True)

    for label, data in all_data.items():
        p = Path(label)
        # Show meaningful dir name: ecutwfc_40, kspacing_0.140, etc.
        if p.parent.name.startswith("OUT.") and p.parent.parent.name:
            short_label = p.parent.parent.name
        else:
            short_label = p.parent.name if "OUT" in label else p.name
        conv = "✓" if data["converged"] else "✗"
        final_de = f"{data['ediffs'][-1]:.2e}" if data["ediffs"] else "N/A"
        from abacuscopilot.core.constants import RY_TO_EV
        e_ry = data.get("final_energy", 0.0)
        e_ev = e_ry * RY_TO_EV
        table.add_row(
            short_label,
            str(data["nsteps"]),
            f"{e_ry:.6f}" if e_ry else "N/A",
            f"{e_ev:.4f}" if e_ry else "N/A",
            conv,
            final_de,
        )

    console.print()
    console.print(table)
    console.print()


# =============================================================================
# Task 713 (also accessible as task 9908): Calculation Status Check
# =============================================================================


def _find_abacus_output_dir() -> Path | None:
    """Find the ABACUS output directory.

    Prefers OUT.* subdirectories (e.g., OUT.ABACUS).  Falls back to the
    current directory only when it contains running log files directly.
    """
    cwd = Path(".")
    # OUT.* subdirectory takes priority — NSCF calcs may have
    # ABACUS-CHARGE-DENSITY.restart in cwd, which is NOT an output dir.
    out_dirs = sorted(cwd.glob("OUT.*"))
    if out_dirs:
        return max(out_dirs, key=lambda p: p.stat().st_mtime)
    # Fallback: cwd *is* the output dir (user cd'd into OUT.ABACUS)
    if list(cwd.glob("running_*.log")):
        return cwd
    return None


def _parse_calculation_status(out_dir: Path) -> dict:
    """Parse an ABACUS output directory for calculation status.

    Returns a dict with completion status, convergence, energies, forces, stress.
    All energies are in Ry (ABACUS native).
    """
    result = {
        "completed": False,
        "converged": False,
        "final_energy_ry": None,
        "n_scf_steps": 0,
        "n_ionic_steps": 0,
        "max_force": None,
        "max_stress": None,
        "wall_time": None,
        "natom": 0,
        "errors": [],
        "output_dir": str(out_dir),
    }

    # Find the main log file
    log_files = sorted(out_dir.glob("running_*.log")) + sorted(out_dir.glob("output*")) + sorted(out_dir.glob("*.log"))
    if not log_files:
        result["errors"].append("No log file found")
        return result

    log_path = log_files[0]
    try:
        content = log_path.read_text(errors="ignore")
    except Exception:
        result["errors"].append(f"Cannot read {log_path}")
        return result

    # Check for normal completion
    completed_re = re.compile(
        r"(?:REACH|JOB\s+DONE|calculation\s+finished|PROGRAM\s+ENDS|END\s+OF\s+ABACUS"
        r"|FINISH\s+Time|Finish\s+Time)",
        re.IGNORECASE,
    )
    result["completed"] = bool(completed_re.search(content))

    # SCF convergence. Markers differ between ABACUS builds:
    #   CPU build: "charge density convergence is achieved"
    #   GPU build: "#SCF IS CONVERGED#"
    # (The old "converge[nc]ed" was a typo matching "convergned"/"convergced",
    #  never the real "converged", so GPU logs read as not converged.)
    converged_re = re.compile(
        r"(?:converged|convergence\s+is\s+achieved|SCF\s+is\s+converged"
        r"|SCF\s+done|reach\s+convergence)",
        re.IGNORECASE,
    )
    result["converged"] = bool(converged_re.search(content))

    # Final total energy — v3.10 outputs eV; older versions use Ry
    # Patterns in priority order (most specific first)
    energy_patterns = [
        r"final\s+etot\s+is\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*eV",  # v3.10 PW
        r"!FINAL_ETOT_IS\s*[=:]*\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*eV",  # v3.x LCAO
        r"FINAL\s+ENERGY\s*[=:]\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
        r"(?:total energy|ETOT|ENERGY)\s*[=:]\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
    ]
    for pat in energy_patterns:
        matches = re.findall(pat, content, re.IGNORECASE)
        if matches:
            e_val = float(matches[-1])
            # Patterns with "eV" suffix are already in eV; convert to Ry for storage
            if "eV" in pat:
                from abacuscopilot.core.constants import RY_TO_EV
                result["final_energy_ry"] = e_val / RY_TO_EV
            else:
                result["final_energy_ry"] = e_val
            break

    # SCF steps — ABACUS v3.10 counts by E_KohnSham lines
    scf_re = re.compile(r"(?:SCF|ITER)\s*=\s*(\d+)", re.IGNORECASE)
    scf_matches = scf_re.findall(content)
    if scf_matches:
        result["n_scf_steps"] = int(scf_matches[-1])
    else:
        # v3.10+: count E_KohnSham lines (one per SCF step)
        n_kohnsham = len(re.findall(r"^\s*E_KohnSham", content, re.MULTILINE))
        if n_kohnsham > 0:
            result["n_scf_steps"] = n_kohnsham

    # Wall time
    time_re = re.compile(r"Total\s+Time\s*:\s*(.+)", re.IGNORECASE)
    time_match = time_re.search(content)
    if time_match:
        result["wall_time"] = time_match.group(1).strip()
    else:
        result["wall_time"] = None

    # Ionic steps (for relax/MD)
    ionic_re = re.compile(r"(?:ION|MD)\s*STEP\s*[=:]*\s*(\d+)", re.IGNORECASE)
    ionic_matches = ionic_re.findall(content)
    if ionic_matches:
        result["n_ionic_steps"] = int(ionic_matches[-1])

    # Max force (eV/Å)
    force_re = re.compile(
        r"(?:TOTAL-FORCE|max\s*force|MAX\s*FORCE)\s*[=:]*\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    force_matches = force_re.findall(content)
    if force_matches:
        result["max_force"] = float(force_matches[-1])

    # Max stress (kbar)
    stress_re = re.compile(
        r"(?:TOTAL-STRESS|max\s*stress|MAX\s*STRESS|PRESSURE)\s*[=:]*\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    stress_matches = stress_re.findall(content)
    if stress_matches:
        result["max_stress"] = float(stress_matches[-1])

    # Atom count from ABACUS log (authoritative)
    natom_re = re.compile(r"TOTAL\s+ATOM\s+NUMBER\s*[=:]\s*(\d+)", re.IGNORECASE)
    natom_match = natom_re.search(content)
    if natom_match:
        result["natom"] = int(natom_match.group(1))

    return result


@task(713, category="SCF Analysis", name="Calc Status",
      description="Check ABACUS calculation completion, convergence, and final energy")
@task(9908, category="System", name="Calc Status",
      description="Quick check: is the calculation done? Converged? Final energy?")
def task_calc_status(args: list[str] | None = None, interactive: bool = True) -> None:
    """Check the status of an ABACUS calculation in the current directory.

    Reports:
    - Whether the calculation completed normally
    - SCF convergence status
    - Final total energy (Ry, Ha, eV)
    - Forces and stress (if available)
    """
    console = _get_console()
    from abacuscopilot.core.constants import RY_TO_EV

    console.print()
    console.print("[bold cyan]=== Calculation Status ===[/bold cyan]")
    console.print()

    # Find output directory
    out_dir = None
    if args:
        for arg in args:
            p = Path(arg)
            if p.is_dir():
                out_dir = p
                break
            if p.exists():
                # A specific log file
                console.print(f"  [dim]Reading: {p}[/dim]")
                out_dir = p.parent
                # Quick parse
                status = _parse_calculation_status(p.parent)
                break

    if out_dir is None:
        out_dir = _find_abacus_output_dir()

    if out_dir is None:
        console.print("[yellow]No OUT.* directory found.[/yellow]")
        console.print("[dim]Run an ABACUS calculation first. Expected directory: OUT.ABACUS/[/dim]")
        return

    status = _parse_calculation_status(out_dir)
    console.print(f"  [dim]Output: {status['output_dir']}[/dim]")
    console.print()

    # --- Completion ---
    icon = "[green]✓[/green]" if status["completed"] else "[yellow]…[/yellow]"
    console.print(f"  {icon} Completed: {'[green]Yes[/green]' if status['completed'] else '[yellow]Not yet / still running[/yellow]'}")

    # --- Convergence ---
    icon = "[green]✓[/green]" if status["converged"] else "[red]✗[/red]" if status["completed"] else "[dim]—[/dim]"
    console.print(f"  {icon} SCF Converged: {'[green]Yes[/green]' if status['converged'] else '[red]No[/red]'}")

    # --- Energy ---
    if status["final_energy_ry"] is not None:
        e_ry = status["final_energy_ry"]
        e_ha = e_ry / 2.0        # 1 Ha = 2 Ry
        e_ev = e_ry * RY_TO_EV
        console.print()
        console.print("  [bold]Final Energy:[/bold]")
        console.print(f"    {e_ry:.8f} Ry   (ABACUS native)")
        console.print(f"    {e_ha:.8f} Ha   (Hartree)")
        console.print(f"    {e_ev:.6f} eV")
    else:
        console.print()
        console.print("  [yellow]Final energy not found in log.[/yellow]")

    # --- Steps ---
    if status["n_scf_steps"]:
        console.print(f"\n  [bold]SCF Steps:[/bold] {status['n_scf_steps']}")
    if status["n_ionic_steps"]:
        console.print(f"  [bold]Ionic Steps:[/bold] {status['n_ionic_steps']}")

    # --- Forces ---
    if status["max_force"] is not None:
        console.print(f"\n  [bold]Max Force:[/bold] {status['max_force']:.6f} eV/Å")

    # --- Stress ---
    if status["max_stress"] is not None:
        from abacuscopilot.core.constants import KBAR_TO_GPA
        stress_gpa = status["max_stress"] * KBAR_TO_GPA
        console.print(f"  [bold]Max Stress:[/bold] {status['max_stress']:.4f} kbar = {stress_gpa:.4f} GPa")

    # --- Errors ---
    if status["errors"]:
        console.print()
        for e in status["errors"]:
            console.print(f"  [red]![/red] {e}")

    console.print()
