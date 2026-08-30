"""SCF convergence analysis tasks for ABACUS output.

Task IDs 701-709

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
# Task 701: SCF convergence check
# =============================================================================

@task(701, category="SCF Analysis", name="SCF Convergence",
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
# Task 702: Compare SCF convergence
# =============================================================================

def _extract_wall_time(log_path: Path) -> float | None:
    """Extract wall time in seconds from an ABACUS log file.

    Looks for a ``Total Time`` line first (the authoritative total).
    Falls back to summing per-ion-step ``DONE : INIT SCF Time`` lines.
    Handles three ABACUS formats:

    * ``TOTAL  Time  : 15``  (integer seconds, often in ``*.out``)
    * ``Total  Time  : 0 h 0 mins 15 secs``  (human-readable)
    * ``DONE : INIT SCF Time : 1.23 (SEC)``  (per-ion-step)
    """
    _re_total_int = re.compile(r"TOTAL\s+Time\s*:\s*(\d+)", re.IGNORECASE)
    _re_total_hms = re.compile(
        r"Total\s+Time\s*:\s*(?:(\d+)\s*h\s*)?\s*(?:(\d+)\s*mins?\s*)?\s*(\d+)\s*secs?",
        re.IGNORECASE,
    )
    _re_per_step = re.compile(
        r"DONE\s*:\s*INIT\s+SCF\s+Time\s*:\s*(\d+\.?\d*)", re.IGNORECASE,
    )

    try:
        text = log_path.read_text(errors="ignore")
    except OSError:
        return None

    # 1) Try "Total  Time  : 0 h 0 mins 15 secs" (human-readable)
    m = _re_total_hms.search(text)
    if m:
        h = int(m.group(1) or 0)
        mi = int(m.group(2) or 0)
        s = int(m.group(3) or 0)
        return float(h * 3600 + mi * 60 + s)

    # 2) Try "TOTAL  Time  : 15" (simple integer seconds)
    m = _re_total_int.search(text)
    if m:
        return float(m.group(1))

    # 3) Sum per-ion-step SCF times
    step_times = [float(t) for t in _re_per_step.findall(text)]
    if step_times:
        return round(sum(step_times), 1)

    return None


def _format_wall_time(seconds: float) -> str:
    """Format seconds as human-readable wall time."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        h = int(seconds // 3600)
        m = int((seconds % 3600) // 60)
        return f"{h}h{m:02d}m"


@task(702, category="SCF Analysis", name="SCF Compare",
      description="Compare SCF convergence between multiple calculations")
def task_scf_compare(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compare SCF convergence across multiple runs.

    Auto-discovers subdirectories under the current working directory
    (``*/OUT.*/running*``), or accepts explicit file/directory paths.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Compare SCF Convergence ===[/bold cyan]")
    console.print()

    log_files: list[Path] = []
    if args:
        for arg in args:
            p = Path(arg)
            if p.is_dir():
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
    all_data: dict[str, dict] = {}
    for f in log_files:
        data = parse_scf_log(f)
        data["wall_time_s"] = _extract_wall_time(f)
        all_data[str(f)] = data

    # ---- Comparison table ----
    from rich.table import Table

    table = Table(title="SCF Convergence Comparison")
    table.add_column("Run", style="cyan", min_width=18, no_wrap=True)
    table.add_column("Steps", justify="right", min_width=7)
    table.add_column("Final E (Ry)", justify="right", min_width=16, no_wrap=True)
    table.add_column("Final E (eV)", justify="right", min_width=16, no_wrap=True)
    table.add_column("Conv.", justify="center", min_width=5)
    table.add_column("Final |dE|", justify="right", min_width=8, no_wrap=True)
    table.add_column("Wall Time", justify="right", min_width=8)

    from abacuscopilot.core.constants import RY_TO_EV

    for label, data in all_data.items():
        p = Path(label)
        if p.parent.name.startswith("OUT.") and p.parent.parent.name:
            short_label = p.parent.parent.name
        else:
            short_label = p.parent.name if "OUT" in label else p.name

        conv = "[green]✓[/green]" if data["converged"] else "[red]✗[/red]"
        final_de = f"{data['ediffs'][-1]:.2e}" if data["ediffs"] else "N/A"
        e_ry = data.get("final_energy", 0.0) or 0.0
        e_ev = e_ry * RY_TO_EV

        wt = data.get("wall_time_s")
        wall_str = _format_wall_time(wt) if wt is not None else "[dim]—[/dim]"

        table.add_row(
            short_label,
            str(data["nsteps"]),
            f"{e_ry:.6f}" if e_ry else "N/A",
            f"{e_ev:.4f}" if e_ry else "N/A",
            conv,
            final_de,
            wall_str,
        )

    console.print()
    console.print(table)
    console.print()


# =============================================================================
# Task 703: Per-ionic-step summary table for relax, MD, and SCF runs
# =============================================================================


# ---------------------------------------------------------------------------
# Task 703 helpers — per-ionic-step summary
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
        r"!FINAL_ETOT_IS\s*[=:]*\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*eV",  # v3.x LCAO (older)
        r"#TOTAL\s+ENERGY#\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*eV",  # v3.x LCAO (newer)
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


# ---------------------------------------------------------------------------

# Patterns for ABACUS log parsing (compiled once at import time)
_RE_STEP_RELAX = re.compile(r"STEP\s+OF\s+RELAXATION\s*:\s*(\d+)", re.IGNORECASE)
_RE_STEP_MD = re.compile(r"STEP\s+OF\s+MOLECULAR\s+DYNAMICS\s*:\s*(\d+)", re.IGNORECASE)
_RE_ION_ELEC_LEGACY = re.compile(r"ION=\s*(\d+)\s+ELEC=\s*(\d+)")
_RE_ION_ELEC_V3 = re.compile(r"#ION\s+MOVE#\s+(\d+)\s+#ELEC\s+ITER#\s+(\d+)")
_RE_FINAL_ETOT = re.compile(r"final\s+etot\s+is\s+(-?\d+\.?\d+(?:[eE][+-]?\d+)?)\s*eV", re.IGNORECASE)
_RE_FINAL_ETOT_BANG = re.compile(r"!FINAL_ETOT_IS\s*(-?\d+\.?\d+(?:[eE][+-]?\d+)?)\s*eV", re.IGNORECASE)
_RE_TOTAL_ENERGY = re.compile(r"#TOTAL\s+ENERGY#\s*(-?\d+\.?\d+(?:[eE][+-]?\d+)?)\s*eV", re.IGNORECASE)
_RE_CONVERGED = re.compile(
    r"(?:charge\s+density\s+convergence\s+is\s+achieved|#SCF\s+IS\s+CONVERGED#"
    r"|converged|SCF\s+done|reach\s+convergence)",
    re.IGNORECASE,
)
# Ionic (geometry) convergence — per-step relaxation verdict
_RE_IONIC_CONVERGED = re.compile(
    r"Relaxation\s+is\s+converged", re.IGNORECASE,
)
_RE_IONIC_NOT_CONVERGED = re.compile(
    r"Relaxation\s+is\s+not\s+converged", re.IGNORECASE,
)
_RE_FORCE_HEADER = re.compile(r"TOTAL-FORCE\s+\(eV/Angstrom\)", re.IGNORECASE)
_RE_FORCE_LINE = re.compile(
    r"^\s*\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
)
_RE_WALLTIME = re.compile(r"Total\s+Time\s*:\s*(.+)", re.IGNORECASE)
_RE_COMPLETED = re.compile(
    r"(?:REACH|JOB\s+DONE|calculation\s+finished|PROGRAM\s+ENDS|END\s+OF\s+ABACUS"
    r"|FINISH\s+Time|Finish\s+Time)",
    re.IGNORECASE,
)
_RE_NATOM = re.compile(r"TOTAL\s+ATOM\s+NUMBER\s*[=:]\s*(\d+)", re.IGNORECASE)


def _parse_ionic_steps(out_dir: Path) -> dict:
    """Parse ABACUS log line-by-line extracting per-ionic-step data.

    Returns a dict with ``calc_type``, ``completed``, ``wall_time``,
    ``natom``, and a ``steps`` list of per-ion-step dicts::

        {"ion": int, "elec": int, "converged": bool,
         "energy_ev": float|None, "max_force": float|None}

    Handles relax, cell-relax, MD, and SCF-only runs.  Works with both
    legacy (``ION=N ELEC=M``) and v3.x (``#ION MOVE#``) log formats.
    """
    result: dict = {
        "calc_type": "scf",
        "completed": False,
        "wall_time": None,
        "natom": 0,
        "steps": [],
    }

    # Find the main log file
    log_files = (
        sorted(out_dir.glob("running_*.log"))
        or sorted(out_dir.glob("output*"))
        or sorted(out_dir.glob("*.log"))
    )
    if not log_files:
        return result

    log_path = log_files[0]

    # State-machine variables
    cur_step: dict | None = None  # current ionic step being built
    calc_detected = False
    in_force_block = False
    _force_sep_seen = False  # tracks the "---" separator right after TOTAL-FORCE header

    def _push_step(step: dict):
        """Finalise *step* and append it to the result list."""
        if step is None:
            return
        # Only keep steps that have at least one electronic iteration
        if step.get("elec", 0) > 0:
            result["steps"].append(step)

    try:
        with open(log_path, errors="ignore") as fh:
            for line in fh:
                # --- Detect calculation type from step boundary markers ---
                relax_match = _RE_STEP_RELAX.search(line)
                md_match = _RE_STEP_MD.search(line)
                if relax_match:
                    calc_detected = True
                    if result["calc_type"] != "cell-relax":
                        result["calc_type"] = "relax"
                    _push_step(cur_step)
                    cur_step = {
                        "ion": int(relax_match.group(1)),
                        "elec": 0,
                        "converged": False,
                        "ionic_converged": None,
                        "energy_ev": None,
                        "max_force": None,
                    }
                    in_force_block = False
                    continue
                if md_match:
                    calc_detected = True
                    result["calc_type"] = "md"
                    _push_step(cur_step)
                    cur_step = {
                        "ion": int(md_match.group(1)),
                        "elec": 0,
                        "converged": False,
                        "ionic_converged": None,
                        "energy_ev": None,
                        "max_force": None,
                    }
                    in_force_block = False
                    continue

                # --- Electronic step counters ---
                ion_elec_legacy = _RE_ION_ELEC_LEGACY.search(line)
                ion_elec_v3 = _RE_ION_ELEC_V3.search(line)
                elec_match = ion_elec_legacy or ion_elec_v3
                if elec_match:
                    ion_n = int(elec_match.group(1))
                    elec_n = int(elec_match.group(2))
                    # Detect ionic step boundary: #ION MOVE# increment
                    if cur_step is not None and cur_step.get("elec", 0) > 0 and ion_n > cur_step.get("ion", 0):
                        _push_step(cur_step)
                        cur_step = {
                            "ion": ion_n,
                            "elec": 0,
                            "converged": False,
                            "ionic_converged": None,
                            "energy_ev": None,
                            "max_force": None,
                        }
                    # If we haven't seen a STEP marker yet (SCF-only),
                    # start an implicit first ionic step.
                    if cur_step is None:
                        cur_step = {
                            "ion": ion_n,
                            "elec": 0,
                            "converged": False,
                            "ionic_converged": None,
                            "energy_ev": None,
                            "max_force": None,
                        }
                        calc_detected = True
                    # Ensure the current step tracks the right ion number
                    # (handles rare cases where step boundary markers differ)
                    if cur_step is not None:
                        cur_step["ion"] = max(cur_step["ion"], ion_n)
                        cur_step["elec"] = max(cur_step["elec"], elec_n)
                    in_force_block = False
                    continue

                # --- Ionic (geometry) convergence ---
                # Must be checked BEFORE SCF convergence because
                # "Relaxation is not converged yet!" contains the word
                # "converged" and would be falsely caught by _RE_CONVERGED.
                if cur_step is not None and _RE_IONIC_CONVERGED.search(line):
                    cur_step["ionic_converged"] = True
                    continue
                if cur_step is not None and _RE_IONIC_NOT_CONVERGED.search(line):
                    cur_step["ionic_converged"] = False
                    continue

                # --- SCF convergence ---
                if cur_step is not None and _RE_CONVERGED.search(line):
                    cur_step["converged"] = True
                    continue

                # --- Energy ---
                etot_match = _RE_FINAL_ETOT.search(line) or _RE_FINAL_ETOT_BANG.search(line) or _RE_TOTAL_ENERGY.search(line)
                if etot_match and cur_step is not None:
                    try:
                        cur_step["energy_ev"] = float(etot_match.group(1))
                    except ValueError:
                        pass
                    continue

                # --- TOTAL-FORCE block ---
                # Two formats exist:
                #   PW / MD:  header → --- → force_data  (1 separator)
                #   LCAO GPU: header → --- → column_hdr → --- → force_data  (2 separators)
                if _RE_FORCE_HEADER.search(line):
                    in_force_block = True
                    _force_sep_count = 0
                    continue
                if in_force_block and cur_step is not None:
                    stripped = line.strip()
                    # Separator line
                    if stripped.startswith("---"):
                        _force_sep_count += 1
                        continue
                    if _force_sep_count == 0:
                        continue  # before first separator
                    # After first separator: try parsing as force data (PW format)
                    fmatch = _RE_FORCE_LINE.match(line)
                    if fmatch:
                        fx = abs(float(fmatch.group(1)))
                        fy = abs(float(fmatch.group(2)))
                        fz = abs(float(fmatch.group(3)))
                        fmax = max(fx, fy, fz)
                        if cur_step["max_force"] is None or fmax > cur_step["max_force"]:
                            cur_step["max_force"] = fmax
                        continue
                    # Not a force line after first separator:
                    # LCAO GPU has a column header here ("Atoms  Force_x  …") —
                    # wait for the second separator, then data.
                    if _force_sep_count >= 2:
                        # Force data started but this line isn't force → exit block
                        in_force_block = False
                    # else: _force_sep_count == 1, non-force line → skip (header row)

                # --- Wall time ---
                twall = _RE_WALLTIME.search(line)
                if twall:
                    result["wall_time"] = twall.group(1).strip()
                    continue

                # --- Completion ---
                if not result["completed"] and _RE_COMPLETED.search(line):
                    result["completed"] = True
                    continue

                # --- Atom count ---
                if result["natom"] == 0:
                    natom_match = _RE_NATOM.search(line)
                    if natom_match:
                        result["natom"] = int(natom_match.group(1))

    except OSError:
        pass

    # Push the last step (possibly incomplete / still running)
    _push_step(cur_step)

    # Refine calc_type if the INPUT file is available (cell-relax vs relax)
    try:
        input_path = out_dir.parent / "INPUT"
        if not input_path.exists():
            input_path = out_dir / "INPUT"
        if input_path.exists():
            inp = input_path.read_text(errors="ignore")
            for param_line in inp.split("\n"):
                stripped = param_line.strip()
                if re.match(r"^\s*calculation\s+", stripped, re.IGNORECASE):
                    calc_val = stripped.split(None, 1)[-1].strip().lower()
                    if calc_val in ("cell-relax", "relax", "md", "scf", "nscf"):
                        if calc_val == "cell-relax":
                            result["calc_type"] = "cell-relax"
                        elif calc_val == "relax":
                            result["calc_type"] = "relax"
                        elif calc_val in ("md", "nscf"):
                            result["calc_type"] = calc_val
                    break
    except OSError:
        pass

    # If no STEP markers but we have electronic steps, it's SCF
    if not calc_detected and not result["steps"]:
        result["calc_type"] = "scf"

    return result


@task(703, category="SCF Analysis", name="Ion Steps",
      description="Per-ion-step summary: SCF iters, convergence, energy, forces")
def task_sys_status(args: list[str] | None = None, interactive: bool = True) -> None:
    """Show a per-ionic-step summary table — ideal for monitoring relax/MD runs."""
    console = _get_console()
    from rich.table import Table

    console.print()
    console.print("[bold cyan]=== Calculation Status ===[/bold cyan]")
    console.print()

    out_dir = _find_abacus_output_dir()
    if out_dir is None:
        console.print("[yellow]No OUT.* directory found.[/yellow]")
        return

    data = _parse_ionic_steps(out_dir)
    steps = data["steps"]

    if not steps:
        console.print("[yellow]No ionic/electronic step data found in log.[/yellow]")
        return

    # ---- Header line ----
    calc_labels = {
        "scf": "SCF", "relax": "relax", "cell-relax": "cell-relax",
        "md": "MD", "nscf": "NSCF",
    }
    calc_name = calc_labels.get(data["calc_type"], data["calc_type"])
    completed = data["completed"]
    status_icon = "[green]✓[/green]" if completed else "[yellow]…[/yellow]"
    status_text = "completed" if completed else "running"
    is_relax = data["calc_type"] in ("relax", "cell-relax")
    is_md = data["calc_type"] == "md"

    # Convergence counts depend on calculation type
    if is_relax:
        n_scf_conv = sum(1 for s in steps if s["converged"])
        n_ionic_conv = sum(1 for s in steps if s["ionic_converged"] is True)
        n_total = len(steps)
    else:
        n_scf_conv = sum(1 for s in steps if s["converged"])
        n_total = len(steps)

    header_parts = [
        f"[dim]Output:[/dim] {out_dir}",
        f"[dim]Calc:[/dim] {calc_name}",
        f"[dim]Status:[/dim] {status_icon} {status_text}",
    ]
    if data["wall_time"]:
        header_parts.append(f"[dim]Wall:[/dim] {data['wall_time']}")
    console.print(f"  {'  |  '.join(header_parts)}")
    console.print()

    # ---- Build the table ----
    show_force = is_relax

    # Choose Conv column header based on calc type
    conv_header = "Conv" if is_relax else "SCF"

    table = Table(expand=False, box=None, padding=(0, 1))
    table.add_column("ION", justify="right", style="dim", min_width=4)
    table.add_column("SCF", justify="right", min_width=4)
    table.add_column(conv_header, justify="center", min_width=5)
    table.add_column("Energy(eV)", justify="right", min_width=16)
    table.add_column("ΔE(eV)", justify="right", min_width=10)
    if show_force:
        table.add_column("Max|F|", justify="right", min_width=9)

    prev_energy = None
    for s in steps:
        ion = str(s["ion"])
        elec = str(s["elec"])

        # Convergence cell: for relax, show ionic convergence;
        # fall back to SCF convergence when no ionic verdict (e.g. running step)
        if is_relax:
            if s["ionic_converged"] is True:
                conv = "[green]✓[/green]"
            elif s["ionic_converged"] is False:
                # SCF converged but ionic not yet — normal relax progression
                if s["converged"]:
                    conv = "[yellow]→[/yellow]"
                else:
                    conv = "[red]⚡[/red]"  # SCF failed too
            else:
                # No ionic verdict yet (last step still running)
                conv = "[green]✓[/green]" if s["converged"] else "[red]⚡[/red]"
        else:
            # MD / SCF: show SCF convergence
            conv = "[green]✓[/green]" if s["converged"] else "[red]✗[/red]"

        if s["energy_ev"] is not None:
            energy = f"{s['energy_ev']: .4f}"
            if prev_energy is not None:
                delta = s["energy_ev"] - prev_energy
                delta_str = f"{delta:+.4f}"
            else:
                delta_str = "[dim]—[/dim]"
            prev_energy = s["energy_ev"]
        else:
            energy = "[yellow](no energy)[/yellow]"
            delta_str = "[dim]—[/dim]"

        row = [ion, elec, conv, energy, delta_str]
        if show_force:
            if s["max_force"] is not None:
                row.append(f"{s['max_force']:.4f}")
            else:
                row.append("[dim]—[/dim]")
        table.add_row(*row)

    console.print(table)

    # ---- Footer ----
    if is_relax:
        # Ionic convergence — just say yes or no
        if n_ionic_conv > 0:
            footer_parts = ["[green]✓ Ionic converged[/green]"]
        else:
            footer_parts = ["[yellow]→ Ionic not converged yet[/yellow]"]
        # SCF is secondary info — only mention if any step failed
        if n_scf_conv < n_total:
            failed_indices = [str(s["ion"]) for s in steps if not s["converged"]]
            footer_parts.append(
                f"[dim]SCF: {n_scf_conv}/{n_total} converged"
                + (f" (step{'s' if len(failed_indices) > 1 else ''} {','.join(failed_indices)} failed)" if failed_indices else "")
                + "[/dim]"
            )
    else:
        scf_icon = "[green]✓[/green]" if n_scf_conv == n_total else "[red]✗[/red]"
        footer_parts = [
            f"{scf_icon} [bold]SCF: {n_scf_conv}/{n_total} converged[/bold]"
        ]
    if data["natom"] > 0:
        footer_parts.append(f"{data['natom']} atoms")
    if is_relax:
        if n_total > 0 and steps[-1]["energy_ev"] is not None:
            final_e = steps[-1]["energy_ev"]
            footer_parts.append(f"final {final_e:.4f} eV")
            if data["natom"] > 0:
                footer_parts.append(f"{final_e / data['natom']:.4f} eV/atom")
    if data["wall_time"]:
        footer_parts.append(f"wall: {data['wall_time']}")
    console.print(f"  {'  |  '.join(footer_parts)}")

    # Legend for Conv column symbols
    if is_relax:
        console.print()
        console.print(
            "  [dim]Conv:  [green]✓[/] converged  [yellow]→[/] in progress  [red]⚡[/] SCF failed[/dim]"
        )

    console.print()
