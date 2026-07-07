"""Mechanical properties analysis tasks for ABACUS output.

Task IDs 751-769

Computes elastic constants, bulk/shear/Young's moduli, Poisson ratio
from ABACUS elastic tensor output. Also supports equation-of-state (EOS)
fitting from energy-volume data.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import task

# =============================================================================
# Elastic tensor parser
# =============================================================================


def parse_elastic_tensor(filepath: str | Path) -> dict[str, Any] | None:
    """Parse elastic tensor from ABACUS output.

    ABACUS outputs the 6×6 elastic tensor (Voigt notation) in the
    running log or a dedicated elastic output file.

    Args:
        filepath: Path to the output file.

    Returns:
        Dict with 'C' (6x6 tensor in GPa), 'unit', and raw data,
        or None if not found.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    content = filepath.read_text()

    # Look for "ELASTIC CONSTANTS" section in ABACUS output
    # Format varies by version - try multiple patterns
    patterns = [
        # ABACUS v3.x: "ELASTIC CONSTANTS" header followed by 6 rows
        r"ELASTIC\s+CONSTANTS.*?\n(.*?)(?:\n\s*\n|$)",
        # "ELASTIC TENSOR (GPa)" format
        r"(?:ELASTIC|ELASTICITY).*?(?:TENSOR|CONSTANTS).*?\n(.*?)(?:\n\s*\n|$)",
    ]

    C = np.zeros((6, 6))
    found = False

    # Try parsing the elastic constants section
    for pattern in patterns:
        m = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
        if m:
            section = m.group(1)
            rows = []
            for line in section.strip().split("\n"):
                parts = line.split()
                nums = []
                for p in parts:
                    try:
                        nums.append(float(p))
                    except ValueError:
                        break
                if len(nums) >= 6:
                    rows.append(nums[:6])
            if len(rows) >= 6:
                C = np.array(rows[:6])
                found = True
                break

    # Try simpler: look for "C11 C12 C13 ..." format
    if not found:
        # Look for 6 lines each containing 6 numbers near "Stiffness" or "Elastic"
        lines = content.split("\n")
        for i, line in enumerate(lines):
            if any(kw in line.upper() for kw in ("STIFFNESS", "ELASTIC", "C11")):
                try:
                    row_data = []
                    for j in range(i, min(i + 6, len(lines))):
                        parts = [float(x) for x in lines[j].split() if _is_number(x)]
                        if len(parts) >= 6:
                            row_data.append(parts[:6])
                    if len(row_data) == 6:
                        C = np.array(row_data)
                        found = True
                        break
                except Exception:
                    continue

    if not found:
        return None

    return {
        "C": C,
        "unit": "GPa",
    }


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# =============================================================================
# Mechanical properties from elastic tensor
# =============================================================================


def compute_mechanical_properties(C: np.ndarray) -> dict[str, float]:
    """Compute mechanical properties from the 6×6 elastic tensor (Voigt, GPa).

    Uses the Voigt-Reuss-Hill (VRH) approximation for polycrystalline averages.

    Args:
        C: 6×6 elastic tensor in Voigt notation (GPa).

    Returns:
        Dict with Voigt/Reuss/Hill moduli, Poisson ratio, Pugh ratio.
    """
    C = np.array(C)
    if C.shape != (6, 6):
        raise ValueError(f"Expected 6×6 elastic tensor, got {C.shape}")

    # Voigt average (upper bound)
    c11_22_33 = (C[0, 0] + C[1, 1] + C[2, 2]) / 3
    c12_13_23 = (C[0, 1] + C[0, 2] + C[1, 2]) / 3
    c44_55_66 = (C[3, 3] + C[4, 4] + C[5, 5]) / 3

    B_V = (c11_22_33 + 2 * c12_13_23) / 3
    G_V = (c11_22_33 - c12_13_23 + 3 * c44_55_66) / 5

    # Reuss average (lower bound)
    S = np.linalg.inv(C)  # Compliance tensor
    s11_22_33 = (S[0, 0] + S[1, 1] + S[2, 2]) / 3
    s12_13_23 = (S[0, 1] + S[0, 2] + S[1, 2]) / 3
    s44_55_66 = (S[3, 3] + S[4, 4] + S[5, 5]) / 3

    B_R = 1 / (3 * s11_22_33 + 6 * s12_13_23)
    G_R = 15 / (12 * s11_22_33 - 12 * s12_13_23 + 3 * s44_55_66)

    # Hill average
    B_H = (B_V + B_R) / 2
    G_H = (G_V + G_R) / 2

    # Young's modulus and Poisson ratio (Hill)
    E_H = 9 * B_H * G_H / (3 * B_H + G_H)
    nu_H = (3 * B_H - 2 * G_H) / (2 * (3 * B_H + G_H))

    # Pugh ratio (ductility indicator: >1.75 = ductile, <1.75 = brittle)
    pugh = B_H / G_H

    # Anisotropy factors
    A_Zener = 2 * C[3, 3] / (C[0, 0] - C[0, 1]) if abs(C[0, 0] - C[0, 1]) > 1e-10 else 0.0

    return {
        "B_V": B_V, "B_R": B_R, "B_H": B_H,
        "G_V": G_V, "G_R": G_R, "G_H": G_H,
        "E_H": E_H,
        "nu_H": nu_H,
        "pugh_ratio": pugh,
        "A_Zener": A_Zener,
    }


# =============================================================================
# Equation of State (EOS) fitting
# =============================================================================


def fit_eos_birch_murnaghan(
    volumes: np.ndarray,
    energies: np.ndarray,
) -> dict[str, float]:
    """Fit a 3rd-order Birch-Murnaghan equation of state.

    E(V) = E0 + 9V0*B0/16 * {[(V0/V)^(2/3)-1]^3 * B0' + [(V0/V)^(2/3)-1]^2 * [6-4*(V0/V)^(2/3)]}

    Args:
        volumes: Volume array (any unit).
        energies: Energy array (any unit).

    Returns:
        Dict with 'E0', 'V0', 'B0', 'B0_prime'.
    """
    from scipy.optimize import curve_fit

    def bm3(V, E0, V0, B0, B0p):
        """3rd-order Birch-Murnaghan EOS."""
        eta = (V0 / V) ** (2.0 / 3.0)
        return E0 + 9.0 * V0 * B0 / 16.0 * (
            (eta - 1.0) ** 3 * B0p +
            (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
        )

    # Initial guesses
    imin = np.argmin(energies)
    E0_guess = energies[imin]
    V0_guess = volumes[imin]
    B0_guess = 100.0  # GPa equivalent
    B0p_guess = 4.0

    try:
        popt, _ = curve_fit(
            bm3, volumes, energies,
            p0=[E0_guess, V0_guess, B0_guess, B0p_guess],
            maxfev=10000,
        )
        return {"E0": popt[0], "V0": popt[1], "B0": popt[2], "B0_prime": popt[3]}
    except Exception:
        return {"E0": E0_guess, "V0": V0_guess, "B0": 0.0, "B0_prime": 0.0}


# =============================================================================
# Task 751: Parse elastic constants
# =============================================================================


@task(1201, category="Mechanics", name="Elastic Constants",
      description="Parse elastic tensor from ABACUS output and compute mechanical properties",
      cli_args=[
          {"name": "--file", "type": str, "default": None, "help": "Path to elastic tensor output"},
      ])
def task_elastic_constants(args: list[str] | None = None, interactive: bool = True,
                           parsed_args=None) -> None:
    """Parse elastic constants from ABACUS output and compute moduli."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Elastic Constants & Mechanical Properties ===[/bold cyan]")
    console.print()

    # Find elastic output
    elastic_path = None
    if parsed_args and parsed_args.file:
        elastic_path = Path(parsed_args.file)
    if elastic_path is None and args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                elastic_path = p
                break

    if elastic_path is None:
        # Auto-find
        candidates = (
            list(Path().glob("running*.log")) +
            list(Path().glob("elastic*.txt")) +
            list(Path().glob("OUT.*/running*.log")) +
            list(Path().glob("OUT.*/elastic*"))
        )
        if candidates:
            elastic_path = max(candidates, key=lambda p: p.stat().st_mtime)

    if elastic_path is None:
        console.print("[red]No elastic tensor data found.[/red]")
        console.print("[dim]Run ABACUS cell-relax calculation (cal_stress=1) to obtain elastic constants.[/dim]")
        return

    console.print(f"  [dim]Reading: {elastic_path}[/dim]")

    result = parse_elastic_tensor(elastic_path)
    if result is None:
        console.print("[red]Could not parse elastic tensor from the file.[/red]")
        console.print("[dim]Look for 'ELASTIC CONSTANTS' section in OUT.ABACUS/running_*.log[/dim]")
        return

    C = result["C"]
    console.print()
    console.print("  [bold]Elastic Tensor (Voigt, GPa):[/bold]")
    for i in range(6):
        row_str = "  ".join(f"{C[i, j]:10.3f}" for j in range(6))
        console.print(f"    {row_str}")

    # Compute properties
    props = compute_mechanical_properties(C)

    console.print()
    console.print("  [bold]Mechanical Properties (VRH):[/bold]")
    console.print(f"    Bulk Modulus    B = {props['B_H']:.2f} GPa (Voigt: {props['B_V']:.2f}, Reuss: {props['B_R']:.2f})")
    console.print(f"    Shear Modulus   G = {props['G_H']:.2f} GPa (Voigt: {props['G_V']:.2f}, Reuss: {props['G_R']:.2f})")
    console.print(f"    Young's Modulus E = {props['E_H']:.2f} GPa")
    console.print(f"    Poisson Ratio   ν = {props['nu_H']:.4f}")
    console.print(f"    Pugh Ratio   B/G = {props['pugh_ratio']:.3f} "
                  f"({'[green]ductile[/green]' if props['pugh_ratio'] > 1.75 else '[yellow]brittle[/yellow]'})")
    if abs(props["A_Zener"]) > 1e-10:
        console.print(f"    Zener Anisotropy = {props['A_Zener']:.4f} "
                      f"({'isotropic' if abs(props['A_Zener']-1.0)<0.05 else 'anisotropic'})")

    console.print()


# =============================================================================
# Task 752: EOS fitting
# =============================================================================


@task(1202, category="Mechanics", name="EOS Fitting",
      description="Fit Birch-Murnaghan equation of state from energy-volume data",
      cli_args=[
          {"name": "--file", "type": str, "default": None, "help": "Path to energy-volume data"},
      ])
def task_eos_fitting(args: list[str] | None = None, interactive: bool = True,
                     parsed_args=None) -> None:
    """Fit equation of state from energy-volume data.

    Can read from:
    - A user-provided data file (volume energy per line)
    - OUT.ABACUS energy output at different volumes
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Equation of State Fitting ===[/bold cyan]")
    console.print("[dim]Birch-Murnaghan (3rd order): E(V) → V₀, B₀, B₀'[/dim]")
    console.print()

    # Look for ev.dat style file
    ev_path = None
    if parsed_args and parsed_args.file:
        ev_path = Path(parsed_args.file)
    if ev_path is None and args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                ev_path = p
                break

    if ev_path is None:
        for name in ("ev.dat", "e_vs_v.dat", "energy_volume.dat"):
            p = Path(name)
            if p.exists():
                ev_path = p
                break

    if ev_path is None:
        console.print("[red]No energy-volume data file found.[/red]")
        console.print("[dim]Provide a file with lines: volume energy (per line).[/dim]")
        console.print("[dim]Example: abacuscopilot -task 1202 ev.dat[/dim]")
        return

    console.print(f"  [dim]Data file: {ev_path}[/dim]")

    # Read data
    data = np.loadtxt(ev_path)
    if data.ndim != 2 or data.shape[1] < 2:
        console.print("[red]File must have at least 2 columns: volume energy.[/red]")
        return

    volumes = data[:, 0]
    energies = data[:, 1]

    console.print(f"  Data points: {len(volumes)}")

    # Fit EOS
    eos = fit_eos_birch_murnaghan(volumes, energies)

    console.print()
    console.print("  [bold]Birch-Murnaghan (3rd order) fit:[/bold]")
    console.print(f"    V₀  = {eos['V0']:.4f}  (equilibrium volume)")
    console.print(f"    B₀  = {eos['B0']:.2f}  (bulk modulus, same units as energy/volume)")
    console.print(f"    B₀' = {eos['B0_prime']:.3f}  (pressure derivative)")
    console.print(f"    E₀  = {eos['E0']:.6f}  (equilibrium energy)")

    # Plot
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(volumes, energies, color="#1f77b4", s=30, zorder=5, label="Data")

    # Fitted curve
    V_fine = np.linspace(volumes.min() * 0.95, volumes.max() * 1.05, 200)

    def bm3(V, E0, V0, B0, B0p):
        eta = (V0 / V) ** (2.0 / 3.0)
        return E0 + 9.0 * V0 * B0 / 16.0 * (
            (eta - 1.0) ** 3 * B0p + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
        )

    E_fit = bm3(V_fine, eos["E0"], eos["V0"], eos["B0"], eos["B0_prime"])

    ax.plot(V_fine, E_fit, color="#d62728", linewidth=1.5, label="BM3 fit")
    ax.axvline(x=eos["V0"], color="gray", linestyle="--", linewidth=0.8, label=f"V₀ = {eos['V0']:.2f}")
    ax.set_xlabel("Volume")
    ax.set_ylabel("Energy")
    ax.set_title("Equation of State")
    ax.legend()

    save_name = "eos_fit.png"
    fig.savefig(save_name)
    console.print(f"  [green]✓ Plot saved to {save_name}[/green]")
    plt.close(fig)
    console.print()
