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


# (parse_elastic_tensor removed — replaced by stress–strain fitting in 1201)


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
    # Use the Voigt compliance matrix S = inv(C) directly in the sum
    # formula to avoid losing the shear-coefficient factor (9 vs 3).
    S = np.linalg.inv(C)  # Voigt compliance
    s11_22_33 = (S[0, 0] + S[1, 1] + S[2, 2]) / 3
    s12_13_23 = (S[0, 1] + S[0, 2] + S[1, 2]) / 3
    shear_sum = S[3, 3] + S[4, 4] + S[5, 5]      # sum, not average

    B_R = 1 / (3 * s11_22_33 + 6 * s12_13_23)

    # G_R = 15 / [4(S₁₁+S₂₂+S₃₃) − 4(S₁₂+S₂₃+S₁₃) + 3(S₄₄+S₅₅+S₆₆)]
    G_R = 15 / (4 * (S[0, 0] + S[1, 1] + S[2, 2])
                - 4 * (S[0, 1] + S[0, 2] + S[1, 2])
                + 3 * shear_sum)

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
    # Universal anisotropy index (Ranganathan & Ostoja-Starzewski)
    A_universal = 5 * G_V / G_R + B_V / B_R - 6 if G_R > 1e-10 else 0.0

    # Born stability criteria (cubic)
    born_1 = C[0, 0] - C[0, 1]  # C11 - C12 > 0
    born_2 = C[0, 0] + 2 * C[0, 1]  # C11 + 2C12 > 0
    born_3 = C[3, 3]  # C44 > 0
    born_stable = born_1 > 0 and born_2 > 0 and born_3 > 0

    # Vickers hardness (Tian 2012 model)
    k = G_H / B_H if B_H > 1e-10 else 0.0
    Hv_tian = 0.92 * (k ** 1.137) * (G_H ** 0.708) if k > 1e-10 else 0.0

    return {
        "B_V": B_V, "B_R": B_R, "B_H": B_H,
        "G_V": G_V, "G_R": G_R, "G_H": G_H,
        "E_H": E_H,
        "nu_H": nu_H,
        "pugh_ratio": pugh,
        "A_Zener": A_Zener,
        "A_universal": A_universal,
        "born_stable": born_stable,
        "born_1": born_1, "born_2": born_2, "born_3": born_3,
        "Hv_tian": Hv_tian,
    }


# =============================================================================
# Equation of State (EOS) fitting
# =============================================================================


def _bm3(V, E0, V0, B0, B0p):
    """3rd-order Birch-Murnaghan EOS.

    E(V) = E0 + 9V₀B₀/16 · {[(V₀/V)^(2/3)−1]³·B₀' + [(V₀/V)^(2/3)−1]²·[6−4(V₀/V)^(2/3)]}
    """
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + 9.0 * V0 * B0 / 16.0 * (
        (eta - 1.0) ** 3 * B0p + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def fit_eos_birch_murnaghan(
    volumes: np.ndarray,
    energies: np.ndarray,
) -> dict[str, float]:
    """Fit a 3rd-order Birch-Murnaghan equation of state.

    E(V) = E0 + 9V₀B₀/16 · {[(V₀/V)^(2/3)−1]³·B₀' + [(V₀/V)^(2/3)−1]²·[6−4(V₀/V)^(2/3)]}

    Args:
        volumes: Volume array (any unit).
        energies: Energy array (any unit).

    Returns:
        Dict with 'E0', 'V0', 'B0', 'B0_prime'.
    """
    from scipy.optimize import curve_fit

    imin = np.argmin(energies)
    E0_guess = energies[imin]
    V0_guess = volumes[imin]
    B0_guess = 100.0
    B0p_guess = 4.0

    try:
        popt, _ = curve_fit(
            _bm3, volumes, energies,
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
      description="Fit elastic tensor from stress–strain data (task.000–023 dirs)",
      cli_args=[])
def task_elastic_constants(args: list[str] | None = None, interactive: bool = True,
                           parsed_args=None) -> None:
    """Fit the 6×6 elastic tensor from ABACUS stress–strain data.

    Reads ``task.000/`` through ``task.023/`` directories (generated by
    task 111), extracts the final stress tensor from each
    ``OUT.ABACUS/running_scf.log``, and performs linear regression:

        σ_i = C_ij · ε_j     (Voigt notation)

    to obtain the elastic tensor C (GPa).  Then computes Voigt–Reuss–Hill
    bulk/shear moduli, Young's modulus, Poisson ratio, Pugh ratio, and
    Zener anisotropy.
    """
    import json

    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Elastic Constants (Stress–Strain Fitting) ===[/bold cyan]")
    console.print("[dim]Linear regression: σ = C·ε (Voigt notation)[/dim]")
    console.print()

    # Scan task directories
    task_dirs = sorted(Path().glob("task.[0-9][0-9][0-9]"))
    if not task_dirs:
        console.print("[red]No task.* directories found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 111' first to generate deformed structures.[/dim]")
        return

    strains_list: list[np.ndarray] = []   # each: (6,) Voigt strain
    stresses_list: list[np.ndarray] = []  # each: (6,) Voigt stress (kbar)

    for td in task_dirs:
        log_file = td / "OUT.ABACUS" / "running_relax.log"
        if not log_file.exists():
            log_file = td / "OUT.ABACUS" / "running_scf.log"
        strain_file = td / "strain.json"
        if not log_file.exists():
            console.print(f"  [yellow]![/yellow] {td.name}: no OUT.ABACUS/running_*.log")
            continue
        if not strain_file.exists():
            console.print(f"  [yellow]![/yellow] {td.name}: no strain.json")
            continue

        # Read strain vector
        try:
            strain_info = json.loads(strain_file.read_text())
            eps_voigt = np.array(strain_info["strain_voigt"], dtype=float)  # (6,)
        except Exception:
            console.print(f"  [yellow]![/yellow] {td.name}: bad strain.json")
            continue

        # Parse stress tensor from log
        stress_3x3 = _parse_stress_tensor(log_file)
        if stress_3x3 is None:
            console.print(f"  [yellow]![/yellow] {td.name}: stress not found in log")
            continue

        # Convert 3×3 stress to Voigt: (σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy)
        sigma_voigt = np.array([
            stress_3x3[0, 0], stress_3x3[1, 1], stress_3x3[2, 2],
            stress_3x3[1, 2], stress_3x3[0, 2], stress_3x3[0, 1],
        ])

        strains_list.append(eps_voigt)
        stresses_list.append(sigma_voigt)
        console.print(f"  [dim]{td.name}: ε=({eps_voigt[0]:+.4f},{eps_voigt[1]:+.4f},{eps_voigt[2]:+.4f},"
                      f"{eps_voigt[3]:+.4f},{eps_voigt[4]:+.4f},{eps_voigt[5]:+.4f}) "
                      f"σ=({sigma_voigt[0]:.2f},{sigma_voigt[1]:.2f},{sigma_voigt[2]:.2f},"
                      f"{sigma_voigt[3]:.2f},{sigma_voigt[4]:.2f},{sigma_voigt[5]:.2f}) kbar[/dim]")

    if len(strains_list) < 6:
        console.print(f"[red]Need at least 6 valid data points, got {len(strains_list)}.[/red]")
        return

    console.print(f"  Data points: {len(strains_list)}")

    # Build matrices: Σ (N×6), Ε (N×6)
    E_mat = np.array(strains_list)   # N×6
    # ABACUS stress convention: compressive = positive.  Standard elasticity
    # uses the opposite sign (tensile = positive).  Negate here.
    S_mat = -np.array(stresses_list)  # N×6, kbar, standard sign

    # Fit elastic tensor row-by-row:  σ_i = Σ_j C_ij · ε_j
    # For each row i of C:  C_i = (E^T·E)^(-1)·E^T·s_i
    C = np.zeros((6, 6))
    for i in range(6):
        C[i, :], _res, _rank, _sv = np.linalg.lstsq(E_mat, S_mat[:, i], rcond=None)

    # Convert kbar → GPa  (1 kbar = 0.1 GPa)
    KBAR_TO_GPA = 0.1
    C_gpa = C * KBAR_TO_GPA

    # Symmetrize (elastic tensor should be symmetric)
    C_gpa = (C_gpa + C_gpa.T) / 2.0

    console.print()
    console.print("  [bold]Elastic Tensor (Voigt, GPa):[/bold]")
    for i in range(6):
        row_str = "  ".join(f"{C_gpa[i, j]:10.3f}" for j in range(6))
        console.print(f"    {row_str}")

    # Compute mechanical properties
    props = compute_mechanical_properties(C_gpa)

    # ---- density from STRU (for sound velocity / Debye temperature) ----
    rho_gcm3: float | None = None
    n_atoms: int | None = None
    M_gmol: float | None = None
    try:
        from abacuscopilot.io.stru_file import _ATOMIC_MASSES, read_stru
        stru_path = Path("STRU")
        if not stru_path.exists():
            # Try task.000/STRU
            for sd in sorted(Path().glob("task.*")):
                s = sd / "STRU"
                if s.exists():
                    stru_path = s
                    break
        if stru_path.exists():
            s = read_stru(str(stru_path))
            n_atoms = len(s.atoms)
            vol_a3 = s.lattice.volume_angstrom
            # Atomic masses from IUPAC
            total_mass_amu = 0.0
            for sp in s.species_order:
                count = sum(1 for a in s.atoms if a.species == sp)
                total_mass_amu += count * _ATOMIC_MASSES.get(sp, 0.0)
            # ρ = mass / vol:  amu→g (×1.66054e-24), Å³→cm³ (×1e-24)
            AMU_TO_GRAM = 1.66053906660e-24
            rho_gcm3 = total_mass_amu * AMU_TO_GRAM / (vol_a3 * 1e-24)
            M_gmol = total_mass_amu  # g/mol (amu ≈ g/mol)
    except Exception:
        pass

    console.print()
    console.print("  [bold]Mechanical Properties (VRH):[/bold]")
    console.print(f"    Bulk Modulus    B = {props['B_H']:.2f} GPa (Voigt: {props['B_V']:.2f}, Reuss: {props['B_R']:.2f})")
    console.print(f"    Shear Modulus   G = {props['G_H']:.2f} GPa (Voigt: {props['G_V']:.2f}, Reuss: {props['G_R']:.2f})")
    console.print(f"    Young's Modulus E = {props['E_H']:.2f} GPa")
    console.print(f"    Poisson Ratio   ν = {props['nu_H']:.4f}")
    console.print(f"    Pugh Ratio   B/G = {props['pugh_ratio']:.3f} "
                  f"({'[green]ductile[/green]' if props['pugh_ratio'] > 1.75 else '[yellow]brittle[/yellow]'})")

    # Anisotropy
    console.print()
    console.print("  [bold]Anisotropy:[/bold]")
    if abs(props["A_Zener"]) > 1e-10:
        console.print(f"    Zener Anisotropy     = {props['A_Zener']:.4f} "
                      f"({'isotropic' if abs(props['A_Zener']-1.0)<0.05 else 'anisotropic'})")
    console.print(f"    Universal Anisotropy = {props['A_universal']:.4f} "
                  f"({'isotropic' if abs(props['A_universal'])<0.1 else 'anisotropic'})")

    # Hardness
    console.print()
    console.print("  [bold]Hardness (Tian 2012):[/bold]")
    console.print(f"    Hv = {props['Hv_tian']:.2f} GPa")

    # Born stability
    console.print()
    console.print("  [bold]Born Stability Criteria (cubic):[/bold]")
    console.print(f"    C₁₁ − C₁₂ = {props['born_1']:.2f} > 0  "
                  f"({'[green]✓[/green]' if props['born_1'] > 0 else '[red]✗[/red]'})")
    console.print(f"    C₁₁ + 2C₁₂ = {props['born_2']:.2f} > 0  "
                  f"({'[green]✓[/green]' if props['born_2'] > 0 else '[red]✗[/red]'})")
    console.print(f"    C₄₄ = {props['born_3']:.2f} > 0  "
                  f"({'[green]✓[/green]' if props['born_3'] > 0 else '[red]✗[/red]'})")
    console.print(f"    → {'[green]Mechanically stable[/green]' if props['born_stable'] else '[red]UNSTABLE[/red]'}")

    # Sound velocities & Debye temperature
    H_PLANCK = 6.62607015e-34
    KB = 1.380649e-23
    NA = 6.02214076e23
    _v_l = _v_t = _v_m = _theta_D = None

    if rho_gcm3 is not None and n_atoms and M_gmol:
        B_pa = props["B_H"] * 1e9
        G_pa = props["G_H"] * 1e9
        rho_kgm3 = rho_gcm3 * 1000.0

        _v_l = np.sqrt((B_pa + 4.0 * G_pa / 3.0) / rho_kgm3)
        _v_t = np.sqrt(G_pa / rho_kgm3)
        _v_m = (1.0 / 3.0 * (2.0 / _v_t**3 + 1.0 / _v_l**3)) ** (-1.0 / 3.0)
        debye_prefactor = H_PLANCK / KB * (3.0 * n_atoms * NA * rho_kgm3 / (4.0 * np.pi * M_gmol * 0.001)) ** (1.0 / 3.0)
        _theta_D = debye_prefactor * _v_m

        console.print()
        console.print("  [bold]Sound Velocities & Debye Temperature:[/bold]")
        console.print(f"    Density       ρ = {rho_gcm3:.3f} g/cm³")
        console.print(f"    Longitudinal  v_l = {_v_l:.1f} m/s")
        console.print(f"    Transverse    v_t = {_v_t:.1f} m/s")
        console.print(f"    Average       v_m = {_v_m:.1f} m/s")
        console.print(f"    Debye Temp    Θ_D = {_theta_D:.1f} K")

    # Save elastic tensor
    out_dat = "elastic_tensor.dat"
    with open(out_dat, "w") as f:
        f.write("# Elastic tensor (Voigt, GPa) — 6×6\n")
        for i in range(6):
            f.write("  ".join(f"{C_gpa[i, j]:12.6f}" for j in range(6)) + "\n")
    console.print()
    console.print(f"  [green]✓ {out_dat}[/green]")

    # Save mechanical properties
    prop_dat = "elastic_properties.dat"
    with open(prop_dat, "w") as f:
        f.write("# Mechanical Properties (VRH)\n")
        f.write(f"# Bulk Modulus    B = {props['B_H']:.2f} GPa "
                f"(Voigt: {props['B_V']:.2f}, Reuss: {props['B_R']:.2f})\n")
        f.write(f"# Shear Modulus   G = {props['G_H']:.2f} GPa "
                f"(Voigt: {props['G_V']:.2f}, Reuss: {props['G_R']:.2f})\n")
        f.write(f"# Young's Modulus E = {props['E_H']:.2f} GPa\n")
        f.write(f"# Poisson Ratio   ν = {props['nu_H']:.4f}\n")
        f.write(f"# Pugh Ratio   B/G = {props['pugh_ratio']:.3f} "
                f"({'ductile' if props['pugh_ratio'] > 1.75 else 'brittle'})\n")
        f.write(f"# Zener Anisotropy  = {props['A_Zener']:.4f}\n")
        f.write(f"# Universal Anisotropy = {props['A_universal']:.4f}\n")
        f.write(f"# Hardness (Tian)    = {props['Hv_tian']:.2f} GPa\n")
        f.write(f"# Born stable        = {'Yes' if props['born_stable'] else 'No'}\n")
        if _v_l is not None:
            f.write(f"# Density          ρ = {rho_gcm3:.3f} g/cm³\n")
            f.write(f"# v_longitudinal     = {_v_l:.1f} m/s\n")
            f.write(f"# v_transverse       = {_v_t:.1f} m/s\n")
            f.write(f"# v_average          = {_v_m:.1f} m/s\n")
            f.write(f"# Debye Temp      Θ_D = {_theta_D:.1f} K\n")
    console.print(f"  [green]✓ {prop_dat}[/green]")

    console.print()


def _parse_stress_tensor(filepath: Path) -> np.ndarray | None:
    """Parse the final TOTAL-STRESS tensor from an ABACUS log file.

    Returns 3×3 stress tensor in kbar, or None if not found.
    """
    if not filepath.exists():
        return None
    text = filepath.read_text()
    lines = text.split("\n")

    # Find the last occurrence of TOTAL-STRESS
    stress_blocks: list[list[str]] = []
    in_block = False
    block: list[str] = []
    for line in lines:
        if re.search(r"TOTAL-STRESS\s*\(KBAR\)", line, re.IGNORECASE):
            if block:
                stress_blocks.append(block)
            block = []
            in_block = True
            continue
        if in_block:
            if block and (not line.strip() or re.match(r"^\s*$", line)):
                stress_blocks.append(block)
                in_block = False
                continue
            # Match lines with 3 numbers
            parts = line.split()
            nums = []
            for p in parts:
                try:
                    nums.append(float(p))
                except ValueError:
                    break
            if len(nums) >= 3:
                block.append(nums[:3])

    if block:
        stress_blocks.append(block)

    if not stress_blocks:
        return None

    # Use the last block
    last = stress_blocks[-1]
    if len(last) >= 3:
        return np.array(last[:3])
    return None


# =============================================================================
# Task 752: EOS fitting
# =============================================================================


@task(1202, category="Mechanics", name="EOS Fitting",
      description="Extract E-V data from scale_* dirs, fit Birch-Murnaghan EOS, plot",
      cli_args=[
          {"name": "--file", "type": str, "default": None, "help": "Path to energy-volume data file (ev.dat)"},
          {"name": "--no-plot", "action": "store_true", "default": False,
           "help": "Skip plotting"},
      ])
def task_eos_fitting(args: list[str] | None = None, interactive: bool = True,
                     parsed_args=None) -> None:
    """Fit equation of state from energy-volume data.

    Two-stage workflow:

    1. **Extract** — scan ``scale_*/`` directories, read STRU (volume)
       and ``OUT.ABACUS/running_scf.log`` (!FINAL_ETOT_IS), write ``ev.dat``.

    2. **Fit** — read ``ev.dat``, fit 3rd-order Birch-Murnaghan EOS,
       write ``eos_fit.dat`` (fitted curve), print V₀/B₀/B₀'/E₀.

    If ``--file`` is given, skip stage 1 and use that file directly.
    """
    import matplotlib.pyplot as plt
    from abacuscopilot.plotting.style import load_style_from_config
    from abacuscopilot.io.stru_file import read_stru

    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Equation of State Fitting ===[/bold cyan]")
    console.print("[dim]Birch-Murnaghan (3rd order): E(V) → V₀, B₀, B₀'[/dim]")
    console.print()

    # ---- stage 1: auto-extract from scale_* directories ----
    ev_path = None
    if parsed_args and parsed_args.file:
        ev_path = Path(parsed_args.file)
    if ev_path is None and args:
        for arg in args:
            p = Path(arg)
            if p.exists() and p.suffix in (".dat", ".txt"):
                ev_path = p
                break

    _ev_dirs: list[str] = []  # scale dir names (parallel to volumes after extraction)

    if ev_path is None:
        # Auto-extract from scale_* directories
        scale_dirs = sorted(Path().glob("scale_*"))
        if scale_dirs:
            console.print("[bold]Step 1: Extract E-V from scale directories[/bold]")
            console.print()

            points: list[tuple[float, float, str]] = []  # (vol, energy, dir_name)
            for sd in scale_dirs:
                stru_file = sd / "STRU"
                log_file = sd / "OUT.ABACUS" / "running_scf.log"
                if not stru_file.exists() or not log_file.exists():
                    console.print(f"  [yellow]![/yellow] {sd.name}: missing STRU or OUT.ABACUS/running_scf.log")
                    continue

                # Volume from STRU
                try:
                    s = read_stru(str(stru_file))
                    vol = s.lattice.volume_angstrom
                except Exception:
                    console.print(f"  [yellow]![/yellow] {sd.name}: failed to read STRU")
                    continue

                # Energy from SCF log
                try:
                    text = log_file.read_text()
                    m = re.search(r"!FINAL_ETOT_IS\s+([\-\d\.Ee+]+)", text)
                    if not m:
                        console.print(f"  [yellow]![/yellow] {sd.name}: !FINAL_ETOT_IS not found in log")
                        continue
                    energy_ev = float(m.group(1))
                except Exception:
                    console.print(f"  [yellow]![/yellow] {sd.name}: failed to read energy")
                    continue

                points.append((vol, energy_ev, sd.name))
                console.print(f"  [dim]{sd.name}: V={vol:.4f} Å³, E={energy_ev:.4f} eV[/dim]")

            if not points:
                console.print("[red]No valid E-V data extracted from scale_* directories.[/red]")
                return

            # Sort by volume
            points.sort(key=lambda x: x[0])

            _ev_dirs = [p[2] for p in points]  # dir names in sorted order
            ev_path = Path("ev.dat")
            with open(ev_path, "w") as f:
                f.write("# volume(A^3)  energy(eV)\n")
                for v, e, _d in points:
                    f.write(f"{v:.6f}  {e:.8f}\n")
            console.print(f"  [green]✓ ev.dat[/green] ({len(points)} points)")
            console.print()

    # ---- stage 2: fit EOS ----
    if ev_path is None:
        # Fallback: scan for ev.dat etc.
        for name in ("ev.dat", "e_vs_v.dat", "energy_volume.dat"):
            p = Path(name)
            if p.exists():
                ev_path = p
                break

    if ev_path is None:
        console.print("[red]No energy-volume data file found.[/red]")
        console.print("[dim]Run this task in a directory with scale_*/ subdirs,[/dim]")
        console.print("[dim]or provide: abacuscopilot -task 1202 --file ev.dat[/dim]")
        return

    console.print(f"  [dim]Data file: {ev_path}[/dim]")

    data = np.loadtxt(ev_path)
    if data.ndim != 2 or data.shape[1] < 2:
        console.print("[red]File must have at least 2 columns: volume energy.[/red]")
        return

    volumes = data[:, 0]
    energies = data[:, 1]
    console.print(f"  Data points: {len(volumes)}")

    # Fit
    eos = fit_eos_birch_murnaghan(volumes, energies)
    if eos["B0"] == 0.0:
        console.print("[red]EOS fitting failed.[/red]")
        return

    # Convert B0 from eV/Å³ → GPa
    b0_ev_a3 = eos["B0"]
    b0_gpa = b0_ev_a3 * 160.2177

    # Equilibrium energy per atom (estimate)
    n_atoms: int | None = None
    try:
        s0 = read_stru("STRU")
        n_atoms = len(s0.atoms)
    except Exception:
        for sd in sorted(Path().glob("scale_*")):
            try:
                s0 = read_stru(str(sd / "STRU"))
                n_atoms = len(s0.atoms)
                break
            except Exception:
                continue

    console.print()
    console.print("  [bold]Birch-Murnaghan (3rd order) fit:[/bold]")
    console.print(f"    V₀  = {eos['V0']:.4f} Å³  (equilibrium volume)")
    console.print(f"    B₀  = {b0_gpa:.2f} GPa  (bulk modulus)")
    console.print(f"    B₀' = {eos['B0_prime']:.3f}  (pressure derivative)")
    console.print(f"    E₀  = {eos['E0']:.6f} eV  (equilibrium energy)")
    if n_atoms:
        console.print(f"    E₀/atom = {eos['E0']/n_atoms:.6f} eV")
    console.print()

    # Locate the data point closest to V0
    idx_min = np.argmin(np.abs(volumes - eos["V0"]))
    console.print(f"  Minimum-energy data point: V = {volumes[idx_min]:.4f} Å³, "
                  f"E = {energies[idx_min]:.6f} eV")

    # ---- fitted curve dat ----
    V_fine = np.linspace(volumes.min() * 0.95, volumes.max() * 1.05, 200)
    E_fit = np.array([_bm3(v, eos["E0"], eos["V0"], eos["B0"], eos["B0_prime"])
                      for v in V_fine])

    dat_fit = Path("eos_fit.dat")
    with open(dat_fit, "w") as f:
        f.write("# volume(A^3)  energy_fit(eV)\n")
        for v, e in zip(V_fine, E_fit):
            f.write(f"{v:.6f}  {e:.8f}\n")
    console.print(f"  [green]✓ eos_fit.dat[/green] (fitted curve, {len(V_fine)} points)")

    # ---- min.STRU: structure closest to equilibrium volume ----
    import shutil
    idx_min = np.argmin(np.abs(volumes - eos["V0"]))
    min_stru_src = None
    if _ev_dirs and idx_min < len(_ev_dirs):
        min_stru_src = Path(_ev_dirs[idx_min]) / "STRU"
    if min_stru_src is None or not min_stru_src.exists():
        # Fallback: scan scale_* directories for matching volume
        for sd in sorted(Path().glob("scale_*")):
            try:
                s_test = read_stru(str(sd / "STRU"))
                if abs(s_test.lattice.volume_angstrom - volumes[idx_min]) < 0.01:
                    min_stru_src = sd / "STRU"
                    break
            except Exception:
                continue
    if min_stru_src and min_stru_src.exists():
        shutil.copy2(min_stru_src, "min.STRU")
        console.print(f"  [green]✓ min.STRU[/green] (from {min_stru_src.parent.name}/)")
    else:
        console.print("  [yellow]![/yellow] min.STRU: source not found")

    # ---- plot ----
    do_plot = True
    if parsed_args and parsed_args.no_plot:
        do_plot = False
    elif interactive:
        from abacuscopilot.console_utils import _prompt_choice
        answer = _prompt_choice(console, "Generate EOS plot?", ["Yes", "No"], "Yes")
        do_plot = "Yes" in answer

    if do_plot:
        load_style_from_config()
        fig, ax = plt.subplots(figsize=(8, 6))

        # Data: scatter
        ax.scatter(volumes, energies, color="#1f77b4", s=40, zorder=5,
                   edgecolors="white", linewidths=0.5, label="Data")

        # Fitted curve: line
        ax.plot(V_fine, E_fit, color="#d62728", linewidth=1.5, zorder=4,
                label="BM3 fit")

        # Equilibrium marker (open circle)
        ax.scatter([eos["V0"]], [eos["E0"]], marker="o", facecolors="none",
                   edgecolors="#d62728", s=120, zorder=6, linewidths=1.5)
        ax.annotate(
            f"V0 = {eos['V0']:.2f} A^3\nE0 = {eos['E0']:.2f} eV",
            xy=(eos["V0"], eos["E0"]),
            xytext=(15, -25), textcoords="offset points",
            fontsize=8, color="#d62728",
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.8),
        )

        ax.set_xlabel("Volume (A^3)")
        ax.set_ylabel("Energy (eV)")
        ax.set_title("Equation of State — Birch-Murnaghan Fit")
        ax.legend(loc="upper left")

        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
            spine.set_visible(True)
        ax.tick_params(axis="both", direction="out")

        fig.tight_layout(pad=1.2)
        out_png = "eos_fit.png"
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        console.print(f"  [green]✓ Plot: {out_png}[/green]")

    console.print()
