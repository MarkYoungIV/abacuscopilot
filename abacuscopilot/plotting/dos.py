"""Density of states plotting for ABACUS output.

Reads DOS_*.dat and PDOS_*.dat files and produces publication-quality
total DOS, projected DOS, and combined DOS+band plots.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Default color cycle for PDOS
PDOS_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
]


def read_dos_dat(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Read an ABACUS DOS_*.dat file.

    The DOS file format is:
        line 1: nedos e_fermi
        Then: nedos lines of (energy, dos_tot, dos_s, dos_p, dos_d, ...)

    Or for older versions:
        line 1: nedos
        Then: nedos lines of (energy, dos_tot)

    Args:
        filepath: Path to the DOS file.

    Returns:
        Tuple of (energies, dos, e_fermi).
        energies: shape (nedos,) in eV.
        dos: shape (nedos,) total DOS in states/eV.
        e_fermi: Fermi energy in eV.
    """
    filepath = Path(filepath)

    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip()]

    # Auto-detect: no header (data starts with energy values) vs
    # old header format (first line: nedos [e_fermi]).
    first = lines[0].split()
    try:
        n_check = int(first[0])
        # Valid header line has a positive integer line count (nedos > 0).
        # Negative values are energy data (TDOS.dat format), not a header.
        if n_check > 0:
            if len(first) >= 2:
                nedos = n_check
                e_fermi = float(first[1])
            else:
                nedos = n_check
                e_fermi = 0.0
            data_start = 1
        else:
            raise ValueError("not a header")
    except (ValueError, IndexError):
        # No header — first line is data
        nedos = len(lines)
        e_fermi = 0.0
        data_start = 0

    data = []
    for line in lines[data_start:]:
        data.extend(float(x) for x in line.split())

    data = np.array(data)
    if nedos <= 0 or len(data) == 0:
        return np.array([]), np.array([]), 0.0

    ncols = max(len(data) // nedos, 1)  # at least energy column

    # Column 0 = energy, column 1 = DOS (total), column 2 = integrated (optional)
    energies = data[0:nedos * ncols:ncols]
    dos_tot = data[1:nedos * ncols:ncols] if ncols >= 2 else np.zeros(nedos)

    return energies, dos_tot, e_fermi


def _read_pdos_xml(filepath: Path) -> dict[str, np.ndarray]:
    """Parse ABACUS XML PDOS file (v3.7+).

    Returns dict with: energies, e_fermi, and per-species keys like
    'Li_s', 'Li_p', 'Li_total', ... and 'total'.
    """
    import xml.etree.ElementTree as ET
    from io import StringIO

    tree = ET.parse(str(filepath))
    root = tree.getroot()

    # Read energy grid
    ev = root.find("energy_values")
    energies = np.loadtxt(StringIO(ev.text))

    # Get e_fermi from running_scf.log or set to 0
    e_fermi = 0.0

    # Sum per species+orbital (l)
    l_names = ["s", "p", "d", "f", "g"]
    species_total: dict[str, np.ndarray] = {}
    species_shells: dict[str, dict[str, np.ndarray]] = {}

    for orb in root.findall("orbital"):
        sp = orb.get("species", "unknown")
        l_val = int(orb.get("l", 0))
        data_text = orb.find("data")
        if data_text is None or data_text.text is None:
            continue
        arr = np.loadtxt(StringIO(data_text.text))  # (nedos,)

        key = f"{sp}_{l_names[l_val]}" if l_val < len(l_names) else f"{sp}_l{l_val}"
        if key not in species_shells.setdefault(sp, {}):
            species_shells[sp][key] = arr
        else:
            species_shells[sp][key] += arr

        if sp not in species_total:
            species_total[sp] = arr.copy()
        else:
            species_total[sp] += arr

    total = np.zeros_like(energies)
    result = {"energies": energies, "e_fermi": e_fermi}
    for sp in species_total:
        result[f"{sp}_total"] = species_total[sp]
        total += species_total[sp]
    result["total"] = total

    # Also add per-shell keys
    for sp, shells in species_shells.items():
        for key, arr in shells.items():
            result[key] = arr

    return result


def read_pdos_dat(filepath: str | Path) -> dict[str, np.ndarray]:
    """Read an ABACUS PDOS file (old .dat or new XML format).

    Args:
        filepath: Path to the PDOS file.

    Returns:
        Dict with keys: 'energies', 'e_fermi', and per-orbital keys like
        's', 'p', 'd', 'f', 'total', and per-species like 'Si_s', 'Si_p', etc.
        Each value is a 1D numpy array.
    """
    filepath = Path(filepath)

    # Try XML format first (ABACUS v3.7+)
    try:
        with open(filepath) as f:
            first = f.read(100)
        if first.strip().startswith("<pdos"):
            return _read_pdos_xml(filepath)
    except Exception:
        pass

    # Legacy text format
    with open(filepath) as f:
        lines = f.readlines()

    header = lines[0].split()
    nedos = int(header[0])
    e_fermi = float(header[1]) if len(header) >= 2 else 0.0

    data = []
    for line in lines[1:]:
        data.extend(float(x) for x in line.split())

    data = np.array(data)
    ncols = len(data) // nedos if nedos > 0 else 1
    reshaped = data[:nedos * ncols].reshape(nedos, ncols)

    result = {"energies": reshaped[:, 0], "e_fermi": e_fermi}

    col_names = ["energy", "s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2",
                 "fy3x2", "fxyz", "fz3", "fzx2", "fx3", "tot"]
    for i in range(1, min(ncols, len(col_names))):
        if i < reshaped.shape[1]:
            result[col_names[i]] = reshaped[:, i]

    return result


def plot_dos(
    dos_path: str | Path,
    *,
    e_range: tuple[float, float] | None = None,
    color: str = "#1f77b4",
    fill: bool = True,
    fill_alpha: float = 0.3,
    linewidth: float = 1.5,
    label: str = "Total DOS",
    title: str | None = None,
    show_fermi: bool = True,
    save: str | Path | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot total density of states from ABACUS output.

    Args:
        dos_path: Path to DOS_*.dat file.
        e_range: (emin, emax) energy range in eV relative to E_Fermi.
        color: Line/fill color.
        fill: Whether to fill under the DOS curve.
        fill_alpha: Transparency of fill.
        linewidth: Line width.
        label: Legend label.
        title: Plot title.
        show_fermi: Whether to show Fermi level.
        save: Path to save figure.
        ax: Existing Axes to plot on.

    Returns:
        matplotlib Axes.
    """
    from abacuscopilot.plotting.style import load_style_from_config

    load_style_from_config()

    if ax is None:
        _, ax = plt.subplots()

    energies, dos, e_fermi = read_dos_dat(dos_path)

    if energies.size == 0:
        ax.text(0.5, 0.5, "No DOS data found", transform=ax.transAxes, ha="center")
        return ax

    energies_shifted = energies - e_fermi

    ax.plot(energies_shifted, dos, color=color, linewidth=linewidth, label=label)

    if fill:
        ax.fill_between(energies_shifted, 0, dos, color=color, alpha=fill_alpha)

    if show_fermi:
        ax.axvline(x=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7,
                   label=f"E$_F$ = {e_fermi:.3f} eV")

    ax.set_xlabel("E − E$_F$ (eV)")
    ax.set_ylabel("DOS (states/eV)")

    if e_range:
        ax.set_xlim(*e_range)
    else:
        margin = 1.0
        ax.set_xlim(energies_shifted.min() - margin, energies_shifted.max() + margin)

    ax.set_ylim(bottom=0)

    if label:
        ax.legend(loc="upper right")

    if title:
        ax.set_title(title)

    _style_ax(ax)

    if save:
        ax.figure.tight_layout(pad=1.2)
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")

    return ax


def _style_ax(ax):
    """Apply consistent style: full frame + outward ticks."""
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")


def plot_pdos(
    pdos_path: str | Path,
    *,
    orbitals: list[str] | None = None,
    colors: list[str] | None = None,
    e_range: tuple[float, float] | None = None,
    fill: bool = True,
    fill_alpha: float = 0.3,
    linewidth: float = 1.2,
    title: str | None = None,
    show_fermi: bool = True,
    save: str | Path | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot projected density of states.

    Args:
        pdos_path: Path to PDOS_*.dat file.
        orbitals: List of orbital names to plot (e.g., ['s', 'p', 'd']).
                 If None, plots all available orbitals.
        colors: Colors for each orbital.
        e_range: Energy range.
        fill: Whether to fill under curves.
        fill_alpha: Fill transparency.
        linewidth: Line width.
        title: Plot title.
        show_fermi: Show Fermi level.
        save: Path to save figure.
        ax: Existing Axes to plot on.

    Returns:
        matplotlib Axes.
    """
    from abacuscopilot.plotting.style import load_style_from_config

    load_style_from_config()

    if ax is None:
        _, ax = plt.subplots()

    data = read_pdos_dat(pdos_path)

    if "energies" not in data or data["energies"].size == 0:
        ax.text(0.5, 0.5, "No PDOS data found", transform=ax.transAxes, ha="center")
        return ax

    energies = data["energies"]
    e_fermi = data.get("e_fermi", 0.0)
    energies_shifted = energies - e_fermi

    if orbitals is None:
        # Auto-detect orbitals
        orbitals = [k for k in data if k not in ("energies", "e_fermi")
                    and isinstance(data[k], np.ndarray)
                    and data[k].shape == energies.shape]

    if colors is None:
        colors = PDOS_COLORS[:len(orbitals)]

    for i, orb in enumerate(orbitals[:len(colors)]):
        color = colors[i]
        ax.plot(energies_shifted, data[orb], color=color, linewidth=linewidth,
                label=orb)
        if fill:
            ax.fill_between(energies_shifted, 0, data[orb], color=color,
                           alpha=fill_alpha)

    if show_fermi:
        ax.axvline(x=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

    ax.set_xlabel("E − E$_F$ (eV)")
    ax.set_ylabel("PDOS (states/eV)")
    ax.set_ylim(bottom=0)

    if e_range:
        ax.set_xlim(*e_range)
    else:
        margin = 1.0
        ax.set_xlim(energies_shifted.min() - margin, energies_shifted.max() + margin)

    if orbitals:
        ax.legend(loc="upper right")

    if title:
        ax.set_title(title)

    _style_ax(ax)

    if save:
        ax.figure.tight_layout(pad=1.2)
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")

    return ax


def plot_dos_bands_combined(
    bands_path: str | Path,
    dos_path: str | Path,
    kpt: object | str | Path | None = None,
    *,
    e_range: tuple[float, float] | None = None,
    band_color: str = "#1f77b4",
    dos_color: str = "#ff7f0e",
    band_alpha: float = 0.7,
    band_lw: float = 1.0,
    dos_lw: float = 1.5,
    title: str | None = None,
    save: str | Path | None = None,
) -> plt.Figure:
    """Combine band structure and DOS side by side.

    Args:
        bands_path: Path to BANDS_*.dat.
        dos_path: Path to DOS_*.dat.
        kpt: KPT information for band path labels.
        e_range: Shared energy range.
        band_color: Band line color.
        dos_color: DOS line color.
        band_alpha: Band transparency.
        band_lw: Band linewidth.
        dos_lw: DOS linewidth.
        title: Overall figure title.
        save: Path to save figure.

    Returns:
        matplotlib Figure.
    """
    from abacuscopilot.plotting.bands import read_bands_with_kpt
    from abacuscopilot.plotting.style import load_style_from_config

    load_style_from_config()

    fig, (ax_bands, ax_dos) = plt.subplots(1, 2, figsize=(12, 6),
                                            gridspec_kw={"width_ratios": [2, 1]})

    # === Band structure (left) ===
    k_dists, energies, e_fermi_b, labels, label_positions = read_bands_with_kpt(
        bands_path, kpt
    )

    if energies.size > 0:
        energies_shifted = energies - e_fermi_b
        nbands = energies.shape[0]
        for ib in range(nbands):
            ax_bands.plot(k_dists, energies_shifted[ib], color=band_color,
                         alpha=band_alpha, linewidth=band_lw)

        ax_bands.axhline(y=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

        if labels and label_positions:
            ax_bands.set_xticks(label_positions)
            ax_bands.set_xticklabels(labels)
            for pos in label_positions:
                ax_bands.axvline(x=pos, color="gray", linestyle=":", linewidth=0.5,
                                alpha=0.5)

    ax_bands.set_ylabel("E − E$_F$ (eV)")
    ax_bands.set_title("Band Structure")

    # === DOS (right) ===
    energies_dos, dos, e_fermi_d = read_dos_dat(dos_path)

    if energies_dos.size > 0:
        energies_d_shifted = energies_dos - e_fermi_d
        ax_dos.plot(energies_d_shifted, dos, color=dos_color, linewidth=dos_lw)
        ax_dos.fill_between(energies_d_shifted, 0, dos, color=dos_color, alpha=0.2)
        ax_dos.axhline(y=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

    ax_dos.set_xlabel("DOS (states/eV)")
    ax_dos.set_title("Density of States")
    ax_dos.set_yticklabels([])
    ax_dos.set_ylim(ax_bands.get_ylim())

    # Shared energy range
    e_fermi = e_fermi_b if e_fermi_b != 0.0 else e_fermi_d
    if e_range:
        ax_bands.set_ylim(*e_range)
        ax_dos.set_ylim(*e_range)

    _style_ax(ax_dos)
    _style_ax(ax_bands)

    fig.tight_layout(pad=1.2)

    if title:
        fig.suptitle(title, y=1.02)

    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    return fig
