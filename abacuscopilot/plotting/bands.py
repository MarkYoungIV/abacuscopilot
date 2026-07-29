"""Band structure plotting for ABACUS output.

Reads BANDS_*.dat files and produces publication-quality band structure
plots with optional fat-band and projected-band capabilities.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from abacuscopilot.core.models import KPoints


def read_bands_dat(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Read an ABACUS BANDS_*.dat file.

    The BANDS file format is:
        line 1: nbands nkpoints
        line 2: e_fermi (optional, may be absent or on line 1)
        Then: nkpoints blocks, each with nbands energies

    Args:
        filepath: Path to the BANDS file.

    Returns:
        Tuple of (k_distances, energies, e_fermi).
        k_distances: shape (nkpoints,) cumulative k-distance.
        energies: shape (nbands, nkpoints) in eV.
        e_fermi: Fermi energy in eV (0.0 if not found).
    """
    filepath = Path(filepath)

    with open(filepath) as f:
        lines = f.readlines()

    # Parse header — if the first line has > 10 values it's actually data (ABACUS v3.10+)
    first_line = lines[0].split()
    ncols_first = len(first_line)
    header = [float(x) for x in first_line]
    e_fermi = 0.0
    nbands = 0
    data_start = 1

    if ncols_first <= 10:
        # Looks like a header: (nbands, nkpoints) or (e_fermi, nkpts) or (nbands, nkpts, e_fermi)
        if len(header) == 1:
            nkpts = int(header[0])
        elif len(header) == 2:
            v1, v2 = header[0], header[1]
            if v1 > 50 and v2 > 1:
                nbands, nkpts = int(v1), int(v2)
            elif v2 < 100:
                e_fermi, nkpts = v1, int(v2)
            else:
                nbands, nkpts = int(v1), int(v2)
        elif len(header) >= 3:
            nbands, nkpts = int(header[0]), int(header[1])
            e_fermi = header[2]

        # Older ABACUS: e_fermi may be on a separate line before data
        if len(header) <= 2 and len(lines) > 1:
            line2 = lines[1].split()
            if len(line2) == 1:
                try:
                    e_fermi = float(line2[0])
                    data_start = 2
                except ValueError:
                    pass
    else:
        # First line is data (no header) — include it in the data
        data_start = 0

    # Parse data — ABACUS v3.10+ writes one row per k-point:
    #   kpt_index  k_distance  E(band1)  E(band2)  ...  E(bandN)
    raw = []
    for line in lines[data_start:]:
        parts = line.split()
        if parts:
            raw.append([float(x) for x in parts])

    if not raw:
        return np.array([]), np.array([]), e_fermi

    data = np.array(raw)           # (nkpts, 2 + nbands)
    # Column 0 = k-point index (1-based), column 1 = k-distance, columns 2+ = band energies
    k_distances = data[:, 1]
    if nbands == 0:
        nbands = data.shape[1] - 2
    energies = data[:, 2:2 + nbands].T   # (nbands, nkpoints)

    # Auto-detect e_fermi from logs if not found in header
    if e_fermi == 0.0:
        e_fermi = _find_e_fermi(filepath)

    return k_distances, energies, e_fermi


def read_bands_with_kpt(
    bands_path: str | Path,
    kpt: KPoints | str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, float, list[str], list[float]]:
    """Read BANDS file and combine with KPT for high-symmetry point labels.

    Args:
        bands_path: Path to BANDS_*.dat file.
        kpt: KPoints object or path to KPT file for line-mode path info.
             If None, attempts to find KPT in the same directory.

    Returns:
        Tuple of (k_distances, energies, e_fermi, labels, label_positions).
    """
    k_distances_raw, energies, e_fermi = read_bands_dat(bands_path)

    if energies.size == 0:
        return k_distances_raw, energies, e_fermi, [], []

    nkpts = energies.shape[1]
    labels = []
    label_positions = []

    # Load KPT if provided
    if isinstance(kpt, (str, Path)):
        from abacuscopilot.io.kpt_file import read_kpt
        kpt = read_kpt(kpt)
    elif kpt is None:
        # Search BANDS directory and parent (KPT is typically in the working dir)
        bands_dir = Path(bands_path).parent
        search_dirs = [bands_dir, bands_dir.parent]
        for search_dir in search_dirs:
            for kpt_name in ("KPT", "KPT_LINE", "KPT_BAND"):
                kpt_path = search_dir / kpt_name
                if kpt_path.exists():
                    try:
                        from abacuscopilot.io.kpt_file import read_kpt
                        kpt = read_kpt(kpt_path)
                        break
                    except Exception:
                        pass
            if kpt is not None:
                break

    if kpt is not None and kpt.mode in ("line", "line_cartesian"):
        # Use k-distances from the BANDS file (column 1 = ABACUS-computed path length)
        k_distances = k_distances_raw

        # Map labels to k-distances based on cumulative npoints per segment.
        # Connector segments (npoints == 1) mark path discontinuities
        # (seekpath "Z|X", "R|M" breaks).  Merge them with the following
        # label and place the merged label at the discontinuity position.
        cum_kpt = 0
        pending_break_label = None
        pending_break_pos = None
        for i, seg in enumerate(kpt.line_path):
            npts = seg.get("npoints", 20)
            label = seg.get("label", "")
            if npts == 1 and label:
                # Connector segment — defer, will merge with next label
                pending_break_label = label
                pending_break_pos = float(k_distances[cum_kpt]) if cum_kpt < len(k_distances) else None
            else:
                if label:
                    if pending_break_label and pending_break_pos is not None:
                        labels.append(f"{pending_break_label} ∥ {label}")
                        label_positions.append(pending_break_pos)
                        pending_break_label = None
                        pending_break_pos = None
                    else:
                        labels.append(label)
                        label_positions.append(float(k_distances[cum_kpt]))
            cum_kpt += npts

        # Add final label (not a connector)
        if kpt.line_path:
            end_label = kpt.line_path[-1].get("end_label", "")
            if end_label and cum_kpt < len(k_distances):
                labels.append(end_label)
                label_positions.append(float(k_distances[min(cum_kpt, len(k_distances) - 1)]))
    else:
        k_distances = k_distances_raw if k_distances_raw.size > 0 else np.arange(nkpts, dtype=float)

    return k_distances, energies, e_fermi, labels, label_positions


def _insert_path_breaks(k_dists: np.ndarray, energies: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Insert NaN breaks at k-path discontinuities (seekpath ``Z|X``-style jumps).

    When consecutive k-points share the same k-distance, the second one is a
    path-discontinuity connector.  We *drop* that duplicate point entirely and
    replace it with NaN so matplotlib breaks the line both before AND after
    the jump, preventing a meaningless straight line across the gap (e.g. the
    Z→X or R→M segments that have zero intermediate k-points).
    """
    diffs = np.diff(k_dists)
    breaks = np.where(diffs < 1e-10)[0]  # indices i where k_dists[i]==k_dists[i+1]
    if len(breaks) == 0:
        return k_dists, energies

    nkpts = len(k_dists)
    nbands = energies.shape[0]
    # Each break removes the duplicate point (brk+1) and inserts a NaN in its place.
    # Net length change = 0 (remove 1, add 1).
    new_n = nkpts
    new_k = np.empty(new_n, dtype=k_dists.dtype)
    new_e = np.empty((nbands, new_n), dtype=energies.dtype)

    src = 0
    dst = 0
    for brk in breaks:
        # Copy all points up to and including brk (the last point before the jump)
        ncopy = brk - src + 1
        new_k[dst:dst + ncopy] = k_dists[src:src + ncopy]
        new_e[:, dst:dst + ncopy] = energies[:, src:src + ncopy]
        dst += ncopy
        # Replace the duplicate point (src + ncopy = brk+1) with NaN
        new_k[dst] = np.nan
        new_e[:, dst] = np.nan
        dst += 1
        # Skip the duplicate point in the source arrays
        src = brk + 2

    # Copy remaining points (if any)
    remaining = nkpts - src
    if remaining > 0:
        new_k[dst:dst + remaining] = k_dists[src:]
        new_e[:, dst:dst + remaining] = energies[:, src:]

    return new_k, new_e


# Gap width between k-path segments at a discontinuity (in normalised units).
_PATH_BREAK_GAP = 0.03


def _compress_path_breaks(
    k_dists: np.ndarray,
    energies: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[int, float]]:
    """Compress NaN-separated k-path segments so large gaps collapse to a
    small visual break, while preserving relative distances within each segment.

    Returns ``(new_k, new_e, label_map)`` where *label_map* maps old
    k-distance positions to new compressed positions.
    """
    # Identify continuous segments (runs of non-NaN values)
    mask = ~np.isnan(k_dists)
    if np.all(mask):
        return k_dists, energies, {}

    segments = []
    start = 0
    n = len(k_dists)
    while start < n:
        while start < n and not mask[start]:
            start += 1
        if start >= n:
            break
        end = start
        while end < n and mask[end]:
            end += 1
        segments.append((start, end))
        start = end + 1

    if not segments:
        return k_dists, energies, {}

    nbands = energies.shape[0]

    # Build compressed arrays: each segment's original x-range is normalised
    # to its original span, then shifted so segments are adjacent.
    new_n = sum(end - start for start, end in segments)
    new_k = np.empty(new_n, dtype=k_dists.dtype)
    new_e = np.empty((nbands, new_n), dtype=energies.dtype)
    label_map: dict[int, float] = {}

    dst = 0
    prev_k_end = 0.0
    for seg_idx, (seg_start, seg_end) in enumerate(segments):
        seg_len = seg_end - seg_start
        seg_k = k_dists[seg_start:seg_end]        # original x-values
        seg_e = energies[:, seg_start:seg_end]     # (nbands, seg_len)

        # Normalise: map [seg_k[0], seg_k[-1]] → [0, 1] then shift
        seg_k_min = seg_k[0]
        seg_k_max = seg_k[-1]
        seg_span = seg_k_max - seg_k_min
        if seg_span > 1e-12:
            seg_norm = (seg_k - seg_k_min) / seg_span
        else:
            seg_norm = np.zeros(seg_len)

        if seg_idx > 0:
            prev_k_end += _PATH_BREAK_GAP
        new_seg = seg_norm * seg_span + prev_k_end

        new_k[dst:dst + seg_len] = new_seg
        new_e[:, dst:dst + seg_len] = seg_e
        label_map[seg_k_min] = float(new_seg[0])

        dst += seg_len
        prev_k_end = float(new_seg[-1])

    return new_k, new_e, label_map


def plot_bands(
    bands_path: str | Path,
    kpt: KPoints | str | Path | None = None,
    *,
    e_range: tuple[float, float] | None = None,
    color: str = "#1f77b4",
    alpha: float = 0.8,
    linewidth: float = 1.2,
    title: str | None = None,
    show_fermi: bool = True,
    save: str | Path | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot a band structure from ABACUS output.

    Args:
        bands_path: Path to BANDS_*.dat file.
        kpt: KPoints object or KPT file path for symmetry labels.
        e_range: (emin, emax) energy range in eV relative to E_Fermi.
        color: Band line color.
        alpha: Band line transparency.
        linewidth: Band line width.
        title: Plot title.
        show_fermi: Whether to show the Fermi level as a dashed line.
        save: Path to save the figure (PNG/PDF).
        ax: Existing matplotlib Axes to plot on.

    Returns:
        matplotlib Axes object.
    """
    from abacuscopilot.plotting.style import load_style_from_config

    load_style_from_config()

    if ax is None:
        _, ax = plt.subplots()

    k_dists, energies, e_fermi, labels, label_positions = read_bands_with_kpt(
        bands_path, kpt
    )

    if energies.size == 0:
        ax.text(0.5, 0.5, "No band data found", transform=ax.transAxes, ha="center")
        return ax

    nbands = energies.shape[0]

    # Shift energies to E_Fermi = 0
    energies_shifted = energies - e_fermi

    # Break lines at path discontinuities (seekpath Z|X / R|M jumps),
    # then compress the x-axis so the empty gap collapses to a small visual break.
    k_dists_broken, energies_broken = _insert_path_breaks(k_dists, energies_shifted)
    k_dists_plot, energies_plot, label_map = _compress_path_breaks(k_dists_broken, energies_broken)

    # Remap label positions to compressed x-axis
    label_positions_plot = [label_map.get(pos, pos) for pos in label_positions]
    nbands_plot = energies_plot.shape[0]

    # Plot each band
    for ib in range(nbands_plot):
        ax.plot(k_dists_plot, energies_plot[ib], color=color, alpha=alpha,
                linewidth=linewidth)

    # Fermi level
    if show_fermi:
        ax.axhline(y=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

    # High-symmetry labels — replace GAMMA with Γ
    if labels and label_positions_plot:
        display_labels = [lab.replace("GAMMA", "Γ") for lab in labels]
        ax.set_xticks(label_positions_plot)
        ax.set_xticklabels(display_labels)
        for pos in label_positions_plot:
            ax.axvline(x=pos, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)
        # Trim x-axis to data range (no empty space on sides)
        ax.set_xlim(k_dists_plot[0], k_dists_plot[-1])
    else:
        ax.set_xticks([])

    # Full frame (all four spines)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    # Ticks outward
    ax.tick_params(axis="both", direction="out")

    # Axis labels
    ax.set_ylabel("E − E$_F$ (eV)")

    # Energy range
    if e_range:
        ax.set_ylim(*e_range)
    else:
        ymin = energies_shifted.min()
        ymax = energies_shifted.max()
        margin = 1.0
        ax.set_ylim(ymin - margin, ymax + margin)

    if title:
        ax.set_title(title)

    if save:
        ax.figure.tight_layout(pad=1.2)
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")

    return ax


def plot_fatbands(
    bands_path: str | Path,
    proj_path: str | Path,
    kpt: KPoints | str | Path | None = None,
    *,
    species_colors: dict[str, str] | None = None,
    e_range: tuple[float, float] | None = None,
    linewidth: float = 2.5,
    title: str | None = None,
    save: str | Path | None = None,
) -> plt.Figure:
    """Plot a fat-band (projected band) structure.

    Args:
        bands_path: Path to BANDS_*.dat file.
        proj_path: Path to projected band data (ProjBands_*.dat or similar).
        kpt: KPoints object or KPT file path.
        species_colors: Mapping from species to color.
        e_range: Energy range (emin, emax) relative to E_Fermi.
        linewidth: Maximum linewidth for fat bands.
        title: Plot title.
        save: Path to save the figure.

    Returns:
        matplotlib Figure.
    """
    from abacuscopilot.plotting.style import load_style_from_config

    load_style_from_config()

    k_dists, energies, e_fermi, labels, label_positions = read_bands_with_kpt(
        bands_path, kpt
    )

    if energies.size == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No band data found", transform=ax.transAxes, ha="center")
        return fig

    nbands, nkpts = energies.shape
    energies_shifted = energies - e_fermi

    # Try to read projected weights
    proj_data = _read_proj_bands(proj_path)

    if species_colors is None:
        species_colors = {}

    if proj_data is not None:
        # proj_data: dict[species -> (nbands, nkpts) weights]
        fig, ax = plt.subplots()

        # Break lines at path discontinuities and compress x-axis gaps
        k_dists_broken, energies_broken = _insert_path_breaks(k_dists, energies_shifted)
        k_dists_plot, energies_plot, label_map = _compress_path_breaks(k_dists_broken, energies_broken)
        label_positions_plot = [label_map.get(pos, pos) for pos in label_positions]
        nbands_plot = energies_plot.shape[0]

        # Identify which original k-point indices survive compression
        # (connector kpts are removed; the rest are kept in order)
        diffs_orig = np.diff(k_dists)
        keep_mask = np.ones(nkpts, dtype=bool)
        for brk in np.where(diffs_orig < 1e-10)[0]:
            keep_mask[brk + 1] = False  # remove the duplicate connector
        kept_indices = np.where(keep_mask)[0]  # maps compressed_idx → original_idx

        # Draw thin gray reference lines for all bands
        for ib in range(nbands_plot):
            ax.plot(k_dists_plot, energies_plot[ib], color="gray", linewidth=0.3,
                    alpha=0.4, zorder=1)

        # Compress scatter x-positions and weights by removing connector kpts.
        # k_dists_plot has the compressed positions (same length as keep_mask.sum()).
        energies_compressed = energies_shifted[:, keep_mask]
        for species, weights in proj_data.items():
            color = species_colors.get(species)
            w_compressed = weights[:, keep_mask]  # (nbands, n_kept)
            for ib in range(nbands):
                if np.any(w_compressed[ib] > 0.01):
                    widths = np.maximum(w_compressed[ib] * linewidth, 0.0)
                    ax.scatter(k_dists_plot, energies_compressed[ib], s=widths * 10,
                              c=color if color else None, alpha=0.6, linewidths=0,
                              label=species if ib == 0 else "", zorder=2)

        # Deduplicate legend
        handles, labels_ = ax.get_legend_handles_labels()
        by_label = dict(zip(labels_, handles))
        if by_label:
            ax.legend(by_label.values(), by_label.keys(), loc="upper right")
    else:
        # Fallback: simple line plot with wider lines
        fig, ax = plt.subplots()
        ax.plot(k_dists, energies_shifted.T, linewidth=linewidth, alpha=0.7)
        ax = fig.axes[0]

    # Fermi level
    ax.axhline(y=0.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

    # Symmetry labels — replace GAMMA with Γ
    if labels and label_positions_plot:
        display_labels = [lab.replace("GAMMA", "Γ") for lab in labels]
        ax.set_xticks(label_positions_plot)
        ax.set_xticklabels(display_labels)
        for pos in label_positions_plot:
            ax.axvline(x=pos, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)
        ax.set_xlim(k_dists_plot[0], k_dists_plot[-1])
    else:
        ax.set_xticks([])

    # Full frame and ticks
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")

    ax.set_ylabel("E − E$_F$ (eV)")

    if e_range:
        ax.set_ylim(*e_range)
    else:
        ymin = energies_shifted.min()
        ymax = energies_shifted.max()
        margin = 1.0
        ax.set_ylim(ymin - margin, ymax + margin)

    if title:
        ax.set_title(title)

    if save:
        fig.tight_layout(pad=1.2)
        fig.savefig(save, dpi=300, bbox_inches="tight")

    return fig


def _find_e_fermi(bands_path: Path) -> float:
    """Try to extract the Fermi energy from SCF/NSCF logs near the BANDS file.

    Uses the shared dos.py parser which handles all ABACUS log formats
    (PW / LCAO GPU / LCAO CPU).  Falls back to generic pattern matching.
    Returns 0.0 if not found.
    """
    import re

    # Prefer the unified parser from dos.py (handles all known ABACUS formats)
    from abacuscopilot.plotting.dos import _find_fermi_from_log
    e_fermi = _find_fermi_from_log(bands_path)
    if e_fermi != 0.0:
        return e_fermi

    # Fallback: generic search in sibling log files
    out_dir = bands_path.parent  # OUT.ABACUS/
    for log_name in ("running_scf.log", "running_nscf.log", "running_relax.log"):
        log_path = out_dir / log_name
        if log_path.exists():
            try:
                content = log_path.read_text(errors="ignore")
                # ABACUS v3.10: "E_Fermi  0.193797 Ry  2.636743 eV" — take the LAST number (eV)
                m_line = re.search(
                    r"(?:E_Fermi|E-fermi|EFERMI|Fermi\s*(?:energy|level)).*",
                    content, re.IGNORECASE,
                )
                if m_line:
                    nums = re.findall(r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)", m_line.group())
                    if nums:
                        return float(nums[-1])  # last number on the line = eV
            except Exception:
                pass
    return 0.0


def _read_proj_bands(filepath: str | Path) -> dict[str, np.ndarray] | None:
    """Parse an ABACUS projected band file (old .dat or new XML PBANDS).

    Returns:
        Dict mapping species label -> (nbands, nkpts) weight array, or None.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    # Try XML format first (ABACUS v3.7+ PBANDS_1, PBANDS_2)
    try:
        return _read_proj_bands_xml(filepath)
    except Exception:
        pass

    # Fall back to old text format (ProjBands*.dat)
    try:
        with open(filepath) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]

        if not lines:
            return None

        header = lines[0].split()
        if len(header) < 2:
            return None

        nkpts = int(header[0])
        nbands = int(header[1])
        nspecies = int(header[2]) if len(header) >= 3 else 1

        species_names = []
        if len(header) > 3:
            species_names = header[3:3 + nspecies]

        all_data = []
        for line in lines[1:]:
            all_data.extend(float(x) for x in line.split())

        total = nkpts * nbands * nspecies
        if len(all_data) < total:
            return None

        data = np.array(all_data[:total]).reshape(nspecies, nkpts, nbands)
        data = data.transpose(0, 2, 1)  # (nspecies, nbands, nkpts)

        result = {}
        for i in range(nspecies):
            name = species_names[i] if i < len(species_names) else f"sp{i+1}"
            result[name] = data[i]

        return result

    except Exception:
        return None


def _read_proj_bands_xml(filepath: Path) -> dict[str, np.ndarray]:
    """Parse ABACUS XML PBANDS file (v3.7+).

    Each <orbital> element has species, l, m info and <data> with (nkpts, nbands)
    weights.  Returns dict[species -> (nbands, nkpts)].
    """
    import xml.etree.ElementTree as ET
    from io import StringIO

    tree = ET.parse(str(filepath))
    root = tree.getroot()

    # Sum weights by species
    species_weights: dict[str, list] = {}
    for orb in root.findall("orbital"):
        sp = orb.get("species", "unknown")
        data_text = orb.find("data")
        if data_text is None or data_text.text is None:
            continue
        arr = np.loadtxt(StringIO(data_text.text))
        # arr shape: (nkpts, nbands)
        if sp not in species_weights:
            species_weights[sp] = []
        species_weights[sp].append(arr)

    # Sum all orbitals of same species, then average
    result = {}
    for sp, arrays in species_weights.items():
        stacked = np.stack(arrays)             # (n_orbitals, nkpts, nbands)
        summed = np.sum(stacked, axis=0)       # (nkpts, nbands)
        result[sp] = summed.T                  # (nbands, nkpts)

    return result
