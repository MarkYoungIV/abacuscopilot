"""Reader and writer for ABACUS KPT files.

The KPT file controls Brillouin zone sampling for ABACUS calculations.
Supports three modes:
- Automatic Monkhorst-Pack mesh generation (mode=0)
- Explicit k-point list (mode=Nkpoints)
- Line mode for band structure calculations
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.core.exceptions import FileFormatError, FileNotFoundError_
from abacuscopilot.core.models import KPoints, Lattice


def read_kpt(filepath: str | Path) -> KPoints:
    """Read and parse an ABACUS KPT file.

    Args:
        filepath: Path to the KPT file.

    Returns:
        KPoints object with parsed data.

    Raises:
        FileNotFoundError_: If file doesn't exist.
        FileFormatError: If the file format is invalid.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError_(str(filepath), "KPT file not found")

    with open(filepath) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    if not lines:
        raise FileFormatError(str(filepath), "KPT file is empty")

    # First line should be K_POINTS, KPOINTS, or K
    header = lines[0].upper()
    if not (header.startswith("K_POINTS") or header.startswith("KPOINTS") or header == "K"):
        raise FileFormatError(
            str(filepath),
            f"Expected 'K_POINTS' on first line, got '{lines[0]}'"
        )

    # Second line determines mode
    second = int(lines[1])

    if second == 0:
        # Auto MP mesh
        kpts = KPoints()
        kpts.mode = "gamma"

        # Third line: Gamma or MP
        if len(lines) > 2:
            center = lines[2].strip().lower()
            kpts.gamma_centered = center in ("gamma", "g", "1")

        # Fourth line: N1 N2 N3 S1 S2 S3
        if len(lines) > 3:
            parts = lines[3].split()
            n_vals = [int(x) for x in parts[:3]]
            s_vals = [float(x) for x in parts[3:6]] if len(parts) >= 6 else [0.0, 0.0, 0.0]
            kpts.grid = tuple(n_vals)
            kpts.shift = tuple(s_vals)
        else:
            kpts.grid = (1, 1, 1)

        return kpts

    elif second > 0:
        # Check third line for mode
        if len(lines) < 3:
            raise FileFormatError(str(filepath), "Missing coordinate type line")

        third = lines[2].strip()
        third_upper = third.upper()

        if third_upper.startswith("LINE"):
            # Line mode for band structure
            kpts = KPoints()
            kpts.mode = "line_cartesian" if "CARTESIAN" in third_upper else "line"

            n_endpoints = second
            # Each line after the third: kx ky kz npoints [label]
            # The very last entry is an endpoint only (npoints typically = 1)
            endpoint_lines = lines[3:3 + n_endpoints]
            for i, line in enumerate(endpoint_lines):
                parts = line.split()
                if len(parts) < 4:
                    continue

                kx, ky, kz = float(parts[0]), float(parts[1]), float(parts[2])
                npts = int(parts[3])
                # Label: handle both "GAMMA" and "# GAMMA" formats
                label = ""
                if len(parts) >= 5:
                    label = parts[4]
                    if label == "#" and len(parts) >= 6:
                        label = parts[5]

                # Set the end of the previous segment to this point
                if kpts.line_path:
                    kpts.line_path[-1]["end"] = (kx, ky, kz)
                    kpts.line_path[-1]["end_label"] = label

                # Start a new segment unless this is the last entry
                is_last = (i == len(endpoint_lines) - 1)
                if not is_last and npts > 0:
                    kpts.line_path.append({
                        "start": (kx, ky, kz),
                        "end": (0.0, 0.0, 0.0),  # placeholder
                        "npoints": npts,
                        "label": label,
                        "end_label": "",
                    })

            return kpts

        else:
            # Explicit k-point list
            kpts = KPoints()
            kpts.mode = "direct"

            nkpts = second
            for i in range(3, min(3 + nkpts, len(lines))):
                parts = lines[i].split()
                if len(parts) >= 4:
                    kpts.explicit_kpoints.append(
                        (float(parts[0]), float(parts[1]),
                         float(parts[2]), float(parts[3]))
                    )
                elif len(parts) >= 3:
                    kpts.explicit_kpoints.append(
                        (float(parts[0]), float(parts[1]),
                         float(parts[2]), 1.0)
                    )

            return kpts

    raise FileFormatError(str(filepath), f"Unrecognized KPT format (second line = {second})")


def write_kpt(kpts: KPoints, filepath: str | Path = "KPT") -> None:
    """Write a KPoints object to an ABACUS KPT file.

    Args:
        kpts: KPoints object to write.
        filepath: Output file path.
    """
    content = kpts.to_string()
    with open(filepath, "w") as f:
        f.write(content)


def auto_mp_kpts(
    lattice: Lattice,
    kspacing: float = 0.20,
    gamma_centered: bool = True,
) -> KPoints:
    """Generate auto Monkhorst-Pack k-point mesh from k-spacing.

    The grid dimensions are computed from the reciprocal lattice vectors
    using the formula: N_i = max(1, ceil(|b_i| / kspacing)).

    kspacing is in units of 2π/Å (ABACUS v3.x convention).
    This is equivalent to the VASP rule-of-thumb a·N ≈ 30:
        kspacing ≈ 2π / 30 ≈ 0.21

    Args:
        lattice: Lattice object.
        kspacing: Maximum k-point spacing in 2π/Å (default 0.20).
        gamma_centered: Whether to use Gamma-centered mesh.

    Returns:
        KPoints object with auto MP mesh.
    """
    # Reciprocal lattice in 2π/Bohr
    recp = lattice.reciprocal_cell  # rows = b1, b2, b3
    recp_lengths = np.linalg.norm(recp, axis=1)  # |b_i| in 2π/Bohr

    # Convert |b_i| from 2π/Bohr to 2π/Å to match kspacing units
    #   |b| (2π/Å) = 2π / a_Å = 2π / (a_Bohr * BOHR_TO_ANGSTROM)
    #   |b| (2π/Å) = |b| (2π/Bohr) * ANGSTROM_TO_BOHR
    from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
    recp_lengths_ang = recp_lengths * ANGSTROM_TO_BOHR  # → 2π/Å

    # N_i = max(1, ceil(|b_i| / kspacing))
    grid = tuple(
        max(1, int(np.ceil(length / kspacing)))
        for length in recp_lengths_ang
    )

    return KPoints(
        mode="gamma" if gamma_centered else "mp",
        grid=grid,
        shift=(0.0, 0.0, 0.0),
        gamma_centered=gamma_centered,
    )


def line_mode_kpts_from_path(
    path: list[tuple[list[float], list[float], int]],
    labels: list[str] | None = None,
    mode: str = "line",
) -> KPoints:
    """Generate line-mode KPT for band structure from a custom path.

    Args:
        path: List of (start_xyz, end_xyz, npoints) tuples.
        labels: Optional list of high-symmetry point labels.
        mode: 'line' (fractional) or 'line_cartesian'.

    Returns:
        KPoints object configured for band structure calculation.
    """
    kpts = KPoints(mode=mode)

    for i, (start, end, npts) in enumerate(path):
        seg = {
            "start": tuple(start),
            "end": tuple(end),
            "npoints": npts,
            "label": labels[i] if labels and i < len(labels) else "",
            "end_label": labels[i + 1] if labels and i + 1 < len(labels) else "",
        }
        kpts.line_path.append(seg)

    return kpts
