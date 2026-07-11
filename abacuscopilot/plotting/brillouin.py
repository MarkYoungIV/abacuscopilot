"""3D Brillouin-zone visualization for high-symmetry k-paths.

Draws the first Brillouin zone (via seekpath's BZ geometry), the
high-symmetry points, and the k-path connecting them, as a static PNG.
Used by tasks 302 (band path) and 304 (phonon q-path).
"""

from __future__ import annotations

import numpy as np


def plot_brillouin_zone(
    recip_lattice,
    segments: list,
    out_png: str = "brillouin_zone.png",
    title: str | None = None,
) -> str | None:
    """Draw the Brillouin zone with the actual k-path segments on it.

    Args:
        recip_lattice: 3x3 reciprocal primitive lattice; rows are b1, b2, b3.
        segments: list of segment dicts, each with keys
            "start" / "end" (fractional coords in reciprocal basis) and
            "label" / "end_label" (point names). This is exactly what gets
            written to the KPT file (kpts.line_path), so the plot matches the
            KPT 1:1 — important for lattices (e.g. FCC) where seekpath breaks
            the path at equivalent points like U/K.
        out_png: output PNG path.
        title: optional figure title.

    Returns the output path on success, or None if plotting failed (plotting is
    an optional add-on and must never break the KPT/band.conf that was written).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        from seekpath.brillouinzone.brillouinzone import BZ

        from abacuscopilot.plotting.style import load_style_from_config
        load_style_from_config()

        b = np.asarray(recip_lattice, dtype=float)   # rows b1,b2,b3
        b1, b2, b3 = b[0], b[1], b[2]

        # --- BZ polyhedron faces (cartesian vertices) ---
        zone = BZ(b1, b2, b3)
        faces = zone.faces

        # fractional (reciprocal-basis) -> cartesian
        def to_cart(frac):
            f = np.asarray(frac, dtype=float)
            return f[0] * b1 + f[1] * b2 + f[2] * b3

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")

        # BZ shell: translucent faces + solid edges
        poly = Poly3DCollection(
            faces, alpha=0.10, facecolor="#4c72b0", edgecolor="#333333", linewidths=0.8
        )
        ax.add_collection3d(poly)

        # --- draw path segments exactly as written to KPT ---
        # Collect the labelled points actually visited (dedupe by cartesian pos).
        labelled = {}   # display_label -> cartesian xyz
        for seg in segments:
            p1 = to_cart(seg["start"])
            p2 = to_cart(seg["end"])
            ax.plot(
                [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                color="#d62728", linewidth=2.0, zorder=4,
            )
            for lbl, xyz in ((seg.get("label", ""), p1), (seg.get("end_label", ""), p2)):
                if lbl:
                    disp = "Γ" if lbl.upper() == "GAMMA" else lbl
                    labelled[disp] = xyz

        # high-symmetry points + labels (only those on the path)
        for disp, xyz in labelled.items():
            ax.scatter(*xyz, color="#d62728", s=35, depthshade=False, zorder=5)
            ax.text(xyz[0], xyz[1], xyz[2], f"  {disp}", fontsize=11, zorder=6)

        # cosmetics: equal aspect, clean background, sensible view
        _set_equal_3d(ax, faces)
        ax.set_xlabel("$k_x$")
        ax.set_ylabel("$k_y$")
        ax.set_zlabel("$k_z$")
        ax.view_init(elev=18, azim=30)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        # de-clutter: hide panes/grid so the BZ stands out
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_alpha(0.0)
        if title:
            ax.set_title(title)

        fig.tight_layout(pad=1.2)
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return out_png
    except Exception as exc:  # never propagate — plotting is optional
        try:
            from abacuscopilot.console_utils import _get_console
            _get_console().print(
                f"  [yellow]Brillouin-zone plot skipped: {exc}[/yellow]"
            )
        except Exception:
            pass
        return None


def _set_equal_3d(ax, faces) -> None:
    """Give the 3D axes equal scale so the BZ isn't distorted."""
    pts = np.array([v for face in faces for v in face], dtype=float)
    if pts.size == 0:
        return
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = (mins + maxs) / 2.0
    half = (maxs - mins).max() / 2.0
    if half <= 0:
        return
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
