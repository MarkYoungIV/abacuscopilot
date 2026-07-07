"""Matplotlib style configuration for abacuscopilot.

Provides a clean, publication-quality default style and helper functions
for applying it to matplotlib figures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib as mpl

# Default abacuscopilot style
ABACUSKIT_STYLE: dict[str, Any] = {
    "figure.figsize": (8, 6),
    "figure.dpi": 300,
    "figure.facecolor": "white",
    "figure.edgecolor": "white",
    "font.size": 12,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "axes.linewidth": 1.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.format": "png",
}


def apply_style(style_name: str = "abacuscopilot") -> None:
    """Apply a named style to matplotlib.

    Args:
        style_name: 'abacuscopilot' (default), or a path to a custom style file.
    """
    if style_name == "abacuscopilot":
        _apply_dict_style(ABACUSKIT_STYLE)
    else:
        # Try loading from ~/.abacuscopilot/ or as a direct path
        custom_path = Path(style_name)
        if not custom_path.exists():
            custom_path = Path.home() / ".abacuscopilot" / f"{style_name}.mplstyle"
        if custom_path.exists():
            mpl.style.use(str(custom_path))
        else:
            _apply_dict_style(ABACUSKIT_STYLE)


def _apply_dict_style(style: dict[str, Any]) -> None:
    """Apply style from a dictionary."""
    for key, value in style.items():
        try:
            mpl.rcParams[key] = value
        except KeyError:
            pass  # Skip unknown rcParams


def load_style_from_config() -> None:
    """Load plotting preferences from the user's abacuscopilot config."""
    try:
        from abacuscopilot.config import load_config
        config = load_config()
        plot_config = config.get("plotting", {})

        if plot_config.get("style", "abacuscopilot") != "abacuscopilot":
            apply_style(plot_config["style"])
        else:
            apply_style("abacuscopilot")

        # Override with specific config values
        for key in ("dpi", "font_size", "figure_format"):
            if key in plot_config:
                rc_key = {
                    "dpi": "figure.dpi",
                    "font_size": "font.size",
                    "figure_format": "savefig.format",
                }.get(key)
                if rc_key:
                    mpl.rcParams[rc_key] = plot_config[key]

        if "figure_size" in plot_config:
            mpl.rcParams["figure.figsize"] = plot_config["figure_size"]

    except Exception:
        apply_style("abacuscopilot")  # Fallback on error


# --- Color schemes for band structure plots ---

COLOR_SCHEMES: dict[str, dict] = {
    "deep": {
        "name": "深色稳重",
        "colors": ["#1a1a2e", "#16213e", "#0f3460", "#533483", "#e94560"],
        "type": "discrete",
    },
    "nature": {
        "name": "Nature 风格",
        "colors": ["#2d5f8b", "#3d7ea6", "#6baed6", "#fd8d3c", "#31a354"],
        "type": "discrete",
    },
    "flat": {
        "name": "现代扁平",
        "colors": ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
                   "#937860", "#8ca5c8", "#ccb974"],
        "type": "discrete",
    },
    "dark": {
        "name": "暗色背景",
        "colors": ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"],
        "type": "discrete",
    },
    "viridis": {
        "name": "Viridis",
        "colors": None,
        "type": "cmap",
        "cmap": "viridis",
    },
    "plasma": {
        "name": "Plasma",
        "colors": None,
        "type": "cmap",
        "cmap": "plasma",
    },
    "warm": {
        "name": "Warm",
        "colors": None,
        "type": "cmap",
        "cmap": "Warm",
    },
    "coolwarm": {
        "name": "Coolwarm",
        "colors": None,
        "type": "cmap",
        "cmap": "coolwarm",
    },
    "tab10": {
        "name": "Tab10 (色盲友好)",
        "colors": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
        "type": "discrete",
    },
    "set2": {
        "name": "Set2 (柔和)",
        "colors": ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
                   "#ffd92f", "#e5c494", "#b3b3b3"],
        "type": "discrete",
    },
}

COLOR_SCHEME_ORDER = ["deep", "nature", "flat", "dark",
                      "viridis", "plasma", "warm", "coolwarm",
                      "tab10", "set2"]


def get_colors(scheme_key: str, n: int) -> list:
    """Return *n* colors sampled from the named scheme.

    For discrete palettes the colors cycle.  For continuous colormaps the
    first *n* evenly-spaced values in [0,1] are sampled.
    """
    scheme = COLOR_SCHEMES.get(scheme_key)
    if scheme is None:
        scheme = COLOR_SCHEMES["deep"]

    if scheme["type"] == "cmap":
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(scheme["cmap"])
        return [cmap(i / max(n - 1, 1)) for i in range(n)]

    colors = scheme["colors"]
    return [colors[i % len(colors)] for i in range(n)]


def get_base_color(scheme_key: str) -> str:
    """Return the primary line colour for a scheme."""
    return get_colors(scheme_key, 1)[0]


def configure_latex(use_latex: bool = False) -> None:
    """Configure matplotlib to use LaTeX for text rendering.

    Args:
        use_latex: If True, enable LaTeX rendering.
    """
    if use_latex:
        mpl.rcParams.update({
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "font.family": "serif",
        })
