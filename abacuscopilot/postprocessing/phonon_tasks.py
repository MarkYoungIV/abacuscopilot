"""Phonon post-processing tasks for ABACUS + Phonopy workflow.

Task 1501: Extract forces, compute phonon bands, and plot.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from abacuscopilot.config import load_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.tasks import task


def _find_phonopy() -> str | None:
    """Return path to the phonopy executable, or None."""
    import shutil as _shutil
    path = _shutil.which("phonopy")
    if path:
        return path
    try:
        subprocess.run(
            [sys.executable, "-m", "phonopy", "--version"],
            capture_output=True, timeout=10,
        )
        return sys.executable
    except Exception:
        return None


def _run_phonopy(args: list[str], cwd: str | Path = ".", timeout: int = 120) -> subprocess.CompletedProcess:
    """Run phonopy as a subprocess (for -p, band plotting, etc.)."""
    import shutil as _shutil, os
    for candidate in [
        _shutil.which("phonopy"),
        Path(sys.executable).parent / "phonopy",
        Path(sys.prefix) / "bin" / "phonopy",
    ]:
        if candidate and Path(str(candidate)).exists():
            cmd = [str(candidate)] + args
            break
    else:
        cmd = [sys.executable, "-m", "phonopy"] + args
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", str(min(os.cpu_count() or 4, 16)))
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=timeout, env=env,
    )


def _run_phonopy_init(args: list[str], cwd: str | Path = ".", timeout: int = 120) -> subprocess.CompletedProcess:
    """Run phonopy-init (v4+) for setup operations: -d, -f."""
    import shutil as _shutil
    for candidate in [
        _shutil.which("phonopy-init"),
        Path(sys.executable).parent / "phonopy-init",
        Path(sys.prefix) / "bin" / "phonopy-init",
    ]:
        if candidate and Path(str(candidate)).exists():
            cmd = [str(candidate)] + args
            break
    else:
        cmd = [sys.executable, "-m", "phonopy", "--init"] + args
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=timeout,
    )


@task(1501, category="Lattice Dynamics", name="Phonon Analysis",
      description="Extract FORCE_SETS, compute phonon bands, and plot dispersion",
      cli_args=[
          {"name": "--dim", "type": str, "default": "2 2 2",
           "help": "Supercell dimensions (must match setup)"},
          {"name": "--mesh", "type": str, "default": "8 8 8",
           "help": "q-point mesh for DOS (e.g. '8 8 8')"},
          {"name": "--fmin", "type": float, "default": None,
           "help": "Minimum frequency for plot y-axis (auto if unset)"},
          {"name": "--fmax", "type": float, "default": None,
           "help": "Maximum frequency for plot y-axis (auto if unset)"},
          {"name": "--no-plot", "action": "store_true", "default": False,
           "help": "Skip plotting"},
      ])
def task_phonon_analysis(args: list[str] | None = None, interactive: bool = True,
                         parsed_args=None) -> None:
    """Post-process phonon calculation: extract forces, compute bands, plot.

    Requires phonopy ≥ 2.19.1 and completed ABACUS runs in all disp-*/ directories.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Phonon Analysis (Phonopy) ===[/bold cyan]")
    console.print("[dim]Force extraction → band structure → plot[/dim]")
    console.print()

    # --- Parameters (auto-read from phonopy_setup.json if available) ---
    setup_info: dict = {}
    setup_file = Path("phonopy_setup.json")
    if setup_file.exists():
        try:
            import json as _json
            setup_info = _json.loads(setup_file.read_text())
            console.print(f"  [dim]Read setup from phonopy_setup.json[/dim]")
        except Exception:
            pass

    if setup_info.get("dim"):
        dim_s = setup_info["dim"]
    elif interactive:
        dim_s = _prompt(console, "Supercell dimensions (must match setup)", "2 2 2")
    elif parsed_args:
        dim_s = parsed_args.dim
    else:
        dim_s = "2 2 2"

    # Y-axis range
    fmin_plot = None
    fmax_plot = None
    if interactive:
        mesh_s = _prompt(console, "q-point mesh for DOS/thermal", "8 8 8")
    elif parsed_args:
        mesh_s = parsed_args.mesh
        if getattr(parsed_args, "fmin", None) is not None:
            fmin_plot = float(parsed_args.fmin)
        if getattr(parsed_args, "fmax", None) is not None:
            fmax_plot = float(parsed_args.fmax)
    else:
        mesh_s = "8 8 8"

    # --- Step 1: Extract FORCE_SETS ---
    console.print("[bold]Step 1: Extract FORCE_SETS[/bold]")

    # Collect running_scf.log from all disp-* directories
    log_paths: list[Path] = []
    disp_dirs = sorted(Path(".").glob("disp-*"))
    if not disp_dirs:
        console.print("[red]No disp-*/ directories found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 112' first.[/dim]")
        return

    for dd in disp_dirs:
        candidates = list(dd.glob("OUT*/running_scf.log")) + list(dd.glob("OUT*/running_relax.log"))
        if candidates:
            log_paths.append(candidates[0])
        else:
            console.print(f"  [yellow]![/yellow] {dd.name}: no OUT*/running_scf.log")

    if not log_paths:
        console.print("[red]No ABACUS output logs found in disp-*/ directories.[/red]")
        return

    console.print(f"  Logs found: {len(log_paths)}/{len(disp_dirs)}")

    # Run phonopy-init -f (v4+ moved -f from phonopy to phonopy-init)
    force_args = ["-f"] + [str(p) for p in log_paths]
    result = _run_phonopy_init(force_args, timeout=30)
    if result.returncode != 0:
        console.print("[red]phonopy-init -f failed:[/red]")
        console.print(result.stderr)
        return
    if result.stderr:
        console.print(f"  [dim]{result.stderr.strip()}[/dim]")

    # Check FORCE_SETS was created
    if not Path("FORCE_SETS").exists():
        console.print("[red]FORCE_SETS was not created.[/red]")
        return
    console.print("  [green]✓ FORCE_SETS[/green]")

    # --- Step 2: Generate band.conf ---
    console.print()
    console.print("[bold]Step 2: Generate band.conf[/bold]")

    # Read species from STRU
    species_list: list[str] = []
    try:
        from abacuscopilot.io.stru_file import read_stru
        s = read_stru("STRU")
        species_list = list(s.species_order)
    except Exception:
        pass

    # Try to get k-path from seekpath
    kpath_str = ""
    try:
        from seekpath import get_path
        from abacuscopilot.io.stru_file import read_stru
        from ase.data import atomic_numbers
        s = read_stru("STRU")
        # seekpath expects fractional coords; convert to plain lists
        positions = np.array([a.position for a in s.atoms])
        if s.coordinate_type != "Direct":
            cell_inv = np.linalg.inv(s.lattice.cell_angstrom)
            positions = positions @ cell_inv
        numbers = [atomic_numbers.get(sp, 0) for sp in [a.species for a in s.atoms]]
        sp = get_path((s.lattice.cell_angstrom.tolist(), positions.tolist(), numbers),
                      with_time_reversal=True)
        point_coords = sp["point_coords"]  # {label: [k1, k2, k3]}
        band_parts: list[str] = []
        band_labels_parts: list[str] = []
        for seg in sp["path"]:  # [('GAMMA', 'X'), ('X', 'U'), ...]
            start_label, end_label = seg[0], seg[1]
            for label in (start_label, end_label):
                c = point_coords[label]
                band_parts.append(f"{c[0]:.6f} {c[1]:.6f} {c[2]:.6f}")
                band_labels_parts.append(label)
        kpath_str = "  ".join(band_parts)
        klabels_str = " ".join(band_labels_parts)
        console.print("  [dim]k-path from seekpath[/dim]")
    except Exception:
        kpath_str = "0 0 0  0.5 0.5 0  0.5 0.5 0.5  0 0 0"
        klabels_str = "GM X X UK GM GM L L W W X"
        console.print("  [yellow]![/yellow] cannot determine k-path, using default")

    # PRIMITIVE_AXES — default identity for now
    prim_axes = "1 0 0  0 1 0  0 0 1"

    band_conf = f"""DIM = {dim_s}
MESH = {mesh_s}
PRIMITIVE_AXES = {prim_axes}
"""
    with open("band.conf", "w") as f:
        f.write(band_conf)
    console.print("  [green]✓ band.conf[/green]")

    # --- Step 2.5: Compute FORCE_CONSTANTS (expensive, one-time) ---
    fc_file = Path("FORCE_CONSTANTS")
    if not fc_file.exists():
        console.print()
        console.print("[bold]Step 2.5: Compute FORCE_CONSTANTS (one-time)[/bold]")
        console.print("  [dim]This may take minutes for large systems. Subsequent runs reuse this file.[/dim]")
        _run_phonopy_with_progress(console, ["--writefc", "band.conf"], timeout=1200)
        if not fc_file.exists():
            console.print("[red]FORCE_CONSTANTS was not created.[/red]")
            return
        console.print("  [green]✓ FORCE_CONSTANTS[/green]")
    else:
        console.print(f"  [dim]FORCE_CONSTANTS exists ({fc_file.stat().st_size // 1024} KB), reusing[/dim]")

    # --- Step 3: Compute phonon bands (fast with --readfc) ---
    console.print()
    console.print("[bold]Step 3: Compute phonon bands[/bold]")
    console.print("  Running: [dim]phonopy --readfc --band ... band.conf[/dim]")

    kpoints = kpath_str.split()
    klabels = klabels_str.split()
    band_args = (
        ["--band"] + kpoints +
        ["--band-labels"] + klabels +
        ["--band-points", "101", "--band-connection", "band.conf"]
    )

    result = _run_phonopy(band_args, timeout=600)
    if result.returncode != 0:
        console.print("[red]phonopy band calculation failed:[/red]")
        console.print(result.stderr)
        return
    if result.stderr:
        for line in result.stderr.strip().split("\n")[:5]:
            console.print(f"  [dim]{line}[/dim]")

    if not Path("band.yaml").exists():
        console.print("[red]band.yaml was not created.[/red]")
        return
    console.print("  [green]✓ band.yaml[/green]")

    # --- Scan band.yaml for frequency range (light: grep, don't parse YAML) ---
    f_auto_min = 0.0
    f_auto_max = 20.0
    try:
        import re
        freq_vals = []
        with open("band.yaml") as f:
            for line in f:
                m = re.search(r"frequency:\s*([-\d.]+)", line)
                if m:
                    freq_vals.append(float(m.group(1)))
        if freq_vals:
            f_auto_min = min(freq_vals) * 1.05
            f_auto_max = max(freq_vals) * 1.05
    except Exception:
        pass

    # --- Step 4: Plot ---
    do_plot = True
    if parsed_args and parsed_args.no_plot:
        do_plot = False
    elif interactive:
        answer = _prompt_choice(console, "Generate phonon band plot?", ["Yes", "No"], "Yes")
        do_plot = "Yes" in answer
        if do_plot:
            fmin_s = _prompt(console,
                f"Frequency min (auto: {f_auto_min:.1f} THz, enter for auto)", "")
            fmax_s = _prompt(console,
                f"Frequency max (auto: {f_auto_max:.1f} THz, enter for auto)", "")
            fmin_plot = float(fmin_s) if fmin_s.strip() else None
            fmax_plot = float(fmax_s) if fmax_s.strip() else None

    if do_plot:
        console.print()
        console.print("[bold]Step 4: Plot phonon dispersion[/bold]")
        try:
            _plot_phonon_bands(console, fmin_plot, fmax_plot)
            # Combined three-panel figure if DOS data exists
            tdos_path = Path("total_dos.dat")
            pdos_path = Path("partial_dos.dat")
            if tdos_path.exists():
                _plot_phonon_combined(console, tdos_path, pdos_path, fmin_plot, fmax_plot)
        except Exception as e:
            console.print(f"[red]Plot failed: {e}[/red]")

    console.print()


def _parse_band_yaml(band_yaml: Path):
    """Parse band.yaml by streaming (avoids loading huge YAML into memory).

    Returns (xs, freqs, tick_pos, tick_labels) where:
      xs[:]      — continuous x-axis (distance)
      freqs[:,:] — shape (n_qpts, n_bands)
      tick_pos   — list of x positions for segment boundaries
      tick_labels — corresponding high-symmetry labels
    """
    import re as _re

    # Stream-parse: extract distances, frequencies, labels, segment_nqpoint
    distances: list[float] = []
    all_freqs: list[list[float]] = []
    current_bands: list[float] = []
    in_band = False
    seg_nq: list[int] = []
    raw_labels: list[list[str]] = []

    with open(band_yaml) as f:
        for line in f:
            # Track q-point entries
            if _re.match(r"^- q-position:", line):
                if current_bands:
                    all_freqs.append(current_bands)
                    current_bands = []
                in_band = False
            elif _re.match(r"^\s+distance:", line):
                m = _re.search(r"([-\d.]+(?:[eE][+-]?\d+)?)", line)
                if m:
                    distances.append(float(m.group(1)))
            elif _re.match(r"^\s+frequency:", line):
                m = _re.search(r"([-\d.]+(?:[eE][+-]?\d+)?)", line)
                if m:
                    current_bands.append(float(m.group(1)))
            elif _re.match(r"^\s+segment_nqpoint:", line):
                for m in _re.finditer(r"\d+", line):
                    seg_nq.append(int(m.group()))
            # Parse labels section
            elif _re.match(r"^labels:", line):
                pass  # handled below via yaml for simplicity

    if current_bands:
        all_freqs.append(current_bands)

    if not all_freqs:
        raise ValueError("band.yaml has no phonon data")

    distances = np.array(distances, dtype=float)
    freqs_raw = np.array(all_freqs, dtype=float)

    # Parse labels with yaml (small section, fast)
    import yaml
    labels_text = []
    with open(band_yaml) as f:
        in_labels = False
        for line in f:
            if line.startswith("labels:"):
                in_labels = True
                continue
            if in_labels:
                if line.startswith("reciprocal_lattice:") or line.startswith("natom:"):
                    break
                labels_text.append(line)
    try:
        raw_labels = yaml.safe_load("".join(labels_text))
    except Exception:
        raw_labels = []

    if not seg_nq:
        # Estimate: each segment has same number of q-points
        if raw_labels:
            seg_nq = [len(distances) // len(raw_labels)] * len(raw_labels)
        else:
            seg_nq = [len(distances)]

    # ... rest of the function below, same as before

    # Find segment breaks where distance is not monotonically increasing
    breaks = [0]
    for i in range(1, len(distances)):
        if distances[i] < distances[i - 1] + 1e-10:
            breaks.append(i)
    breaks.append(len(distances))

    # Build continuous x-axis per segment
    xs_list, freq_list = [], []
    for k in range(len(breaks) - 1):
        s, e = breaks[k], breaks[k + 1]
        seg_dist = distances[s:e] - distances[s]
        if k > 0:
            cursor = xs_list[-1][-1]
        else:
            cursor = 0.0
        xs_list.append(seg_dist + cursor)
        freq_list.append(freqs_raw[s:e])

    xs_all = np.concatenate(xs_list)
    freq_all = np.concatenate(freq_list, axis=0)

    # Build tick positions from segment_nqpoint + labels (already parsed)
    tick_pos = []
    tick_labels_raw: list[list[str]] = []

    if raw_labels and seg_nq:
        q_idx = 0
        for i, seg in enumerate(raw_labels):
            sl, el = seg[0], seg[1]
            n_pts = seg_nq[i] if i < len(seg_nq) else 101
            tick_labels_raw.append([sl, el])
            if i == 0:
                tick_pos.append(distances[q_idx] if q_idx < len(distances) else 0.0)
            q_idx += n_pts
            if q_idx - 1 < len(distances):
                tick_pos.append(distances[q_idx - 1])
    else:
        # Fallback: no explicit labels, use segment boundaries
        for k in range(len(breaks) - 1):
            s, e = breaks[k], breaks[k + 1]
            sl = phonon[s].get("label", "").strip()
            el = phonon[e - 1].get("label", "").strip()
            tick_labels_raw.append([sl, el])
            if k == 0:
                tick_pos.append(distances[s])
            tick_pos.append(distances[e - 1])

    tick_labels: list[str] = []
    for k, (seg_xs, (s_lbl, e_lbl)) in enumerate(zip(xs_list, tick_labels_raw)):
        if k == 0:
            tick_labels.append(_fmt_label(s_lbl))
        else:
            prev = tick_labels[-1]
            cur = _fmt_label(s_lbl)
            if prev != cur and cur and prev:
                tick_labels[-1] = f"{prev}|{cur}"
            elif cur and not prev:
                tick_labels[-1] = cur
        tick_labels.append(_fmt_label(e_lbl))

    # Handle case where tick_labels might be shorter than tick_pos
    while len(tick_labels) < len(tick_pos):
        tick_labels.append("")

    return xs_all, freq_all, tick_pos, tick_labels


def _fmt_label(lbl: str) -> str:
    if not lbl:
        return ""
    u = lbl.upper()
    if u in ("GAMMA", "GM", "G"):
        return r"$\Gamma$"
    if lbl in (r"\Gamma", r"$\Gamma$"):
        return r"$\Gamma$"
    if lbl.startswith("$") and lbl.endswith("$"):
        return lbl
    return f"${lbl}$"


def _run_phonopy_with_progress(console, args, timeout=1200):
    """Run phonopy with real-time progress from stderr."""
    import shutil as _shutil
    for candidate in [
        _shutil.which("phonopy"),
        Path(sys.executable).parent / "phonopy",
        Path(sys.prefix) / "bin" / "phonopy",
    ]:
        if candidate and Path(str(candidate)).exists():
            cmd = [str(candidate)] + args
            break
    else:
        cmd = [sys.executable, "-m", "phonopy"] + args

    import time
    start = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                            stderr=subprocess.PIPE, text=True,
                            bufsize=1)
    last_line = ""
    while proc.poll() is None:
        if proc.stderr:
            line = proc.stderr.readline()
            if line:
                stripped = line.strip()
                if stripped and stripped != last_line:
                    elapsed = int(time.time() - start)
                    console.print(f"  [dim][{elapsed}s] {stripped[:120]}[/dim]")
                    last_line = stripped
        if time.time() - start > timeout:
            proc.kill()
            console.print(f"[red]Timed out after {int(time.time()-start)}s.[/red]")
            return

    if proc.returncode != 0:
        remaining = proc.stderr.read() if proc.stderr else ""
        console.print(f"[red]phonopy failed (code {proc.returncode}):[/red]")
        for line in remaining.strip().split("\n")[-5:]:
            if line.strip():
                console.print(f"  [dim]{line.strip()[:120]}[/dim]")


def _parse_total_dos(path: Path):
    """Return (freq, dos) from total_dos.dat."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def _parse_partial_dos(path: Path):
    """Return (freq, pdos_per_element, labels) from partial_dos.dat."""
    # Read header to get column labels
    labels: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                parts = line.lstrip("#").strip().split()
                labels = parts[1:]  # skip "frequency(THz)"
                break
    data = np.loadtxt(path)
    freq = data[:, 0]
    pdos = data[:, 1:]
    # Filter out "total" column for plotting
    plot_cols = [i for i, lbl in enumerate(labels) if lbl.lower() != "total"]
    plot_labels = [labels[i] for i in plot_cols]
    return freq, pdos[:, plot_cols], plot_labels


def _plot_phonon_bands(console, fmin: float | None = None, fmax: float | None = None) -> None:
    """Parse band.yaml and plot phonon dispersion."""
    xs, freqs, tick_pos, tick_labels = _parse_band_yaml(Path("band.yaml"))

    f_lo = fmin if fmin is not None else float(np.min(freqs)) * 1.05
    f_hi = fmax if fmax is not None else float(np.max(freqs)) * 1.05

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    n_bands = freqs.shape[1]
    for b in range(n_bands):
        ax.plot(xs, freqs[:, b], color="#2c7bb6", lw=0.8, alpha=0.85)

    ax.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.5)

    # Segment boundary lines
    for xp in tick_pos:
        ax.axvline(xp, color="black", lw=0.8)

    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(f_lo, f_hi)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=11)
    ax.set_ylabel("Frequency (THz)", fontsize=12)
    ax.set_title("Phonon Dispersion", fontsize=13)
    ax.yaxis.grid(True, ls=":", alpha=0.4)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)
    ax.tick_params(axis="both", direction="out")

    fig.tight_layout(pad=1.2)
    out_png = "phonon_band.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    console.print(f"  [green]✓ {out_png}[/green]")

    # Save pho.dat (gnuplot-compatible)
    with open("pho.dat", "w") as f:
        for iq in range(len(xs)):
            freqs_str = "  ".join(f"{freqs[iq, ib]:.6f}" for ib in range(n_bands))
            f.write(f"{xs[iq]:.8f}  {freqs_str}\n")
    console.print("  [green]✓ pho.dat[/green]")


def _plot_phonon_combined(console, tdos_path: Path, pdos_path: Path,
                          fmin: float | None = None, fmax: float | None = None) -> None:
    """Combined figure: band + PDOS (if available), or band + TDOS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    xs, freqs, tick_pos, tick_labels = _parse_band_yaml(Path("band.yaml"))
    has_pdos = pdos_path.exists()
    has_tdos = tdos_path.exists()

    if not has_pdos and not has_tdos:
        raise FileNotFoundError("Neither partial_dos.dat nor total_dos.dat found")

    freq_t, dos_t = (None, None)
    if has_tdos:
        freq_t, dos_t = _parse_total_dos(tdos_path)

    freq_p, pdos_data, pdos_labels = (None, None, None)
    if has_pdos:
        freq_p, pdos_data, pdos_labels = _parse_partial_dos(pdos_path)

    f_lo = fmin if fmin is not None else float(np.min(freqs)) * 1.05
    f_hi = fmax if fmax is not None else float(np.max(freqs)) * 1.05

    # Layout: 2 panels (band + DOS/PDOS)
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1.5], wspace=0.06)
    ax_band = fig.add_subplot(gs[0])
    ax_dos = fig.add_subplot(gs[1])

    # ── Band dispersion ──
    n_bands = freqs.shape[1]
    for b in range(n_bands):
        ax_band.plot(xs, freqs[:, b], color="#2c7bb6", lw=0.8, alpha=0.85)
    ax_band.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.5)
    for xp in tick_pos:
        ax_band.axvline(xp, color="black", lw=0.8)
    ax_band.set_xlim(xs[0], xs[-1])
    ax_band.set_ylim(f_lo, f_hi)
    ax_band.set_xticks(tick_pos)
    ax_band.set_xticklabels(tick_labels, fontsize=11)
    ax_band.set_ylabel("Frequency (THz)", fontsize=12)
    ax_band.set_title("Phonon Dispersion", fontsize=13, pad=8)
    for spine in ax_band.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)
    ax_band.tick_params(axis="both", direction="out")
    ax_band.yaxis.grid(True, ls=":", alpha=0.4)

    # ── PDOS (preferred) or TDOS ──
    pdos_colors = ["#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                   "#F4A261", "#264653", "#A8DADC", "#6A0572"]
    if has_pdos:
        max_val = max(np.max(pdos_data[:, i]) for i in range(pdos_data.shape[1]))
        max_val = max_val if max_val > 0 else 1.0
        for i in range(pdos_data.shape[1]):
            label = pdos_labels[i] if i < len(pdos_labels) else f"atom {i + 1}"
            c = pdos_colors[i % len(pdos_colors)]
            dos_norm = pdos_data[:, i] / max_val
            ax_dos.fill_betweenx(freq_p, 0, dos_norm, color=c, alpha=0.30, lw=0)
            ax_dos.plot(dos_norm, freq_p, color=c, lw=1.3, label=label)
        ax_dos.legend(loc="upper right", fontsize=9, framealpha=0.7)
        ax_dos.set_xlabel("pDOS (norm.)", fontsize=11)
        ax_dos.set_title("Projected DOS", fontsize=13, pad=8)
    else:
        ax_dos.fill_betweenx(freq_t, 0, dos_t, color="#2c7bb6", alpha=0.35, lw=0)
        ax_dos.plot(dos_t, freq_t, color="#2c7bb6", lw=1.2)
        ax_dos.set_xlabel("DOS", fontsize=11)
        ax_dos.set_title("Total DOS", fontsize=13, pad=8)
    ax_dos.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.5)
    ax_dos.set_ylim(f_lo, f_hi)
    ax_dos.set_yticklabels([])
    ax_dos.yaxis.grid(True, ls=":", alpha=0.4)
    for spine in ax_dos.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)
    ax_dos.tick_params(axis="both", direction="out")

    fig.suptitle("Phonon Properties", fontsize=15, y=1.01)
    fig.tight_layout()
    out_png = "phonon_combined.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    console.print(f"  [green]✓ {out_png}[/green] (combined view)")


@task(1502, category="Lattice Dynamics", name="Phonon DOS",
      description="Compute total phonon DOS from FORCE_SETS",
      cli_args=[
          {"name": "--mesh", "type": str, "default": "8 8 8",
           "help": "q-point mesh for DOS (e.g. '8 8 8')"},
          {"name": "--fmin", "type": float, "default": None,
           "help": "Minimum frequency for plot x-axis (auto if unset)"},
          {"name": "--fmax", "type": float, "default": None,
           "help": "Maximum frequency for plot x-axis (auto if unset)"},
          {"name": "--no-plot", "action": "store_true", "default": False,
           "help": "Skip plotting"},
      ])
def task_phonon_dos(args: list[str] | None = None, interactive: bool = True,
                    parsed_args=None) -> None:
    """Compute total phonon DOS from FORCE_SETS using phonopy Python API."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Phonon DOS ===[/bold cyan]")
    console.print()

    if not Path("FORCE_SETS").exists():
        console.print("[red]FORCE_SETS not found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1501' first to extract forces.[/dim]")
        return

    fmin_plot = getattr(parsed_args, "fmin", None) if parsed_args else None
    fmax_plot = getattr(parsed_args, "fmax", None) if parsed_args else None

    # MESH
    setup_info: dict = {}
    setup_file = Path("phonopy_setup.json")
    if setup_file.exists():
        try:
            import json as _json
            setup_info = _json.loads(setup_file.read_text())
        except Exception:
            pass
    dim_default = setup_info.get("dim", "2 2 2")

    # Read default mesh from existing band.conf if available
    _mesh_default = "8 8 8"
    _bc = Path("band.conf")
    if _bc.exists():
        import re as _re
        _txt = _bc.read_text()
        _m = _re.search(r"MESH\s*=\s*(\d+\s+\d+\s+\d+)", _txt)
        if _m:
            _mesh_default = _m.group(1)
    if interactive:
        mesh_s = _prompt(console, "q-point mesh for DOS", _mesh_default)
    elif parsed_args:
        mesh_s = parsed_args.mesh
    else:
        mesh_s = "8 8 8"

    try:
        mesh = [int(x) for x in mesh_s.split()]
    except ValueError:
        console.print("[red]Mesh must be 3 integers, e.g. '8 8 8'[/red]")
        return

    console.print(f"  MESH: {mesh_s}")
    console.print("  Computing via phonopy Python API ...")

    try:
        from phonopy import load
        fc_kwargs = {"calculator": "abacus"}
        if Path("FORCE_CONSTANTS").exists():
            fc_kwargs["force_constants_filename"] = "FORCE_CONSTANTS"
            fc_kwargs["produce_fc"] = False
        else:
            fc_kwargs["force_sets_filename"] = "FORCE_SETS"
            fc_kwargs["produce_fc"] = True
        phonon = load(
            supercell_matrix=[int(x) for x in dim_default.split()],
            primitive_matrix="auto",
            unitcell_filename="STRU",
            is_nac=False,
            **fc_kwargs,
        )
        phonon.run_mesh(mesh)
        phonon.run_total_dos()
        td = phonon._total_dos
        freq = td.frequency_points
        dos = td.dos

        np.savetxt("total_dos.dat", np.column_stack([freq, dos]),
                   header="frequency(THz)  DOS")
        console.print("  [green]✓ total_dos.dat[/green]")
    except Exception as e:
        console.print(f"[red]DOS computation failed: {e}[/red]")
        return

    # Plot
    do_plot = not (parsed_args and parsed_args.no_plot)
    if interactive and do_plot:
        answer = _prompt_choice(console, "Generate DOS plot?", ["Yes", "No"], "Yes")
        do_plot = "Yes" in answer
        if do_plot:
            fmin_s = _prompt(console,
                f"Frequency min (range: {freq[0]:.1f} ~ {freq[-1]:.1f} THz, enter for auto)", "")
            fmax_s = _prompt(console,
                f"Frequency max (range: {freq[0]:.1f} ~ {freq[-1]:.1f} THz, enter for auto)", "")
            fmin_plot = float(fmin_s) if fmin_s.strip() else None
            fmax_plot = float(fmax_s) if fmax_s.strip() else None

    if do_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from abacuscopilot.plotting.style import load_style_from_config
            load_style_from_config()

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.plot(freq, dos, color="#2c7bb6", lw=1.5, label="Total")
            ax.fill_between(freq, 0, dos, color="#2c7bb6", alpha=0.15)
            ax.axvline(0, color="gray", lw=0.7, ls="--", alpha=0.6)
            f_lo_dos = fmin_plot if fmin_plot is not None else float(np.min(freq))
            f_hi_dos = fmax_plot if fmax_plot is not None else float(np.max(freq)) * 1.05
            ax.set_xlim(f_lo_dos, f_hi_dos)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("Frequency (THz)", fontsize=12)
            ax.set_ylabel("DOS (states/THz)", fontsize=12)
            ax.set_title("Phonon Density of States", fontsize=13)
            ax.grid(True, ls=":", alpha=0.4)
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_visible(True)
            ax.tick_params(axis="both", direction="out")
            fig.tight_layout(pad=1.2)
            out_png = "phonon_dos.png"
            fig.savefig(out_png, dpi=300, bbox_inches="tight")
            plt.close(fig)
            console.print(f"  [green]✓ {out_png}[/green]")
        except Exception as e:
            console.print(f"  [yellow]! Plot failed: {e}[/yellow]")

    console.print()


@task(1503, category="Lattice Dynamics", name="Phonon PDOS",
      description="Compute projected phonon DOS (atom-resolved) from FORCE_SETS",
      cli_args=[
          {"name": "--mesh", "type": str, "default": "8 8 8",
           "help": "q-point mesh for PDOS (e.g. '8 8 8')"},
          {"name": "--fmin", "type": float, "default": None,
           "help": "Minimum frequency for plot x-axis (auto if unset)"},
          {"name": "--fmax", "type": float, "default": None,
           "help": "Maximum frequency for plot x-axis (auto if unset)"},
          {"name": "--no-plot", "action": "store_true", "default": False,
           "help": "Skip plotting"},
      ])
def task_phonon_pdos(args: list[str] | None = None, interactive: bool = True,
                     parsed_args=None) -> None:
    """Compute projected phonon DOS from FORCE_SETS using phonopy Python API."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Phonon Projected DOS ===[/bold cyan]")
    console.print()

    if not Path("FORCE_SETS").exists():
        console.print("[red]FORCE_SETS not found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1501' first to extract forces.[/dim]")
        return

    fmin_plot = getattr(parsed_args, "fmin", None) if parsed_args else None
    fmax_plot = getattr(parsed_args, "fmax", None) if parsed_args else None

    # Read STRU for species info
    species_list: list[str] = []
    species_indices: dict[str, list[int]] = {}  # {species: [atom indices in STRU]}
    try:
        from abacuscopilot.io.stru_file import read_stru as _read_stru2
        s = _read_stru2("STRU")
        species_list = list(s.species_order)
        for i, atom in enumerate(s.atoms):
            sp = atom.species
            species_indices.setdefault(sp, []).append(i + 1)
    except Exception:
        pass

    # MESH + DIM
    setup_info: dict = {}
    setup_file = Path("phonopy_setup.json")
    if setup_file.exists():
        try:
            import json as _json
            setup_info = _json.loads(setup_file.read_text())
        except Exception:
            pass
    dim_default = setup_info.get("dim", "2 2 2")

    _mesh_default = "8 8 8"
    _bc = Path("band.conf")
    if _bc.exists():
        import re as _re
        _txt = _bc.read_text()
        _m = _re.search(r"MESH\s*=\s*(\d+\s+\d+\s+\d+)", _txt)
        if _m:
            _mesh_default = _m.group(1)
    if interactive:
        mesh_s = _prompt(console, "q-point mesh for PDOS", _mesh_default)
    elif parsed_args:
        mesh_s = parsed_args.mesh
    else:
        mesh_s = "8 8 8"

    try:
        mesh = [int(x) for x in mesh_s.split()]
    except ValueError:
        console.print("[red]Mesh must be 3 integers, e.g. '8 8 8'[/red]")
        return

    # Species selection
    selected_labels: list[str] = []
    if interactive and len(species_list) > 1:
        opts = ["All"] + species_list
        answer = _prompt_choice(console, "Project onto which species?", opts, "All")
        if "All" in answer:
            selected_labels = species_list
        else:
            selected_labels = [sp for sp in species_list if sp in answer]
    else:
        selected_labels = species_list  # single species or non-interactive

    console.print(f"  MESH: {mesh_s}")
    console.print(f"  Projecting: {', '.join(selected_labels)}")
    console.print("  Computing via phonopy Python API ...")

    try:
        from phonopy import load
        fc_kwargs = {"calculator": "abacus"}
        if Path("FORCE_CONSTANTS").exists():
            fc_kwargs["force_constants_filename"] = "FORCE_CONSTANTS"
            fc_kwargs["produce_fc"] = False
        else:
            fc_kwargs["force_sets_filename"] = "FORCE_SETS"
            fc_kwargs["produce_fc"] = True
        phonon = load(
            supercell_matrix=[int(x) for x in dim_default.split()],
            primitive_matrix="auto",
            unitcell_filename="STRU",
            is_nac=False,
            **fc_kwargs,
        )
        phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        phonon.run_projected_dos()

        pd = phonon._pdos
        freq = pd.frequency_points
        pdos_data = pd.projected_dos  # shape (n_atoms_prim, n_freq)

        # Get primitive cell element mapping from phonopy
        prim_cell = phonon.primitive
        prim_symbols = prim_cell.symbols  # list of element symbols per atom
        # Build per-element grouping: {symbol: [atom_indices]}
        elem_groups: dict[str, list[int]] = {}
        for i, sym in enumerate(prim_symbols):
            elem_groups.setdefault(sym, []).append(i)

        # Save: freq + per-element PDOS (grouped) + total
        total_pdos = pdos_data.sum(axis=0)
        col_labels: list[str] = []
        col_arrays: list[np.ndarray] = [freq]
        for sym, indices in elem_groups.items():
            grouped = pdos_data[indices, :].sum(axis=0)
            col_arrays.append(grouped)
            col_labels.append(sym)
        col_arrays.append(total_pdos)
        col_labels.append("total")
        header = "frequency(THz)  " + "  ".join(col_labels)
        np.savetxt("partial_dos.dat", np.column_stack(col_arrays), header=header)
        console.print("  [green]✓ partial_dos.dat[/green]")
    except Exception as e:
        console.print(f"[red]PDOS computation failed: {e}[/red]")
        return

    # Plot
    do_plot = not (parsed_args and parsed_args.no_plot)
    if interactive and do_plot:
        answer = _prompt_choice(console, "Generate PDOS plot?", ["Yes", "No"], "Yes")
        do_plot = "Yes" in answer
        if do_plot:
            fmin_s = _prompt(console,
                f"Frequency min (range: {freq[0]:.1f} ~ {freq[-1]:.1f} THz, enter for auto)", "")
            fmax_s = _prompt(console,
                f"Frequency max (range: {freq[0]:.1f} ~ {freq[-1]:.1f} THz, enter for auto)", "")
            fmin_plot = float(fmin_s) if fmin_s.strip() else None
            fmax_plot = float(fmax_s) if fmax_s.strip() else None

    if do_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from abacuscopilot.plotting.style import load_style_from_config
            load_style_from_config()

            fig, ax = plt.subplots(figsize=(7, 5))
            pdos_colors = ["#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                           "#F4A261", "#264653", "#A8DADC", "#6A0572"]

            # Plot by element group (sum over atoms of same element)
            color_idx = 0
            for sym, indices in elem_groups.items():
                if selected_labels and sym not in selected_labels:
                    continue
                grouped = pdos_data[indices, :].sum(axis=0)
                c = pdos_colors[color_idx % len(pdos_colors)]
                color_idx += 1
                ax.plot(freq, grouped, color=c, lw=1.3, label=sym)
                ax.fill_between(freq, 0, grouped, color=c, alpha=0.12)

            # Total PDOS (black)
            ax.plot(freq, total_pdos, color="black", lw=1.8, label="Total", zorder=5)

            ax.axvline(0, color="gray", lw=0.7, ls="--", alpha=0.6)
            f_lo_pdos = fmin_plot if fmin_plot is not None else float(np.min(freq))
            f_hi_pdos = fmax_plot if fmax_plot is not None else float(np.max(freq)) * 1.05
            ax.set_xlim(f_lo_pdos, f_hi_pdos)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("Frequency (THz)", fontsize=12)
            ax.set_ylabel("PDOS (states/THz)", fontsize=12)
            ax.set_title("Phonon Projected Density of States", fontsize=13)
            ax.legend(fontsize=10, framealpha=0.8)
            ax.grid(True, ls=":", alpha=0.4)
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_visible(True)
            ax.tick_params(axis="both", direction="out")
            fig.tight_layout(pad=1.2)
            out_png = "phonon_pdos.png"
            fig.savefig(out_png, dpi=300, bbox_inches="tight")
            plt.close(fig)
            console.print(f"  [green]✓ {out_png}[/green]")
        except Exception as e:
            console.print(f"  [yellow]! Plot failed: {e}[/yellow]")

    console.print()


@task(1504, category="Lattice Dynamics", name="Phonon Combined",
      description="Combined phonon band + DOS + PDOS figure (side-by-side)",
      cli_args=[
          {"name": "--fmin", "type": float, "default": None,
           "help": "Minimum frequency for plot y-axis (auto if unset)"},
          {"name": "--fmax", "type": float, "default": None,
           "help": "Maximum frequency for plot y-axis (auto if unset)"},
      ])
def task_phonon_combined(args: list[str] | None = None, interactive: bool = True,
                         parsed_args=None) -> None:
    """Generate combined phonon band + DOS + PDOS figure from existing data."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Phonon Combined Plot ===[/bold cyan]")
    console.print()

    band_yaml = Path("band.yaml")
    tdos_dat = Path("total_dos.dat")
    pdos_dat = Path("partial_dos.dat")

    if not band_yaml.exists():
        console.print("[red]band.yaml not found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1501' first.[/dim]")
        return
    if not tdos_dat.exists() and not pdos_dat.exists():
        console.print("[red]Neither total_dos.dat nor partial_dos.dat found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1502' or '-task 1503' first.[/dim]")
        return

    fmin_plot = getattr(parsed_args, "fmin", None) if parsed_args else None
    fmax_plot = getattr(parsed_args, "fmax", None) if parsed_args else None

    if interactive and fmin_plot is None:
        try:
            _, freqs_tmp, _, _ = _parse_band_yaml(band_yaml)
            fmin_s = _prompt(console,
                f"Frequency min (auto: {float(np.min(freqs_tmp))*1.05:.1f} THz, enter for auto)", "")
            fmax_s = _prompt(console,
                f"Frequency max (auto: {float(np.max(freqs_tmp))*1.05:.1f} THz, enter for auto)", "")
            fmin_plot = float(fmin_s) if fmin_s.strip() else None
            fmax_plot = float(fmax_s) if fmax_s.strip() else None
        except Exception:
            pass

    try:
        _plot_phonon_combined(console, tdos_dat, pdos_dat, fmin_plot, fmax_plot)
    except Exception as e:
        console.print(f"[red]Combined plot failed: {e}[/red]")

    console.print()


# =============================================================================
# MLP-SSCHA: temperature-dependent phonons via machine-learned potential
# =============================================================================

@task(1505, category="Lattice Dynamics", name="MLP-SSCHA Setup",
      description="Generate random displacements for MLP training (SSCHA workflow)",
      cli_args=[
          {"name": "--dim", "type": str, "default": "2 2 2",
           "help": "Supercell dimensions (e.g. '2 2 2')"},
          {"name": "--rd", "type": int, "default": 1000,
           "help": "Number of random displacements (default 1000)"},
          {"name": "--amin", "type": float, "default": 0.03,
           "help": "Minimum displacement amplitude (default 0.03)"},
          {"name": "--amax", "type": float, "default": 1.5,
           "help": "Maximum displacement amplitude (default 1.5)"},
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw, dp"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_mlp_sscha_setup(args=None, interactive=True, parsed_args=None):
    """Generate random displacements for MLP-SSCHA phonon training.

    Uses ``phonopy-init --rd N`` to create many supercells with random
    atomic displacements (not symmetry-adapted).  These are used as
    training data for ``pypolymlp`` to fit a machine-learned potential,
    which then feeds the SSCHA (self-consistent harmonic approximation)
    to compute temperature-dependent phonons.
    """
    import shutil as _shutil

    console = _get_console()

    console.print()
    console.print("[bold cyan]=== MLP-SSCHA Setup ===[/bold cyan]")
    console.print("[dim]Random displacements → MLP training → SSCHA phonons[/dim]")
    console.print()

    # --- Parameters ---
    if interactive:
        dim_s = _prompt(console, "Supercell dimensions", "2 2 2")
    elif parsed_args:
        dim_s = parsed_args.dim
    else:
        dim_s = "2 2 2"

    if interactive:
        rd_s = _prompt(console, "Number of random displacements", "300")
    elif parsed_args:
        rd_s = str(parsed_args.rd)
    else:
        rd_s = "300"
    n_rd = int(rd_s)

    if interactive:
        amin_s = _prompt(console, "Min displacement amplitude", "0.03")
        amax_s = _prompt(console, "Max displacement amplitude", "1.5")
    elif parsed_args:
        amin_s = str(parsed_args.amin)
        amax_s = str(parsed_args.amax)
    else:
        amin_s, amax_s = "0.03", "1.5"

    console.print(f"  Supercell: {dim_s},  displacements: {n_rd}")
    console.print(f"  Amplitude: {amin_s} ~ {amax_s}")
    if n_rd > 100:
        console.print("  [bold yellow]⚠  MLP-SSCHA works best for small primitive cells (Al, Si, GaN, MgO).[/bold yellow]")
        console.print("  [bold yellow]    Large/complex systems (30+ atoms) may run out of memory.[/bold yellow]")

    # --- Basis type ---
    from abacuscopilot.preprocessing.input_tasks import (
        _prompt_choice as _pc, _get_template, _apply_template,
        _ask_lcao_solver, _apply_solver_override,
        _auto_prepare_files,
    )
    if interactive:
        basis = _pc(console, "Basis type", ["lcao", "pw", "dp"], "lcao")
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw", "dp") else "lcao"
    else:
        basis = "lcao"

    # --- INPUT params ---
    from abacuscopilot.core.models import InputParams
    params = InputParams()
    params.suffix = "ABACUS"
    params.calculation = "scf"
    params.cal_force = 1
    # Point to parent dir for pseudo/orbital to avoid copying to every subdir
    params.pseudo_dir = "../"
    if basis != "dp":
        params.orbital_dir = "../"

    if basis == "dp":
        params.esolver_type = "dp"
        pot_file = "graph.pb"
        if interactive:
            pot_file = _prompt(console, "DP model file", "graph.pb")
        if not Path(pot_file).exists():
            console.print(f"[red]DP model not found: {pot_file}[/red]")
            return
        params.pot_file = pot_file
    else:
        template = _get_template(basis, "scf")
        if template:
            _apply_template(params, template)
        if interactive and basis == "lcao":
            _ask_lcao_solver(console, params)

        # Force cal_force into INPUT (default=1, template overwrites template_keys)
        if "cal_force" not in params.extras.get("_template_keys", []):
            params.extras.setdefault("_template_keys", []).append("cal_force")

        # --- Functional + D3 (DFT only, skip for DP) ---
        if interactive and basis != "dp":
            use_func = _pc(console, "Exchange-correlation functional", ["PBEsol", "PBE"], "PBEsol")
            if "PBEsol" in use_func:
                params.dft_functional = "pbesol"

            use_d3 = _pc(console, "D3 dispersion correction", ["No", "d3_0 (zero-damping)", "d3_bj (Becke-Johnson)"], "No")
            if "d3_0" in use_d3:
                params.vdw_method = "d3_0"
            elif "d3_bj" in use_d3:
                params.vdw_method = "d3_bj"
            if params.vdw_method != "none" and params.dft_functional == "pbesol":
                console.print("  [bold yellow]⚠  ABACUS D3 does not support dft_functional=pbesol.[/bold yellow]")
                console.print("  [dim]    Setting dft_functional to 'pbe' for D3 compatibility.[/dim]")
                console.print("  [dim]    PBEsol reuses PBE D3 parameters — this is standard practice.[/dim]")
                params.dft_functional = "pbe"

            # Force only non-default values into INPUT
            for _k, _def in (("dft_functional", "pbe"), ("vdw_method", "none")):
                if getattr(params, _k) != _def:
                    params.extras.setdefault("_template_keys", []).append(_k)
                    # Remove redundant comment hint (already active)
                    _hints = params.extras.get("_comment_hints", {})
                    if isinstance(_hints, dict):
                        _hints.pop(_k, None)

    # Remove _comment_hints for cal_force/cal_stress (already active)
    _hints = params.extras.get("_comment_hints", {})
    for _k in ("cal_force", "cal_stress"):
        if isinstance(_hints, dict):
            _hints.pop(_k, None)

    # --- Read STRU ---
    from abacuscopilot.io.stru_file import read_stru, write_stru
    stru_path = Path("STRU")
    if not stru_path.exists():
        console.print("[red]No STRU file found.[/red]")
        return

    structure = read_stru(stru_path)
    if structure.coordinate_type != "Direct":
        if "Cartesian" in structure.coordinate_type:
            cell_bohr = structure.lattice.cell
            cell_inv = np.linalg.inv(cell_bohr)
            for atom in structure.atoms:
                pos = atom.position.copy()
                if "angstrom" in structure.coordinate_type.lower():
                    from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
                    pos = pos * ANGSTROM_TO_BOHR
                atom.position = pos @ cell_inv
        structure.coordinate_type = "Direct"

    # Copy pseudopotentials/orbitals (DFT only)
    pseudo_files: list[Path] = []
    orbital_files: list[Path] = []
    if basis != "dp":
        _auto_prepare_files(console, params, interactive)
        structure = read_stru(stru_path)
        for species in structure.species_order:
            for pat in Path(".").glob(f"{species}_*.upf"):
                if pat not in pseudo_files:
                    pseudo_files.append(pat)
            if basis == "lcao":
                for pat in Path(".").glob(f"{species}_*.orb"):
                    if pat not in orbital_files:
                        orbital_files.append(pat)

    # --- Run phonopy-init --rd ---
    console.print()
    console.print(f"  Running: [dim]phonopy-init --rd {n_rd} --dim {dim_s} --amin {amin_s} --amax {amax_s} --abacus[/dim]")
    rd_result = _run_phonopy_init(
        ["--rd", str(n_rd), "--dim", dim_s,
         "--amin", amin_s, "--amax", amax_s, "--abacus"],
        timeout=300,
    )
    if rd_result.returncode != 0:
        console.print("[red]phonopy-init --rd failed:[/red]")
        console.print(rd_result.stderr)
        return
    if rd_result.stderr:
        for line in rd_result.stderr.strip().split("\n")[:5]:
            console.print(f"  [dim]{line}[/dim]")

    # Find generated STRU-XXX files
    stru_files = sorted(Path(".").glob("STRU-*"))
    if not stru_files:
        console.print("[red]No STRU-* files generated.[/red]")
        return
    console.print(f"  Generated: {len(stru_files)} structures")

    # --- Handle sub.abacus (sub.abacus-dp for DP) ---
    config = load_config()
    sub_src: Path | None = None
    if basis == "dp":
        dp_path = config.get("paths", {}).get("sub_script_dp", "")
        if dp_path:
            p = Path(dp_path).expanduser()
            if p.exists():
                sub_src = p
    if sub_src is None:
        sub_path = config.get("paths", {}).get("sub_script", "")
        if sub_path:
            p = Path(sub_path).expanduser()
            if p.exists():
                sub_src = p
    if sub_src:
        console.print(f"  sub.abacus: [dim]{sub_src}[/dim]")
    else:
        console.print("  sub.abacus: [dim]skipped[/dim]")

    # --- Create disp-XXX directories ---
    console.print()
    console.print("[bold]Creating directories:[/bold]")

    from abacuscopilot.io.input_file import write_input
    for i, sf in enumerate(stru_files, 1):
        dir_name = f"disp-{i:04d}"
        dir_path = Path(dir_name)
        dir_path.mkdir(exist_ok=True)
        _shutil.copy2(sf, dir_path / "STRU")

        if basis == "dp":
            with open(dir_path / "INPUT", "w") as f:
                f.write("INPUT_PARAMETERS\n")
                f.write("calculation          scf\n")
                f.write("esolver_type         dp\n")
                pot_ref = f"../{params.pot_file}" if not params.pot_file.startswith("/") else params.pot_file
                f.write(f"pot_file             {pot_ref}\n")
                f.write("cal_force            1\n")
                f.write("symmetry             1\n")
        else:
            write_input(params, dir_path / "INPUT")

        if sub_src and sub_src.exists():
            _shutil.copy2(sub_src, dir_path / "sub.abacus")

        if i % 100 == 0:
            console.print(f"  [dim]... {i}/{len(stru_files)}[/dim]")

    # --- Cleanup parent sub.abacus ---
    local_sub = Path("sub.abacus")
    if local_sub.exists() and sub_src:
        local_sub.unlink()

    # --- Save setup ---
    import json as _json
    _setup = {"dim": dim_s, "basis": basis, "n_rd": n_rd,
              "amin": float(amin_s), "amax": float(amax_s),
              "n_structures": len(stru_files)}
    with open("phonopy_setup.json", "w") as f:
        _json.dump(_setup, f, indent=2)

    # Keep UPF/ORB in parent dir — all disp dirs reference them via ../pseudo_dir
    console.print()
    console.print(f"[green]✓ MLP-SSCHA setup complete: {len(stru_files)} structures[/green]")
    console.print(f"  Run ABACUS in all disp-*/ dirs, then: [bold]abacuscopilot -task 1506[/bold]")
    console.print()


@task(1506, category="Lattice Dynamics", name="MLP Training",
      description="Extract forces from disp dirs and train MLP (pypolymlp) for SSCHA",
      cli_args=[
          {"name": "--ntrain", "type": int, "default": 100,
           "help": "Number of structures for training (default 100)"},
          {"name": "--ntest", "type": int, "default": 20,
           "help": "Number of structures for testing (default 20)"},
          {"name": "--range", "type": str, "default": "",
           "help": "Range of disp dirs to use, e.g. '1-120'"},
      ])
def task_mlp_training(args=None, interactive=True, parsed_args=None):
    """Extract forces from ABACUS outputs and train an MLP via pypolymlp.

    Step 1: ``phonopy --sp -f disp-*/OUT*/running_scf.log``
    Step 2: Compress ``phonopy_params.yaml`` → ``.xz``
    Step 3: ``phonopy-load ... --pypolymlp`` → ``polymlp.yaml``
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== MLP Training (pypolymlp) ===[/bold cyan]")
    console.print()

    # --- Collect log files ---
    disp_dirs = sorted(Path(".").glob("disp-*"))
    if not disp_dirs:
        console.print("[red]No disp-*/ directories found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1505' first.[/dim]")
        return

    log_paths = []
    range_str = getattr(parsed_args, "range", "") if parsed_args else ""
    if range_str:
        import re as _re
        m = _re.match(r"(\d+)-(\d+)", range_str)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            disp_dirs = [d for d in disp_dirs if lo <= int(d.name.split("-")[-1]) <= hi]

    for dd in disp_dirs:
        candidates = list(dd.glob("OUT*/running_scf.log")) + list(dd.glob("OUT*/running_relax.log"))
        if candidates:
            log_paths.append(candidates[0])
    if not log_paths:
        console.print("[red]No ABACUS output logs found.[/red]")
        return
    console.print(f"  Using {len(log_paths)}/{len(disp_dirs)} completed calculations")

    # --- Step 1: phonopy --sp -f ---
    console.print()
    console.print("[bold]Step 1: Extract forces (phonopy-init --sp -f)[/bold]")
    sp_args = ["--sp", "-f"] + [str(p) for p in log_paths]
    result = _run_phonopy_init(sp_args, timeout=300)
    if result.returncode != 0:
        console.print("[red]Force extraction failed:[/red]")
        console.print(result.stderr[-500:])
        return
    params_file = Path("phonopy_params.yaml")
    if not params_file.exists():
        console.print("[red]phonopy_params.yaml was not created.[/red]")
        return
    console.print("  [green]✓ phonopy_params.yaml[/green]")

    # --- Step 1.5: Inject supercell energies (ABACUS doesn't output them) ---
    console.print()
    console.print("[bold]Step 1.5: Inject supercell energies[/bold]")
    try:
        import yaml as _yaml
        with open("phonopy_params.yaml") as f:
            _pp = _yaml.safe_load(f)
        _ds = _pp.get("dataset", {})
        _n_disp = len(_ds.get("forces", []) if isinstance(_ds, dict) else [])
        _energies = []
        for _i in range(_n_disp):
            _dn = f"disp-{_i + 1:04d}"
            _log = Path(_dn) / "OUT.ABACUS" / "running_scf.log"
            if not _log.exists():
                _log = Path(_dn) / "OUT.ABACUS" / "running_relax.log"
            _ev = 0.0
            if _log.exists():
                with open(_log) as _lf:
                    for _line in _lf:
                        _m = re.search(r"!FINAL_ETOT_IS\s+([\-\d\.Ee+]+)", _line)
                        if _m:
                            _ev = float(_m.group(1))
                            break
            _energies.append(_ev)
        if _energies and isinstance(_ds, dict):
            _ds["supercell_energies"] = _energies
            with open("phonopy_params.yaml", "w") as f:
                _yaml.safe_dump(_pp, f, default_flow_style=False)
            console.print(f"  [green]✓ Injected energies for {_n_disp} supercells[/green]")
        else:
            console.print("  [yellow]! Could not inject energies[/yellow]")
    except Exception as _e:
        console.print(f"  [yellow]! Energy injection failed: {_e}[/yellow]")

    # --- Step 2: Compress ---
    console.print()
    console.print("[bold]Step 2: Compress dataset[/bold]")
    import subprocess as _sp
    xz_file = Path("phonopy_params.yaml.xz")
    _sp.run(["xz", "-f", "phonopy_params.yaml"], check=False)
    if xz_file.exists():
        console.print("  [green]✓ phonopy_params.yaml.xz[/green]")
    else:
        console.print("  [yellow]! xz compression may have failed, trying without[/yellow]")

    # --- Step 3: MLP training ---
    console.print()
    console.print("[bold]Step 3: Train MLP (pypolymlp)[/bold]")
    ntrain = getattr(parsed_args, "ntrain", 100) if parsed_args else 100
    ntest = getattr(parsed_args, "ntest", 20) if parsed_args else 20
    if interactive:
        ntrain_s = _prompt(console, "Training set size", str(ntrain))
        ntrain = int(ntrain_s) if ntrain_s.strip() else ntrain
        ntest_s = _prompt(console, "Test set size", str(ntest))
        ntest = int(ntest_s) if ntest_s.strip() else ntest

    mlp_args = ["--pypolymlp", f"--mlp-params=ntrain={ntrain}, ntest={ntest}, gtinv_order=2"]
    if xz_file.exists():
        mlp_args = [str(xz_file)] + mlp_args
    else:
        mlp_args = [str(params_file)] + mlp_args

    console.print(f"  Running: [dim]phonopy-load ... --pypolymlp ntrain={ntrain}, ntest={ntest}[/dim]")
    result = _run_phonopy(mlp_args, timeout=600)
    if result.returncode != 0:
        console.print("[red]MLP training failed:[/red]")
        console.print(result.stderr[-500:])
        return
    if result.stderr:
        for line in result.stderr.strip().split("\n")[:5]:
            console.print(f"  [dim]{line}[/dim]")
    if Path("polymlp.yaml").exists():
        console.print("  [green]✓ polymlp.yaml[/green]")
    else:
        console.print("  [yellow]! polymlp.yaml not found, check logs[/yellow]")

    console.print()
    console.print("[green]✓ MLP training complete[/green]")
    console.print("  [dim]Next: abacuscopilot -task 1507 (SSCHA)[/dim]")
    console.print()


@task(1507, category="Lattice Dynamics", name="SSCHA Run",
      description="Run SSCHA iterations for temperature-dependent force constants",
      cli_args=[
          {"name": "--temperature", "type": float, "default": 300,
           "help": "Target temperature in K (default 300)"},
          {"name": "--iterations", "type": int, "default": 10,
           "help": "Number of SSCHA iterations (default 10)"},
          {"name": "--rd", "type": int, "default": 1000,
           "help": "Number of random displacements (default 1000)"},
      ])
def task_sscha_run(args=None, interactive=True, parsed_args=None):
    """Run SSCHA self-consistent iterations for temperature-dependent phonons.

    Requires ``polymlp.yaml`` and ``phonopy_params.yaml.xz`` from task 1506.
    Outputs ``phonopy_sscha_fc_N.yaml.xz`` for each iteration.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== SSCHA Run ===[/bold cyan]")
    console.print("[dim]Self-consistent harmonic approximation[/dim]")
    console.print()

    params_file = Path("phonopy_params.yaml.xz")
    if not params_file.exists():
        params_file = Path("phonopy_params.yaml")
    if not params_file.exists():
        console.print("[red]phonopy_params.yaml[.xz] not found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1506' first.[/dim]")
        return
    if not Path("polymlp.yaml").exists():
        console.print("[red]polymlp.yaml not found.[/red]")
        return

    if interactive:
        t_s = _prompt(console, "Target temperature (K)", "300")
        temperature = float(t_s) if t_s.strip() else 300.0
        n_iter_s = _prompt(console, "SSCHA iterations", "10")
        n_iter = int(n_iter_s) if n_iter_s.strip() else 10
        n_rd_s = _prompt(console, "Random displacements per iteration", "1000")
        n_rd = int(n_rd_s) if n_rd_s.strip() else 1000
    elif parsed_args:
        temperature = float(parsed_args.temperature)
        n_iter = int(parsed_args.iterations)
        n_rd = int(parsed_args.rd)
    else:
        temperature, n_iter, n_rd = 300.0, 10, 1000

    console.print(f"  Temperature: {temperature} K,  iterations: {n_iter},  rd: {n_rd}")
    console.print(f"  Running: [dim]phonopy-load ... --pypolymlp --sscha {n_iter} --rd-temperature {temperature} --rd {n_rd}[/dim]")
    console.print("  [dim]This may take a while for large systems...[/dim]")

    sscha_args = [
        str(params_file), "--pypolymlp",
        "--sscha", str(n_iter),
        "--rd-temperature", str(int(temperature)),
        "--rd", str(n_rd),
    ]
    # Stream stderr in real-time for progress
    import shutil as _sh, os as _os
    for _cand in [_sh.which("phonopy"), _sh.which("phonopy-load"),
                 Path(sys.executable).parent / "phonopy",
                 Path(sys.prefix) / "bin" / "phonopy",
                 Path(sys.executable).parent / "phonopy-load",
                 Path(sys.prefix) / "bin" / "phonopy-load"]:
        if _cand and Path(str(_cand)).exists():
            _exe = str(_cand)
            break
    else:
        _exe = sys.executable
        sscha_args = ["-m", "phonopy"] + sscha_args
    _env = _os.environ.copy()
    _env.setdefault("OMP_NUM_THREADS", str(_os.cpu_count() or 4))
    _proc = subprocess.Popen([_exe] + sscha_args, stderr=subprocess.PIPE,
                              text=True, bufsize=1, env=_env)
    _last = ""
    while _proc.poll() is None:
        _line = _proc.stderr.readline() if _proc.stderr else ""
        if _line:
            _s = _line.strip()
            if _s and _s != _last:
                console.print(f"  [dim]{_s[:140]}[/dim]")
                _last = _s
    _rc = _proc.returncode
    if _rc != 0:
        _rem = _proc.stderr.read() if _proc.stderr else ""
        console.print(f"[red]SSCHA failed (code {_rc}):[/red]")
        console.print(_rem[-500:])
        return

    fc_files = sorted(Path(".").glob("phonopy_sscha_fc_*.yaml.xz"))
    if fc_files:
        console.print(f"  [green]✓ {len(fc_files)} force constant files generated[/green]")
        for fc in fc_files:
            console.print(f"    [dim]{fc.name}[/dim]")
    else:
        console.print("  [yellow]! No SSCHA force constant files found[/yellow]")

    console.print()
    console.print("[green]✓ SSCHA complete[/green]")
    console.print(f"  Next: [bold]abacuscopilot -task 1508[/bold] (plot convergence)")
    console.print()


@task(1508, category="Lattice Dynamics", name="SSCHA Plot",
      description="Plot SSCHA convergence: band structures across iterations",
      cli_args=[
          {"name": "--fmin", "type": float, "default": None,
           "help": "Minimum frequency for plot y-axis (auto if unset)"},
          {"name": "--fmax", "type": float, "default": None,
           "help": "Maximum frequency for plot y-axis (auto if unset)"},
      ])
def task_sscha_plot(args=None, interactive=True, parsed_args=None):
    """Plot phonon band convergence across SSCHA iterations.

    Reads ``phonopy_sscha_fc_N.yaml.xz`` for N=1..max_iter,
    computes band structure for each, and overlays them on one plot.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== SSCHA Convergence Plot ===[/bold cyan]")
    console.print()

    fc_files = sorted(Path(".").glob("phonopy_sscha_fc_*.yaml.xz"))
    if not fc_files:
        console.print("[red]No phonopy_sscha_fc_*.yaml.xz files found.[/red]")
        console.print("[dim]Run 'abacuscopilot -task 1507' first.[/dim]")
        return

    console.print(f"  Found {len(fc_files)} SSCHA iteration files")

    # --- Compute band for each iteration ---
    band_files = []
    for fc in fc_files:
        band_out = Path(f"band_sscha_{fc.stem}.yaml")
        if not band_out.exists():
            console.print(f"  Computing band for {fc.name} ...")
            result = _run_phonopy(
                ["--band", "auto", "--band-points", "101", "-s", str(fc)],
                timeout=300,
            )
            if result.returncode == 0 and Path("band.yaml").exists():
                import shutil as _shutil
                _shutil.move("band.yaml", str(band_out))
        if band_out.exists():
            band_files.append(band_out)

    if not band_files:
        console.print("[red]No band structures computed.[/red]")
        return

    # Read auto frequency range from last iteration
    try:
        _last_freqs = []
        with open(band_files[-1]) as f:
            for line in f:
                m = re.search(r"frequency:\s*([-\d.]+)", line)
                if m:
                    _last_freqs.append(float(m.group(1)))
        _auto_min = min(_last_freqs) * 1.05 if _last_freqs else 0.0
        _auto_max = max(_last_freqs) * 1.05 if _last_freqs else 20.0
    except Exception:
        _auto_min, _auto_max = 0.0, 20.0

    fmin_plot = getattr(parsed_args, "fmin", None) if parsed_args else None
    fmax_plot = getattr(parsed_args, "fmax", None) if parsed_args else None

    if interactive and fmin_plot is None:
        fmin_s = _prompt(console,
            f"Frequency min (auto: {_auto_min:.1f} THz, enter for auto)", "")
        fmax_s = _prompt(console,
            f"Frequency max (auto: {_auto_max:.1f} THz, enter for auto)", "")
        fmin_plot = float(fmin_s) if fmin_s.strip() else None
        fmax_plot = float(fmax_s) if fmax_s.strip() else None

    # Plot overlay
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(band_files)))

    for i, bf in enumerate(band_files):
        try:
            xs, freqs, _, _ = _parse_band_yaml(bf)
            n_bands = freqs.shape[1]
            alpha = 0.3 if i < len(band_files) - 1 else 0.9
            lw = 0.5 if i < len(band_files) - 1 else 1.5
            for b in range(n_bands):
                ax.plot(xs, freqs[:, b], color=colors[i], lw=lw, alpha=alpha)
        except Exception:
            pass

    # Last iteration labels
    try:
        _, _, tick_pos, tick_labels = _parse_band_yaml(band_files[-1])
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_labels, fontsize=10)
        ax.set_xlim(0, tick_pos[-1])
    except Exception:
        pass

    ax.set_ylabel("Frequency (THz)", fontsize=12)
    ax.set_title("SSCHA Convergence", fontsize=13)
    f_lo = fmin_plot if fmin_plot is not None else _auto_min
    f_hi = fmax_plot if fmax_plot is not None else _auto_max
    ax.set_ylim(f_lo, f_hi)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    out_png = "sscha_convergence.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    console.print(f"  [green]✓ {out_png}[/green] ({len(band_files)} iterations)")

    console.print()


