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
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=timeout,
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
      description="Extract FORCE_SETS, compute phonon bands, and plot dispersion + DOS",
      cli_args=[
          {"name": "--dim", "type": str, "default": "2 2 2",
           "help": "Supercell dimensions (must match setup)"},
          {"name": "--mesh", "type": str, "default": "8 8 8",
           "help": "q-point mesh for DOS (e.g. '8 8 8')"},
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

    if interactive:
        mesh_s = _prompt(console, "q-point mesh for DOS/thermal", "8 8 8")
    elif parsed_args:
        mesh_s = parsed_args.mesh
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

    # --- Step 3: Compute phonon bands ---
    console.print()
    console.print("[bold]Step 3: Compute phonon bands[/bold]")
    console.print("  Running: [dim]phonopy --band ... band.conf[/dim]")

    # Split k-path + labels into individual arguments
    kpoints = kpath_str.split()
    klabels = klabels_str.split()
    band_args = (
        ["--band"] + kpoints +
        ["--band-labels"] + klabels +
        ["--band-points", "101", "--band-connection", "band.conf"]
    )

    result = _run_phonopy(band_args, timeout=120)
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

    # --- Step 4: Plot ---
    do_plot = True
    if parsed_args and parsed_args.no_plot:
        do_plot = False
    elif interactive:
        answer = _prompt_choice(console, "Generate phonon band plot?", ["Yes", "No"], "Yes")
        do_plot = "Yes" in answer

    if do_plot:
        console.print()
        console.print("[bold]Step 4: Plot phonon dispersion[/bold]")
        try:
            _plot_phonon_bands(console)
            # Combined three-panel figure if DOS data exists
            tdos_path = Path("total_dos.dat")
            pdos_path = Path("partial_dos.dat")
            if tdos_path.exists():
                _plot_phonon_combined(console, tdos_path, pdos_path)
        except Exception as e:
            console.print(f"[red]Plot failed: {e}[/red]")

    console.print()


def _parse_band_yaml(band_yaml: Path):
    """Parse band.yaml with segment detection by distance breaks.

    Returns (xs, freqs, tick_pos, tick_labels) where:
      xs[:]      — continuous x-axis (distance)
      freqs[:,:] — shape (n_qpts, n_bands)
      tick_pos   — list of x positions for segment boundaries
      tick_labels — corresponding high-symmetry labels
    """
    import yaml

    data = yaml.safe_load(band_yaml.read_text(encoding="utf-8"))
    phonon = data.get("phonon", [])
    if not phonon:
        raise ValueError("band.yaml has no phonon data")

    distances = np.array([q["distance"] for q in phonon], dtype=float)
    freqs_raw = np.array([[b["frequency"] for b in q["band"]] for q in phonon], dtype=float)

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

    # Tick positions and labels from per-segment start/end
    tick_pos = []
    tick_labels_raw: list[list[str]] = []
    for k in range(len(breaks) - 1):
        s, e = breaks[k], breaks[k + 1]
        sl = phonon[s].get("label", "").strip()
        el = phonon[e - 1].get("label", "").strip()
        tick_labels_raw.append([sl, el])

    tick_labels: list[str] = []
    for k, (seg_xs, (s_lbl, e_lbl)) in enumerate(zip(xs_list, tick_labels_raw)):
        if k == 0:
            tick_pos.append(seg_xs[0])
            tick_labels.append(_fmt_label(s_lbl))
        else:
            prev = tick_labels[-1]
            cur = _fmt_label(s_lbl)
            if prev != cur and cur:
                tick_labels[-1] = f"{prev}|{cur}"
            elif cur and not prev:
                tick_labels[-1] = cur
        tick_pos.append(seg_xs[-1])
        tick_labels.append(_fmt_label(e_lbl))

    return xs_all, freq_all, tick_pos, tick_labels


def _fmt_label(lbl: str) -> str:
    if not lbl:
        return ""
    u = lbl.upper()
    if u in ("GAMMA", "GM", "G", "\\GAMMA"):
        return r"$\Gamma$"
    if lbl == r"\Gamma":
        return r"$\Gamma$"
    return f"${lbl}$"


def _parse_total_dos(path: Path):
    """Return (freq, dos) from total_dos.dat."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def _parse_partial_dos(path: Path):
    """Return (freq, pdos_per_atom) from partial_dos.dat."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1:]


def _plot_phonon_bands(console) -> None:
    """Parse band.yaml and plot phonon dispersion."""
    xs, freqs, tick_pos, tick_labels = _parse_band_yaml(Path("band.yaml"))

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
    ax.set_ylim(bottom=0)
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


def _plot_phonon_combined(console, tdos_path: Path, pdos_path: Path) -> None:
    """Combined three-panel figure: band + TDOS + PDOS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    xs, freqs, tick_pos, tick_labels = _parse_band_yaml(Path("band.yaml"))
    freq_t, dos_t = _parse_total_dos(tdos_path)

    has_pdos = pdos_path.exists()
    if has_pdos:
        freq_p, pdos_data = _parse_partial_dos(pdos_path)
        # Determine species from STRU
        species_list: list[str] = []
        try:
            from abacuscopilot.io.stru_file import read_stru
            s = read_stru("STRU")
            species_list = list(s.species_order)
        except Exception:
            species_list = [f"atom {i + 1}" for i in range(pdos_data.shape[1])]

    f_lo = 0.0
    f_hi = float(np.max(freqs)) * 1.05

    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[3, 1.2, 1.2], wspace=0.08)
    ax_band = fig.add_subplot(gs[0])
    ax_tdos = fig.add_subplot(gs[1])
    ax_pdos = fig.add_subplot(gs[2])

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
    ax_band.yaxis.grid(True, ls=":", alpha=0.4)

    # ── Total DOS ──
    ax_tdos.fill_betweenx(freq_t, 0, dos_t, color="#2c7bb6", alpha=0.35, lw=0)
    ax_tdos.plot(dos_t, freq_t, color="#2c7bb6", lw=1.2)
    ax_tdos.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.5)
    ax_tdos.set_ylim(f_lo, f_hi)
    ax_tdos.set_xlabel("DOS", fontsize=11)
    ax_tdos.set_title("Total DOS", fontsize=13, pad=8)
    ax_tdos.set_yticklabels([])
    ax_tdos.yaxis.grid(True, ls=":", alpha=0.4)

    # ── PDOS ──
    pdos_colors = ["#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                   "#F4A261", "#264653", "#A8DADC", "#6A0572"]
    if has_pdos:
        max_val = max(np.max(pdos_data[:, i]) for i in range(pdos_data.shape[1]))
        max_val = max_val if max_val > 0 else 1.0
        for i in range(pdos_data.shape[1]):
            label = species_list[i] if i < len(species_list) else f"atom {i + 1}"
            c = pdos_colors[i % len(pdos_colors)]
            dos_norm = pdos_data[:, i] / max_val
            ax_pdos.fill_betweenx(freq_p, 0, dos_norm, color=c, alpha=0.30, lw=0)
            ax_pdos.plot(dos_norm, freq_p, color=c, lw=1.3, label=label)
        ax_pdos.legend(loc="upper right", fontsize=9, framealpha=0.7)
    else:
        ax_pdos.text(0.5, 0.5, "PDOS not available", ha="center", va="center",
                     transform=ax_pdos.transAxes, fontsize=11, color="gray")
    ax_pdos.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.5)
    ax_pdos.set_ylim(f_lo, f_hi)
    ax_pdos.set_xlabel("pDOS (norm.)", fontsize=11)
    ax_pdos.set_title("Projected DOS", fontsize=13, pad=8)
    ax_pdos.set_yticklabels([])
    ax_pdos.yaxis.grid(True, ls=":", alpha=0.4)

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

    if interactive:
        mesh_s = _prompt(console, "q-point mesh for DOS", "8 8 8")
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
        phonon = load(
            supercell_matrix=[int(x) for x in dim_default.split()],
            primitive_matrix="auto",
            unitcell_filename="STRU",
            force_sets_filename="FORCE_SETS",
            calculator="abacus",
            is_nac=False,
            produce_fc=True,
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
            ax.set_xlim(left=0)
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

    if interactive:
        mesh_s = _prompt(console, "q-point mesh for PDOS", "8 8 8")
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
        phonon = load(
            supercell_matrix=[int(x) for x in dim_default.split()],
            primitive_matrix="auto",
            unitcell_filename="STRU",
            force_sets_filename="FORCE_SETS",
            calculator="abacus",
            is_nac=False,
            produce_fc=True,
        )
        phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        phonon.run_projected_dos()

        pd = phonon._pdos
        freq = pd.frequency_points
        pdos_data = pd.projected_dos  # shape (n_atoms_prim, n_freq)
        # Sum PDOS over all atoms
        total_pdos = pdos_data.sum(axis=0)

        # Save: freq + per-atom PDOS + total
        cols = [freq] + [pdos_data[i, :] for i in range(pdos_data.shape[0])] + [total_pdos]
        header = "frequency(THz)  " + "  ".join(
            [f"atom_{i + 1}" for i in range(pdos_data.shape[0])] + ["total"]
        )
        np.savetxt("partial_dos.dat", np.column_stack(cols), header=header)
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
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from abacuscopilot.plotting.style import load_style_from_config
            load_style_from_config()

            fig, ax = plt.subplots(figsize=(7, 5))
            pdos_colors = ["#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                           "#F4A261", "#264653", "#A8DADC", "#6A0572"]
            for i in range(pdos_data.shape[0]):
                if i < len(species_list):
                    label = species_list[i]
                else:
                    label = f"atom {i + 1}"
                if selected_labels and i < len(species_list) and species_list[i] not in selected_labels:
                    continue
                c = pdos_colors[i % len(pdos_colors)]
                ax.plot(freq, pdos_data[i, :], color=c, lw=1.3, label=label)
                ax.fill_between(freq, 0, pdos_data[i, :], color=c, alpha=0.12)

            ax.axvline(0, color="gray", lw=0.7, ls="--", alpha=0.6)
            ax.set_xlim(left=0)
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
