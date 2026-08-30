"""Molecular dynamics postprocessing tasks for ABACUS.

Task IDs 781-799

Reads MD_dump trajectory files from ABACUS MD output and provides:
- PDB export for VMD visualisation
- Frame extraction
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.tasks import task

# =============================================================================
# MD_dump parser
# =============================================================================


def parse_md_dump(filepath: str | Path) -> list[dict]:
    """Parse an ABACUS MD_dump trajectory file.

    Each frame starts with ``MDSTEP: N``, followed by lattice info and
    a table of atomic positions / forces / velocities.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"MD_dump not found: {filepath}")

    frames: list[dict] = []
    current_frame: dict | None = None
    in_table = False
    lat_vecs: list[list[float]] | None = None  # None = not collecting, [] = collecting
    # NVT MD_dump writes lattice info only on the first frame (volume is
    # constant).  NPT writes it on every frame.  Carry the last-seen lattice
    # forward so NVT frames also have it (needed for POSCAR export etc.).
    last_lc: float | None = None
    last_lv: np.ndarray | None = None

    with open(filepath, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("MDSTEP:"):
                if current_frame is not None and current_frame.get("atoms"):
                    if "lattice_vectors" not in current_frame and last_lv is not None:
                        current_frame["lattice_constant"] = last_lc
                        current_frame["lattice_vectors"] = last_lv
                    frames.append(current_frame)
                # Fast int parse: "MDSTEP:  123" → skip 7 chars, strip, int
                step_str = s[7:].strip()
                current_frame = {"step": int(step_str), "atoms": []}
                in_table = False
                lat_vecs = None
                continue

            if current_frame is None:
                continue

            if s.startswith("LATTICE_CONSTANT:"):
                # "LATTICE_CONSTANT: 1.0 Angstrom" → extract float
                lc_str = s.split(":")[1].strip().split()[0]
                last_lc = current_frame["lattice_constant"] = float(lc_str)
                continue

            if s == "LATTICE_VECTORS":
                lat_vecs = []
                continue

            # Still collecting lattice vectors (list, not yet converted to array)
            if isinstance(current_frame.get("lattice_vectors"), list) or lat_vecs is not None:
                if s.startswith("INDEX") or s.startswith("MDSTEP"):
                    if len(lat_vecs) == 3:
                        current_frame["lattice_vectors"] = np.array(lat_vecs)
                    lat_vecs = None
                    continue
                parts = s.split()
                if len(parts) >= 3:
                    try:
                        lat_vecs.append([float(x) for x in parts[:3]])
                    except ValueError:
                        if len(lat_vecs) == 3:
                            current_frame["lattice_vectors"] = np.array(lat_vecs)
                        lat_vecs = None
                        continue
                    if len(lat_vecs) == 3:
                        current_frame["lattice_vectors"] = np.array(lat_vecs)
                        last_lv = current_frame["lattice_vectors"]
                        lat_vecs = None  # done collecting
                    continue

            if s.startswith("INDEX"):
                in_table = True
                continue

            if in_table:
                parts = s.split()
                if len(parts) < 10:
                    in_table = False
                    continue
                try:
                    current_frame["atoms"].append({
                        "index": int(parts[0]),
                        "label": parts[1],
                        "xyz": np.array([float(parts[2]), float(parts[3]), float(parts[4])]),
                        "force": np.array([float(parts[5]), float(parts[6]), float(parts[7])]),
                        "velocity": np.array([float(parts[8]), float(parts[9]), float(parts[10])]),
                    })
                except (ValueError, IndexError):
                    in_table = False

    if current_frame is not None and current_frame.get("atoms"):
        if "lattice_vectors" not in current_frame and last_lv is not None:
            current_frame["lattice_constant"] = last_lc
            current_frame["lattice_vectors"] = last_lv
        frames.append(current_frame)

    return frames


# =============================================================================
# PDB writer
# =============================================================================

_PDB_ATOM_FMT = (
    "ATOM  {serial:>5d} {name:>4s} {res:>3s} {chain:1s}{res_seq:>4d}    "
    "{x:8.3f}{y:8.3f}{z:8.3f}  {occ:4.2f}  {bfactor:5.2f}           {elem:>2s}\n"
)


def write_pdb(
    frame: dict,
    filepath: str | Path,
    *,
    model_num: int | None = None,
) -> None:
    """Write a single MD frame as a PDB file.

    Args:
        frame: Frame dict from ``parse_md_dump``.
        filepath: Output .pdb path.
        model_num: Optional MODEL record number (for multi-model PDBs).
    """
    filepath = Path(filepath)
    mode = "a" if model_num is not None and model_num > 1 else "w"

    with open(filepath, mode) as f:
        if model_num is not None:
            f.write(f"MODEL     {model_num:>4d}\n")

        # CRYST1 record (lattice info — optional but helps VMD)
        lc = frame.get("lattice_constant", 1.0)
        lv = frame.get("lattice_vectors")
        if lv is not None and lv.shape == (3, 3):
            cell = lc * np.array(lv)
            a = np.linalg.norm(cell[0])
            b = np.linalg.norm(cell[1])
            c = np.linalg.norm(cell[2])
            # Compute angles
            alpha = np.degrees(
                np.arccos(np.dot(cell[1], cell[2]) / (b * c))
            )
            beta = np.degrees(
                np.arccos(np.dot(cell[0], cell[2]) / (a * c))
            )
            gamma = np.degrees(
                np.arccos(np.dot(cell[0], cell[1]) / (a * b))
            )
            f.write(
                f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
                f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n"
            )

        for atom in frame["atoms"]:
            elem = atom["label"]
            x, y, z = atom["xyz"]
            f.write(
                _PDB_ATOM_FMT.format(
                    serial=atom["index"] + 1,
                    name=elem,
                    res=elem,
                    chain="A",
                    res_seq=1,
                    x=x,
                    y=y,
                    z=z,
                    occ=1.00,
                    bfactor=0.00,
                    elem=elem,
                )
            )

        if model_num is not None:
            f.write("ENDMDL\n")


def _find_md_dump(args: list[str] | None = None) -> str:
    """Auto-locate the MD_dump trajectory file.

    Search order: OUT.ABACUS/MD_dump → MD_dump → OUT.*/MD_dump.
    CLI arguments override the auto-detected path.
    """
    md_path = "OUT.ABACUS/MD_dump"
    for candidate in ("OUT.ABACUS/MD_dump", "MD_dump"):
        if Path(candidate).exists():
            md_path = candidate
            break
    else:
        for d in sorted(Path().glob("OUT.*")):
            p = d / "MD_dump"
            if p.exists():
                md_path = str(p)
                break
    if args:
        for arg in args:
            p = Path(arg)
            if p.exists() and p.is_file():
                md_path = str(p)
                break
    return md_path


# =============================================================================
# Task 781: MD_dump → PDB
# =============================================================================


@task(3101, category="MD Analysis", name="MD Trajectory → PDB",
      description="Convert ABACUS MD_dump trajectory to PDB format for VMD")
def task_md_to_pdb(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert ABACUS MD_dump trajectory file to PDB format.

    PDB files can be loaded directly into VMD/Chimera/PyMOL for
    visualisation and analysis of the MD trajectory.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== MD Trajectory → PDB ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    # Parse
    console.print(f"  [dim]Reading: {md_path}[/dim]")
    try:
        frames = parse_md_dump(md_path)
    except Exception as e:
        console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
        return

    if not frames:
        console.print("[red]No frames found in MD_dump.[/red]")
        return

    n_frames = len(frames)
    n_atoms = len(frames[0]["atoms"])
    console.print(f"  Frames: {n_frames}, Atoms per frame: {n_atoms}")
    console.print(f"  Steps: {frames[0]['step']} → {frames[-1]['step']}")
    console.print()

    # Stride (frame sampling)
    if interactive:
        stride = int(_prompt(console, "Frame stride (1 = all frames)", "1"))
    else:
        stride = 1
    stride = max(1, stride)
    if stride > 1:
        frames = frames[::stride]
        n_frames = len(frames)
        console.print(f"  [dim]After stride {stride}: {n_frames} frames[/dim]")
        console.print()

    # Choose output
    if interactive:
        mode = _prompt_choice(
            console,
            "Output mode",
            [
                "Single PDB (all frames in one file, MODEL/ENDMDL)",
                "First frame only",
                "Last frame only",
                "Range of frames (start-end)",
            ],
            "Single PDB (all frames in one file, MODEL/ENDMDL)",
        )
    else:
        mode = "Single PDB (all frames in one file, MODEL/ENDMDL)"

    # Determine output path
    out_path = Path(md_path).with_suffix(".pdb")

    if "all frames" in mode:
        # Multi-model PDB — keep original MD step numbers
        from rich.progress import Progress
        with Progress() as progress:
            task = progress.add_task("[cyan]Writing PDB...", total=len(frames))
            for frame in frames:
                write_pdb(frame, out_path, model_num=frame["step"])
                progress.update(task, advance=1)
        console.print(f"  [green]✓ {n_frames} frames → {out_path}[/green]")

    elif "First" in mode:
        write_pdb(frames[0], out_path)
        console.print(f"  [green]✓ Frame {frames[0]['step']} → {out_path}[/green]")

    elif "Last" in mode:
        write_pdb(frames[-1], out_path)
        console.print(f"  [green]✓ Frame {frames[-1]['step']} → {out_path}[/green]")

    elif "Range" in mode:
        rng = _prompt(console, f"  Frame range (1-{n_frames})", f"1-{min(10, n_frames)}")
        try:
            parts = [int(x.strip()) for x in rng.split("-")]
            start, end = max(1, parts[0]), min(n_frames, parts[-1])
        except (ValueError, IndexError):
            start, end = 1, min(10, n_frames)
        for i in range(start - 1, end):
            write_pdb(frames[i], out_path, model_num=frames[i]["step"])
        console.print(f"  [green]✓ Frames {start}→{end} → {out_path}[/green]")

    console.print(f"  [dim]Open with: vmd {out_path}[/dim]")
    console.print()


# =============================================================================
# MD_dump writer
# =============================================================================


def write_md_dump(frames: list[dict], filepath: str | Path) -> None:
    """Write selected frames to an ABACUS MD_dump file.

    Produces a valid MD_dump that can be re-read by ``parse_md_dump``
    or visualised after PDB conversion.
    """
    filepath = Path(filepath)
    with open(filepath, "w") as f:
        for frame in frames:
            f.write(f"MDSTEP:  {frame['step']}\n")
            lc = frame.get("lattice_constant", 1.0)
            f.write(f"LATTICE_CONSTANT: {lc:.12f} Angstrom\n")
            f.write("LATTICE_VECTORS\n")
            lv = frame.get("lattice_vectors")
            if lv is not None:
                for row in lv:
                    f.write(f"  {row[0]:.12f}  {row[1]:.12f}  {row[2]:.12f}\n")
            f.write(
                "INDEX    LABEL    POSITION (Angstrom)    "
                "FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)\n"
            )
            for atom in frame["atoms"]:
                x, y, z = atom["xyz"]
                fx, fy, fz = atom["force"]
                vx, vy, vz = atom["velocity"]
                f.write(
                    f"  {atom['index']}  {atom['label']}"
                    f"  {x:.12f}  {y:.12f}  {z:.12f}"
                    f"  {fx:.12f}  {fy:.12f}  {fz:.12f}"
                    f"  {vx:.12f}  {vy:.12f}  {vz:.12f}\n"
                )


# =============================================================================
# Task 783: Extract frames (stride)
# =============================================================================


@task(3102, category="MD Analysis", name="Extract Frames",
      description="Extract every Nth frame from MD_dump to create a lighter trajectory")
def task_extract_frames(args: list[str] | None = None, interactive: bool = True) -> None:
    """Extract frames at a regular stride, writing a smaller MD_dump."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Extract Frames (Stride) ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    console.print(f"  [dim]Reading: {md_path}[/dim]")
    try:
        frames = parse_md_dump(md_path)
    except Exception as e:
        console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
        return

    n_frames = len(frames)
    console.print(f"  Frames: {n_frames}")
    console.print()

    # Stride
    if interactive:
        stride = int(_prompt(console, "Extract every N frames (stride)", "10"))
    else:
        stride = 10
    stride = max(1, stride)

    selected = frames[::stride]
    out_path = Path(md_path).parent / "MD_sampled.dump"
    write_md_dump(selected, out_path)

    console.print(f"  [green]✓ {len(selected)} frames → {out_path}[/green]")
    console.print(f"  Stride: {stride}, original: {n_frames} → reduced: {len(selected)}")
    console.print()


# =============================================================================
# MSD calculation
# =============================================================================


def _read_md_dt(md_path: str) -> float:
    """Read md_dt (fs) from INPUT or running_md.log.

    Looks in order:
    1. INPUT file next to MD_dump
    2. OUT.ABACUS/INPUT.info (DP-MD convention)
    3. OUT.ABACUS/running_md.log (fallback — md_dt is logged there)
    """
    md_dir = Path(md_path).parent
    for candidate in (md_dir / "INPUT", md_dir / "INPUT.info",
                      md_dir.parent / "INPUT", md_dir / "running_md.log"):
        if not candidate.exists():
            continue
        content = candidate.read_text()
        m = re.search(r"md_dt\s+([\d.]+)", content)
        if m:
            return float(m.group(1))
    # Also try OUT.ABACUS/running_md.log if md_dir is not already OUT.ABACUS
    out_log = md_dir / "OUT.ABACUS" / "running_md.log"
    if out_log.exists() and out_log != md_dir / "running_md.log":
        content = out_log.read_text()
        m = re.search(r"md_dt\s+([\d.]+)", content)
        if m:
            return float(m.group(1))
    raise FileNotFoundError(f"md_dt not found — no INPUT or running_md.log near {md_path}")


def _read_md_dumpfreq(md_path: str) -> int:
    """Read md_dumpfreq from INPUT — the trajectory output stride in steps.

    MD_dump frames are written every *md_dumpfreq* MD steps, so the real
    time between two consecutive trajectory frames is:

        frame_dt (fs) = md_dt * md_dumpfreq
    """
    md_dir = Path(md_path).parent
    for candidate in (md_dir / "INPUT", md_dir / "INPUT.info",
                      md_dir.parent / "INPUT", md_dir / "running_md.log"):
        if not candidate.exists():
            continue
        content = candidate.read_text()
        m = re.search(r"md_dumpfreq\s+(\d+)", content)
        if m:
            return int(m.group(1))
    # Try OUT.ABACUS/running_md.log
    out_log = md_dir / "OUT.ABACUS" / "running_md.log"
    if out_log.exists():
        content = out_log.read_text()
        m = re.search(r"md_dumpfreq\s+(\d+)", content)
        if m:
            return int(m.group(1))
    return 1  # default: every step dumped


def _resolve_frame_dt(md_path: str | Path, console, interactive: bool = True) -> float:
    """Real time between consecutive trajectory frames, in fs.

    For ABACUS-native MD it is md_dt × md_dumpfreq (read from INPUT /
    running_md.log).  For trajectories converted from VASP XDATCAR there is no
    ABACUS INPUT, so in interactive mode we ask the user for the frame spacing
    (e.g. VASP POTIM); non-interactively we warn and assume 1.0 fs.
    """
    try:
        md_dt = _read_md_dt(str(md_path))
        dumpfreq = _read_md_dumpfreq(str(md_path))
        frame_dt = md_dt * dumpfreq
        console.print(f"  [dim]md_dt = {md_dt:.1f} fs, dumpfreq = {dumpfreq} → frame_dt = {frame_dt:.1f} fs[/dim]")
        return frame_dt
    except FileNotFoundError:
        if interactive:
            frame_dt = float(_prompt(
                console,
                "No ABACUS INPUT found — real time between trajectory frames (fs)",
                "1.0",
            ))
            console.print(f"  [dim]frame_dt = {frame_dt:.1f} fs (manual)[/dim]")
            return frame_dt
        console.print("  [yellow]! No ABACUS INPUT found — assuming frame_dt = 1.0 fs.[/yellow]")
        console.print("  [yellow]  If this trajectory was converted from VASP XDATCAR, the real frame")
        console.print("  [yellow]  spacing may differ — re-run interactively to set it.[/yellow]")
        return 1.0


def _prompt_range(console, label: str, prompt_text: str) -> tuple[bool, tuple | None]:
    """Ask for a 'min,max' range interactively.

    Returns (False, None) when left empty (keep current), or (True, (lo, hi))
    for a parsed range.  Re-asks on unparseable input.
    """
    while True:
        val = _prompt(console, prompt_text, "")
        if not val.strip():
            return False, None
        try:
            lo, hi = (float(x) for x in val.replace(" ", "").split(","))
            return True, (lo, hi)
        except (ValueError, TypeError):
            console.print(
                f"  [yellow]! Could not parse {label} range — use 'min,max' (e.g. '0,12').[/yellow]"
            )


def _plot_msd(data_file: str, species: str, dt_fs: float,
              console=None, interactive: bool = True) -> None:
    """Plot MSD curve(s).

    Interactive: after each plot asks "Adjust plot ranges?"; choosing Yes lets
    the user set the X and Y axis ranges (empty keeps the current value) and
    redraws immediately.  Choosing No keeps the plot and returns.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    data = np.loadtxt(data_file)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    time_frames = data[:, 0]
    n_curves = data.shape[1] - 1

    # --- Colour map ---
    _COLORS: dict[str, str] = {
        "total": "#cb5c55",  # warm coral
        "x":     "#0060a6",  # deep blue
        "y":     "#969fb5",  # muted purple-grey
        "z":     "#4b8fc1",  # medium blue
    }
    _FALLBACK = ["#cb5c55", "#0060a6", "#969fb5", "#4b8fc1"]

    time_ps = time_frames * dt_fs / 1000.0

    comp_names = {1: ["total"], 3: ["x", "y", "z"], 4: ["total", "x", "y", "z"]}
    names = comp_names.get(n_curves, [str(i) for i in range(n_curves)])

    ylim = None
    xlim = None
    while True:
        fig, ax = plt.subplots(figsize=(8, 6))
        for ic in range(n_curves):
            msd_vals = data[:, ic + 1]
            label = names[ic]
            c = _COLORS.get(label, _FALLBACK[ic % len(_FALLBACK)])
            ax.plot(time_ps, msd_vals, "-", color=c, linewidth=1.2, label=label)

        # --- Labels (journal style: no in-figure title) ---
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("MSD (Å²)")
        ax.legend(frameon=False, handlelength=2.0, handletextpad=0.6)

        # --- Spine styling ---
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        ax.tick_params(axis="both", direction="out")

        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)

        fig.tight_layout(pad=1.2)
        save_name = f"MSD_{species}.png"
        fig.savefig(save_name, dpi=300, bbox_inches="tight")
        plt.close(fig)
        if console is not None:
            console.print(f"  [green]✓ MSD plot: {save_name}[/green]")

        if not interactive:
            break

        # Explicit gate after each plot: keep it, or adjust X/Y ranges.
        adjust = _prompt_choice(
            console,
            "Adjust plot ranges?",
            ["No — keep this plot", "Yes — adjust X/Y ranges"],
            "No — keep this plot",
        )
        if adjust.startswith("No"):
            break

        # Ask X and Y; empty keeps the current value.  Loop back to redraw.
        ch, val = _prompt_range(console, "Y", "Y-axis range (e.g. '0,12', empty=keep): ")
        if ch:
            ylim = val
        ch, val = _prompt_range(console, "X", "X-axis range (e.g. '0,50', empty=keep): ")
        if ch:
            xlim = val


def _compute_msd_fft(
    traj: np.ndarray,
    max_lag: int,
    axes: dict[str, list[int]],
) -> dict[str, np.ndarray]:
    """MSD via FFT autocorrelation — O(N log N) instead of the naive O(N²).

    Uses the identity
        MSD(τ) = [ S_r(τ) + S_f(N−τ) − 2·C(τ) ] / (N−τ) / n_atoms
    where S_r/S_f are cumulative sums of |r|² (O(N)) and C(τ) is the
    dot-product autocorrelation Σ_t r(t+τ)·r(t), computed with the FFT
    (Wiener–Khinchin) in O(N log N).  Same result as the naive pair loop
    (verified to ~1e-13), returns the per-atom mean MSD.

    Args:
        traj: (n_frames, n_atoms, 3) coordinates in Å.
        max_lag: Number of lags to compute (≤ n_frames − 1).
        axes: {label: [component indices]}, e.g. {"xyz": [0, 1, 2]}.

    Returns:
        {label: msd[0..max_lag]}.
    """
    n_frames, n_atoms, _ = traj.shape
    n_fft = 1 << (2 * n_frames - 1).bit_length()  # ≥ 2N → linear autocorr

    results: dict[str, np.ndarray] = {}
    for label, comps in axes.items():
        # Dot-product autocorrelation summed over atoms & components, via FFT.
        x = traj[:, :, comps].reshape(n_frames, -1)  # (N, n_atoms*ncomp)
        X = np.fft.rfft(x, n=n_fft, axis=0)
        C = np.fft.irfft((X.real ** 2 + X.imag ** 2).sum(axis=1), n=n_fft)[:n_frames]

        # Cumulative sums of |r|² over time (per atom).
        sum2 = (traj[:, :, comps] ** 2).sum(axis=2)                 # (N, n_atoms)
        S_f = np.zeros((n_frames + 1, n_atoms))
        np.cumsum(sum2, axis=0, out=S_f[1:])                        # S_f[k] = Σ_{t<k}
        S_r = np.flip(np.cumsum(np.flip(sum2, axis=0), axis=0), axis=0)  # S_r[k] = Σ_{t≥k}

        msd = np.zeros(max_lag + 1)
        for lag in range(1, max_lag + 1):
            num = S_r[lag].sum() + S_f[n_frames - lag].sum() - 2.0 * C[lag]
            msd[lag] = num / (n_atoms * (n_frames - lag))
        results[label] = msd

    return results


def compute_msd(
    frames: list[dict],
    atom_indices: list[int],
    directions: str = "xyz",
    *,
    unwrap: bool = True,
) -> dict[str, np.ndarray]:
    """Compute mean square displacement (MSD) for selected atoms.

    MSD_d(t) = <|r_d(t + τ) − r_d(τ)|²>

    where *d* is one of x, y, z or the total magnitude, and <...> is the mean
    over starting times *and* over the selected atoms (per-atom mean MSD).

    Uses the FFT autocorrelation (O(N log N)) — same numbers as the naive
    pair loop but orders of magnitude faster on long trajectories.

    Args:
        frames: List of frame dicts from ``parse_md_dump``.
        atom_indices: Indices of atoms to track (same across all frames).
        directions: Which directions to compute:
            ``"x"`` / ``"y"`` / ``"z"`` — single component,
            ``"xyz"`` — total (sum over all 3),
            ``"all"`` — x, y, z separately + total.
        unwrap: If True (default), apply minimum-image PBC unwrapping
            via ``_unwrap_trajectory`` before computing MSD.  This prevents
            artificial jumps when atoms cross periodic boundaries.
            Set to False only for debugging / non-periodic systems.

    Returns:
        Dict mapping label → 1-D MSD array.  Also includes a ``"lag"`` key
        with the frame-lag axis.
    """
    n_frames = len(frames)
    n_atoms = len(atom_indices)

    # Build trajectory array: (n_frames, n_atoms, 3) — vectorized extraction
    traj = np.array([[frame["atoms"][i]["xyz"] for i in atom_indices]
                      for frame in frames])

    # ------------------------------------------------------------------
    # PBC unwrapping — remove artificial jumps across periodic boundaries
    # before computing displacements.  Without this, an atom crossing
    # from x≈L to x≈0 registers a jump of ~−L, inflating MSD.
    # ------------------------------------------------------------------
    if unwrap:
        unwrapped = _unwrap_trajectory(frames, atom_indices)  # (n_atoms, 3, n_frames)
        traj = unwrapped.transpose(2, 0, 1)                   # → (n_frames, n_atoms, 3)

    max_lag = n_frames - 1
    result = {"lag": np.arange(max_lag + 1, dtype=float)}

    # Which axes to compute
    axes: dict[str, list[int]] = {}
    if directions == "all":
        axes = {"total": [0, 1, 2], "x": [0], "y": [1], "z": [2]}
    elif directions in ("x", "y", "z"):
        idx = {"x": 0, "y": 1, "z": 2}[directions]
        axes = {directions: [idx]}
    else:  # "xyz" (default)
        axes = {"xyz": [0, 1, 2]}

    # FFT-based MSD — O(N log N) instead of O(N²). Same result as the naive
    # pair loop (verified to ~1e-13), returns the per-atom mean MSD.
    console = _get_console()
    console.print("  [dim]Using FFT autocorrelation (O(N log N))...[/dim]")
    result.update(_compute_msd_fft(traj, max_lag, axes))
    return result


# =============================================================================
# Task 782: MSD
# =============================================================================


@task(3103, category="MD Analysis", name="MSD",
      description="Compute Mean Square Displacement from MD_dump trajectory")
def task_md_msd(args: list[str] | None = None, interactive: bool = True) -> None:
    """Calculate mean square displacement (MSD) for selected atomic species.

    Writes ``MSD_{species}.dat`` with columns:  frame_lag  MSD_Å²
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Mean Square Displacement (MSD) ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    console.print(f"  [dim]Reading: {md_path}[/dim]")
    try:
        frames = parse_md_dump(md_path)
    except Exception as e:
        console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
        return

    if not frames:
        console.print("[red]No frames found.[/red]")
        return

    n_frames = len(frames)
    n_atoms = len(frames[0]["atoms"])

    # Collect species in trajectory order (not alphabetical)
    species_seen = []
    for a in frames[0]["atoms"]:
        if a["label"] not in species_seen:
            species_seen.append(a["label"])
    species_counts = {s: sum(1 for a in frames[0]["atoms"] if a["label"] == s)
                      for s in species_seen}

    console.print(f"  Frames: {n_frames}, Atoms: {n_atoms}")
    console.print(f"  Species: {', '.join(f'{s}({species_counts[s]})' for s in species_seen)}")
    console.print(f"  Steps: {frames[0]['step']} → {frames[-1]['step']}")
    console.print()

    # --- Choose species ---
    if interactive:
        species_choices = [
            f"{s} ({species_counts[s]} atoms)" for s in species_seen
        ]
        chosen = _prompt_choice(console, "Select species for MSD", species_choices,
                                species_choices[0])
        # Extract element symbol from "Li (24 atoms)"
        species = chosen.split()[0]
    else:
        species = species_seen[0]

    # --- Frame range ---
    last_frame = n_frames - 1
    if interactive:
        start = int(_prompt(console, "Start frame", "500"))
        end = int(_prompt(console, "End frame", str(last_frame)))
        step = int(_prompt(console, "Frame interval", "1"))
    else:
        start, end, step = 500, last_frame, 1

    # Clamp and validate
    start = max(0, start)
    end = min(last_frame, end)
    step = max(1, step)
    if end <= start:
        console.print(f"[red]End frame ({end}) must be greater than start ({start}).[/red]")
        return

    selected_frames = frames[start:end:step]
    n_selected = len(selected_frames)
    console.print(f"\n  Selected frames: {start} → {end} (step {step}) = {n_selected} frames")
    console.print()

    # Get atom indices for chosen species
    atom_indices = [
        i for i, a in enumerate(selected_frames[0]["atoms"])
        if a["label"] == species
    ]
    console.print(f"  {len(atom_indices)} {species} atoms selected for MSD")
    console.print()

    # --- Choose direction ---
    if interactive:
        dir_choice = _prompt_choice(
            console,
            "MSD direction",
            ["xyz (total)", "x only", "y only", "z only", "x, y, z + total"],
            "xyz (total)",
        )
    else:
        dir_choice = "xyz (total)"

    if "x, y, z" in dir_choice:
        directions = "all"
    elif "x only" in dir_choice:
        directions = "x"
    elif "y only" in dir_choice:
        directions = "y"
    elif "z only" in dir_choice:
        directions = "z"
    else:
        directions = "xyz"

    # Compute MSD
    console.print("  Computing MSD...")
    result = compute_msd(selected_frames, atom_indices, directions=directions)
    lag = result.pop("lag")

    # Build data columns: lag + each MSD column
    cols = [lag]
    col_labels = []
    for label, msd_arr in result.items():
        cols.append(msd_arr)
        col_labels.append(label)

    # Write output
    out_file = f"MSD_{species}.dat"
    header = f"# MSD for {species} — frames {start}→{end} step {step}\n"
    header += "# per-atom mean MSD (Å²); D = slope/(2·d·frame_dt) with d = dims\n"
    header += "# lag(frames)  " + "  ".join(f"MSD_{lab}(A^2)" for lab in col_labels)
    np.savetxt(out_file, np.column_stack(cols),
               fmt="%8d  " + "  ".join("%.8f" for _ in col_labels),
               header=header)
    console.print(f"  [green]✓ MSD saved to {out_file}[/green]")

    # Quick summary
    for lab in col_labels:
        arr = result[lab]
        if len(arr) > 50:
            console.print(f"  MSD_{lab}(10) = {arr[10]:.4f} Å²,  MSD_{lab}(50) = {arr[50]:.4f} Å²")

    # Auto-plot and compute diffusion coefficient
    frame_dt = _resolve_frame_dt(md_path, console, interactive)
    _plot_msd(out_file, species, frame_dt, console=console, interactive=interactive)
    console.print(f"  [green]✓ MSD plot saved to MSD_{species}.png[/green]")
    for lab in col_labels:
        arr = result[lab]
        if len(arr) > 1:
            n = len(arr)
            fit_start = max(1, n // 10)
            fit_end = int(n * 0.8)
            slope, _ = np.polyfit(lag[fit_start:fit_end], arr[fit_start:fit_end], 1)
            # Convert to Å²/fs then to cm²/s
            # D_total = slope / (2d * frame_dt), d=3 → 1/(6*frame_dt); per-direction → 1/(2*frame_dt)
            diff_cm2_s = abs(slope) / (6.0 * frame_dt) * 1e-1 if lab == "total" else abs(slope) / (2.0 * frame_dt) * 1e-1
            console.print(f"  D_{lab}({species}) = {diff_cm2_s:.4e} cm²/s")
    console.print()


# =============================================================================
# RDF calculation
# =============================================================================


def _minimum_image(dr: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Apply minimum-image convention to a displacement vector.

    Args:
        dr: (3,) or (N,3) displacement vectors (Angstrom).
        cell: (3,3) cell vectors (Angstrom).

    Returns:
        Corrected displacement vectors (Angstrom).
    """
    # Convert to fractional, wrap to [-0.5, 0.5), convert back
    cell_inv = np.linalg.inv(cell.T)
    frac = dr @ cell_inv
    frac -= np.floor(frac + 0.5)
    return frac @ cell.T


def compute_rdf(
    frames: list[dict],
    species_a: str,
    species_b: str,
    r_max: float = 10.0,
    n_bins: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute partial radial distribution function g_{ab}(r).

    g(r) = (V / (4π r² Δr N_a N_b)) Σ δ(r - r_ij)

    Averages over all selected frames.

    Args:
        frames: List of frame dicts.
        species_a: Central atom species.
        species_b: Neighbour atom species.
        r_max: Maximum distance (Å).
        n_bins: Number of histogram bins.

    Returns:
        (r, g) — bin centres and g(r).
    """
    # Find indices
    a_idx = [i for i, a in enumerate(frames[0]["atoms"]) if a["label"] == species_a]
    b_idx = [i for i, a in enumerate(frames[0]["atoms"]) if a["label"] == species_b]
    n_a = len(a_idx)
    n_b = len(b_idx)

    if n_a == 0 or n_b == 0:
        raise ValueError(f"No atoms found for {species_a}→{species_b}")

    # Cell volume (same for all frames with NPT, approximate otherwise)
    lv = frames[0].get("lattice_vectors")
    lc = frames[0].get("lattice_constant", 1.0)
    cell = lc * np.array(lv) if lv is not None else None
    if cell is not None:
        vol = abs(np.linalg.det(cell))
    else:
        vol = 1.0

    dr_bin = r_max / n_bins
    edges = np.linspace(0, r_max, n_bins + 1)
    r = 0.5 * (edges[:-1] + edges[1:])
    hist = np.zeros(n_bins)

    from rich.progress import Progress
    with Progress() as progress:
        task = progress.add_task("[cyan]Computing RDF...", total=len(frames))
        for frame in frames:
            pos_a = np.array([frame["atoms"][i]["xyz"] for i in a_idx])
            pos_b = np.array([frame["atoms"][i]["xyz"] for i in b_idx])
            lv_f = frame.get("lattice_vectors", lv)
            lc_f = frame.get("lattice_constant", lc)
            cell_frame = lc_f * np.array(lv_f) if lv_f is not None else None

            # Broadcasting: (n_a, 1, 3) - (1, n_b, 3) → (n_a, n_b, 3)
            d = pos_a[:, None, :] - pos_b[None, :, :]
            if cell_frame is not None:
                d = _minimum_image(d.reshape(-1, 3), cell_frame).reshape(n_a, n_b, 3)
            dist = np.linalg.norm(d, axis=2)  # (n_a, n_b)
            if species_a == species_b:
                dist[dist < 1e-6] = np.inf
            h, _ = np.histogram(dist, bins=edges)
            hist += h
            progress.update(task, advance=1)

    n_frames = len(frames)
    rho_b = n_b / vol
    for i in range(n_bins):
        shell_vol = 4.0 * np.pi * r[i] ** 2 * dr_bin
        hist[i] /= (n_frames * n_a * rho_b * shell_vol)

    return r, hist


def _plot_rdf(r: np.ndarray, g: np.ndarray, species_a: str, species_b: str) -> None:
    """Plot RDF curve."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(r, g, "b-", linewidth=1.2)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("g(r)")
    ax.set_title(f"Radial Distribution Function — {species_a}–{species_b}")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

    fig.tight_layout(pad=1.2)
    save_name = f"RDF_{species_a}-{species_b}.png"
    fig.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_cn(r: np.ndarray, cn: np.ndarray, species_a: str, species_b: str) -> None:
    """Plot coordination number curve."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(r, cn, "r-", linewidth=1.2)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("CN(r)")
    ax.set_title(f"Coordination Number — {species_a}–{species_b}")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")

    fig.tight_layout(pad=1.2)
    save_name = f"CN_{species_a}-{species_b}.png"
    fig.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_rdf_cn_combined(
    r: np.ndarray, g: np.ndarray, cn: np.ndarray,
    species_a: str, species_b: str,
) -> None:
    """Plot RDF and CN on dual Y-axes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax1 = plt.subplots(figsize=(8, 6))

    # g(r) on left Y-axis
    ax1.plot(r, g, "b-", linewidth=1.2, label="g(r)")
    ax1.set_xlabel("r (Å)")
    ax1.set_ylabel("g(r)", color="b")
    ax1.tick_params(axis="y", labelcolor="b")
    ax1.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.tick_params(axis="both", direction="out")

    # CN(r) on right Y-axis
    ax2 = ax1.twinx()
    ax2.plot(r, cn, "r-", linewidth=1.2, label="CN(r)")
    ax2.set_ylabel("CN(r)", color="r")
    ax2.tick_params(axis="y", labelcolor="r", direction="out")

    ax1.set_title(f"RDF + Coordination Number — {species_a}–{species_b}")
    for spine in ax1.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.5)
        spine.set_visible(True)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    fig.tight_layout(pad=1.2)
    save_name = f"RDF_CN_{species_a}-{species_b}.png"
    fig.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Task 784: RDF
# =============================================================================


@task(3104, category="MD Analysis", name="RDF",
      description="Compute Radial Distribution Function from MD_dump trajectory")
def task_md_rdf(args: list[str] | None = None, interactive: bool = True) -> None:
    """Calculate partial radial distribution function g(r)."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Radial Distribution Function (RDF) ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    console.print(f"  [dim]Reading: {md_path}[/dim]")
    try:
        frames = parse_md_dump(md_path)
    except Exception as e:
        console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
        return

    n_frames = len(frames)
    # Collect species in trajectory order
    species_order = []
    for a in frames[0]["atoms"]:
        if a["label"] not in species_order:
            species_order.append(a["label"])
    species_counts = {s: sum(1 for a in frames[0]["atoms"] if a["label"] == s)
                      for s in species_order}

    console.print(f"  Frames: {n_frames}, Atoms: {len(frames[0]['atoms'])}")
    console.print(f"  Species: {', '.join(f'{s}({species_counts[s]})' for s in species_order)}")
    console.print()

    # --- Choose atom pairs ---
    if interactive:
        central_choices = [f"{s} ({species_counts[s]} atoms)" for s in species_order]
        cent = _prompt_choice(console, "Central (source) species", central_choices,
                              central_choices[0])
        species_a = cent.split()[0]

        # For neighbour, show all species again
        neig = _prompt_choice(console, "Neighbour (target) species", central_choices,
                              cent)  # default = same as central
        species_b = neig.split()[0]

        # Parameters
        r_max = float(_prompt(console, "Max distance r_max (Å)", "10.0"))
        n_bins = int(_prompt(console, "Number of bins", "200"))
    else:
        species_a = species_order[0]
        species_b = species_order[0]
        r_max = 10.0
        n_bins = 200

    # --- Frame range ---
    last_frame = n_frames - 1
    if interactive:
        start = int(_prompt(console, "Start frame", "500"))
        end = int(_prompt(console, "End frame", str(last_frame)))
        step = int(_prompt(console, "Frame interval", "1"))
    else:
        start, end, step = 500, last_frame, 1

    start = max(0, start)
    end = min(last_frame, end)
    step = max(1, step)
    if end <= start:
        console.print(f"[red]End frame ({end}) must be greater than start ({start}).[/red]")
        return

    selected_frames = frames[start:end:step]
    console.print(f"\n  Frames: {start}→{end} step {step} = {len(selected_frames)}")
    console.print(f"  Computing RDF {species_a}–{species_b}...")

    r, g = compute_rdf(selected_frames, species_a, species_b, r_max=r_max, n_bins=n_bins)

    # Save RDF
    out_file = f"RDF_{species_a}-{species_b}.dat"
    np.savetxt(out_file, np.column_stack([r, g]),
               fmt="%.8f  %.8f", header="# r(Ang)  g(r)")
    console.print(f"  [green]✓ RDF saved to {out_file}[/green]")

    # First peak
    peak_idx = int(np.argmax(g[1:]) + 1)
    console.print(f"  First peak: r = {r[peak_idx]:.3f} Å, g = {g[peak_idx]:.3f}")

    # --- Coordination number ---
    n_b_species = species_counts.get(species_b, 0)
    lv0 = frames[0].get("lattice_vectors")
    lc0 = frames[0].get("lattice_constant", 1.0)
    cell0 = lc0 * np.array(lv0) if lv0 is not None else None
    vol = abs(np.linalg.det(cell0)) if cell0 is not None else 1.0
    rho_b = n_b_species / vol if vol > 0 else 0.0

    compute_cn = False
    if interactive:
        want_cn = _prompt_choice(console, "Compute coordination number?",
                                  ["Yes", "No"], "Yes")
        compute_cn = "Yes" in want_cn

    if compute_cn:
        dr_bin = r[1] - r[0]
        cn = np.zeros(len(r))
        for i in range(1, len(r)):
            shell_vol = 4.0 * np.pi * r[i] ** 2 * dr_bin
            cn[i] = cn[i - 1] + g[i] * rho_b * shell_vol

        cn_file = f"CN_{species_a}-{species_b}.dat"
        np.savetxt(cn_file, np.column_stack([r, cn]),
                   fmt="%.8f  %.8f", header="# r(Ang)  CN(r)")
        console.print(f"  [green]✓ Coordination number saved to {cn_file}[/green]")

        if interactive:
            plot_mode = _prompt_choice(console, "Plot mode",
                                        ["Dual Y-axis (g(r) + CN in one figure)",
                                         "Separate figures"],
                                        "Dual Y-axis (g(r) + CN in one figure)")
            combined = "Dual" in plot_mode
        else:
            combined = True

        if combined:
            _plot_rdf_cn_combined(r, g, cn, species_a, species_b)
            console.print(f"  [green]✓ Combined plot → RDF_CN_{species_a}-{species_b}.png[/green]")
        else:
            _plot_rdf(r, g, species_a, species_b)
            _plot_cn(r, cn, species_a, species_b)
            console.print(f"  [green]✓ RDF plot → RDF_{species_a}-{species_b}.png[/green]")
            console.print(f"  [green]✓ CN plot → CN_{species_a}-{species_b}.png[/green]")
    else:
        _plot_rdf(r, g, species_a, species_b)
        console.print(f"  [green]✓ RDF plot → RDF_{species_a}-{species_b}.png[/green]")

    console.print()


# =============================================================================
# Probability density → CHGCAR + POSCAR export
# =============================================================================


def write_poscar(
    frame: dict,
    filepath: str | Path,
) -> None:
    """Write a single MD frame as a VASP POSCAR file.

    Uses the lattice vectors and atomic positions from the frame.
    """
    lv = frame.get("lattice_vectors")
    lc = frame.get("lattice_constant", 1.0)
    if lv is None:
        raise ValueError("Frame has no lattice vectors")

    # MD_dump already stores lattice_constant in Angstrom, so this is Å.
    cell = lc * np.array(lv)

    # Gather species in order
    species_order = []
    species_atoms = {}
    for a in frame["atoms"]:
        s = a["label"]
        if s not in species_order:
            species_order.append(s)
            species_atoms[s] = []
        species_atoms[s].append(a)

    filepath = Path(filepath)
    with open(filepath, "w") as f:
        f.write(f"MD step {frame['step']}\n")
        f.write("1.0\n")
        for row in cell:
            f.write(f"  {row[0]:.16f}  {row[1]:.16f}  {row[2]:.16f}\n")
        f.write("  " + "  ".join(species_order) + "\n")
        f.write("  " + "  ".join(str(len(species_atoms[s])) for s in species_order) + "\n")
        f.write("Direct\n")
        # Convert Cartesian (Å) → fractional (cell is in Å, so units match)
        cell_inv = np.linalg.inv(cell.T)
        for s in species_order:
            for a in species_atoms[s]:
                frac = a["xyz"] @ cell_inv
                f.write(f"  {frac[0]:.16f}  {frac[1]:.16f}  {frac[2]:.16f}\n")


def compute_probability_density(
    frames: list[dict],
    species: str,
    grid_size: tuple[int, int, int] = (80, 80, 80),
) -> tuple[np.ndarray, np.ndarray]:
    """Compute 3-D probability density for a species over all frames.

    Args:
        frames: List of frame dicts.
        species: Atom species to track.
        grid_size: (nx, ny, nz).

    Returns:
        (density_3d, lattice_vectors) — density as probability [0,1],
        lattice in Angstrom.
    """
    lv = frames[0].get("lattice_vectors")
    lc = frames[0].get("lattice_constant", 1.0)
    if lv is None:
        raise ValueError("No lattice vectors in trajectory")

    # MD_dump stores lattice_constant in Angstrom → cell is in Å
    cell = lc * np.array(lv)

    nx, ny, nz = grid_size
    density = np.zeros((nx, ny, nz), dtype=np.float64)

    # Pre-compute cell inverse once
    cell_inv = np.linalg.inv(cell.T)
    total_count = 0
    from rich.progress import Progress
    with Progress() as progress:
        task = progress.add_task("[cyan]Computing density...", total=len(frames))
        for frame in frames:
            for a in frame["atoms"]:
                if a["label"] != species:
                    continue
                frac = a["xyz"] @ cell_inv
                ix = int(frac[0] * nx) % nx
                iy = int(frac[1] * ny) % ny
                iz = int(frac[2] * nz) % nz
                density[ix, iy, iz] += 1
                total_count += 1
            progress.update(task, advance=1)

    if total_count > 0:
        density /= total_count

    return density, cell


def write_chgcar(
    density: np.ndarray,
    lv: np.ndarray,
    species: str,
    filepath: str | Path,
) -> None:
    """Write probability density as VASP CHGCAR for VESTA."""
    nx, ny, nz = density.shape
    vol = abs(np.linalg.det(lv))

    filepath = Path(filepath)
    with open(filepath, "w") as f:
        f.write(f"Probability density for {species}\n")
        f.write("   1.00000000000000\n")
        for row in lv:
            f.write(f"  {row[0]:.16f}  {row[1]:.16f}  {row[2]:.16f}\n")
        f.write(f"  {species}\n")
        f.write("  1\n")
        f.write("Direct\n")
        f.write("  0.0000000000000000  0.0000000000000000  0.0000000000000000\n")
        f.write("\n")
        f.write(f"  {nx}  {ny}  {nz}\n")

        # VASP CHGCAR convention: x fastest (innermost), z slowest (outermost)
        #   write(((density[x,y,z], x=1,nx), y=1,ny), z=1,nz)
        count = 0
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    val = density[ix, iy, iz] / vol * 1000.0
                    f.write(f"  {val:21.14E}")
                    count += 1
                    if count % 5 == 0:
                        f.write("\n")
        if count % 5 != 0:
            f.write("\n")


# =============================================================================
# Task 785: Probability density
# =============================================================================


@task(3105, category="MD Analysis", name="Probability Density",
      description="Compute atomic probability density from MD trajectory, export as CHGCAR + POSCAR")
def task_probability_density(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compute 3-D probability density for a species and export CHGCAR + POSCAR."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Probability Density → CHGCAR + POSCAR ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    console.print(f"  [dim]Reading: {md_path}[/dim]")
    try:
        frames = parse_md_dump(md_path)
    except Exception as e:
        console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
        return

    n_frames = len(frames)
    # Species order from trajectory
    species_order = []
    for a in frames[0]["atoms"]:
        if a["label"] not in species_order:
            species_order.append(a["label"])
    species_counts = {s: sum(1 for a in frames[0]["atoms"] if a["label"] == s)
                      for s in species_order}

    console.print(f"  Frames: {n_frames}, Atoms: {len(frames[0]['atoms'])}")
    console.print(f"  Species: {', '.join(f'{s}({species_counts[s]})' for s in species_order)}")
    console.print()

    # --- POSCAR export (first frame) ---
    try:
        write_poscar(frames[0], "STRU.vasp")
        console.print("  [green]✓ First frame → STRU.vasp[/green]")
    except Exception as e:
        console.print(f"  [yellow]! POSCAR failed: {e}[/yellow]")

    # --- Choose species ---
    if interactive:
        species_choices = [f"{s} ({species_counts[s]} atoms)" for s in species_order]
        chosen = _prompt_choice(console, "Select species for density", species_choices,
                                species_choices[0])
        species = chosen.split()[0]
        grid_in = _prompt(console, "Grid size (nx ny nz)", "80 80 80")
        try:
            g = [int(x) for x in grid_in.split()]
            grid_size = (g[0], g[1], g[2])
        except (ValueError, IndexError):
            grid_size = (80, 80, 80)
    else:
        species = species_order[0]
        grid_size = (80, 80, 80)

    # --- Frame range ---
    last_frame = n_frames - 1
    if interactive:
        start = int(_prompt(console, "Start frame", "500"))
        end = int(_prompt(console, "End frame", str(last_frame)))
        step = int(_prompt(console, "Frame interval", "1"))
    else:
        start, end, step = 500, last_frame, 1

    start = max(0, start)
    end = min(last_frame, end)
    step = max(1, step)
    if end <= start:
        console.print(f"[red]End frame ({end}) must be greater than start ({start}).[/red]")
        return

    selected = frames[start:end:step]
    console.print(f"\n  Frames: {start}→{end} step {step} = {len(selected)}")
    console.print(f"  Grid: {grid_size[0]}×{grid_size[1]}×{grid_size[2]}")
    console.print(f"  Computing density for {species}...")

    density, lv = compute_probability_density(selected, species, grid_size=grid_size)

    out_file = f"probability_{species}.vasp"
    write_chgcar(density, lv, species, out_file)
    console.print(f"  [green]✓ Density saved to {out_file}[/green]")
    console.print("  [dim]Open with VESTA: Properties → Isosurfaces[/dim]")

    console.print()


# =============================================================================
# van Hove correlation functions
# =============================================================================


def _pairwise_dist_mic(
    pos_a: np.ndarray, pos_b: np.ndarray,
    cell: np.ndarray | None, cell_inv: np.ndarray | None,
) -> np.ndarray:
    """Pairwise distances with minimum-image convention.

    Args:
        pos_a: (n_a, 3) positions.
        pos_b: (n_b, 3) positions.
        cell: (3, 3) cell matrix (Å), or None for no PBC.
        cell_inv: (3, 3) inverse cell.T, precomputed.

    Returns:
        (n_a, n_b) distance matrix.
    """
    # Try PyTorch MPS on Apple Silicon for large arrays (>10K pairs)
    if pos_a.shape[0] * pos_b.shape[0] > 10000:
        try:
            import torch
            if hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                device = torch.device("mps")
                # torch.cdist is the fastest way to get pairwise distances on GPU
                da = torch.from_numpy(pos_a.astype(np.float32)).to(device).unsqueeze(0)
                db = torch.from_numpy(pos_b.astype(np.float32)).to(device).unsqueeze(0)
                if cell is not None:
                    # Apply minimum-image in batch on GPU
                    L = torch.from_numpy(
                        np.linalg.norm(cell, axis=1).astype(np.float32)
                    ).to(device)
                    delta = da.squeeze(0)[:, None, :] - db.squeeze(0)[None, :, :]
                    delta -= torch.round(delta / L) * L
                    dist = torch.norm(delta, dim=2)
                else:
                    dist = torch.cdist(da, db).squeeze(0)
                return dist.cpu().numpy().astype(np.float64)
        except (ImportError, RuntimeError):
            pass

    # NumPy path
    delta = pos_a[:, None, :] - pos_b[None, :, :]
    if cell is not None:
        frac = delta.reshape(-1, 3) @ cell_inv
        frac -= np.floor(frac + 0.5)
        delta = (frac @ cell.T).reshape(pos_a.shape[0], pos_b.shape[0], 3)
    return np.linalg.norm(delta, axis=2)


def _unwrap_trajectory(frames: list[dict], atom_indices: list[int]) -> np.ndarray:
    """Remove periodic boundary jumps from a trajectory.

    Uses the same algorithm as the reference GPU script:
        delta = coords[k] - coords[k-1]
        delta -= round(delta / L) * L
    where L = cell box lengths (orthorhombic) or full matrix (triclinic).

    Returns (n_atoms, 3, n_frames) unwrapped coordinates in Angstrom.
    """
    n_atoms = len(atom_indices)
    n_frames = len(frames)

    coords = np.zeros((n_atoms, 3, n_frames), dtype=np.float64)
    for k, frame in enumerate(frames):
        for i, idx in enumerate(atom_indices):
            coords[i, :, k] = frame["atoms"][idx]["xyz"]

    unwrapped = np.zeros_like(coords)
    unwrapped[:, :, 0] = coords[:, :, 0]

    # Get box lengths from first frame — use diagonal or full matrix
    lv0 = frames[0].get("lattice_vectors")
    lc0 = frames[0].get("lattice_constant", 1.0)
    use_orthorhombic = lv0 is not None and np.allclose(lv0 - np.diag(np.diag(lv0)), 0)

    for k in range(1, n_frames):
        lv = frames[k - 1].get("lattice_vectors")
        lc = frames[k - 1].get("lattice_constant", 1.0)
        delta = coords[:, :, k] - coords[:, :, k - 1]

        if lv is not None:
            if use_orthorhombic:
                L = lc * np.diag(lv)  # fast: just box lengths
                delta -= np.round(delta / L) * L
            else:
                cell = lc * np.array(lv)
                cell_inv = np.linalg.inv(cell.T)
                frac = delta @ cell_inv
                frac -= np.round(frac)
                delta = frac @ cell.T

        unwrapped[:, :, k] = unwrapped[:, :, k - 1] + delta

    return unwrapped


def compute_van_hove(
    frames: list[dict],
    species_a: str,
    species_b: str | None = None,
    *,
    r_max: float = 8.0,
    dr: float = 0.05,
    frame_stride: int = 10,
) -> dict:
    """Compute van Hove correlation functions (Gs, Gd) and NGP.

    MSD is intentionally excluded — use task 3103 for proper ensemble-averaged MSD.

    Directional self-part Gs is included: P(|Δx|,t), P(|Δy|,t), P(|Δz|,t),
    i.e. 1-D histograms of the unwrapped displacement components. These are
    normalized by the atom count only (no 4πr² shell factor) so they are
    probability distributions per bin, matching the reference directional
    van Hove script.

    Args:
        frames: List of frame dicts.
        species_a: Central atom species (for Gs and Gd).
        species_b: Neighbour species for Gd (None → same as species_a).
        r_max: Max distance for histogram (Å).
        dr: Bin width (Å).
        frame_stride: Only compute every N-th lag to save time.

    Returns:
        Dict with keys: r, time_lags, gs_matrix, gsx_matrix, gsy_matrix,
        gsz_matrix, gd_matrix, ngp.
    """
    # Find atom indices
    a_idx = sorted(i for i, a in enumerate(frames[0]["atoms"])
                   if a["label"] == species_a)
    n_a = len(a_idx)

    if species_b and species_b != species_a:
        b_idx = sorted(i for i, a in enumerate(frames[0]["atoms"])
                       if a["label"] == species_b)
    else:
        b_idx = a_idx
        species_b = species_a
    n_b = len(b_idx)

    # Volume
    lv = frames[0].get("lattice_vectors")
    lc0 = frames[0].get("lattice_constant", 1.0)
    cell0 = lc0 * np.array(lv) if lv is not None else None
    vol = abs(np.linalg.det(cell0)) if cell0 is not None else 1.0
    rho = n_b / vol  # number density of species_b

    # Histogram bins
    r_edges = np.arange(0, r_max + dr, dr)
    r_centers = (r_edges[:-1] + r_edges[1:]) / 2
    n_r = len(r_centers)
    n_frames = len(frames)

    # Unwrap species_a
    print("  Unwrapping trajectory...")
    coords_unwrapped = _unwrap_trajectory(frames, a_idx)  # (n_a, 3, n_frames)

    # Reference positions at t=0
    r0_unwrapped = coords_unwrapped[:, :, 0]  # (n_a, 3)

    # Non-Gaussian parameter: α₂(t) = 3<r⁴>/(5<r²>²) - 1
    # Pre-compute all species_a coords: (n_a, 3, n_frames)
    coords_a = np.zeros((n_a, 3, n_frames), dtype=np.float64)
    for k in range(n_frames):
        for ia, idx in enumerate(a_idx):
            coords_a[ia, :, k] = frames[k]["atoms"][idx]["xyz"]

    ngp = np.zeros(n_frames, dtype=np.float64)
    lv = frames[0].get("lattice_vectors")
    cell_inv = np.linalg.inv(cell0.T) if cell0 is not None else None
    for k in range(1, n_frames):
        delta = coords_a[:, :, k] - coords_a[:, :, 0]
        if cell0 is not None:
            frac = delta @ cell_inv
            frac -= np.floor(frac + 0.5)
            delta = frac @ cell0.T
        r2_per_atom = np.sum(delta ** 2, axis=1)
        r2 = np.mean(r2_per_atom)
        if r2 > 1e-12:
            r4 = np.mean(r2_per_atom ** 2)
            ngp[k] = (3.0 * r4) / (5.0 * r2 ** 2) - 1.0

    # --- Gs(r,t): self-part ---
    # Sample time lags
    lag_indices = list(range(0, n_frames, frame_stride))
    n_lags = len(lag_indices)
    gs_matrix = np.zeros((n_r, n_lags), dtype=np.float64)
    # Directional self-part: P(|Δx|,t), P(|Δy|,t), P(|Δz|,t)
    gsx_matrix = np.zeros((n_r, n_lags), dtype=np.float64)
    gsy_matrix = np.zeros((n_r, n_lags), dtype=np.float64)
    gsz_matrix = np.zeros((n_r, n_lags), dtype=np.float64)
    shell_vol = 4.0 * np.pi * r_centers ** 2 * dr

    print(f"  Computing Gs(r,t) + directional Gs(x/y/z) at {n_lags} time lags...")
    for li, lag in enumerate(lag_indices):
        dr_vec = coords_unwrapped[:, :, lag] - r0_unwrapped  # (n_a, 3)
        dr_mag = np.linalg.norm(dr_vec, axis=1)             # (n_a,)
        hist, _ = np.histogram(dr_mag, bins=r_edges)
        gs_matrix[:, li] = hist / (n_a * shell_vol)
        # Directional: 1-D histograms of |Δx|, |Δy|, |Δz|, normalized by
        # atom count only (no 4πr² shell factor) — probability per bin.
        for axis, out in enumerate((gsx_matrix, gsy_matrix, gsz_matrix)):
            hist_axis, _ = np.histogram(np.abs(dr_vec[:, axis]), bins=r_edges)
            out[:, li] = hist_axis / n_a

    # --- Gd(r,t): distinct-part ---
    # Pre-compute species_b trajectory for all frames: (n_b, 3, n_frames)
    b_coords_all = np.zeros((n_b, 3, n_frames), dtype=np.float64)
    for k in range(n_frames):
        for i, idx in enumerate(b_idx):
            b_coords_all[i, :, k] = frames[k]["atoms"][idx]["xyz"]
    b_coords_0 = b_coords_all[:, :, 0]  # (n_b, 3)

    gd_matrix = np.zeros((n_r, n_lags), dtype=np.float64)
    print(f"  Computing Gd(r,t) at {n_lags} time lags...")
    from rich.progress import Progress
    with Progress() as progress:
        task = progress.add_task("[cyan]Computing Gd(r,t)...", total=n_lags)
        for li, lag in enumerate(lag_indices):
            b_coords_t = b_coords_all[:, :, lag]  # (n_b, 3) — no dict access

            dist = _pairwise_dist_mic(b_coords_0, b_coords_t, cell0, cell_inv)
            if species_a == species_b:
                np.fill_diagonal(dist, np.inf)
            hist, _ = np.histogram(dist, bins=r_edges)
            gd_matrix[:, li] = hist / (n_b * n_b * shell_vol * rho)
            progress.update(task, advance=1)

    result = {
        "r": r_centers,
        "time_lags": np.array(lag_indices, dtype=int),
        "time_frames": np.arange(n_frames),
        "gs_matrix": gs_matrix,
        "gsx_matrix": gsx_matrix,
        "gsy_matrix": gsy_matrix,
        "gsz_matrix": gsz_matrix,
        "gd_matrix": gd_matrix,
        "ngp": ngp,
        "species_a": species_a,
        "species_b": species_b,
    }
    return result


def _plot_ngp_only(result: dict, prefix: str, dt_fs: float = 1.0) -> None:
    """Plot NGP curve."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    time_ps = result["time_frames"] * dt_fs / 1000.0
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(time_ps, result["ngp"], "r-", linewidth=1.2)
    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Time lag (ps)")
    ax.set_ylabel(r"$\alpha_2(t)$")
    ax.set_title(f"Non-Gaussian Parameter — {result['species_a']}")
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="x", scilimits=(-2, 3))
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    fig.savefig(f"{prefix}_NGP.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_gs_heatmap(r, gs_matrix, lags_ps, sa, prefix, vmax):
    """Plot Gs heatmap with user-specified color max."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    gs_shell = 4.0 * np.pi * r[:, np.newaxis] ** 2 * gs_matrix
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(gs_shell, extent=[lags_ps[0], lags_ps[-1], r[0], r[-1]],
                   aspect="auto", origin="lower", cmap="turbo",
                   vmin=0, vmax=vmax)
    ax.set_xlabel("Time lag (ps)")
    ax.set_ylabel(r"$r$ (Å)")
    ax.set_title(rf"$4\pi r^2 \cdot G_s(r,t)$ — {sa}")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$4\pi r^2 \cdot G_s(r,t)$")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    fig.savefig(f"{prefix}_Gs_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_gs_directional_heatmap(r, gs_axis_matrix, lags_ps, sa, prefix, axis, vmax):
    """Plot a directional self-part van Hove heatmap P(|d_axis|,t).

    axis: 0 → x, 1 → y, 2 → z. Plots the per-bin probability distribution
    directly (no 4πr² shell factor), matching the reference script.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    axis_label = {0: "Δx", 1: "Δy", 2: "Δz"}[axis]
    out_tag = {0: "X", 1: "Y", 2: "Z"}[axis]

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(gs_axis_matrix, extent=[lags_ps[0], lags_ps[-1], r[0], r[-1]],
                   aspect="auto", origin="lower", cmap="turbo",
                   vmin=0, vmax=vmax)
    ax.set_xlabel("Time lag (ps)")
    ax.set_ylabel(f"|{axis_label}| (Å)")
    ax.set_title(f"P(|{axis_label}|,t) — {sa}")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(f"P(|{axis_label}|,t)")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    fig.savefig(f"{prefix}_Gs_{out_tag}_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_gd_heatmap(r, gd_matrix, lags_ps, sa, sb, prefix, vmax):
    """Plot Gd heatmap with user-specified color max."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(gd_matrix, extent=[lags_ps[0], lags_ps[-1], r[0], r[-1]],
                   aspect="auto", origin="lower", cmap="turbo",
                   vmin=0, vmax=vmax)
    ax.set_xlabel("Time lag (ps)")
    ax.set_ylabel(r"$r$ (Å)")
    ax.set_title(rf"$G_d(r,t)$ — {sa}–{sb}")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$G_d(r,t)$")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    fig.savefig(f"{prefix}_Gd_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_gs_slices(r, windows, tick_labels, species, prefix, smooth=0):
    """Plot r^2*Gs(r) curves for multiple time windows with a colorbar.

    Args:
        tick_labels: two formatted time strings — the first window's start
            time and the last window's end time (ps), shown at the two ends
            of the colorbar.
        smooth: moving-average window size (0 = no smoothing).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import turbo
    from matplotlib.colors import Normalize

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    data = np.array(windows)
    if smooth > 1:
        kernel = np.ones(smooth) / smooth
        for i in range(len(data)):
            data[i] = np.convolve(data[i], kernel, mode="same")

    n_curves = len(data)
    cmap = turbo
    norm = Normalize(vmin=0, vmax=n_curves - 1)

    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(n_curves):
        ax.plot(r, data[i], color=cmap(norm(i)), linewidth=1.0, alpha=0.85)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=n_curves - 1))
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Time (ps)")
    cbar.set_ticks([0, n_curves - 1])
    cbar.set_ticklabels(tick_labels)

    ax.set_xlabel(r"$r$ (Å)")
    ax.set_ylabel(r"$r^2 \cdot G_s(r)$")
    ax.set_title(f"Self-part van Hove — {species}")
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")
    fig.tight_layout(pad=1.2)
    out = f"{prefix}_Gs_slices.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# van Hove cache helpers (task 3106)
# =============================================================================

# Keys a cached van Hove .npz must contain to be reusable.  Files saved by
# older versions (before directional-Gs / metadata were added) are rejected so
# they are silently recomputed instead of producing empty directional plots.
_VH_CACHE_REQUIRED = {"r", "time_lags", "gs", "gsx", "gsy", "gsz", "gd", "ngp"}


def _load_van_hove_cache(npz_path: Path) -> dict | None:
    """Load a previously saved van Hove .npz cache into a ``compute_van_hove``-style dict.

    Returns None if the file is unreadable or predates the directional-Gs
    matrices, in which case the caller should recompute from the trajectory.
    """
    try:
        data = np.load(npz_path)
    except Exception:
        return None

    if not _VH_CACHE_REQUIRED.issubset(data.files):
        return None

    species_a = (str(data["species_a"].item())
                 if "species_a" in data else npz_path.stem.split("_", 1)[1])
    species_b = (str(data["species_b"].item())
                 if "species_b" in data else species_a)

    meta_keys = {"r_max", "dr", "frame_stride", "start_frame", "end_frame",
                 "frame_step", "frame_dt"}
    return {
        "r": data["r"],
        "time_lags": np.asarray(data["time_lags"], dtype=int),
        "time_frames": (data["time_frames"] if "time_frames" in data
                        else np.arange(len(data["ngp"]))),
        "gs_matrix": data["gs"],
        "gsx_matrix": data["gsx"],
        "gsy_matrix": data["gsy"],
        "gsz_matrix": data["gsz"],
        "gd_matrix": data["gd"],
        "ngp": data["ngp"],
        "species_a": species_a,
        "species_b": species_b,
        "meta": ({k: data[k] for k in meta_keys} if meta_keys.issubset(data.files) else None),
    }


# =============================================================================
# Task 786: van Hove
# =============================================================================


@task(3106, category="MD Analysis", name="van Hove",
      description="Compute van Hove functions Gs(r,t), directional Gs(x/y/z), Gd(r,t) and NGP from MD trajectory")
def task_van_hove(args: list[str] | None = None, interactive: bool = True) -> None:
    """Calculate self-part, directional self-part and distinct-part van Hove functions."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== van Hove Correlation Functions ===[/bold cyan]")
    console.print()

    md_path = _find_md_dump(args)
    if interactive:
        inp = _prompt(console, "MD_dump file path", md_path)
        if inp:
            md_path = inp

    if not Path(md_path).exists():
        console.print(f"[red]File not found: {md_path}[/red]")
        return

    # --- Cache reuse: skip the slow trajectory read when a valid .npz exists ---
    result = None
    prefix = None
    if interactive:
        reusable = []
        for c in sorted(Path.cwd().glob("vanHove_*.npz")):
            r = _load_van_hove_cache(c)
            if r is not None:
                reusable.append((c, r))
            else:
                console.print(f"  [dim]Ignoring outdated cache {c.name} (will recompute).[/dim]")

        if reusable:
            console.print("  [dim]Found cached van Hove results:[/dim]")
            for c, r in reusable:
                m = r["meta"]
                if m:
                    console.print(f"    [dim]{c.name} — {r['species_a']}, r_max={float(m['r_max']):.2f}, "
                                  f"dr={float(m['dr']):.2f}, stride={int(m['frame_stride'])}, "
                                  f"frames={int(m['start_frame'])}→{int(m['end_frame'])} "
                                  f"step {int(m['frame_step'])}[/dim]")
                else:
                    console.print(f"    [dim]{c.name} — {r['species_a']} (old format)[/dim]")
            reuse_map = {f"Reuse {c.name}": (c, r) for c, r in reusable}
            choices = list(reuse_map) + ["No — recompute"]
            reuse_choice = _prompt_choice(console, "Reuse cached data (skip reading trajectory)?",
                                          choices, choices[0])
            if reuse_choice in reuse_map:
                cache_path, result = reuse_map[reuse_choice]
                prefix = cache_path.stem
                console.print(f"  [green]✓ Loaded cached result for {result['species_a']} "
                              f"from {cache_path.name}[/green]")
                _meta = result.get("meta") or {}
                if _meta and "frame_dt" in _meta:
                    frame_dt = float(_meta["frame_dt"])
                    console.print(f"  [dim]frame_dt = {frame_dt:.1f} fs (from cache)[/dim]")
                else:
                    frame_dt = _resolve_frame_dt(md_path, console, interactive)

    if result is None:
        console.print(f"  [dim]Reading: {md_path}[/dim]")
        try:
            frames = parse_md_dump(md_path)
        except Exception as e:
            console.print(f"[red]Failed to parse MD_dump: {e}[/red]")
            return

        n_frames = len(frames)
        species_order = []
        for a in frames[0]["atoms"]:
            if a["label"] not in species_order:
                species_order.append(a["label"])
        species_counts = {s: sum(1 for a in frames[0]["atoms"] if a["label"] == s)
                          for s in species_order}

        console.print(f"  Frames: {n_frames}, Atoms: {len(frames[0]['atoms'])}")
        console.print(f"  Species: {', '.join(f'{s}({species_counts[s]})' for s in species_order)}")
        console.print()

        # --- Choose species ---
        if interactive:
            species_choices = [f"{s} ({species_counts[s]} atoms)" for s in species_order]
            cent = _prompt_choice(console, "Select species", species_choices,
                                  species_choices[0])
            species_a = cent.split()[0]

            r_max = float(_prompt(console, "Max distance r_max (Å)", "8.0"))
            dr = float(_prompt(console, "Bin width dr (Å)", "0.05"))
            f_stride = int(_prompt(console, "Time-lag stride (1=all frames)", "10"))
        else:
            species_a = species_order[0]
            r_max, dr, f_stride = 8.0, 0.05, 10

        # --- Frame range ---
        last_frame = n_frames - 1
        if interactive:
            start = int(_prompt(console, "Start frame", "0"))
            end = int(_prompt(console, "End frame", str(last_frame)))
            step = int(_prompt(console, "Frame interval", "1"))
        else:
            start, end, step = 0, last_frame, 1

        start = max(0, start)
        end = min(last_frame, end)
        step = max(1, step)
        if end <= start:
            console.print("[red]End frame must be greater than start.[/red]")
            return

        selected = frames[start:end:step]
        console.print(f"\n  Frames: {start}→{end} step {step} = {len(selected)}")

        # Real time between frames — needed for every time axis; ask BEFORE the
        # (slow) computation, not after.
        frame_dt = _resolve_frame_dt(md_path, console, interactive)

        console.print(f"  Computing van Hove for {species_a}...")
        console.print(f"  r_max={r_max} Å, dr={dr} Å, time stride={f_stride}")

        result = compute_van_hove(
            selected, species_a,
            r_max=r_max, dr=dr, frame_stride=f_stride,
        )

        # Save data (incl. parameters so the .npz can be reused as a cache later)
        prefix = f"vanHove_{species_a}"
        np.savez(f"{prefix}.npz",
                 r=result["r"], time_lags=result["time_lags"],
                 time_frames=result["time_frames"],
                 gs=result["gs_matrix"],
                 gsx=result["gsx_matrix"],
                 gsy=result["gsy_matrix"],
                 gsz=result["gsz_matrix"],
                 gd=result["gd_matrix"], ngp=result["ngp"],
                 species_a=result["species_a"], species_b=result["species_b"],
                 r_max=r_max, dr=dr, frame_stride=f_stride,
                 start_frame=start, end_frame=end, frame_step=step,
                 frame_dt=frame_dt)
        console.print(f"  [green]✓ Data saved to {prefix}.npz[/green]")

        np.savetxt(f"{prefix}_NGP.dat",
                   np.column_stack([result["time_frames"], result["ngp"]]),
                   fmt="%d  %.8f", header="# frame  NGP")

    # species_a must be defined for the shared plotting section (fresh or cached)
    species_a = result["species_a"]

    # --- NGP plot (always made) ---
    _plot_ngp_only(result, prefix, frame_dt)
    console.print(f"  [green]✓ NGP plot: {prefix}_NGP.png")

    # NGP peak info
    ngp_peak_idx = np.argmax(result["ngp"])
    console.print(f"  NGP peak: α₂ = {result['ngp'][ngp_peak_idx]:.4f} at frame {ngp_peak_idx}")

    # --- Gs heatmap — interactive color range ---
    gs_shell = 4.0 * np.pi * result["r"][:, np.newaxis] ** 2 * result["gs_matrix"]
    gs_pos = gs_shell[gs_shell > 0]
    lags_ps = result["time_lags"] * frame_dt / 1000.0

    if interactive:
        while True:
            console.print()
            console.print(f"  [dim]Gs data range: [{np.min(gs_pos):.4e}, {np.max(gs_pos):.4e}][/dim]")
            vmax_in = _prompt(console, "  Color max for Gs heatmap (0=auto)", "0")
            try:
                vmax_gs = float(vmax_in) if float(vmax_in) > 0 else float(np.max(gs_pos)) / 40.0
            except ValueError:
                vmax_gs = float(np.max(gs_pos)) / 40.0
            _plot_gs_heatmap(result["r"], result["gs_matrix"], lags_ps,
                             result["species_a"], prefix, vmax_gs)
            console.print(f"  [green]✓ Gs heatmap → {prefix}_Gs_heatmap.png[/green]")
            redo = _prompt_choice(console, "Redo Gs?", ["No", "Yes"], "No")
            if "No" in redo:
                break
    else:
        _plot_gs_heatmap(result["r"], result["gs_matrix"], lags_ps,
                         result["species_a"], prefix, float(np.max(gs_pos)) / 40.0)
        console.print(f"  [green]✓ Gs heatmap → {prefix}_Gs_heatmap.png[/green]")

    # --- Directional Gs (x/y/z) heatmaps — P(|Δx|,t), P(|Δy|,t), P(|Δz|,t) ---
    gs_dir_matrices = (result["gsx_matrix"], result["gsy_matrix"], result["gsz_matrix"])
    gs_dir_axes = "XYZ"

    if interactive:
        console.print()
        do_dir = _prompt_choice(console, "Plot directional Gs (x/y/z) heatmaps?",
                                ["Yes", "No"], "Yes")
        if "Yes" in do_dir:
            dir_max = max(float(np.max(m[m > 0])) for m in gs_dir_matrices)
            while True:
                console.print(f"  [dim]Directional Gs range: [0, {dir_max:.4e}][/dim]")
                vmax_in = _prompt(console, "  Color max for directional Gs (0=auto)", "0")
                try:
                    vmax_dir = float(vmax_in) if float(vmax_in) > 0 else dir_max / 40.0
                except ValueError:
                    vmax_dir = dir_max / 40.0
                for axis, m in enumerate(gs_dir_matrices):
                    _plot_gs_directional_heatmap(result["r"], m, lags_ps,
                                                 result["species_a"], prefix, axis, vmax_dir)
                    console.print(f"  [green]✓ Directional Gs → "
                                  f"{prefix}_Gs_{gs_dir_axes[axis]}_heatmap.png[/green]")
                redo = _prompt_choice(console, "Redo directional Gs?", ["No", "Yes"], "No")
                if "No" in redo:
                    break
    else:
        dir_max = max(float(np.max(m[m > 0])) for m in gs_dir_matrices)
        for axis, m in enumerate(gs_dir_matrices):
            _plot_gs_directional_heatmap(result["r"], m, lags_ps,
                                         result["species_a"], prefix, axis, dir_max / 40.0)
        console.print("  [green]✓ Directional Gs heatmaps → "
                      f"{prefix}_Gs_X/Y/Z_heatmap.png[/green]")

    # --- Gd heatmap — interactive color range ---
    gd_pos = result["gd_matrix"][result["gd_matrix"] > 0]

    if interactive:
        while True:
            console.print()
            console.print(f"  [dim]Gd data range: [{np.min(gd_pos):.4e}, {np.max(gd_pos):.4e}][/dim]")
            vmax_in = _prompt(console, "  Color max for Gd heatmap (0=auto)", "0")
            try:
                vmax_gd = float(vmax_in) if float(vmax_in) > 0 else float(np.max(gd_pos)) / 40.0
            except ValueError:
                vmax_gd = float(np.max(gd_pos)) / 40.0
            _plot_gd_heatmap(result["r"], result["gd_matrix"], lags_ps,
                             result["species_a"], result["species_b"], prefix, vmax_gd)
            console.print(f"  [green]✓ Gd heatmap → {prefix}_Gd_heatmap.png[/green]")
            redo = _prompt_choice(console, "Redo Gd?", ["No", "Yes"], "No")
            if "No" in redo:
                break
    else:
        _plot_gd_heatmap(result["r"], result["gd_matrix"], lags_ps,
                         result["species_a"], result["species_b"], prefix,
                         float(np.max(gd_pos)) / 40.0)
        console.print(f"  [green]✓ Gd heatmap → {prefix}_Gd_heatmap.png[/green]")

    # --- Time-window sliced Gs(r) curves ---
    if interactive:
        console.print()
        do_window = _prompt_choice(console, "Plot time-sliced r^2*Gs(r) curves?",
                                   ["Yes", "No"], "Yes")
        if "Yes" in do_window:
            lags = result["time_lags"]
            n_total = len(result["time_frames"])

            t_start = int(_prompt(console, f"  Start frame (0–{n_total - 1})",
                                 "0"))
            t_end = int(_prompt(console, f"  End frame (0–{n_total - 1})",
                               str(n_total - 1)))
            default_intv = max(1, (n_total - 1) // 10)
            interval = int(_prompt(console, "  Window width (frames)",
                                  str(default_intv)))
            n_slices = int(_prompt(console, "  Number of windows", "10"))

            t_start = max(lags[0], t_start)
            t_end = min(lags[-1], t_end)
            interval = max(1, interval)
            n_slices = max(2, min(n_slices, 20))

            if t_end <= t_start or t_end - t_start < interval:
                console.print("[red]Range too small for the given window width[/red]")
            else:
                r = result["r"]
                # Evenly spaced window starts across [t_start, t_end-interval]
                total_span = t_end - t_start
                n_slices = min(n_slices, total_span // interval)

                def _fmt_ps(v: float) -> str:
                    """Round to a 'nice' 1-2-5 value and format compactly.

                    Display-friendly for non-experts: 199.8 -> '200', 19.98 -> '20',
                    0.5 -> '0.5'. The colorbar is only a qualitative guide.
                    """
                    if v <= 0:
                        v = 0.0
                    else:
                        base = 10.0 ** np.floor(np.log10(v))
                        for nice in (1, 2, 5, 10):
                            if v <= nice * base:
                                v = nice * base
                                break
                    return f"{v:.3f}".rstrip("0").rstrip(".")

                windows = []
                win_starts = []  # formatted start time (ps) of each window
                win_ends = []    # formatted end time (ps) of each window
                for i in range(n_slices):
                    w0 = t_start + i * (total_span // n_slices)
                    w1 = w0 + interval
                    if w1 > t_end:
                        w1 = t_end
                    idx0 = np.searchsorted(lags, w0)
                    idx1 = np.searchsorted(lags, w1)
                    if idx1 <= idx0:
                        continue
                    gs_avg = np.mean(result["gs_matrix"][:, idx0:idx1 + 1], axis=1)
                    r2_gs = r ** 2 * gs_avg
                    windows.append(r2_gs)
                    # Time in ps = frame × frame_dt(fs) / 1000
                    win_starts.append(_fmt_ps(w0 * frame_dt / 1000.0))
                    win_ends.append(_fmt_ps(w1 * frame_dt / 1000.0))

                if not windows:
                    console.print("[red]No valid windows (interval too small?)[/red]")
                else:
                    windows = np.array(windows)  # (n_windows, n_r)
                    # Colorbar ticks: first window's start, last window's end
                    tick_labels = [win_starts[0], win_ends[-1]]
                    _plot_gs_slices(r, windows, tick_labels, species_a, prefix)
                    console.print(f"  [green]✓ Plot: {prefix}_Gs_slices.png[/green]")
                    # --- Smoothing (re-prompt after the first pass) ---
                    smoothed = False
                    while True:
                        q_smooth = "Re-apply smoothing?" if smoothed else "Apply smoothing?"
                        want = _prompt_choice(console, q_smooth, ["No", "Yes"], "No")
                        if "No" in want:
                            break
                        sw_in = _prompt(console, "  Smoothing window (odd number, 3–21)", "5")
                        try:
                            sw = int(sw_in)
                            if sw < 2:
                                break
                        except ValueError:
                            break
                        _plot_gs_slices(r, windows, tick_labels, species_a, prefix, smooth=sw)
                        console.print(f"  [green]✓ Plot: {prefix}_Gs_slices.png (smooth={sw})[/green]")
                        smoothed = True

    console.print()


# =============================================================================
# Task 3107: LAMMPS dump → ABACUS MD_dump
# =============================================================================

@task(3107, category="MD Analysis", name="LAMMPS→MD_dump",
      description="Convert LAMMPS trajectory dump to ABACUS MD_dump format")
def task_lammps_to_md_dump(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert a LAMMPS dump trajectory to ABACUS MD_dump format.

    LAMMPS "real" units: distance = Å, time = fs.
    ABACUS output: LATTICE_CONSTANT = 1.0 Å, coordinates in Cartesian Å.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== LAMMPS Dump → ABACUS MD_dump ===[/bold cyan]")
    console.print()

    lammps_path = None
    if args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                lammps_path = p
                break
    if interactive and lammps_path is None:
        inp = _prompt(console, "LAMMPS dump file path", "")
        if inp:
            lammps_path = Path(inp)

    if lammps_path is None or not Path(lammps_path).exists():
        console.print("[red]LAMMPS dump file not found.[/red]")
        return

    console.print(f"  [dim]Scanning: {lammps_path}[/dim]")
    type_set: set[int] = set()
    n_atoms = 0
    with open(lammps_path) as f:
        in_atoms = False
        for line in f:
            line = line.strip()
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                continue
            if line.startswith("ITEM: ATOMS"):
                in_atoms = True
                continue
            if in_atoms:
                if not line or line.startswith("ITEM:"):
                    break
                parts = line.split()
                type_set.add(int(parts[1]))
                n_atoms += 1

    console.print(f"  Atoms per frame: {n_atoms}")
    console.print(f"  Atom types: {sorted(type_set)}")
    console.print()

    type_map: dict[int, str] = {}
    if interactive:
        console.print("[bold]Map LAMMPS types to element symbols:[/bold]")
        for t in sorted(type_set):
            elem = _prompt(console, f"  Type {t} → element symbol", "")
            if not elem:
                console.print("[red]Element symbol required for each type.[/red]")
                return
            type_map[t] = elem.strip()
    else:
        for t in sorted(type_set):
            type_map[t] = f"X{t}"

    console.print()
    console.print(f"  Mapping: {type_map}")
    console.print()

    from rich.progress import Progress

    out_path = Path(str(lammps_path) + ".ABACUS.dump")
    n_frames = 0
    n_written = 0

    with open(lammps_path) as fin, open(out_path, "w") as fout:
        frame_lines: list[str] = []
        in_atoms = False
        box_bounds: list[tuple[float, float]] = []
        timestep = 0
        expect_timestep = False
        expect_bounds = 0

        with Progress() as progress:
            task = progress.add_task("[cyan]Converting...", total=None)

            for line in fin:
                s = line.strip()

                if s.startswith("ITEM: TIMESTEP"):
                    if frame_lines and len(box_bounds) == 3:
                        _write_abacus_md_frame(fout, frame_lines, box_bounds,
                                                type_map, timestep)
                        n_written += 1
                    frame_lines = []
                    in_atoms = False
                    box_bounds = []
                    expect_timestep = True
                    expect_bounds = 0
                    n_frames += 1
                    progress.update(task, advance=1,
                                    description=f"[cyan]Converting... frame {n_frames}")
                    continue

                if expect_timestep:
                    try:
                        timestep = int(s)
                    except ValueError:
                        pass
                    expect_timestep = False
                    continue

                if s.startswith("ITEM: BOX BOUNDS"):
                    expect_bounds = 3
                    continue

                if expect_bounds > 0:
                    parts = s.split()
                    if len(parts) >= 2:
                        box_bounds.append((float(parts[0]), float(parts[1])))
                    expect_bounds -= 1
                    continue

                if s.startswith("ITEM: ATOMS"):
                    in_atoms = True
                    continue

                if s.startswith("ITEM:"):
                    continue

                if in_atoms:
                    if s:
                        frame_lines.append(s)
                    else:
                        in_atoms = False

        if frame_lines and len(box_bounds) == 3:
            _write_abacus_md_frame(fout, frame_lines, box_bounds, type_map, timestep)
            n_written += 1

    console.print()
    console.print(f"[green]✓ Converted {n_written} frames → {out_path}[/green]")
    console.print(f"  Cell: {box_bounds[0][1]-box_bounds[0][0]:.4f} × "
                  f"{box_bounds[1][1]-box_bounds[1][0]:.4f} × "
                  f"{box_bounds[2][1]-box_bounds[2][0]:.4f} Å³")
    console.print()


def _write_abacus_md_frame(
    fout,
    atom_lines: list[str],
    box_bounds: list[tuple[float, float]],
    type_map: dict[int, str],
    timestep: int,
) -> None:
    """Write a single frame in ABACUS MD_dump format."""
    lx = box_bounds[0][1] - box_bounds[0][0]
    ly = box_bounds[1][1] - box_bounds[1][0]
    lz = box_bounds[2][1] - box_bounds[2][0]

    fout.write(f"MDSTEP: {timestep}\n")
    fout.write("LATTICE_CONSTANT: 1.0 Angstrom\n")
    fout.write("LATTICE_VECTORS\n")
    fout.write(f"   {lx:15.10f}  {0:15.10f}  {0:15.10f}\n")
    fout.write(f"   {0:15.10f}  {ly:15.10f}  {0:15.10f}\n")
    fout.write(f"   {0:15.10f}  {0:15.10f}  {lz:15.10f}\n")
    fout.write("\n")
    fout.write("INDEX   LABEL   X           Y           Z           "
               "FX   FY   FZ   VX   VY   VZ\n")

    # Sort by original LAMMPS ID — ensures consistent atom order across frames
    parsed = []
    for line in atom_lines:
        parts = line.split()
        parsed.append((int(parts[0]), int(parts[1]),
                        float(parts[2]), float(parts[3]), float(parts[4])))
    parsed.sort(key=lambda t: t[0])

    for idx, (atom_id, lmp_type, x, y, z) in enumerate(parsed, start=1):
        elem = type_map.get(lmp_type, f"X{lmp_type}")
        fout.write(f"{idx:>4d}  {elem:>4s}  "
                   f"{x:12.8f}  {y:12.8f}  {z:12.8f}  "
                   f"0.0  0.0  0.0  0.0  0.0  0.0\n")

    fout.write("\n")


# =============================================================================
# Task 3108: XDATCAR → ABACUS MD_dump
# =============================================================================

def _parse_xdatcar_header(filepath: Path) -> dict:
    """Read XDATCAR header: lattice, species, atom counts, coord type."""
    with open(filepath) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])

    # Line 5: species names or atom counts
    line5_parts = lines[5].split()
    try:
        counts = [int(x) for x in line5_parts]
        species = [f"X{i+1}" for i in range(len(counts))]
        header_end = 7
    except ValueError:
        species = line5_parts
        counts = [int(x) for x in lines[6].split()]
        header_end = 8

    # Handle both "Direct" and "Direct configuration= 1"
    coord_raw = lines[header_end - 1].lower()
    coord_type = "Direct" if coord_raw.startswith("direct") else "Cartesian"
    total_atoms = sum(counts)

    return {
        "scale": scale,
        "cell": cell * scale,
        "species": species,
        "counts": counts,
        "coord_type": coord_type,
        "total_atoms": total_atoms,
        "header_lines": header_end,
        "all_lines": lines,
    }


@task(3108, category="MD Analysis", name="XDATCAR→MD_dump",
      description="Convert VASP XDATCAR trajectory to ABACUS MD_dump format")
def task_xdatcar_to_md_dump(args: list[str] | None = None, interactive: bool = True) -> None:
    """Convert VASP XDATCAR to ABACUS MD_dump format.

    Supports both NVT (fixed cell) and NPT (per-config cell) XDATCAR formats.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== XDATCAR → ABACUS MD_dump ===[/bold cyan]")
    console.print()

    xdatcar_path = None
    if args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                xdatcar_path = p
                break
    if interactive and xdatcar_path is None:
        inp = _prompt(console, "XDATCAR file path", "XDATCAR")
        if inp:
            xdatcar_path = Path(inp)

    if xdatcar_path is None or not Path(xdatcar_path).exists():
        console.print("[red]XDATCAR not found.[/red]")
        return

    hdr = _parse_xdatcar_header(xdatcar_path)
    console.print(f"  Atoms: {hdr['total_atoms']}")
    console.print(f"  Species: {', '.join(f'{s}({c})' for s, c in zip(hdr['species'], hdr['counts']))}")
    console.print(f"  Coord type: {hdr['coord_type']}")

    # Detect NPT: look for scale + lattice lines before second "Direct" header
    config_lines = hdr["all_lines"][hdr["header_lines"]:]
    atoms_per_config = hdr["total_atoms"]
    all_lines = hdr["all_lines"]

    # Detect NPT: after first config's atoms, if the next "Direct/Cartesian"
    # header is more than 1 line away, it's NPT with per-config sub-headers.
    # NPT sub-header: comment, scale, lattice×3, species, counts = 7 lines.
    is_npt = False
    probe_start = hdr["header_lines"] + atoms_per_config
    for i in range(probe_start, min(len(all_lines), probe_start + 10)):
        if all_lines[i].lower().startswith(("direct", "cartesian")):
            is_npt = (i - probe_start > 1)  # >1 gap means sub-header exists
            break

    n_configs = count_xdatcar_configs(all_lines, hdr["header_lines"], atoms_per_config, is_npt)
    console.print(f"  Configurations: {n_configs}")
    console.print(f"  Ensemble: {'NPT (per-config cell)' if is_npt else 'NVT (fixed cell)'}")
    console.print()

    if n_configs == 0:
        console.print("[red]No configurations found in XDATCAR.[/red]")
        return

    atom_labels = []
    for s, c in zip(hdr["species"], hdr["counts"]):
        atom_labels.extend([s] * c)

    out_path = Path(str(xdatcar_path) + ".ABACUS.dump")

    if is_npt:
        _convert_xdatcar_npt(all_lines, hdr, atom_labels, atoms_per_config,
                             n_configs, out_path, console)
    else:
        _convert_xdatcar_nvt(config_lines, hdr, atom_labels, atoms_per_config,
                             n_configs, out_path, console)


def count_xdatcar_configs(all_lines, header_end, atoms_per_config, is_npt):
    """Count configurations in XDATCAR."""
    count = 1  # first config is always in the header
    remaining = all_lines[header_end + atoms_per_config:] if is_npt else all_lines[header_end:]
    count += sum(1 for l in remaining
                 if l.lower().startswith(("direct", "cartesian")))
    return count


def _convert_xdatcar_nvt(config_lines, hdr, atom_labels, atoms_per_config,
                          n_configs, out_path, console):
    """Convert NVT XDATCAR (fixed cell for all configs)."""
    from rich.progress import Progress

    line_idx = 0
    config_count = 0
    cell = hdr["cell"]

    with open(out_path, "w") as fout, Progress() as progress:
        task = progress.add_task("[cyan]Converting...", total=n_configs)

        while config_count < n_configs and line_idx < len(config_lines):
            if line_idx < len(config_lines) and config_lines[line_idx].lower().startswith(
                ("direct", "cartesian")):
                line_idx += 1

            if line_idx + atoms_per_config > len(config_lines):
                break

            _write_md_frame_nvt(fout, config_count + 1, cell, config_lines,
                                line_idx, atoms_per_config, hdr, atom_labels)
            line_idx += atoms_per_config
            config_count += 1
            progress.update(task, advance=1,
                            description=f"[cyan]Converting... config {config_count}")

    console.print()
    console.print(f"[green]✓ Converted {config_count} configurations → {out_path}[/green]")
    console.print(f"  Cell: {cell[0,0]:.4f} × {cell[1,1]:.4f} × {cell[2,2]:.4f} Å³")
    console.print()


def _convert_xdatcar_npt(all_lines, hdr, atom_labels, atoms_per_config,
                          n_configs, out_path, console):
    """Convert NPT XDATCAR (per-config cell vectors)."""
    from rich.progress import Progress

    config_count = 0
    i = hdr["header_lines"]

    with open(out_path, "w") as fout, Progress() as progress:
        task = progress.add_task("[cyan]Converting...", total=n_configs)

        # First config: positions are at header_end, cell from first header
        _write_md_frame_nvt(fout, 1, hdr["cell"], all_lines,
                            i, atoms_per_config, hdr, atom_labels)
        i += atoms_per_config
        config_count = 1
        progress.update(task, advance=1, description="[cyan]Converting... config 1")

        # Remaining configs: each has a NPT sub-header
        while config_count < n_configs and i < len(all_lines):
            # Skip NPT sub-header (7 lines) to find next "Direct/Cartesian"
            while i < len(all_lines) and not all_lines[i].lower().startswith(
                ("direct", "cartesian")):
                i += 1
            if i >= len(all_lines):
                break

            # Extract per-config cell — the sub-header ends at line i
            # Layout: i-7=comment, i-6=scale, i-5..i-3=lattice, i-2=species, i-1=counts, i=Direct
            if i >= 7:
                try:
                    scale = float(all_lines[i - 6])
                    lv = np.array([[float(x) for x in all_lines[i - 5 + j].split()]
                                    for j in range(3)])
                    cell = lv * scale
                except (ValueError, IndexError):
                    cell = hdr["cell"]
            else:
                cell = hdr["cell"]

            i += 1  # skip coord header

            if i + atoms_per_config > len(all_lines):
                break

            _write_md_frame_nvt(fout, config_count + 1, cell, all_lines,
                                i, atoms_per_config, hdr, atom_labels)
            i += atoms_per_config
            config_count += 1
            progress.update(task, advance=1,
                            description=f"[cyan]Converting... config {config_count}")

    console.print()
    console.print(f"[green]✓ Converted {config_count} configurations → {out_path}[/green]")
    console.print()


def _write_md_frame_nvt(fout, step, cell, lines, start, natoms, hdr, labels):
    """Write one MD_dump frame from XDATCAR atom data."""
    fout.write(f"MDSTEP: {step}\n")
    fout.write("LATTICE_CONSTANT: 1.0 Angstrom\n")
    fout.write("LATTICE_VECTORS\n")
    for row in cell:
        fout.write(f"   {row[0]:15.10f}  {row[1]:15.10f}  {row[2]:15.10f}\n")
    fout.write("\n")
    fout.write("INDEX   LABEL   X           Y           Z           "
               "FX   FY   FZ   VX   VY   VZ\n")

    for j in range(natoms):
        if start + j >= len(lines):
            break
        parts = lines[start + j].split()
        if len(parts) < 3:
            continue
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])

        if hdr["coord_type"].lower().startswith("direct"):
            frac = np.array([x, y, z])
            cart = frac @ cell
            x, y, z = cart[0], cart[1], cart[2]

        label = labels[j] if j < len(labels) else "X"
        fout.write(f"{j+1:>4d}  {label:>4s}  "
                   f"{x:12.8f}  {y:12.8f}  {z:12.8f}  "
                   f"0.0  0.0  0.0  0.0  0.0  0.0\n")

    fout.write("\n")


# =============================================================================
# Task 3109: Energy & Temperature vs Time
# =============================================================================

@task(3109, category="MD Analysis", name="E-T vs Time",
      description="Extract energy and temperature vs time from running_md.log, plot and save")
def task_energy_temp_vs_time(args: list[str] | None = None, interactive: bool = True) -> None:
    """Read MD step summaries from running_md.log, output energy & temperature
    as a function of time, and produce publication-quality plots.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Energy & Temperature vs Time ===[/bold cyan]")
    console.print()

    from abacuscopilot.preprocessing.system_tasks import _find_md_log, _parse_md_progress
    log_path = _find_md_log()
    if log_path is None:
        console.print("[red]No running_md.log found (OUT.ABACUS/running_md.log).[/red]")
        return

    console.print(f"  [dim]Log: {log_path}[/dim]")

    data = _parse_md_progress(log_path, tail=0)
    if not data:
        console.print("[red]No MD steps found in log.[/red]")
        return

    md_dt = _read_md_dt(str(log_path.parent))
    try:
        md_dt = float(md_dt)
    except (ValueError, TypeError):
        console.print("[red]Cannot read md_dt from INPUT.[/red]")
        return

    steps = np.array([d["step"] for d in data])
    times = steps * md_dt / 1000.0  # fs -> ps
    from abacuscopilot.core.constants import RY_TO_EV
    energies_ry = np.array([d["energy_ry"] for d in data])
    energies = energies_ry * RY_TO_EV  # Ry → eV
    temperatures = np.array([d["temperature_k"] for d in data])

    console.print(f"  Steps: {steps[0]} -> {steps[-1]} ({len(steps)} points)")
    console.print(f"  Time:  {times[0]:.2f} -> {times[-1]:.2f} ps  (md_dt={md_dt} fs)")
    console.print(f"  T:     {temperatures.min():.1f} -> {temperatures.max():.1f} K")

    out_dir = log_path.parent  # OUT.ABACUS
    out_dat = str(out_dir / "energy_temperature.dat")
    with open(out_dat, "w") as f:
        f.write("# time(ps)  energy(eV)  temperature(K)\n")
        for t, e, temp in zip(times, energies, temperatures):
            f.write(f"{t:.6f}  {e:.8f}  {temp:.4f}\n")
    console.print(f"  [green]... Data: {out_dat}[/green]")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    if interactive:
        mode = _prompt_choice(
            console, "Plot style",
            ["Separate panels (E and T)", "Dual-Y axes (E+T on one plot)"],
            "Separate panels (E and T)",
        )
    else:
        mode = "Separate"

    out_png = str(out_dir / "energy_temperature.png")

    # --- Y-axis limits ---
    t_max_lim = max(temperatures[0], temperatures[-1]) * 1.8
    e_min, e_max = energies.min(), energies.max()
    e_range = e_max - e_min if e_max > e_min else abs(e_max) * 0.01
    e_pad = e_range * 0.15
    e_lo, e_hi = e_min - e_pad, e_max + e_pad

    if "Dual" in mode:
        fig, ax1 = plt.subplots(figsize=(8, 6))
        ax2 = ax1.twinx()
        ax1.plot(times, energies, "-", color="#1f77b4", linewidth=1.2, label="Energy (eV)")
        ax2.plot(times, temperatures, "-", color="#d62728", linewidth=1.2, label="Temperature (K)")
        ax1.set_xlabel("Time (ps)")
        ax1.set_ylabel("Energy (eV)", color="#1f77b4")
        ax2.set_ylabel("Temperature (K)", color="#d62728")
        ax1.set_ylim(e_lo, e_hi)
        ax2.set_ylim(0, t_max_lim)
        ax1.tick_params(axis="y", labelcolor="#1f77b4", direction="out")
        ax2.tick_params(axis="y", labelcolor="#d62728", direction="out")
        ax1.tick_params(axis="x", direction="out")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
        for spine in ax1.spines.values():
            spine.set_linewidth(0.5)
        ax1.spines["top"].set_visible(True)
        ax1.spines["right"].set_visible(True)
        for spine in ax2.spines.values():
            spine.set_linewidth(0.5)
    else:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
        ax1.plot(times, energies, "-", color="#1f77b4", linewidth=1.2)
        ax1.set_ylabel("Energy (eV)")
        ax1.set_title("MD - Energy & Temperature vs Time")
        ax1.set_ylim(e_lo, e_hi)
        ax2.plot(times, temperatures, "-", color="#d62728", linewidth=1.2)
        ax2.set_xlabel("Time (ps)")
        ax2.set_ylabel("Temperature (K)")
        ax2.set_ylim(0, t_max_lim)
        for ax in (ax1, ax2):
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
            ax.tick_params(axis="both", direction="out")

    fig.tight_layout(pad=1.2)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    console.print(f"  [green]... Plot: {out_png}[/green]")
    console.print()
