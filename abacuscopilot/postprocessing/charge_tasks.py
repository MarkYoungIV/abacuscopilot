"""Charge density analysis tasks for ABACUS output.

Task IDs 1001-1003 (Charge Density) and 1304-1305 (Population, Bader/Hirshfeld)

Reads charge density files (CHG*, SPIN*) from ABACUS output, computes
1D planar averages, provides Cube/XSF export, difference charge density,
and atomic charge analysis (Bader / Hirshfeld live here because they read
the charge density; they are registered in the Population menu).
"""

from __future__ import annotations

import re
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console
from abacuscopilot.core.constants import ANGSTROM_TO_BOHR
from abacuscopilot.io.stru_file import read_stru
from abacuscopilot.tasks import task

# =============================================================================
# Charge density file reader
# =============================================================================


def _find_charge_file(pattern: str = "CHG") -> Path | None:
    """Find a charge density file in the working directory."""
    candidates = (list(Path().glob(f"{pattern}*")) +
                  list(Path().glob(f"OUT.*/{pattern}*")) +
                  list(Path().glob(f"OUT.*/*/{pattern}*")))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _read_str_u() -> tuple[np.ndarray, int, tuple, np.ndarray] | None:
    """Read STRU file to get cell and atomic information needed by charge routines."""
    try:
        from abacuscopilot.io.stru_file import read_stru
    except ImportError:
        return None
    stru_path = Path("STRU")
    if not stru_path.exists():
        return None
    s = read_stru(stru_path)
    return s.lattice.cell_angstrom, s.num_atoms, s.species_order, s.positions


# -----------------------------------------------------------------------------
# ABACUS binary rhog restart reader (e.g. `OUT.*/ABACUS-CHARGE-DENSITY.restart`,
# commonly copied out and renamed `CHG` / `chg1` etc.)
#
# Format (ModuleIO::write_rhog, source/module_io/rhog_io.cpp), all little-endian:
#   [int 3][int gamma_only][int ngm][int nspin][int 3]
#   [int 9][double GT[9]][int 9]        # GT = reciprocal lattice in 1/Angstrom
#   [int 3*ngm][int miller[3*ngm]][int 3*ngm]   # G-vectors as Miller indices
#   per spin: [int ngm][complex<double> rho[ngm]][int ngm]
# The FFT grid is not stored; ABACUS picks FFT-friendly (2,3,5,7-only factor)
# sizes, which we reconstruct from the miller-index extents.
# -----------------------------------------------------------------------------


def _next_smooth(n: int) -> int:
    """Smallest integer >= n whose prime factors are only 2, 3, 5, 7 (FFTW-friendly)."""
    n = int(n)
    while True:
        m = n
        for p in (2, 3, 5, 7):
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def read_charge_rhog(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read an ABACUS binary rhog (charge density) restart file.

    Returns (data_3d, cell_bohr, origin) with data ordered (nz, ny, nx),
    matching the Gaussian-Cube convention, or None if the file is not a valid
    rhog restart. `cell_bohr` is the 3x3 real-space lattice (rows = a, b, c).
    """
    filepath = Path(filepath)
    raw = filepath.read_bytes()

    def rd_int(pos):
        if pos + 4 > len(raw):
            return None, pos
        return struct.unpack_from("<i", raw, pos)[0], pos + 4

    def rd_dbl(pos):
        if pos + 8 > len(raw):
            return None, pos
        return struct.unpack_from("<d", raw, pos)[0], pos + 8

    pos = 0
    # Header: [3, gamma_only, ngm, nspin, 3]
    hdr = []
    for _ in range(5):
        v, pos = rd_int(pos)
        if v is None:
            return None
        hdr.append(v)
    if hdr[0] != 3 or hdr[4] != 3 or hdr[2] <= 0:
        return None
    ngm, nspin = hdr[2], hdr[3]

    # Reciprocal lattice GT (1/Angstrom). GT = (latvec_Bohr)^{-1} * ang_to_bohr.
    marker, pos = rd_int(pos)
    if marker != 9:
        return None
    gt = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            gt[i, j], pos = rd_dbl(pos)
    marker, pos = rd_int(pos)
    if marker != 9:
        return None
    cell_bohr = np.linalg.inv(gt) * ANGSTROM_TO_BOHR  # real-space lattice, Bohr

    # Miller indices (one int per G-vector component)
    nm, pos = rd_int(pos)
    if nm <= 0 or nm % 3 != 0 or pos + 4 * nm > len(raw):
        return None
    miller = np.frombuffer(raw, dtype="<i4", count=nm, offset=pos).reshape(-1, 3)
    pos += 4 * nm
    marker, pos = rd_int(pos)

    # Reconstruct the FFT grid (ABACUS uses smooth sizes; grid must hold the
    # full planewave set without aliasing).
    ix, iy, iz = miller[:, 0], miller[:, 1], miller[:, 2]
    nx = _next_smooth(2 * max(abs(ix.min()), abs(ix.max())) + 1)
    ny = _next_smooth(2 * max(abs(iy.min()), abs(iy.max())) + 1)
    nz = _next_smooth(2 * max(abs(iz.min()), abs(iz.max())) + 1)

    # rho(G) per spin channel (complex<double>), summed to total real-space rho
    rho_r = np.zeros((nx, ny, nz), dtype=np.float64)
    nvec = len(miller)
    gx = np.where(ix < 0, ix + nx, ix)
    gy = np.where(iy < 0, iy + ny, iy)
    gz = np.where(iz < 0, iz + nz, iz)
    for _ in range(nspin):
        nm, pos = rd_int(pos)
        if nm != nvec or pos + 16 * nvec > len(raw):
            return None
        rhog = np.frombuffer(raw, dtype=np.complex128, count=nvec, offset=pos)
        pos += 16 * nvec
        marker, pos = rd_int(pos)
        arr = np.zeros((nx, ny, nz), dtype=np.complex128)
        arr[gx, gy, gz] = rhog
        # rho(r) = sum_G rho(G) exp(i G.r); numpy ifftn already carries the 1/N
        rho_r += np.fft.ifftn(arr).real * (nx * ny * nz)

    # Transpose (nx, ny, nz) -> (nz, ny, nx) to match the Cube convention
    return rho_r.transpose(2, 1, 0), cell_bohr, np.zeros(3)


def read_charge_cube(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read a charge density file — Gaussian Cube text or ABACUS binary rhog.

    Returns (data_3d, cell_bohr, origin) or None if the file is unreadable.
    data_3d is ordered (nz, ny, nx); cell_bohr rows are the real-space lattice.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None
    try:
        with open(filepath, encoding="utf-8") as f:
            header = f.readline().strip()
    except UnicodeDecodeError:
        # Binary file — try the ABACUS rhog restart format
        return read_charge_rhog(filepath)
    if "CUBE" not in header.upper() and "CUBE" not in filepath.name.upper():
        return None

    # Simple Cube reader for ABACUS-generated files
    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip()]

    # Skip comment lines
    data_start = 2
    # atoms + origin line
    parts = lines[data_start].split()
    natom = int(parts[0])
    origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])

    # Three grid-spacing lines follow (line 3 = nx, then ny, nz).
    data_start += 1
    dims = []
    cell = np.zeros((3, 3))
    for i in range(3):
        parts = lines[data_start + i].split()
        n = int(parts[0])
        dims.append(n)
        cell[i] = np.array([float(parts[1]) * n, float(parts[2]) * n, float(parts[3]) * n])
    nx, ny, nz = dims

    # Skip the natom atom-position lines, then read the volumetric data.
    data_start += 3 + natom
    raw = []
    for line in lines[data_start:]:
        raw.extend(float(x) for x in line.split())

    expected = nx * ny * nz
    data = np.array(raw[:expected])
    # ABACUS CUBE writes z-fastest ("Inner loop is z ..."); standard Gaussian is
    # x-fastest.  Transpose so data is (nz, ny, nx) like the copilot's own cubes.
    if "inner loop is z" in header.lower():
        data = data.reshape(nx, ny, nz).transpose(2, 1, 0)
    else:
        data = data.reshape(nz, ny, nx)

    # ABACUS CUBE is already in Bohr / e-Bohr³ (matches the rhog reader).
    return data, cell, origin


# =============================================================================
# 1D planar average
# =============================================================================


def planar_average_1d(data: np.ndarray, axis: int = 2) -> np.ndarray:
    """Compute 1D planar average along a given axis.

    Args:
        data: 3D numpy array.
        axis: Axis to average over (0=x, 1=y, 2=z).

    Returns:
        1D array of averaged values.
    """
    axes = tuple(i for i in range(3) if i != axis)
    return data.mean(axis=axes)


# =============================================================================
# Export to Cube / XSF
# =============================================================================


def write_cube(
    filepath: str | Path,
    data: np.ndarray,
    cell: np.ndarray,
    origin: tuple = (0.0, 0.0, 0.0),
    title: str = "abacuscopilot charge density",
) -> None:
    """Write a 3D scalar field in Gaussian Cube format.

    Args:
        filepath: Output path.
        data: 3D array (nz, ny, nx) as in Cube convention.
        cell: 3x3 cell matrix in Bohr.
        origin: Origin of the grid.
        title: Header title line.
    """
    nz, ny, nx = data.shape
    with open(filepath, "w") as f:
        f.write(f"{title}\n")
        f.write("abacuscopilot Cube export\n")
        f.write(f"  1  {origin[0]:.6f}  {origin[1]:.6f}  {origin[2]:.6f}\n")
        # ABACUS CUBE convention: cell axes in Bohr, density in e/Bohr³.
        for i, n in enumerate((nx, ny, nz)):
            f.write(f"  {n}  {cell[i,0]/n:.10f}"
                    f"  {cell[i,1]/n:.10f}"
                    f"  {cell[i,2]/n:.10f}\n")
        f.write("  1  1.0  0.0  0.0  0.0\n")  # dummy atom

        for iz in range(nz):
            for iy in range(ny):
                for ix in range(0, nx, 6):
                    chunk = data[iz, iy, ix:min(ix+6, nx)]
                    f.write("  " + "  ".join(f"{v:.6e}" for v in chunk) + "\n")


def write_xsf(
    filepath: str | Path,
    data: np.ndarray,
    structure,
    title: str = "abacuscopilot charge density",
) -> None:
    """Write a 3D scalar field in XCrySDen XSF format.

    Args:
        filepath: Output path.
        data: 3D array (nx, ny, nz) in XSF convention (x-fastest).
        structure: Structure object for atomic positions.
        title: Header title.
    """
    nx, ny, nz = data.shape
    cell = structure.lattice.cell_angstrom  # Å for XSF

    with open(filepath, "w") as f:
        f.write("CRYSTAL\nPRIMVEC\n")
        for v in cell:
            f.write(f"  {v[0]:.10f}  {v[1]:.10f}  {v[2]:.10f}\n")
        f.write(f"PRIMCOORD\n  {structure.num_atoms}  1\n")
        for atom in structure.atoms:
            f.write(f"  {atom.species:>3s}  {atom.position[0]:.10f}  "
                    f"{atom.position[1]:.10f}  {atom.position[2]:.10f}\n")

        f.write(f"BEGIN_BLOCK_DATAGRID_3D\n  {title}\n")
        f.write("  BEGIN_DATAGRID_3D\n")
        f.write(f"    {nx} {ny} {nz}\n")
        f.write("    0.0 0.0 0.0\n")
        for v in cell:
            f.write(f"    {v[0]:.10f}  {v[1]:.10f}  {v[2]:.10f}\n")

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(0, nz, 6):
                    chunk = data[ix, iy, iz:min(iz+6, nz)]
                    f.write("  " + "  ".join(f"{v:.6e}" for v in chunk) + "\n")

        f.write("  END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D\n")


# =============================================================================
# Task 741: 1D planar average charge density
# =============================================================================


@task(1001, category="Charge Density", name="1D Planar Avg Charge",
      description="Compute 1D planar-averaged charge density and save plot",
      cli_args=[
          {"name": "--file", "type": str, "default": None, "help": "Path to CHG cube file"},
      ])
def task_planar_avg_charge(args: list[str] | None = None, interactive: bool = True,
                           parsed_args=None) -> None:
    """Compute and plot the 1D planar-averaged charge density."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== 1D Planar-Averaged Charge Density ===[/bold cyan]")
    console.print()

    # Find charge file
    chg_path = None
    if parsed_args and parsed_args.file:
        chg_path = Path(parsed_args.file)
    if chg_path is None and args:
        for arg in args:
            if Path(arg).exists():
                chg_path = Path(arg)
                break
    if chg_path is None:
        chg_path = _find_charge_file("CHG")

    if chg_path is None:
        console.print("[red]No CHG file found. Run an ABACUS calculation with out_chg=True first.[/red]")
        return

    console.print(f"  [dim]Charge file: {chg_path}[/dim]")

    # Try reading as Cube
    result = read_charge_cube(chg_path)
    if result is None:
        console.print("[yellow]Could not read charge density as Cube format.[/yellow]")
        console.print("[dim]ABACUS cube-format CHG files are in OUT.ABACUS/ with .cube extension.[/dim]")
        return

    data, cell, origin = result
    nz, ny, nx = data.shape
    console.print(f"  [bold]Grid:[/bold] {nx} × {ny} × {nz}")

    # Default: average along z (axis 2 in z,y,x → axis 0 in nz,ny,nx)
    axis = 0  # z-direction in (nz, ny, nx) ordering
    avg = planar_average_1d(data, axis=axis).flatten()
    npts = avg.shape[0]
    z_axis = origin[2] + np.linspace(0, cell[2, 2], npts, endpoint=False)  # in Bohr
    from abacuscopilot.core.constants import BOHR_TO_ANGSTROM
    z_axis_ang = z_axis * BOHR_TO_ANGSTROM

    console.print(f"  [bold]Average:[/bold] {avg.mean():.4f} e/Å³")

    # Plot
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_axis_ang, avg, color="#1f77b4", linewidth=1.5)
    ax.set_xlabel("z (Å)")
    ax.set_ylabel("ρ (e/Å³)")
    ax.set_title("1D Planar-Averaged Charge Density")
    ax.set_xlim(z_axis_ang.min(), z_axis_ang.max())

    # Save data
    data_file = "planar_avg_chg.dat"
    np.savetxt(data_file, np.column_stack([z_axis_ang, avg]),
               fmt="%.8f", header="z_Angstrom  rho_e_per_A3")
    console.print(f"  [dim]Data saved to {data_file}[/dim]")

    save_name = "planar_avg_chg.png"
    fig.savefig(save_name)
    console.print(f"  [green]✓ Plot saved to {save_name}[/green]")
    plt.close(fig)
    console.print()


# =============================================================================
# Task 742: Export charge density to Cube/XSF
# =============================================================================


@task(1002, category="Charge Density", name="Export Cube/XSF",
      description="Export charge density to Cube or XCrySDen XSF format for VESTA")
def task_export_chg(args: list[str] | None = None, interactive: bool = True) -> None:
    """Export charge density to Cube or XSF format for visualization in VESTA."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Export Charge Density ===[/bold cyan]")
    console.print()

    chg_path = None
    if args:
        for arg in args:
            if Path(arg).exists():
                chg_path = Path(arg)
                break
    if chg_path is None:
        chg_path = _find_charge_file("CHG")

    if chg_path is None:
        console.print("[red]No CHG file found.[/red]")
        return

    console.print(f"  [dim]Charge file: {chg_path}[/dim]")

    result = read_charge_cube(chg_path)
    if result is None:
        console.print("[red]Could not read charge density. Ensure ABACUS cube output is enabled.[/red]")
        return

    data, cell, origin = result

    # Load STRU for XSF export
    from abacuscopilot.io.stru_file import read_stru
    structure = None
    if Path("STRU").exists():
        try:
            structure = read_stru("STRU")
        except Exception:
            pass

    # Ask format
    fmt = "cube"
    if interactive:
        from rich.prompt import Prompt
        fmt = Prompt.ask("  Export format", choices=["cube", "xsf"], default="cube")
    else:
        fmt = "xsf" if structure else "cube"

    if fmt == "cube":
        write_cube("charge_density.cube", data, cell)
        console.print("  [green]✓ Cube file written: charge_density.cube[/green]")
    else:
        if structure is None:
            console.print("[yellow]STRU file not found — needed for XSF format. Falling back to Cube.[/yellow]")
            write_cube("charge_density.cube", data, cell)
            console.print("  [green]✓ Cube file written: charge_density.cube[/green]")
        else:
            # Transpose to (nx, ny, nz) for XSF
            data_xsf = data.transpose(2, 1, 0)
            write_xsf("charge_density.xsf", data_xsf, structure)
            console.print("  [green]✓ XSF file written: charge_density.xsf[/green]")

    console.print()


# =============================================================================
# Task 743: Difference charge density
# =============================================================================


@task(1003, category="Charge Density", name="Diff. Charge Density",
      description="Compute difference charge density Δρ = ρ_AB − ρ_A − ρ_B")
def task_diff_charge(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compute difference charge density from three charge density files.

    Δρ = ρ(AB) − ρ(A) − ρ(B)

    Useful for analyzing charge transfer at interfaces/adsorbates.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Difference Charge Density ===[/bold cyan]")
    console.print()
    console.print("[dim]Δρ = ρ(AB) − ρ(A) − ρ(B)[/dim]")
    console.print()

    if interactive:
        chg_ab = _prompt_path(console, "Path to ρ(AB) (total system)")
        chg_a = _prompt_path(console, "Path to ρ(A) (component A)")
        chg_b = _prompt_path(console, "Path to ρ(B) (component B)")
    else:
        console.print("[yellow]Please provide 3 CHG files as arguments: chg_ab chg_a chg_b[/yellow]")
        return

    # Read all three. Each file may be a Gaussian Cube text file or an ABACUS
    # binary rhog restart (e.g. OUT.*/ABACUS-CHARGE-DENSITY.restart, or copies
    # commonly renamed `chg1`/`chg2`/...). Auto-detected by read_charge_cube.
    datasets = []
    cell = None
    for label, fp in [("AB", chg_ab), ("A", chg_a), ("B", chg_b)]:
        if not fp or not Path(fp).exists():
            console.print(f"[red]File not found for ρ({label}): {fp}[/red]")
            return
        r = read_charge_cube(fp)
        if r is None:
            console.print(f"[red]Could not read ρ({label}) from {fp}[/red]")
            console.print("[dim]Expected a Gaussian Cube file or an ABACUS binary CHG restart file.[/dim]")
            return
        datasets.append(r[0])  # data only
        if cell is None:
            cell = r[1]  # real-space lattice (Bohr) for the cube header

    d_ab, d_a, d_b = datasets

    # Grids must match for element-wise subtraction. If they differ (e.g. the
    # three runs used different cells), resample onto the ρ(AB) grid.
    if d_ab.shape != d_a.shape or d_ab.shape != d_b.shape:
        console.print("[yellow]Grid dimensions differ between the three files.[/yellow]")
        console.print(f"[yellow]  ρ(AB) {d_ab.shape}  ρ(A) {d_a.shape}  ρ(B) {d_b.shape}[/yellow]")
        console.print("[yellow]  Resampling ρ(A) and ρ(B) onto the ρ(AB) grid...[/yellow]")
        from scipy.ndimage import zoom
        factors = [t / s for t, s in zip(d_ab.shape, d_a.shape)]
        d_a = zoom(d_a, factors, order=1)
        factors = [t / s for t, s in zip(d_ab.shape, d_b.shape)]
        d_b = zoom(d_b, factors, order=1)

    diff = d_ab - d_a - d_b

    write_cube("diff_charge_density.cube", diff, cell)

    console.print()
    console.print("[green]✓ Difference charge density written: diff_charge_density.cube[/green]")
    console.print("  Open with VESTA/Jmol for visualization.")
    console.print()


def _prompt_path(console, question: str) -> str:
    """Prompt for a file path."""
    from rich.prompt import Prompt
    return Prompt.ask(f"  {question}")


# =============================================================================
# Atomic charge analysis: Mulliken / Bader / Hirshfeld
# =============================================================================


def _atom_labels(structure) -> list[str]:
    """Per-atom labels in STRU order, GLOBALLY numbered (C1..C24, H25..H36).

    Matches the Mulliken task's numbering (atom index = STRU position), so
    the Bader / RESP tables line up atom-for-atom with the 1301 population
    table.  E.g. 24 C then 12 H → C1..C24, H25..H36.
    """
    return [f"{a.species}{i + 1}" for i, a in enumerate(structure.atoms)]


def _locate_bader() -> Path | None:
    """Locate the bundled (or system) bader.x binary."""
    pkg_root = Path(__file__).resolve().parent.parent.parent
    candidates = [
        pkg_root / "scripts" / "bader" / "bader.x",
        Path("scripts") / "bader" / "bader.x",
    ]
    for name in ("bader.x", "bader"):
        p = shutil.which(name)
        if p:
            return Path(p)
    for c in candidates:
        if c.exists():
            return c
    return None


def _read_upf_zval(pp_name: str) -> float | None:
    """Read z_valence (valence electron count) from a UPF file.

    Searches the current directory, OUT.* and the configured pseudo library.
    """
    if not pp_name:
        return None
    from abacuscopilot.config import load_config
    libs = load_config().get("libraries", {})
    search = [Path(".")] + [Path(d) for d in (libs.get("pseudo_library") or [])]
    for d in search:
        for p in [d / pp_name] + list(d.glob(f"**/{pp_name}")):
            if p.is_file():
                try:
                    content = p.read_text(errors="ignore")
                except Exception:
                    continue
                m = re.search(r'z_valence\s*=\s*"\s*([\d.]+)', content)
                if m:
                    return float(m.group(1))
    return None


_DEFAULT_ZVAL = {
    "H": 1, "He": 2, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8,
    "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12, "Ga": 3, "Ge": 4, "As": 5,
    "Se": 6, "Br": 7, "Kr": 8,
    "Rb": 1, "Sr": 2, "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8,
    "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12, "In": 3, "Sn": 4, "Sb": 5,
    "Te": 6, "I": 7, "Xe": 8,
    "Cs": 1, "Ba": 2, "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7,
    "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12, "Tl": 3, "Pb": 4,
    "Bi": 5, "Sm": 8, "Ce": 8, "Yb": 8, "Lu": 9,
}


def _zval_for(species: str, structure) -> float:
    """Valence electrons for a species — from the actual UPF when possible."""
    pp = structure.pseudo_files.get(species, "")
    z = _read_upf_zval(pp)
    if z is not None:
        return z
    return float(_DEFAULT_ZVAL.get(species, 0.0))


def _write_bader_cube(filepath, data, cell_bohr, origin, structure) -> None:
    """Write a Gaussian-Cube density with real atom positions (STRU order).

    The Bader program reports ACF.dat rows in this exact atom order, which is
    how we map Bader charges back to atoms.  Atom positions are in Bohr.
    """
    try:
        from ase.data import atomic_numbers
    except ImportError:
        atomic_numbers = {}
    nz, ny, nx = data.shape
    dims = [nx, ny, nz]
    atom_lines = []
    cart = _atom_cartesian_angstrom(structure)
    for i, a in enumerate(structure.atoms):
        pos = cart[i] * ANGSTROM_TO_BOHR  # ABACUS CUBE atom positions are in Bohr
        znum = int(atomic_numbers.get(a.species, 1))
        atom_lines.append(
            f"  {znum:3d}  {0.0:10.6f}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}"
        )
    with open(filepath, "w") as f:
        f.write("abacuscopilot Bader input\n")
        f.write("charge density\n")
        f.write(f"  {len(atom_lines):5d}  {origin[0]:12.6f}  {origin[1]:12.6f}  {origin[2]:12.6f}\n")
        # ABACUS CUBE convention: cell axes in Bohr, density in e/Bohr³.
        for i, n in enumerate(dims):
            f.write(f"  {n:5d}  {cell_bohr[i, 0] / n:12.6f}"
                    f"  {cell_bohr[i, 1] / n:12.6f}"
                    f"  {cell_bohr[i, 2] / n:12.6f}\n")
        for line in atom_lines:
            f.write(line + "\n")
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(0, nx, 6):
                    chunk = data[iz, iy, ix:ix + 6]
                    f.write("  " + "  ".join(f"{v:.6e}" for v in chunk) + "\n")


def _read_cube_standard(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read a standard Gaussian CUBE file (the format VESTA/Multiwfn/bader use).

    Returns (data_3d (nz, ny, nx), cell_bohr, origin) or None on failure.
    """
    try:
        with open(filepath) as f:
            lines = [l.strip() for l in f if l.strip()]
    except Exception:
        return None
    try:
        parts = lines[2].split()
        natom = int(parts[0])
        origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        dims = []
        cell = np.zeros((3, 3))
        for i in range(3):
            parts = lines[3 + i].split()
            n = int(parts[0])
            dims.append(n)
            cell[i] = np.array([float(parts[1]) * n, float(parts[2]) * n, float(parts[3]) * n])
        nx, ny, nz = dims
        start = 3 + 3 + natom
        raw = []
        for line in lines[start:]:
            raw.extend(float(x) for x in line.split())
        data = np.array(raw[: nx * ny * nz])
        # ABACUS writes CUBE with the INNER (fastest) loop = z — the header says
        # "Inner loop is z, followed by y and x" — which is the OPPOSITE of the
        # standard Gaussian x-fastest order.  Detect it and transpose so the
        # returned array is (nz, ny, nx) with x-fastest like the copilot's own
        # cubes.  (A mis-read here silently scrambles Bader/RESP/ESP — water
        # Bader gave O=8/H=0 before this fix.)
        if "inner loop is z" in lines[0].lower():
            data = data.reshape(nx, ny, nz).transpose(2, 1, 0)
        else:
            data = data.reshape(nz, ny, nx)
        # ABACUS writes SPIN1_CHG.cube in BOHR units (cell axes, atom positions
        # in Bohr, density in e/Bohr³) — already the code's internal convention,
        # matching the rhog reader.  No unit conversion.
        return data, cell, origin
    except Exception:
        return None


def _read_density_any(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read a charge density in standard-cube or ABACUS-binary-rhog format."""
    try:
        with open(filepath, encoding="utf-8") as f:
            f.readline()
    except UnicodeDecodeError:
        return read_charge_rhog(filepath)
    return _read_cube_standard(filepath)


def _parse_acf(acf_path: Path) -> list[tuple[float, float, float, float]]:
    """Parse ACF.dat into [(x, y, z, charge), ...] in atom order."""
    rows: list[tuple[float, float, float, float]] = []
    for line in acf_path.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) >= 7 and parts[0].lstrip("-").isdigit():
            rows.append((float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
    return rows


def _find_restart_file() -> Path | None:
    """Find ABACUS's always-written binary charge-density restart file."""
    for cand in (Path("ABACUS-CHARGE-DENSITY.restart"),
                 Path("OUT.ABACUS/ABACUS-CHARGE-DENSITY.restart")):
        if cand.exists():
            return cand
    for d in sorted(Path(".").glob("OUT.*")):
        p = d / "ABACUS-CHARGE-DENSITY.restart"
        if p.exists():
            return p
    return None


def _find_charge_density() -> tuple[np.ndarray | None, np.ndarray, np.ndarray, str]:
    """Find + read the total (spin-summed) charge density.

    Tries, in order: SPIN1/SPIN2 cubes, CHG/chg/CHGCAR copies, and finally the
    ABACUS-CHARGE-DENSITY.restart binary rhog that is written on every SCF even
    without ``out_chg``.

    Returns (data, cell_bohr, origin, description); data is None if not found.
    """
    spin1 = _find_charge_file("SPIN1_CHG")
    spin2 = _find_charge_file("SPIN2_CHG")
    if spin1 and spin2:
        d1, cell, origin = _read_density_any(spin1)
        d2, _, _ = _read_density_any(spin2)
        if d1 is not None and d2 is not None:
            return d1 + d2, cell, origin, f"SPIN1 + SPIN2 ({spin1.name}, {spin2.name})"
    if spin1:
        data, cell, origin = _read_density_any(spin1)
        if data is not None:
            return data, cell, origin, str(spin1)
    for pat in ("CHG", "chg", "CHGCAR"):
        fp = _find_charge_file(pat)
        if fp:
            data, cell, origin = _read_density_any(fp)
            if data is not None:
                return data, cell, origin, str(fp)
    rp = _find_restart_file()
    if rp:
        data, cell, origin = _read_density_any(rp)
        if data is not None:
            return data, cell, origin, str(rp)
    return None, None, None, ""


# HIDDEN (2026-08-26): Bader gives wrong charges — ABACUS SPIN1_CHG.cube is
# z-fastest ("Inner loop is z") but readers assumed x-fastest, scrambling the
# density; even after the transpose fix, bader.x still can't resolve the small
# H basins on the 192³ grid (water H got 0). Under development.
# @task(1304, category="Population", name="Bader Charge",
#       description="Bader atomic charges from the charge density (bundled bader.x)")
def task_bader_charge(args: list[str] | None = None, interactive: bool = True) -> None:
    """Compute Bader atomic charges by partitioning the total charge density.

    Reads the ABACUS charge density (cube or binary rhog), sums spins if
    spin-polarized, runs the bundled `bader.x`, parses ACF.dat and reports
    per-atom Bader charges (net charge = ZVAL − Bader electrons).
    """
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== Bader Charge ===[/bold cyan]")
    console.print()

    data, cell, origin, desc = _find_charge_density()
    if data is None:
        console.print("[red]Could not read a charge density file.[/red]")
        console.print("[dim]Looked for SPIN1_CHG.cube / SPIN2_CHG.cube, CHG, chg1, "
                      "ABACUS-CHARGE-DENSITY.restart.[/dim]")
        return
    console.print(f"  [dim]Density: {desc} ({data.shape[0]}×{data.shape[1]}×{data.shape[2]})[/dim]")

    if not Path("STRU").exists():
        console.print("[red]No STRU file — cannot map Bader charges to atoms.[/red]")
        return
    structure = read_stru("STRU")
    if not structure.atoms:
        console.print("[red]STRU contains no atoms.[/red]")
        return

    cube_path = "bader_input.cube"
    _write_bader_cube(cube_path, data, cell, origin, structure)

    bader = _locate_bader()
    if bader is None:
        console.print("[red]Bader program not found.[/red]")
        console.print("[dim]Run ./setup.sh to build it, or install the Henkelman bader code.[/dim]")
        return

    console.print("  Computing Bader partition...")
    try:
        subprocess.run([str(bader), cube_path], check=True, cwd=".",
                       capture_output=True, text=True, timeout=900)
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Bader failed: {(e.stderr or '').strip() or e}[/red]")
        return
    except FileNotFoundError:
        console.print(f"[red]Could not execute {bader}[/red]")
        return

    acf = Path("ACF.dat")
    if not acf.exists():
        console.print("[red]Bader did not produce ACF.dat.[/red]")
        return
    charges = _parse_acf(acf)
    n_atoms = len(structure.atoms)
    if len(charges) != n_atoms:
        console.print(f"[yellow]ACF.dat has {len(charges)} rows, STRU has {n_atoms} atoms — "
                      "checking the atom order/mapping.[/yellow]")

    labels = _atom_labels(structure)
    from rich.table import Table
    table = Table(title=f"Bader Charges ({n_atoms} atoms)")
    table.add_column("Atom")
    table.add_column("Element")
    table.add_column("Bader (e)")
    table.add_column("Net Charge (e)")

    rows: list[tuple[str, str, float, float]] = []
    total_net = 0.0
    for i, a in enumerate(structure.atoms):
        bader_e = charges[i][3] if i < len(charges) else float("nan")
        zval = _zval_for(a.species, structure)
        net = zval - bader_e
        total_net += net
        rows.append((labels[i], a.species, bader_e, net))
        table.add_row(labels[i], a.species, f"{bader_e:.4f}", f"{net:+.4f}")
    table.add_section()
    table.add_row("[bold]Total[/bold]", "", "", f"[bold]{total_net:+.4f}[/bold]")

    console.print()
    console.print(table)

    with open("Bader.dat", "w") as f:
        f.write("# Atom  Element  Bader(e)  NetCharge(e)\n")
        for label, sp, be, net in rows:
            f.write(f"{label:>6s}  {sp:>3s}  {be:12.4f}  {net:12.4f}\n")
    console.print("  [green]✓ Bader charges saved to Bader.dat[/green]")
    console.print()


# HIDDEN (2026-08-26): ABACUS v3.10 does not support out_hirshfeld (the keyword
# is fatal), so this task can never find data on the current stack.
# @task(1305, category="Population", name="Hirshfeld Charge",
#       description="Hirshfeld atomic charges from ABACUS LCAO output (out_hirshfeld 1)")
def task_hirshfeld_charge(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display Hirshfeld atomic charges, parsed from the ABACUS running log."""
    from abacuscopilot.postprocessing.population_tasks import (
        _display_population_table,
        _find_running_log,
        parse_hirshfeld_from_log,
    )
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== Hirshfeld Charge ===[/bold cyan]")
    console.print()

    log_path = None
    if args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                log_path = p
                break
    if log_path is None:
        log_path = _find_running_log()
    if log_path is None:
        console.print("[red]No ABACUS running log found.[/red]")
        console.print("[dim]Run an ABACUS LCAO SCF with out_hirshfeld 1 first.[/dim]")
        return

    data = parse_hirshfeld_from_log(log_path)
    if data is None:
        console.print("[yellow]No Hirshfeld data found in the log.[/yellow]")
        console.print("[dim]Set out_hirshfeld 1 in INPUT and re-run the SCF.[/dim]")
        return

    console.print(f"  [dim]Reading: {log_path}[/dim]")
    _display_population_table(console, data)
    console.print()


# =============================================================================
# RESP (Restrained ElectroStatic Potential) charges — for isolated molecules
# =============================================================================


def _downsample_for_esp(data: np.ndarray, cell_bohr: np.ndarray,
                        max_spacing_bohr: float = 0.4) -> np.ndarray:
    """Decimate the density grid so the diagonal spacing ≤ *max_spacing_bohr*.

    The ESP at vdW shells is a smooth 1/r average, so a ~0.2 Å grid is plenty;
    this cuts the convolution cost ~(spacing/0.2)³ for fine-grid densities.
    """
    nz, ny, nx = data.shape
    step = [max(1, int(round(max_spacing_bohr / (abs(cell_bohr[i, i]) / n))))
            for i, n in enumerate((nz, ny, nx))]
    if max(step) > 1:
        return data[:: step[0], :: step[1], :: step[2]]
    return data


def _atom_cartesian_angstrom(structure) -> np.ndarray:
    """Atomic positions in Cartesian Å, converting Direct (fractional) if needed.

    STRU files may use Direct coordinates; ESP/Bader need real-space positions.
    """
    cell_ang = structure.lattice.cell_angstrom  # rows = a, b, c (Å)
    pos = np.array([a.position for a in structure.atoms], dtype=float)
    if getattr(structure, "coordinate_type", "Direct").lower().startswith("direct"):
        pos = pos @ cell_ang
    return pos


def _esp_points_for_atoms(structure, shell_factors=(1.4, 1.6, 1.8, 2.0),
                          n_per_shell: int = 60) -> np.ndarray:
    """Generate ESP sampling points on vdW shells around each atom (Bohr).

    ``n_per_shell`` defaults to 60: fewer points (e.g. 20) make the fit
    sensitive to the (non-mirror-symmetric) golden-angle sampling, so
    symmetry-equivalent atoms (e.g. the two H of H2O) can get unequal charges.
    """
    from ase.data import atomic_numbers, vdw_radii  # vdw_radii indexed by Z
    cart = _atom_cartesian_angstrom(structure)
    pts: list[np.ndarray] = []
    ga = np.pi * (3 - np.sqrt(5))  # golden angle → ~equidistributed sphere
    for i, a in enumerate(structure.atoms):
        pos = cart[i] * ANGSTROM_TO_BOHR
        z = atomic_numbers.get(a.species, 1)
        r_vdw = float(vdw_radii[z]) * ANGSTROM_TO_BOHR if z < len(vdw_radii) else 1.7
        for f in shell_factors:
            r = r_vdw * f
            for m in range(n_per_shell):
                y = 1 - 2 * (m + 0.5) / n_per_shell
                rp = np.sqrt(max(0.0, 1 - y * y))
                th = ga * m
                pts.append(pos + r * np.array([rp * np.cos(th), rp * np.sin(th), y]))
    return np.array(pts)


def _compute_esp(points: np.ndarray, structure, rho: np.ndarray,
                 cell_bohr: np.ndarray, origin: np.ndarray) -> np.ndarray:
    """Electrostatic potential (a.u.) at *points* from the density.

    ESP(r) = Σ_A Z_A/|r−R_A| − ∫ρ(r')/|r−r'| dr'.  Density ρ in e/Bohr³.

    The nuclear charge used is the VALENCE charge Z_val (from the UPF): the
    density from a pseudopotential calculation is valence-only, and treating
    the atom as Z_val + valence density reproduces the molecular ESP outside
    the (localized) core region.
    """
    nz, ny, nx = rho.shape
    k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    grid = (origin[None, None, None, :]
            + (i[..., None] / nx) * cell_bohr[0]
            + (j[..., None] / ny) * cell_bohr[1]
            + (k[..., None] / nz) * cell_bohr[2]).reshape(-1, 3)
    dv = abs(np.linalg.det(cell_bohr)) / (nx * ny * nz)

    # Mask to the molecule's occupied region: the 1/r sum over a molecule in a
    # big vacuum box would otherwise scale with the BOX volume.  Dropping the
    # negligible low-density tail (ρ < 1e-6·max) makes it scale with the
    # molecular volume — ~100× faster for water, crucial for bigger molecules.
    rho_flat = rho.ravel()
    keep = rho_flat > rho_flat.max() * 1e-6
    grid = grid[keep].astype(np.float32)
    rho_dv = rho_flat[keep].astype(np.float32) * np.float32(dv)

    atom_pos = _atom_cartesian_angstrom(structure) * ANGSTROM_TO_BOHR
    atom_z = np.array([_zval_for(a.species, structure) for a in structure.atoms])

    # nuclear term (vectorized over points × atoms)
    d_nuc = np.linalg.norm(points[:, None, :] - atom_pos[None, :, :], axis=-1)
    esp = np.sum(atom_z[None, :] / d_nuc, axis=1)

    # electronic term (batched over points to bound memory; batch size shrinks
    # for big molecules so the per-batch distance array stays small)
    batch = max(1, min(20, 4_000_000 // max(1, len(grid))))
    for lo in range(0, len(points), batch):
        hi = min(lo + batch, len(points))
        d = np.linalg.norm(points[lo:hi, None, :] - grid[None, :, :], axis=-1)
        esp[lo:hi] -= np.sum(rho_dv[None, :] / d, axis=1)
    return esp


def _detect_equivalent_atoms(structure, tol: float = 0.05) -> list[list[int]]:
    """Group symmetry-equivalent atoms (same element + identical distance signature).

    Used for the RESP stage-2 constraint so equivalent atoms get equal charges
    (e.g. the two H of H2O).  Compares, for each same-element pair, the sorted
    list of distances to all other atoms.
    """
    from scipy.spatial.distance import cdist
    cart = _atom_cartesian_angstrom(structure)  # Å
    n = len(structure.atoms)
    dist = cdist(cart, cart)
    groups: list[list[int]] = []
    used: set[int] = set()
    for i in range(n):
        if i in used:
            continue
        grp = [i]
        for j in range(i + 1, n):
            if j in used:
                continue
            if structure.atoms[i].species != structure.atoms[j].species:
                continue
            si = np.sort(np.delete(dist[i], [i, j]))
            sj = np.sort(np.delete(dist[j], [i, j]))
            if np.allclose(si, sj, atol=tol):
                grp.append(j)
        if len(grp) >= 2:
            groups.append(grp)
            used.update(grp)
    return groups


def _resp_fit(points: np.ndarray, esp: np.ndarray, structure,
              k: float = 0.01, a: float = 0.1,
              total_charge: float = 0.0,
              equiv_groups: list[list[int]] | None = None,
              maxiter: int = 60, tol: float = 1e-9) -> np.ndarray:
    """Restrained ESP fit (RESP): minimize χ² + Σ_A k[(q_A²+a²)^{1/2}−a].

    Enforces total-charge conservation (Σ q = *total_charge*, default 0 for a
    neutral molecule) and, when *equiv_groups* is given, equal charges within
    each group (RESP stage-2 symmetry constraints).  Constraints are enforced
    via Lagrange multipliers.  The hyperbolic restraint is applied by
    iteratively reweighting the diagonal (weight k/sqrt(q²+a²)).
    Returns per-atom charges.
    """
    atom_pos = _atom_cartesian_angstrom(structure) * ANGSTROM_TO_BOHR
    A = 1.0 / np.linalg.norm(points[:, None, :] - atom_pos[None, :, :], axis=-1)
    n = A.shape[1]

    # Build constraint rows C q = c : total charge, then q_i − q_k = 0 per group.
    c_rows: list[np.ndarray] = [np.ones(n)]
    c_vals: list[float] = [total_charge]
    if equiv_groups:
        for grp in equiv_groups:
            g0 = grp[0]
            for gi in grp[1:]:
                row = np.zeros(n)
                row[g0] = 1.0
                row[gi] = -1.0
                c_rows.append(row)
                c_vals.append(0.0)
    ncon = len(c_rows)
    C = np.array(c_rows)
    cv = np.array(c_vals)

    M = np.zeros((n + ncon, n + ncon))
    M[:n, :n] = A.T @ A
    M[:n, n:] = C.T
    M[n:, :n] = C
    b = np.zeros(n + ncon)
    b[:n] = A.T @ esp
    b[n:] = cv
    q = np.linalg.solve(M, b)[:n]
    for _ in range(maxiter):
        M2 = M.copy()
        M2[:n, :n] = A.T @ A + np.diag(k / np.sqrt(q ** 2 + a ** 2))
        q_new = np.linalg.solve(M2, b)[:n]
        if np.max(np.abs(q_new - q)) < tol:
            q = q_new
            break
        q = q_new
    return q


# HIDDEN (2026-08-26): RESP fit is unstable for large molecules (ill-conditioned
# LSQ — 36-atom coronene gives ±8 charges; restraint far too weak vs A^T A).
# Under development; re-enable by un-commenting the @task decorator.
# @task(1306, category="Population", name="RESP Charge",
#       description="Restrained ElectroStatic Potential (RESP) charges — for isolated molecules")
def task_resp_charge(args: list[str] | None = None, interactive: bool = True) -> None:
    """Fit RESP atomic charges to the molecular electrostatic potential.

    For ISOLATED MOLECULES in a large box (e.g. ABACUS LCAO/PW with a big
    vacuum box), not periodic solids.  Reads the charge density, computes the
    ESP on vdW shells, and fits charges with the RESP hyperbolic restraint.
    """
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== RESP Charge (ESP fitting) ===[/bold cyan]")
    console.print()
    console.print("  [yellow]Note: RESP fits atomic charges to a molecular electrostatic potential —")
    console.print("  [yellow]for ISOLATED MOLECULES (large vacuum box), not periodic solids.[/yellow]")
    console.print()

    data, cell, origin, desc = _find_charge_density()
    if data is None:
        console.print("[red]Could not read a charge density file.[/red]")
        console.print("[dim]Looked for SPIN1_CHG.cube / CHG / ABACUS-CHARGE-DENSITY.restart.[/dim]")
        return
    if not Path("STRU").exists():
        console.print("[red]No STRU file — cannot map RESP charges to atoms.[/red]")
        return
    structure = read_stru("STRU")
    if not structure.atoms:
        console.print("[red]STRU contains no atoms.[/red]")
        return

    # Density unit sanity check (ABACUS cube may be e/Bohr³ or e/Å³).
    nz, ny, nx = data.shape
    dv = abs(np.linalg.det(cell)) / (nx * ny * nz)
    z_total = sum(_zval_for(a.species, structure) for a in structure.atoms)
    total = float(np.sum(data) * dv)
    if abs(total - z_total) > 0.1 * z_total:
        data_conv = data * ANGSTROM_TO_BOHR ** 3  # e/Å³ → e/Bohr³
        total_conv = float(np.sum(data_conv) * dv)
        if abs(total_conv - z_total) < abs(total - z_total):
            data = data_conv
            total = total_conv
    console.print(f"  [dim]Density: {desc} ({nx}×{ny}×{nz})[/dim]")
    console.print(f"  [dim]Valence electrons: {z_total:.0f}, density integral: {total:.2f}[/dim]")

    console.print("  Computing ESP on vdW shells...")
    # Scale sampling density by molecule size: small molecules keep 60 pts/shell,
    # big ones cap the total (~8000 ESP points) so the fit stays fast.
    n_shell = max(10, min(60, 8000 // max(1, len(structure.atoms) * 4)))
    pts = _esp_points_for_atoms(structure, n_per_shell=n_shell)
    # Downsample the density for the ESP convolution (a ~0.2 Å grid suffices).
    data_esp = _downsample_for_esp(data, cell)
    esp = _compute_esp(pts, structure, data_esp, cell, origin)

    # Symmetry-equivalent atoms → RESP stage-2: force equal charges so
    # equivalent atoms (e.g. the two H of H2O) don't split numerically.
    equiv_groups = _detect_equivalent_atoms(structure)
    labels = _atom_labels(structure)
    if equiv_groups:
        for g in equiv_groups:
            console.print(f"  [dim]Equivalent atoms (constrained equal): "
                          f"{', '.join(labels[i] for i in g)}[/dim]")

    console.print("  Fitting RESP charges (two-stage: restraint + equivalent-atom constraints)...")
    # Stage 1: all atoms free, standard restraint weight.
    _resp_fit(pts, esp, structure, k=0.01, total_charge=0.0)
    # Stage 2: smaller restraint + equal-charge constraints on equivalent atoms.
    charges = _resp_fit(pts, esp, structure, k=0.001, total_charge=0.0,
                        equiv_groups=equiv_groups)
    from rich.table import Table
    table = Table(title=f"RESP Charges ({len(structure.atoms)} atoms)")
    table.add_column("Atom")
    table.add_column("Element")
    table.add_column("RESP Charge (e)")
    rows = []
    for i, a in enumerate(structure.atoms):
        rows.append((labels[i], a.species, charges[i]))
        table.add_row(labels[i], a.species, f"{charges[i]:+.6f}")
    table.add_section()
    table.add_row("[bold]Sum[/bold]", "", f"[bold]{sum(charges):+.6f}[/bold]")
    console.print()
    console.print(table)

    with open("RESP.dat", "w") as f:
        f.write("# Atom  Element  RESPCharge(e)\n")
        for label, sp, q in rows:
            f.write(f"{label:>6s}  {sp:>3s}  {q:12.6f}\n")
    console.print("  [green]✓ RESP charges saved to RESP.dat[/green]")
    console.print()
