"""Charge density analysis tasks for ABACUS output.

Task IDs 741-759

Reads charge density files (CHG*, SPIN*) from ABACUS output,
computes 1D planar averages, provides Cube/XSF export, and
handles difference charge density analysis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.console_utils import _get_console
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


def read_charge_cube(filepath: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Read a Cube-format charge density file.

    Returns (data_3d, cell_angstrom, origin) or None if not a Cube file.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    with open(filepath) as f:
        header = f.readline().strip()
    if "CUBE" not in header.upper() and "CUBE" not in f.name.upper():
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

    data_start += 1
    # Grid dimensions
    parts = lines[data_start].split()
    nx = int(parts[0])
    data_start += 1
    # Read grid spacing vectors
    cell = np.zeros((3, 3))
    for i in range(3):
        parts = lines[data_start + i].split()
        n = int(parts[0])
        cell[i] = np.array([float(parts[1]) * n, float(parts[2]) * n, float(parts[3]) * n])

    data_start += 3
    # Skip atom positions
    data_start += natom

    # Read nx, ny, nz
    ny = int(lines[data_start].split()[0]) if data_start < len(lines) else nx
    data_start += 1
    nz = int(lines[data_start].split()[0]) if data_start < len(lines) else nx
    data_start += 1

    # Read volumetric data
    raw = []
    for line in lines[data_start:]:
        raw.extend(float(x) for x in line.split())

    expected = nx * ny * nz
    data = np.array(raw[:expected]).reshape(nz, ny, nx)  # z, y, x ordering for Cube

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
        f.write(f"  {nx}  {cell[0,0]/nx:.10f}  {cell[0,1]/nx:.10f}  {cell[0,2]/nx:.10f}\n")
        f.write(f"  {ny}  {cell[1,0]/ny:.10f}  {cell[1,1]/ny:.10f}  {cell[1,2]/ny:.10f}\n")
        f.write(f"  {nz}  {cell[2,0]/nz:.10f}  {cell[2,1]/nz:.10f}  {cell[2,2]/nz:.10f}\n")
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

    # Read all three
    datasets = []
    for label, fp in [("AB", chg_ab), ("A", chg_a), ("B", chg_b)]:
        if not fp or not Path(fp).exists():
            console.print(f"[red]File not found for ρ({label}): {fp}[/red]")
            return
        r = read_charge_cube(fp)
        if r is None:
            console.print(f"[red]Could not read ρ({label}) from {fp}[/red]")
            return
        datasets.append(r[0])  # data only

    d_ab, d_a, d_b = datasets

    # Check shapes match
    if d_ab.shape != d_a.shape or d_ab.shape != d_b.shape:
        console.print("[red]Grid dimensions must match for all three files.[/red]")
        return

    diff = d_ab - d_a - d_b

    write_cube("diff_charge_density.cube", diff, datasets[0])

    console.print()
    console.print("[green]✓ Difference charge density written: diff_charge_density.cube[/green]")
    console.print("  Open with VESTA/Jmol for visualization.")
    console.print()


def _prompt_path(console, question: str) -> str:
    """Prompt for a file path."""
    from rich.prompt import Prompt
    return Prompt.ask(f"  {question}")
