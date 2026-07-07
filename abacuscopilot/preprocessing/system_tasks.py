"""System configuration tasks for abacuscopilot.

Task IDs 001-099

System setup, configuration editing, pseudopotential path management,
and environment checking.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from abacuscopilot.config import load_config, save_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.tasks import task


def _prompt(console, question: str, default: Any = None) -> str:
    if default is not None:
        console.print(f"  {question} [", end="")
        console.print(str(default), style="dim", end="")
        result = console.input("]: ")
        return result.strip() if result.strip() else str(default)
    return console.input(f"  {question}: ").strip()


# =============================================================================
# Public API: prepare calculation files
# =============================================================================


def _find_file_for_element(library_dir: str, element: str, suffix: str) -> str | None:
    """Find a file in library_dir whose name starts with '{Element}_' and ends with suffix.

    Returns the filename (not full path) if found, or None.
    """
    lib = Path(library_dir)
    if not lib.is_dir():
        return None
    prefix = f"{element}_"
    for f in sorted(lib.iterdir()):
        if f.is_file() and f.name.startswith(prefix) and f.name.endswith(suffix):
            return f.name
    return None


def read_species_from_stru(stru_path: str | Path = "STRU") -> list[str]:
    """Extract element species from a STRU file.

    Reads the ATOMIC_SPECIES section to get unique element labels.
    Falls back to scanning the ATOMIC_POSITIONS block if needed.

    Returns a list of element symbols (e.g., ['Si', 'O']).
    """
    stru_path = Path(stru_path)
    if not stru_path.exists():
        return []

    with open(stru_path) as f:
        lines = f.readlines()

    species = []
    in_species = False
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("!"):
            continue
        upper = stripped.upper()
        if upper.startswith("ATOMIC_SPECIES"):
            in_species = True
            continue
        if in_species:
            # Next section keyword ends ATOMIC_SPECIES
            if any(upper.startswith(kw) for kw in (
                "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS",
                "LATTICE_PARAMETERS", "ATOMIC_POSITIONS",
            )):
                break
            # Each line: Element mass pseudo_file [pseudo_type]
            parts = stripped.split()
            if parts and parts[0] not in species:
                species.append(parts[0])

    return species


def prepare_calculation_files(
    species: list[str],
    basis_type: str,
    pseudo_library: str,
    orbital_library: str = "",
    target_dir: str | Path = ".",
    dry_run: bool = False,
) -> dict:
    """Copy pseudopotential (and orbital if LCAO) files to the target directory.

    Args:
        species: List of element symbols (e.g., ['Si', 'O']).
        basis_type: 'pw', 'lcao', or 'lcao_in_pw'.
        pseudo_library: Directory containing .upf pseudopotential files.
        orbital_library: Directory containing .orb numerical orbital files.
        target_dir: Where to copy files (default: current directory).
        dry_run: If True, only report what would be done without copying.

    Returns:
        Dict with keys 'pseudo_files', 'orbital_files', 'errors'.
    """
    target = Path(target_dir)
    result = {"pseudo_files": [], "orbital_files": [], "errors": []}

    if not pseudo_library or not Path(pseudo_library).is_dir():
        if not dry_run:
            result["errors"].append(f"Pseudopotential library not found: {pseudo_library}")
        return result

    is_lcao = basis_type.startswith("lcao")

    for elem in species:
        # --- Pseudopotential ---
        pp_file = _find_file_for_element(pseudo_library, elem, ".upf")
        if pp_file:
            src = Path(pseudo_library) / pp_file
            dst = target / pp_file
            if not dry_run:
                if not dst.exists() or src.stat().st_mtime > dst.stat().st_mtime:
                    shutil.copy2(src, dst)
            result["pseudo_files"].append(pp_file)
        else:
            msg = f"No pseudopotential found for {elem}"
            result["errors"].append(msg)

        # --- Orbital (LCAO only) ---
        if is_lcao and orbital_library:
            orb_file = _find_file_for_element(orbital_library, elem, ".orb")
            if orb_file:
                src = Path(orbital_library) / orb_file
                dst = target / orb_file
                if not dry_run:
                    if not dst.exists() or src.stat().st_mtime > dst.stat().st_mtime:
                        shutil.copy2(src, dst)
                result["orbital_files"].append(orb_file)
            else:
                msg = f"No orbital found for {elem}"
                result["errors"].append(msg)

    return result


# =============================================================================
# Task 010: Prepare Calculation Files
# =============================================================================

@task(9907, category="System", name="Prepare Files",
      description="Auto-copy pseudopotential & orbital files from library to current directory")
def task_prepare_files(args: list[str] | None = None, interactive: bool = True) -> None:
    """Copy pseudopotential (and orbital if LCAO) files for the current structure.

    Reads species from STRU file, then copies matching pseudopotential/orbital
    files from the configured library directories to the current directory.
    """
    console = _get_console()
    config = load_config()

    console.print()
    console.print("[bold cyan]=== Prepare Calculation Files ===[/bold cyan]")
    console.print()

    # Determine species
    species = read_species_from_stru("STRU")

    if not species:
        # Try reading from a CIF if available
        cif_files = list(Path(".").glob("*.cif"))
        if cif_files and interactive:
            console.print("[yellow]No STRU file found.[/yellow]")
            try:
                from ase.io import read as ase_read
                atoms = ase_read(str(cif_files[0]))
                species = list(set(atoms.get_chemical_symbols()))
                console.print(f"  Read {len(species)} species from {cif_files[0].name}")
            except ImportError:
                pass
            except Exception as e:
                console.print(f"[red]Failed to read CIF: {e}[/red]")
                return

    if not species:
        console.print("[red]No STRU or CIF file found. Cannot determine species.[/red]")
        console.print("[dim]Run this task after generating STRU (task 201).[/dim]")
        return

    console.print(f"[bold]Elements:[/bold] {', '.join(species)}")
    console.print()

    # Determine basis type from existing INPUT file
    basis_type = "lcao"
    input_path = Path("INPUT")
    if input_path.exists():
        from abacuscopilot.io.input_file import read_input
        try:
            params = read_input(input_path)
            basis_type = params.basis_type
            console.print(f"  Basis type (from INPUT): [green]{basis_type}[/green]")
        except Exception:
            console.print(f"  Basis type: [dim]{basis_type}[/dim] (default)")
    elif interactive:
        from .input_tasks import _prompt_choice as _pc
        basis_type = _pc(console, "Basis type", ["lcao", "pw"], "lcao")
    console.print()

    # Get library paths
    libraries = config.get("libraries", {})
    pseudo_lib = libraries.get("pseudo_library", "")
    orbital_lib = libraries.get("orbital_library", "")

    if not pseudo_lib:
        console.print("[yellow]Pseudopotential library path not configured.[/yellow]")
        console.print("[dim]Run System Setup (task 1501) to configure paths, or set them in ~/.abacuscopilot/config.yaml[/dim]")
        return

    # Confirm
    console.print(f"  Pseudopotential library: [dim]{pseudo_lib}[/dim]")
    if basis_type.startswith("lcao"):
        console.print(f"  Orbital library:         [dim]{orbital_lib}[/dim]")
    console.print()

    # Dry run first
    result = prepare_calculation_files(
        species, basis_type, pseudo_lib, orbital_lib, ".", dry_run=True
    )

    if result["pseudo_files"]:
        console.print("[bold]Pseudopotential files to copy:[/bold]")
        for f in result["pseudo_files"]:
            console.print(f"  [green]✓[/green] {f}")
    if result["orbital_files"]:
        console.print("[bold]Orbital files to copy:[/bold]")
        for f in result["orbital_files"]:
            console.print(f"  [green]✓[/green] {f}")
    if result["errors"]:
        console.print("[bold red]Missing:[/bold red]")
        for e in result["errors"]:
            console.print(f"  [red]✗[/red] {e}")

    if not result["pseudo_files"] and not result["orbital_files"]:
        console.print("[red]No matching files found![/red]")
        console.print("[dim]Check that pseudo_library and orbital_library paths are correct.[/dim]")
        return

    # Execute
    console.print()
    result = prepare_calculation_files(
        species, basis_type, pseudo_lib, orbital_lib, ".", dry_run=False
    )

    copied = len(result["pseudo_files"]) + len(result["orbital_files"])
    console.print(f"[bold green]✓ {copied} file(s) copied to current directory.[/bold green]")
    console.print()


# =============================================================================
# Task 001: System Setup Wizard
# =============================================================================

@task(9901, category="System", name="System Setup",
      description="Configure abacuscopilot: pseudopotential paths, orbital paths, and defaults")
def task_system_setup(args: list[str] | None = None, interactive: bool = True) -> None:
    """Interactive setup wizard for abacuscopilot configuration.

    Configures:
    - Pseudopotential directory paths (by functional: PBE, LDA, etc.)
    - Numerical orbital directory paths (LCAO calculations)
    - ABACUS binary path
    - Default calculation parameters
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]====================================[/bold cyan]")
    console.print("[bold cyan]  AbacusKit System Setup Wizard[/bold cyan]")
    console.print("[bold cyan]====================================[/bold cyan]")
    console.print()
    console.print("[dim]This wizard will help you configure abacuscopilot.[/dim]")
    console.print("[dim]Configuration is saved to ~/.abacuscopilot/config.yaml[/dim]")
    console.print()

    config = load_config()

    # Global pseudo dir
    current = config["defaults"].get("pseudo_dir", "./")
    path = _prompt(console, "Default pseudopotential directory", current)
    config["defaults"]["pseudo_dir"] = path if path else "./"

    console.print()

    # === Orbital paths (LCAO) ===
    console.print("[bold yellow]--- Numerical Orbital Paths (LCAO) ---[/bold yellow]")
    console.print("[dim]Numerical atomic orbital files for LCAO basis[/dim]")
    console.print()

    current = config["defaults"].get("orbital_dir", "./")
    path = _prompt(console, "Default orbital directory", current)
    config["defaults"]["orbital_dir"] = path if path else "./"

    console.print()

    # === ABACUS binary ===
    console.print("[bold yellow]--- ABACUS Binary ---[/bold yellow]")
    console.print()

    current = config["paths"].get("abacus_binary", "abacus")
    path = _prompt(console, "ABACUS binary name or path", current)
    config["paths"]["abacus_binary"] = path if path else "abacus"

    current = config["paths"].get("mpirun", "mpirun")
    path = _prompt(console, "MPI launcher", current)
    config["paths"]["mpirun"] = path if path else "mpirun"

    console.print()

    # === Default calculation parameters ===
    console.print("[bold yellow]--- Default Calculation Parameters ---[/bold yellow]")
    console.print()

    current = config["defaults"].get("kspacing", 0.04)
    val = _prompt(console, "Default k-spacing (2π/Å)", str(current))
    config["defaults"]["kspacing"] = float(val) if val else 0.04

    current = config["defaults"].get("ecutwfc", 100.0)
    val = _prompt(console, "Default ecutwfc (Ry)", str(current))
    config["defaults"]["ecutwfc"] = float(val) if val else 100.0

    current = config["defaults"].get("scf_thr", 1e-7)
    val = _prompt(console, "Default SCF convergence threshold (Ry)", str(current))
    config["defaults"]["scf_thr"] = float(val) if val else 1e-7

    current = config["defaults"].get("force_thr", 0.001)
    val = _prompt(console, "Default force convergence threshold (eV/Å)", str(current))
    config["defaults"]["force_thr"] = float(val) if val else 0.001

    current = config["defaults"].get("basis_type", "pw")
    val = _prompt(console, "Default basis type (pw / lcao)", current)
    config["defaults"]["basis_type"] = val if val else "pw"

    current = config["defaults"].get("dft_functional", "pbe")
    val = _prompt(console, "Default DFT functional", current)
    config["defaults"]["dft_functional"] = val if val else "pbe"

    console.print()

    # === Save ===
    save_config(config)

    console.print("[bold green]✓ Configuration saved to ~/.abacuscopilot/config.yaml[/bold green]")
    console.print()
    console.print("[dim]You can edit this file manually at any time.[/dim]")
    console.print("[dim]Run this wizard again with: abacuscopilot -task 1501[/dim]")
    console.print()


# =============================================================================
# Task 002: Show current configuration
# =============================================================================

@task(9902, category="System", name="Show Config",
      description="Display current abacuscopilot configuration")
def task_show_config(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display the current configuration."""
    console = _get_console()
    config = load_config()

    console.print()
    console.print("[bold cyan]=== Current Configuration ===[/bold cyan]")
    console.print()

    # Pseudo-lib directory already shown above

    console.print()

    # Defaults
    console.print("[bold]Calculation Defaults:[/bold]")
    defaults = config.get("defaults", {})
    for key in ("pseudo_dir", "orbital_dir", "basis_type", "ecutwfc",
                "kspacing", "dft_functional", "scf_thr", "force_thr",
                "calculation"):
        val = defaults.get(key, "—")
        console.print(f"  {key}: {val}")

    console.print()

    # Paths
    console.print("[bold]Binary Paths:[/bold]")
    paths = config.get("paths", {})
    for key, val in paths.items():
        console.print(f"  {key}: {val}")

    console.print()

    # Plotting
    console.print("[bold]Plotting:[/bold]")
    plotting = config.get("plotting", {})
    for key in ("style", "dpi", "figure_format", "font_size"):
        val = plotting.get(key, "—")
        console.print(f"  {key}: {val}")

    console.print()
    console.print("[dim]Config file: ~/.abacuscopilot/config.yaml[/dim]")
    console.print()


# =============================================================================
# Task 003: Check environment
# =============================================================================

@task(9903, category="System", name="Check Environment",
      description="Check if required tools and files are available")
def task_check_env(args: list[str] | None = None, interactive: bool = True) -> None:
    """Check the computing environment."""
    import shutil

    console = _get_console()
    config = load_config()

    console.print()
    console.print("[bold cyan]=== Environment Check ===[/bold cyan]")
    console.print()

    checks = []

    # Python version
    import sys
    py_ver = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    checks.append(("Python", py_ver, True))

    # ABACUS binary
    abacus_bin = config["paths"].get("abacus_binary", "abacus")
    abacus_path = shutil.which(abacus_bin)
    checks.append(("ABACUS binary", abacus_path or f"{abacus_bin} (not found in PATH)",
                   abacus_path is not None))

    # MPI
    mpirun = config["paths"].get("mpirun", "mpirun")
    mpi_path = shutil.which(mpirun)
    checks.append(("MPI launcher", mpi_path or f"{mpirun} (not found)",
                   mpi_path is not None))

    # Pseudopotential paths
    pseudo_dir = config["defaults"].get("pseudo_dir", "./")
    pseudo_ok = Path(str(pseudo_dir)).exists()
    checks.append(("Pseudopotential files", "found" if pseudo_ok else "not found",
                   pseudo_ok))

    # ASE (optional)
    try:
        import ase
        checks.append(("ASE", f"v{ase.__version__}", True))
    except ImportError:
        checks.append(("ASE", "not installed (optional)", False))

    # spglib (optional)
    try:
        import spglib
        checks.append(("spglib", "installed", True))
    except ImportError:
        checks.append(("spglib", "not installed (optional)", False))

    # NumPy
    import numpy
    checks.append(("NumPy", numpy.__version__, True))

    # SciPy
    import scipy
    checks.append(("SciPy", scipy.__version__, True))

    # Matplotlib
    import matplotlib
    checks.append(("Matplotlib", matplotlib.__version__, True))

    # Display
    for name, value, ok in checks:
        icon = "[green]✓[/green]" if ok else "[yellow]![/yellow]"
        console.print(f"  {icon} {name}: {value}")

    console.print()


# =============================================================================
# Task 004: Clean working directory
# =============================================================================

_CLEAN_KEEP_PATTERNS = [
    "INPUT", "KPT", "STRU",
    "*.upf", "*.UPF",
    "*.orb", "*.ORB",
    "sub*", "Sub*", "SUB*",
]


@task(9904, category="System", name="Clean Directory",
      description="Remove all files except INPUT, KPT, STRU, *.{upf,orb}, sub* — directories untouched")
def task_clean_directory(args: list[str] | None = None, interactive: bool = True) -> None:
    """Remove temporary/result files, keeping only essential inputs.

    Preserves: INPUT, KPT, STRU, *.upf, *.orb, sub* (submission scripts).
    Directories (including OUT.*/) are never touched.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Clean Directory ===[/bold cyan]")
    console.print()

    cwd = Path(".")
    keep: set[Path] = set()

    # Collect files to keep
    for pattern in _CLEAN_KEEP_PATTERNS:
        for p in cwd.glob(pattern):
            if p.is_file():
                keep.add(p)

    # Find all files in current directory (not recursive, no dirs)
    all_files = sorted(p for p in cwd.iterdir() if p.is_file())

    to_delete = [p for p in all_files if p not in keep]

    if not to_delete:
        console.print("[green]No stray files to clean.[/green]")
    else:
        console.print(f"  Keeping ({len(keep)} files):")
        for p in sorted(keep, key=lambda x: x.name):
            console.print(f"    [green]✓[/green] {p.name}")
        console.print()
        console.print(f"  To delete ({len(to_delete)} files):")
        for p in to_delete:
            console.print(f"    [red]✗[/red] {p.name}")
    for p in sorted(keep, key=lambda x: x.name):
        console.print(f"    [green]✓[/green] {p.name}")
    console.print()
    console.print(f"  To delete ({len(to_delete)} files):")
    for p in to_delete:
        console.print(f"    [red]✗[/red] {p.name}")

    if to_delete:
        if interactive:
            console.print()
            confirm = _prompt(console, "Confirm deletion? (yes/no)", "no")
            if confirm.lower() not in ("yes", "y"):
                console.print("[yellow]File cleanup cancelled.[/yellow]")
            else:
                deleted = 0
                for p in to_delete:
                    try:
                        p.unlink()
                        deleted += 1
                    except OSError as e:
                        console.print(f"  [red]Failed to delete {p.name}: {e}[/red]")
                console.print(f"[green]✓ {deleted} file(s) deleted, {len(keep)} kept.[/green]")

    console.print()

    # --- OUT.* directories ---
    out_dirs = sorted(Path(".").glob("OUT.*"))
    out_dirs = [d for d in out_dirs if d.is_dir()]
    if out_dirs:
        if interactive:
            rm_out = _prompt_choice(console, "Delete OUT.* directories?",
                                    ["Yes, delete them", "No, keep them"],
                                    "No, keep them")
        else:
            rm_out = "No"
        if "Yes" in rm_out:
            import shutil
            for d in out_dirs:
                try:
                    shutil.rmtree(d)
                    console.print(f"  [dim]Deleted {d.name}/[/dim]")
                except OSError as e:
                    console.print(f"  [red]Failed to delete {d.name}/: {e}[/red]")
        else:
            console.print("  [dim]OUT.* directories kept.[/dim]")

    console.print()


# =============================================================================
# Task 005: MD progress monitor
# =============================================================================


def _parse_md_progress(log_path: Path) -> list[dict]:
    """Extract MD step summaries from a running_md.log file."""
    import re as _re

    content = log_path.read_text(errors="ignore")
    results = []
    step_pat = _re.compile(r"STEP OF MOLECULAR DYNAMICS\s*:\s*(\d+)", _re.IGNORECASE)
    num_pat = _re.compile(r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)")

    lines = content.split("\n")
    # Pattern for the temperature table header
    temp_header = _re.compile(r"Energy\s*\(Ry\)\s+Potential\s*\(Ry\)\s+Kinetic\s*\(Ry\)\s+Temperature", _re.IGNORECASE)

    i = 0
    while i < len(lines):
        m = step_pat.search(lines[i])
        if m:
            step = int(m.group(1))
            # Search forward for the temperature header, then take the next numeric line
            for j in range(i + 1, len(lines)):
                if step_pat.search(lines[j]):
                    break  # next step reached
                if temp_header.search(lines[j]):
                    # Next meaningful line should have the 4 values
                    for k in range(j + 1, min(j + 5, len(lines))):
                        nums = num_pat.findall(lines[k])
                        if len(nums) >= 4:
                            try:
                                results.append({
                                    "step": step,
                                    "energy_ry": float(nums[0]),
                                    "potential_ry": float(nums[1]),
                                    "kinetic_ry": float(nums[2]),
                                    "temperature_k": float(nums[3]),
                                })
                            except ValueError:
                                pass
                            break
                    break
        i += 1
    return results


def _find_md_log() -> Path | None:
    """Find the MD log file."""
    for p in Path().glob("running_md.log"):
        return p
    for p in Path().glob("OUT.*/running_md.log"):
        return p
    return None


@task(9905, category="System", name="MD Monitor",
      description="Live monitor of MD simulation: step, energy, temperature")
def task_md_monitor(args: list[str] | None = None, interactive: bool = True) -> None:
    """Continuously display MD simulation progress from running_md.log.

    Press q / Esc / Ctrl+C to exit.
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== MD Progress Monitor ===[/bold cyan]")

    log_path = _find_md_log()
    if log_path is None:
        console.print("[red]No running_md.log found.[/red]")
        return


    data = _parse_md_progress(log_path)
    total_steps = len(data)

    # Print all existing steps
    for d in data:
        console.print(
            f"  [bold]{d['step']:>6d}[/bold]  "
            f"E={d['energy_ry']:12.6f}  "
            f"V={d['potential_ry']:12.6f}  "
            f"K={d['kinetic_ry']:10.6f}  "
            f"T=[green]{d['temperature_k']:8.2f}[/green] K"
        )

    try:
        while True:
            data = _parse_md_progress(log_path)
            new_count = len(data) - total_steps

            # Print new steps
            if new_count > 0:
                for d in data[-new_count:]:
                    console.print(
                        f"  [bold]{d['step']:>6d}[/bold]  "
                        f"E={d['energy_ry']:12.6f}  "
                        f"V={d['potential_ry']:12.6f}  "
                        f"K={d['kinetic_ry']:10.6f}  "
                        f"T=[green]{d['temperature_k']:8.2f}[/green] K"
                    )
                total_steps = len(data)

            # Status bar — only on step change or first run, avoid flooding
            if data:
                last = data[-1]
                current_status = (
                    f"[dim]Step {last['step']}  |  "
                    f"T = {last['temperature_k']:.1f} K  |  "
                    f"E = {last['energy_ry']:.4f} Ry  |  "
                    f"q / Ctrl+C to exit  |  {log_path}[/dim]"
                )
                if new_count > 0:
                    console.print(current_status)
            else:
                if total_steps == 0:
                    console.print(f"[dim]Waiting for data...  |  q / Ctrl+C to exit  |  {log_path}[/dim]")

            # Check for exit key
            try:
                import select
                import sys as _sys

                r, _, _ = select.select([_sys.stdin], [], [], 1.5)
                if r:
                    c = _sys.stdin.read(1)
                    if c in ("q", "Q", "\x03"):
                        break
            except (OSError, ValueError):
                pass

    except KeyboardInterrupt:
        pass

    console.print()
    console.print("[yellow]Monitor stopped.[/yellow]")
    console.print()


# =============================================================================
# Task 006: Set job submission script path
# =============================================================================


@task(9906, category="System", name="Set Submit Script",
      description="Configure path to a SLURM/PBS submission script for batch test tasks")
def task_set_sub_script(args: list[str] | None = None, interactive: bool = True) -> None:
    """Set the path to a user's job submission script.

    This script will be copied into convergence test directories by tasks
    107/108.  If not configured, those tasks skip script generation and
    print a reminder.
    """
    console = _get_console()
    config = load_config()

    console.print()
    console.print("[bold cyan]=== Set Submission Script ===[/bold cyan]")
    console.print()

    if "paths" not in config:
        config["paths"] = {}

    current = config["paths"].get("sub_script", "")
    if current:
        expanded = str(Path(current).expanduser())
        console.print("  Current script: ", end="")
        console.print(expanded, style="green")
        if Path(expanded).exists():
            console.print("  [green]✓ File exists[/green]")
        else:
            console.print("  [yellow]! File not found[/yellow]")
        console.print()

    if interactive:
        new_path = _prompt(console, "Path to submission script (blank=keep, \"clear\"=remove)", current)
        if new_path and new_path.lower().strip() == "clear":
            config["paths"]["sub_script"] = ""
            save_config(config)
            console.print("  [dim]Script path cleared.[/dim]")
        elif new_path and new_path != current:
            expanded = str(Path(new_path).expanduser())
            config["paths"]["sub_script"] = expanded
            if Path(expanded).exists():
                console.print("  [green]✓ Script set:[/green] ", end="")
                console.print(expanded)
            else:
                console.print("  [yellow]Warning: file does not exist — check path[/yellow]")
            save_config(config)
        # else: blank → keep current, do nothing
        console.print()
