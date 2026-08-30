"""System configuration tasks for abacuscopilot.

Task IDs 001-099

System setup, configuration editing, pseudopotential path management,
and environment checking.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

from abacuscopilot.config import load_config, save_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.standards import is_lcao_basis
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


# Canonical DZP basis used as the default when a library ships several
# orbitals per element (e.g. the APNS lanthanide bundles offer 2s1p1d /
# 4s2p2d1f / 6s3p3d2f at rcut 6-10 au).  Kept in sync with the SG15 default.
_DEFAULT_ORB_BASIS = "4s2p2d1f"
_DEFAULT_ORB_RCUT = 7


def _candidate_rank(name: str, suffix: str) -> tuple:
    """Deterministic preference key for library files of one element.

    Used only when several files match an element; the smallest key wins.
    Prefers the canonical DZP orbital at 7 au (matches the SG15 convention),
    then falls back to rcut closest to 7 au, then alphabetical.
    """
    base = name.lower()
    rank = [0, 0]
    if _DEFAULT_ORB_BASIS not in base:
        rank[0] = 1
    if suffix == ".orb":
        m = re.search(r"(\d+)au", base)
        rcut = int(m.group(1)) if m else None
        rank[1] = abs(rcut - _DEFAULT_ORB_RCUT) if rcut is not None else 99
    return tuple(rank) + (name,)


def _as_dir_list(value: Any) -> list[str]:
    """Normalize a config library value (str or list) to a list of strings."""
    if not value:
        return []
    return [value] if isinstance(value, str) else list(value)


def _find_file_in_libraries(library_dirs, element: str, suffix: str) -> Path | None:
    """Locate the best file for *element* among the library dirs.

    *library_dirs* may be a single path string or a list of paths.  Each dir
    is searched recursively (the APNS lanthanide bundles nest files under
    element/basis subfolders).  Filenames are matched case-insensitively on a
    prefix of ``{Element}`` followed by ``_`` or ``.`` (optionally with a
    charge state such as ``Sm3+_f--core-icmod1.PD04.PBE.UPF``) and the given
    suffix (also case-insensitive).  Returns the best-matching Path or None.
    """
    pattern = re.compile(rf"^{re.escape(element)}(?:\d+\+)?[_.]", re.IGNORECASE)
    best_path: Path | None = None
    best_key: tuple | None = None
    for d in _as_dir_list(library_dirs):
        root = Path(d)
        if not root.is_dir():
            continue
        for f in sorted(root.rglob("*")):
            if not f.is_file():
                continue
            if not pattern.match(f.name):
                continue
            if not f.name.lower().endswith(suffix):
                continue
            key = _candidate_rank(f.name, suffix)
            if best_key is None or key < best_key:
                best_key = key
                best_path = f
    return best_path


def _find_file_for_element(library_dir: Any, element: str, suffix: str) -> str | None:
    """Return the best-matching filename (not path) for *element*, or None.

    Thin wrapper over :func:`_find_file_in_libraries` kept for callers that
    only need the file name.
    """
    p = _find_file_in_libraries(library_dir, element, suffix)
    return p.name if p else None


# =============================================================================
# Large-core (f-electron-pseudized) lanthanide awareness
# =============================================================================

# Filename markers of the APNS lanthanide bundle: f-electrons in the core,
# icmod1 model-core-charge correction (e.g. 'Sm3+_f--core-icmod1.PD04.PBE.UPF').
_F_CORE_MARKERS = ("f--core", "icmod1")

# Orbital energy cutoffs above this (Ry) are atypical for the 100 Ry SG15
# orbitals and trigger the ecutwfc warning (APNS lanthanides use 300 Ry).
_F_CORE_ECUT_ALERT = 100.0


def _f_core_info_from_structure(structure) -> dict:
    """Identify large-core f-electron-pseudized species from resolved files.

    A species is flagged when its resolved pseudopotential filename carries an
    f-core marker (e.g. ``Sm3+_f--core-icmod1.PD04.PBE.UPF``).  Orbital energy
    cutoffs (Ry) are parsed from the resolved orbital filenames.

    Returns:
        dict with ``f_core_species`` [(species, pp_name), ...], ``ecut_entries``
        [(species, orb_name, ecut), ...] and ``max_orb_ecut`` (float, 0 if none).
    """
    f_core_species: list[tuple[str, str]] = []
    ecut_entries: list[tuple[str, str, float]] = []
    for sp in structure.species_order:
        pp = structure.pseudo_files.get(sp, "")
        if pp and any(m in pp.lower() for m in _F_CORE_MARKERS):
            f_core_species.append((sp, pp))
        orb = structure.orbital_files.get(sp, "")
        m = re.search(r"(\d+(?:\.\d+)?)ry", orb.lower()) if orb else None
        if m:
            ecut_entries.append((sp, orb, float(m.group(1))))
    max_ecut = max((e[2] for e in ecut_entries), default=0.0)
    return {"f_core_species": f_core_species, "ecut_entries": ecut_entries,
            "max_orb_ecut": max_ecut}


def analyze_f_core(structure_or_species, libraries: dict | None = None) -> dict:
    """Like :func:`_f_core_info_from_structure`, but resolves filenames from libraries.

    Useful before the STRU is written, or for structures/species that don't yet
    carry resolved pseudo/orbital filenames.  Accepts a :class:`Structure` or a
    plain iterable of element symbols (e.g. from reading an existing STRU).

    Args:
        structure_or_species: A Structure, or an iterable of element symbols.
        libraries: Config ``libraries`` dict (auto-loaded if None).
    """
    if libraries is None:
        from abacuscopilot.config import load_config
        libraries = load_config().get("libraries", {})
    species = (structure_or_species.species_order
               if hasattr(structure_or_species, "species_order")
               else list(structure_or_species))
    pseudo_lib = libraries.get("pseudo_library", "")
    orb_lib = libraries.get("orbital_library", "")

    f_core_species: list[tuple[str, str]] = []
    ecut_entries: list[tuple[str, str, float]] = []
    for sp in species:
        pp = _find_file_for_element(pseudo_lib, sp, ".upf")
        if pp and any(m in pp.lower() for m in _F_CORE_MARKERS):
            f_core_species.append((sp, pp))
        orb = _find_file_for_element(orb_lib, sp, ".orb")
        m = re.search(r"(\d+(?:\.\d+)?)ry", orb.lower()) if orb else None
        if m:
            ecut_entries.append((sp, orb, float(m.group(1))))
    max_ecut = max((e[2] for e in ecut_entries), default=0.0)
    return {"f_core_species": f_core_species, "ecut_entries": ecut_entries,
            "max_orb_ecut": max_ecut}


def warn_f_core(console, info: dict) -> None:
    """Print a prominent notice about large-core (f-electron-pseudized) species.

    Covers what these PPs are, what they are suited for, what they are NOT
    suited for, and that they do not apply to unaries.
    """
    if not info.get("f_core_species"):
        return
    console.print()
    console.print("[bold yellow]⚠ 大核赝势提醒 (large-core / f-electron-pseudized PP)[/bold yellow]")
    for sp, pp in info["f_core_species"]:
        console.print(f"  [yellow]{sp}[/yellow]  →  {pp}")
    console.print("  [yellow]• f 电子赝化进芯、按 +3 价生成的大核赝势。[/yellow]")
    console.print("  [yellow]• 适合:结构优化 / MD / 离子输运等占据态性质;Li3MCl6 型卤化物(+3 金属)是官方目标体系。[/yellow]")
    console.print("  [yellow]• 不适合:需要 f 电子参与的计算(磁性、光谱、含 f 的非占据态)。[/yellow]")
    console.print("  [yellow]• 不适用于单质 / 金属体系(通常无法收敛)。[/yellow]")
    max_ecut = info.get("max_orb_ecut", 0.0)
    if max_ecut > _F_CORE_ECUT_ALERT:
        console.print(
            f"  [yellow]• 根据官方说明,本体系轨道能量截断最高 {max_ecut:.0f} Ry,ecutwfc 需 ≥ {max_ecut:.0f} Ry;"
            f"但此值可能过于保守,强烈推荐做截断能收敛测试(任务 107 生成 / 任务 109 分析)。[/yellow]"
        )


def adjust_ecutwfc_for_f_core(console, params, info: dict,
                              interactive: bool = True) -> bool:
    """Offer to raise *params.ecutwfc* to the f-core orbital cutoff (e.g. 300 Ry).

    The APNS lanthanide NAOs carry a 300 Ry energy cutoff; leaving the default
    100 Ry would silently truncate them.  In interactive mode, prompts the user;
    if they confirm, sets ``ecutwfc`` to the highest orbital cutoff in the
    system.  In non-interactive mode, raises it automatically (running a task
    via CLI is taken as intent).

    Returns True if *params.ecutwfc* was changed.
    """
    if not info.get("f_core_species"):
        return False
    max_ecut = info.get("max_orb_ecut", 0.0)
    current = params.ecutwfc or 0
    if max_ecut <= current:
        return False
    if not interactive:
        params.ecutwfc = int(round(max_ecut))
        return True
    species_str = ", ".join(sp for sp, _ in info["f_core_species"])
    set_label = f"Set ecutwfc = {int(round(max_ecut))} Ry"
    answer = _prompt_choice(
        console,
        f"体系含 f 电子进芯大核赝势({species_str}),"
        f"其 NAO 轨道截断最高 {max_ecut:.0f} Ry"
        f"(当前 ecutwfc = {current:.0f} Ry)。建议先跑收敛测试(107/109)"
        f"确认低截断是否够用。统一提高截断能?",
        [set_label, "Keep current"],
        set_label,
    )
    if answer.startswith("Set"):
        params.ecutwfc = int(round(max_ecut))
        return True
    return False


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


def resolve_basis_type(structure=None, interactive: bool = True,
                       default: str = "lcao") -> str:
    """Authoritative determination of basis_type for a STRU-writing task.

    Standard规范: geometry-editing tasks (supercell, slab, fix atoms, ...) have
    no intrinsic basis_type, so resolve it consistently:
      1. A nearby INPUT file (its basis_type is authoritative).
      2. Whether the source structure carries orbital files (loaded from an
         LCAO STRU) — implies lcao.
      3. Interactive prompt (lcao / pw) when running interactively.
      4. Fall back to `default`.
    """
    # 1. nearby INPUT is authoritative
    input_path = Path("INPUT")
    if input_path.exists():
        try:
            from abacuscopilot.io.input_file import read_input
            bt = read_input(input_path).basis_type
            if bt:
                return bt
        except Exception:
            pass

    # 2. source structure already has orbital files -> lcao
    if structure is not None and getattr(structure, "orbital_files", None):
        return "lcao"

    # 3. ask the user
    if interactive:
        return _prompt_choice(_get_console(), "Basis type", ["lcao", "pw"], default)

    # 4. fallback
    return default


def prepare_calculation_files(
    species: list[str],
    basis_type: str,
    pseudo_library: Any,
    orbital_library: Any = "",
    target_dir: str | Path = ".",
    dry_run: bool = False,
) -> dict:
    """Copy pseudopotential (and orbital if LCAO) files to the target directory.

    Args:
        species: List of element symbols (e.g., ['Si', 'O']).
        basis_type: 'pw', 'lcao', or 'lcao_in_pw'.
        pseudo_library: Directory, or list of directories (searched in order),
            containing .upf pseudopotential files.
        orbital_library: Directory, or list of directories, containing .orb
            numerical orbital files.
        target_dir: Where to copy files (default: current directory).
        dry_run: If True, only report what would be done without copying.

    Returns:
        Dict with keys 'pseudo_files', 'orbital_files', 'errors'.
    """
    target = Path(target_dir)
    result = {"pseudo_files": [], "orbital_files": [], "errors": []}

    pseudo_dirs = _as_dir_list(pseudo_library)
    if not pseudo_dirs or not any(Path(d).is_dir() for d in pseudo_dirs):
        if not dry_run:
            result["errors"].append(
                f"Pseudopotential library not found: {pseudo_library}")
        return result

    is_lcao = is_lcao_basis(basis_type)

    for elem in species:
        # --- Pseudopotential ---
        pp_path = _find_file_in_libraries(pseudo_dirs, elem, ".upf")
        if pp_path:
            dst = target / pp_path.name
            if not dry_run:
                if not dst.exists() or pp_path.stat().st_mtime > dst.stat().st_mtime:
                    shutil.copy2(pp_path, dst)
            result["pseudo_files"].append(pp_path.name)
        else:
            msg = f"No pseudopotential found for {elem}"
            result["errors"].append(msg)

        # --- Orbital (LCAO only) ---
        if is_lcao and _as_dir_list(orbital_library):
            orb_path = _find_file_in_libraries(orbital_library, elem, ".orb")
            if orb_path:
                dst = target / orb_path.name
                if not dry_run:
                    if not dst.exists() or orb_path.stat().st_mtime > dst.stat().st_mtime:
                        shutil.copy2(orb_path, dst)
                result["orbital_files"].append(orb_path.name)
            else:
                msg = f"No orbital found for {elem}"
                result["errors"].append(msg)

    return result


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
    console.print("[bold cyan]  AbacusCopilot System Setup Wizard[/bold cyan]")
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
    console.print("[dim]Run this wizard again with: abacuscopilot -task 9901[/dim]")
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

    # Library directories (auto-detected)
    libs = config.get("libraries", {})
    console.print("[bold]PP/Orbital Libraries:[/bold]")
    console.print("  Pseudopotential:")
    for d in _as_dir_list(libs.get("pseudo_library", "")):
        console.print(f"    [dim]{d}[/dim]")
    console.print("  Orbital:")
    for d in _as_dir_list(libs.get("orbital_library", "")):
        console.print(f"    [dim]{d}[/dim]")
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


def _parse_md_progress(log_path: Path, tail: int = 0) -> list[dict]:
    """Extract MD step summaries from a running_md.log file.

    For large logs (1+ GB / 10M+ lines) we only parse the tail of the
    file to avoid reading the entire thing into memory.  ABACUS writes
    ~2700 lines per MD step, so we read ``tail`` × 3000 lines from the
    end, which is orders of magnitude faster for live monitoring.

    Args:
        log_path: path to running_md.log.
        tail: if > 0, read at most this many lines from the end of the file.
              Use tail=0 to parse the entire file (may be slow on large logs).
    """
    import os
    import re as _re

    step_pat = _re.compile(r"STEP OF MOLECULAR DYNAMICS\s*:\s*(\d+)", _re.IGNORECASE)
    num_pat = _re.compile(r"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)")
    temp_header = _re.compile(
        r"Energy\s*\(Ry\)\s+Potential\s*\(Ry\)\s+Kinetic\s*\(Ry\)\s+Temperature",
        _re.IGNORECASE,
    )

    if tail > 0:
        # Read only the tail — fast for huge logs.
        # DP logs are ~27k lines/step; LCAO/PW ~2.7k.  Use 30k to be safe.
        chunk_lines = tail * 30000
        with open(log_path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            if size > chunk_lines * 120:  # ~120 bytes per line avg
                f.seek(max(0, size - chunk_lines * 120))
                raw = f.read().decode("utf-8", errors="ignore")
                # Skip partial first line (broken by the seek)
                nl = raw.find("\n")
                content = raw[nl + 1:] if nl >= 0 else raw
            else:
                content = log_path.read_text(errors="ignore")
    else:
        content = log_path.read_text(errors="ignore")

    lines = content.split("\n")
    results = []
    i = 0
    while i < len(lines):
        m = step_pat.search(lines[i])
        if m:
            step = int(m.group(1))
            for j in range(i + 1, len(lines)):
                if step_pat.search(lines[j]):
                    break
                if temp_header.search(lines[j]):
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


@task(704, category="SCF Analysis", name="MD Monitor",
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


    import os
    import time

    # Fast tail-scan — works even on 1+ GB logs (reads only the last ~10K lines).
    data = _parse_md_progress(log_path, tail=3)
    if not data:
        console.print("[yellow]No MD steps found in log yet.[/yellow]")
        return

    last = data[-1]
    size_gb = os.path.getsize(log_path) / 1e9
    console.print(f"  [dim]Log: {log_path.name} ({size_gb:.1f} GB)[/dim]")
    console.print(
        f"  Latest: step [bold]{last['step']:>6d}[/bold]  "
        f"E={last['energy_ry']:12.6f} Ry  "
        f"T=[green]{last['temperature_k']:8.2f}[/green] K"
    )

    known_step = last["step"]
    try:
        while True:
            new_data = _parse_md_progress(log_path, tail=3)
            for d in new_data:
                if d["step"] > known_step:
                    console.print(
                        f"  [bold]{d['step']:>6d}[/bold]  "
                        f"E={d['energy_ry']:12.6f}  "
                        f"T=[green]{d['temperature_k']:8.2f}[/green] K"
                    )
                    known_step = d["step"]

            # Check for exit key
            try:
                import select
                import sys as _sys
                r, _, _ = select.select([_sys.stdin], [], [], 2.0)
                if r:
                    c = _sys.stdin.read(1)
                    if c in ("q", "Q", "\x03"):
                        break
            except (OSError, ValueError):
                pass
            time.sleep(0.5)

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
