"""INPUT file generation tasks for ABACUS calculations.

Task IDs 101-199

Provides interactive INPUT file generation with sensible defaults for
common calculation types. Each (calculation, basis_type) pair has a
well-tested template — the interactive mode only asks the key questions
and applies proven defaults for the rest.

Templates are defined as dicts below and applied to InputParams before
writing. Users can override any parameter via config or CLI flags.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from abacuscopilot.config import load_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.models import InputParams
from abacuscopilot.tasks import task

# =============================================================================
# Helper utilities
# =============================================================================


def _prompt_choice(console, question: str, options: list[str], default: str = "") -> str:
    """Prompt user to choose from a list of options."""
    console.print(f"  {question}")
    for i, opt in enumerate(options, 1):
        marker = " [dim](default)[/dim]" if opt == default else ""
        console.print(f"    {i}. {opt}{marker}")

    answer = console.input("  Choice: ").strip()

    if not answer and default:
        return default

    try:
        idx = int(answer) - 1
        if 0 <= idx < len(options):
            return options[idx]
    except ValueError:
        if answer in options:
            return answer

    return default or options[0]


def _apply_template(params: InputParams, template: dict) -> InputParams:
    """Apply a template dict to an InputParams object.

    Template keys are Python field names (lowercase).
    Values can be plain (``"scf"``) or a ``(value, "inline comment")`` tuple.
    Use ``None`` to mean "use the InputParams default".

    Special meta-keys (prefixed with ``_``):
    - ``_comment_hints``: dict[str, str] — bottom-of-file ``#param  value`` hints.
    - ``_section_hints``: dict[str, list[str]] — commented lines within sections.
    - ``_inline_comments``: dict[str, str] — inline ``# comment`` after a param.
      (Can also be specified as ``(value, comment)`` tuples on regular params.)
    """
    touched: list[str] = []
    inline: dict[str, str] = {}
    for key, entry in template.items():
        if key.startswith("_"):
            if key == "_comment_hints" and isinstance(entry, dict):
                params.extras["_comment_hints"] = entry
            elif key == "_section_hints" and isinstance(entry, dict):
                params.extras["_section_hints"] = entry
            elif key == "_inline_comments" and isinstance(entry, dict):
                inline.update(entry)
            continue

        # Unpack (value, comment) tuples
        if isinstance(entry, tuple) and len(entry) == 2:
            value, comment = entry
            if value is not None:
                params.set_param(key, value)
                touched.append(key)
            if comment:
                inline[key] = comment
        elif entry is not None:
            params.set_param(key, entry)
            touched.append(key)

    params.extras["_template_keys"] = touched
    if inline:
        # Merge with any existing inline comments
        existing = params.extras.get("_inline_comments", {})
        existing.update(inline)
        params.extras["_inline_comments"] = existing
    return params


def _check_charge_density(console) -> None:
    """Check whether ABACUS-CHARGE-DENSITY.restart is present.

    NSCF/DOS calculations need the charge density from a preceding SCF run.
    If the file is missing, print a clear reminder.
    """
    if Path("ABACUS-CHARGE-DENSITY.restart").exists():
        console.print("  [green]✓[/green] ABACUS-CHARGE-DENSITY.restart found")
    elif Path("OUT.ABACUS/ABACUS-CHARGE-DENSITY.restart").exists():
        console.print("  [green]✓[/green] OUT.ABACUS/ABACUS-CHARGE-DENSITY.restart found")
    else:
        console.print("  [yellow]![/yellow] ABACUS-CHARGE-DENSITY.restart not found in current directory")
        console.print("    [dim]Copy it from your SCF directory before running NSCF:[/dim]")
        console.print("    [dim]  cp ../000000.1/OUT.ABACUS/ABACUS-CHARGE-DENSITY.restart ./[/dim]")


# =============================================================================
# INPUT templates — proven defaults for each (calculation, basis_type) pair.
# Uses base+override pattern: one change here propagates to all templates.
# =============================================================================

# --- Base templates (shared fields) ------------------------------------------

_LCAO_BASE = {
    "symmetry": 1, "kspacing": 0.14, "precision": "double",
    "ecutwfc": 100.0, "basis_type": "lcao", "ks_solver": "genelpa",
    "smearing_method": "gauss", "smearing_sigma": 0.01,
    "mixing_type": "broyden", "mixing_beta": 0.8, "scf_nmax": 100, "scf_thr": 1e-7,
}

_PW_BASE = {
    "symmetry": 1, "kspacing": 0.14, "precision": "double",
    "ecutwfc": 80.0, "pw_diag_nmax": 20, "pw_diag_ndim": 2,
    "basis_type": "pw", "ks_solver": "dav_subspace",
    "smearing_method": "gauss", "smearing_sigma": 0.01,
    "mixing_type": "broyden", "mixing_beta": 0.8, "scf_nmax": 100, "scf_thr": 1e-8,
}

_RELAX_FIELDS = {
    "relax_method": "cg", "relax_nmax": 60,
    "cal_force": 1, "force_thr_ev": 0.01, "stress_thr": 0.5, "fixed_axes": "None",
}

_NSCF_READ = {"read_file_dir": "./", "init_chg": "file"}

# --- LCAO templates ----------------------------------------------------------

TEMPLATE_LCAO_CELL_RELAX = {**_LCAO_BASE, "calculation": "cell-relax", **_RELAX_FIELDS, "cal_stress": 1}
TEMPLATE_LCAO_RELAX      = {**_LCAO_BASE, "calculation": "relax",      **_RELAX_FIELDS, "cal_stress": 0}
TEMPLATE_LCAO_SCF         = {**_LCAO_BASE, "calculation": "scf"}
TEMPLATE_LCAO_BAND        = {**_LCAO_BASE, "calculation": "nscf", "symmetry": 0, "kspacing": None, **_NSCF_READ, "out_band": True}
TEMPLATE_LCAO_DOS         = {**_LCAO_BASE, "calculation": "nscf", **_NSCF_READ, "out_dos": 1, "dos_emin_ev": -15.0, "dos_emax_ev": 15.0, "dos_edelta_ev": 0.01, "dos_sigma": 0.03}

# --- PW templates ------------------------------------------------------------

TEMPLATE_PW_CELL_RELAX = {**_PW_BASE, "calculation": "cell-relax", **_RELAX_FIELDS, "cal_stress": 1}
TEMPLATE_PW_RELAX      = {**_PW_BASE, "calculation": "relax",      **_RELAX_FIELDS, "cal_stress": 0}
TEMPLATE_PW_SCF         = {**_PW_BASE, "calculation": "scf",
    "_comment_hints": {"cal_force": "1", "cal_stress": "1"}}
TEMPLATE_PW_BAND        = {**_PW_BASE, "calculation": "nscf", "symmetry": 0, "kspacing": None, **_NSCF_READ, "out_band": True}
TEMPLATE_PW_DOS         = {**_PW_BASE, "calculation": "nscf", **_NSCF_READ, "out_dos": 1, "dos_emin_ev": -15.0, "dos_emax_ev": 15.0, "dos_edelta_ev": 0.01, "dos_sigma": 0.03}

# --- MD templates (distinct enough to keep inline) ---------------------------

TEMPLATE_LCAO_MD = {
    "calculation": "md", "symmetry": 0, "gamma_only": 1,
    "precision": ("double", "or single"),
    "ecutwfc": 100.0, "basis_type": "lcao",
    "ks_solver": ("genelpa", "cusolver for single GPU task or cusolvermp for multi GPUs task"),
    "smearing_method": "gauss", "smearing_sigma": 0.01,
    "mixing_type": "broyden", "mixing_beta": 0.8, "scf_nmax": 100, "scf_thr": 1e-7,
    "md_type": ("nvt", "or npt, nve, langevin, fire, msst"),
    "md_nstep": ("10000", "number of steps"), "md_dt": ("1", "unit in fs"),
    "md_tfirst": ("300", "unit in K"), "md_tlast": ("300", "unit in K"),
    "md_dumpfreq": 1, "md_restartfreq": 1, "out_level": "h",
    "_section_hints": {"System variables": ["kspacing            0.14 # unit in 1/bohr"]},
    "_comment_hints": {"cal_force": "1", "cal_stress": "1"},
}

TEMPLATE_PW_MD = {
    "calculation": "md", "symmetry": 0, "gamma_only": 1,
    "precision": ("double", "or single"),
    "ecutwfc": 80.0, "pw_diag_nmax": 20, "pw_diag_ndim": 2,
    "basis_type": ("pw", "or lcao"),
    "ks_solver": ("dav_subspace", "or genelpa for lcao basis"),
    "smearing_method": "gauss", "smearing_sigma": 0.01,
    "mixing_type": "broyden", "mixing_beta": 0.8, "scf_nmax": 100, "scf_thr": 1e-8,
    "md_type": ("nvt", "or npt, nve, langevin, fire, msst"),
    "md_nstep": ("10000", "number of steps"), "md_dt": ("1", "unit in fs"),
    "md_tfirst": ("300", "unit in K"), "md_tlast": ("300", "unit in K"),
    "md_dumpfreq": 1, "md_restartfreq": 1, "out_level": "h",
    "_section_hints": {"System variables": ["kspacing            0.14 # unit in 1/bohr"]},
    "_comment_hints": {"cal_force": "1", "cal_stress": "1"},
}

TEMPLATE_DP_MD = {
    "calculation": "md", "esolver_type": "dp",
    "symmetry": 0,
    "cal_force": 1, "cal_stress": 1,
    "pot_file": ("graph-compress.pb", "DeePMD-kit model, v1=.pb, v2=.pt"),
    "md_type": ("nvt", "or npt, nve, langevin, fire, msst"),
    "md_nstep": ("10000", "number of steps"), "md_dt": ("1", "unit in fs"),
    "md_tfirst": ("300", "unit in K"), "md_tlast": ("300", "unit in K"),
    "md_restartfreq": ("md_nstep/10", "write RESTART every N steps"),
    "out_level": "m",
    "_comment_hints": {"dump_force": "1", "dump_vel": "1"},
}


# =============================================================================
# Shared interactive entrypoint
# =============================================================================


def _ask_basis_and_calc(
    console,
    calc_choices: list[tuple[str, str]],
    default_calc_index: int = 0,
) -> tuple[str, str]:
    """Ask the user for basis type and calculation type.

    Returns (basis_type, calculation) — e.g. ("lcao", "cell-relax").
    """
    basis = _prompt_choice(console, "Basis type",
                           ["lcao", "pw"],
                           "lcao")
    labels = [c[0] for c in calc_choices]
    default_label = labels[default_calc_index]
    chosen = _prompt_choice(console, "Calculation type", labels, default_label)
    for label, calc in calc_choices:
        if label == chosen:
            return basis, calc
    return basis, calc_choices[default_calc_index][1]


def _ask_lcao_solver(console, params: InputParams) -> None:
    """For LCAO basis, ask whether the target server is CPU or GPU.

    GPU servers should use cusolver (CUDA-optimised).  CPU servers use
    genelpa (default).  Only called when basis_type is lcao / lcao_in_pw.
    """
    from abacuscopilot.core.standards import is_lcao_basis, solver_for
    if not is_lcao_basis(params.basis_type):
        return
    choice = _prompt_choice(
        console,
        "Target server",
        ["CPU (genelpa)", "GPU (cusolver)"],
        "CPU (genelpa)",
    )
    device = "gpu" if "GPU" in choice else "cpu"
    params.device = device
    params.ks_solver = solver_for(params.basis_type, device)


def _apply_solver_override(console, params: InputParams, solver: str) -> None:
    """Apply a non-interactive --solver override, validated against the标准规范.

    Invalid combinations (e.g. --basis pw --solver genelpa) are auto-corrected
    to the规范 default with a warning, so output is always consistent.
    """
    from abacuscopilot.core.standards import validate_solver
    corrected, warning = validate_solver(params.basis_type, solver, params.device)
    if warning:
        console.print(f"  [yellow]![/yellow] {warning}")
    params.ks_solver = corrected


def _print_summary(console, params: InputParams):
    """Print a compact summary of the key parameters."""
    console.print()
    console.print("[green]✓ INPUT file written successfully.[/green]")
    console.print(f"  Calculation: {params.calculation}, Basis: {params.basis_type}")
    console.print(f"  ecutwfc: {params.ecutwfc} Ry, ks_solver: {params.ks_solver}")
    console.print(f"  SCF thr: {params.scf_thr}, mixing_beta: {params.mixing_beta}")
    if params.calculation in ("relax", "cell-relax"):
        console.print(f"  Relax method: {params.relax_method}, steps: {params.relax_nmax}")
        console.print(f"  Force thr: {params.force_thr_ev} eV/A, Stress thr: {params.stress_thr} kbar")
        if params.calculation == "cell-relax":
            console.print(f"  fixed_axes: {params.fixed_axes}")
    console.print()


def _auto_prepare_files(console, params: InputParams, interactive: bool = True) -> None:
    """Optionally copy pseudopotential/orbital files for the current structure.

    Called after INPUT generation. Reads species from STRU, looks up library
    paths from config, and copies matching files to the current directory.
    After copying, rewrites STRU so its ATOMIC_SPECIES upf (and NUMERICAL_ORBITAL
    orb) filenames match the real library filenames — otherwise STRU may still
    reference e.g. "K.upf" while the copied file is "K_ONCV_PBE-1.0.upf", and
    ABACUS fails to find the pseudopotential (标准规范).
    """
    from abacuscopilot.preprocessing.system_tasks import (
        prepare_calculation_files,
        read_species_from_stru,
    )

    config = load_config()
    libraries = config.get("libraries", {})
    pseudo_lib = libraries.get("pseudo_library", "")
    orbital_lib = libraries.get("orbital_library", "")

    if not pseudo_lib:
        return  # not configured, skip silently

    species = read_species_from_stru("STRU")
    if not species:
        # Try reading from CIF files in current directory
        cif_files = sorted(Path(".").glob("*.cif"))
        if cif_files:
            try:
                from ase.io import read as ase_read
                atoms = ase_read(str(cif_files[0]))
                species = list(set(atoms.get_chemical_symbols()))
            except Exception:
                pass
    if not species:
        return  # nothing to work with, skip

    # Check if files are already present — if so, skip silently
    all_present = True
    for elem in species:
        if not list(Path(".").glob(f"{elem}_*.upf")):
            all_present = False
            break
        if params.basis_type.startswith("lcao") and orbital_lib:
            if not list(Path(".").glob(f"{elem}_*.orb")):
                all_present = False
                break

    if all_present:
        return  # files already copied

    console.print("[bold]Pseudopotential / orbital files:[/bold]")
    console.print(f"  Library: [dim]{pseudo_lib}[/dim]")
    console.print(f"  Elements: {', '.join(species)}")
    console.print()

    result = prepare_calculation_files(
        species, params.basis_type, pseudo_lib, orbital_lib, ".", dry_run=False
    )

    if result["pseudo_files"]:
        console.print("  [green]✓ Pseudopotential:[/green] " + ", ".join(result["pseudo_files"]))
    if result["orbital_files"]:
        console.print("  [green]✓ Orbitals:[/green] " + ", ".join(result["orbital_files"]))
    if result["errors"]:
        for e in result["errors"]:
            console.print(f"  [yellow]![/yellow] {e}")
    console.print()

    # Rewrite STRU so its upf/orb filenames match the real (copied) library
    # files. Only when a STRU exists here (species may have come from a CIF).
    _sync_stru_filenames(console, params, interactive)


def _sync_stru_filenames(console, params: InputParams, interactive: bool) -> None:
    """Update STRU's ATOMIC_SPECIES/NUMERICAL_ORBITAL filenames to real library
    names, per 标准规范. Preserves all other STRU content (m-flags, mag, velocity).

    Interactive: ask whether to overwrite STRU (default yes); if no, save as
    STRU-new. Non-interactive: overwrite in place.
    """
    stru_path = Path("STRU")
    if not stru_path.exists():
        return  # nothing to sync (e.g. species came from a CIF)

    from abacuscopilot.core.standards import is_lcao_basis
    from abacuscopilot.io.stru_file import read_stru
    from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare

    try:
        structure = read_stru("STRU")
    except Exception:
        return  # unreadable STRU — leave it untouched

    is_lcao = is_lcao_basis(params.basis_type)

    # Decide output target
    target = "STRU"
    if interactive:
        choice = _prompt_choice(
            console,
            "Update STRU with resolved pseudopotential/orbital filenames?",
            ["Overwrite STRU", "Save as STRU-new"],
            "Overwrite STRU",
        )
        if "new" in choice.lower():
            target = "STRU-new"

    _write_stru_bare(structure, is_lcao=is_lcao, filepath=target)
    if target == "STRU":
        console.print("  [green]✓ STRU updated[/green] (filenames match library)")
    else:
        console.print(f"  [green]✓ Written {target}[/green] (original STRU unchanged)")
    console.print()


# =============================================================================
# Template lookup
# =============================================================================

_TEMPLATES = {
    ("lcao", "cell-relax"): TEMPLATE_LCAO_CELL_RELAX,
    ("lcao", "relax"): TEMPLATE_LCAO_RELAX,
    ("lcao", "scf"): TEMPLATE_LCAO_SCF,
    ("lcao", "nscf"): TEMPLATE_LCAO_BAND,
    ("lcao", "dos"): TEMPLATE_LCAO_DOS,
    ("lcao", "md"): TEMPLATE_LCAO_MD,
    ("pw", "cell-relax"): TEMPLATE_PW_CELL_RELAX,
    ("pw", "relax"): TEMPLATE_PW_RELAX,
    ("pw", "scf"): TEMPLATE_PW_SCF,
    ("pw", "md"): TEMPLATE_PW_MD,
    ("dp", "md"): TEMPLATE_DP_MD,
    ("pw", "nscf"): TEMPLATE_PW_BAND,
    ("pw", "dos"): TEMPLATE_PW_DOS,
}


def _get_template(basis: str, calc: str) -> dict | None:
    """Look up a template by (basis, calculation)."""
    return _TEMPLATES.get((basis, calc))


def _check_cli_basis_misuse(args, console) -> bool:
    """Return True if user used --pw/--lcao/--dp instead of --basis X,
    or --server instead of --solver."""
    if args:
        for a in args:
            if a in ("--pw", "--lcao", "--dp") and not any(
                    x.startswith("--basis") for x in args):
                console.print(f"[yellow]'{a}' is not a valid flag. Use '--basis {a[2:]}'.[/yellow]")
                console.print(f"[dim]Example: abacuscopilot -task 101 --basis {a[2:]}[/dim]")
                return True
            if a == "--server":
                console.print("[yellow]'--server' is not a valid flag. Use '--solver'.[/yellow]")
                return True
    return False


# =============================================================================
# Task 101: SCF INPUT generation
# =============================================================================

@task(101, category="INPUT", name="SCF INPUT",
      description="Generate INPUT file for self-consistent field calculation",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_scf_input(args: list[str] | None = None, interactive: bool = True,
                   parsed_args=None) -> None:
    """Generate an INPUT file for an SCF calculation using proven defaults."""
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate SCF INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis, _ = _ask_basis_and_calc(console, [("SCF", "scf")])
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
    else:
        basis = "lcao"
    template = _get_template(basis, "scf")

    if template:
        _apply_template(params, template)

    if interactive:
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    _print_summary(console, params)
    _auto_prepare_files(console, params, interactive)


# =============================================================================
# Task 102: Relax INPUT generation
# =============================================================================

@task(102, category="INPUT", name="Relax INPUT",
      description="Generate INPUT file for geometry / cell relaxation",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--relax-type", "type": str, "default": "cell-relax",
           "help": "Relaxation type: cell-relax, relax"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_relax_input(args: list[str] | None = None, interactive: bool = True,
                     parsed_args=None) -> None:
    """Generate an INPUT file for relaxation with proven defaults."""
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate Relax INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis, calc = _ask_basis_and_calc(console, [
            ("cell-relax (atoms + cell)", "cell-relax"),
            ("relax (atoms only)", "relax"),
        ])
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
        calc = getattr(parsed_args, "relax_type", "cell-relax")
    else:
        basis, calc = "lcao", "cell-relax"
    template = _get_template(basis, calc)

    if template:
        _apply_template(params, template)

    if interactive:
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    _print_summary(console, params)
    _auto_prepare_files(console, params, interactive)


# =============================================================================
# Task 103: MD INPUT generation
# =============================================================================

@task(103, category="INPUT", name="MD INPUT",
      description="Generate INPUT file for molecular dynamics",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis: lcao, pw, dp"},
          {"name": "--ensemble", "type": str, "default": "nvt",
           "help": "MD ensemble: nvt, npt, nve, langevin, fire, msst"},
          {"name": "--nstep", "type": int, "default": 10000,
           "help": "Number of MD steps"},
          {"name": "--dt", "type": float, "default": 1.0,
           "help": "Time step (fs)"},
          {"name": "--tfirst", "type": float, "default": 300.0,
           "help": "Initial temperature (K)"},
          {"name": "--tlast", "type": float, "default": 300.0,
           "help": "Final temperature (K)"},
          {"name": "--press", "type": float, "default": 0.0,
           "help": "Target pressure in kbar (NPT only)"},
          {"name": "--pmode", "type": str, "default": "iso",
           "help": "Pressure control mode: iso, aniso, tri"},
          {"name": "--dumpfreq", "type": int, "default": 10,
           "help": "MD_dump output frequency (steps)"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_md_input(args: list[str] | None = None, interactive: bool = True,
                  parsed_args=None) -> None:
    """Generate an INPUT file for MD simulation with proven defaults.
    - scf_thr=1e-5 (looser convergence)
    - symmetry=0, gamma_only=1 (faster MD)
    - chg_extrap=second-order
    """
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate MD INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis = _prompt_choice(console, "Basis type", ["lcao", "pw", "dp (Deep Potential)"],
                               "lcao")
        if basis.startswith("dp"):
            basis = "dp"
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw", "dp") else "lcao"
    else:
        basis = "lcao"
    template = _get_template(basis, "md")
    if template:
        _apply_template(params, template)

    if interactive:
        params.md_type = _prompt_choice(console, "MD ensemble",
                                        ["nvt", "npt", "nve", "langevin", "fire", "msst"],
                                        "nvt")
        if params.md_type == "npt":
            pmode = _prompt_choice(console, "Pressure control mode",
                                   ["iso (isotropic)", "aniso (anisotropic)", "tri (full)"],
                                   "iso (isotropic)")
            params.md_pmode = pmode.split()[0]
            if params.md_pmode == "iso":
                p = float(_prompt(console, "Target pressure (kbar)", "0.0"))
                params.press1 = params.press2 = params.press3 = p
            else:
                p1 = float(_prompt(console, "  press1 (kbar, x-direction)", "0.0"))
                p2 = float(_prompt(console, "  press2 (kbar, y-direction)", "0.0"))
                p3 = float(_prompt(console, "  press3 (kbar, z-direction)", "0.0"))
                params.press1, params.press2, params.press3 = p1, p2, p3
        params.md_nstep = int(_prompt(console, "Number of MD steps", 10000))
        params.md_dt = float(_prompt(console, "Time step (fs)", 1.0))
        params.md_tfirst = float(_prompt(console, "Initial temperature (K)", 300.0))
        params.md_tlast = float(_prompt(console, "Final temperature (K)", 300.0))
        if basis == "dp":
            freq = int(_prompt(console, "MD_dump output frequency (steps)", "100"))
            params.md_dumpfreq = freq
    elif parsed_args:
        params.md_type = parsed_args.ensemble
        params.md_nstep = parsed_args.nstep
        params.md_dt = parsed_args.dt
        params.md_tfirst = parsed_args.tfirst
        params.md_tlast = parsed_args.tlast
        if params.md_type == "npt":
            if parsed_args.press > 0:
                params.press1 = params.press2 = params.press3 = parsed_args.press
            params.md_pmode = parsed_args.pmode
        if basis == "dp":
            params.md_dumpfreq = parsed_args.dumpfreq
    else:
        params.md_type = "nvt"
        params.md_nstep = 10000
        params.md_dt = 1.0
        params.md_tfirst = 300.0
        params.md_tlast = 300.0

    # DP: restartfreq = nstep/10; LCAO/PW: keep template default (1)
    if params.esolver_type == "dp":
        params.md_restartfreq = max(1, params.md_nstep // 10)

    if interactive and params.esolver_type != "dp":
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ MD INPUT file written successfully.[/green]")
    console.print(f"  Ensemble: {params.md_type}, {params.md_nstep} steps, dt={params.md_dt} fs")
    if params.md_type == "npt":
        if params.md_pmode == "iso":
            console.print(f"  Pressure: {params.press1:.1f} kbar, mode: iso")
        else:
            console.print(f"  Pressure: ({params.press1:.1f}, {params.press2:.1f}, {params.press3:.1f}) kbar, mode: {params.md_pmode}")
    console.print(f"  T: {params.md_tfirst}→{params.md_tlast} K")
    if params.esolver_type == "dp":
        model = params.get_param("pot_file")
        console.print(f"  Deep Potential model: {model}")
        if not Path(model).exists():
            console.print(f"  [bold red]! Model file '{model}' not found in current directory![/bold red]")
            console.print("  [bold red]  This file is REQUIRED for DP-MD — place it here before running.[/bold red]")
        # Write STRU-dp: same structure, no pseudo/orbital info needed
        from abacuscopilot.io.stru_file import read_stru
        from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
        try:
            structure = read_stru("STRU")
            _write_stru_bare(structure, is_lcao=False, filepath="STRU-dp", is_dp=True)
            console.print("  [green]✓ STRU-dp written (no pseudo/orbital files needed)[/green]")
            console.print(f"  [dim]Initial velocities: auto-generated from Maxwell-Boltzmann at {params.md_tfirst} K[/dim]")
            removed = []
            for pat in ("*.upf", "*.orb", "KPT"):
                for f in sorted(Path(".").glob(pat)):
                    f.unlink()
                    removed.append(str(f))
            if removed:
                console.print(f"  [dim]Removed: {', '.join(removed)}[/dim]")
        except Exception:
            console.print("  [yellow]! STRU-dp not written (no STRU found)[/yellow]")
    else:
        console.print(f"  Basis: {params.basis_type}")
        console.print(f"  mixing: {params.mixing_type}/beta={params.mixing_beta}, scf_thr={params.scf_thr}")
        # Check STRU has proper upf/orb info; fix if missing
        from abacuscopilot.io.stru_file import read_stru
        from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
        try:
            structure = read_stru("STRU")
            need_fix = False
            for sp in structure.species_order:
                pf = structure.pseudo_files.get(sp, "")
                if not pf or pf == f"{sp}.upf" or "/" not in str(pf):
                    need_fix = True
                    break
            if params.basis_type == "lcao" and not structure.orbital_files:
                need_fix = True
            if need_fix:
                is_lcao = params.basis_type == "lcao"
                _write_stru_bare(structure, is_lcao=is_lcao, filepath="STRU")
                console.print(f"  [green]✓ STRU updated with {'LCAO' if is_lcao else 'PW'} file info[/green]")
            else:
                console.print("  [dim]STRU already has upf/orb info[/dim]")
        except Exception:
            console.print("  [yellow]! STRU not found, skipping file check[/yellow]")
        _auto_prepare_files(console, params, interactive)
    console.print()


# =============================================================================
# Task 104: Band Structure INPUT
# =============================================================================

@task(104, category="INPUT", name="Band Structure INPUT",
      description="Generate INPUT file for band structure (NSCF) calculation",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--proj", "action": "store_true",
           "help": "Output projected bands"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_band_input(args: list[str] | None = None, interactive: bool = True,
                    parsed_args=None) -> None:
    """Generate an INPUT file for band structure calculation with proven defaults."""
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate Band Structure INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis, _ = _ask_basis_and_calc(console, [("Band (NSCF)", "nscf")])
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
    else:
        basis = "lcao"
    template = _get_template(basis, "nscf")

    if template:
        _apply_template(params, template)

    if interactive:
        _ask_lcao_solver(console, params)
        want_proj = _prompt_choice(console, "Output projected bands?",
                                   ["No (band only)", "Yes (band + projection)"],
                                   "No (band only)")
        params.out_proj_band = "Yes" in want_proj
    elif parsed_args:
        if parsed_args.solver:
            _apply_solver_override(console, params, parsed_args.solver)
        if parsed_args.proj:
            params.out_proj_band = True

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ Band structure INPUT file written successfully.[/green]")
    console.print(f"  Calculation: nscf, Basis: {params.basis_type}, out_band: 1"
                  + (", out_proj_band: 1" if params.out_proj_band else ""))
    console.print(f"  symmetry: {params.symmetry} (off for band path)")
    console.print()

    _auto_prepare_files(console, params, interactive)

    # Check for charge density file — show last so it's not buried
    _check_charge_density(console)
    console.print()


# =============================================================================
# Task 105: DOS INPUT
# =============================================================================

@task(105, category="INPUT", name="DOS INPUT",
      description="Generate INPUT file for DOS / PDOS calculation",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--pdos", "action": "store_true",
           "help": "Projected DOS (PDOS) instead of total DOS"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_dos_input(args: list[str] | None = None, interactive: bool = True,
                   parsed_args=None) -> None:
    """Generate an INPUT file for DOS calculation with proven defaults."""
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate DOS INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis, _ = _ask_basis_and_calc(console, [("DOS (NSCF)", "dos")])
        template = _get_template(basis, "dos")
        if template:
            _apply_template(params, template)

        dos_type = _prompt_choice(console, "DOS type",
                                  ["Total DOS", "Projected DOS (PDOS)"],
                                  "Total DOS")
        params.out_dos = 2 if "Projected" in dos_type else 1
    else:
        template = _get_template("lcao", "dos")
        if template:
            _apply_template(params, template)

    if interactive:
        _ask_lcao_solver(console, params)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ DOS INPUT file written successfully.[/green]")
    console.print(f"  DOS type: {'PDOS' if params.out_dos == 2 else 'Total DOS'}")
    console.print(f"  Energy range: [{params.dos_emin_ev}, {params.dos_emax_ev}] eV")
    console.print(f"  Basis: {params.basis_type}, ks_solver: {params.ks_solver}")
    console.print()
    _auto_prepare_files(console, params, interactive)

    _check_charge_density(console)
    console.print()


# =============================================================================
# Task 106: Work Function INPUT
# =============================================================================

@task(106, category="INPUT", name="Work Function INPUT",
      description="Generate INPUT for work function calculation (outputs electrostatic potential)",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_wf_input(args: list[str] | None = None, interactive: bool = True,
                  parsed_args=None) -> None:
    """Generate INPUT for work function calculation.

    Uses the SCF template as baseline, with out_pot enabled.
    """
    console = _get_console()

    if _check_cli_basis_misuse(args, console):
        return

    console.print()
    console.print("[bold cyan]=== Generate Work Function INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis = _prompt_choice(console, "Basis type", ["lcao", "pw"], "lcao")
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
    else:
        basis = "lcao"
    template = _get_template(basis, "scf")

    if template:
        _apply_template(params, template)

    if interactive:
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    params.out_pot = 2  # electrostatic potential (needed for work function)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ Work function INPUT written successfully.[/green]")
    console.print("  out_pot=True (needed for work function analysis)")
    console.print(f"  Basis: {params.basis_type}, ecutwfc: {params.ecutwfc} Ry")
    console.print()
    _auto_prepare_files(console, params, interactive)


# =============================================================================
# Convergence test helpers
# =============================================================================


def _setup_convergence_subdir(
    console,
    subdir: str,
    params: InputParams,
    kpts,
    dry_run: bool = False,
) -> Path:
    """Create a sub-directory with INPUT, STRU, KPT, and support files."""
    import shutil

    dir_path = Path(subdir)
    if not dry_run:
        dir_path.mkdir(parents=True, exist_ok=True)

    for f in ("STRU", "stru"):
        if Path(f).exists() and not dry_run:
            shutil.copy(f, dir_path / "STRU")
            # Always rewrite STRU with properly resolved upf/orb filenames
            from abacuscopilot.io.stru_file import read_stru
            from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
            structure = read_stru(dir_path / "STRU")
            _write_stru_bare(structure, is_lcao=params.basis_type == "lcao",
                             filepath=str(dir_path / "STRU"))

    if not dry_run:
        from abacuscopilot.io.input_file import write_input
        write_input(params, filepath=dir_path / "INPUT")

    if not dry_run and kpts is not None:
        from abacuscopilot.io.kpt_file import write_kpt
        write_kpt(kpts, filepath=dir_path / "KPT")

    if not dry_run:
        from abacuscopilot.config import load_config
        from abacuscopilot.preprocessing.system_tasks import (
            prepare_calculation_files,
            read_species_from_stru,
        )
        species = read_species_from_stru("STRU")
        if species:
            config = load_config()
            libs = config.get("libraries", {})
            pseudo_lib = libs.get("pseudo_library", "")
            orbital_lib = libs.get("orbital_library", "")
            if pseudo_lib:
                prepare_calculation_files(
                    species, params.basis_type, pseudo_lib, orbital_lib,
                    target_dir=str(dir_path), dry_run=False,
                )

    return dir_path


def _get_sub_script_path() -> str | None:
    """Return the user-configured submission script path, or None."""
    config = load_config()
    path = config.get("paths", {}).get("sub_script", "")
    if path:
        expanded = str(Path(path).expanduser())
        if Path(expanded).exists():
            return expanded
    return None


def _copy_sub_script(target_dir: Path, script_path: str | None, console) -> None:
    """Copy submission script to target_dir if configured, otherwise remind."""
    import shutil
    if script_path and Path(script_path).exists():
        dest = target_dir / Path(script_path).name
        shutil.copy(script_path, dest)
        (target_dir / dest.name).chmod(0o755)
    else:
        console.print("  [yellow]![/yellow] No submission script configured (task 1506 to set)")


def _build_sub_script(name: str) -> str:
    """Minimal SLURM submission script — kept as fallback (no longer used by default)."""
    return f"""#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source ~/softwares/abacus-develop-LTSv3.10.0/toolchain/abacus_env.sh
export OMP_NUM_THREADS=1

srun abacus > {name}.log 2>&1
"""


# =============================================================================
# Task 107: Ecutwfc convergence test
# =============================================================================

@task(107, category="INPUT", name="Ecutwfc Test",
      description="Generate INPUTs for ecutwfc convergence test (PW only)")
def task_ecutwfc_test(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a series of INPUT files sweeping ecutwfc (PW basis only)."""
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== Ecutwfc Convergence Test ===[/bold cyan]")
    console.print()

    # Basis type — LCAO makes no sense for this test
    basis = _prompt_choice(console, "Basis type", ["pw", "lcao"], "pw") if interactive else "pw"
    if basis == "lcao":
        console.print()
        console.print("[yellow]For LCAO, ecutwfc should match the orbital file's cutoff (e.g. 100Ry).[/yellow]")
        console.print("[yellow]Convergence testing is not meaningful without multiple orbital libraries.[/yellow]")
        console.print("[dim]Use Kspacing Test for LCAO convergence testing instead.[/dim]")
        console.print()
        return

    params = InputParams()
    params.suffix = "ABACUS"
    if Path("INPUT").exists():
        try:
            from abacuscopilot.io.input_file import read_input
            params = read_input("INPUT")
            console.print("  [dim]Loaded base INPUT[/dim]")
        except Exception:
            _apply_template(params, TEMPLATE_PW_SCF)
            console.print("  [dim]Using PW SCF defaults[/dim]")
    else:
        _apply_template(params, TEMPLATE_PW_SCF)
        console.print("  [dim]Using PW SCF defaults[/dim]")

    kpts = None
    if Path("STRU").exists():
        from abacuscopilot.io.stru_file import read_stru
        structure = read_stru("STRU")
        from abacuscopilot.io.kpt_file import auto_mp_kpts
        kpts = auto_mp_kpts(structure.lattice, kspacing=0.14)

    e_start = float(_prompt(console, "Start ecutwfc (Ry)", "40"))
    e_end = float(_prompt(console, "End ecutwfc (Ry)", "100"))
    e_step = float(_prompt(console, "Step (Ry)", "10"))
    ecut_values = np.arange(e_start, e_end + 0.1, e_step)

    console.print(f"\n  {len(ecut_values)} directories: ecutwfc = {', '.join(f'{v:.0f}' for v in ecut_values)} Ry")
    console.print()

    for v in ecut_values:
        params.ecutwfc = v
        _setup_convergence_subdir(console, f"ecutwfc_{v:.0f}", params, kpts, dry_run=True)
        console.print(f"  [dim]ecutwfc_{v:.0f}/[/dim]")

    confirm = _prompt(console, "Proceed? (y/n)", "y")
    if confirm.lower() not in ("y", "yes"):
        console.print("[yellow]Cancelled.[/yellow]")
        return

    sub_script = _get_sub_script_path()
    for v in ecut_values:
        params.ecutwfc = v
        d = _setup_convergence_subdir(console, f"ecutwfc_{v:.0f}", params, kpts)
        _copy_sub_script(d, sub_script, console)
        console.print(f"  [green]✓[/green] ecutwfc_{v:.0f}/")

    console.print()
    console.print(f"[green]✓ {len(ecut_values)} directories ready.[/green]")
    if sub_script:
        console.print("  [dim]Submit: for d in ecutwfc_*/; do cd $d && sbatch *.sh && cd ..; done[/dim]")
    console.print()


# =============================================================================
# Task 108: Kspacing convergence test
# =============================================================================

@task(108, category="INPUT", name="Kspacing Test",
      description="Generate INPUTs for kspacing convergence test (PW & LCAO)")
def task_kspacing_test(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a series of INPUT files sweeping kspacing."""
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== Kspacing Convergence Test ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"
    if Path("INPUT").exists():
        try:
            from abacuscopilot.io.input_file import read_input
            params = read_input("INPUT")
            console.print("  [dim]Loaded base INPUT[/dim]")
        except Exception:
            _apply_template(params, TEMPLATE_LCAO_SCF)
    else:
        _apply_template(params, TEMPLATE_LCAO_SCF)

    if not Path("STRU").exists():
        console.print("[red]No STRU found.[/red]")
        return
    from abacuscopilot.io.stru_file import read_stru
    structure = read_stru("STRU")

    if interactive:
        _ask_lcao_solver(console, params)

    ks_start = float(_prompt(console, "Start kspacing (2pi/A)", "0.30"))
    ks_end = float(_prompt(console, "End kspacing (2pi/A)", "0.06"))
    ks_step = float(_prompt(console, "Step (2pi/A)", "-0.02"))

    ks_values = []
    v = ks_start
    while (ks_step < 0 and v >= ks_end - 0.0001) or (ks_step > 0 and v <= ks_end + 0.0001):
        ks_values.append(round(v, 4))
        v += ks_step

    # Pre-compute grids and detect duplicates
    from abacuscopilot.io.kpt_file import auto_mp_kpts
    ks_grids = []
    for v in ks_values:
        kts = auto_mp_kpts(structure.lattice, kspacing=v)
        ks_grids.append((v, kts))

    # Check for duplicate grids
    seen_grids = {}
    duplicates = set()
    for v, kts in ks_grids:
        grid_key = kts.grid
        if grid_key in seen_grids:
            duplicates.add(grid_key)
        seen_grids[grid_key] = v  # last (tightest) kspacing for this grid

    filtered = ks_grids
    if duplicates:
        console.print(f"\n  [yellow]{len(duplicates)} duplicate grid(s) detected:[/yellow]")
        for g in duplicates:
            dup_vals = [f"{v:.3f}" for v, kts in ks_grids if kts.grid == g]
            console.print(f"    {g[0]}x{g[1]}x{g[2]}: kspacing = {', '.join(dup_vals)} → keep {seen_grids[g]:.3f}")

        if interactive:
            skip = _prompt(console, "Skip duplicates? (y/n)", "y")
            if skip.lower() not in ("n", "no"):
                filtered = [(v, kts) for v, kts in ks_grids
                            if kts.grid not in duplicates or abs(v - seen_grids[kts.grid]) < 0.0001]
                console.print(f"  [dim]Kept {len(filtered)}/{len(ks_values)} unique grids[/dim]")

    console.print()
    for v, kts in filtered:
        _setup_convergence_subdir(console, f"kspacing_{v:.3f}", params, kts, dry_run=True)
        console.print(f"  [dim]kspacing_{v:.3f}/  {kts.grid[0]}x{kts.grid[1]}x{kts.grid[2]}[/dim]")
    console.print()

    confirm = _prompt(console, "Proceed? (y/n)", "y")
    if confirm.lower() not in ("y", "yes"):
        console.print("[yellow]Cancelled.[/yellow]")
        return

    sub_script = _get_sub_script_path()
    for v, kts in filtered:
        params.kspacing = v
        d = _setup_convergence_subdir(console, f"kspacing_{v:.3f}", params, kts)
        _copy_sub_script(d, sub_script, console)
        console.print(f"  [green]✓[/green] kspacing_{v:.3f}/  {kts.grid[0]}x{kts.grid[1]}x{kts.grid[2]}")

    console.print()
    console.print(f"[green]✓ {len(filtered)} directories ready.[/green]")
    if sub_script:
        console.print("  [dim]Submit: for d in kspacing_*/; do cd $d && sbatch *.sh && cd ..; done[/dim]")
    console.print()


# =============================================================================
# Task 109: Convergence analysis
# =============================================================================


def _detect_conv_dirs() -> tuple[str, list[tuple[float, Path]]]:
    """Find ecutwfc_*/ or kspacing_*/ convergence test directories.

    Returns (param_name, [(value, dir_path), ...]) sorted by value.
    """
    cwd = Path(".")
    for prefix in ("ecutwfc_", "kspacing_"):
        dirs = sorted(cwd.glob(f"{prefix}*"))
        if dirs:
            pairs = []
            for d in dirs:
                try:
                    val = float(d.name.replace(prefix, "").replace("_", "."))
                except ValueError:
                    continue
                pairs.append((val, d))
            if pairs:
                return prefix.rstrip("_"), sorted(pairs)
    return "", []


def _check_convergence(out_dir: Path) -> dict:
    """Check convergence of a single calculation directory."""
    from abacuscopilot.postprocessing.scf_tasks import _parse_calculation_status
    return _parse_calculation_status(out_dir)


@task(109, category="INPUT", name="Convergence Analysis",
      description="Check convergence test results and plot energy vs parameter")
def task_conv_analysis(args: list[str] | None = None, interactive: bool = True) -> None:
    """Scan ecutwfc_*/ or kspacing_*/ directories and plot convergence."""
    console = _get_console()
    console.print()
    console.print("[bold cyan]=== Convergence Analysis ===[/bold cyan]")
    console.print()

    param_name, dirs = _detect_conv_dirs()
    if not dirs:
        console.print("[red]No ecutwfc_*/ or kspacing_*/ directories found.[/red]")
        console.print("[dim]Run task 107 or 108 first to generate test directories.[/dim]")
        return

    console.print(f"  Found {len(dirs)} {param_name} test directories")

    # Check each directory
    from abacuscopilot.core.constants import RY_TO_EV

    # Read atom count from log (authoritative — what ABACUS actually used)
    natom = 0

    results = []
    for val, d in dirs:
        out_dir = d / "OUT.ABACUS"
        status = _check_convergence(out_dir)
        e_ry = status.get("final_energy_ry")
        e_ev = e_ry * RY_TO_EV if e_ry is not None else None
        # Use natom from first valid log
        if natom == 0 and status.get("natom"):
            natom = status["natom"]
        e_per_atom = e_ev / natom if e_ev is not None and natom > 0 else None
        results.append({
            "value": val,
            "dir": d,
            "completed": status["completed"],
            "converged": status["converged"],
            "energy_ry": e_ry,
            "energy_ev": e_per_atom,
            "n_scf_steps": status["n_scf_steps"],
            "wall_time": status.get("wall_time", "—") or "—",
        })

    if natom > 0:
        console.print(f"  [dim]Atoms: {natom} (from log, energies per atom)[/dim]")

    # Summary table
    from rich.table import Table
    table = Table(title=f"{param_name} Convergence Summary")
    table.add_column(f"{param_name}", justify="right")
    table.add_column("SCF Steps", justify="right")
    table.add_column("Wall Time", justify="right")
    table.add_column("Energy (eV/atom)", justify="right")
    table.add_column("Status")
    for r in results:
        icon = "[green]✓[/green]" if r["converged"] else "[red]✗[/red]"
        e_str = f"{r['energy_ev']:.6f}" if r["energy_ev"] is not None else "—"
        wt_str = r["wall_time"].strip() if r["wall_time"] else "—"
        table.add_row(f"{r['value']:.4f}", str(r["n_scf_steps"]), wt_str, e_str, icon)
    console.print()
    console.print(table)
    console.print()

    # Check if all converged
    converged_count = sum(1 for r in results if r["converged"])
    if converged_count == 0:
        console.print("[red]No calculations converged. Cannot plot.[/red]")
        return

    if interactive:
        do_plot = _prompt(console, "Plot convergence curve? (y/n)", "y")
        if do_plot.lower() not in ("y", "yes"):
            return

    # Filter converged results
    conv_results = [r for r in results if r["converged"] and r["energy_ev"] is not None]
    if len(conv_results) < 2:
        console.print("[red]Need at least 2 converged points to plot.[/red]")
        return

    xs = [r["value"] for r in conv_results]
    ys = [r["energy_ev"] for r in conv_results]

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from abacuscopilot.plotting.style import load_style_from_config
    load_style_from_config()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(xs, ys, "o-", color="#1f77b4", linewidth=1.2, markersize=5)
    ax.set_xlabel(param_name)
    ax.set_ylabel("Energy (eV/atom)")
    ax.set_title(f"{param_name} Convergence")
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.6f}"))
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis="both", direction="out")

    save_name = f"{param_name}_convergence.png"
    fig.tight_layout(pad=1.2)
    fig.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Also save data
    dat_name = f"{param_name}_convergence.dat"
    with open(dat_name, "w") as f:
        f.write(f"# {param_name}  Energy(eV/atom)\n")
        for v, e in zip(xs, ys):
            f.write(f"{v:.6f}  {e:.8f}\n")

    console.print(f"  [green]✓ Plot: {save_name}[/green]")
    console.print(f"  [green]✓ Data: {dat_name}[/green]")

    # Estimate converged value (minimum energy)
    imin = np.argmin(ys)
    console.print(f"  Lowest energy: {ys[imin]:.6f} eV/atom at {param_name} = {xs[imin]:.4f}")
    console.print()


# =============================================================================
# Task 110: NEB INPUT
# =============================================================================

TEMPLATE_LCAO_NEB = {
    "calculation": "nscf", "symmetry": 1,
    "kspacing": 0.14, "precision": "double",
    "ecutwfc": 100.0, "basis_type": "lcao",
    "ks_solver": "genelpa", "smearing_method": "gauss",
    "smearing_sigma": 0.01, "mixing_type": "broyden",
    "mixing_beta": 0.8, "scf_nmax": 100, "scf_thr": 1e-7,
    "nbands": 256,
}

TEMPLATE_PW_NEB = {
    "calculation": "nscf", "symmetry": 1,
    "kspacing": 0.14, "precision": "double",
    "ecutwfc": 80.0, "pw_diag_nmax": 20, "pw_diag_ndim": 2,
    "basis_type": "pw", "ks_solver": "dav_subspace",
    "smearing_method": "gauss", "smearing_sigma": 0.01,
    "mixing_type": "broyden", "mixing_beta": 0.8,
    "scf_nmax": 100, "scf_thr": 1e-8, "nbands": 256,
}

# Register NEB templates
_TEMPLATES[("lcao", "neb")] = TEMPLATE_LCAO_NEB
_TEMPLATES[("pw", "neb")] = TEMPLATE_PW_NEB


@task(110, category="INPUT", name="NEB INPUT",
      description="Generate INPUT file for Nudged Elastic Band calculation",
      cli_args=[
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--images", "type": int, "default": 5,
           "help": "Number of NEB images"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_neb_input(args: list[str] | None = None, interactive: bool = True,
                   parsed_args=None) -> None:
    """Generate an INPUT file for NEB calculation."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Generate NEB INPUT ===[/bold cyan]")
    console.print()

    params = InputParams()
    params.suffix = "ABACUS"

    if interactive:
        basis = _prompt_choice(console, "Basis type", ["lcao", "pw"], "lcao")
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
    else:
        basis = "lcao"
    template = _get_template(basis, "neb")
    if template:
        _apply_template(params, template)

    if interactive:
        n_images = int(_prompt(console, "Number of NEB images", "5"))
    elif parsed_args:
        n_images = parsed_args.images
    else:
        n_images = 5

    if interactive:
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    from abacuscopilot.io.input_file import write_input
    write_input(params)

    console.print()
    console.print("[green]✓ NEB INPUT file written successfully.[/green]")
    console.print(f"  Calculation: nscf, Basis: {params.basis_type}")
    console.print(f"  NEB images: {n_images}")
    console.print("  [dim]Place 00/ → NN/ subdirectories with POSCAR files[/dim]")
    console.print()
    _auto_prepare_files(console, params, interactive)
