"""Phonon setup tasks for ABACUS + Phonopy workflow.

Task 112: Generate displaced supercells and INPUT files for phonon calculation.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from abacuscopilot.config import load_config
from abacuscopilot.console_utils import _get_console, _prompt, _prompt_choice
from abacuscopilot.core.models import InputParams
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


def _find_phonopy_init() -> str | None:
    """Return path to phonopy-init (for v4+), or None."""
    import shutil as _shutil
    path = _shutil.which("phonopy-init")
    if path:
        return path
    # phonopy v4+ also supports python -m phonopy --init
    return None


def _run_phonopy_cmd(args: list[str], cwd: str | Path = ".") -> subprocess.CompletedProcess:
    """Run phonopy with the given arguments."""
    exe = _find_phonopy()
    if exe is None:
        raise FileNotFoundError("phonopy not found")
    if exe == sys.executable:
        cmd = [sys.executable, "-m", "phonopy"] + args
    else:
        cmd = [exe] + args
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=120,
    )


def _run_phonopy_init(args: list[str], cwd: str | Path = ".") -> subprocess.CompletedProcess:
    """Run phonopy-init -d ... (phonopy v4+)."""
    import shutil as _shutil
    # Try all possible locations for phonopy-init
    for candidate in [
        _shutil.which("phonopy-init"),
        Path(sys.executable).parent / "phonopy-init",
        Path(sys.prefix) / "bin" / "phonopy-init",
    ]:
        if candidate and Path(str(candidate)).exists():
            cmd = [str(candidate)] + args
            break
    else:
        # Last resort: python -m phonopy --init
        cmd = [sys.executable, "-m", "phonopy", "--init"] + args
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(cwd), timeout=120,
    )


@task(112, category="INPUT", name="Phonon Setup",
      description="Generate displaced supercells (phonopy -d) and SCF INPUT files for phonon calculation",
      cli_args=[
          {"name": "--dim", "type": str, "default": "2 2 2",
           "help": "Supercell dimensions (e.g. '2 2 2')"},
          {"name": "--basis", "type": str, "default": "lcao",
           "help": "Basis type: lcao, pw"},
          {"name": "--solver", "type": str, "default": "",
           "help": "LCAO solver: genelpa (CPU) or cusolver (GPU)"},
      ])
def task_phonon_setup(args: list[str] | None = None, interactive: bool = True,
                      parsed_args=None) -> None:
    """Generate displaced supercells via phonopy and create disp-XXX directories.

    Requires phonopy ≥ 2.19.1 (``pip install phonopy``).
    """
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Phonon Setup (Phonopy) ===[/bold cyan]")
    console.print("[dim]Finite-displacement method: phonopy -d --abacus[/dim]")
    console.print()

    # --- Check phonopy ---
    if _find_phonopy() is None:
        console.print("[red]Phonopy is not installed.[/red]")
        console.print("[dim]Install: pip install phonopy[/dim]")
        return

    # --- Supercell dimensions ---
    if interactive:
        dim_s = _prompt(console, "Supercell dimensions", "2 2 2")
    elif parsed_args:
        dim_s = parsed_args.dim
    else:
        dim_s = "2 2 2"
    dim_parts = dim_s.strip().split()
    if len(dim_parts) != 3:
        console.print("[red]Dimensions must be 3 integers, e.g. '2 2 2'[/red]")
        return
    console.print(f"  Supercell: {dim_parts[0]}×{dim_parts[1]}×{dim_parts[2]}")

    # --- Basis type ---
    from abacuscopilot.preprocessing.input_tasks import (
        _apply_template, _ask_lcao_solver, _apply_solver_override,
        _auto_prepare_files, _get_template,
    )

    if interactive:
        basis = _prompt_choice(console, "Basis type", ["lcao", "pw"], "lcao")
    elif parsed_args:
        basis = parsed_args.basis if parsed_args.basis in ("lcao", "pw") else "lcao"
    else:
        basis = "lcao"

    # --- INPUT params from SCF template ---
    params = InputParams()
    params.suffix = "ABACUS"
    template = _get_template(basis, "scf")
    if template:
        _apply_template(params, template)
    # Phonon: SCF with forces, no relaxation
    params.calculation = "scf"
    params.cal_force = 1
    if "_template_keys" not in params.extras:
        params.extras["_template_keys"] = []
    if "cal_force" not in params.extras["_template_keys"]:
        params.extras["_template_keys"].append("cal_force")

    # --- LCAO solver ---
    if interactive and basis == "lcao":
        _ask_lcao_solver(console, params)
    elif parsed_args and parsed_args.solver:
        _apply_solver_override(console, params, parsed_args.solver)

    # --- D3 / functional ---
    if interactive:
        use_d3 = _prompt_choice(console, "D3 dispersion correction", ["No", "d3_0 (zero-damping)", "d3_bj (Becke-Johnson)"], "No")
        if "d3_0" in use_d3:
            params.vdw_method = "d3_0"
        elif "d3_bj" in use_d3:
            params.vdw_method = "d3_bj"

    if interactive:
        use_func = _prompt_choice(console, "Exchange-correlation functional", ["PBEsol", "PBE"], "PBEsol")
        if "PBEsol" in use_func:
            params.dft_functional = "pbesol"

    # --- Read STRU and prepare files ---
    from abacuscopilot.io.stru_file import read_stru, write_stru

    stru_path = Path("STRU")
    if not stru_path.exists():
        console.print("[red]No STRU file found in current directory.[/red]")
        return

    structure = read_stru(stru_path)
    _auto_prepare_files(console, params, interactive)

    # Gather UPF/ORB files
    pseudo_files: list[Path] = []
    orbital_files: list[Path] = []
    for species in structure.species_order:
        for pat in Path(".").glob(f"{species}_*.upf"):
            if pat not in pseudo_files:
                pseudo_files.append(pat)
        if basis == "lcao":
            for pat in Path(".").glob(f"{species}_*.orb"):
                if pat not in orbital_files:
                    orbital_files.append(pat)

    # --- Run phonopy -d ---
    console.print()
    console.print(f"  Running: [dim]phonopy-init -d --dim=\"{dim_s}\" --abacus[/dim]")
    result = _run_phonopy_init(["-d", f"--dim={dim_s}", "--abacus"])
    if result.returncode != 0:
        console.print("[red]phonopy-init -d failed:[/red]")
        console.print(result.stderr)
        return
    if result.stderr:
        # phonopy prints info to stderr
        for line in result.stderr.strip().split("\n"):
            console.print(f"  [dim]{line}[/dim]")

    # Delete STRU.in (phonopy supercell template, redundant with STRU-001)
    stru_in = Path("STRU.in")
    if stru_in.exists():
        stru_in.unlink()

    # Find generated STRU-XXX files
    stru_files = sorted(Path(".").glob("STRU-*"))
    if not stru_files:
        console.print("[red]No STRU-* files generated by phonopy.[/red]")
        console.print("[dim]Check that STRU is valid and phonopy ≥ 2.19.1 is installed.[/dim]")
        return

    console.print(f"  Displacements: {len(stru_files)}")

    # --- Create disp-XXX directories ---
    console.print()
    console.print("[bold]Creating displacement directories:[/bold]")

    from abacuscopilot.io.input_file import write_input

    for i, sf in enumerate(stru_files, 1):
        dir_name = f"disp-{i:03d}"
        dir_path = Path(dir_name)
        dir_path.mkdir(exist_ok=True)

        # Copy displaced STRU
        shutil.copy2(sf, dir_path / "STRU")

        # Write INPUT
        write_input(params, dir_path / "INPUT")

        # Copy pseudopotential + orbital files
        for pf in pseudo_files:
            shutil.copy2(pf, dir_path / pf.name)
        for of in orbital_files:
            shutil.copy2(of, dir_path / of.name)

        # Copy sub.abacus if present
        sub_path = Path("sub.abacus")
        if sub_path.exists():
            shutil.copy2(sub_path, dir_path / "sub.abacus")

        console.print(f"  [green]✓[/green] {dir_name}/  ({sf.name})")

    # --- Cleanup ---
    for pf in pseudo_files:
        pf.unlink(missing_ok=True)
    for of in orbital_files:
        of.unlink(missing_ok=True)
    local_sub = Path("sub.abacus")
    if local_sub.exists():
        local_sub.unlink()

    # Save setup info for task 1501
    import json as _json
    _setup_info = {
        "dim": dim_s,
        "basis": basis,
        "ks_solver": params.ks_solver,
        "n_displacements": len(stru_files),
    }
    with open("phonopy_setup.json", "w") as _f:
        _json.dump(_setup_info, _f, indent=2)
    console.print("  [green]✓ phonopy_setup.json[/green]")

    # Keep STRU-* files and phonopy_disp.yaml (needed by phonopy later)
    console.print()
    console.print("  [bold yellow]⚠  For accurate phonon results, use a well-relaxed[/bold yellow]")
    console.print("  [bold yellow]    primitive cell as the starting STRU.[/bold yellow]")
    console.print("  [dim]    Run 'abacuscopilot -task 502' to convert to primitive cell first.[/dim]")
    console.print()
    console.print(f"[green]✓ Phonon setup complete: {len(stru_files)} displacement directories[/green]")
    console.print(f"  Basis: {basis}, ks_solver: {params.ks_solver}")
    console.print(f"  Run all disp-*/ tasks, then: [bold]abacuscopilot -task 1501[/bold]")
    console.print("  [dim]STRU-* and phonopy_disp.yaml kept for phonopy post-processing.[/dim]")
    console.print()
