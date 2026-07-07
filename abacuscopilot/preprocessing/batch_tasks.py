"""Batch job submission tasks for ABACUS calculations.

Task IDs 801-899

Generates PBS/Torque and SLURM job scripts from templates, validates
INPUT+STRU+KPT consistency before submission, and supports dry-run mode.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

from abacuscopilot.console_utils import _get_console, _prompt
from abacuscopilot.tasks import task

# =============================================================================
# Job script templates
# =============================================================================

PBS_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -l walltime={walltime}
#PBS -q {queue}
#PBS -j oe
#PBS -o {job_name}.out
#PBS -V

cd $PBS_O_WORKDIR

{env_setup}

mpirun -np {total_cores} {abacus_binary} > {job_name}.log 2>&1
"""

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --time={walltime}
#SBATCH --partition={partition}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err

cd $SLURM_SUBMIT_DIR

{env_setup}

srun {abacus_binary} > {job_name}.log 2>&1
"""


def _load_job_config() -> dict:
    """Load job submission settings from config."""
    from abacuscopilot.config import load_config
    config = load_config()
    return config.get("job", {})


def _get_abacus_binary() -> str:
    """Get the ABACUS binary path from config."""
    from abacuscopilot.config import load_config
    config = load_config()
    return config.get("paths", {}).get("abacus_binary", "abacus")


# =============================================================================
# Task 801: Generate PBS job script
# =============================================================================


@task(601, category="Batch", name="PBS Script",
      description="Generate PBS/Torque job submission script from template")
def task_pbs_script(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a PBS/Torque job script for ABACUS."""
    console = _get_console()
    job_config = _load_job_config()
    abacus_bin = _get_abacus_binary()

    console.print()
    console.print("[bold cyan]=== Generate PBS Job Script ===[/bold cyan]")
    console.print()

    # Collect parameters
    if interactive:
        job_name = _prompt(console, "Job name", "ABACUS")
        default_nodes = job_config.get("pbs_nodes", "1")
        nodes = _prompt(console, "Number of nodes", str(default_nodes))
        default_ppn = job_config.get("pbs_ppn", "32")
        ppn = _prompt(console, "Processors per node", str(default_ppn))
        walltime = _prompt(console, "Walltime (HH:MM:SS)", "48:00:00")
        queue = _prompt(console, "Queue name",
                        job_config.get("pbs_queue", "batch"))
        env_setup = _prompt(console, "Environment setup commands (optional)",
                            job_config.get("env_setup", ""))
    else:
        job_name = job_config.get("job_name", "ABACUS")
        nodes = str(job_config.get("pbs_nodes", 1))
        ppn = str(job_config.get("pbs_ppn", 32))
        walltime = job_config.get("walltime", "48:00:00")
        queue = job_config.get("pbs_queue", "batch")
        env_setup = job_config.get("env_setup", "")

    total_cores = str(int(nodes) * int(ppn))

    # Generate script
    script = PBS_TEMPLATE.format(
        job_name=job_name,
        nodes=nodes,
        ppn=ppn,
        walltime=walltime,
        queue=queue,
        total_cores=total_cores,
        abacus_binary=abacus_bin,
        env_setup=env_setup,
    )

    script_path = f"{job_name}.pbs"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)

    console.print()
    console.print(f"[green]✓ PBS script written to {script_path}[/green]")
    console.print(f"  Nodes: {nodes}, PPN: {ppn}, Cores: {total_cores}")
    console.print(f"  Walltime: {walltime}, Queue: {queue}")
    console.print()
    console.print(f"  [dim]Submit with: qsub {script_path}[/dim]")
    console.print()


# =============================================================================
# Task 802: Generate SLURM job script
# =============================================================================


@task(602, category="Batch", name="SLURM Script",
      description="Generate SLURM job submission script from template")
def task_slurm_script(args: list[str] | None = None, interactive: bool = True) -> None:
    """Generate a SLURM job script for ABACUS."""
    console = _get_console()
    job_config = _load_job_config()
    abacus_bin = _get_abacus_binary()

    console.print()
    console.print("[bold cyan]=== Generate SLURM Job Script ===[/bold cyan]")
    console.print()

    if interactive:
        job_name = _prompt(console, "Job name", "ABACUS")
        default_nodes = job_config.get("slurm_nodes", "1")
        nodes = _prompt(console, "Number of nodes", str(default_nodes))
        default_nt = job_config.get("slurm_ntasks_per_node", "32")
        ntasks = _prompt(console, "Tasks per node", str(default_nt))
        default_cpt = job_config.get("slurm_cpus_per_task", "1")
        cpus = _prompt(console, "CPUs per task", str(default_cpt))
        walltime = _prompt(console, "Walltime (HH:MM:SS)", "48:00:00")
        partition = _prompt(console, "Partition",
                            job_config.get("slurm_partition", "compute"))
        env_setup = _prompt(console, "Environment setup (optional)",
                            job_config.get("env_setup", ""))
    else:
        job_name = job_config.get("job_name", "ABACUS")
        nodes = str(job_config.get("slurm_nodes", 1))
        ntasks = str(job_config.get("slurm_ntasks_per_node", 32))
        cpus = str(job_config.get("slurm_cpus_per_task", 1))
        walltime = job_config.get("walltime", "48:00:00")
        partition = job_config.get("slurm_partition", "compute")
        env_setup = job_config.get("env_setup", "")

    script = SLURM_TEMPLATE.format(
        job_name=job_name,
        nodes=nodes,
        ntasks_per_node=ntasks,
        cpus_per_task=cpus,
        walltime=walltime,
        partition=partition,
        abacus_binary=abacus_bin,
        env_setup=env_setup,
    )

    script_path = f"{job_name}.slurm"
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)

    console.print()
    console.print(f"[green]✓ SLURM script written to {script_path}[/green]")
    console.print(f"  Nodes: {nodes}, Tasks/node: {ntasks}, CPUs/task: {cpus}")
    console.print(f"  Walltime: {walltime}, Partition: {partition}")
    console.print()
    console.print(f"  [dim]Submit with: sbatch {script_path}[/dim]")
    console.print()


# =============================================================================
# Task 803: Validate calculation inputs
# =============================================================================


def _check_required_files() -> dict[str, Any]:
    """Check that required input files exist and are consistent."""
    result = {
        "ok": True,
        "issues": [],
        "warnings": [],
    }

    # Check INPUT exists
    if not Path("INPUT").exists():
        result["ok"] = False
        result["issues"].append("INPUT file is missing")
    else:
        try:
            from abacuscopilot.io.input_file import read_input, validate_input
            params = read_input("INPUT")
            warnings = validate_input(params)
            for w in warnings:
                result["warnings"].append(w)
        except Exception as e:
            result["warnings"].append(f"Could not validate INPUT: {e}")

    # Check STRU exists
    if not Path("STRU").exists():
        result["ok"] = False
        result["issues"].append("STRU file is missing")
    else:
        try:
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru("STRU")
            # Check species match
            species_in_stru = set(structure.species_order)
            # Check pseudo files exist
            missing_pp = []
            for sp, pp_file in structure.pseudo_files.items():
                if not Path(pp_file).exists():
                    missing_pp.append(f"{sp}: {pp_file}")
            if missing_pp:
                result["warnings"].append(
                    f"Missing pseudopotential files: {', '.join(missing_pp)}"
                )
        except Exception as e:
            result["warnings"].append(f"Could not validate STRU: {e}")

    # Check KPT exists
    if not Path("KPT").exists():
        result["warnings"].append("KPT file is missing — ABACUS will use gamma-point only")
    else:
        try:
            from abacuscopilot.io.kpt_file import read_kpt
            kpts = read_kpt("KPT")
            if kpts.mode in ("gamma", "mp") and kpts.grid:
                nk = kpts.grid[0] * kpts.grid[1] * kpts.grid[2]
                if nk == 1:
                    result["warnings"].append("KPT mesh is 1×1×1 (gamma-point only)")
        except Exception as e:
            result["warnings"].append(f"Could not validate KPT: {e}")

    return result


@task(603, category="Batch", name="Validate Inputs",
      description="Validate INPUT+STRU+KPT consistency before job submission")
def task_validate_inputs(args: list[str] | None = None, interactive: bool = True) -> None:
    """Check input files before submitting a calculation."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Validate Calculation Inputs ===[/bold cyan]")
    console.print()

    result = _check_required_files()

    if result["ok"] and not result["issues"]:
        console.print("[green]✓ All required input files found.[/green]")
    else:
        for issue in result["issues"]:
            console.print(f"[red]✗[/red] {issue}")

    if result["warnings"]:
        console.print()
        for w in result["warnings"]:
            console.print(f"[yellow]![/yellow] {w}")

    if not result["issues"] and not result["warnings"]:
        console.print("[green]✓ Ready for submission.[/green]")

    console.print()


# =============================================================================
# Task 804: Batch submit sub-directories
# =============================================================================


def _detect_scheduler(script_path: Path) -> str | None:
    """Detect SLURM or PBS from script content."""
    content = script_path.read_text()
    if "SBATCH" in content.upper() or "srun" in content:
        return "sbatch"
    if "PBS" in content.upper() or "qsub" in content:
        return "qsub"
    return None


@task(604, category="Batch", name="Batch Submit",
      description="Submit all sub-directories with the configured submission script")
def task_batch_submit(args: list[str] | None = None, interactive: bool = True) -> None:
    """Scan subdirectories for submission scripts and submit all."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Batch Submit ===[/bold cyan]")
    console.print()

    from abacuscopilot.config import load_config
    config = load_config()
    sub_script = config.get("paths", {}).get("sub_script", "")
    if not sub_script:
        console.print("[red]No submission script configured.[/red]")
        console.print("[dim]Use System → Set Submit Script to configure one.[/dim]")
        return

    script_name = Path(sub_script).name
    console.print(f"  Looking for [dim]{script_name}[/dim] in subdirectories")

    jobs = []
    for d in sorted(Path(".").iterdir()):
        if not d.is_dir():
            continue
        s = d / script_name
        if s.exists():
            jobs.append((d, s, _detect_scheduler(s)))

    if not jobs:
        console.print(f"\n[yellow]No subdirectories contain {script_name}.[/yellow]")
        return

    console.print(f"\n  Found {len(jobs)} job(s):")
    for d, s, sch in jobs:
        console.print(f"    {d.name}/  →  {sch or 'unknown'}")

    if interactive:
        confirm = _prompt(console, "\n  Submit all? (y/n)", "y")
        if confirm.lower() not in ("y", "yes"):
            console.print("[yellow]Cancelled.[/yellow]")
            return

    import subprocess

    # Check if scheduler is available
    dry_run = False
    first_scheduler = jobs[0][2] or "sbatch"
    if not shutil.which(first_scheduler):
        console.print(f"  [yellow]! {first_scheduler} not found — dry-run mode[/yellow]")
        console.print()
        dry_run = True

    submitted = 0
    for d, s, scheduler in jobs:
        cmd = scheduler or "sbatch"
        if dry_run:
            console.print(f"  [dim]→[/dim] cd {d.name} && {cmd} {s.name}")
            submitted += 1
        else:
            try:
                result = subprocess.run([cmd, str(s.name)], cwd=str(d),
                                        capture_output=True, text=True, timeout=10)
                jid = result.stdout.strip().split()[-1] if result.stdout.strip() else "OK"
                console.print(f"  [green]✓[/green] {d.name}/ → {jid}")
                submitted += 1
            except Exception as e:
                console.print(f"  [red]✗[/red] {d.name}/ → {e}")

    console.print()
    if dry_run:
        console.print(f"[green]✓ {submitted} jobs (dry-run).[/green]")
    else:
        console.print(f"[green]✓ {submitted}/{len(jobs)} submitted.[/green]")
    console.print()
