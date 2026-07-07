"""Population analysis tasks for ABACUS output.

Task IDs 761-769

Mulliken and Lowdin population analysis, bond order analysis,
and charge decomposition from ABACUS LCAO output.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import task

# =============================================================================
# Population data parsers
# =============================================================================


def _find_population_file() -> Path | None:
    """Find the Mulliken/population output file."""
    candidates = (
        list(Path().glob("mulliken*")) +
        list(Path().glob("MULLIKEN*")) +
        list(Path().glob("OUT.*/mulliken*")) +
        list(Path().glob("OUT.*/MULLIKEN*"))
    )
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _find_running_log() -> Path | None:
    """Find the ABACUS running log file."""
    candidates = (
        list(Path().glob("running_*.log")) +
        list(Path().glob("OUT.*/running_scf.log")) +
        list(Path().glob("OUT.*/running*.log"))
    )
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def parse_mulliken_from_log(filepath: str | Path) -> dict[str, Any] | None:
    """Extract Mulliken population analysis from ABACUS output.

    ABACUS (LCAO mode) prints Mulliken charges in the running log
    with a section like:

        Mulliken population analysis:
         Atom    Charge    ...
         Si1     3.8452    ...
         O1      6.1548    ...

    Args:
        filepath: Path to the running log file.

    Returns:
        Dict with 'charges' (dict[atom_label -> charge]), 'total_charge',
        or None if not found.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    content = filepath.read_text(errors="ignore")

    # Look for Mulliken population section
    mulliken_patterns = [
        r"Mulliken\s+population.*?\n(.*?)(?:\n\s*\n|\n\s*(?:Lowdin|Hirshfeld|END))",
        r"MULLIKEN\s+CHARGES.*?\n(.*?)(?:\n\s*\n|\n\s*(?:LOWDIN|HIRSHFELD))",
        r"(?:Atom|Element)\s+Charge.*?\n(.*?)(?:\n\s*\n|\Z)",
    ]

    for pattern in mulliken_patterns:
        m = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
        if m:
            section = m.group(1).strip()
            lines = section.split("\n")
            charges = {}
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    label = parts[0]
                    try:
                        charge = float(parts[1])
                        charges[label] = charge
                    except ValueError:
                        continue

            if charges:
                total = sum(charges.values())
                return {
                    "charges": charges,
                    "total_charge": total,
                    "method": "Mulliken",
                }

    # Also try to find total charge per atom from final SCF output
    charge_re = re.findall(
        r"(?:charge|Charge)\s+(?:of|on)\s+atom\s+(\w+)\s*[=:]\s*(-?\d+\.?\d*)",
        content, re.IGNORECASE
    )
    if charge_re:
        charges = {label: float(c) for label, c in charge_re}
        if charges:
            return {
                "charges": charges,
                "total_charge": sum(charges.values()),
                "method": "Mulliken",
            }

    return None


def parse_lowdin_from_log(filepath: str | Path) -> dict[str, Any] | None:
    """Extract Lowdin population analysis from ABACUS output.

    Similar format to Mulliken but with "Lowdin" header.

    Args:
        filepath: Path to the running log file.

    Returns:
        Dict with 'charges', 'total_charge', or None if not found.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    content = filepath.read_text(errors="ignore")

    lowdin_patterns = [
        r"L[oö]wdin\s+population.*?\n(.*?)(?:\n\s*\n|\n\s*(?:Mulliken|Hirshfeld|END))",
        r"LOWDIN\s+CHARGES.*?\n(.*?)(?:\n\s*\n|\n\s*(?:MULLIKEN|HIRSHFELD))",
    ]

    for pattern in lowdin_patterns:
        m = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
        if m:
            section = m.group(1).strip()
            lines = section.split("\n")
            charges = {}
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    label = parts[0]
                    try:
                        charge = float(parts[1])
                        charges[label] = charge
                    except ValueError:
                        continue

            if charges:
                total = sum(charges.values())
                return {
                    "charges": charges,
                    "total_charge": total,
                    "method": "Lowdin",
                }

    return None


def _display_population_table(console, data: dict[str, Any]) -> None:
    """Display population analysis results as a rich table."""
    from rich.table import Table

    table = Table(title=f"{data['method']} Population Analysis")
    table.add_column("Atom", style="cyan")
    table.add_column("Charge (e)", justify="right")
    table.add_column("Net Charge (e)", justify="right")

    charges = data["charges"]
    for label, charge in charges.items():
        # Net charge = valence - Mulliken charge (approximate based on label)
        # For now, just show the raw charge
        table.add_row(label, f"{charge:.6f}", "—")

    table.add_section()
    table.add_row(
        "[bold]Total[/bold]",
        f"[bold]{data['total_charge']:.6f}[/bold]",
        "",
    )

    console.print()
    console.print(table)


# =============================================================================
# Task 761: Mulliken population analysis
# =============================================================================


@task(1301, category="Population", name="Mulliken Analysis",
      description="Parse Mulliken population charges from ABACUS output")
def task_mulliken(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display Mulliken population charges from ABACUS LCAO output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Mulliken Population Analysis ===[/bold cyan]")
    console.print()

    # Find population file
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
        console.print("[dim]Run an ABACUS LCAO SCF calculation first.[/dim]")
        return

    console.print(f"  [dim]Reading: {log_path}[/dim]")

    data = parse_mulliken_from_log(log_path)

    if data is None:
        console.print("[yellow]Mulliken population data not found in the log.[/yellow]")
        console.print("[dim]ABACUS LCAO mode with Mulliken output enabled is required.[/dim]")
        console.print("[dim]Check that your ABACUS version supports Mulliken analysis.[/dim]")
        return

    _display_population_table(console, data)
    console.print()


# =============================================================================
# Task 762: Lowdin population analysis
# =============================================================================


@task(1302, category="Population", name="Lowdin Analysis",
      description="Parse Lowdin population charges from ABACUS output")
def task_lowdin(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display Lowdin population charges from ABACUS LCAO output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Lowdin Population Analysis ===[/bold cyan]")
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
        return

    console.print(f"  [dim]Reading: {log_path}[/dim]")

    # Try Lowdin first, fall back to Mulliken
    data = parse_lowdin_from_log(log_path)
    if data is None:
        data = parse_mulliken_from_log(log_path)
        if data:
            console.print("[yellow]Lowdin data not found. Displaying Mulliken instead.[/yellow]")

    if data is None:
        console.print("[yellow]No population analysis data found in the log.[/yellow]")
        return

    _display_population_table(console, data)
    console.print()


# =============================================================================
# Task 763: Bond order analysis
# =============================================================================


def _parse_bond_orders(filepath: str | Path) -> dict[str, Any] | None:
    """Parse bond order / bond population data from ABACUS output.

    ABACUS may output bond populations in various formats.
    This function tries multiple patterns.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    content = filepath.read_text(errors="ignore")

    # Look for bond populations
    # Pattern: "Bond population between atom X and atom Y: value"
    bond_pat = re.compile(
        r"(?:bond|Bond)\s+(?:population|order)\s+between\s+atom\s+(\w+)\s+and\s+atom\s+(\w+)\s*[=:]\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )
    bonds = bond_pat.findall(content)

    if bonds:
        return {
            "bonds": [
                {"atom1": a1, "atom2": a2, "population": float(pop)}
                for a1, a2, pop in bonds
            ],
        }

    # Try simpler format: table with columns
    # "Overlap population matrix" section
    ovlp_pat = re.compile(
        r"(?:Overlap\s+population|Bond\s+population).*?\n(.*?)(?:\n\s*\n|\Z)",
        re.DOTALL | re.IGNORECASE,
    )
    m = ovlp_pat.search(content)
    if m:
        # Parse matrix rows
        section = m.group(1)
        matrix_data = []
        for line in section.strip().split("\n"):
            parts = line.split()
            nums = []
            for p in parts:
                try:
                    nums.append(float(p))
                except ValueError:
                    pass
            if nums:
                matrix_data.append(nums)

        if matrix_data:
            # Convert to bond list (upper triangle)
            bonds = []
            n = len(matrix_data)
            for i in range(n):
                for j in range(i + 1, min(len(matrix_data[i]), n)):
                    if abs(matrix_data[i][j]) > 0.001:
                        bonds.append({
                            "atom1": str(i + 1),
                            "atom2": str(j + 1),
                            "population": matrix_data[i][j],
                        })
            if bonds:
                return {"bonds": bonds}

    return None


@task(1303, category="Population", name="Bond Order",
      description="Extract bond populations / bond orders from ABACUS output")
def task_bond_order(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display bond populations from ABACUS LCAO output."""
    console = _get_console()

    console.print()
    console.print("[bold cyan]=== Bond Order Analysis ===[/bold cyan]")
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
        return

    console.print(f"  [dim]Reading: {log_path}[/dim]")

    data = _parse_bond_orders(log_path)

    if data is None or not data.get("bonds"):
        console.print("[yellow]No bond population data found in the log.[/yellow]")
        console.print("[dim]Bond population analysis is available in ABACUS LCAO mode.[/dim]")
        return

    from rich.table import Table

    table = Table(title="Bond Populations")
    table.add_column("Atom 1", style="cyan")
    table.add_column("Atom 2", style="cyan")
    table.add_column("Population", justify="right")
    table.add_column("Type", justify="center")

    for bond in data["bonds"]:
        pop = bond["population"]
        bond_type = (
            "[green]bonding[/green]" if pop > 0.1
            else "[yellow]weak[/yellow]" if abs(pop) < 0.1
            else "[red]antibonding[/red]"
        )
        table.add_row(
            bond["atom1"],
            bond["atom2"],
            f"{pop:.6f}",
            bond_type,
        )

    console.print()
    console.print(table)
    console.print()
