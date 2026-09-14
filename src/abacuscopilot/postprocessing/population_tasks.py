"""Population analysis tasks for ABACUS output.

Task IDs 1301-1401

Mulliken and Hirshfeld population analysis, bond order analysis, and charge
decomposition from ABACUS LCAO output.
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


def parse_mulliken_from_file(filepath: str | Path) -> dict[str, Any] | None:
    """Parse ABACUS's ``mulliken.txt`` file (produced by ``out_mul 1``).

    ABACUS v3.10+ writes the Mulliken populations to a separate
    ``OUT.ABACUS/mulliken.txt`` file (not inline in the running log) with per
    atom blocks::

        0   Zeta of C   Spin 1
        ...
        Total Charge on atom:  C  4.0506

    Returns dict with 'charges' (atom_label → population) in atom order.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None
    charges: dict[str, float] = {}
    cur_idx: int | None = None
    for line in filepath.read_text(errors="ignore").splitlines():
        m = re.match(r"^\s*(\d+)\s+Zeta of\s+(\w+)", line)
        if m:
            cur_idx = int(m.group(1))
            continue
        m = re.match(r"^\s*Total Charge on atom:\s*(\w+)\s+(-?\d+\.?\d*)", line)
        if m and cur_idx is not None:
            charges[f"{m.group(1)}{cur_idx + 1}"] = float(m.group(2))
    if not charges:
        return None
    return {"charges": charges, "total_charge": sum(charges.values()),
            "method": "Mulliken"}


def parse_hirshfeld_from_log(filepath: str | Path) -> dict[str, Any] | None:
    """Extract Hirshfeld charges from ABACUS output.

    Requires the calculation to have run with ``out_hirshfeld 1`` in INPUT.
    ABACUS prints a "Hirshfeld charges" section in the running log (one charge
    per atom), e.g.::

        Hirshfeld charges of Atom 1 (Li)   =  0.876
        Hirshfeld charges of Atom 2 (Cl)   = -0.292
        ...

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

    patterns = [
        # Section table like Mulliken/Lowdin
        r"Hirshfeld\s+charges.*?\n(.*?)(?:\n\s*\n|\n\s*(?:Mulliken|Lowdin|END)|\Z)",
        r"HIRSHFELD\s+CHARGES.*?\n(.*?)(?:\n\s*\n|\n\s*(?:MULLIKEN|LOWDIN)|\Z)",
        # Per-atom lines: "Hirshfeld charges of Atom N (El) = value"
        r"(?:Hirshfeld|hirshfeld)\s+charges?\s+of\s+Atom\s+\d+\s*\((\w+)\)\s*=\s*(-?\d+\.?\d*)",
    ]

    # 1) per-atom lines first (more specific — the table regex would otherwise
    #    swallow "Hirshfeld charges of Atom N (El) = v" as a bogus table)
    found = re.findall(patterns[2], content, re.IGNORECASE)
    if found:
        charges = {el: float(v) for el, v in found}
        return {"charges": charges,
                "total_charge": sum(charges.values()),
                "method": "Hirshfeld"}

    # 2) table-style section
    for pattern in patterns[:2]:
        m = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
        if m:
            section = m.group(1).strip()
            charges = {}
            for line in section.split("\n"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        charge = float(parts[-1])
                        label = parts[0]
                        charges[label] = charge
                    except ValueError:
                        continue
            if charges:
                return {"charges": charges,
                        "total_charge": sum(charges.values()),
                        "method": "Hirshfeld"}

    return None


def _display_population_table(console, data: dict[str, Any]) -> None:
    """Display population analysis results as a rich table.

    For Mulliken/Lowdin the raw value is the electron POPULATION on the atom;
    the Net Charge column reports Z_val − population, i.e. how many electrons
    the atom gained (negative) or lost (positive) relative to neutral.  For
    Hirshfeld the raw value is already the (net) charge.
    """
    from rich.table import Table

    method = data.get("method", "")
    is_population = method.lower() in ("mulliken", "lowdin")
    is_hirshfeld = method.lower() == "hirshfeld"

    # Valence electrons per element (from the STRU's UPF when available).
    structure = None
    if Path("STRU").exists():
        try:
            from abacuscopilot.io.stru_file import read_stru
            structure = read_stru("STRU")
        except Exception:
            structure = None

    def _net(label: str, value: float) -> float:
        if is_hirshfeld:
            return value
        m = re.match(r"([A-Za-z]+)", label)
        el = m.group(1) if m else label
        if structure is not None:
            from abacuscopilot.postprocessing.charge_tasks import _zval_for
            zval = _zval_for(el, structure)
        else:
            from abacuscopilot.postprocessing.charge_tasks import _DEFAULT_ZVAL
            zval = float(_DEFAULT_ZVAL.get(el, 0.0))
        return zval - value  # positive = lost e⁻, negative = gained e⁻

    table = Table(title=f"{method} Population Analysis")
    table.add_column("Atom", style="cyan")
    table.add_column("Population (e)" if is_population else "Charge (e)", justify="right")
    table.add_column("Net Charge (e)", justify="right")

    charges = data["charges"]
    net_total = 0.0
    for label, value in charges.items():
        net = _net(label, value)
        net_total += net
        table.add_row(label, f"{value:.6f}", f"{net:+.6f}")

    table.add_section()
    table.add_row("[bold]Total[/bold]",
                  f"[bold]{data['total_charge']:.6f}[/bold]",
                  f"[bold]{net_total:+.6f}[/bold]")

    console.print()
    console.print(table)
    if is_population:
        console.print("  [dim]Net Charge (e) = valence − population: + means the atom lost electrons "
                      "(cation-like), − means it gained (anion-like).[/dim]")


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

    # Find the population file — ABACUS v3.10+ writes out_mul data to a
    # separate mulliken.txt; older versions printed it inline in the log.
    src = None
    if args:
        for arg in args:
            p = Path(arg)
            if p.exists():
                src = p
                break
    data = None
    if src is None:
        src = _find_population_file()
    if src is not None:
        data = parse_mulliken_from_file(src)
        console.print(f"  [dim]Reading: {src}[/dim]")
    if data is None:
        log_path = _find_running_log()
        if log_path:
            console.print(f"  [dim]Reading: {log_path}[/dim]")
            data = parse_mulliken_from_log(log_path)

    if data is None:
        console.print("[yellow]Mulliken population data not found.[/yellow]")
        console.print("[dim]ABACUS LCAO mode with out_mul 1 (Mulliken output) is required.[/dim]")
        console.print("[dim]Looked for OUT.ABACUS/mulliken.txt and the running log.[/dim]")
        return

    _display_population_table(console, data)
    console.print()


# =============================================================================
# Task 1401: Bond order analysis
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


@task(1401, category="Bond Order", name="Mulliken Bond Order",
      description="Mulliken bond order / overlap population from ABACUS LCAO output (out_mul 1)")
def task_bond_order(args: list[str] | None = None, interactive: bool = True) -> None:
    """Display Mulliken bond orders (overlap populations) from ABACUS LCAO output.

    Note: this is the Mulliken-type bond order (overlap population), derived
    from the off-diagonal (PS) terms ABACUS prints with ``out_mul 1``.  It is
    the simplest, most basis-dependent bond-order type — not the Mayer/Wiberg
    bond orders used in most modern analyses.
    """
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
