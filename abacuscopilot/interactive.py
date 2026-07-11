"""Interactive TUI for abacuscopilot — VASPKIT-style hierarchical menus.

Provides a menu-driven terminal interface similar to VASPKIT's
interactive mode. Users navigate numbered two-level menus to
select and execute tasks.

Default view is the main menu showing task categories grouped
into sections. Type a category number to enter its sub-menu.
Type '0' to quit, '9' to go back, or a direct task ID to run.

From the main menu, numbers 1..N always navigate to categories.
To run a task directly, type its numeric ID — it works from any
menu level (including the main menu, for IDs above N).

KeyboardInterrupt (Ctrl‑C) stays in the loop; EOF / 0 / quit exits.
"""

from __future__ import annotations

from dataclasses import dataclass

from abacuscopilot.console_utils import _get_console
from abacuscopilot.tasks import TaskRegistry

# ---------------------------------------------------------------------------
# Category → section mapping
# ---------------------------------------------------------------------------

# Human-readable labels shown in the main menu for each category key.
_CATEGORY_LABELS: dict[str, str] = {
    "INPUT":               "INPUT File Generation",
    "STRU":                "STRU File Generation",
    "KPT":                 "K-Point Generation",
    "Structure Editing":   "Structure Editor",
    "Symmetry":            "Symmetry Analysis",
    "SCF Analysis":        "SCF Convergence Analysis",
    "Band Structure":      "Band Structure Visualization",
    "DOS/PDOS":            "Density of States (DOS/PDOS)",
    "Charge Density":      "Charge Density Analysis",
    "Work Function":       "Work Function Analysis",
    "Mechanics":           "Mechanical Properties",
    "System":              "System & Configuration",
    "Batch":               "Batch Job Submission",
    "Population":          "Population Analysis",
    "MD Analysis":         "MD Trajectory Analysis",
    "Lattice Dynamics":    "Lattice Dynamics",
    "Reaction Dynamics":   "Reaction Dynamics",
    "User Extensions":     "User Extensions",
}

# Categories grouped into VASPKIT-style sections (display order).
_STRUCTURAL_CATEGORIES = [
    "INPUT", "STRU", "KPT", "Structure Editing", "Symmetry", "Batch",
]

_ELECTRONIC_CATEGORIES = [
    "SCF Analysis", "Band Structure", "DOS/PDOS",
    "Charge Density", "Work Function", "Mechanics", "Population",
]

_DYNAMICS_CATEGORIES = [
    "MD Analysis", "Lattice Dynamics", "Reaction Dynamics",
]

_MISC_CATEGORIES = [
    "System", "User Extensions",
]

# Categories whose menu numbers are fixed (not sequential)
_FIXED_MENU_NUMBERS: dict[str, int] = {
    "System": 99,
}

_ALL_KNOWN_CATEGORIES = (
    _STRUCTURAL_CATEGORIES + _ELECTRONIC_CATEGORIES +
    _DYNAMICS_CATEGORIES + _MISC_CATEGORIES
)


@dataclass
class Section:
    """A group of related categories shown under a section divider."""
    title: str
    categories: list[tuple[str, str, list]]  # [(key, display_label, tasks)]


def _build_sections(registry: TaskRegistry) -> list[Section]:
    """Build the section list from registered tasks.

    Only categories that actually have tasks are included.
    Each category appears in exactly one section — predefined
    categories go to their assigned section; leftovers land in
    Misc Utilities.
    """
    by_cat = registry.get_tasks_by_category()
    assigned: set[str] = set()

    def _pick(cat_keys: list[str]) -> list[tuple[str, str, list]]:
        """Collect (key, label, tasks) for each category in *cat_keys*.
        Empty categories get a placeholder label.
        """
        result: list[tuple[str, str, list]] = []
        for ck in cat_keys:
            if ck in assigned:
                continue
            label = _CATEGORY_LABELS.get(ck, ck)
            tasks = sorted(by_cat.get(ck, []), key=lambda t: t.task_id)
            if not tasks and ck not in _ALL_KNOWN_CATEGORIES:
                continue  # skip truly unknown categories
            if not tasks:
                label = f"{label} (coming soon)"
            result.append((ck, label, tasks))
            assigned.add(ck)
        return result

    sections: list[Section] = []

    structural = _pick(_STRUCTURAL_CATEGORIES)
    if structural:
        sections.append(Section(title="Structural Utilities", categories=structural))

    electronic = _pick(_ELECTRONIC_CATEGORIES)
    if electronic:
        sections.append(Section(title="Electronic Utilities", categories=electronic))

    dynamics = _pick(_DYNAMICS_CATEGORIES)
    if dynamics:
        sections.append(Section(title="Dynamics Utilities", categories=dynamics))

    # Misc: predefined list + any remaining unassigned categories
    misc = _pick(_MISC_CATEGORIES)
    leftovers = _pick(sorted(by_cat.keys()))
    misc.extend(leftovers)
    if misc:
        sections.append(Section(title="Misc Utilities", categories=misc))

    return sections


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------

def _show_welcome(console) -> None:
    """Print the welcome banner and version line."""
    from abacuscopilot.menus import WELCOME_BANNER, get_developer_info, get_version_info
    console.print(f"[bold cyan]{WELCOME_BANNER}[/bold cyan]")
    console.print(f"[dim]{get_version_info()}[/dim]")
    console.print(f"[dim]{get_developer_info()}[/dim]")


def _section_header(title: str, width: int = 60) -> str:
    """Build a VASPKIT-style section divider.

    Example:  ==================== Structural Utilities ====================
    """
    inner = f" {title} "
    side = (width - len(inner)) // 2
    return f"[bold cyan]{'=' * side}{inner}{'=' * (width - len(inner) - side)}[/bold cyan]"


def _show_main_menu(console, sections: list[Section]) -> None:
    """Render the main menu with VASPKIT-style section dividers.

    Categories are numbered sequentially (1‑based) across all
    sections. Within each section they are laid out in two columns.
    """
    console.print()

    # Build flat index: sequential number → category key
    # Categories in _FIXED_MENU_NUMBERS use their hard-coded number
    flat: list[tuple[int, str, int]] = []  # [(index, label, n_tasks)]
    for section in sections:
        for _key, label, tasks in section.categories:
            idx = _FIXED_MENU_NUMBERS.get(_key, len(flat) + 1)
            flat.append((idx, label, len(tasks)))

    for section in sections:
        console.print(" " + _section_header(section.title))
        console.print()

        # Build the list of (index, label) pairs for this section
        items: list[tuple[int, str]] = []
        for _key, label, _tasks in section.categories:
            for idx, lbl, _ in flat:
                if lbl == label:
                    items.append((idx, label))
                    break

        col_width = 36
        for i in range(0, len(items), 2):
            left_idx, left_label = items[i]
            left_str = f" {left_idx:>2})  {left_label}"

            if i + 1 < len(items):
                right_idx, right_label = items[i + 1]
                right_str = f" {right_idx:>2})  {right_label}"
                console.print(left_str.ljust(col_width) + right_str)
            else:
                console.print(left_str)

        console.print()

    console.print(" [dim]0)[/dim]  Quit")
    console.print()


def _show_category_submenu(console, category_key: str, tasks: list) -> None:
    """Render the sub-menu for a specific category.

    Tasks are listed with their real (3‑digit) task IDs.  Direct task ID
    entry works from any menu level.
    """
    console.print()

    label = _CATEGORY_LABELS.get(category_key, category_key)
    console.print(" " + _section_header(label))
    console.print()

    for t in tasks:
        desc = f" -- {t.description}" if t.description else ""
        console.print(f" [green]{t.task_id:>4d})[/green]  {t.name}{desc}")

    console.print()
    console.print(" [dim]0)[/dim]  Quit")
    console.print(" [dim]9)[/dim]  Back")
    console.print()


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------

def _show_help(console) -> None:
    """Display help information."""
    console.print()
    console.print("[bold]AbacusCopilot Help[/bold]")
    console.print()
    console.print("  [green]abacuscopilot[/green]           Launch interactive mode (this)")
    console.print("  [green]abacuscopilot -task N[/green]   Run specific task N")
    console.print()
    console.print("[bold]Interactive commands:[/bold]")
    console.print("  [cyan]<number>[/cyan]    Select category (main menu) or run task (sub-menu)")
    console.print("  [cyan]0[/cyan]           Quit from anywhere")
    console.print("  [cyan]9[/cyan]           Back to main menu")
    console.print("  [cyan]menu[/cyan]        Jump to main menu")
    console.print("  [cyan]help[/cyan]        Show this help")
    console.print("  [cyan]version[/cyan]     Show version")
    console.print()
    console.print("[bold]Direct task entry:[/bold]")
    console.print("  Type a task ID (e.g. [green]301[/green], [green]711[/green]) from any menu level to run it directly.")
    console.print()


# ---------------------------------------------------------------------------
# Task dispatch helper
# ---------------------------------------------------------------------------

def _dispatch_task(registry: TaskRegistry, task_id: int, console) -> bool:
    """Run a task by ID.  Returns True if the task completed normally."""
    try:
        registry.dispatch(task_id, args=[], interactive=True)
        console.print()
        return True
    except Exception as e:
        msg = str(e).replace("[", "\\[").replace("]", "\\]")
        console.print(f"[red]Error running task {task_id}: {msg}[/red]")
        return False


# ---------------------------------------------------------------------------
# Main interactive loop
# ---------------------------------------------------------------------------

def run_interactive():
    """Launch the interactive TUI loop — VASPKIT-style hierarchical menus.

    Flow
    ----
    1. Welcome banner + version
    2. Main menu (categories grouped in sections, numbered 1…N)
    3. User picks a category → sub-menu with tasks
    4. From sub-menu: task ID runs a task, 0 quits, 9 goes back
    5. Direct task IDs work from any level
    """
    console = _get_console()

    _show_welcome(console)

    registry = TaskRegistry()
    registry.discover_modules()

    sections = _build_sections(registry)

    # --- Build lookup tables ---
    # Sequential index (1‑based) → category key
    cat_index: dict[int, str] = {}
    for section in sections:
        for cat_key, _label, _tasks in section.categories:
            idx = _FIXED_MENU_NUMBERS.get(cat_key, len(cat_index) + 1)
            cat_index[idx] = cat_key

    # Category key → sorted task list
    cat_tasks: dict[str, list] = {}
    for section in sections:
        for cat_key, _label, tasks in section.categories:
            cat_tasks[cat_key] = tasks

    # Set of all known task IDs (for direct dispatch)
    all_task_ids: set[int] = set()
    for tasks in cat_tasks.values():
        for t in tasks:
            all_task_ids.add(t.task_id)

    n_categories = len(cat_index)
    current: str = "main"  # "main" | category_key

    while True:
        # --- Render current view ---
        if current == "main":
            _show_main_menu(console, sections)
            prompt = "[bold green]abacuscopilot>[/bold green] "
        else:
            tasks = cat_tasks.get(current, [])
            _show_category_submenu(console, current, tasks)
            prompt = f"[bold green]abacuscopilot[/bold green][green]/[/green][bold yellow]{current}[/bold yellow][bold green]>[/bold green] "

        # --- Read input ---
        try:
            user_input = console.input(prompt).strip()
        except KeyboardInterrupt:
            console.print("\n[yellow]Interrupted. Type '0' to quit, '9' for main menu.[/yellow]")
            continue
        except EOFError:
            from abacuscopilot.quotes import farewell
            farewell(console)
            break

        if not user_input:
            continue

        cmd = user_input.lower()

        # --- Global text commands ---
        if cmd in ("quit", "exit", "q"):
            from abacuscopilot.quotes import farewell
            farewell(console)
            break

        if cmd == "help":
            _show_help(console)
            continue

        if cmd == "version":
            from abacuscopilot import __version__
            console.print(f"AbacusCopilot v{__version__}")
            continue

        # --- Navigation ---
        if cmd == "0":
            from abacuscopilot.quotes import farewell
            farewell(console)
            break

        if cmd == "9":
            if current != "main":
                current = "main"
                continue
            # 9 from main — let it fall through to numeric processing below

        if cmd in ("menu", "m"):
            current = "main"
            continue

        # --- Numeric input ---
        try:
            num = int(user_input)
        except ValueError:
            console.print(
                f"[red]Unknown command: '{user_input}'. "
                f"Type 'help' for assistance.[/red]"
            )
            continue

        if current == "main":
            # Main menu: 1..N picks a category, anything higher is a direct task ID
            if num in cat_index:
                current = cat_index[num]
            elif num in all_task_ids:
                if _dispatch_task(registry, num, console):
                    from abacuscopilot.quotes import farewell
                    farewell(console)
                    break
            else:
                console.print(
                    f"[red]Invalid selection: {num}. "
                    f"Enter a menu number or task ID.[/red]"
                )
        else:
            # Sub-menu: any number is treated as a task ID
            if num in all_task_ids:
                if _dispatch_task(registry, num, console):
                    from abacuscopilot.quotes import farewell
                    farewell(console)
                    break
            else:
                console.print(
                    f"[red]Unknown task ID: {num}. "
                    f"Type '9' to go back to the main menu.[/red]"
                )
