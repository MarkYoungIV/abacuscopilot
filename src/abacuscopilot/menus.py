"""Menu tree builder — legacy data model preserved for future use.

The interactive TUI (interactive.py) now renders menus directly without
these classes.  Menu / MenuItem / build_main_menu() are kept in case a
richer tree-based UI (e.g. rich Tree) is desired later.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from abacuscopilot.tasks import TaskRegistry


@dataclass
class MenuItem:
    """A single item in the interactive menu.

    Can represent either a task (has task_id) or a sub-menu entry (has sub_menu).
    """

    label: str
    task_id: int | None = None
    description: str = ""
    sub_menu: "Menu | None" = None


@dataclass
class Menu:
    """A menu level in the interactive hierarchy."""

    title: str
    items: list[MenuItem] = field(default_factory=list)
    sub_menus: dict[str, "Menu"] = field(default_factory=dict)


# =============================================================================
# Hard-coded welcome / main menu content
# =============================================================================

WELCOME_BANNER = """\
    ୧(๑•̀⌄•́๑)૭                               o(╥﹏╥)o
╭────────────────────────────────────────────────────╮
│                                                    │
│            A B A C U S - C O P I L O T             │
│    A pre- & post-processing copilot for ABACUS     │
│                                                    │
╰────────────────────────────────────────────────────╯"""


def get_version_info() -> str:
    """Get version string with date."""
    from abacuscopilot import __version__, __version_date__
    return f"AbacusCopilot v{__version__} ({__version_date__})"


def get_developer_info() -> str:
    """Get developer credit lines."""
    return "Developer: Xu Yang (xuyangmark@foxmail.com)\n           Rong-yu Zhang"


def build_main_menu() -> Menu:
    """Build the complete interactive menu tree from registered tasks.

    Tasks are grouped by their category into sub-menus.
    Categories are roughly ordered by the task numbering scheme
    (pre-processing first, then post-processing).

    Returns:
        Root Menu object with all task categories as sub-items.
    """
    registry = TaskRegistry()
    registry.discover_modules()

    main_menu = Menu("Main Menu")
    tasks_by_cat = registry.get_tasks_by_category()

    # Define display order for categories (pre-processing first)
    category_order = [
        "System",
        "INPUT",
        "STRU",
        "KPT",
        "Structure Editing",
        "Symmetry",
        "Batch",
        "SCF Analysis",
        "Band Structure",
        "DOS/PDOS",
        "Charge Density",
        "Mechanics",
        "Population",
        "Bond Order",
        "Work Function",
        "User Extensions",
    ]

    category_labels = {
        "System": "System & Configuration",
        "INPUT": "INPUT File Generation",
        "STRU": "STRU File Generation",
        "KPT": "K-Point Generation",
        "Structure Editing": "Structure Editor",
        "Symmetry": "Symmetry Analysis",
        "Batch": "Batch & Validation",
        "SCF Analysis": "SCF Convergence Analysis",
        "Band Structure": "Band Structure",
        "DOS/PDOS": "Density of States",
        "Charge Density": "Charge Density Analysis",
        "Mechanics": "Mechanical Properties",
        "Population": "Population Analysis",
        "Bond Order": "Bond Order",
        "Work Function": "Work Function",
        "User Extensions": "User Extensions",
    }

    # Display categories in a fixed order
    for cat_key in category_order:
        if cat_key in tasks_by_cat:
            tasks = tasks_by_cat[cat_key]
            sub = Menu(category_labels.get(cat_key, cat_key))

            for t in sorted(tasks, key=lambda x: x.task_id):
                label = f"[{t.task_id}] {t.name}"
                sub.items.append(
                    MenuItem(
                        label=label,
                        task_id=t.task_id,
                        description=t.description,
                    )
                )

            main_menu.items.append(
                MenuItem(
                    label=category_labels.get(cat_key, cat_key),
                    description=f"{len(tasks)} tasks",
                    sub_menu=sub,
                )
            )

    # Add any categories not in the predefined order
    for cat, tasks in sorted(tasks_by_cat.items()):
        if cat not in category_order:
            sub = Menu(cat)
            for t in sorted(tasks, key=lambda x: x.task_id):
                sub.items.append(
                    MenuItem(
                        label=f"[{t.task_id}] {t.name}",
                        task_id=t.task_id,
                        description=t.description,
                    )
                )
            main_menu.items.append(
                MenuItem(label=cat, sub_menu=sub)
            )

    # Add utility entries
    main_menu.items.append(
        MenuItem(
            label="[config] Edit Configuration",
            description="Edit ~/.abacuscopilot/config.yaml",
            task_id=2,
        )
    )

    return main_menu


def get_task_by_menu_selection(menu: Menu, selection: str) -> int | None:
    """Resolve a user input string to a task ID.

    Handles both direct task ID entry (e.g., "101") and menu navigation.
    Returns None if the selection doesn't match a task.
    """
    # Try direct numeric task ID
    try:
        tid = int(selection)
        registry = TaskRegistry()
        if registry.get_task(tid):
            return tid
    except ValueError:
        pass
    return None
