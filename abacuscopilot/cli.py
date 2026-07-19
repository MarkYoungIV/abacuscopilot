"""Command-line interface for abacuscopilot.

Supports both interactive mode (abacuscopilot with no arguments) and
task-driven mode (abacuscopilot -task <id> with template defaults).

Usage:
    abacuscopilot                          # Launch interactive TUI
    abacuscopilot -task 101                # Run task 101 non-interactively (template defaults)
    abacuscopilot -task 721 --bands OUT.Si/BANDS_1.dat  # Task with args
    abacuscopilot --list-tasks             # List all registered tasks
    abacuscopilot --version                # Print version
"""

from __future__ import annotations

import argparse
import sys

from abacuscopilot import __version__
from abacuscopilot.tasks import TaskRegistry


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the main argument parser."""
    parser = argparse.ArgumentParser(
        prog="abacuscopilot",
        description=f"A pre- and post-processing toolkit for the ABACUS DFT software (v{__version__})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  abacuscopilot                          Launch interactive mode
  abacuscopilot -task 101                Generate SCF INPUT interactively
  abacuscopilot -task 721 --bands OUT/BANDS_1.dat
  abacuscopilot --list-tasks             Show all available tasks
        """,
    )

    parser.add_argument(
        "-task", "--task",
        type=int,
        metavar="ID",
        help="Task number to execute",
    )
    parser.add_argument(
        "--version", "-V",
        action="store_true",
        help="Print version and exit",
    )
    parser.add_argument(
        "--list-tasks", "-l",
        action="store_true",
        help="List all registered tasks",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        metavar="DIR",
        default=".",
        help="Output directory for generated files (default: current directory)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clean working directory (keep INPUT, KPT, STRU, *.{upf,orb}, sub*)",
    )
    parser.add_argument(
        "--state",
        action="store_true",
        help="Quick calculation status: done? converged? energy? (same as task 9908)",
    )
    parser.add_argument(
        "--md",
        action="store_true",
        help="Monitor a running MD calculation (auto-find OUT.ABACUS/MD_dump)",
    )

    return parser


def print_version():
    """Print abacuscopilot version."""
    from abacuscopilot import __version__
    print(f"AbacusCopilot v{__version__} (2026-07-14)")
    print("Developer: Xu Yang (xuyangmark@foxmail.com)")
    print("           Rong-yu Zhang")


def print_task_list():
    """Print all registered tasks in menu order with section headers."""
    from abacuscopilot.interactive import (
        _CATEGORY_LABELS,
        _DYNAMICS_CATEGORIES,
        _ELECTRONIC_CATEGORIES,
        _FIXED_MENU_NUMBERS,
        _MISC_CATEGORIES,
        _STRUCTURAL_CATEGORIES,
    )
    from abacuscopilot.tasks import TaskRegistry
    registry = TaskRegistry()
    registry.discover_modules()
    by_category = registry.get_tasks_by_category()

    sections = [
        ("Structural Utilities", _STRUCTURAL_CATEGORIES),
        ("Electronic Utilities", _ELECTRONIC_CATEGORIES),
        ("Dynamics Utilities", _DYNAMICS_CATEGORIES),
        ("Misc Utilities", _MISC_CATEGORIES),
    ]

    n = 1
    print()
    for sec_title, cat_keys in sections:
        print(f"  ==================== {sec_title} ====================")
        print()
        for ck in cat_keys:
            tasks = by_category.get(ck, [])
            label = _CATEGORY_LABELS.get(ck, ck)
            menu_num = _FIXED_MENU_NUMBERS.get(ck, n)
            if tasks:
                print(f"  {menu_num:>2})  {label}")
                if menu_num == n:
                    n += 1
                for t in tasks:
                    desc = f"  -- {t.description}" if t.description else ""
                    print(f"        {t.task_id:>4d})  {t.name}{desc}")
                print()
            elif ck in _CATEGORY_LABELS:
                # Show empty categories as coming-soon placeholders
                print(f"  {menu_num:>2})  {label} (coming soon)")
                if menu_num == n:
                    n += 1
                print()
    print()


def parse_task_args(argv: list[str]) -> tuple[int, list[str]]:
    """Parse task-specific arguments from the command line.

    All arguments after --task are passed to the task handler.
    This allows tasks to define their own argument handling.

    Args:
        argv: Full command-line arguments.

    Returns:
        Tuple of (task_id, remaining_args).
    """
    # Find --task or -task
    task_id = None
    remaining = []

    i = 0
    while i < len(argv):
        if argv[i] in ("--task", "-task"):
            if i + 1 < len(argv):
                try:
                    task_id = int(argv[i + 1])
                except ValueError:
                    print(f"Error: Invalid task ID: {argv[i + 1]}", file=sys.stderr)
                    sys.exit(1)
                i += 2
                continue
        remaining.append(argv[i])
        i += 1

    return task_id, remaining


def main(argv: list[str] | None = None):
    """Main entry point for abacuscopilot.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:]).
    """
    if argv is None:
        argv = sys.argv[1:]

    # Single parse_known_args() call: --task/<task-id> is a known argument;
    # everything else (e.g. --bands, --erange, positional file paths) lands
    # in `remaining` and is forwarded to the task handler as-is.
    parser = build_argument_parser()
    args, remaining = parser.parse_known_args(argv)

    # Handle special flags
    if args.version:
        print_version()
        return

    if args.list_tasks:
        print_task_list()
        return

    if args.clean:
        registry = TaskRegistry()
        registry.discover_modules()
        try:
            registry.dispatch(9904, args=remaining, interactive=True,
                            output_dir=args.output_dir)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        return

    if args.state:
        registry = TaskRegistry()
        registry.discover_modules()
        try:
            registry.dispatch(9908, args=remaining, interactive=True,
                            output_dir=args.output_dir)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        return

    if args.md:
        registry = TaskRegistry()
        registry.discover_modules()
        try:
            registry.dispatch(9905, args=remaining, interactive=True,
                            output_dir=args.output_dir)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        return

    # Dispatch task or enter interactive mode
    if args.task is not None:
        task_id = args.task
        registry = TaskRegistry()
        registry.discover_modules()
        try:
            registry.dispatch(task_id, args=remaining, interactive=False,
                            output_dir=args.output_dir)
        except Exception as e:
            print(f"Error executing task {task_id}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # No task specified, launch interactive mode
        from abacuscopilot.interactive import run_interactive
        run_interactive()


if __name__ == "__main__":
    main()
