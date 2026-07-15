"""Shared console utilities for interactive prompts.

All task modules import from here instead of defining their own copies.
"""

from __future__ import annotations

from typing import Any

_console = None  # cached singleton — only one Console created per session


def _get_console():
    """Return the cached Rich Console singleton."""
    global _console
    if _console is None:
        from rich.console import Console
        _console = Console()
    return _console


def _prompt(console, question: str, default: Any = None) -> str:
    """Prompt the user for input with an optional default value."""
    if default is not None:
        result = console.input(f"  {question} [{default}]: ")
        return result.strip() if result.strip() else str(default)
    return console.input(f"  {question}: ").strip()


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
