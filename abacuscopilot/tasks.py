"""Task registry and decorator for abacuscopilot.

Tasks are registered with a numeric ID and metadata via the @task decorator.
The central TaskRegistry singleton handles dispatch, discovery, and listing.

Usage:
    from abacuscopilot.tasks import task

    @task(101, category="INPUT", name="SCF INPUT",
          description="Generate INPUT file for SCF calculation")
    def generate_scf_input(args=None, interactive=True):
        ...
"""

from __future__ import annotations

import importlib
import importlib.util
import pkgutil
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable


@dataclass
class TaskMetadata:
    """Metadata for a registered task."""

    task_id: int
    category: str
    name: str
    handler: Callable
    description: str = ""
    requires_files: list[str] = field(default_factory=list)
    cli_args: list[dict] = field(default_factory=list)


class TaskRegistry:
    """Singleton registry mapping task IDs to handler metadata.

    Supports:
    - @task decorator for registration
    - dispatch by task_id
    - listing by category
    - auto-discovery of task modules
    """

    _instance: "TaskRegistry | None" = None
    _tasks: dict[int, TaskMetadata]
    _discovered: bool

    def __new__(cls) -> "TaskRegistry":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._tasks = {}
            cls._instance._discovered = False
        return cls._instance

    def register(
        self,
        task_id: int,
        category: str,
        name: str,
        *,
        description: str = "",
        requires_files: list[str] | None = None,
        cli_args: list[dict] | None = None,
    ) -> Callable:
        """Decorator that registers a function as a task handler.

        Args:
            task_id: Unique numeric task identifier.
            category: Task category for menu grouping.
            name: Short display name.
            description: Longer description shown in help.
            requires_files: Files that must exist in the working directory.
            cli_args: Additional CLI argument specs for this task.

        Returns:
            Decorator function.

        Raises:
            ValueError: If task_id is already registered.
        """
        if task_id in self._tasks:
            raise ValueError(f"Task {task_id} already registered: {self._tasks[task_id].name}")

        def decorator(func: Callable) -> Callable:
            self._tasks[task_id] = TaskMetadata(
                task_id=task_id,
                category=category,
                name=name,
                handler=func,
                description=description,
                requires_files=requires_files or [],
                cli_args=cli_args or [],
            )
            return func

        return decorator

    def dispatch(
        self, task_id: int, args: list[str] | None = None, interactive: bool = True,
        output_dir: str = ".",
    ) -> Any:
        """Execute a task by its numeric ID.

        Args:
            task_id: The task number to run.
            args: Positional arguments for CLI mode.
            interactive: Whether to run in interactive mode.
            output_dir: Output directory for generated files.

        Returns:
            Return value of the task handler.

        Raises:
            TaskNotFoundError: If task_id is not registered.
        """
        import argparse as _argparse

        from abacuscopilot.core.exceptions import TaskNotFoundError

        if task_id not in self._tasks:
            self.discover_modules()
            if task_id not in self._tasks:
                raise TaskNotFoundError(task_id)

        metadata = self._tasks[task_id]

        # Build a per-task argument parser when running non-interactively
        # and the task has declared cli_args metadata.
        parsed_args = None
        if not interactive and metadata.cli_args and args:
            parser = _argparse.ArgumentParser(
                prog=f"abacuscopilot -task {task_id}",
                description=metadata.name,
                add_help=True,
            )
            for spec in metadata.cli_args:
                name = spec["name"]
                kwargs = {k: v for k, v in spec.items() if k != "name"}
                parser.add_argument(name, **kwargs)
            parsed_args, unknown = parser.parse_known_args(args)
            if unknown:
                import sys as _sys
                print(f"Warning: unrecognized arguments: {' '.join(unknown)}",
                      file=_sys.stderr)
                print(f"  Run 'abacuscopilot -task {task_id} --help' for supported options.",
                      file=_sys.stderr)

        # Only pass parsed_args / output_dir if the handler accepts them (backward compat)
        import inspect as _inspect
        sig = _inspect.signature(metadata.handler)
        handler_kwargs = {"args": args, "interactive": interactive}
        has_var_kwargs = any(
            p.kind == _inspect.Parameter.VAR_KEYWORD
            for p in sig.parameters.values()
        )
        if "parsed_args" in sig.parameters or has_var_kwargs:
            handler_kwargs["parsed_args"] = parsed_args
        if "output_dir" in sig.parameters or has_var_kwargs:
            handler_kwargs["output_dir"] = output_dir
        return metadata.handler(**handler_kwargs)

    def get_task(self, task_id: int) -> TaskMetadata | None:
        """Get task metadata by ID."""
        return self._tasks.get(task_id)

    def get_categories(self) -> list[str]:
        """Return sorted list of unique categories."""
        return sorted(set(t.category for t in self._tasks.values()))

    def list_tasks(self, category: str | None = None) -> list[TaskMetadata]:
        """List tasks, optionally filtered by category.

        Returns tasks sorted by task_id.
        """
        tasks = self._tasks.values()
        if category:
            tasks = [t for t in tasks if t.category == category]
        return sorted(tasks, key=lambda t: t.task_id)

    def get_tasks_by_category(self) -> dict[str, list[TaskMetadata]]:
        """Group all tasks by category."""
        groups = defaultdict(list)
        for t in self._tasks.values():
            groups[t.category].append(t)
        return {cat: sorted(ts, key=lambda t: t.task_id) for cat, ts in groups.items()}

    def discover_modules(self):
        """Auto-discover task modules by scanning subpackages.

        This imports all modules under abacuscopilot.preprocessing and
        abacuscopilot.postprocessing to trigger @task registrations.

        In PyInstaller-frozen mode, pkgutil.walk_packages cannot walk
        the virtual filesystem, so we use an explicit module list instead.
        """
        if self._discovered:
            return

        import sys as _sys

        if getattr(_sys, "frozen", False):
            # PyInstaller frozen: explicit imports (walk_packages won't work)
            _frozen_modules = [
                "abacuscopilot.preprocessing.system_tasks",
                "abacuscopilot.preprocessing.input_tasks",
                "abacuscopilot.preprocessing.kpt_tasks",
                "abacuscopilot.preprocessing.stru_tasks",
                "abacuscopilot.preprocessing.structure_editing_tasks",
                "abacuscopilot.preprocessing.symmetry_tasks",
                "abacuscopilot.preprocessing.batch_tasks",
                "abacuscopilot.postprocessing.scf_tasks",
                "abacuscopilot.postprocessing.band_tasks",
                "abacuscopilot.postprocessing.dos_tasks",
                "abacuscopilot.postprocessing.charge_tasks",
                "abacuscopilot.postprocessing.workfunc_tasks",
                "abacuscopilot.postprocessing.mechanics_tasks",
                "abacuscopilot.postprocessing.population_tasks",
                "abacuscopilot.postprocessing.md_tasks",
            ]
            for mod_name in _frozen_modules:
                try:
                    importlib.import_module(mod_name)
                except ImportError:
                    pass
        else:
            packages = [
                "abacuscopilot.preprocessing",
                "abacuscopilot.postprocessing",
            ]

            for package_name in packages:
                try:
                    package = importlib.import_module(package_name)
                    for _, mod_name, _ in pkgutil.walk_packages(
                        package.__path__, prefix=package_name + "."
                    ):
                        try:
                            importlib.import_module(mod_name)
                        except ImportError:
                            pass  # Skip modules with missing optional dependencies
                except ImportError:
                    pass  # Package may not exist yet

        # Discover user extensions (~/.abacuscopilot/extensions/*.py)
        self._discover_user_extensions()

        self._discovered = True

    def _discover_user_extensions(self):
        """Scan ~/.abacuscopilot/extensions/ for user-defined task modules.

        Each .py file in the extensions directory is imported.  Tasks
        registered with IDs >= 9000 automatically appear under the
        "User Extensions" category in the interactive menu.
        """
        from pathlib import Path as _Path
        ext_dir = _Path.home() / ".abacuscopilot" / "extensions"
        if not ext_dir.is_dir():
            return

        for py_file in sorted(ext_dir.glob("*.py")):
            if py_file.name.startswith("_"):
                continue
            try:
                spec = importlib.util.spec_from_file_location(
                    f"abacuscopilot_user_ext.{py_file.stem}",
                    str(py_file),
                )
                if spec and spec.loader:
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
            except Exception:
                pass  # Skip broken user modules


# Convenience function for the decorator pattern
def task(
    task_id: int,
    category: str,
    name: str,
    *,
    description: str = "",
    requires_files: list[str] | None = None,
    cli_args: list[dict] | None = None,
) -> Callable:
    """Register a function as a task handler.

    Usage:
        @task(101, category="INPUT", name="SCF INPUT",
              description="Generate INPUT file for SCF calculation")
        def generate_scf_input(args=None, interactive=True):
            ...

    Args:
        task_id: Unique numeric task identifier.
        category: Task category for menu grouping (e.g., "INPUT", "Band Structure").
        name: Short display name.
        description: Longer description for help output.
        requires_files: Files that should exist in the working directory.
        cli_args: Additional argparse argument specs.

    Returns:
        Decorator that registers the function.
    """
    return TaskRegistry().register(
        task_id,
        category,
        name,
        description=description,
        requires_files=requires_files,
        cli_args=cli_args,
    )
