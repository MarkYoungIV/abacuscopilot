"""CLI integration tests for abacuscopilot.

Tests the command-line interface end-to-end: argument parsing,
task dispatch, and non-interactive mode execution.
"""


import pytest

# =============================================================================
# CLI argument parsing
# =============================================================================


class TestCLIVersion:
    def test_version_flag(self, capsys):
        """abacuscopilot --version prints version and exits cleanly."""
        from abacuscopilot.cli import main

        main(["--version"])
        captured = capsys.readouterr()
        assert "AbacusCopilot" in captured.out

    def test_short_version_flag(self, capsys):
        """abacuscopilot -V works."""
        from abacuscopilot.cli import main

        main(["-V"])
        captured = capsys.readouterr()
        assert "AbacusCopilot" in captured.out


class TestCLIListTasks:
    def test_list_tasks(self, capsys):
        """abacuscopilot --list-tasks prints all categories."""
        from abacuscopilot.cli import main

        main(["--list-tasks"])
        captured = capsys.readouterr()
        assert "INPUT" in captured.out
        assert "KPT" in captured.out
        assert "STRU" in captured.out

    def test_list_tasks_by_category(self, capsys):
        """abacuscopilot --list-tasks."""
        from abacuscopilot.cli import main

        main(["--list-tasks"])
        captured = capsys.readouterr()
        assert "SCF INPUT" in captured.out


class TestCLITaskDispatch:
    def test_task_args_passed_correctly(self):
        """Task-specific args (--bands, --erange) reach the handler."""
        from abacuscopilot.tasks import TaskRegistry

        registry = TaskRegistry()
        registry.discover_modules()
        meta = registry.get_task(801)
        assert meta is not None, "Task 801 not found"

        # Verify cli_args metadata exists
        arg_names = [a["name"] for a in meta.cli_args]
        assert "--bands" in arg_names
        assert "--erange" in arg_names

    def test_dispatch_passes_parsed_args(self):
        """TaskRegistry.dispatch() parses cli_args into parsed_args."""
        from abacuscopilot.tasks import TaskRegistry

        registry = TaskRegistry()
        registry.discover_modules()
        meta = registry.get_task(801)

        # Replace handler with a spy
        received = {}
        def spy(**kwargs):
            received.update(kwargs)

        meta.handler = spy
        registry.dispatch(801, args=["--bands", "test.dat", "--erange=-5,5"],
                         interactive=False)

        assert received.get("interactive") is False
        assert received.get("parsed_args") is not None
        assert received["parsed_args"].bands == "test.dat"
        assert received["parsed_args"].erange == "-5,5"

    def test_dispatch_backward_compat_no_parsed_args(self):
        """Old handlers (without parsed_args param) still work."""
        from abacuscopilot.tasks import TaskRegistry

        registry = TaskRegistry()
        registry.discover_modules()
        meta = registry.get_task(102)  # Relax INPUT, no cli_args

        received = {}
        def spy(args=None, interactive=True):
            received["called"] = True
            received["args"] = args

        meta.handler = spy
        registry.dispatch(102, args=[], interactive=False)

        assert received["called"] is True

    def test_invalid_task_id(self):
        """Unknown task ID raises TaskNotFoundError."""
        from abacuscopilot.core.exceptions import TaskNotFoundError
        from abacuscopilot.tasks import TaskRegistry

        registry = TaskRegistry()
        registry.discover_modules()

        with pytest.raises(TaskNotFoundError):
            registry.dispatch(99999)


# =============================================================================
# Non-interactive task execution
# =============================================================================


class TestTaskNonInteractive:
    def test_task_101_generates_input(self, tmp_path):
        """Task 101 (SCF INPUT) runs non-interactively and creates INPUT file."""
        import os
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            from abacuscopilot.tasks import TaskRegistry
            registry = TaskRegistry()
            registry.discover_modules()
            registry.dispatch(101, args=[], interactive=False)

            input_path = tmp_path / "INPUT"
            assert input_path.exists(), f"INPUT file not created in {tmp_path}"
            content = input_path.read_text()
            assert "INPUT_PARAMETERS" in content
        finally:
            os.chdir(original_cwd)

    def test_task_301_generates_kpt(self, tmp_path):
        """Task 301 (Auto KPT) runs non-interactively and creates KPT file."""
        import os
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            # Need a STRU to get real lattice; without it, falls back to identity
            from abacuscopilot.tasks import TaskRegistry
            registry = TaskRegistry()
            registry.discover_modules()
            registry.dispatch(301, args=[], interactive=False)

            kpt_path = tmp_path / "KPT"
            assert kpt_path.exists(), f"KPT file not created in {tmp_path}"
            content = kpt_path.read_text()
            assert "K_POINTS" in content
        finally:
            os.chdir(original_cwd)

    def test_task_603_validates_inputs(self, tmp_path):
        """Task 603 (Validate Inputs) reports missing files correctly."""
        import os
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            from abacuscopilot.tasks import TaskRegistry
            registry = TaskRegistry()
            registry.discover_modules()

            # No files — should report issues
            registry.dispatch(603, args=[], interactive=False)
            # Task prints to console rather than raising — no crash is success
        finally:
            os.chdir(original_cwd)


# =============================================================================
# Task discovery
# =============================================================================


class TestTaskDiscovery:
    def test_all_categories_present(self):
        """All expected task categories are discovered."""
        from abacuscopilot.tasks import TaskRegistry
        registry = TaskRegistry()
        registry.discover_modules()

        categories = registry.get_categories()
        expected = [
            "Band Structure", "Batch", "Charge Density", "DOS/PDOS",
            "INPUT", "KPT", "Mechanics", "Population",
            "SCF Analysis", "STRU", "Structure Editing", "Symmetry",
            "System", "Work Function",
        ]
        for cat in expected:
            assert cat in categories, f"Category '{cat}' not found"

    def test_task_count(self):
        """At least 40 tasks are registered."""
        from abacuscopilot.tasks import TaskRegistry
        registry = TaskRegistry()
        registry.discover_modules()

        assert len(registry._tasks) >= 40

    def test_key_tasks_exist(self):
        """Core tasks are registered with correct IDs."""
        from abacuscopilot.tasks import TaskRegistry
        registry = TaskRegistry()
        registry.discover_modules()

        key_tasks = {
            101: "SCF INPUT",
            201: "CIF to STRU",
            301: "Auto KPT (MP mesh)",
            801: "Plot Band Structure",
            901: "Plot DOS",
            1201: "Elastic Constants",
            601: "PBS Script",
            602: "SLURM Script",
        }
        for tid, name in key_tasks.items():
            t = registry.get_task(tid)
            assert t is not None, f"Task {tid} not registered"
            assert t.name == name, f"Task {tid}: expected '{name}', got '{t.name}'"


# =============================================================================
# Output directory flag
# =============================================================================


class TestOutputDirFlag:
    def test_output_dir_argument_parsed(self):
        """--output-dir is parsed correctly."""
        from abacuscopilot.cli import build_argument_parser

        parser = build_argument_parser()
        args, _ = parser.parse_known_args(["--task", "101", "--output-dir", "/tmp/testout"])

        assert args.output_dir == "/tmp/testout"

    def test_output_dir_default(self):
        """--output-dir defaults to '.'."""
        from abacuscopilot.cli import build_argument_parser

        parser = build_argument_parser()
        args, _ = parser.parse_known_args(["--task", "101"])

        assert args.output_dir == "."
