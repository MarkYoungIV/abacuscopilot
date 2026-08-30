"""Tests for the task registry, decorator, and discovery."""

import pytest

from abacuscopilot.core.exceptions import TaskNotFoundError
from abacuscopilot.tasks import TaskMetadata, TaskRegistry, task


class TestTaskDecorator:
    def test_register_task(self):
        registry = TaskRegistry()
        # Use unique ID to avoid collision with real tasks
        tid = 99991

        @registry.register(tid, "Test", "Test Task",
                          description="A test task for unit tests")
        def dummy_handler(args=None, interactive=True):
            return "done"

        assert tid in registry._tasks
        assert registry._tasks[tid].name == "Test Task"
        assert registry._tasks[tid].category == "Test"
        assert registry._tasks[tid].description == "A test task for unit tests"

        # Clean up
        del registry._tasks[tid]

    def test_duplicate_registration(self):
        registry = TaskRegistry()
        tid = 88881

        @registry.register(tid, "Test", "First")
        def first(args=None, interactive=True):
            pass

        with pytest.raises(ValueError, match="already registered"):
            @registry.register(tid, "Test", "Second")
            def second(args=None, interactive=True):
                pass

        # Clean up
        del registry._tasks[tid]

    def test_convenience_decorator(self):
        registry = TaskRegistry()
        tid = 77771

        @task(tid, "Test", "Convenience Test")
        def handler(args=None, interactive=True):
            return 42

        result = registry.dispatch(tid)
        assert result == 42

        # Clean up
        del registry._tasks[tid]


class TestTaskRegistry:
    def test_singleton(self):
        r1 = TaskRegistry()
        r2 = TaskRegistry()
        assert r1 is r2

    def test_dispatch_unknown(self):
        registry = TaskRegistry()
        # Skip discovery to speed up — the task doesn't exist
        old_discovered = registry._discovered
        registry._discovered = True
        try:
            with pytest.raises(TaskNotFoundError):
                registry.dispatch(99999)
        finally:
            registry._discovered = old_discovered

    def test_discover_modules(self):
        registry = TaskRegistry()
        # Ensure discovery has happened (may be no-op if already done)
        registry.discover_modules()

        assert len(registry._tasks) > 0
        task_ids = list(registry._tasks.keys())
        assert 101 in task_ids  # SCF INPUT
        assert 301 in task_ids  # Auto KPT
        assert 701 in task_ids  # SCF convergence
        assert 801 in task_ids  # Plot bands
        assert 901 in task_ids  # Plot DOS

    def test_get_categories(self):
        registry = TaskRegistry()
        registry.discover_modules()

        categories = registry.get_categories()
        assert "INPUT" in categories
        assert "KPT" in categories
        assert "Band Structure" in categories
        assert "DOS/PDOS" in categories
        assert "SCF Analysis" in categories

    def test_list_tasks_by_category(self):
        registry = TaskRegistry()
        registry.discover_modules()

        input_tasks = registry.list_tasks("INPUT")
        assert len(input_tasks) >= 6
        input_ids = [t.task_id for t in input_tasks]
        assert 101 in input_ids
        assert 102 in input_ids

    def test_get_tasks_by_category(self):
        registry = TaskRegistry()
        registry.discover_modules()

        by_cat = registry.get_tasks_by_category()
        assert "INPUT" in by_cat
        assert len(by_cat["INPUT"]) >= 6

    def test_get_task(self):
        registry = TaskRegistry()
        registry.discover_modules()

        t = registry.get_task(101)
        assert t is not None
        assert t.name == "SCF INPUT"
        assert t.category == "INPUT"

        t_none = registry.get_task(99999)
        assert t_none is None


class TestTaskMetadata:
    def test_creation(self):
        def dummy(args=None, interactive=True):
            pass

        meta = TaskMetadata(
            task_id=100,
            category="Test",
            name="Test Task",
            handler=dummy,
            description="desc",
            requires_files=["INPUT"],
        )
        assert meta.task_id == 100
        assert meta.name == "Test Task"
        assert meta.description == "desc"
        assert meta.requires_files == ["INPUT"]
        assert meta.handler is dummy
