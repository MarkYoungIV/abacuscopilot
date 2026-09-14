"""Custom exception hierarchy for abacuscopilot."""


class AbacusCopilotError(Exception):
    """Base exception for all abacuscopilot errors."""


# Backward-compatible alias (software was formerly named "AbacusKit").
AbacusKitError = AbacusCopilotError


class FileFormatError(AbacusCopilotError):
    """Raised when an input file has invalid or unexpected format."""

    def __init__(self, filepath: str, message: str = ""):
        self.filepath = filepath
        super().__init__(f"Format error in '{filepath}': {message}")


class FileNotFoundError_(AbacusCopilotError):
    """Raised when a required input file is missing."""

    def __init__(self, filepath: str, message: str = ""):
        self.filepath = filepath
        super().__init__(f"File not found: '{filepath}'. {message}")


class MissingSectionError(FileFormatError):
    """Raised when a required section is missing from a file."""

    def __init__(self, filepath: str, section: str):
        super().__init__(filepath, f"Missing required section: '{section}'")


class ParameterError(AbacusCopilotError):
    """Raised for invalid parameter values."""

    def __init__(self, param: str, value, expected: str = ""):
        self.param = param
        self.value = value
        msg = f"Invalid value for '{param}': {value}"
        if expected:
            msg += f". Expected: {expected}"
        super().__init__(msg)


class TaskNotFoundError(AbacusCopilotError):
    """Raised when a requested task ID is not registered."""

    def __init__(self, task_id: int):
        self.task_id = task_id
        super().__init__(f"Task {task_id} not found. Use --list-tasks to see available tasks.")


class ConvergenceError(AbacusCopilotError):
    """Raised when SCF or geometry optimization fails to converge."""


class SymmetryError(AbacusCopilotError):
    """Raised when symmetry analysis fails (e.g., spglib not available)."""


class ASEImportError(AbacusCopilotError):
    """Raised when ASE-dependent functionality is used but ASE is not installed."""

    def __init__(self, feature: str = ""):
        msg = "This feature requires ASE. Install with: pip install ase"
        if feature:
            msg = f"'{feature}' requires ASE. Install with: pip install ase"
        super().__init__(msg)


class LibraryFamilyError(AbacusCopilotError):
    """Raised when the active user-configured library family cannot provide a
    required element — the auto-copy / STRU-sync is stopped instead of mixing
    files across series."""

    def __init__(self, family_label: str, missing_elements: list[str]):
        self.family_label = family_label
        self.missing_elements = sorted(set(missing_elements))
        elem_list = ", ".join(self.missing_elements)
        super().__init__(
            f"{family_label} has no files for element(s): {elem_list}. "
            "The auto-copy/STRU-sync was stopped to avoid mixing different library "
            "series.  Switch back to SG15 or add the missing files and retry."
        )
