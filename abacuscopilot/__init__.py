"""
AbacusCopilot: A pre- and post-processing copilot for the ABACUS DFT software.

It provides both an interactive menu-driven interface and a command-line task
mode for generating ABACUS input files, analyzing calculation outputs, and
producing publication-quality figures.
"""

__version__ = "0.1.34b"
__version_date__ = "2026-09-06"
__author__ = "AbacusCopilot Developers"
__license__ = "GPL-3.0"

from abacuscopilot.core.models import Atom, InputParams, KPoints, Lattice, Structure
from abacuscopilot.tasks import TaskRegistry, task

__all__ = [
    "Atom",
    "InputParams",
    "KPoints",
    "Lattice",
    "Structure",
    "TaskRegistry",
    "__version__",
    "__version_date__",
    "task",
]
