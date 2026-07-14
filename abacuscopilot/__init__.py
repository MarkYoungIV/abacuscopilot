"""
AbacusCopilot: A pre- and post-processing toolkit for the ABACUS DFT software.

Inspired by VASPKIT, AbacusCopilot provides both an interactive menu-driven
interface and a command-line task mode for generating ABACUS input files,
analyzing calculation outputs, and producing publication-quality figures.
"""

__version__ = "0.1.7"
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
    "task",
]
