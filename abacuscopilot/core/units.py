"""Unit conversion utilities for abacuscopilot."""

from __future__ import annotations

import numpy as np

from abacuscopilot.core.constants import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_RY,
    RY_TO_EV,
)


def bohr_to_angstrom(value: float | np.ndarray) -> float | np.ndarray:
    """Convert Bohr to Angstrom."""
    return value * BOHR_TO_ANGSTROM


def angstrom_to_bohr(value: float | np.ndarray) -> float | np.ndarray:
    """Convert Angstrom to Bohr."""
    return value * ANGSTROM_TO_BOHR


def ry_to_ev(value: float | np.ndarray) -> float | np.ndarray:
    """Convert Rydberg to eV."""
    return value * RY_TO_EV


def ev_to_ry(value: float | np.ndarray) -> float | np.ndarray:
    """Convert eV to Rydberg."""
    return value * EV_TO_RY


def ha_to_ev(value: float | np.ndarray) -> float | np.ndarray:
    """Convert Hartree to eV."""
    from abacuscopilot.core.constants import HA_TO_EV
    return value * HA_TO_EV


def ev_to_ha(value: float | np.ndarray) -> float | np.ndarray:
    """Convert eV to Hartree."""
    from abacuscopilot.core.constants import EV_TO_HA
    return value * EV_TO_HA


def kbar_to_gpa(value: float | np.ndarray) -> float | np.ndarray:
    """Convert kbar to GPa."""
    from abacuscopilot.core.constants import KBAR_TO_GPA
    return value * KBAR_TO_GPA


def gpa_to_kbar(value: float | np.ndarray) -> float | np.ndarray:
    """Convert GPa to kbar."""
    from abacuscopilot.core.constants import GPA_TO_KBAR
    return value * GPA_TO_KBAR
