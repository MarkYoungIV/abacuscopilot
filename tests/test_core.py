"""Tests for core models, constants, units, and exceptions."""

import numpy as np
import pytest

from abacuscopilot.core.constants import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HA,
    EV_TO_RY,
    HA_TO_EV,
    RY_TO_EV,
)
from abacuscopilot.core.exceptions import (
    AbacusCopilotError,
    ASEImportError,
    FileFormatError,
    FileNotFoundError_,
    MissingSectionError,
    ParameterError,
    TaskNotFoundError,
)
from abacuscopilot.core.models import Atom, InputParams, KPoints, Lattice, Structure
from abacuscopilot.core.units import (
    angstrom_to_bohr,
    bohr_to_angstrom,
    ev_to_ry,
    gpa_to_kbar,
    ha_to_ev,
    kbar_to_gpa,
    ry_to_ev,
)

# =============================================================================
# Constants
# =============================================================================

class TestConstants:
    def test_energy_conversions(self):
        assert abs(RY_TO_EV - 13.605693122994) < 1e-10
        assert abs(EV_TO_RY * RY_TO_EV - 1.0) < 1e-10
        assert abs(HA_TO_EV - 27.211386245988) < 1e-10
        assert abs(EV_TO_HA * HA_TO_EV - 1.0) < 1e-10

    def test_length_conversions(self):
        assert abs(BOHR_TO_ANGSTROM - 0.529177210903) < 1e-10
        assert abs(ANGSTROM_TO_BOHR * BOHR_TO_ANGSTROM - 1.0) < 1e-10


# =============================================================================
# Unit conversions
# =============================================================================

class TestUnits:
    def test_bohr_angstrom(self):
        assert abs(bohr_to_angstrom(1.0) - BOHR_TO_ANGSTROM) < 1e-10
        assert abs(angstrom_to_bohr(1.0) - ANGSTROM_TO_BOHR) < 1e-10
        assert abs(angstrom_to_bohr(bohr_to_angstrom(5.0)) - 5.0) < 1e-10

    def test_energy_units(self):
        assert abs(ry_to_ev(1.0) - RY_TO_EV) < 1e-10
        assert abs(ev_to_ry(ry_to_ev(10.0)) - 10.0) < 1e-10
        assert abs(ha_to_ev(1.0) - HA_TO_EV) < 1e-10

    def test_pressure_units(self):
        assert abs(kbar_to_gpa(10.0) - 1.0) < 1e-10
        assert abs(gpa_to_kbar(1.0) - 10.0) < 1e-10

    def test_array_conversions(self):
        arr = np.array([1.0, 2.0, 3.0])
        result = bohr_to_angstrom(arr)
        assert isinstance(result, np.ndarray)
        assert result.shape == arr.shape
        assert abs(result[0] - BOHR_TO_ANGSTROM) < 1e-10


# =============================================================================
# Lattice
# =============================================================================

class TestLattice:
    def test_default_lattice(self):
        lat = Lattice()
        assert lat.constant == 1.0
        np.testing.assert_array_equal(lat.vectors, np.eye(3))

    def test_cell_property(self):
        lat = Lattice(constant=10.0, vectors=np.eye(3))
        np.testing.assert_array_equal(lat.cell, 10.0 * np.eye(3))

    def test_volume_cubic(self):
        lat = Lattice(constant=5.0, vectors=np.eye(3))
        assert abs(lat.volume - 125.0) < 1e-10

    def test_reciprocal_cell(self):
        lat = Lattice(constant=1.0, vectors=np.eye(3))
        expected = 2.0 * np.pi * np.eye(3)
        np.testing.assert_array_almost_equal(lat.reciprocal_cell, expected)

    def test_from_cell(self):
        cell = 5.0 * np.eye(3)
        lat = Lattice.from_cell(cell, constant=5.0)
        np.testing.assert_array_almost_equal(lat.vectors, np.eye(3))
        np.testing.assert_array_almost_equal(lat.cell, cell)

    def test_from_cell_angstrom(self):
        cell_ang = 5.0 * np.eye(3)
        lat = Lattice.from_cell_angstrom(cell_ang)
        expected_bohr = 5.0 * ANGSTROM_TO_BOHR
        np.testing.assert_array_almost_equal(lat.cell, expected_bohr * np.eye(3))

    def test_invalid_shape(self):
        with pytest.raises(ValueError):
            Lattice(vectors=np.zeros((2, 3)))


# =============================================================================
# Atom
# =============================================================================

class TestAtom:
    def test_defaults(self):
        atom = Atom(species="Si", position=np.array([0.0, 0.0, 0.0]))
        assert atom.species == "Si"
        assert atom.fix == (True, True, True)
        assert atom.magmom == 0.0
        assert atom.velocity is None
        assert atom.angle1 is None

    def test_custom(self):
        atom = Atom(
            species="Fe",
            position=np.array([0.5, 0.5, 0.5]),
            fix=(True, True, False),
            magmom=2.0,
        )
        assert atom.species == "Fe"
        assert atom.fix == (True, True, False)
        assert atom.magmom == 2.0

    def test_invalid_position(self):
        with pytest.raises(ValueError):
            Atom(species="H", position=np.array([0.0, 0.0]))


# =============================================================================
# Structure
# =============================================================================

class TestStructure:
    def test_empty(self):
        s = Structure()
        assert s.num_atoms == 0
        assert s.num_species == 0
        assert s.coordinate_type == "Direct"

    def test_basic(self):
        s = Structure()
        s.atoms = [
            Atom(species="Si", position=np.array([0.0, 0.0, 0.0])),
            Atom(species="Si", position=np.array([0.25, 0.25, 0.25])),
            Atom(species="O", position=np.array([0.5, 0.5, 0.5])),
        ]
        s.species_order = ["Si", "O"]
        assert s.num_atoms == 3
        assert s.num_species == 2
        assert s.count_species("Si") == 2
        assert s.count_species("O") == 1

    def test_get_atoms_by_species(self):
        s = Structure()
        s.atoms = [
            Atom(species="Si", position=np.array([0, 0, 0])),
            Atom(species="O", position=np.array([1, 1, 1])),
            Atom(species="Si", position=np.array([0.5, 0.5, 0.5])),
        ]
        s.species_order = ["Si", "O"]
        si_indices = s.get_atoms_by_species("Si")
        assert si_indices == [0, 2]

    def test_positions_array(self):
        s = Structure()
        s.atoms = [
            Atom(species="H", position=np.array([1.0, 2.0, 3.0])),
            Atom(species="H", position=np.array([4.0, 5.0, 6.0])),
        ]
        expected = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        np.testing.assert_array_equal(s.positions, expected)


# =============================================================================
# KPoints
# =============================================================================

class TestKPoints:
    def test_default(self):
        kp = KPoints()
        assert kp.mode == "gamma"
        assert kp.grid is None

    def test_auto_mp_string(self):
        kp = KPoints(mode="gamma", grid=(4, 4, 4), gamma_centered=True)
        s = kp.to_string()
        assert "K_POINTS" in s
        assert "Gamma" in s
        assert "4 4 4" in s

    def test_mp_string(self):
        kp = KPoints(mode="mp", grid=(6, 6, 4), gamma_centered=False,
                      shift=(0.5, 0.5, 0.5))
        s = kp.to_string()
        assert "MP" in s
        assert "0.5 0.5 0.5" in s

    def test_line_mode_string(self):
        kp = KPoints(mode="line")
        kp.line_path = [
            {"start": (0, 0, 0), "end": (0.5, 0, 0), "npoints": 20, "label": "G"},
        ]
        s = kp.to_string()
        assert "Line" in s or "LINE" in s

    def test_explicit_string(self):
        kp = KPoints(mode="direct", explicit_kpoints=[
            (0.0, 0.0, 0.0, 1.0),
            (0.5, 0.5, 0.5, 1.0),
        ])
        s = kp.to_string()
        assert "Direct" in s

    def test_invalid_mode(self):
        kp = KPoints(mode="unknown")
        with pytest.raises(ValueError):
            kp.to_string()

    def test_line_cartesian_cumulative_distances(self):
        """Regression test: line_cartesian mode must use Cartesian coords directly.

        This catches the self-assignment bug where ``end_cart = end_cart``
        was used instead of ``end_cart = end``.
        """
        kp = KPoints(mode="line_cartesian")
        kp.line_path = [
            {"start": (0.0, 0.0, 0.0), "end": (1.0, 0.0, 0.0), "npoints": 10,
             "label": "G", "end_label": "X"},
            {"start": (1.0, 0.0, 0.0), "end": (1.0, 1.0, 0.0), "npoints": 10,
             "label": "X", "end_label": "M"},
        ]
        # Identity reciprocal cell — distances are directly in Cartesian
        recip = np.eye(3)
        dists = kp.get_cumulative_distances(recip)
        assert len(dists) == 20  # 10 + 10 points
        # First segment: from (0,0,0) to (1,0,0) — length 1.0
        assert dists[0] == pytest.approx(0.0, abs=1e-10)
        assert dists[9] == pytest.approx(0.9, abs=1e-6)
        # Labels: first label + end_labels
        assert kp.labels == ["G", "X", "M"]
        assert kp.label_positions[0] == pytest.approx(0.0)
        assert kp.label_positions[1] == pytest.approx(1.0)
        assert kp.label_positions[2] == pytest.approx(2.0)

    def test_line_mode_cumulative_distances(self):
        """Sanity check: line mode with fractional coords converts correctly."""
        kp = KPoints(mode="line")
        kp.line_path = [
            {"start": (0.0, 0.0, 0.0), "end": (0.5, 0.0, 0.0), "npoints": 10,
             "label": "G", "end_label": "X"},
        ]
        recip = np.eye(3) * 2.0 * np.pi  # b1 = 2π x̂
        dists = kp.get_cumulative_distances(recip)
        assert len(dists) == 10
        assert kp.labels == ["G", "X"]


# =============================================================================
# InputParams
# =============================================================================

class TestInputParams:
    def test_defaults(self):
        p = InputParams()
        assert p.calculation == "scf"
        assert p.ecutwfc == 100.0
        assert p.scf_thr == 1e-7
        assert p.basis_type == "pw"

    def test_get_param(self):
        p = InputParams(ecutwfc=60.0)
        assert p.get_param("ecutwfc") == 60.0
        assert p.get_param("nonexistent", "default") == "default"

    def test_set_param(self):
        p = InputParams()
        p.set_param("ecutwfc", 80.0)
        assert p.ecutwfc == 80.0
        p.set_param("custom_param", "value")
        assert p.extras["custom_param"] == "value"

    def test_as_dict(self):
        p = InputParams(ecutwfc=80.0, calculation="relax", nspin=2)
        d = p.as_dict()
        assert d["ecutwfc"] == 80.0
        assert d["calculation"] == "relax"
        assert d["nspin"] == 2

    def test_extras(self):
        p = InputParams()
        p.extras["custom_key"] = "custom_value"
        assert p.get_param("custom_key") == "custom_value"
        d = p.as_dict()
        assert d["custom_key"] == "custom_value"


# =============================================================================
# Exceptions
# =============================================================================

class TestExceptions:
    def test_base(self):
        with pytest.raises(AbacusCopilotError):
            raise AbacusCopilotError("test")

    def test_file_format_error(self):
        e = FileFormatError("test.txt", "bad format")
        assert "test.txt" in str(e)
        assert "bad format" in str(e)

    def test_file_not_found(self):
        e = FileNotFoundError_("missing.txt", "gone")
        assert "missing.txt" in str(e)

    def test_missing_section(self):
        e = MissingSectionError("stru", "ATOMIC_SPECIES")
        assert "ATOMIC_SPECIES" in str(e)

    def test_parameter_error(self):
        e = ParameterError("ecutwfc", -1, "> 0")
        assert "ecutwfc" in str(e)
        assert str(-1) in str(e)

    def test_task_not_found(self):
        e = TaskNotFoundError(999)
        assert "999" in str(e)

    def test_ase_import_error(self):
        e = ASEImportError("to_ase")
        assert "to_ase" in str(e)
        assert "ASE" in str(e)
