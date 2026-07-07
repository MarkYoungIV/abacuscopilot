"""Tests for INPUT, STRU, and KPT reader/writer modules."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from abacuscopilot.core.models import Atom, InputParams, KPoints, Lattice, Structure
from abacuscopilot.io.input_file import (
    _format_value,
    _parse_value,
    read_input,
    validate_input,
    write_input,
)
from abacuscopilot.io.kpt_file import (
    auto_mp_kpts,
    line_mode_kpts_from_path,
    read_kpt,
    write_kpt,
)
from abacuscopilot.io.stru_file import read_stru, write_stru

# =============================================================================
# INPUT file tests
# =============================================================================

class TestInputReadWrite:
    def test_parse_value_bool(self):
        assert _parse_value("True") is True
        assert _parse_value("false") is False
        assert _parse_value("1") is True
        assert _parse_value("0") is False

    def test_parse_value_int(self):
        assert _parse_value("42") == 42
        assert _parse_value("-10") == -10

    def test_parse_value_float(self):
        assert abs(_parse_value("1.5") - 1.5) < 1e-10
        assert abs(_parse_value("1e-7") - 1e-7) < 1e-10

    def test_parse_value_string(self):
        assert _parse_value("hello") == "hello"
        assert _parse_value('"quoted"') == "quoted"

    def test_format_value(self):
        assert _format_value(True) == "1"
        assert _format_value(False) == "0"
        assert _format_value(42) == "42"
        assert _format_value(1.5) == "1.5"

    def test_write_and_read_roundtrip(self):
        """Write an INPUT file and read it back."""
        params = InputParams(
            calculation="scf",
            ecutwfc=80.0,
            nspin=2,
            scf_thr=1e-6,
            out_chg=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "INPUT"
            write_input(params, fpath)

            # Read back
            loaded = read_input(fpath)
            assert loaded.calculation == "scf"
            assert loaded.ecutwfc == 80.0
            assert loaded.nspin == 2
            assert loaded.scf_thr == 1e-6
            assert loaded.out_chg is True

    def test_write_with_extras(self):
        params = InputParams()
        params.extras["custom_param"] = "value"
        params.extras["another"] = 123

        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "INPUT"
            write_input(params, fpath)

            loaded = read_input(fpath)
            assert loaded.extras.get("custom_param") == "value"
            assert loaded.extras.get("another") == 123

    def test_write_preserves_comment(self):
        params = InputParams()
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "INPUT"
            write_input(params, fpath, comment="Test calculation")
            content = fpath.read_text()
            assert "Test calculation" in content

    def test_read_missing_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "INPUT"
            fpath.write_text("wrong_header\n")
            with pytest.raises(Exception):
                read_input(fpath)


class TestInputValidation:
    def test_valid_params(self):
        params = InputParams(ntype=2, ecutwfc=100.0, scf_thr=1e-7)
        warnings = validate_input(params)
        assert warnings == []

    def test_negative_ntype(self):
        params = InputParams(ntype=0)
        warnings = validate_input(params)
        assert any("ntype" in w for w in warnings)

    def test_unknown_calculation(self):
        params = InputParams(calculation="impossible")
        warnings = validate_input(params)
        assert any("calculation" in w for w in warnings)

    def test_unknown_basis(self):
        params = InputParams(basis_type="unknown")
        warnings = validate_input(params)
        assert any("basis_type" in w for w in warnings)

    def test_zero_ecutwfc(self):
        params = InputParams(ecutwfc=0)
        warnings = validate_input(params)
        assert any("ecutwfc" in w for w in warnings)

    def test_lcao_solver_check(self):
        params = InputParams(basis_type="lcao", ks_solver="unknown_solver")
        warnings = validate_input(params)
        assert any("ks_solver" in w for w in warnings)


# =============================================================================
# STRU file tests
# =============================================================================

SAMPLE_STRU = """ATOMIC_SPECIES
Si  28.085  Si.upf

NUMERICAL_ORBITAL
Si.orb

LATTICE_CONSTANT
  10.263

LATTICE_VECTORS
  0.0000000000  0.5000000000  0.5000000000
  0.5000000000  0.0000000000  0.5000000000
  0.5000000000  0.5000000000  0.0000000000

ATOMIC_POSITIONS
Direct

Si
0.0
2
  0.0000000000  0.0000000000  0.0000000000  0 0 0
  0.2500000000  0.2500000000  0.2500000000  1 1 1
"""


class TestStruReadWrite:
    def test_read_sample(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "STRU"
            fpath.write_text(SAMPLE_STRU)

            structure = read_stru(fpath)
            assert structure.num_atoms == 2
            assert structure.num_species == 1
            assert structure.species_order == ["Si"]
            assert structure.coordinate_type == "Direct"
            assert "Si" in structure.pseudo_files
            assert structure.pseudo_files["Si"] == "Si.upf"
            assert abs(structure.lattice.constant - 10.263) < 1e-6

    def test_roundtrip(self):
        """Write and read back a Structure."""
        structure = Structure()
        structure.lattice = Lattice(constant=10.0, vectors=np.eye(3))
        structure.atoms = [
            Atom(species="Si", position=np.array([0.0, 0.0, 0.0])),
            Atom(species="Si", position=np.array([0.25, 0.25, 0.25])),
        ]
        structure.species_order = ["Si"]
        structure.pseudo_files = {"Si": "Si.upf"}
        structure.coordinate_type = "Direct"

        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "STRU"
            write_stru(structure, fpath, is_lcao=False)

            loaded = read_stru(fpath)
            assert loaded.num_atoms == 2
            assert loaded.coordinate_type == "Direct"
            assert "Si" in loaded.pseudo_files

    def test_read_missing_file(self):
        with pytest.raises(Exception):
            read_stru("nonexistent_file_12345.str")

    def test_read_missing_section(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "STRU"
            fpath.write_text("ATOMIC_SPECIES\nSi 0 Si.upf\n")
            with pytest.raises(Exception):
                read_stru(fpath)


# =============================================================================
# KPT file tests
# =============================================================================

class TestKptReadWrite:
    def test_read_auto_mp(self):
        content = "K_POINTS\n0\nGamma\n4 4 4  0 0 0\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "KPT"
            fpath.write_text(content)

            kpts = read_kpt(fpath)
            assert kpts.mode == "gamma"
            assert kpts.grid == (4, 4, 4)
            assert kpts.gamma_centered is True

    def test_read_standard_mp(self):
        content = "KPOINTS\n0\nMP\n6 6 4  0.5 0.5 0.5\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "KPT"
            fpath.write_text(content)

            kpts = read_kpt(fpath)
            assert kpts.grid == (6, 6, 4)
            assert kpts.gamma_centered is False
            assert kpts.shift == (0.5, 0.5, 0.5)

    def test_read_explicit(self):
        content = "K_POINTS\n2\nDirect\n  0.0  0.0  0.0  1.0\n  0.5  0.5  0.5  1.0\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "KPT"
            fpath.write_text(content)

            kpts = read_kpt(fpath)
            assert kpts.mode == "direct"
            assert len(kpts.explicit_kpoints) == 2

    def test_read_line_mode(self):
        content = """K_POINTS
5
Line
  0.000  0.000  0.000  20  G
  0.500  0.000  0.000  20  X
  0.500  0.500  0.000  20  M
  0.500  0.500  0.500  20  R
  0.000  0.000  0.000  1   G
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "KPT"
            fpath.write_text(content)

            kpts = read_kpt(fpath)
            assert kpts.mode == "line"
            assert len(kpts.line_path) == 4

    def test_roundtrip_auto(self):
        kpts = KPoints(mode="gamma", grid=(4, 4, 4), gamma_centered=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "KPT"
            write_kpt(kpts, fpath)

            loaded = read_kpt(fpath)
            assert loaded.grid == (4, 4, 4)
            assert loaded.gamma_centered is True

    def test_auto_mp_kpts(self):
        lat = Lattice(constant=1.0, vectors=np.eye(3))
        kpts = auto_mp_kpts(lat, kspacing=0.04)
        assert kpts.mode == "gamma"
        assert kpts.grid is not None
        assert all(n >= 1 for n in kpts.grid)

    def test_line_mode_from_path(self):
        path = [
            ([0, 0, 0], [0.5, 0, 0], 20),
            ([0.5, 0, 0], [0.5, 0.5, 0], 20),
        ]
        labels = ["G", "X", "M"]
        kpts = line_mode_kpts_from_path(path, labels=labels)
        assert kpts.mode == "line"
        assert len(kpts.line_path) == 2
