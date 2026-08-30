"""Tests for atomic charge analysis tasks (1004 Mulliken / 1005 Bader / 1006 Hirshfeld)."""

from pathlib import Path

import numpy as np
import pytest

from abacuscopilot.core.models import Atom, Lattice, Structure
from abacuscopilot.postprocessing import charge_tasks as ct
from abacuscopilot.postprocessing.population_tasks import parse_hirshfeld_from_log


def _small_structure() -> Structure:
    """Two Li + two Cl in a 10 Bohr cubic box."""
    s = Structure()
    s.lattice = Lattice(constant=10.0, vectors=np.eye(3))
    bohr = 0.529177
    s.atoms = [
        Atom(species="Li", position=np.array([2, 2, 2]) * bohr),
        Atom(species="Li", position=np.array([8, 8, 8]) * bohr),
        Atom(species="Cl", position=np.array([2, 8, 8]) * bohr),
        Atom(species="Cl", position=np.array([8, 2, 2]) * bohr),
    ]
    s.species_order = ["Li", "Cl"]
    return s


class TestParseHirshfeld:
    def test_table_style(self, tmp_path):
        log = tmp_path / "running.log"
        log.write_text("""
some output
 Mulliken population:
  Li1  2.8421
 Hirshfeld charges:
  Li1  0.8765
  Cl1  -0.8765
""")
        data = parse_hirshfeld_from_log(log)
        assert data is not None
        assert data["method"] == "Hirshfeld"
        assert data["charges"] == {"Li1": 0.8765, "Cl1": -0.8765}

    def test_per_atom_lines(self, tmp_path):
        log = tmp_path / "running.log"
        log.write_text("""
 Hirshfeld charges of Atom 1 (Li)   =  0.876
 Hirshfeld charges of Atom 2 (Cl)   = -0.292
""")
        data = parse_hirshfeld_from_log(log)
        assert data["charges"] == {"Li": 0.876, "Cl": -0.292}

    def test_missing(self, tmp_path):
        log = tmp_path / "running.log"
        log.write_text("no hirshfeld here\n")
        assert parse_hirshfeld_from_log(log) is None


class TestReadCubeStandard:
    def test_roundtrip(self, tmp_path):
        rng = np.random.default_rng(0)
        data = rng.random((12, 14, 16)) + 0.5
        cell = np.eye(3) * 10.0
        origin = (0.0, 0.0, 0.0)
        path = tmp_path / "test.cube"
        ct._write_bader_cube(path, data, cell, origin, _small_structure())
        out, out_cell, out_origin = ct._read_cube_standard(path)
        assert out is not None
        # cube writes with %.6e → ~1e-6 relative precision
        np.testing.assert_allclose(out, data, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(out_cell, cell, rtol=1e-4, atol=1e-4)

    def test_abacus_z_fastest_transposed(self, tmp_path):
        """ABACUS CUBE header says 'Inner loop is z' → data is z-fastest, which is
        the OPPOSITE of the standard Gaussian x-fastest order.  Readers must
        transpose, else Bader/RESP/ESP silently read a scrambled density (water
        Bader gave O=8 / H=0 before this fix)."""
        nx, ny, nz = 2, 3, 4
        vals = [float(1000 * x + 100 * y + z)
                for x in range(nx) for y in range(ny) for z in range(nz)]
        lines = [
            "STEP: 0  Cubefile created from ABACUS. Inner loop is z, followed by y and x",
            "1 (nspin) 0.0",
            "1 0 0 0",
            f"{nx} 1.0 0 0", f"{ny} 0 1.0 0", f"{nz} 0 0 1.0",
            "1 1.0 0 0 0",  # dummy atom
        ]
        for i in range(0, len(vals), 6):
            lines.append("  " + "  ".join(f"{v:.1f}" for v in vals[i:i + 6]))
        path = tmp_path / "SPIN1_CHG.cube"
        path.write_text("\n".join(lines) + "\n")

        data, cell, _ = ct._read_cube_standard(path)
        assert data.shape == (nz, ny, nx)
        # data[z, y, x] must equal 1000*x + 100*y + z
        assert data[3, 2, 1] == 1000 * 1 + 100 * 2 + 3
        assert data[0, 1, 1] == 1000 * 1 + 100 * 1 + 0
        # legacy reader too
        data2, _, _ = ct.read_charge_cube(path)
        assert data2.shape == (nz, ny, nx)
        assert data2[3, 2, 1] == 1000 * 1 + 100 * 2 + 3

    def test_standard_cube_untouched(self, tmp_path):
        """A cube WITHOUT the ABACUS z-fastest marker reads as standard x-fastest."""
        nx, ny, nz = 2, 3, 4
        vals = [float(1000 * z + 100 * y + x)
                for z in range(nz) for y in range(ny) for x in range(nx)]
        lines = [
            "title", "subtitle", "1 0 0 0",
            f"{nx} 1.0 0 0", f"{ny} 0 1.0 0", f"{nz} 0 0 1.0",
            "1 1.0 0 0 0",
        ]
        for i in range(0, len(vals), 6):
            lines.append("  " + "  ".join(f"{v:.1f}" for v in vals[i:i + 6]))
        path = tmp_path / "std.cube"
        path.write_text("\n".join(lines) + "\n")
        data, _, _ = ct._read_cube_standard(path)
        assert data.shape == (nz, ny, nx)
        assert data[3, 2, 1] == 1000 * 3 + 100 * 2 + 1


class TestParseAcf:
    def test_parses_rows_in_order(self, tmp_path):
        acf = tmp_path / "ACF.dat"
        acf.write_text("""    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL
 --------------------------------------------------------------------------------
    1    0.000000    0.000000    0.000000    3.100000     0.000000   10.000000
    2    2.000000    2.000000    2.000000    7.200000     0.500000   12.000000
""")
        rows = ct._parse_acf(acf)
        assert len(rows) == 2
        assert rows[0][3] == pytest.approx(3.1)
        assert rows[1][3] == pytest.approx(7.2)


class TestAtomLabels:
    def test_global_numbering_matches_mulliken(self):
        """Bader/RESP atom labels must use STRU-order global numbering
        (C1..C24, H25..H36) so they line up with the 1301 Mulliken table."""
        from abacuscopilot.postprocessing.charge_tasks import _atom_labels
        # _small_structure() is [Li, Li, Cl, Cl]
        assert _atom_labels(_small_structure()) == ["Li1", "Li2", "Cl3", "Cl4"]


class TestZval:
    def test_reads_from_upf(self, tmp_path, monkeypatch):
        upf = tmp_path / "Sm3+_f--core-icmod1.PD04.PBE.UPF"
        upf.write_text('<UPF>\n z_valence="   11.00"\n</UPF>\n')
        monkeypatch.setattr(
            "abacuscopilot.config.load_config",
            lambda: {"libraries": {"pseudo_library": [str(tmp_path)]}},
        )
        assert ct._read_upf_zval(upf.name) == 11.0

    def test_fallback_table(self):
        assert ct._zval_for("Cl", _small_structure()) == 7.0


class TestTaskRegistration:
    def test_charge_tasks_in_population_menu(self):
        from abacuscopilot.tasks import TaskRegistry
        r = TaskRegistry()
        r.discover_modules()
        # Only Mulliken (1301) is active in Population.  Bader/Hirshfeld/RESP are
        # hidden (1304/1305/1306 unregistered) — see the notes above their tasks.
        assert sorted(t for t in r._tasks
                      if r._tasks[t].category == "Population") == [1301]
        assert 1301 in r._tasks
        for hidden in (1302, 1304, 1305, 1306):
            assert hidden not in r._tasks
        # No duplicate Mulliken under 1004 (1301 already covers it).
        assert 1004 not in r._tasks


class TestBaderE2E:
    def test_task_bader_charge(self, tmp_path, monkeypatch):
        if ct._locate_bader() is None:
            pytest.skip("bader.x not available (run setup.sh)")
        monkeypatch.chdir(tmp_path)

        # Write STRU with resolved PP filenames
        from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare
        structure = _small_structure()
        _write_stru_bare(structure, is_lcao=False, filepath="STRU")

        # Synthetic density: Gaussian blobs at each atom (Bohr positions)
        n = 40
        ax = np.linspace(0, 10, n)
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        density = np.zeros_like(X)
        for a in structure.atoms:
            p = a.position * ct.ANGSTROM_TO_BOHR
            density += np.exp(-((X - p[0]) ** 2 + (Y - p[1]) ** 2 + (Z - p[2]) ** 2) / 0.5)
        density = density.transpose(2, 1, 0)
        ct._write_bader_cube("SPIN1_CHG.cube", density, np.eye(3) * 10.0, (0.0, 0.0, 0.0), structure)

        ct.task_bader_charge(interactive=False)

        assert Path("Bader.dat").exists()
        lines = Path("Bader.dat").read_text().strip().splitlines()
        # header + 4 atoms
        assert len(lines) == 5


class TestResp:
    def test_esp_points_count(self):
        s = _small_structure()  # 4 atoms
        pts = ct._esp_points_for_atoms(s)
        # 4 atoms × 4 vdW shells × 60 points each
        assert len(pts) == 4 * 4 * 60
        assert np.isfinite(pts).all()

    def test_compute_esp_finite(self):
        s = _small_structure()
        rng = np.random.default_rng(4)
        rho = np.abs(rng.random((16, 16, 16))) + 0.01
        cell = np.eye(3) * 10.0
        origin = np.zeros(3)
        pts = ct._esp_points_for_atoms(s, n_per_shell=4)
        esp = ct._compute_esp(pts, s, rho, cell, origin)
        assert np.isfinite(esp).all()
        assert esp.shape == (len(pts),)

    def test_resp_fit_conserves_charge(self):
        s = _small_structure()
        rng = np.random.default_rng(5)
        pts = ct._esp_points_for_atoms(s, n_per_shell=6)
        esp = np.linalg.norm(pts - np.array([5.0, 5.0, 5.0]), axis=1)
        q = ct._resp_fit(pts, esp, s, total_charge=0.0)
        assert q.shape == (len(s.atoms),)
        assert abs(sum(q)) < 1e-6
        assert np.isfinite(q).all()


class TestAtomCartesian:
    def test_direct_to_cartesian(self):
        """Direct (fractional) STRU positions must convert to Cartesian Å."""
        s = Structure()
        s.lattice = Lattice(constant=ct.ANGSTROM_TO_BOHR, vectors=np.eye(3) * 10.0)  # 10 Å cube
        s.atoms = [Atom(species="O", position=np.array([0.5, 0.5, 0.5]))]
        s.species_order = ["O"]
        s.coordinate_type = "Direct"
        cart = ct._atom_cartesian_angstrom(s)
        np.testing.assert_allclose(cart, [[5.0, 5.0, 5.0]])

    def test_cartesian_passthrough(self):
        s = Structure()
        s.lattice = Lattice(constant=ct.ANGSTROM_TO_BOHR, vectors=np.eye(3) * 10.0)
        s.atoms = [Atom(species="O", position=np.array([2.0, 3.0, 4.0]))]
        s.species_order = ["O"]
        s.coordinate_type = "Cartesian_angstrom"
        cart = ct._atom_cartesian_angstrom(s)
        np.testing.assert_allclose(cart, [[2.0, 3.0, 4.0]])


class TestRespEquivalence:
    def _h2o(self):
        s = Structure()
        s.lattice = Lattice(constant=ct.ANGSTROM_TO_BOHR, vectors=np.eye(3) * 15.0)
        ang = 0.529177
        o = np.array([7.5, 7.5, 7.5])
        half = 104.5 / 2 * np.pi / 180
        s.atoms = [
            Atom(species="O", position=o * ang),
            Atom(species="H", position=(o + 0.96 * np.array([np.cos(half), np.sin(half), 0])) * ang),
            Atom(species="H", position=(o + 0.96 * np.array([np.cos(half), -np.sin(half), 0])) * ang),
        ]
        s.species_order = ["O", "H"]
        s.coordinate_type = "Cartesian_angstrom"
        return s

    def test_detect_equivalent_h(self):
        groups = ct._detect_equivalent_atoms(self._h2o())
        assert groups == [[1, 2]]

    def test_equiv_constraint_forces_equal(self):
        s = self._h2o()
        rng = np.random.default_rng(7)
        pts = ct._esp_points_for_atoms(s, n_per_shell=20)
        esp = np.linalg.norm(pts - np.array([7.5, 7.5, 7.5]) * ct.ANGSTROM_TO_BOHR, axis=1)
        q = ct._resp_fit(pts, esp, s, equiv_groups=[[1, 2]])
        assert abs(q[1] - q[2]) < 1e-8
        assert abs(sum(q)) < 1e-6


class TestParseMullikenFile:
    def test_parses_atom_blocks(self, tmp_path):
        from abacuscopilot.postprocessing.population_tasks import parse_mulliken_from_file
        f = tmp_path / "mulliken.txt"
        f.write_text("""STEP: 0
CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM
0 Zeta of C  Spin 1
Total Charge on atom:  C  4.0506
1 Zeta of C  Spin 1
Total Charge on atom:  C  4.0165
2 Zeta of H  Spin 1
Total Charge on atom:  H  0.9160
""")
        d = parse_mulliken_from_file(f)
        assert d is not None
        assert d["charges"] == {"C1": 4.0506, "C2": 4.0165, "H3": 0.9160}
        assert abs(d["total_charge"] - 8.9831) < 1e-6


class TestMullikenNetCharge:
    """Regression: 1301 Net Charge column must be computed, not '—'.

    Old code hard-coded the Net Charge cell to an em-dash.  For Mulliken/Lowdin
    the net charge is valence − population; for Hirshfeld the raw value already
    is the (net) charge, so it passes through unchanged.
    """

    def _render(self, monkeypatch, tmp_path, data) -> str:
        monkeypatch.chdir(tmp_path)
        from io import StringIO

        from rich.console import Console

        import abacuscopilot.postprocessing.population_tasks as pt

        buf = StringIO()
        console = Console(file=buf, width=200)
        monkeypatch.setattr(pt, "_get_console", lambda: console)
        pt._display_population_table(console, data)
        return buf.getvalue()

    def test_mulliken_net_charge_computed(self, monkeypatch, tmp_path):
        data = {
            "charges": {"C1": 4.0506, "H2": 0.9160},
            "total_charge": 4.9666,
            "method": "Mulliken",
        }
        out = self._render(monkeypatch, tmp_path, data)
        assert "—" not in out
        assert "Population (e)" in out
        assert "-0.050600" in out  # C: 4 − 4.0506
        assert "+0.084000" in out  # H: 1 − 0.9160

    def test_hirshfeld_net_charge_passthrough(self, monkeypatch, tmp_path):
        data = {
            "charges": {"Li1": 0.8765, "Cl1": -0.8765},
            "total_charge": 0.0,
            "method": "Hirshfeld",
        }
        out = self._render(monkeypatch, tmp_path, data)
        assert "—" not in out
        assert "Charge (e)" in out
        assert "+0.876500" in out
        assert "-0.876500" in out
