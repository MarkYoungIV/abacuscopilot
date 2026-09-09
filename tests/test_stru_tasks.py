"""Tests for STRU-writing helpers: library resolution and missing-species warnings."""

import numpy as np

from abacuscopilot.core.models import Atom, Structure
from abacuscopilot.preprocessing import input_tasks as it
from abacuscopilot.preprocessing import stru_tasks as st
from abacuscopilot.preprocessing.stru_tasks import _write_stru_bare


def _make_structure(species: list[str]) -> Structure:
    s = Structure()
    s.atoms = [Atom(species=sp, position=np.array([0.0, 0.0, 0.0])) for sp in species]
    s.species_order = list(species)
    return s


def _fake_config(upf_dir, orb_dir) -> dict:
    return {"libraries": {"pseudo_library": str(upf_dir), "orbital_library": str(orb_dir)}}


class TestWriteStruLibraryResolution:
    def test_resolves_filenames_from_library(self, tmp_path, monkeypatch):
        upf_dir = tmp_path / "upf"
        orb_dir = tmp_path / "orb"
        upf_dir.mkdir()
        orb_dir.mkdir()
        (upf_dir / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        (orb_dir / "Li_gga_7au_100Ry_2s2p1d.orb").write_text("dummy")

        monkeypatch.setattr(
            "abacuscopilot.config.load_config",
            lambda: _fake_config(upf_dir, orb_dir),
        )

        structure = _make_structure(["Li"])
        _write_stru_bare(structure, is_lcao=True, filepath=str(tmp_path / "STRU"))

        assert structure.pseudo_files["Li"] == "Li_ONCV_PBE-1.0.upf"
        assert structure.orbital_files["Li"] == "Li_gga_7au_100Ry_2s2p1d.orb"

    def test_warns_when_species_missing_from_library(self, tmp_path, monkeypatch, capsys):
        upf_dir = tmp_path / "upf"
        orb_dir = tmp_path / "orb"
        upf_dir.mkdir()
        orb_dir.mkdir()
        (upf_dir / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        (orb_dir / "Li_gga_7au_100Ry_2s2p1d.orb").write_text("dummy")

        monkeypatch.setattr(
            "abacuscopilot.config.load_config",
            lambda: _fake_config(upf_dir, orb_dir),
        )

        structure = _make_structure(["Li", "Sm"])
        _write_stru_bare(structure, is_lcao=True, filepath=str(tmp_path / "STRU"))

        # Li resolved from library; Sm falls back to placeholder filenames
        assert structure.pseudo_files["Li"].startswith("Li_")
        assert structure.orbital_files["Li"].startswith("Li_")
        assert structure.pseudo_files["Sm"] == "Sm.upf"
        assert structure.orbital_files["Sm"] == "Sm.orb"

        out = capsys.readouterr().out
        assert "Sm" in out
        assert "not found" in out
        assert "pseudopotential" in out
        assert "orbital" in out

    def test_no_warning_when_all_present(self, tmp_path, monkeypatch, capsys):
        upf_dir = tmp_path / "upf"
        orb_dir = tmp_path / "orb"
        upf_dir.mkdir()
        orb_dir.mkdir()
        (upf_dir / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        (upf_dir / "Sm_ONCV_PBE-1.0.upf").write_text("dummy")
        (orb_dir / "Li_gga_7au_100Ry_2s2p1d.orb").write_text("dummy")
        (orb_dir / "Sm_gga_10au_100Ry_2s2p1d.orb").write_text("dummy")

        monkeypatch.setattr(
            "abacuscopilot.config.load_config",
            lambda: _fake_config(upf_dir, orb_dir),
        )

        structure = _make_structure(["Li", "Sm"])
        _write_stru_bare(structure, is_lcao=True, filepath=str(tmp_path / "STRU"))

        assert structure.pseudo_files["Sm"].startswith("Sm_")
        assert structure.orbital_files["Sm"].startswith("Sm_")
        assert "not found" not in capsys.readouterr().out


class TestPdbToStru:
    WATER = """HETATM    1  O           0      -0.538   3.372   0.000                       O
HETATM    2  H           0       0.422   3.372   0.000                       H
HETATM    3  H           0      -0.858   4.277   0.000                       H
END
"""

    def test_read_pdb_no_box(self, tmp_path):
        from abacuscopilot.preprocessing.stru_tasks import _read_pdb
        pdb = tmp_path / "water.pdb"
        pdb.write_text(self.WATER)
        atoms = _read_pdb(pdb)
        assert atoms is not None
        assert len(atoms) == 3
        assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]
        assert atoms.cell.volume < 1e-6  # no CRYST1

    def test_read_pdb_with_cryst1(self, tmp_path):
        from abacuscopilot.preprocessing.stru_tasks import _read_pdb
        pdb = tmp_path / "box.pdb"
        pdb.write_text(
            "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
            + self.WATER
        )
        atoms = _read_pdb(pdb)
        assert abs(atoms.cell.volume - 8000.0) < 1e-6

    def test_center_molecule_in_box(self, tmp_path):
        from abacuscopilot.preprocessing.stru_tasks import _center_molecule_in_box, _read_pdb
        pdb = tmp_path / "water.pdb"
        pdb.write_text(self.WATER)
        atoms = _read_pdb(pdb)
        _center_molecule_in_box(atoms, 15.0)
        assert atoms.cell.volume > 1e-6
        pos = atoms.get_positions()
        # bounding box centered at the box center (7.5, 7.5, 7.5)
        assert np.allclose((pos.min(axis=0) + pos.max(axis=0)) / 2.0, 7.5, atol=1e-6)


class _FullCalcConsole:
    def __init__(self):
        self.out: list[str] = []

    def print(self, *args, **kwargs):
        self.out.append(" ".join(str(a) for a in args))


class _FullCalcAnswers:
    """Pop answers for _prompt_choice in call order; else use default."""

    def __init__(self, choices):
        self.choices = list(choices)

    def choice(self, console, question, options, default):
        return self.choices.pop(0) if self.choices else default


class TestFullCalcSharedInput:
    """Full-calc (201/202/207) must route per-calculation questions AND the
    post-INPUT pipeline through the shared input_tasks machinery, so its INPUT
    is consistent with running the matching INPUT-module task."""

    def test_relax_default_pbe_calls_auto_prepare(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        console = _FullCalcConsole()
        # basis, calc, server; XC/D3 fall back to defaults → PBE + No.
        a = _FullCalcAnswers(
            ["lcao", "relax (atoms only)", "CPU (genelpa)"]
        )
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(
            it, "_prompt",
            lambda c, q, default=None: default,
        )

        captured: dict = {}

        def _fake_auto_prepare(console_, params, interactive=True):
            captured["params"] = params

        monkeypatch.setattr(it, "_auto_prepare_files", _fake_auto_prepare)

        st._run_full_calculation_setup(console, interactive=True)

        assert "params" in captured
        params = captured["params"]
        assert params.calculation == "relax"
        assert params.basis_type == "lcao"
        # Default PBE → no functional override; vdw stays none.
        assert params.dft_functional == "pbe"
        assert params.vdw_method == "none"
        assert "dft_functional" not in params.extras.get("_template_keys", [])

        # INPUT must NOT carry an active dft_functional line (ABACUS PBE default),
        # but must keep the pbesol hint comment.
        input_text = (tmp_path / "INPUT").read_text()
        active = [ln for ln in input_text.splitlines() if ln.startswith("dft_functional")]
        assert active == []
        assert "#dft_functional" in input_text

    def test_relax_pbesol_writes_override(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        console = _FullCalcConsole()
        # Explicit PBEsol, D3 default No.
        a = _FullCalcAnswers(
            ["lcao", "relax (atoms only)", "CPU (genelpa)", "PBEsol", "No"]
        )
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(
            it, "_prompt",
            lambda c, q, default=None: default,
        )

        captured: dict = {}

        def _fake_auto_prepare(console_, params, interactive=True):
            captured["params"] = params

        monkeypatch.setattr(it, "_auto_prepare_files", _fake_auto_prepare)

        st._run_full_calculation_setup(console, interactive=True)

        params = captured["params"]
        assert params.dft_functional == "pbesol"
        assert "dft_functional" in params.extras.get("_template_keys", [])

        # INPUT carries an active pbesol line; the pbesol hint is dropped.
        input_text = (tmp_path / "INPUT").read_text()
        active = [ln for ln in input_text.splitlines() if ln.startswith("dft_functional")]
        assert active and active[0].endswith("pbesol")
        assert "#dft_functional" not in input_text

    def test_md_branch_dispatches_to_md_thermo(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        console = _FullCalcConsole()
        # basis, calc=MD, server, then MD ensemble default nvt.
        a = _FullCalcAnswers(
            ["lcao", "MD", "CPU (genelpa)"]
        )
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(
            it, "_prompt",
            lambda c, q, default=None: default,
        )

        captured: dict = {}

        def _fake_auto_prepare(console_, params, interactive=True):
            captured["params"] = params

        monkeypatch.setattr(it, "_auto_prepare_files", _fake_auto_prepare)

        st._run_full_calculation_setup(console, interactive=True)

        params = captured["params"]
        assert params.calculation == "md"
        assert params.md_type == "nvt"
        assert params.md_nstep == 10000
        assert params.md_dt == 1.0
