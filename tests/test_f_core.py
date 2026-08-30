"""Tests for large-core (f-electron-pseudized) lanthanide awareness."""

from types import SimpleNamespace

import numpy as np

from abacuscopilot.core.models import Atom, Structure
from abacuscopilot.preprocessing.system_tasks import (
    _f_core_info_from_structure,
    adjust_ecutwfc_for_f_core,
    analyze_f_core,
    warn_f_core,
)


def _make_structure(species: list[str]) -> Structure:
    s = Structure()
    s.atoms = [Atom(species=sp, position=np.array([0.0, 0.0, 0.0])) for sp in species]
    s.species_order = list(species)
    return s


def _lanthanide_lib(tmp_path):
    """Build a mini library dir: an f-core Sm PP + a 300 Ry Sm orbital."""
    sm = tmp_path / "Sm"
    sm.mkdir()
    (sm / "Sm3+_f--core-icmod1.PD04.PBE.UPF").write_text("dummy")
    (sm / "Sm_gga_7au_300.0Ry_4s2p2d1f.orb").write_text("dummy")
    return str(tmp_path)


class TestAnalyzeFCore:
    def test_detects_f_core_and_orbital_cutoff(self, tmp_path):
        lib = _lanthanide_lib(tmp_path)
        info = analyze_f_core(_make_structure(["Sm"]), {
            "pseudo_library": lib, "orbital_library": lib,
        })
        assert [sp for sp, _ in info["f_core_species"]] == ["Sm"]
        assert info["max_orb_ecut"] == 300.0

    def test_normal_species_not_flagged(self, tmp_path):
        lib = _lanthanide_lib(tmp_path)
        (tmp_path / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        (tmp_path / "Li_gga_7au_100Ry_4s1p.orb").write_text("dummy")
        info = analyze_f_core(_make_structure(["Li"]), {
            "pseudo_library": lib, "orbital_library": lib,
        })
        assert info["f_core_species"] == []
        assert info["max_orb_ecut"] == 100.0

    def test_from_resolved_files(self):
        s = _make_structure(["Sm"])
        s.pseudo_files["Sm"] = "Sm3+_f--core-icmod1.PD04.PBE.UPF"
        s.orbital_files["Sm"] = "Sm_gga_7au_300.0Ry_4s2p2d1f.orb"
        info = _f_core_info_from_structure(s)
        assert info["f_core_species"] == [("Sm", "Sm3+_f--core-icmod1.PD04.PBE.UPF")]
        assert info["max_orb_ecut"] == 300.0

    def test_accepts_plain_species_list(self, tmp_path):
        """analyze_f_core works with a list of element symbols (not just Structure)."""
        lib = _lanthanide_lib(tmp_path)
        info = analyze_f_core(["Sm"], {"pseudo_library": lib, "orbital_library": lib})
        assert [sp for sp, _ in info["f_core_species"]] == ["Sm"]
        assert info["max_orb_ecut"] == 300.0


class TestWarnFCore:
    def test_prints_suitability_notice(self, capsys):
        from abacuscopilot.console_utils import _get_console
        info = {"f_core_species": [("Sm", "Sm3+_f--core-icmod1.PD04.PBE.UPF")],
                "ecut_entries": [], "max_orb_ecut": 300.0}
        warn_f_core(_get_console(), info)
        out = capsys.readouterr().out
        assert "大核赝势" in out
        assert "f 电子赝化进芯" in out
        assert "不适合" in out
        assert "单质" in out
        assert "300" in out

    def test_silent_without_f_core(self, capsys):
        from abacuscopilot.console_utils import _get_console
        warn_f_core(_get_console(), {"f_core_species": [], "max_orb_ecut": 0.0})
        assert capsys.readouterr().out == ""


class TestAdjustEcutwfc:
    def test_raises_to_orbital_cutoff(self, monkeypatch):
        from abacuscopilot.preprocessing import system_tasks
        monkeypatch.setattr(
            system_tasks, "_prompt_choice",
            lambda console, q, options, default: options[0],  # "Set ecutwfc = 300 Ry"
        )
        params = SimpleNamespace(ecutwfc=100.0)
        info = {"f_core_species": [("Sm", "Sm3+_...UPF")],
                "max_orb_ecut": 300.0}
        changed = adjust_ecutwfc_for_f_core(object(), params, info)
        assert changed is True
        assert params.ecutwfc == 300

    def test_keeps_current_on_decline(self, monkeypatch):
        from abacuscopilot.preprocessing import system_tasks
        monkeypatch.setattr(
            system_tasks, "_prompt_choice",
            lambda console, q, options, default: "Keep current",
        )
        params = SimpleNamespace(ecutwfc=100.0)
        info = {"f_core_species": [("Sm", "Sm3+_...UPF")],
                "max_orb_ecut": 300.0}
        changed = adjust_ecutwfc_for_f_core(object(), params, info)
        assert changed is False
        assert params.ecutwfc == 100.0

    def test_auto_raises_when_noninteractive(self, monkeypatch):
        """CLI/non-interactive runs raise ecutwfc without prompting."""
        from abacuscopilot.preprocessing import system_tasks
        called = []

        def fake_prompt(console, q, options, default):
            called.append(1)
            return options[0]

        monkeypatch.setattr(system_tasks, "_prompt_choice", fake_prompt)
        params = SimpleNamespace(ecutwfc=100.0)
        info = {"f_core_species": [("Sm", "Sm3+_...UPF")],
                "max_orb_ecut": 300.0}
        changed = adjust_ecutwfc_for_f_core(object(), params, info, interactive=False)
        assert changed is True
        assert params.ecutwfc == 300
        assert not called  # no prompt in non-interactive mode

    def test_noop_when_cutoff_already_high_enough(self, monkeypatch):
        from abacuscopilot.preprocessing import system_tasks
        called = []

        def fake_prompt(console, q, options, default):
            called.append(1)
            return options[0]

        monkeypatch.setattr(system_tasks, "_prompt_choice", fake_prompt)
        params = SimpleNamespace(ecutwfc=300.0)
        info = {"f_core_species": [("Sm", "Sm3+_...UPF")],
                "max_orb_ecut": 300.0}
        assert adjust_ecutwfc_for_f_core(object(), params, info) is False
        assert not called  # no prompt needed

    def test_noop_without_f_core(self, monkeypatch):
        params = SimpleNamespace(ecutwfc=100.0)
        info = {"f_core_species": [], "max_orb_ecut": 0.0}
        assert adjust_ecutwfc_for_f_core(object(), params, info) is False
        assert params.ecutwfc == 100.0
