"""Tests for the library "family" abstraction (ABACUS-APNS support).

Covers the config schema additions, deterministic per-element defaults for the
APNS efficiency/precision orbital sets, materializing a family into the flat
pseudo_library/orbital_library lists, and the interactive picker.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import abacuscopilot.library_families as lf
from abacuscopilot.config import (
    DEFAULT_CONFIG,
    _path_within,
    _registered_family_dirs,
)
from abacuscopilot.preprocessing.system_tasks import (
    _candidate_rank,
    _find_file_for_element,
    _find_file_in_libraries,
    missing_element_blockers,
    missing_library_files,
    orbital_rank_mode,
)


def _mkdir_with(tmp_path, *names: str) -> Path:
    d = tmp_path / "lib"
    d.mkdir(exist_ok=True)
    for n in names:
        (d / n).write_text("dummy")
    return d


class FakeConsole:
    """Minimal Rich-console stand-in: answers fed to input(), output captured."""

    def __init__(self, answers=()):
        self.answers = list(answers)
        self.out: list[str] = []

    def print(self, *args, **kwargs):
        self.out.append(" ".join(str(a) for a in args))

    def input(self, prompt: str = "") -> str:
        self.out.append(str(prompt))
        return self.answers.pop(0) if self.answers else ""


def _dojo_config(pseudo_dir: str, orb_root: str, family: str = "sg15") -> dict:
    """Dojo-NC-FR config: all three tiers share the Orbitals_v2.0 root dir."""
    return {
        "libraries": {
            "pseudo_library": [],
            "orbital_library": [],
            "family": family,
            "families": {
                "dojoncfr": {
                    "pseudo_dir": pseudo_dir,
                    "orbital_dirs": {"sz": orb_root, "dzp": orb_root,
                                     "tzdp": orb_root},
                }
            },
        }
    }


def _apns_config(pseudo_dir: str, eff_dir: str, family: str = "sg15") -> dict:
    return {
        "libraries": {
            "pseudo_library": [],
            "orbital_library": [],
            "family": family,
            "families": {
                "apns": {
                    "pseudo_dir": pseudo_dir,
                    "orbital_dirs": {"efficiency": eff_dir},
                }
            },
        }
    }


class TestConfigSchema:
    def test_default_family_is_sg15(self):
        libs = DEFAULT_CONFIG["libraries"]
        assert libs["family"] == "sg15"
        assert libs["families"] == {}
        # The legacy flat lists are still present (source of truth).
        assert "pseudo_library" in libs
        assert "orbital_library" in libs
        # rcut-copy policy defaults to the canonical SG15 7 au.
        assert libs["rcut_policy"] == "7"

    def test_family_accessors(self):
        assert lf.current_family({}) == "sg15"
        assert lf.current_family({"libraries": {"family": "apns/precision"}}) == "apns/precision"
        assert lf.is_user_family("sg15") is False
        assert lf.is_user_family("custom") is False
        assert lf.is_user_family("apns") is True
        assert lf.is_user_family("apns/precision") is True


    def test_known_families_include_dojoncfr(self):
        spec = lf.KNOWN_FAMILIES["dojoncfr"]
        assert spec.label.startswith("Dojo-NC-FR")
        # dzp is the recommended default and must be the blank/Enter default.
        assert list(spec.variants) == ["dzp", "sz", "tzdp"]
        assert spec.pseudo_key[1:] == ("families", "dojoncfr", "pseudo_dir")

    def test_family_accessors_dojo(self):
        assert lf.current_family({"libraries": {"family": "dojoncfr/dzp"}}) == "dojoncfr/dzp"
        assert lf.is_user_family("dojoncfr") is True
        assert lf.is_user_family("dojoncfr/tzdp") is True
        assert "Dojo-NC-FR" in lf.family_label("dojoncfr/dzp")


class TestRankModes:
    def test_orbital_rank_mode_mapping(self):
        assert orbital_rank_mode(None) == "sg15"
        assert orbital_rank_mode("sg15") == "sg15"
        assert orbital_rank_mode("custom") == "sg15"
        assert orbital_rank_mode("apns") == "sg15"
        assert orbital_rank_mode("apns/efficiency") == "apns-efficiency"
        assert orbital_rank_mode("apns/precision") == "apns-precision"
        assert orbital_rank_mode("dojoncfr") == "sg15"
        assert orbital_rank_mode("dojoncfr/sz") == "dojo-sz"
        assert orbital_rank_mode("dojoncfr/dzp") == "dojo-dzp"
        assert orbital_rank_mode("dojoncfr/tzdp") == "dojo-tzdp"

    def test_orbital_rank_mode_rcut_policy(self):
        """libraries.rcut_policy encodes -min/-max onto dojo modes only."""
        pol = lambda p: {"libraries": {"rcut_policy": p}}  # noqa: E731
        assert orbital_rank_mode("dojoncfr/dzp", pol("7")) == "dojo-dzp"
        assert orbital_rank_mode("dojoncfr/dzp") == "dojo-dzp"          # default
        assert orbital_rank_mode("dojoncfr/dzp", pol("min")) == "dojo-dzp-min"
        assert orbital_rank_mode("dojoncfr/dzp", pol("max")) == "dojo-dzp-max"
        assert orbital_rank_mode("dojoncfr/sz", pol("max")) == "dojo-sz-max"
        # bogus policy falls back to the default (no suffix)
        assert orbital_rank_mode("dojoncfr/tzdp", pol("bogus")) == "dojo-tzdp"
        # apns / sg15 ignore the policy entirely.
        assert orbital_rank_mode("apns/efficiency", pol("max")) == "apns-efficiency"
        assert orbital_rank_mode("sg15", pol("min")) == "sg15"
        assert orbital_rank_mode("dojoncfr", pol("min")) == "sg15"      # no tier -> sg15

    def test_precision_prefers_most_complete_basis(self):
        big = _candidate_rank("B_gga_4s4p3d2f_1.2au_120Ry.orb", ".orb", "apns-precision")
        small = _candidate_rank("B_gga_3s3p2d1f_1.2au_120Ry.orb", ".orb", "apns-precision")
        assert big < small

    def test_efficiency_prefers_smaller_rcut(self):
        ten = _candidate_rank("Cs_gga_10au_80Ry_4s2p1d.orb", ".orb", "apns-efficiency")
        twelve = _candidate_rank("Cs_gga_12au_80Ry_4s2p1d.orb", ".orb", "apns-efficiency")
        assert ten < twelve

    def test_precision_resolves_most_complete_file(self, tmp_path):
        d = _mkdir_with(
            tmp_path,
            "B_gga_3s3p2d1f_1.2au_120Ry.orb",
            "B_gga_4s4p3d2f_1.2au_120Ry.orb",
        )
        assert _find_file_for_element(str(d), "B", ".orb", "apns-precision") == \
            "B_gga_4s4p3d2f_1.2au_120Ry.orb"

    def test_efficiency_resolves_cs_10au(self, tmp_path):
        d = _mkdir_with(
            tmp_path,
            "Cs_gga_10au_80Ry_4s2p1d.orb",
            "Cs_gga_12au_80Ry_4s2p1d.orb",
        )
        assert _find_file_for_element(str(d), "Cs", ".orb", "apns-efficiency") == \
            "Cs_gga_10au_80Ry_4s2p1d.orb"

    def test_sg15_mode_unaffected(self, tmp_path):
        """Default sg15 ranking is untouched by the new modes."""
        d = _mkdir_with(tmp_path, "Sm_gga_7au_300.0Ry_4s2p2d1f.orb",
                        "Sm_gga_7au_300.0Ry_2s1p1d.orb")
        assert _find_file_for_element(str(d), "Sm", ".orb") == \
            "Sm_gga_7au_300.0Ry_4s2p2d1f.orb"


class TestDojoTiers:
    """Tier-aware Dojo-NC-FR orbital resolution.

    The Dojo package ships one ``{El}_{SZ,DZP,TZDP}`` folder per element and
    tier (under the same Orbitals_v2.0 root), each holding several rcut copies
    of the same basis.  A ``dojo-<tier>`` rank restricts ``.orb`` matches to the
    requested tier folder — so a missing DZP never silently falls back onto an
    SZ/TZDP file (mixing tiers) — and picks the rcut closest to the canonical
    7 au.  Pseudopotential lookups are never tier-gated (Dojo PPs are a flat
    ``Pseudopotential/`` dir).
    """

    @staticmethod
    def _write_orb(folder: Path, el: str, rcut: int, comp: str) -> None:
        (folder / f"{el}_gga_{rcut}au_100Ry_{comp}.orb").write_text("dummy")

    @staticmethod
    def _orb_root(tmp_path) -> Path:
        return tmp_path / "Orbitals_v2.0"

    def test_dzp_picks_dzp_folder_7au_over_same_file_in_other_tiers(self, tmp_path):
        root = self._orb_root(tmp_path)
        for tier in ("SZ", "DZP", "TZDP"):
            (root / f"Hf_{tier}").mkdir(parents=True)
            self._write_orb(root / f"Hf_{tier}", "Hf", 7, "2s2p1d")
        got = _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-dzp")
        assert got is not None
        assert got.parent.name == "Hf_DZP"
        assert got.name == "Hf_gga_7au_100Ry_2s2p1d.orb"

    def test_dzp_prefers_rcut_closest_to_7_within_tier(self, tmp_path):
        d = self._orb_root(tmp_path) / "Hf_DZP"
        d.mkdir(parents=True)
        for rcut in (6, 9, 10):          # |6-7| is the smallest distance
            self._write_orb(d, "Hf", rcut, "2s2p1d")
        got = _find_file_in_libraries(str(d.parent), "Hf", ".orb", "dojo-dzp")
        assert got.name == "Hf_gga_6au_100Ry_2s2p1d.orb"

    def test_sz_and_tzdp_each_pick_own_tier(self, tmp_path):
        root = self._orb_root(tmp_path)
        for tier in ("SZ", "DZP", "TZDP"):
            (root / f"O_{tier}").mkdir(parents=True)
        for tier, comp in (("SZ", "1s1p"), ("DZP", "2s2p1d"),
                           ("TZDP", "3s3p2d1f")):
            self._write_orb(root / f"O_{tier}", "O", 7, comp)
        assert _find_file_in_libraries(str(root), "O", ".orb",
                                       "dojo-sz").parent.name == "O_SZ"
        assert _find_file_in_libraries(str(root), "O", ".orb",
                                       "dojo-tzdp").parent.name == "O_TZDP"

    def test_dzp_never_leaks_across_tiers(self, tmp_path):
        """SZ/TZDP @7au present, DZP only @8au: dojo-dzp must take the DZP
        file (gated by folder), never the other tiers' closer-to-7 file."""
        root = self._orb_root(tmp_path)
        (root / "Hf_SZ").mkdir(parents=True)
        self._write_orb(root / "Hf_SZ", "Hf", 7, "1s1p")
        (root / "Hf_DZP").mkdir(parents=True)
        self._write_orb(root / "Hf_DZP", "Hf", 8, "2s2p1d")
        (root / "Hf_TZDP").mkdir(parents=True)
        self._write_orb(root / "Hf_TZDP", "Hf", 7, "3s3p2d1f")
        got = _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-dzp")
        assert got.parent.name == "Hf_DZP"
        assert got.name == "Hf_gga_8au_100Ry_2s2p1d.orb"

    def test_missing_tier_returns_none(self, tmp_path):
        """No DZP folder -> dojo-dzp resolves to None (hard-error trigger),
        while the present tiers still resolve under their own modes."""
        root = self._orb_root(tmp_path)
        (root / "Hf_SZ").mkdir(parents=True)
        self._write_orb(root / "Hf_SZ", "Hf", 7, "1s1p")
        (root / "Hf_TZDP").mkdir(parents=True)
        self._write_orb(root / "Hf_TZDP", "Hf", 7, "3s3p2d1f")
        assert _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-dzp") is None
        assert _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-sz") is not None

    def test_rcut_policy_min_max_within_tier(self, tmp_path):
        """dojo-<tier>-min / -max pick the smallest / largest rcut copy of that
        tier — and stay inside it (a smaller/larger SZ/TZDP file must not leak
        in, which would mix tiers)."""
        root = self._orb_root(tmp_path)
        (root / "Hf_SZ").mkdir(parents=True)
        self._write_orb(root / "Hf_SZ", "Hf", 5, "1s1p")         # globally smallest
        (root / "Hf_DZP").mkdir(parents=True)
        for rcut in (6, 7, 8, 9, 10):
            self._write_orb(root / "Hf_DZP", "Hf", rcut, "2s2p1d")
        (root / "Hf_TZDP").mkdir(parents=True)
        self._write_orb(root / "Hf_TZDP", "Hf", 12, "3s3p2d1f")  # globally largest
        got = _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-dzp-min")
        assert got.parent.name == "Hf_DZP"
        assert got.name == "Hf_gga_6au_100Ry_2s2p1d.orb"
        got = _find_file_in_libraries(str(root), "Hf", ".orb", "dojo-dzp-max")
        assert got.parent.name == "Hf_DZP"
        assert got.name == "Hf_gga_10au_100Ry_2s2p1d.orb"
        # End-to-end: config policy -> mode -> resolved file.
        mode = orbital_rank_mode("dojoncfr/dzp",
                                 {"libraries": {"rcut_policy": "max"}})
        assert mode == "dojo-dzp-max"
        got = _find_file_in_libraries(str(root), "Hf", ".orb", mode)
        assert got.name == "Hf_gga_10au_100Ry_2s2p1d.orb"

    def test_pseudo_lookup_never_tier_gated(self, tmp_path):
        """Flat .upf under a plain folder resolves even with a dojo-* rank."""
        root = tmp_path / "Pseudopotential"
        root.mkdir(parents=True)
        (root / "Hf.upf").write_text("dummy")
        assert _find_file_for_element(str(root), "Hf", ".upf", "dojo-dzp") == "Hf.upf"


class TestMaterialize:
    def test_apns_efficiency_sets_lists(self, tmp_path):
        pseudo = tmp_path / "apns-pseudo"
        eff = tmp_path / "apns-eff"
        pseudo.mkdir()
        eff.mkdir()
        cfg = _apns_config(str(pseudo), str(eff))
        assert lf.materialize_family(cfg, "apns", "efficiency") is True
        assert cfg["libraries"]["family"] == "apns/efficiency"
        assert cfg["libraries"]["pseudo_library"] == [str(pseudo.resolve())]
        assert cfg["libraries"]["orbital_library"] == [str(eff.resolve())]

    def test_dojo_dzp_sets_lists(self, tmp_path):
        pseudo = tmp_path / "Pseudopotential"
        orb = tmp_path / "Orbitals_v2.0"
        pseudo.mkdir()
        orb.mkdir()
        cfg = _dojo_config(str(pseudo), str(orb))
        assert lf.materialize_family(cfg, "dojoncfr", "dzp") is True
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        assert cfg["libraries"]["pseudo_library"] == [str(pseudo.resolve())]
        assert cfg["libraries"]["orbital_library"] == [str(orb.resolve())]

    def test_dojo_sz_sets_lists(self, tmp_path):
        pseudo = tmp_path / "Pseudopotential"
        orb = tmp_path / "Orbitals_v2.0"
        pseudo.mkdir()
        orb.mkdir()
        cfg = _dojo_config(str(pseudo), str(orb))
        assert lf.materialize_family(cfg, "dojoncfr", "sz") is True
        assert cfg["libraries"]["family"] == "dojoncfr/sz"
        assert cfg["libraries"]["orbital_library"] == [str(orb.resolve())]

    def test_dojo_pseudo_only_missing_returns_false(self, tmp_path):
        orb = tmp_path / "Orbitals_v2.0"
        orb.mkdir()
        cfg = _dojo_config(str(tmp_path / "no-pseudo"), str(orb))
        assert lf.materialize_family(cfg, "dojoncfr", "dzp") is False
        assert cfg["libraries"]["family"] == "sg15"

    def test_sg15_reverts_to_auto_detection(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "abacuscopilot.config._detect_library_dirs",
            lambda ext, exclude_containing=(): [f"/x/{ext}"],
        )
        cfg = {"libraries": {"pseudo_library": ["stale"], "orbital_library": ["stale"],
                             "family": "apns/precision", "families": {}}}
        assert lf.materialize_family(cfg, "sg15") is True
        assert cfg["libraries"]["family"] == "sg15"
        assert cfg["libraries"]["pseudo_library"] == ["/x/.upf"]
        assert cfg["libraries"]["orbital_library"] == ["/x/.orb"]

    def test_missing_path_returns_false_and_leaves_config(self, tmp_path):
        pseudo = tmp_path / "apns-pseudo"
        pseudo.mkdir()
        cfg = _apns_config(str(pseudo), str(tmp_path / "does-not-exist"))
        cfg["libraries"]["family"] = "apns/precision"
        before = dict(cfg["libraries"])
        assert lf.materialize_family(cfg, "apns", "precision") is False
        assert cfg["libraries"]["family"] == "apns/precision"
        assert cfg["libraries"] == before

    def test_unknown_family_returns_false(self):
        assert lf.materialize_family({"libraries": {}}, "nope") is False


class TestSg15ExcludesRegisteredFamilies:
    """A registered external family (e.g. an ABACUS-APNS copy placed under
    PP-Orb/) must never be absorbed into the auto-detected SG15 lists, or the
    two series would silently mix under the default SG15 family."""

    def _cfg_with_apns(self, root: Path) -> dict:
        pseudo = root / "apns-pseudopotentials-v1"
        eff = root / "apns-orbitals-efficiency-v1"
        pre = root / "apns-orbitals-precision-v1"
        for d in (pseudo, eff, pre):
            d.mkdir(parents=True)
        return {
            "libraries": {
                "pseudo_library": [],
                "orbital_library": [],
                "family": "apns/efficiency",
                "families": {
                    "apns": {
                        "pseudo_dir": str(pseudo),
                        "orbital_dirs": {
                            "efficiency": str(eff),
                            "precision": str(pre),
                        },
                    }
                },
            }
        }

    def test_path_within(self, tmp_path):
        root = tmp_path / "lib"
        assert _path_within(root, root) is True           # equal counts
        assert _path_within(root / "a" / "b", root) is True
        assert _path_within(root, root / "a") is False    # sibling/child of? no
        assert _path_within(tmp_path, root) is False      # outside

    def test_registered_family_dirs_collected(self, tmp_path):
        cfg = self._cfg_with_apns(tmp_path)
        got = sorted(_registered_family_dirs(cfg))
        assert got == sorted([
            Path(str(tmp_path / "apns-pseudopotentials-v1")),
            Path(str(tmp_path / "apns-orbitals-efficiency-v1")),
            Path(str(tmp_path / "apns-orbitals-precision-v1")),
        ])
        assert _registered_family_dirs({"libraries": {}}) == []
        assert _registered_family_dirs({}) == []

    def test_materialize_sg15_passes_registered_dirs_as_exclusion(
            self, tmp_path, monkeypatch):
        cfg = self._cfg_with_apns(tmp_path)
        captured: list = []

        def fake_detect(ext, exclude_containing=()):
            captured.append((ext, set(exclude_containing)))
            return [f"/bundled/{ext}"]

        monkeypatch.setattr(
            "abacuscopilot.config._detect_library_dirs", fake_detect)
        assert lf.materialize_family(cfg, "sg15") is True
        assert cfg["libraries"]["family"] == "sg15"
        assert cfg["libraries"]["pseudo_library"] == ["/bundled/.upf"]
        assert cfg["libraries"]["orbital_library"] == ["/bundled/.orb"]
        # Both extensions asked with the SAME exclusion set (all registered dirs).
        assert len(captured) == 2
        expected = {
            Path(str(tmp_path / d))
            for d in ("apns-pseudopotentials-v1", "apns-orbitals-efficiency-v1",
                      "apns-orbitals-precision-v1")
        }
        for ext, excl in captured:
            assert excl == expected


class TestMissingElements:
    def test_missing_library_files_reports_uncovered(self, tmp_path):
        d = _mkdir_with(tmp_path, "Si_ONCV_PBE-1.0.upf")
        species = ["Si", "Os"]
        missing = missing_library_files(species, "pw", str(d), "")
        assert ("pseudopotential", "Os") in missing
        assert ("pseudopotential", "Si") not in missing

    def test_cwd_file_covers_a_missing_element(self, tmp_path):
        lib = _mkdir_with(tmp_path, "Si_ONCV_PBE-1.0.upf")
        empty = tmp_path / "work"
        empty.mkdir()
        # Library lacks Os, and cwd has none -> blocker.
        missing = missing_library_files(["Si", "Os"], "pw", str(lib), "")
        assert missing_element_blockers(missing, cwd=empty) == [
            ("pseudopotential", "Os")
        ]
        # User-provided Os.upf in cwd removes the blocker.
        (empty / "Os.upf").write_text("dummy")
        assert missing_element_blockers(missing, cwd=empty) == []


class TestPicker:
    def test_noninteractive_noop(self, tmp_path):
        cfg = _apns_config(str(tmp_path), str(tmp_path))
        console = FakeConsole()
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=False) is False

    def _empty_config(self) -> dict:
        return {"libraries": {"pseudo_library": [], "orbital_library": [],
                              "family": "sg15", "families": {}}}

    def test_empty_config_prints_no_library_hint(self, tmp_path, monkeypatch):
        """Nothing usable under PP-Orb/ and no registered series -> guide the
        user to place their own *.upf / *.orb libraries there."""
        cfg = self._empty_config()
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        console = FakeConsole(answers=[""])  # keep current (SG15)
        assert lf.pick_library_family(console, cfg, basis_type="pw",
                                      interactive=True) is False
        assert any("No pseudopotential/orbital libraries" in line
                   for line in console.out)
        assert any("PP-Orb" in line for line in console.out)

    def test_registered_family_suppresses_hint(self, tmp_path, monkeypatch):
        """Empty lists but an external series registered -> the picker can
        activate it, so the PP-Orb hint is not needed."""
        pseudo = tmp_path / "p"
        pseudo.mkdir()
        cfg = _apns_config(str(pseudo), str(pseudo))
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        console = FakeConsole(answers=[""])  # keep current (SG15)
        assert lf.pick_library_family(console, cfg, basis_type="pw",
                                      interactive=True) is False
        assert not any("No pseudopotential/orbital libraries" in line
                       for line in console.out)

    def test_blank_keeps_current(self, tmp_path, monkeypatch):
        pseudo = tmp_path / "p"
        pseudo.mkdir()
        cfg = _apns_config(str(pseudo), str(pseudo))
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        console = FakeConsole(answers=[""])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is False
        assert saved == []

    def test_switch_to_apns_efficiency(self, tmp_path, monkeypatch):
        pseudo = tmp_path / "p"
        eff = tmp_path / "e"
        pseudo.mkdir()
        eff.mkdir()
        cfg = _apns_config(str(pseudo), str(eff))
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        console = FakeConsole(answers=["2", "1"])  # APNS, then efficiency
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "apns/efficiency"
        assert cfg["libraries"]["pseudo_library"] == [str(pseudo.resolve())]
        assert cfg["libraries"]["orbital_library"] == [str(eff.resolve())]
        assert saved == [cfg]

    def test_switch_to_dojoncfr_dzp(self, tmp_path, monkeypatch):
        """Picking Dojo-NC-FR with no paths configured prompts for pseudo +
        orbital root once, then materializes dojoncfr/dzp and prints the SOC note."""
        pseudo = tmp_path / "Pseudopotential"
        orb = tmp_path / "Orbitals_v2.0"
        pseudo.mkdir()
        orb.mkdir()
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # No PP-Orb drop-in to auto-detect (dev PP-Orb physically holds Dojo —
        # point detection elsewhere so the manual path prompts are exercised).
        monkeypatch.setattr(lf, "_pp_orb_root", lambda: tmp_path / "no-pp-orb")
        # options: 1=SG15, 2=apns, 3=dojoncfr; then 1=dzp variant.
        # pseudo_dir then dzp orbital_dir both get prompted (empty family).
        console = FakeConsole(answers=["3", "1", str(pseudo), str(orb)])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        fam = cfg["libraries"]["families"]["dojoncfr"]
        assert fam["pseudo_dir"] == str(pseudo.resolve())
        assert fam["orbital_dirs"]["dzp"] == str(orb.resolve())
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        assert cfg["libraries"]["pseudo_library"] == [str(pseudo.resolve())]
        assert cfg["libraries"]["orbital_library"] == [str(orb.resolve())]
        assert saved == [cfg]
        # SOC context is spelled out so the user knows what this family is for.
        assert any("spin-orbit coupling" in line for line in console.out)

    def _dojo_pp_orb(self, tmp_path) -> Path:
        """Build a fake PP-Orb root holding the Dojo drop-in folder."""
        root = tmp_path / "pp-orb"
        (root / "Dojo-NC-FR" / "Pseudopotential").mkdir(parents=True, exist_ok=True)
        (root / "Dojo-NC-FR" / "Orbitals_v2.0").mkdir(parents=True, exist_ok=True)
        return root

    def test_dojo_auto_detected_from_pp_orb(self, tmp_path, monkeypatch):
        """Picking Dojo-NC-FR when its folder already sits under PP-Orb/
        auto-registers the paths — no manual prompt needed."""
        monkeypatch.setattr(lf, "_pp_orb_root",
                            lambda: self._dojo_pp_orb(tmp_path))
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # Only the family + variant (+ rcut) prompts: 1=SG15 2=apns 3=dojoncfr.
        console = FakeConsole(answers=["3", "1"])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        root = tmp_path / "pp-orb" / "Dojo-NC-FR"
        fam = cfg["libraries"]["families"]["dojoncfr"]
        assert fam["pseudo_dir"] == str((root / "Pseudopotential").resolve())
        # All three tiers share the single Orbitals_v2.0 root.
        for tier in ("sz", "dzp", "tzdp"):
            assert fam["orbital_dirs"][tier] == str((root / "Orbitals_v2.0").resolve())
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        # No path prompt was ever shown.
        assert not any("Path to dojoncfr" in line for line in console.out)
        assert any("auto-registered" in line for line in console.out)

    def test_dojo_auto_detected_partial_falls_back_to_prompt(self, tmp_path, monkeypatch):
        """Auto-detect only fills what PP-Orb provides; a missing piece (here:
        a variant whose orbital folder is absent) is still asked for."""
        root = tmp_path / "pp-orb"
        (root / "Dojo-NC-FR" / "Pseudopotential").mkdir(parents=True)
        monkeypatch.setattr(lf, "_pp_orb_root", lambda: root)
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # family -> variant -> prompt for the (missing) orbital root dir.
        orb = tmp_path / "my-orbs"
        orb.mkdir()
        console = FakeConsole(answers=["3", "1", str(orb)])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        fam = cfg["libraries"]["families"]["dojoncfr"]
        assert fam["pseudo_dir"] == str((root / "Dojo-NC-FR" / "Pseudopotential").resolve())
        assert fam["orbital_dirs"]["dzp"] == str(orb.resolve())
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"

    def test_apns_auto_detected_from_pp_orb(self, tmp_path, monkeypatch):
        """Same drop-in auto-registration for the APNS series layout."""
        root = tmp_path / "pp-orb"
        top = root / "ABACUS-APNS-PPORBs-v1"
        (top / "apns-pseudopotentials-v1").mkdir(parents=True)
        (top / "apns-orbitals-efficiency-v1").mkdir(parents=True)
        (top / "apns-orbitals-precision-v1").mkdir(parents=True)
        monkeypatch.setattr(lf, "_pp_orb_root", lambda: root)
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # family (2=apns) -> efficiency variant.
        console = FakeConsole(answers=["2", "1"])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        fam = cfg["libraries"]["families"]["apns"]
        assert fam["pseudo_dir"] == str((top / "apns-pseudopotentials-v1").resolve())
        assert (fam["orbital_dirs"]["efficiency"]
                == str((top / "apns-orbitals-efficiency-v1").resolve()))
        assert (fam["orbital_dirs"]["precision"]
                == str((top / "apns-orbitals-precision-v1").resolve()))
        assert cfg["libraries"]["family"] == "apns/efficiency"
        assert not any("Path to apns" in line for line in console.out)

    def test_option_text_shows_auto_register_when_dropped_in(self, tmp_path, monkeypatch):
        """The picker menu says a present drop-in will auto-register instead of
        claiming it will prompt."""
        monkeypatch.setattr(lf, "_pp_orb_root",
                            lambda: self._dojo_pp_orb(tmp_path))
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        console = FakeConsole(answers=[""])  # keep current; just render menu
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        lf.pick_library_family(console, cfg, basis_type="pw", interactive=True)
        text = "\n".join(console.out)
        assert "Dojo-NC-FR" in text
        assert "found in PP-Orb" in text
        # The Dojo option specifically no longer claims it will prompt.
        dojo_line = next(ln for ln in console.out if "Dojo-NC-FR" in ln)
        assert "found in PP-Orb — will auto-register" in dojo_line

    def test_reselect_dojo_blank_keeps_current_tier(self, tmp_path, monkeypatch):
        """Re-picking active dojoncfr re-asks the tier; Enter keeps dzp (never
        silently drops to a smaller tier)."""
        pseudo = tmp_path / "p"
        orb = tmp_path / "o"
        pseudo.mkdir()
        orb.mkdir()
        cfg = _dojo_config(str(pseudo), str(orb), family="dojoncfr/dzp")
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        console = FakeConsole(answers=["3", ""])  # dojoncfr, blank tier prompt
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        assert any("(current)" in line for line in console.out)

    def test_dojo_pick_sets_rcut_policy(self, tmp_path, monkeypatch):
        """Picking Dojo asks which rcut copy to use; choosing one persists it as
        libraries.rcut_policy for every later resolution."""
        pseudo = tmp_path / "Pseudopotential"
        orb = tmp_path / "Orbitals_v2.0"
        pseudo.mkdir()
        orb.mkdir()
        cfg = _dojo_config(str(pseudo), str(orb))
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # dojoncfr -> dzp -> rcut policy option 2 (smallest).
        console = FakeConsole(answers=["3", "1", "2"])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        assert cfg["libraries"]["rcut_policy"] == "min"
        assert saved == [cfg]

    def test_dojo_reselect_blank_keeps_rcut_policy(self, tmp_path, monkeypatch):
        """Re-picking the active Dojo re-asks the rcut copy, but a blank Enter
        keeps the persisted policy (never silently re-defaults to 7 au)."""
        pseudo = tmp_path / "Pseudopotential"
        orb = tmp_path / "Orbitals_v2.0"
        pseudo.mkdir()
        orb.mkdir()
        cfg = _dojo_config(str(pseudo), str(orb), family="dojoncfr/dzp")
        cfg["libraries"]["rcut_policy"] = "max"
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        # dojoncfr, blank tier (=dzp), blank rcut (=keep max).
        console = FakeConsole(answers=["3", "", ""])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "dojoncfr/dzp"
        assert cfg["libraries"]["rcut_policy"] == "max"
        assert any("(current)" in line for line in console.out)

    def test_dojo_path_prompt_abort_leaves_config(self, tmp_path, monkeypatch):
        """Aborting a path prompt must not mutate the config."""
        orb = tmp_path / "o"
        orb.mkdir()
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        # No PP-Orb drop-in to auto-detect, so the prompt actually fires.
        monkeypatch.setattr(lf, "_pp_orb_root", lambda: tmp_path / "no-pp-orb")
        console = FakeConsole(answers=["3", "1", ""])  # abort on pseudo dir
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is False
        assert cfg["libraries"]["family"] == "sg15"
        # Aborting mid-prompt leaves the family stub without any stored path.
        assert "pseudo_dir" not in cfg["libraries"]["families"]["dojoncfr"]

    def _cfg_on_precision(self, tmp_path):
        pseudo = tmp_path / "p"
        eff = tmp_path / "e"
        pre = tmp_path / "r"
        for d in (pseudo, eff, pre):
            d.mkdir()
        cfg = _apns_config(str(pseudo), str(eff))
        cfg["libraries"]["families"]["apns"]["orbital_dirs"]["precision"] = str(pre)
        cfg["libraries"]["family"] = "apns/precision"
        return cfg

    def test_reselect_blank_keeps_current_variant(self, tmp_path, monkeypatch):
        """Re-picking the active APNS re-asks the variant, but Enter keeps the
        current one — precision is never silently downgraded to efficiency."""
        cfg = self._cfg_on_precision(tmp_path)
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        # "2" = APNS, then "" (blank) on the variant prompt = keep current.
        console = FakeConsole(answers=["2", ""])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "apns/precision"
        assert any("(current)" in line for line in console.out)

    def test_reselect_can_switch_variant_explicitly(self, tmp_path, monkeypatch):
        """The variant ask stays reachable: picking 1 flips precision→efficiency."""
        cfg = self._cfg_on_precision(tmp_path)
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        console = FakeConsole(answers=["2", "1"])  # APNS, then efficiency
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "apns/efficiency"
        eff = str((tmp_path / "e").resolve())
        assert cfg["libraries"]["orbital_library"] == [eff]

    def test_prompts_for_missing_paths_first_time(self, tmp_path, monkeypatch):
        """Unconfigured APNS is offered; picking it guides path entry."""
        pseudo = tmp_path / "apns-pseudo"
        eff = tmp_path / "apns-eff"
        pseudo.mkdir()
        eff.mkdir()
        cfg = {"libraries": {"pseudo_library": [], "orbital_library": [],
                             "family": "sg15", "families": {}}}
        saved = []
        monkeypatch.setattr(lf, "save_config", lambda c: saved.append(c))
        # No PP-Orb drop-in to auto-detect (dev PP-Orb physically holds APNS —
        # point detection elsewhere so the manual path prompts are exercised).
        monkeypatch.setattr(lf, "_pp_orb_root", lambda: tmp_path / "no-pp-orb")
        console = FakeConsole(answers=["2", "1", str(pseudo), str(eff)])
        assert lf.pick_library_family(console, cfg, basis_type="lcao",
                                      interactive=True) is True
        fam = cfg["libraries"]["families"]["apns"]
        assert fam["pseudo_dir"] == str(pseudo.resolve())
        assert fam["orbital_dirs"]["efficiency"] == str(eff.resolve())
        assert cfg["libraries"]["family"] == "apns/efficiency"
        assert cfg["libraries"]["pseudo_library"] == [str(pseudo.resolve())]
        assert cfg["libraries"]["orbital_library"] == [str(eff.resolve())]

    def test_switch_back_to_sg15(self, tmp_path, monkeypatch):
        pseudo = tmp_path / "p"
        pseudo.mkdir()
        cfg = _apns_config(str(pseudo), str(pseudo))
        cfg["libraries"]["family"] = "apns/precision"
        monkeypatch.setattr(
            "abacuscopilot.config._detect_library_dirs",
            lambda ext, exclude_containing=(): [],
        )
        # pick_library_family persists with save_config — keep it out of the
        # developer's real ~/.abacuscopilot/config.yaml (regression: an earlier
        # version of this test leaked its tmp_path into the live config).
        monkeypatch.setattr(lf, "save_config", lambda c: None)
        console = FakeConsole(answers=["1"])  # SG15
        assert lf.pick_library_family(console, cfg, basis_type="pw",
                                      interactive=True) is True
        assert cfg["libraries"]["family"] == "sg15"
        assert cfg["libraries"]["pseudo_library"] == []


class TestAutoPrepareHardError:
    """Missing element under a user family (APNS) → hard error, no mixing;
    under the bundled SG15 family the legacy soft-warning path is unchanged."""

    def _monkeypatch_flow(self, tmp_path, monkeypatch, cfg):
        import abacuscopilot.preprocessing.input_tasks as it

        monkeypatch.chdir(tmp_path)
        # Keep the family picker silent; the flow uses the persisted family.
        monkeypatch.setattr("abacuscopilot.library_families.pick_library_family",
                            lambda *a, **k: False)
        # Feed species directly (no STRU needed for this code path).
        monkeypatch.setattr(
            "abacuscopilot.preprocessing.system_tasks.read_species_from_stru",
            lambda path: ["Si", "Os"],
        )
        monkeypatch.setattr(it, "load_config", lambda: cfg)
        from abacuscopilot.core.exceptions import LibraryFamilyError
        from abacuscopilot.core.models import InputParams

        params = InputParams()
        params.basis_type = "lcao"
        return it, params, LibraryFamilyError

    def test_user_family_missing_element_raises(self, tmp_path, monkeypatch):
        from abacuscopilot.preprocessing.input_tasks import _auto_prepare_files

        pseudo = tmp_path / "p"
        eff = tmp_path / "e"
        pseudo.mkdir()
        eff.mkdir()
        (pseudo / "Si_ONCV_PBE-1.0.upf").write_text("dummy")
        (eff / "Si_gga_7au_100.0Ry_4s2p2d1f.orb").write_text("dummy")
        cfg = {"libraries": {"pseudo_library": [str(pseudo)],
                             "orbital_library": [str(eff)],
                             "family": "apns/efficiency",
                             "families": {}}}
        it, params, LibraryFamilyError = self._monkeypatch_flow(
            tmp_path, monkeypatch, cfg)
        console = FakeConsole()
        with pytest.raises(LibraryFamilyError) as ei:
            _auto_prepare_files(console, params, interactive=True)
        assert "Os" in str(ei.value)
        assert "ABACUS-APNS-PPORBs-v1" in str(ei.value)

    def test_sg15_missing_element_stays_soft(self, tmp_path, monkeypatch):
        from abacuscopilot.preprocessing.input_tasks import _auto_prepare_files

        pseudo = tmp_path / "p"
        pseudo.mkdir()
        (pseudo / "Si_ONCV_PBE-1.0.upf").write_text("dummy")
        cfg = {"libraries": {"pseudo_library": [str(pseudo)],
                             "orbital_library": [],
                             "family": "sg15",
                             "families": {}}}
        it, params, LibraryFamilyError = self._monkeypatch_flow(
            tmp_path, monkeypatch, cfg)
        console = FakeConsole()
        # Must NOT raise; legacy soft-warn path (Os missing is only a yellow note).
        _auto_prepare_files(console, params, interactive=True)

    def test_dojo_missing_tier_element_raises(self, tmp_path, monkeypatch):
        """Under dojoncfr/dzp an element with no DZP tier is a hard error
        naming the Dojo family — never a silent SZ/TZDP substitution."""
        from abacuscopilot.preprocessing.input_tasks import _auto_prepare_files

        pseudo = tmp_path / "Pseudopotential"
        pseudo.mkdir()
        (pseudo / "Si.upf").write_text("dummy")
        orb = tmp_path / "Orbitals_v2.0"
        (orb / "Si_DZP").mkdir(parents=True)
        (orb / "Si_DZP" / "Si_gga_7au_100Ry_2s2p1d.orb").write_text("dummy")
        cfg = {"libraries": {"pseudo_library": [str(pseudo)],
                             "orbital_library": [str(orb)],
                             "family": "dojoncfr/dzp",
                             "families": {}}}
        it, params, LibraryFamilyError = self._monkeypatch_flow(
            tmp_path, monkeypatch, cfg)
        console = FakeConsole()
        with pytest.raises(LibraryFamilyError) as ei:
            _auto_prepare_files(console, params, interactive=True)
        assert "Os" in str(ei.value)
        assert "Dojo-NC-FR" in str(ei.value)
