"""Tests for PP/orbital library resolution: multi-dir, recursive, tolerant naming."""


from abacuscopilot.config import _valid_library_dirs
from abacuscopilot.preprocessing.system_tasks import (
    _candidate_rank,
    _find_file_for_element,
    _find_file_in_libraries,
)


class TestMatcherNaming:
    def test_lanthanide_upf(self, tmp_path):
        """'Sm3+_f--core-icmod1.PD04.PBE.UPF' must resolve for element Sm."""
        (tmp_path / "Sm3+_f--core-icmod1.PD04.PBE.UPF").write_text("dummy")
        assert _find_file_for_element(str(tmp_path), "Sm", ".upf") == \
            "Sm3+_f--core-icmod1.PD04.PBE.UPF"

    def test_plain_upf(self, tmp_path):
        """'Ag.upf' must resolve for element Ag."""
        (tmp_path / "Ag.upf").write_text("dummy")
        assert _find_file_for_element(str(tmp_path), "Ag", ".upf") == "Ag.upf"

    def test_dot_family_upf(self, tmp_path):
        """'As.PD04.PBE.UPF' (dot-prefix + uppercase) must resolve for As."""
        (tmp_path / "As.PD04.PBE.UPF").write_text("dummy")
        assert _find_file_for_element(str(tmp_path), "As", ".upf") == "As.PD04.PBE.UPF"

    def test_sulfur_not_samarium(self, tmp_path):
        """Element 'S' must NOT match a Sm3+ file."""
        (tmp_path / "Sm3+_f--core-icmod1.PD04.PBE.UPF").write_text("dummy")
        (tmp_path / "S_ONCV_PBE-1.0.upf").write_text("dummy")
        assert _find_file_for_element(str(tmp_path), "S", ".upf") == "S_ONCV_PBE-1.0.upf"

    def test_missing_returns_none(self, tmp_path):
        assert _find_file_for_element(str(tmp_path), "Sm", ".upf") is None


class TestMatcherMultipleDirs:
    def test_searches_dir_list_in_order(self, tmp_path):
        d1 = tmp_path / "sg15"
        d2 = tmp_path / "lanthanides"
        d1.mkdir()
        d2.mkdir()
        (d1 / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        # Sm only in the second dir
        (d2 / "Sm3+_f--core-icmod1.PD04.PBE.UPF").write_text("dummy")
        dirs = [str(d1), str(d2)]
        assert _find_file_for_element(dirs, "Li", ".upf") == "Li_ONCV_PBE-1.0.upf"
        assert _find_file_for_element(dirs, "Sm", ".upf") == \
            "Sm3+_f--core-icmod1.PD04.PBE.UPF"
        assert _find_file_in_libraries(dirs, "Sm", ".upf").parent == d2

    def test_recursive_nested_lanthanide_layout(self, tmp_path):
        """Files nested under element/basis folders must be found."""
        el = tmp_path / "Sm" / "Sm_4s2p2d1f"
        el.mkdir(parents=True)
        (el / "Sm_gga_7au_300.0Ry_4s2p2d1f.orb").write_text("dummy")
        assert _find_file_for_element(str(tmp_path), "Sm", ".orb") == \
            "Sm_gga_7au_300.0Ry_4s2p2d1f.orb"


class TestMatcherOrbitalPreference:
    def test_prefers_dzp_at_7au(self, tmp_path):
        """Among 15 lanthanide orbitals, pick 4s2p2d1f at 7 au."""
        el = tmp_path / "Sm"
        el.mkdir()
        made = []
        for basis in ("2s1p1d", "4s2p2d1f", "6s3p3d2f"):
            for rcut in (6, 7, 8, 9, 10):
                name = f"Sm_gga_{rcut}au_300.0Ry_{basis}.orb"
                (el / name).write_text("dummy")
                made.append(name)
        assert _find_file_for_element(str(tmp_path), "Sm", ".orb") == \
            "Sm_gga_7au_300.0Ry_4s2p2d1f.orb"

    def test_candidate_rank_ordering(self):
        low = _candidate_rank("Sm_gga_7au_300.0Ry_4s2p2d1f.orb", ".orb")
        other = _candidate_rank("Sm_gga_7au_300.0Ry_2s1p1d.orb", ".orb")
        far = _candidate_rank("Sm_gga_10au_300.0Ry_4s2p2d1f.orb", ".orb")
        assert low < other
        assert low < far


class TestConfigValidation:
    def test_valid_library_dirs_filters_stale(self, tmp_path):
        real = tmp_path / "lib"
        real.mkdir()
        (real / "Li_ONCV_PBE-1.0.upf").write_text("dummy")
        stale = tmp_path / "gone"
        assert _valid_library_dirs(str(real), ".upf") == [str(real.resolve())]
        assert _valid_library_dirs([str(real), str(stale)], ".upf") == [str(real.resolve())]
        assert _valid_library_dirs([str(stale)], ".upf") == []
        assert _valid_library_dirs("", ".upf") == []
