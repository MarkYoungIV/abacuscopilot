"""What reaches ~/.abacuscopilot/config.yaml, and in what shape."""

from pathlib import Path

import yaml

from abacuscopilot import config as cfg


def _use_tmp_home(tmp_path, monkeypatch):
    monkeypatch.setattr(cfg, "_get_config_dir", lambda: tmp_path)


def _rendered(config):
    return [line for line in cfg.render_config(config).splitlines() if line.strip()]


def _key_lines(config):
    """The lines holding a ``key:`` of some section.

    Continuation lines start deeper than the key indent (a block sequence's
    ``- item``, a nested mapping's children) and take no note of their own.
    """
    return [
        line for line in _rendered(config)
        if line.startswith(cfg._INDENT) and not line[len(cfg._INDENT):].startswith(" ")
    ]


class TestNoteLayout:
    """Notes live beside their value, in the style of VASPKIT's ~/.vaspkit.

    Not above it.  A revision that put explanations on their own lines roughly
    doubled the file: comments became the majority of it, and the values -- the
    only reason to open the file -- were what got buried.
    """

    def test_a_note_never_takes_a_line_of_its_own(self):
        """Every comment is a trailing one, attached to a value.

        Section headings are the sole exception, and they are recognised by
        sitting at column 0 directly above the section they introduce.
        """
        lines = _rendered(cfg.DEFAULT_CONFIG)
        standalone = [line for line in lines
                      if line.lstrip().startswith("#") and not line.startswith("#")]
        assert standalone == []

    def test_every_key_carries_its_note_on_the_same_line(self):
        uncommented = [line for line in _key_lines(cfg.DEFAULT_CONFIG) if "#" not in line]
        assert uncommented == []

    def test_notes_line_up_in_a_column(self):
        """A ragged right margin is what makes a file like this unskimmable."""
        columns = {line.index("#") for line in _key_lines(cfg.DEFAULT_CONFIG) if "#" in line}
        assert len(columns) == 1, f"notes at columns {sorted(columns)}"

    def test_comments_are_a_minority_of_the_file(self):
        """The regression guard for the original complaint, in one number.

        The layout this replaced put a paragraph above every key: comments were
        then most of the file (40-odd comment-only lines against 36 keys), and
        the values were what you had to hunt for.  Written BESIDE the value, a
        note costs no line at all -- only the banner and the four section
        headings stand alone.
        """
        lines = _rendered(cfg.DEFAULT_CONFIG)
        comment_only = [line for line in lines if line.lstrip().startswith("#")]
        assert len(comment_only) < sum(len(v) for v in cfg.DEFAULT_CONFIG.values()) / 2


class TestRenderConfig:
    def test_every_default_key_is_documented(self):
        """A new key must not be able to slip in uncommented."""
        undocumented = [
            f"{section}.{key}"
            for section, values in cfg.DEFAULT_CONFIG.items()
            if isinstance(values, dict)
            for key in values
            if f"{section}.{key}" not in cfg._FIELD_HELP
        ]
        assert undocumented == []

    def test_no_note_is_left_for_a_key_that_no_longer_exists(self):
        """A stale entry is a note the user reads and then cannot find."""
        stale = [
            path for path in cfg._FIELD_HELP
            if path.split(".", 1)[-1] not in cfg.DEFAULT_CONFIG.get(path.split(".")[0], {})
        ]
        assert stale == []

    def test_every_section_is_introduced(self):
        assert set(cfg.DEFAULT_CONFIG) <= set(cfg._SECTION_HELP)

    def test_round_trips_the_defaults(self):
        assert yaml.safe_load(cfg.render_config(cfg.DEFAULT_CONFIG)) == cfg.DEFAULT_CONFIG

    def test_round_trips_arbitrary_nested_values(self):
        """libraries.families holds user-defined trees of any shape."""
        config = {
            "libraries": {
                "families": {"apns": {"pseudo_dir": "/x/PP", "orbital_dirs": {"d": "/x/O"}}},
                "family": "apns",
            }
        }
        assert yaml.safe_load(cfg.render_config(config)) == config

    def test_round_trips_an_empty_section(self):
        """A bare ``paths:`` would come back as None, not {}."""
        assert yaml.safe_load(cfg.render_config({"paths": {}})) == {"paths": {}}

    def test_quotes_a_version_string_that_looks_numeric(self):
        """'2.10' unquoted loads as the float 2.1 — a silent wrong version."""
        text = cfg.render_config({"paths": {"abacus_dp_version": "2.10"}})
        assert yaml.safe_load(text)["paths"]["abacus_dp_version"] == "2.10"

    def test_keeps_a_key_the_user_added(self):
        config = {"my_thing": {"a": 1}}
        assert yaml.safe_load(cfg.render_config(config)) == config


class TestSaveConfig:
    def test_writes_what_render_config_produces(self, tmp_path, monkeypatch):
        _use_tmp_home(tmp_path, monkeypatch)
        cfg.save_config(cfg.DEFAULT_CONFIG)
        assert (tmp_path / "config.yaml").read_text() == cfg.render_config(cfg.DEFAULT_CONFIG)


class TestLoadConfig:
    def test_user_values_survive_a_load(self, tmp_path, monkeypatch):
        """Reading is not writing: nothing the user typed may come back changed."""
        _use_tmp_home(tmp_path, monkeypatch)
        path = tmp_path / "config.yaml"
        path.write_text(yaml.safe_dump({"paths": {"abacus_binary": "/opt/abacus"}}))
        before = path.read_text()

        merged = cfg.load_config()

        assert merged["paths"]["abacus_binary"] == "/opt/abacus"
        assert merged["paths"]["abacus_dp_binary"] == ""  # filled in memory only
        assert path.read_text() == before

    def test_a_note_the_user_left_in_the_file_is_left_alone(self, tmp_path, monkeypatch):
        """Their edits and their comments are theirs; loading never rewrites."""
        _use_tmp_home(tmp_path, monkeypatch)
        path = tmp_path / "config.yaml"
        path.write_text("# my own note\npaths:\n  abacus_binary: /opt/abacus\n")
        before = path.read_text()

        cfg.load_config()

        assert path.read_text() == before


class TestReleaseRoot:
    """The data root — PP-Orb/ and scripts/ live there, not next to the package.

    These pin the src layout deliberately.  Every one of them passed under the
    old flat layout too, which is the point: they describe the *contract*
    (release_root() is the tree holding the data dirs), not a hop count, so they
    would have caught PP-Orb silently resolving to <root>/src.
    """

    def test_is_the_tree_holding_pyproject(self):
        root = cfg.release_root()
        assert (root / "pyproject.toml").is_file()

    def test_is_not_the_src_dir(self):
        """The failure mode: counting .parent hops lands one level short."""
        assert cfg.release_root() != Path(cfg.__file__).resolve().parent.parent

    def test_is_an_ancestor_of_the_package(self):
        Path(cfg.__file__).resolve().parent.relative_to(cfg.release_root())

    def test_the_data_dirs_hang_off_it(self):
        """scripts/ is tracked in git; PP-Orb/ may be empty in a bare clone."""
        root = cfg.release_root()
        assert (root / "scripts").is_dir()


def _build_pp_orb(tmp_path):
    """A PP-Orb/ tree in the shipped shape: one top-level folder per series.

    Carries both bundled series plus the two external ones users drop in, since
    it is exactly the *coexistence* of those that has to be handled correctly.
    """
    pp = tmp_path / "PP-Orb"
    sg15_pp = pp / "SG15-Version1p0" / "SG15-Version1p0_Pseudopotential"
    sg15_orb = pp / "SG15-Version1p0" / "SG15-Version1p0__StandardOrbitals-Version2p0"
    lan = pp / "lanthanides-f--core.icmod1" / "PD04.3+f--core.icmod1"
    _write(sg15_pp / "Li_ONCV_PBE-1.0.upf")
    _write(sg15_orb / "Li_gga_7au_100Ry_4s1p.orb")
    _write(lan / "Sm" / "Sm3+_f--core-icmod1.PD04.PBE.UPF")
    _write(lan / "Sm" / "Sm_gga_7au_300.0Ry_4s2p2d1f.orb")
    # The external series, physically present under PP-Orb/ exactly as a user
    # who downloaded them would leave them.  Dojo's bare "Li.upf" sorts AHEAD of
    # SG15's "Li_ONCV_PBE-1.0.upf", which is why its presence in the list used
    # to win the resolution.
    _write(pp / "Dojo-NC-FR" / "Pseudopotential" / "Li.upf")
    _write(pp / "Dojo-NC-FR" / "Orbitals_v2.0" / "Li_DZP" / "Li_gga_7au_100Ry_4s2p2d1f.orb")
    _write(pp / "ABACUS-APNS-PPORBs-v1" / "apns-pseudopotentials-v1" / "Si.upf")
    return pp


def _write(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("dummy")
    return path


class TestBundledSeriesDetection:
    """One top-level folder under PP-Orb/ is one *series*, not one library.

    The rule this replaces ("every top-level folder holding the extension is a
    bundled library root") read one folder per library, which stopped being true
    once a series became a single folder holding both — and it admitted the
    external series folders into the default SG15 lists, so an SG15 job silently
    got Dojo pseudopotentials.  Nothing warned; the file names just were not
    SG15's.
    """

    def test_nested_series_yields_explicit_subroots(self, tmp_path, monkeypatch):
        pp = _build_pp_orb(tmp_path)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)

        pseudo = cfg._detect_library_dirs(".upf")
        orbital = cfg._detect_library_dirs(".orb")

        assert pseudo == [
            str((pp / "SG15-Version1p0" / "SG15-Version1p0_Pseudopotential").resolve()),
            str((pp / "lanthanides-f--core.icmod1" / "PD04.3+f--core.icmod1").resolve()),
        ]
        assert orbital == [
            str((pp / "SG15-Version1p0" / "SG15-Version1p0__StandardOrbitals-Version2p0").resolve()),
            str((pp / "lanthanides-f--core.icmod1" / "PD04.3+f--core.icmod1").resolve()),
        ]
        # The series folder holds both extensions, so taking it as *the* root
        # would make the two lists identical.  They must not be.
        assert pseudo != orbital

    def test_external_series_present_but_not_registered_is_excluded(self, tmp_path, monkeypatch):
        """The regression test for the reported bug.

        Dojo-NC-FR sits under PP-Orb/ and is NOT registered under
        libraries.families, so nothing in the config marks it as foreign.  Its
        files must still never be picked for an SG15 job.
        """
        from abacuscopilot.preprocessing.system_tasks import _find_file_for_element

        _build_pp_orb(tmp_path)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        assert cfg._registered_family_dirs({"libraries": {"families": {}}}) == []

        pseudo = cfg._detect_library_dirs(".upf")
        orbital = cfg._detect_library_dirs(".orb")

        assert not any("Dojo-NC-FR" in d or "ABACUS-APNS-PPORBs-v1" in d
                       for d in pseudo + orbital)
        assert _find_file_for_element(pseudo, "Li", ".upf") == "Li_ONCV_PBE-1.0.upf"
        assert _find_file_for_element(orbital, "Li", ".orb") == "Li_gga_7au_100Ry_4s1p.orb"
        # The lanthanide supplement is still reached for a 4f element.
        assert _find_file_for_element(pseudo, "Sm", ".upf") == "Sm3+_f--core-icmod1.PD04.PBE.UPF"

    def test_a_registered_family_is_excluded_by_path_too(self, tmp_path, monkeypatch):
        """A family the user registered is kept out, even under an unknown name.

        The name-based skip only knows the series shipped in the layout tables;
        a user's own library folder is recognised by *where it is*, which is what
        ``_registered_family_dirs`` feeds in.
        """
        pp = _build_pp_orb(tmp_path)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        mine = _write(pp / "my-series" / "Si.upf").parent.resolve()

        assert mine in [Path(d) for d in cfg._detect_library_dirs(".upf")]
        assert mine not in [Path(d) for d in
                            cfg._detect_library_dirs(".upf", exclude_containing=[mine])]

    def test_unknown_dropin_folder_is_bundled(self, tmp_path, monkeypatch):
        """The drop-in promise: a folder the user adds is detected, no config."""
        pp = _build_pp_orb(tmp_path)
        _write(pp / "my-own-libs" / "Na.upf")
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)

        assert str((pp / "my-own-libs").resolve()) in cfg._detect_library_dirs(".upf")

    def test_a_family_nested_under_a_series_container_keeps_the_bundled_roots(
            self, tmp_path, monkeypatch):
        """The over-broad skip the nesting used to trigger.

        The series container is no longer a candidate root, so a registered
        family sitting *beside* the libraries inside it can no longer take the
        whole series down with it.
        """
        pp = _build_pp_orb(tmp_path)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        nested = pp / "SG15-Version1p0" / "my-family"
        _write(nested / "Si.upf")

        assert str((pp / "SG15-Version1p0" / "SG15-Version1p0_Pseudopotential").resolve()) \
            in cfg._detect_library_dirs(".upf", exclude_containing=[nested])

    def test_bundled_layout_is_declared(self, tmp_path):
        from abacuscopilot.library_families import KNOWN_FAMILIES

        assert set(cfg._BUNDLED_SERIES_LAYOUT) == {"sg15", "lanthanides"}
        for top, subs in cfg._BUNDLED_SERIES_LAYOUT.values():
            assert "/" not in top
            for sub in subs.values():
                assert sub and not sub.startswith("/")
        # The external table and the family registry are two halves of one fact.
        assert set(cfg._EXTERNAL_SERIES_LAYOUT) <= set(KNOWN_FAMILIES)


class TestLibraryMigration:
    """A moved library folder must not silently drop a whole series."""

    def test_a_stale_entry_re_derives_and_keeps_the_rest(self, tmp_path, monkeypatch):
        """The rule this replaces kept the survivors and lost the rest.

        A config written before the PP-Orb reorg lists the flat SG15 path next
        to a still-valid lanthanide path.  Keeping "whatever survived" left only
        the lanthanides and every SG15 element then failed to resolve.
        """
        pp = _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        lan = (pp / "lanthanides-f--core.icmod1" / "PD04.3+f--core.icmod1").resolve()
        (tmp_path / "config.yaml").write_text(yaml.safe_dump({
            "libraries": {
                "pseudo_library": [str(pp / "SG15-Version1p0_Pseudopotential"), str(lan)],
                "orbital_library": [str(pp / "SG15-Version1p0__StandardOrbitals-Version2p0"),
                                    str(lan)],
                "family": "sg15",
            },
        }))

        libs = cfg.load_config()["libraries"]
        pseudo = [Path(d) for d in libs["pseudo_library"]]
        orbital = [Path(d) for d in libs["orbital_library"]]

        assert (pp / "SG15-Version1p0" / "SG15-Version1p0_Pseudopotential").resolve() in pseudo
        assert (pp / "SG15-Version1p0" / "SG15-Version1p0__StandardOrbitals-Version2p0").resolve() \
            in orbital
        assert lan in pseudo and lan in orbital       # the survivor is kept, not replaced
        assert not any("Dojo" in str(d) or "APNS" in str(d) for d in pseudo + orbital)

    def test_a_healthy_list_is_left_exactly_as_configured(self, tmp_path, monkeypatch):
        """A `family: custom` user's own dirs are never unioned with the bundle."""
        _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        mine = _write(tmp_path / "elsewhere" / "my-libs" / "Na.upf").parent
        (tmp_path / "config.yaml").write_text(yaml.safe_dump({
            "libraries": {"pseudo_library": [str(mine)], "family": "custom"},
        }))

        libs = cfg.load_config()["libraries"]

        assert libs["pseudo_library"] == [str(mine.resolve())]

    def test_a_stale_default_dir_falls_back_to_the_job_dir(self, tmp_path, monkeypatch):
        """An absolute path from an older wizard is a directory ABACUS cannot use."""
        _use_tmp_home(tmp_path, monkeypatch)
        (tmp_path / "config.yaml").write_text(yaml.safe_dump({
            "defaults": {"pseudo_dir": str(tmp_path / "gone"),
                         "orbital_dir": str(tmp_path)},   # this one still exists
        }))

        defaults = cfg.load_config()["defaults"]

        assert defaults["pseudo_dir"] == "./"
        assert defaults["orbital_dir"] == str(tmp_path)


class TestFamilyOwnsLibraryLists:
    """The active family decides which series the lists may resolve from.

    Regression for a server whose config said `family: sg15` while
    ``pseudo_library`` still listed Dojo-NC-FR first.  Every entry existed, so
    the stale-entry sanitizer left the list alone, and a Si CIF resolved to
    Dojo's `Si.upf` + SZ `..._1s1p.orb` instead of SG15's
    `Si_ONCV_PBE-1.0.upf` + `2s2p1d`.
    """

    def _write_config(self, tmp_path, libraries):
        (tmp_path / "config.yaml").write_text(yaml.safe_dump({"libraries": libraries}))

    def test_sg15_drops_valid_foreign_series_entries(self, tmp_path, monkeypatch):
        from abacuscopilot.preprocessing.system_tasks import _find_file_for_element

        pp = _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        dojo_p = pp / "Dojo-NC-FR" / "Pseudopotential"
        dojo_o = pp / "Dojo-NC-FR" / "Orbitals_v2.0"
        self._write_config(tmp_path, {
            "family": "sg15",
            "pseudo_library": [str(dojo_p)],
            "orbital_library": [str(dojo_o)],
        })

        libs = cfg.load_config()["libraries"]

        assert not any("Dojo" in d or "APNS" in d
                       for d in libs["pseudo_library"] + libs["orbital_library"])
        assert _find_file_for_element(libs["pseudo_library"], "Li", ".upf") == \
            "Li_ONCV_PBE-1.0.upf"
        assert _find_file_for_element(libs["orbital_library"], "Li", ".orb") == \
            "Li_gga_7au_100Ry_4s1p.orb"

    def test_a_missing_family_key_defaults_to_sg15(self, tmp_path, monkeypatch):
        pp = _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        self._write_config(tmp_path, {
            "pseudo_library": [str(pp / "Dojo-NC-FR" / "Pseudopotential")],
        })

        libs = cfg.load_config()["libraries"]

        assert not any("Dojo" in d for d in libs["pseudo_library"])

    def test_custom_family_keeps_its_hand_configured_list(self, tmp_path, monkeypatch):
        pp = _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        dojo_p = (pp / "Dojo-NC-FR" / "Pseudopotential").resolve()
        self._write_config(tmp_path, {
            "family": "custom",
            "pseudo_library": [str(dojo_p)],
        })

        libs = cfg.load_config()["libraries"]

        assert libs["pseudo_library"] == [str(dojo_p)]

    def test_an_external_family_keeps_its_own_dirs(self, tmp_path, monkeypatch):
        pp = _build_pp_orb(tmp_path)
        _use_tmp_home(tmp_path, monkeypatch)
        monkeypatch.setattr(cfg, "release_root", lambda: tmp_path)
        dojo_p = (pp / "Dojo-NC-FR" / "Pseudopotential").resolve()
        dojo_o = (pp / "Dojo-NC-FR" / "Orbitals_v2.0").resolve()
        self._write_config(tmp_path, {
            "family": "dojoncfr/dzp",
            "pseudo_library": [str(dojo_p)],
            "orbital_library": [str(dojo_o)],
        })

        libs = cfg.load_config()["libraries"]

        assert libs["pseudo_library"] == [str(dojo_p)]
        assert libs["orbital_library"] == [str(dojo_o)]
