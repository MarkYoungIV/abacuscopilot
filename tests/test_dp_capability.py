"""Tests for DP model / ABACUS-DP capability detection.

The regression these guard: a DP-MD INPUT generator used to run
``deepmd.infer.DeepPot(model)`` as a "version mismatch" probe and, on *any*
exception, rename the model to ``*_original`` and convert it in place.  DeepMD
v3 dispatches ``.pb`` to its TensorFlow backend, so that probe raises for every
``.pb`` when the env lacks tensorflow — and working models were silently
replaced by graphs ABACUS could not read (empty type map at run time).
"""

from __future__ import annotations

from pathlib import Path

from abacuscopilot import dp_capability as dp

CHUNK = 1 << 20  # matches dp._scan_markers' default


def _fake_graph(*, tmap: bytes | None = b"Li Ge P S",
                version: bytes | None = b"1.1") -> bytes:
    """Bytes shaped like a DeepMD-kit v2 TensorFlow frozen graph."""
    parts = [
        b"\n;\n\x05t_box\x12\x0bPlaceholder",
        b"descrpt_attr/ntypes\x12\x01i",
        b"fitting_attr/dfparam\x12\x010",
        b"model_attr/model_type\x12\x04ener",
    ]
    if tmap is not None:
        parts.append(b"model_attr/tmap" + tmap)
    if version is not None:
        parts.append(b"model_attr/model_version\x12\x05Constvalue" + version)
    return b"".join(parts)


def _cap(**kw) -> dp.DPCapability:
    base = dict(binary="abacus", resolved="/opt/bin/abacus", has_dp=True,
                deepmd_version="2.2.11", source="binary_probe")
    base.update(kw)
    return dp.DPCapability(**base)


# ---------------------------------------------------------------------------
# Model fingerprint
# ---------------------------------------------------------------------------


class TestFingerprintModel:
    def test_missing_file(self, tmp_path):
        fp = dp.fingerprint_model(tmp_path / "nope.pb")
        assert fp.exists is False
        assert fp.kind == "missing"
        assert fp.is_dp2_tf is False

    def test_directory_is_not_a_model(self, tmp_path):
        assert dp.fingerprint_model(tmp_path).kind == "missing"

    def test_empty_file(self, tmp_path):
        p = tmp_path / "empty.pb"
        p.write_bytes(b"")
        fp = dp.fingerprint_model(p)
        assert (fp.exists, fp.kind, fp.size) == (True, "empty", 0)

    def test_dpmodel_zip(self, tmp_path):
        p = tmp_path / "model.pb"
        p.write_bytes(b"PK\x03\x04" + b"x" * 100)
        assert dp.fingerprint_model(p).kind == "dpmodel_zip"

    def test_healthy_dp2_graph(self, tmp_path):
        p = tmp_path / "graph.pb"
        p.write_bytes(_fake_graph())
        fp = dp.fingerprint_model(p)
        assert fp.kind == "tf_graphdef"
        assert fp.has_tmap is True
        assert fp.model_version == "1.1"
        assert fp.is_dp2_tf is True

    def test_graph_without_tmap_is_still_a_graph(self, tmp_path):
        """The converted-away type map must be reported precisely.

        Treating this as "not a TF graph" would hide the single most common way
        a model goes bad, so the two cases are deliberately distinguished.
        """
        p = tmp_path / "graph.pb"
        p.write_bytes(_fake_graph(tmap=None))
        fp = dp.fingerprint_model(p)
        assert fp.kind == "tf_graphdef"
        assert fp.has_tmap is False
        assert fp.is_dp2_tf is False
        assert "model_attr/tmap" in fp.describe()

    def test_missing_version_node_is_reported_not_hidden(self, tmp_path):
        """Fingerprinting it is fine; *reporting* it as healthy is not.

        DeePMD-kit reads an absent model_attr/model_version as 0.0 and aborts
        on anything below 1.1. Models published in the older uncompressed DP
        format look exactly like this, so the missing node has to be visible.
        """
        p = tmp_path / "graph.pb"
        p.write_bytes(_fake_graph(version=None))
        fp = dp.fingerprint_model(p)
        assert (fp.kind, fp.has_tmap, fp.model_version) == ("tf_graphdef", True, None)
        assert fp.has_version is False
        assert fp.usable is False
        assert "model_version" in fp.describe()

    def test_garbage_is_unknown(self, tmp_path):
        p = tmp_path / "junk.pb"
        p.write_bytes(b"\x00\x01\x02not a graph at all")
        assert dp.fingerprint_model(p).kind == "unknown"


class TestScanMarkers:
    def test_marker_straddling_a_chunk_boundary(self, tmp_path):
        """Markers must survive the chunked read, and offsets must be exact."""
        pad = CHUNK - 3 - len(b"descrpt_attr/")
        p = tmp_path / "big.pb"
        p.write_bytes(b"descrpt_attr/" + b"x" * pad + dp._TMAP_MARKER + b"Li Ge P S")
        assert p.stat().st_size > CHUNK

        found = dp._scan_markers(p, (dp._TMAP_MARKER,))
        assert found[dp._TMAP_MARKER] == CHUNK - 3
        # ...and the version window read stays consistent with that offset.
        assert dp.fingerprint_model(p).has_tmap is True

    def test_all_markers_found_in_one_pass(self, tmp_path):
        p = tmp_path / "g.pb"
        p.write_bytes(_fake_graph())
        found = dp._scan_markers(p, (*dp._TF_GRAPH_MARKERS, dp._TMAP_MARKER))
        assert dp._TMAP_MARKER in found
        assert all(m in found for m in dp._TF_GRAPH_MARKERS)


# ---------------------------------------------------------------------------
# Binary capability
# ---------------------------------------------------------------------------


class TestClassifyDpSupport:
    def test_symbols_mean_dp_is_linked(self):
        assert dp._classify_dp_support("... deepmd::DeepPot ... ProdEnvMatA") is True

    def test_recompile_message_means_no_dp(self):
        blob = "Please recompile with DeePMD-kit (-DUSE_DEEPMD=ON)"
        assert dp._classify_dp_support(blob) is False

    def test_word_deepmd_alone_proves_nothing(self):
        assert dp._classify_dp_support("deepmd mentioned in a doc string") is None

    def test_blank_blob_is_unknown(self):
        assert dp._classify_dp_support("") is None


class TestParseVersionAndLinkage:
    def test_version_only_from_deepmd_lines(self):
        lines = ["/usr/lib/libfoo.so.1.2.3", "/opt/deepmd-2.2.11/lib/libdeepmd.so"]
        assert dp._parse_version(lines) == "2.2.11"

    def test_no_version_when_absent(self):
        assert dp._parse_version(["/usr/lib/libdeepmd.so"]) is None

    def test_build_path_yields_no_version(self):
        """A conda env name contains "deepmd" but names no version.

        This is the real line ABACUS's binary carries: the first X.Y.Z on it is
        ABACUS's own 3.10.0, from abacus-develop-LTSv3.10.0.  Reporting that as
        the DeepMD version flips the verdict to a spurious warning.
        """
        line = (
            "/home/young/miniconda3/envs/deepmd_v2/lib:"
            "/home/young/softwares/abacus-develop-LTSv3.10.0/toolchain/install/"
            "openmpi-5.0.6/lib"
        )
        assert dp._parse_version([line]) is None

    def test_unknown_prefixed_line_does_not_mask_a_later_real_one(self):
        lines = ["/opt/envs/deepmd_v2/lib", "/opt/deepmd-kit-2.2.11/lib/libdeepmd.so"]
        assert dp._parse_version(lines) == "2.2.11"

    def test_version_from_linked_lib_reads_the_owning_env(self, tmp_path):
        libdir = tmp_path / "envs" / "deepmd_v2" / "lib"
        libdir.mkdir(parents=True)
        lib = libdir / "libdeepmd_c.so"
        lib.touch()
        info = libdir / "python3.11" / "site-packages" / "deepmd_kit-2.2.11.dist-info"
        info.mkdir(parents=True)
        (info / "METADATA").write_text("Metadata-Version: 2.1\nVersion: 2.2.11\n")

        assert dp._version_from_linked_lib([str(lib)]) == "2.2.11"

    def test_version_from_linked_lib_ignores_non_deepmd_libs(self, tmp_path):
        assert dp._version_from_linked_lib([str(tmp_path / "libfoo.so")]) is None

    def test_linkage_picks_out_real_paths(self):
        ldd = [
            "libdeepmd.so => /opt/deepmd-2.2.11/lib/libdeepmd.so (0x00007f)",
            "libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007e)",
        ]
        linkage, libs = dp._linkage_info(ldd)
        assert linkage == "dynamic"
        assert libs == ["/opt/deepmd-2.2.11/lib/libdeepmd.so"]

    def test_no_deepmd_libs_is_unknown(self):
        assert dp._linkage_info(["libc.so.6 => /lib/libc.so.6"]) == ("unknown", [])


class TestProbeAbacusDp:
    def _fake_binary(self, tmp_path) -> Path:
        p = tmp_path / "abacus"
        p.write_bytes(b"\x7fELF" + b"\x00" * 64)
        return p

    def test_missing_binary_is_not_an_error(self, tmp_path):
        cap = dp.probe_abacus_dp(str(tmp_path / "nope"), config={"paths": {}},
                                 use_cache=False)
        assert cap.resolved is None
        assert cap.has_dp is None
        assert cap.source == "not_found"

    def test_probes_symbols_and_linkage(self, tmp_path, monkeypatch):
        binary = self._fake_binary(tmp_path)
        monkeypatch.setattr(dp, "save_machine_capability", lambda cap: None)
        blobs = {
            ("strings", str(binary)): "deepmd::DeepPot ProdEnvMatA\n",
            ("ldd", str(binary)): (
                "libdeepmd.so => /opt/deepmd-2.2.11/lib/libdeepmd.so (0x1)\n"),
        }
        monkeypatch.setattr(dp, "_run", lambda cmd, timeout=30: blobs.get(tuple(cmd), ""))

        cap = dp.probe_abacus_dp(str(binary), config={"paths": {}}, use_cache=False)
        assert cap.has_dp is True
        assert cap.source == "binary_probe"
        assert cap.linkage == "dynamic"
        assert cap.deepmd_version == "2.2.11"

    def test_static_binary_has_no_ldd_output(self, tmp_path, monkeypatch):
        binary = self._fake_binary(tmp_path)
        monkeypatch.setattr(dp, "save_machine_capability", lambda cap: None)
        monkeypatch.setattr(dp, "_run", lambda cmd, timeout=30: (
            "deepmd::DeepPot\n" if cmd[0] == "strings" else ""))

        cap = dp.probe_abacus_dp(str(binary), config={"paths": {}}, use_cache=False)
        assert (cap.has_dp, cap.linkage) == (True, "static")

    def test_explicit_config_overrides_detection(self, tmp_path, monkeypatch):
        binary = self._fake_binary(tmp_path)
        monkeypatch.setattr(dp, "save_machine_capability", lambda cap: None)
        monkeypatch.setattr(dp, "_run", lambda cmd, timeout=30: "deepmd::DeepPot\n")

        cap = dp.probe_abacus_dp(str(binary),
                                 config={"paths": {"abacus_dp_version": "2.2.11"}},
                                 use_cache=False)
        assert cap.deepmd_version == "2.2.11"
        assert cap.source.endswith("+explicit")

    def test_dp_binary_preferred_over_generic(self, tmp_path, monkeypatch):
        dp_bin = self._fake_binary(tmp_path)
        monkeypatch.setattr(dp, "save_machine_capability", lambda cap: None)
        monkeypatch.setattr(dp, "_run", lambda cmd, timeout=30: "")
        config = {"paths": {"abacus_dp_binary": str(dp_bin), "abacus_binary": "abacus"}}
        cap = dp.probe_abacus_dp(None, config=config, use_cache=False)
        assert cap.resolved == str(dp_bin)


class TestRuntimeEnvHint:
    def test_dynamic_without_env_file_warns(self):
        assert dp.check_runtime_env_hint(_cap(linkage="dynamic", runtime_env=""))

    def test_dynamic_with_env_file_is_quiet(self):
        assert dp.check_runtime_env_hint(
            _cap(linkage="dynamic", runtime_env="/e.sh")) is None

    def test_static_is_quiet(self):
        assert dp.check_runtime_env_hint(_cap(linkage="static")) is None


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


class TestAssessCompatibility:
    def _fp(self, tmp_path, blob: bytes, name="graph.pb") -> dp.ModelFingerprint:
        p = tmp_path / name
        p.write_bytes(blob)
        return dp.fingerprint_model(p)

    def test_healthy_model_with_matching_binary(self, tmp_path):
        verdict = dp.assess_compatibility(self._fp(tmp_path, _fake_graph()), _cap())
        assert verdict.status == dp.OK
        assert verdict.is_problem is False

    def test_graph_without_version_node_is_a_hard_failure(self, tmp_path):
        """The 192.168.8.27 abort, caught before submission instead of by ABACUS.

        A freshly downloaded, never-compressed model has its type map and
        passes every other check — then DeePMD-kit's C API reads version 0.0
        and throws "incompatable model: version 0.0 in graph, but version 1.1
        supported" (signal 6, core dumped) partway into the run.
        """
        fp = self._fp(tmp_path, _fake_graph(version=None))
        verdict = dp.assess_compatibility(fp, _cap())
        assert verdict.status == dp.INCOMPATIBLE
        assert "model_version" in verdict.message
        assert "dp compress" in verdict.guidance

    def test_converted_away_tmap_is_a_hard_failure(self, tmp_path):
        fp = self._fp(tmp_path, _fake_graph(tmap=None))
        verdict = dp.assess_compatibility(fp, _cap())
        assert verdict.status == dp.INCOMPATIBLE
        assert "type map" in verdict.message
        assert "_original" in verdict.guidance

    def test_empty_file_is_a_hard_failure(self, tmp_path):
        fp = self._fp(tmp_path, b"")
        assert dp.assess_compatibility(fp, _cap()).status == dp.INCOMPATIBLE

    def test_binary_without_dp_support(self, tmp_path):
        fp = self._fp(tmp_path, _fake_graph())
        verdict = dp.assess_compatibility(fp, _cap(has_dp=False))
        assert verdict.status == dp.INCOMPATIBLE
        assert "DeepMD" in verdict.message

    def test_dpmodel_zip_needs_conversion(self, tmp_path):
        fp = self._fp(tmp_path, b"PK\x03\x04zipzip")
        assert dp.assess_compatibility(fp, _cap()).status == dp.NEEDS_CONVERSION

    def test_model_missing(self, tmp_path):
        fp = dp.fingerprint_model(tmp_path / "gone.pb")
        assert dp.assess_compatibility(fp, _cap()).status == dp.MISSING

    def test_unknown_binary_capability_warns_not_fails(self, tmp_path):
        fp = self._fp(tmp_path, _fake_graph())
        verdict = dp.assess_compatibility(fp, _cap(has_dp=None, deepmd_version=None))
        assert verdict.status == dp.WARN
        assert verdict.is_problem is True  # surfaced to the user, but not fatal

    def test_dp3_binary_warns_but_does_not_block(self, tmp_path):
        fp = self._fp(tmp_path, _fake_graph())
        verdict = dp.assess_compatibility(fp, _cap(deepmd_version="3.1.2"))
        assert verdict.status == dp.WARN
        assert "试算" in verdict.guidance


# ---------------------------------------------------------------------------
# Per-host cache
# ---------------------------------------------------------------------------


class TestMachineCache:
    def test_round_trip(self, tmp_path, monkeypatch):
        cache = tmp_path / "machines.yaml"
        monkeypatch.setattr(dp, "_machines_path", lambda: cache)

        dp.save_machine_capability(_cap(linked_libs=["/opt/libdeepmd.so"],
                                        runtime_env="/e.sh", probed_at="2026-09-11"))
        loaded = dp._load_cached_capability("/opt/bin/abacus")
        assert loaded is not None
        assert loaded.deepmd_version == "2.2.11"
        assert loaded.linked_libs == ["/opt/libdeepmd.so"]
        assert loaded.probed_at == "2026-09-11"

    def test_different_binary_invalidates_cache(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dp, "_machines_path", lambda: tmp_path / "machines.yaml")
        dp.save_machine_capability(_cap())
        assert dp._load_cached_capability("/somewhere/else/abacus") is None

    def test_corrupt_cache_is_ignored(self, tmp_path, monkeypatch):
        cache = tmp_path / "machines.yaml"
        monkeypatch.setattr(dp, "_machines_path", lambda: cache)
        cache.write_text("::: not yaml :::\n\t- [")
        assert dp._load_cached_capability("/opt/bin/abacus") is None

    def test_unknown_keys_in_cache_are_dropped(self, tmp_path, monkeypatch):
        """A cache written by a newer version must not crash this one."""
        cache = tmp_path / "machines.yaml"
        monkeypatch.setattr(dp, "_machines_path", lambda: cache)
        dp.save_machine_capability(_cap())
        cache.write_text(cache.read_text().replace(
            "resolved: /opt/bin/abacus",
            "resolved: /opt/bin/abacus\n      future_field: 42"))
        loaded = dp._load_cached_capability("/opt/bin/abacus")
        assert loaded is not None and loaded.binary == "abacus"


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscoverModels:
    def test_finds_models_and_ignores_everything_else(self, tmp_path):
        (tmp_path / "LiGePS.pb").write_bytes(_fake_graph())
        (tmp_path / "weights.pt").write_bytes(b"PK\x03\x04torch")
        (tmp_path / "case.md").write_bytes(b"# notes")
        (tmp_path / "traj.pb").mkdir()  # a directory that looks like a model

        names = [Path(fp.path).name for fp in dp.discover_models(tmp_path)]
        assert names == ["LiGePS.pb", "weights.pt"]

    def test_usable_models_rank_first(self, tmp_path):
        """A broken graph must not be the default just because of its name."""
        (tmp_path / "aaa-broken.pb").write_bytes(_fake_graph(tmap=None))
        (tmp_path / "zzz-good.pb").write_bytes(_fake_graph())

        found = dp.discover_models(tmp_path)
        assert Path(found[0].path).name == "zzz-good.pb"
        assert found[0].usable is True

    def test_broken_model_is_still_returned(self, tmp_path):
        """Dropping it would leave the user with an empty list and no clue."""
        (tmp_path / "graph-compress.pb").write_bytes(_fake_graph(tmap=None))
        found = dp.discover_models(tmp_path)
        assert [Path(fp.path).name for fp in found] == ["graph-compress.pb"]
        assert found[0].has_tmap is False

    def test_loadable_model_outranks_an_uncompressed_download(self, tmp_path):
        """Both have a type map; only one has the version node ABACUS needs."""
        (tmp_path / "LiGePS-SSE-PBEsol-model.pb").write_bytes(
            _fake_graph(version=None))
        (tmp_path / "graph-compress.pb").write_bytes(_fake_graph())

        found = dp.discover_models(tmp_path)
        assert [Path(fp.path).name for fp in found] == [
            "graph-compress.pb", "LiGePS-SSE-PBEsol-model.pb"]
        assert [fp.usable for fp in found] == [True, False]

    def test_backup_siblings_are_not_models(self, tmp_path):
        """*_original / *.broken are recovery copies, not `pot_file` candidates."""
        (tmp_path / "graph-compress.pb").write_bytes(_fake_graph())
        (tmp_path / "graph-compress.pb_original").write_bytes(_fake_graph())
        (tmp_path / "graph-compress.pb.broken").write_bytes(_fake_graph())

        names = [Path(fp.path).name for fp in dp.discover_models(tmp_path)]
        assert names == ["graph-compress.pb"]

    def test_missing_directory_is_empty(self, tmp_path):
        assert dp.discover_models(tmp_path / "nope") == []


class TestTorchModelFingerprint:
    def test_pt_is_a_torch_model_not_a_v3_archive(self, tmp_path):
        """Both are zips; only the suffix tells them apart.

        Reporting a `.pt` as a "v3 backend archive" would send the user off to
        convert a model that was already in the format their build wants.
        """
        p = tmp_path / "model.pt"
        p.write_bytes(b"PK\x03\x04" + b"x" * 100)
        fp = dp.fingerprint_model(p)
        assert fp.kind == "torch_model"
        assert fp.is_torch is True
        assert fp.is_dp2_tf is False
        assert "PyTorch" in fp.describe()

    def test_pb_that_is_really_a_zip_stays_a_v3_archive(self, tmp_path):
        p = tmp_path / "model.pb"
        p.write_bytes(b"PK\x03\x04" + b"x" * 100)
        assert dp.fingerprint_model(p).kind == "dpmodel_zip"

    def test_torch_model_warns_rather_than_blocks(self, tmp_path):
        p = tmp_path / "model.pt"
        p.write_bytes(b"PK\x03\x04" + b"x" * 100)
        verdict = dp.assess_compatibility(dp.fingerprint_model(p), _cap())
        assert verdict.status == dp.WARN
        assert "convert-backend" in verdict.guidance
