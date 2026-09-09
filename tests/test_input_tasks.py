"""Tests for INPUT-generation tasks (task 107 ecutwfc sweep LCAO support)."""

from pathlib import Path

from abacuscopilot.preprocessing import input_tasks as it


class TestEcutwfcTestLCAO:
    def test_lcao_supported_with_reminder(self, tmp_path, monkeypatch, capsys):
        """Task 107 must NOT block LCAO: show a reminder and sweep from the
        orbital energy cutoff (e.g. 300 Ry for the APNS lanthanide NAOs)."""
        monkeypatch.chdir(tmp_path)

        # Force the LCAO branch; skip the solver prompt.
        monkeypatch.setattr(it, "_prompt_choice",
                            lambda c, q, o, d: "lcao")
        monkeypatch.setattr(it, "_ask_lcao_solver", lambda c, p: None)
        # Simulate the resolved orbital cutoff of the species in STRU.
        monkeypatch.setattr(it, "_max_orbital_ecut_from_stru", lambda: 300.0)

        prompts: dict = {}
        monkeypatch.setattr(
            it, "_prompt",
            lambda c, q, default=None: (prompts.update({q: default}) or default),
        )
        # No-op the sub-directory creation; return a usable path.
        monkeypatch.setattr(it, "_setup_convergence_subdir",
                            lambda *a, **k: Path(tmp_path))
        monkeypatch.setattr(it, "_get_sub_script_path", lambda: None)
        # Task 107 offers the pseudopotential/orbital library family
        # picker; keep the current family (blank = no change).
        monkeypatch.setattr("abacuscopilot.library_families._ask_index",
                            lambda *a, **k: "")

        it.task_ecutwfc_test(interactive=True)

        out = capsys.readouterr().out
        # Reminder about the LCAO grid-cutoff meaning is printed (no hard block).
        assert "real-space FFT grid" in out
        assert "truncated" in out
        # Sweep defaults anchored on the orbital cutoff, not 40 Ry.
        assert prompts.get("Start ecutwfc (Ry)") == "300"
        assert prompts.get("End ecutwfc (Ry)") == "400"
        assert prompts.get("Step (Ry)") == "20"

    def test_lcao_defaults_without_stru_cutoff(self, tmp_path, monkeypatch):
        """Fall back to 100→200 when the orbital cutoff can't be determined."""
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(it, "_prompt_choice",
                            lambda c, q, o, d: "lcao")
        monkeypatch.setattr(it, "_ask_lcao_solver", lambda c, p: None)
        monkeypatch.setattr(it, "_max_orbital_ecut_from_stru", lambda: 0.0)
        prompts: dict = {}
        monkeypatch.setattr(
            it, "_prompt",
            lambda c, q, default=None: (prompts.update({q: default}) or default),
        )
        monkeypatch.setattr(it, "_setup_convergence_subdir",
                            lambda *a, **k: Path(tmp_path))
        monkeypatch.setattr(it, "_get_sub_script_path", lambda: None)
        # Task 107 offers the pseudopotential/orbital library family
        # picker; keep the current family (blank = no change).
        monkeypatch.setattr("abacuscopilot.library_families._ask_index",
                            lambda *a, **k: "")

        it.task_ecutwfc_test(interactive=True)
        assert prompts.get("Start ecutwfc (Ry)") == "100"
        assert prompts.get("End ecutwfc (Ry)") == "200"


class _FakeConsole:
    """Minimal console: print() captured so tests can assert warnings."""

    def __init__(self):
        self.out: list[str] = []

    def print(self, *args, **kwargs):
        self.out.append(" ".join(str(a) for a in args))


class _Answers:
    """Feeds answers to _prompt_choice / _prompt in call order.

    Calls with an exhausted queue fall back to their *default* argument, so an
    empty queue exercises the default path of every prompt.
    """

    def __init__(self, choices=(), prompts=()):
        self.choices = list(choices)
        self.prompts = list(prompts)

    def choice(self, console, question, options, default):
        return self.choices.pop(0) if self.choices else default

    def prompt(self, console, question, default=None):
        return self.prompts.pop(0) if self.prompts else default


class TestSharedCalcHelpers:
    """Shared per-calculation helpers (XC/D3, analysis outputs, band/DOS/MD)
    used by BOTH the INPUT module tasks and the STRU full-calc branch."""

    def test_ask_xc_functional_default_is_pbe(self, monkeypatch):
        console, a = _FakeConsole(), _Answers()          # no answers → defaults
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_xc_functional(console, p)
        assert p.dft_functional == "pbe"

    def test_ask_xc_functional_pbesol_selected(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["PBEsol"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_xc_functional(console, p)
        assert p.dft_functional == "pbesol"

    def test_ask_d3_none_keeps_default(self, monkeypatch):
        console, a = _FakeConsole(), _Answers()          # default "No"
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_d3(console, p)
        assert p.vdw_method == "none"

    def test_ask_d3_active_forces_pbe_from_pbesol(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["d3_0 (zero-damping)"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        p.dft_functional = "pbesol"
        it._ask_d3(console, p)
        assert p.vdw_method == "d3_0"
        assert p.dft_functional == "pbe"                 # forced back
        assert "does not support" in "".join(console.out)

    def test_ensure_xc_in_core_only_registers_non_default(self):
        p = it.InputParams()                            # pbe / none = defaults
        it._ensure_xc_in_core(p)
        assert p.extras.get("_template_keys", []) == []

        p2 = it.InputParams()
        p2.dft_functional = "pbesol"                     # non-default
        p2.vdw_method = "d3_0"                           # non-default
        p2.extras["_comment_hints"] = {"vdw_method": "d3_0", "dft_functional": "pbesol"}
        it._ensure_xc_in_core(p2)
        assert sorted(p2.extras["_template_keys"]) == ["dft_functional", "vdw_method"]
        assert "dft_functional" not in p2.extras["_comment_hints"]
        assert "vdw_method" not in p2.extras["_comment_hints"]

    def test_ask_xc_and_d3_default_relax(self, monkeypatch):
        console, a = _FakeConsole(), _Answers()          # PBE + No
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_xc_and_d3(console, p)
        assert p.dft_functional == "pbe"
        assert p.vdw_method == "none"
        assert p.extras.get("_template_keys", []) == []

    def test_ask_xc_and_d3_pbesol_registers_hint_pop(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["PBEsol", "No"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        p.extras["_comment_hints"] = {"vdw_method": "d3_0", "dft_functional": "pbesol"}
        it._ask_xc_and_d3(console, p)
        assert p.dft_functional == "pbesol"
        assert p.extras["_template_keys"] == ["dft_functional"]
        assert "dft_functional" not in p.extras["_comment_hints"]

    def test_ask_analysis_outputs_none(self, monkeypatch):
        console, a = _FakeConsole(), _Answers()          # default "none"
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(it, "_abacus_supports_hirshfeld", lambda: True)
        p = it.InputParams()
        it._ask_analysis_outputs(console, p, "lcao")
        assert p.out_chg is False
        assert p.out_mul is False

    def test_ask_analysis_outputs_all_lcao(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["all (Mulliken + Hirshfeld + CHG)"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(it, "_abacus_supports_hirshfeld", lambda: True)
        p = it.InputParams()
        it._ask_analysis_outputs(console, p, "lcao")
        assert p.out_chg is True
        assert p.out_mul is True
        assert p.get_param("out_hirshfeld") == 1

    def test_ask_band_projection(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["Yes (band + projection)"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_band_projection(console, p)
        assert p.out_proj_band is True

    def test_ask_dos_type(self, monkeypatch):
        console, a = _FakeConsole(), _Answers(choices=["Projected DOS (PDOS)"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        p = it.InputParams()
        it._ask_dos_type(console, p)
        assert p.out_dos == 2

    def test_ask_md_thermo_default_nvt(self, monkeypatch):
        console, a = _FakeConsole(), _Answers()          # nvt default
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(it, "_prompt", a.prompt)
        p = it.InputParams()
        it._ask_md_thermo(console, p, "lcao")
        assert p.md_type == "nvt"
        assert p.md_nstep == 10000
        assert p.md_dt == 1.0
        assert p.md_tfirst == 300.0
        assert p.md_tlast == 300.0

    def test_ask_md_thermo_npt_iso(self, monkeypatch):
        console = _FakeConsole()
        a = _Answers(choices=["npt", "iso (isotropic)"],
                     prompts=["2.5", "10000", "1.0", "400.0", "400.0"])
        monkeypatch.setattr(it, "_prompt_choice", a.choice)
        monkeypatch.setattr(it, "_prompt", a.prompt)
        p = it.InputParams()
        it._ask_md_thermo(console, p, "lcao")
        assert p.md_type == "npt"
        assert p.md_pmode == "iso"
        assert (p.press1, p.press2, p.press3) == (2.5, 2.5, 2.5)
        assert p.md_tfirst == 400.0
