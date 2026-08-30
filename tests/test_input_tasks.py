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

        it.task_ecutwfc_test(interactive=True)
        assert prompts.get("Start ecutwfc (Ry)") == "100"
        assert prompts.get("End ecutwfc (Ry)") == "200"
