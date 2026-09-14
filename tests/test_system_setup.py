"""Wizard 9901 and display 9902 — what they ask, and what they write.

``load_config`` is patched in every test here, but that alone is not enough:
the wizard *saves*, and on a developer machine ``~/.abacuscopilot/config.yaml``
is a live file.  Patching only the read would have the tests overwrite it, so
``save_config`` is replaced alongside and the written dict is inspected instead.
``load_config``'s replacement returns a deep copy for the same reason — the
module-level ``DEFAULT_CONFIG`` must not be edited by a test.
"""

import copy

from abacuscopilot import config as cfg
from abacuscopilot.preprocessing import system_tasks as st


class _Wizard:
    """Runs a task with a recording ``_prompt`` and a captured ``save_config``."""

    def __init__(self, monkeypatch, start=None):
        self.questions: dict[str, object] = {}
        self.saved: dict = {}
        self.config = copy.deepcopy(start if start is not None else cfg.DEFAULT_CONFIG)
        monkeypatch.setattr(st, "load_config", lambda: self.config)
        monkeypatch.setattr(st, "save_config", self._save)
        monkeypatch.setattr(st, "_prompt", self._prompt)

    def _save(self, config):
        self.saved = copy.deepcopy(config)

    def _prompt(self, console, question, default=None):
        """Answer with the default — the wizard's "just press Enter" path."""
        self.questions[question] = default
        return str(default) if default is not None else ""


def test_9901_does_not_ask_for_pseudo_or_orbital_dirs(tmp_path, monkeypatch, capsys):
    """The wizard asks for binaries and numbers, never for library paths.

    Asking was the old behaviour, and the answer it collected was written into
    ``defaults.pseudo_dir`` as an absolute path that then went stale on every
    upgrade.  The library is chosen once, in the INPUT flow, as a *series*.
    """
    wizard = _Wizard(monkeypatch)

    st.task_system_setup(interactive=True)

    asked = " ".join(wizard.questions).lower()
    assert "pseudopotential" not in asked
    assert "orbital" not in asked
    assert "directory" not in asked
    assert "kspacing" in asked or "k-spacing" in asked   # the rest still runs

    # Both stay at the one value ABACUS understands as "the job directory".
    assert wizard.saved["defaults"]["pseudo_dir"] == "./"
    assert wizard.saved["defaults"]["orbital_dir"] == "./"
    # …and the step is still announced, so the wizard does not look truncated.
    assert "auto-detected from PP-Orb/" in capsys.readouterr().out


def test_9901_writes_force_thr_ev_never_force_thr(monkeypatch):
    """The old key was named ``force_thr``, which nothing in the codebase reads.

    Whatever the user typed was stored where no reader would find it, so the
    threshold silently kept its built-in value.
    """
    wizard = _Wizard(monkeypatch)

    st.task_system_setup(interactive=True)

    assert "force_thr_ev" in wizard.saved["defaults"]
    assert "force_thr" not in wizard.saved["defaults"]


def test_9901_prompt_defaults_come_from_the_schema(monkeypatch):
    """Four literals in the wizard disagreed with DEFAULT_CONFIG.

    kspacing was offered as 0.04 against a schema 0.14, basis_type as "pw"
    against "lcao", and the force threshold as 0.001 — so pressing Enter on the
    defaults produced a config that meant something else than the file it was
    loaded from.  This fails on any of those literals.
    """
    wizard = _Wizard(monkeypatch)

    st.task_system_setup(interactive=True)

    assert wizard.questions["Default k-spacing (1/bohr, ABACUS unit)"] == "0.14"
    assert wizard.questions["Default basis type (lcao / pw)"] == "lcao"
    assert wizard.questions["Default force convergence threshold (eV/Å)"] == "0.01"
    # Accepting every default must round-trip the schema unchanged.
    for key in ("kspacing", "ecutwfc", "scf_thr", "force_thr_ev",
                "basis_type", "dft_functional"):
        assert wizard.saved["defaults"][key] == cfg.DEFAULT_CONFIG["defaults"][key]


def test_9902_shows_force_thr_ev(monkeypatch, capsys):
    """A key that exists must not be displayed as the missing-value placeholder."""
    monkeypatch.setattr(st, "load_config",
                        lambda: {**copy.deepcopy(cfg.DEFAULT_CONFIG),
                                 "defaults": {**copy.deepcopy(cfg.DEFAULT_CONFIG["defaults"]),
                                              "force_thr_ev": 0.02}})

    st.task_show_config()

    out = capsys.readouterr().out
    assert "force_thr_ev: 0.02" in out
    assert "force_thr:" not in out
