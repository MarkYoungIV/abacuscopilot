"""Tests for the FFT-based MSD computation (task 3103) and frame timing."""


import numpy as np

from abacuscopilot.postprocessing.md_tasks import (
    _compute_msd_fft,
    _resolve_frame_dt,
)


def _naive_msd(traj, max_lag):
    """Per-atom mean MSD via the direct pair loop (reference)."""
    n_frames = traj.shape[0]
    out = {}
    axes = {"total": [0, 1, 2], "x": [0], "y": [1], "z": [2]}
    for label, comps in axes.items():
        msd = np.zeros(max_lag + 1)
        for lag in range(1, max_lag + 1):
            d = traj[lag:, :, comps] - traj[:-lag, :, comps]
            msd[lag] = np.mean(np.sum(d * d, axis=2))  # sum comps, mean t+atoms
        out[label] = msd
    return out


class TestComputeMsdFft:
    def test_matches_naive(self):
        rng = np.random.default_rng(7)
        traj = np.cumsum(rng.normal(0, 0.3, (300, 8, 3)), axis=0)
        fft = _compute_msd_fft(traj, 299, {"total": [0, 1, 2], "x": [0]})
        ref = _naive_msd(traj, 299)
        assert np.max(np.abs(fft["total"] - ref["total"])) < 1e-9
        assert np.max(np.abs(fft["x"] - ref["x"])) < 1e-9

    def test_zero_lag_is_zero(self):
        rng = np.random.default_rng(1)
        traj = rng.normal(0, 1, (50, 4, 3))
        res = _compute_msd_fft(traj, 49, {"total": [0, 1, 2]})
        assert res["total"][0] == 0.0

    def test_is_per_atom_mean(self):
        """MSD should be the per-atom mean, not the sum over atoms."""
        rng = np.random.default_rng(2)
        traj = rng.normal(0, 1, (40, 6, 3))
        res = _compute_msd_fft(traj, 39, {"total": [0, 1, 2]})
        # Direct check at lag 1: per-atom mean over t of |r(t+1)-r(t)|^2.
        d = traj[1:] - traj[:-1]  # (39, 6, 3)
        expected = np.mean(np.sum(d * d, axis=2))  # sum comps, mean over t+atoms
        assert abs(res["total"][1] - expected) < 1e-9

    def test_scale(self):
        """MSD of a stationary cloud stays ~2·Var·dims; diffusive cloud grows."""
        rng = np.random.default_rng(3)
        stationary = rng.normal(0, 0.1, (200, 5, 3))          # jitter: MSD≈6σ²=0.06
        diffusive = np.cumsum(rng.normal(0, 0.2, (200, 5, 3)), axis=0)
        s = _compute_msd_fft(stationary, 199, {"total": [0, 1, 2]})["total"]
        d = _compute_msd_fft(diffusive, 199, {"total": [0, 1, 2]})["total"]
        assert s[50] < 0.2            # bounded (stationary), ~0.06
        assert d[50] > 1.0            # grows with t


class TestResolveFrameDt:
    def test_auto_detect_from_input(self, tmp_path, monkeypatch, capsys):
        """ABACUS-native: md_dt × dumpfreq read from INPUT."""
        (tmp_path / "INPUT").write_text("md_dt  2.0\nmd_dumpfreq  5\n")
        monkeypatch.chdir(tmp_path)
        from abacuscopilot.console_utils import _get_console
        fd = _resolve_frame_dt(tmp_path / "MD_dump", _get_console(), interactive=True)
        assert fd == 10.0
        assert "frame_dt = 10.0" in capsys.readouterr().out

    def test_prompts_when_no_input(self, tmp_path, monkeypatch, capsys):
        """Converted VASP trajectory (no INPUT): ask the user for the spacing."""
        monkeypatch.chdir(tmp_path)
        from abacuscopilot.postprocessing import md_tasks as mt
        monkeypatch.setattr(mt, "_prompt", lambda c, q, default=None: "2.0")
        from abacuscopilot.console_utils import _get_console
        fd = _resolve_frame_dt(tmp_path / "MD_dump", _get_console(), interactive=True)
        assert fd == 2.0

    def test_noninteractive_defaults_to_1(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        from abacuscopilot.console_utils import _get_console
        fd = _resolve_frame_dt(tmp_path / "MD_dump", _get_console(), interactive=False)
        assert fd == 1.0
        assert "assuming frame_dt = 1.0" in capsys.readouterr().out
