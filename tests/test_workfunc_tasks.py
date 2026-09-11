"""Regression tests for ABACUS work-function tasks 1101 and 1102."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np

from abacuscopilot.core.constants import RY_TO_EV
from abacuscopilot.postprocessing.workfunc_tasks import (
    macroscopic_average,
    read_fermi_energy,
    read_potential_cube,
    read_vacuum_level,
    task_macro_avg_potential,
    task_work_function,
)


def _write_cube(path, profile_ry: np.ndarray) -> None:
    """Write a small x/y/z cube with a known z profile in Rydberg."""
    nx, ny, nz = 2, 2, len(profile_ry)
    values = np.broadcast_to(profile_ry, (nx, ny, nz)).ravel()
    rows = [
        "ABACUS electrostatic potential",
        "synthetic test cube",
        "1 0.0 0.0 0.0",
        f"{nx} 1.0 0.0 0.0",
        f"{ny} 0.0 1.0 0.0",
        f"{nz} 0.0 0.0 1.0",
        "1 1.0 0.0 0.0 0.0",
    ]
    rows.extend(" ".join(f"{value:.10f}" for value in values[i:i + 6])
                for i in range(0, len(values), 6))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_cube_reader_preserves_z_axis_and_units(tmp_path):
    profile = np.linspace(1.0, 2.0, 8)
    cube = tmp_path / "ElecStaticPot.cube"
    _write_cube(cube, profile)

    result = read_potential_cube(cube)
    assert result is not None
    data, cell, origin = result
    assert data.shape == (2, 2, 8)
    assert np.allclose(data.mean(axis=(0, 1)), profile)
    assert np.allclose(cell[2], [0.0, 0.0, 8.0])
    assert np.allclose(origin, [0.0, 0.0, 0.0])


def test_vacuum_and_fermi_parsers_use_final_electronic_values(tmp_path):
    vacuum = tmp_path / "E_vacuum.out"
    vacuum.write_text("E_VACUUM (eV) = 5.75000000000 at z (Angstrom) = 10.0\n")
    log = tmp_path / "running_scf.log"
    log.write_text(
        "E_Fermi        0.1000000000        1.3605693123\n"
        "E_Fermi        0.2000000000 Ry     2.7211386259 eV\n"
    )
    assert read_vacuum_level(vacuum) == 5.75
    assert read_fermi_energy(log) == 2.7211386259


def test_work_function_uses_e_vacuum_and_scf_fermi_without_cube(tmp_path, monkeypatch):
    (tmp_path / "E_vacuum.out").write_text("E_VACUUM (eV) = 5.75\n")
    (tmp_path / "running_scf.log").write_text("E_Fermi 0.2000000000 2.7211386259\n")
    monkeypatch.chdir(tmp_path)

    result = task_work_function(
        args=[],
        interactive=False,
        parsed_args=SimpleNamespace(
            file=None, vacuum_file=None, log=None, fermi=None, no_plot=True
        ),
        output_dir=str(tmp_path / "out"),
    )
    assert result is not None
    assert result["vacuum_level"] == 5.75
    assert np.isclose(result["work_function"], 5.75 - 2.7211386259)
    saved = json.loads((tmp_path / "out" / "work_function.json").read_text())
    assert saved["vacuum_source"].endswith("E_vacuum.out")


def test_macroscopic_average_is_periodic_and_double_filtered():
    values = np.sin(np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)) + 3.0
    result = macroscopic_average(values, 5)
    assert result.shape == values.shape
    assert np.isclose(result.mean(), values.mean())
    assert np.max(np.abs(result - 3.0)) < np.max(np.abs(values - 3.0))


def test_macro_task_converts_ry_and_uses_vacuum_and_fermi(tmp_path, monkeypatch):
    profile = np.r_[np.linspace(0.0, 0.2, 8), np.full(8, 0.5)]
    cube = tmp_path / "ElecStaticPot.cube"
    _write_cube(cube, profile)
    (tmp_path / "E_vacuum.out").write_text("E_VACUUM (eV) = 7.0\n")
    (tmp_path / "running_scf.log").write_text("E_Fermi 0.1000000000 1.3605693123\n")
    monkeypatch.chdir(tmp_path)

    result = task_macro_avg_potential(
        args=[],
        interactive=False,
        parsed_args=SimpleNamespace(
            file=str(cube), vacuum_file=None, log=None, fermi=None,
            period=2.0, no_plot=True,
        ),
        output_dir=str(tmp_path / "out"),
    )
    assert result is not None
    assert result["vacuum_level"] == 7.0
    assert np.isclose(result["work_function"], 7.0 - 1.3605693123)
    macro = np.loadtxt(tmp_path / "out" / "macro_avg_potential.dat")[:, 1]
    assert np.isclose(macro.max(), (profile * RY_TO_EV).max(), atol=1.0)
