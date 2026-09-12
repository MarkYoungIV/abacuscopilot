"""Optional end-to-end check against a completed ABACUS SCF case.

Set ``ABACUS_WORKFUNC_CASE`` to a case directory containing
``OUT.ABACUS/ElecStaticPot.cube`` and ``OUT.ABACUS/running_scf.log`` before
running this file.  The large calculation output is intentionally not stored
in the source repository.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from abacuscopilot.postprocessing.workfunc_tasks import (
    _read_profile,
    estimate_vacuum_level_from_cube,
    macroscopic_average,
    read_fermi_energy,
)


def _real_case_paths() -> tuple[Path, Path]:
    root = os.environ.get("ABACUS_WORKFUNC_CASE")
    if not root:
        pytest.skip("Set ABACUS_WORKFUNC_CASE to run the real ABACUS integration test")
    case = Path(root)
    cube = case / "OUT.ABACUS" / "ElecStaticPot.cube"
    log = case / "OUT.ABACUS" / "running_scf.log"
    if not cube.is_file() or not log.is_file():
        pytest.skip("ABACUS_WORKFUNC_CASE lacks ElecStaticPot.cube or running_scf.log")
    return cube, log


def test_completed_abacus_scf_supports_raw_cube_workfunc_pipeline():
    cube, log = _real_case_paths()
    profile = _read_profile(cube)
    assert profile is not None
    z_ang, planar, dz, shape, atom_z, cell_length = profile
    fermi = read_fermi_energy(log)
    assert fermi is not None
    assert shape[2] > 1 and atom_z.size > 0

    raw_vacuum = estimate_vacuum_level_from_cube(
        planar, z_ang, atom_z, cell_length, exclude_distance=3.0
    )
    assert raw_vacuum is not None
    assert np.isfinite(raw_vacuum["vacuum_level"])
    assert np.isfinite(raw_vacuum["vacuum_level"] - fermi)

    window = max(1, int(round(3.0 / dz)))
    macro = macroscopic_average(planar, window)
    macro_vacuum = estimate_vacuum_level_from_cube(
        macro, z_ang, atom_z, cell_length, exclude_distance=3.0
    )
    assert macro_vacuum is not None
    assert np.isfinite(macro_vacuum["vacuum_level"] - fermi)
