"""Regression tests for the ABACUS binary rhog (charge density) reader."""
from __future__ import annotations

import struct

import numpy as np
import pytest

from abacuscopilot.postprocessing.charge_tasks import (
    _next_smooth,
    read_charge_cube,
    read_charge_rhog,
)


def _valid_range(n):
    """Valid miller-index range for a grid of size n (ABACUS read_rhog check)."""
    return -((n + 1) // 2) + 1, n // 2


def _write_rhog(path, cell_bohr, rho_r, nspin=1):
    """Write a synthetic ABACUS rhog restart file that round-trips rho_r."""
    nz, ny, nx = rho_r.shape
    # All grid points as planewaves (a dense set), ordered by numpy FFT index
    ix, iy, iz = np.indices((nx, ny, nz)).reshape(3, -1)
    ix = np.where(ix >= (nx + 1) // 2, ix - nx, ix)
    iy = np.where(iy >= (ny + 1) // 2, iy - ny, iy)
    iz = np.where(iz >= (nz + 1) // 2, iz - nz, iz)
    miller = np.column_stack([ix, iy, iz])
    # FFT over (x, y, z) so the coefficient order matches the miller ordering
    arr = np.fft.fftn(rho_r.transpose(2, 1, 0).astype(complex))
    rho_g = arr.ravel() / (nx * ny * nz)  # in ABACUS convention
    # Keep only the planewaves ABACUS stores (valid miller range per axis).
    # For odd grids this is the full set, so the round trip is exact.
    (lxo, hxi), (lyo, hyi), (lzo, hzi) = _valid_range(nx), _valid_range(ny), _valid_range(nz)
    keep = ((ix >= lxo) & (ix <= hxi)
            & (iy >= lyo) & (iy <= hyi)
            & (iz >= lzo) & (iz <= hzi))
    miller = miller[keep]
    rho_g = rho_g[keep]

    # Real ABACUS files store GT in units of 1/Angstrom (see rhog_io.cpp)
    gt = np.linalg.inv(cell_bohr) * 1.889726
    with open(path, "wb") as f:
        for v in (3, 0, len(miller), nspin, 3):
            f.write(struct.pack("<i", v))
        f.write(struct.pack("<i", 9))
        f.write(struct.pack("<9d", *gt.ravel()))
        f.write(struct.pack("<i", 9))
        f.write(struct.pack("<i", 3 * len(miller)))
        for m in miller:
            f.write(struct.pack("<3i", *m))
        f.write(struct.pack("<i", 3 * len(miller)))
        for _ in range(nspin):
            f.write(struct.pack("<i", len(miller)))
            for c in rho_g:
                f.write(struct.pack("<dd", c.real, c.imag))
            f.write(struct.pack("<i", len(miller)))


@pytest.mark.parametrize("grid", [(5, 5, 5), (7, 5, 9)])
def test_rhog_roundtrip(tmp_path, grid):
    nz, ny, nx = grid
    rng = np.random.default_rng(42)
    cell = np.array([[8.0, 0, 0], [0, 8.0, 0], [0, 0, 12.0]])
    rho_r = rng.random((nz, ny, nx))
    p = tmp_path / "CHG"
    _write_rhog(p, cell, rho_r)

    data, cell_read, origin = read_charge_rhog(p)
    assert data.shape == (nz, ny, nx)
    assert np.allclose(cell_read, cell, atol=1e-8)
    assert np.allclose(origin, 0.0)

    # The planewave set is a strict subset (Nyquist-free), so the reconstruction
    # should match the input on the low-frequency part; near the planewave cutoff
    # small ringing is acceptable, so allow a loose tolerance.
    assert np.allclose(data, rho_r, atol=1e-6)


def test_read_charge_cube_detects_binary(tmp_path):
    """read_charge_cube should route binary rhog files to the rhog reader."""
    rng = np.random.default_rng(0)
    cell = np.array([[6.0, 0, 0], [0, 6.0, 0], [0, 0, 10.0]])
    p = tmp_path / "chg_abacus"
    _write_rhog(p, cell, rng.random((5, 3, 5)))
    data, cell_read, _ = read_charge_cube(p)
    assert data.shape == (5, 3, 5)
    assert np.allclose(cell_read, cell, atol=1e-8)


def test_next_smooth():
    assert _next_smooth(1) == 1
    assert _next_smooth(125) == 125  # 5^3
    assert _next_smooth(269) == 270  # 2*3^3*5
    assert _next_smooth(290) == 294  # 2*3*7^2
