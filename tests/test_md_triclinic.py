"""Regression tests for triclinic-cell handling in the MD analysis chain.

Two defects fixed in v0.1.38:

* task 3107 (LAMMPS dump -> ABACUS MD_dump) treated the LAMMPS *bounding-box*
  values as true box lengths and dropped the tilt columns, writing an
  orthogonal cell for a triclinic box;
* ``_minimum_image`` used ``inv(cell.T)`` / ``@ cell.T`` even though the whole
  codebase uses the row-vector convention ``cart = frac @ cell`` (see
  ``symmetry_tasks.py`` and ``reaction_tasks.py``), so wrapping happened in
  the wrong space for any non-orthogonal cell.

The MIC is additionally checked against a brute-force nearest-image search,
which is the ground truth used by ASE's ``get_all_distances(mic=True)``.
"""

import numpy as np

from abacuscopilot.postprocessing.md_tasks import (
    _lammps_bounds_to_cell,
    _minimum_image,
    _pairwise_dist_mic,
    task_lammps_to_md_dump,
)

MILD = np.array([[9.9, 0.0, 0.0], [0.3, 11.6, 0.0], [-0.8, -0.4, 12.0]])
SHEARED = np.array([[10.0, 0.0, 0.0], [2.5, 9.0, 0.0], [3.1, 1.7, 11.0]])


def _brute_force_mic(dr: np.ndarray, cell: np.ndarray, r: int = 2) -> np.ndarray:
    """Shortest distance over neighbouring images (ground truth)."""
    best = np.full(len(dr), np.inf)
    for i in range(-r, r + 1):
        for j in range(-r, r + 1):
            for k in range(-r, r + 1):
                shift = i * cell[0] + j * cell[1] + k * cell[2]
                np.minimum(best, np.linalg.norm(dr + shift, axis=1), out=best)
    return best


class TestMinimumImage:
    def test_triclinic_matches_brute_force(self):
        rng = np.random.default_rng(0)
        dr = rng.uniform(-0.5, 0.5, size=(5000, 3)) @ MILD
        got = np.linalg.norm(_minimum_image(dr, MILD), axis=1)
        assert np.allclose(got, _brute_force_mic(dr, MILD), atol=1e-10)

    def test_highly_sheared_matches_brute_force(self):
        rng = np.random.default_rng(1)
        dr = rng.uniform(-0.5, 0.5, size=(5000, 3)) @ SHEARED
        got = np.linalg.norm(_minimum_image(dr, SHEARED), axis=1)
        assert np.allclose(got, _brute_force_mic(dr, SHEARED), atol=1e-10)

    def test_orthogonal_is_exact(self):
        cell = np.diag([10.0, 11.0, 12.0])
        dr = np.array([[4.9, -6.5, 0.3], [-5.4, 6.2, -1.0]])
        got = _minimum_image(dr, cell)
        assert np.allclose(got, [[4.9, 4.5, 0.3], [4.6, -4.8, -1.0]])

    def test_single_vector_shape_preserved(self):
        got = _minimum_image(np.array([4.9, 0.0, 0.0]), np.diag([10.0, 11.0, 12.0]))
        assert got.shape == (3,)

    def test_buggy_transpose_would_differ(self):
        """Guard: the old inv(cell.T)/@cell.T form disagrees on a triclinic cell."""
        rng = np.random.default_rng(2)
        dr = rng.uniform(-0.5, 0.5, size=(2000, 3)) @ MILD
        inv_t = np.linalg.inv(MILD.T)
        frac = dr @ inv_t
        frac -= np.floor(frac + 0.5)
        buggy = frac @ MILD.T
        assert not np.allclose(buggy, _minimum_image(dr, MILD), atol=1e-6)


class TestLammpsBoundsToCell:
    def test_triclinic(self):
        bounds = [(-0.8, 10.2, 0.3), (-0.4, 11.6, -0.8), (0.0, 12.0, -0.4)]
        cell = _lammps_bounds_to_cell(bounds)
        assert np.allclose(cell, MILD)
        assert abs(abs(np.linalg.det(cell)) - 9.9 * 11.6 * 12.0) < 1e-9

    def test_orthorhombic(self):
        bounds = [(0.0, 10.0, 0.0), (0.0, 11.0, 0.0), (0.0, 12.0, 0.0)]
        assert np.allclose(_lammps_bounds_to_cell(bounds), np.diag([10.0, 11.0, 12.0]))

    def test_positive_tilt(self):
        # positive tilt: xlo_bound must be shifted down by xy+xz, not 0
        bounds = [(2.0, 12.0, 0.5), (1.0, 12.0, 0.7), (0.0, 10.0, -0.3)]
        cell = _lammps_bounds_to_cell(bounds)
        assert np.allclose(cell[1, 0], 0.5)   # b_x = xy
        assert np.allclose(cell[2, 0], 0.7)   # c_x = xz
        assert np.allclose(cell[2, 1], -0.3)  # c_y = yz


class TestPairwiseDistMic:
    def test_matches_brute_force(self):
        rng = np.random.default_rng(3)
        pos_a = rng.uniform(0, 1, size=(40, 3)) @ MILD
        pos_b = rng.uniform(0, 1, size=(37, 3)) @ MILD
        got = _pairwise_dist_mic(pos_a, pos_b, MILD, np.linalg.inv(MILD))
        dr = (pos_a[:, None, :] - pos_b[None, :, :]).reshape(-1, 3)
        # (no PBC wrap needed: brute force over images)
        ref = _brute_force_mic(dr, MILD).reshape(40, 37)
        assert np.allclose(got, ref, atol=1e-10)


class TestLammpsDumpConversion:
    DUMP = (
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n4\n"
        "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
        "-0.8 10.2 0.3\n-0.4 11.6 -0.8\n0.0 12.0 -0.4\n"
        "ITEM: ATOMS id type x y z\n"
        "1 1 1.0 1.0 1.0\n2 2 2.0 2.0 2.0\n"
        "3 1 3.0 3.0 3.0\n4 2 4.0 4.0 4.0\n"
    )

    def test_triclinic_cell_written(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "min.dump").write_text(self.DUMP)
        task_lammps_to_md_dump(args=["min.dump"], interactive=False)

        text = (tmp_path / "min.dump.ABACUS.dump").read_text()
        lines = text.splitlines()
        idx = lines.index("LATTICE_VECTORS")
        cell = np.array([[float(x) for x in lines[idx + 1 + j].split()] for j in range(3)])
        assert np.allclose(cell, MILD)

    def test_orthorhombic_dump_unchanged(self, tmp_path, monkeypatch):
        dump = (
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
            "ITEM: BOX BOUNDS pp pp pp\n"
            "0.0 10.0\n0.0 11.0\n0.0 12.0\n"
            "ITEM: ATOMS id type x y z\n"
            "1 1 1.0 1.0 1.0\n2 1 2.0 2.0 2.0\n"
        )
        monkeypatch.chdir(tmp_path)
        (tmp_path / "o.dump").write_text(dump)
        task_lammps_to_md_dump(args=["o.dump"], interactive=False)

        text = (tmp_path / "o.dump.ABACUS.dump").read_text()
        lines = text.splitlines()
        idx = lines.index("LATTICE_VECTORS")
        cell = np.array([[float(x) for x in lines[idx + 1 + j].split()] for j in range(3)])
        assert np.allclose(cell, np.diag([10.0, 11.0, 12.0]))
