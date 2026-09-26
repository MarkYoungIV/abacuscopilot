"""Tests for Structure Editing tasks: 401 Build Supercell and 402 Redefine Lattice.

These guard the matrix-multiplication *order*.  Lattice vectors are the rows
of ``L``, so redefining the cell applies the matrix from the left
(``L' = M @ L``, the VASPKIT TRANSMAT convention) and fractional coordinates
transform as ``f' = f @ M⁻¹``.  With the multiplication on the wrong side the
resulting cell vectors are not integer combinations of the old ones — i.e. not
a supercell at all — while still reporting the correct volume ratio and atom
count, so the mistake is completely silent.
"""

import numpy as np
import pytest

from abacuscopilot.core.models import Atom, Lattice, Structure
from abacuscopilot.preprocessing.structure_editing_tasks import (
    apply_lattice_transform,
    build_supercell,
)

ANG = 1.889726  # Bohr per Angstrom — with this constant `vectors` are in Angstrom

# K3PS4 (Pnma), the cell that exposed the bug.
A, B, C = 9.0336452442, 10.5246413359, 9.1475342167
M_ORTHO = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=float)

# A non-orthogonal (hexagonal-ish) cell, so a diagonal-only shortcut cannot pass.
HEX_CELL = np.array([[5.0, 0.0, 0.0], [-2.5, 4.3301270189, 0.0], [0.0, 0.0, 7.0]], dtype=float)
HEX_POS = [[1.0 / 3.0, 2.0 / 3.0, 0.25]]


def _structure(cell: np.ndarray, positions: list[list[float]]) -> Structure:
    s = Structure()
    s.lattice = Lattice(constant=ANG, vectors=np.array(cell, dtype=float))
    s.coordinate_type = "Direct"
    s.species_order = ["K"]
    s.magnetism = {"K": 0.0}
    s.atoms = [Atom(species="K", position=np.array(p, dtype=float)) for p in positions]
    return s


def _assert_same_crystal(new: Structure, old: Structure, n_image: int) -> None:
    """Every new atom must sit on an original lattice site (modulo the old cell).

    Positions are compared in the OLD cell's fractional basis, so this checks
    the full geometry — cell shape *and* the atom placement — not just the
    volume, which is right either way.
    """
    old_cell = old.lattice.cell_angstrom
    new_cell = new.lattice.cell_angstrom
    expected = [a.position % 1.0 for a in old.atoms]

    assert new.num_atoms == len(expected) * n_image

    matched = [0] * len(expected)
    for atom in new.atoms:
        r = atom.position @ new_cell  # fractional (new) -> Cartesian
        f = (r @ np.linalg.inv(old_cell)) % 1.0
        for i, e in enumerate(expected):
            d = (f - e + 0.5) % 1.0 - 0.5  # minimum image in fractional space
            if np.all(np.abs(d) < 1e-8):
                matched[i] += 1
                break
        else:
            pytest.fail(f"atom at {f} is not on any site of the original crystal")

    assert matched == [n_image] * len(expected)


class TestRedefineLatticeMatrixOrder:
    def test_new_vectors_are_integer_combinations_of_old(self):
        """M's rows give the coefficients: a' = a+b, b' = -a+b, c' = c."""
        old = _structure(np.diag([A, B, C]), [[0.0, 0.0, 0.0]])
        new = apply_lattice_transform(old, M_ORTHO)

        cell_old = old.lattice.cell_angstrom
        cell_new = new.lattice.cell_angstrom
        assert np.allclose(cell_new[0], cell_old[0] + cell_old[1])
        assert np.allclose(cell_new[1], -cell_old[0] + cell_old[1])
        assert np.allclose(cell_new[2], cell_old[2])

        # The coefficients in the old basis must be exactly M — and therefore
        # integers. This is the check the old `vectors @ M` order failed
        # (it produced 0.858 / -1.165).
        coeffs = cell_new @ np.linalg.inv(cell_old)
        assert np.allclose(coeffs, M_ORTHO)
        assert np.allclose(coeffs, np.round(coeffs))

    def test_fractional_coordinates_match_vaspkit(self):
        """Anchor against the VASPKIT SUPERCELL.vasp for the same input."""
        f0 = [0.05273138, 0.54232301, 0.29771140]
        old = _structure(np.diag([A, B, C]), [f0])
        new = apply_lattice_transform(old, M_ORTHO)

        assert new.num_atoms == 2
        assert np.allclose(
            new.atoms[0].position,
            [0.2975271950, 0.2447958150, 0.2977114000],
            atol=1e-8,
        )
        assert np.allclose(
            new.atoms[1].position,
            [0.7975271950, 0.7447958150, 0.2977114000],
            atol=1e-8,
        )

    def test_volume_ratio_and_atom_count(self):
        old = _structure(np.diag([A, B, C]), [[0.1, 0.2, 0.3], [0.6, 0.7, 0.8]])
        new = apply_lattice_transform(old, M_ORTHO)

        assert new.lattice.volume_angstrom == pytest.approx(
            2.0 * old.lattice.volume_angstrom, rel=1e-10
        )
        assert new.num_atoms == 4

    @pytest.mark.parametrize(
        "M",
        [
            M_ORTHO,
            np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
            np.array([[1, 1, 1], [0, 2, 1], [0, 0, 3]], dtype=float),
            np.array([[1, -1, 0], [2, 1, 0], [0, 0, 1]], dtype=float),
        ],
    )
    def test_crystal_is_preserved_for_non_orthogonal_cell(self, M):
        old = _structure(HEX_CELL, HEX_POS)
        new = apply_lattice_transform(old, M)

        n_image = round(float(abs(np.linalg.det(M))))
        _assert_same_crystal(new, old, n_image)

    def test_unimodular_matrix_replicates_nothing(self):
        """A det=±1 change of basis must leave the atom count untouched."""
        M = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=float)  # 90 deg rotation
        old = _structure(HEX_CELL, HEX_POS)
        new = apply_lattice_transform(old, M)

        assert abs(np.linalg.det(M)) == pytest.approx(1.0)
        assert new.num_atoms == old.num_atoms
        _assert_same_crystal(new, old, 1)


class TestBuildSupercell:
    def test_scales_each_lattice_vector(self):
        """Nx x Ny x Nz scales vector i by N_i — not Cartesian component i."""
        old = _structure(HEX_CELL, HEX_POS)
        new = build_supercell(old, 2, 3, 1)

        cell_old = old.lattice.cell_angstrom
        cell_new = new.lattice.cell_angstrom
        assert np.allclose(cell_new[0], 2.0 * cell_old[0])
        assert np.allclose(cell_new[1], 3.0 * cell_old[1])
        assert np.allclose(cell_new[2], 1.0 * cell_old[2])
        assert new.num_atoms == 6
        _assert_same_crystal(new, old, 6)
