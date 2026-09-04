"""Regression tests for reaction pathway tasks (NEB path generation)."""

import numpy as np
from ase import Atoms

from abacuscopilot.preprocessing.reaction_tasks import (
    _compute_max_displacement,
    _linear_interpolate,
)


def test_max_displacement_minimum_image_ignores_image_wrap():
    """A framework atom stored at a different periodic image must not inflate
    the max displacement to ~cell size.

    Regression: two ABACUS endpoints exported to CIF stored an unmoved Mo at a
    periodic image one lattice vector apart, so the raw per-index Cartesian
    difference reported d_max = 11.17 Å (≈ |c|) for a system whose true max
    displacement (migrating Zn) was ~4.36 Å — tripling the suggested NEB image
    count (15 vs 7).
    """
    cell = np.diag([10.0, 10.0, 10.0])
    init = Atoms(
        symbols=["Mo", "Zn"],
        positions=[
            [1.0, 2.0, 3.0],  # framework Mo — physically unmoved
            [1.0, 1.0, 1.0],
        ],  # migrating Zn — real hop
        cell=cell,
        pbc=True,
    )
    final = Atoms(
        symbols=["Mo", "Zn"],
        positions=[
            [11.0, 2.0, 3.0],  # same Mo, stored wrapped by +a
            [5.0, 1.0, 1.0],
        ],  # Zn moved 4 Å
        cell=cell,
        pbc=True,
    )
    d = _compute_max_displacement(init, final)
    assert 3.9 <= d <= 4.1, f"expected ~4 Å physical hop, got {d:.4f}"


def test_max_displacement_nonperiodic_plain_diff():
    """Non-periodic (molecule) input keeps the plain Cartesian difference."""
    init = Atoms(symbols=["C", "O"], positions=[[0, 0, 0], [1.2, 0, 0]], pbc=False)
    final = Atoms(symbols=["C", "O"], positions=[[0, 0, 0], [1.2, 0, 0.5]], pbc=False)
    d = _compute_max_displacement(init, final)
    assert abs(d - 0.5) < 1e-9


def test_linear_interpolate_aligns_wrapped_framework():
    """A framework atom stored one lattice vector apart in the endpoints must
    NOT sweep across the cell during index-wise interpolation.

    Regression: ABACUS relaxations/CIF exports stored an unmoved Mo at a
    different periodic image; linear interpolation then dragged it ~11 Å
    through the framework, giving intermediate NEB images with overlapping
    atoms (Mo–Zn 0.19 Å) that ABACUS rejects as "structure unreasonable".
    """
    cell = np.diag([10.0, 10.0, 10.0])
    init = Atoms(
        symbols=["Mo", "Zn"],
        positions=[
            [1.0, 2.0, 3.0],  # framework Mo
            [1.0, 1.0, 1.0],
        ],  # migrating Zn
        cell=cell,
        pbc=True,
    )
    # same Mo stored wrapped by +a; Zn genuinely hops 4 Å
    final = Atoms(
        symbols=["Mo", "Zn"],
        positions=[[11.0, 2.0, 3.0], [5.0, 1.0, 1.0]],
        cell=cell,
        pbc=True,
    )
    imgs = _linear_interpolate(init, final, 3)
    for i, img in enumerate(imgs):
        mo_x = img.positions[0][0]
        assert abs(mo_x - 1.0) < 1e-6, f"image {i}: framework Mo swept to x={mo_x}"
    # migrating Zn still interpolates linearly 1 → 5
    assert abs(imgs[-1].positions[1][0] - 5.0) < 1e-6
