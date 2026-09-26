"""Regression tests for elastic constants / Born stability (task 1201)."""

import numpy as np

from abacuscopilot.postprocessing.mechanics_tasks import compute_mechanical_properties


def _cubic_tensor(c11: float, c12: float, c44: float) -> np.ndarray:
    C = np.zeros((6, 6))
    C[:3, :3] = c12
    np.fill_diagonal(C[:3, :3], c11)
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


class TestBornStability:
    def test_cubic_stable_is_positive_definite(self):
        props = compute_mechanical_properties(_cubic_tensor(165.0, 64.0, 80.0))
        assert props["born_stable"] is True
        assert props["born_min_eig"] > 0
        assert len(props["born_eigvals"]) == 6
        assert all(e > 0 for e in props["born_eigvals"])

    def test_hexagonal_stable_though_cubic_test_would_fail(self):
        # Negative C12 gives C11 + 2*C12 < 0: the old cubic-only test declared
        # this matrix unstable, but the elastic matrix is positive definite.
        C = np.zeros((6, 6))
        C[0, 0] = C[1, 1] = C[2, 2] = 100.0
        C[0, 1] = C[1, 0] = -60.0
        C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = 30.0
        C[3, 3] = C[4, 4] = C[5, 5] = 30.0
        assert C[0, 0] + 2 * C[0, 1] < 0

        props = compute_mechanical_properties(C)
        assert props["born_stable"] is True
        assert props["born_min_eig"] > 0

    def test_unstable_tensor_is_detected(self):
        props = compute_mechanical_properties(_cubic_tensor(100.0, 120.0, 10.0))
        assert props["born_stable"] is False
        assert props["born_min_eig"] < 0

    def test_soft_mode_is_a_normalized_six_vector(self):
        props = compute_mechanical_properties(_cubic_tensor(100.0, 120.0, 10.0))
        soft = np.array(props["born_soft_mode"])
        assert soft.shape == (6,)
        assert abs(np.linalg.norm(soft) - 1.0) < 1e-9
        assert props["born_eig_tol"] > 0
