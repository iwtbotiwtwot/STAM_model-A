from __future__ import annotations

import numpy as np

from stam_model_a.constants import B_PANTHEON_UNION_STYLE, L_STAM
from stam_model_a.cosmology import (
    A_path_average,
    D_adj,
    D_adj_from_path_integral,
    D_excess,
    D_geo,
    infer_b_from_distance,
)


def test_distance_decomposition() -> None:
    z = np.array([0.1, 0.5, 1.0, 2.0])
    b = B_PANTHEON_UNION_STYLE
    assert np.allclose(D_adj(z, b=b), D_geo(z) + D_excess(z, b=b))


def test_path_integral_identity() -> None:
    z = 1.4
    b = B_PANTHEON_UNION_STYLE
    assert np.isclose(D_adj(z, b=b), D_adj_from_path_integral(z, b=b, n=20001), rtol=1e-10)


def test_average_path_accumulation_identity() -> None:
    z = np.array([0.1, 0.5, 1.0])
    b = B_PANTHEON_UNION_STYLE
    assert np.allclose(A_path_average(z, b=b), D_excess(z, b=b) / D_geo(z))


def test_infer_b_from_distance_round_trip() -> None:
    z = np.array([0.1, 0.5, 1.0])
    b = 461.3626922
    d = D_adj(z, b=b, L=L_STAM)
    assert np.allclose(infer_b_from_distance(z, d, L=L_STAM), b)
