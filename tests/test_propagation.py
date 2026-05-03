from __future__ import annotations

import math

from stam_model_a.propagation import shapiro_delay_straight_path, shapiro_integral_numeric


def test_shapiro_straight_path_numeric_matches_analytic() -> None:
    mass = 1.98847e30
    impact = 6.9634e8
    x1, x2 = -1.0e11, 1.0e11
    analytic = shapiro_delay_straight_path(mass, impact, x1, x2)
    numeric = shapiro_integral_numeric(mass, impact, x1, x2, n=20001)
    assert math.isclose(numeric, analytic, rel_tol=1e-6)
