from __future__ import annotations

import numpy as np

from stam_model_a.local import (
    accumulation_spherical,
    gravity_newtonian_radial,
    gravity_radial_from_A,
    schwarzschild_radius,
)


def test_gravity_identity_spherical() -> None:
    mass = 3.0e30
    radii = np.array([1.0e8, 2.0e9, 4.0e10])
    assert np.allclose(gravity_radial_from_A(mass, radii), gravity_newtonian_radial(mass, radii))


def test_horizon_threshold() -> None:
    mass = 3.0e30
    rs = schwarzschild_radius(mass)
    assert np.isclose(accumulation_spherical(mass, rs), 1.0)
