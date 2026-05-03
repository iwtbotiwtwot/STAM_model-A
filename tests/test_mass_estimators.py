from __future__ import annotations

import math

from stam_model_a.constants import C_SI, G_SI
from stam_model_a.local import schwarzschild_radius
from stam_model_a.mass_estimators import (
    mass_from_acceleration,
    mass_from_gravitational_redshift,
    mass_from_horizon_radius,
    mass_from_lensing_deflection,
    mass_from_orbital_velocity,
    mass_from_shapiro_coefficient,
)


def test_synthetic_mass_estimators_round_trip() -> None:
    mass = 5.0e30
    r = 9.0e10
    rs = schwarzschild_radius(mass)
    g = G_SI * mass / r**2
    v = math.sqrt(G_SI * mass / r)
    K = 2.0 * G_SI * mass / C_SI**3
    z_grav = G_SI * mass / (C_SI**2 * r)
    alpha = 4.0 * G_SI * mass / (C_SI**2 * r)

    estimates = [
        mass_from_horizon_radius(rs),
        mass_from_acceleration(g, r),
        mass_from_orbital_velocity(v, r),
        mass_from_shapiro_coefficient(K),
        mass_from_gravitational_redshift(z_grav, r),
        mass_from_lensing_deflection(alpha, r),
    ]
    for estimate in estimates:
        assert math.isclose(estimate / mass, 1.0, rel_tol=1e-12)
