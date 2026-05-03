"""Propagation/path-integral helpers for STAM_model-A."""

from __future__ import annotations

import numpy as np

from .constants import C_SI, G_SI


def shapiro_delay_straight_path(
    mass_kg: float,
    impact_parameter_m: float,
    x1_m: float,
    x2_m: float,
) -> float:
    """Weak-field accumulation delay for a straight path past a point mass.

    Path parameterization:

        r(x) = sqrt(x^2 + impact_parameter_m^2)

    STAM delay:

        Delta t = (2GM/c^3) int_x1^x2 dx/r(x)
                = (2GM/c^3) [asinh(x/b)]_x1^x2

    This tests the logarithmic/asinh structure only. Coefficient-level comparison
    to Solar System conventions is a separate empirical test.
    """
    if mass_kg <= 0:
        raise ValueError("mass_kg must be positive")
    if impact_parameter_m <= 0:
        raise ValueError("impact_parameter_m must be positive")
    coeff = 2.0 * G_SI * mass_kg / C_SI**3
    integral = np.arcsinh(x2_m / impact_parameter_m) - np.arcsinh(x1_m / impact_parameter_m)
    return float(coeff * integral)


def shapiro_integral_numeric(
    mass_kg: float,
    impact_parameter_m: float,
    x1_m: float,
    x2_m: float,
    n: int = 100_001,
) -> float:
    """Numerically integrate Delta t = (1/c) int A ds for a straight path."""
    if n < 3:
        raise ValueError("n must be at least 3")
    xs = np.linspace(x1_m, x2_m, n)
    r = np.sqrt(xs**2 + impact_parameter_m**2)
    A = 2.0 * G_SI * mass_kg / (C_SI**2 * r)
    return float(np.trapezoid(A, xs) / C_SI)
