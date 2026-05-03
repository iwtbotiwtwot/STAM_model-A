"""Operational mass estimators used in STAM_model-A tests."""

from __future__ import annotations

from .constants import C_SI, G_SI


def mass_from_horizon_radius(r_h_m: float) -> float:
    """M = c^2 r_h / (2G)."""
    return C_SI**2 * r_h_m / (2.0 * G_SI)


def mass_from_acceleration(g_m_s2: float, r_m: float) -> float:
    """M = |g| r^2 / G."""
    return abs(g_m_s2) * r_m**2 / G_SI


def mass_from_orbital_velocity(v_m_s: float, r_m: float) -> float:
    """M = v^2 r / G."""
    return v_m_s**2 * r_m / G_SI


def mass_from_shapiro_coefficient(K_s: float) -> float:
    """M = K c^3 / (2G), where K is the coefficient multiplying int(ds/r)."""
    return K_s * C_SI**3 / (2.0 * G_SI)


def mass_from_gravitational_redshift(z_grav: float, r_m: float) -> float:
    """Weak-field estimate M ~= z_grav c^2 r / G."""
    return z_grav * C_SI**2 * r_m / G_SI


def mass_from_lensing_deflection(alpha_rad: float, impact_parameter_m: float) -> float:
    """Thin-lens weak-field estimate M ~= alpha c^2 b / (4G)."""
    return alpha_rad * C_SI**2 * impact_parameter_m / (4.0 * G_SI)
