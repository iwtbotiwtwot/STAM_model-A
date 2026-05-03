"""Local spherical accumulation functions for STAM_model-A."""

from __future__ import annotations

import numpy as np

from .constants import C_SI, G_SI

ArrayLike = float | np.ndarray


def schwarzschild_radius(mass_kg: ArrayLike) -> ArrayLike:
    """Return Rs = 2GM/c^2 in meters."""
    return 2.0 * G_SI * np.asarray(mass_kg) / C_SI**2


def accumulation_spherical(mass_kg: float, radius_m: ArrayLike) -> ArrayLike:
    """Return A(r) = Rs/r = 2GM/(c^2 r).

    Parameters
    ----------
    mass_kg:
        Source mass in kg.
    radius_m:
        Radius in meters. Must be positive.
    """
    r = np.asarray(radius_m, dtype=float)
    if np.any(r <= 0):
        raise ValueError("radius_m must be positive")
    return schwarzschild_radius(mass_kg) / r


def dA_dr_spherical(mass_kg: float, radius_m: ArrayLike) -> ArrayLike:
    """Return radial derivative dA/dr for A=Rs/r.

    This derivative is negative because A decreases outward.
    """
    r = np.asarray(radius_m, dtype=float)
    if np.any(r <= 0):
        raise ValueError("radius_m must be positive")
    return -schwarzschild_radius(mass_kg) / r**2


def gravity_radial_from_A(mass_kg: float, radius_m: ArrayLike) -> ArrayLike:
    """Return signed radial gravitational acceleration from (c^2/2)dA/dr.

    Negative sign denotes inward acceleration along r_hat.
    """
    return (C_SI**2 / 2.0) * dA_dr_spherical(mass_kg, radius_m)


def gravity_newtonian_radial(mass_kg: float, radius_m: ArrayLike) -> ArrayLike:
    """Return signed Newtonian radial acceleration -GM/r^2."""
    r = np.asarray(radius_m, dtype=float)
    if np.any(r <= 0):
        raise ValueError("radius_m must be positive")
    return -G_SI * mass_kg / r**2


def gravity_magnitude(mass_kg: float, radius_m: ArrayLike) -> ArrayLike:
    """Return |g| = GM/r^2."""
    r = np.asarray(radius_m, dtype=float)
    if np.any(r <= 0):
        raise ValueError("radius_m must be positive")
    return G_SI * mass_kg / r**2
