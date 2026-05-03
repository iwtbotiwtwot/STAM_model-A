"""Cosmological distance and path-accumulation functions for STAM_model-A."""

from __future__ import annotations

import numpy as np

from .constants import C_MLY_PER_MPC, L_STAM

ArrayLike = float | np.ndarray


def _z_array(z: ArrayLike) -> np.ndarray:
    arr = np.asarray(z, dtype=float)
    if np.any(arr < 0):
        raise ValueError("redshift z must be nonnegative")
    return arr


def D_geo(z: ArrayLike, L: float = L_STAM) -> ArrayLike:
    """Geometric-distance proxy D_geo = L z (1 + 0.15 z)."""
    z_arr = _z_array(z)
    return L * z_arr * (1.0 + 0.15 * z_arr)


def D_adj(z: ArrayLike, b: float, L: float = L_STAM) -> ArrayLike:
    """Adjusted/luminosity-distance-equivalent form D_adj = Lz(1+0.5z)+bz."""
    z_arr = _z_array(z)
    return L * z_arr * (1.0 + 0.5 * z_arr) + b * z_arr


def D_excess(z: ArrayLike, b: float, L: float = L_STAM) -> ArrayLike:
    """Traversal excess D_excess = b z + 0.35 L z^2."""
    z_arr = _z_array(z)
    return b * z_arr + 0.35 * L * z_arr**2


def A_path_average(z: ArrayLike, b: float, L: float = L_STAM) -> ArrayLike:
    """Average path accumulation <A_path> = D_excess/D_geo.

    At z=0 this ratio is undefined; the limiting value is b/L.
    This function returns b/L exactly for z=0.
    """
    z_arr = _z_array(z)
    return (b / L + 0.35 * z_arr) / (1.0 + 0.15 * z_arr)


def A_path_local(z: ArrayLike, b: float, L: float = L_STAM) -> ArrayLike:
    """Local differential path accumulation dD_excess/dD_geo."""
    z_arr = _z_array(z)
    return (b / L + 0.70 * z_arr) / (1.0 + 0.30 * z_arr)


def dD_geo_dz(z: ArrayLike, L: float = L_STAM) -> ArrayLike:
    """Derivative of D_geo with respect to redshift."""
    z_arr = _z_array(z)
    return L * (1.0 + 0.30 * z_arr)


def D_adj_from_path_integral(z: float, b: float, L: float = L_STAM, n: int = 10_001) -> float:
    """Numerically evaluate int_0^z [1 + A_path_local(u)] dD_geo/du du."""
    if z < 0:
        raise ValueError("redshift z must be nonnegative")
    if n < 3:
        raise ValueError("n must be at least 3")
    zs = np.linspace(0.0, z, n)
    integrand = (1.0 + A_path_local(zs, b=b, L=L)) * dD_geo_dz(zs, L=L)
    return float(np.trapezoid(integrand, zs))


def infer_b_from_distance(z: ArrayLike, D_obs: ArrayLike, L: float = L_STAM) -> ArrayLike:
    """Infer b from observed adjusted distance D_obs using D_adj = Lz(1+0.5z)+bz.

    z must be positive. This is useful for redshift-bin falsification tests.
    """
    z_arr = np.asarray(z, dtype=float)
    d_arr = np.asarray(D_obs, dtype=float)
    if np.any(z_arr <= 0):
        raise ValueError("z must be positive to infer b")
    return (d_arr - L * z_arr * (1.0 + 0.5 * z_arr)) / z_arr


def distance_modulus_to_mly(mu: ArrayLike) -> ArrayLike:
    """Convert distance modulus to luminosity distance in million light-years.

    Standard relation:

        D_Mpc = 10 ** ((mu - 25)/5)
        D_Mly = D_Mpc * 3.261563776
    """
    mu_arr = np.asarray(mu, dtype=float)
    d_mpc = 10.0 ** ((mu_arr - 25.0) / 5.0)
    return d_mpc * C_MLY_PER_MPC
