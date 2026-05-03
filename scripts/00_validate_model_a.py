#!/usr/bin/env python3
"""Run baseline STAM_model-A identity/smoke checks."""

from __future__ import annotations

from pathlib import Path
import json
import math

import numpy as np

from stam_model_a.constants import C_SI, G_SI, L_STAM, B_PANTHEON_UNION_STYLE
from stam_model_a.local import (
    schwarzschild_radius,
    accumulation_spherical,
    gravity_radial_from_A,
    gravity_newtonian_radial,
)
from stam_model_a.propagation import shapiro_delay_straight_path, shapiro_integral_numeric
from stam_model_a.cosmology import D_adj, D_adj_from_path_integral, D_excess, D_geo
from stam_model_a.mass_estimators import (
    mass_from_acceleration,
    mass_from_gravitational_redshift,
    mass_from_horizon_radius,
    mass_from_lensing_deflection,
    mass_from_orbital_velocity,
    mass_from_shapiro_coefficient,
)


def assert_close(name: str, got: float, expected: float, rtol: float = 1e-10, atol: float = 0.0) -> None:
    if not math.isclose(got, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(f"{name}: got {got}, expected {expected}")


def main() -> None:
    results: dict[str, object] = {}

    # 1. Local gravity identity.
    mass = 1.98847e30  # kg, approximate solar mass
    radii = np.array([1.0e8, 1.0e9, 1.0e10, 1.0e11])
    g_a = gravity_radial_from_A(mass, radii)
    g_n = gravity_newtonian_radial(mass, radii)
    max_rel = float(np.max(np.abs((g_a - g_n) / g_n)))
    if max_rel > 1e-14:
        raise AssertionError(f"local gravity identity failed: max_rel={max_rel}")
    results["local_gravity_identity_max_rel_error"] = max_rel

    # 2. Horizon threshold A=1 at r=Rs.
    rs = float(schwarzschild_radius(mass))
    A_at_rs = float(accumulation_spherical(mass, rs))
    assert_close("A_at_Rs", A_at_rs, 1.0, rtol=1e-14)
    results["A_at_Rs"] = A_at_rs

    # 3. Weak-field Shapiro integral: analytic asinh form matches numerical integral.
    impact = 6.9634e8  # m, rough solar radius
    x1, x2 = -1.0e11, 1.0e11
    analytic_delay = shapiro_delay_straight_path(mass, impact, x1, x2)
    numeric_delay = shapiro_integral_numeric(mass, impact, x1, x2, n=20001)
    rel_delay = abs(analytic_delay - numeric_delay) / analytic_delay
    if rel_delay > 1e-6:
        raise AssertionError(f"Shapiro integral check failed: rel={rel_delay}")
    results["shapiro_analytic_seconds"] = analytic_delay
    results["shapiro_numeric_seconds"] = numeric_delay
    results["shapiro_rel_error"] = rel_delay

    # 4. Mass estimators recover synthetic mass.
    r = 4.2e10
    g = G_SI * mass / r**2
    v = math.sqrt(G_SI * mass / r)
    K = 2.0 * G_SI * mass / C_SI**3
    z_grav = G_SI * mass / (C_SI**2 * r)
    alpha = 4.0 * G_SI * mass / (C_SI**2 * r)
    estimates = {
        "horizon": mass_from_horizon_radius(rs),
        "acceleration": mass_from_acceleration(g, r),
        "orbital_velocity": mass_from_orbital_velocity(v, r),
        "shapiro_coefficient": mass_from_shapiro_coefficient(K),
        "gravitational_redshift": mass_from_gravitational_redshift(z_grav, r),
        "lensing_deflection": mass_from_lensing_deflection(alpha, r),
    }
    ratios = {name: estimate / mass for name, estimate in estimates.items()}
    max_mass_rel = max(abs(ratio - 1.0) for ratio in ratios.values())
    if max_mass_rel > 1e-12:
        raise AssertionError(f"mass estimator check failed: ratios={ratios}")
    results["mass_estimator_ratios"] = ratios

    # 5. Cosmology identities.
    z = 1.25
    b = B_PANTHEON_UNION_STYLE
    d_adj = float(D_adj(z, b=b, L=L_STAM))
    d_sum = float(D_geo(z, L=L_STAM) + D_excess(z, b=b, L=L_STAM))
    d_int = D_adj_from_path_integral(z, b=b, L=L_STAM, n=20001)
    assert_close("D_adj_sum", d_adj, d_sum, rtol=1e-13)
    rel_int = abs(d_adj - d_int) / d_adj
    if rel_int > 1e-10:
        raise AssertionError(f"cosmology path integral failed: rel={rel_int}")
    results["cosmology_D_adj"] = d_adj
    results["cosmology_path_integral_rel_error"] = rel_int

    out_dir = Path("results/identity_checks")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "validate_model_a_results.json"
    out_path.write_text(json.dumps(results, indent=2, sort_keys=True))

    print("PASS: STAM_model-A baseline identity checks passed.")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
