
from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# STAM multi-method mass inference test
# ============================================================
# Purpose:
#   Test whether mass can be inferred consistently through multiple
#   STAM accumulation-field observables.
#
# Core STAM field:
#   A(r) = Rs/r = 2GM/(c^2 r)
#
# Methods tested:
#   1. Horizon threshold:
#        A = 1 -> r = Rs -> M = c^2*r_h/(2G)
#
#   2. Acceleration/gradient:
#        g = (c^2/2)|grad A| = GM/r^2
#        M = g*r^2/G
#
#   3. Orbital velocity:
#        v^2/r = GM/r^2
#        M = v^2*r/G
#
#   4. Shapiro delay coefficient:
#        K = 2GM/c^3
#        M = K*c^3/(2G)
#
#   5. Gravitational redshift, weak-field:
#        z_grav ≈ GM/(c^2 r) = A/2
#        M ≈ z_grav*c^2*r/G
#
#   6. Lensing deflection, weak-field:
#        alpha ≈ 4GM/(c^2 b) = 2A(b)
#        M ≈ alpha*c^2*b/(4G)
#
# Note:
#   This is an internal-consistency test. In synthetic weak-field cases,
#   all routes should recover the same mass. This does not yet prove a
#   new prediction beyond GR/Newtonian weak-field identities.
# ============================================================

BASE = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path(".")
OUT = BASE / "output"
OUT.mkdir(exist_ok=True)

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30

def rs_from_mass(mass_kg):
    return 2.0 * G * mass_kg / C**2

def A_at_r(mass_kg, r_m):
    return rs_from_mass(mass_kg) / r_m

def mass_from_horizon_radius(r_h_m):
    return C**2 * r_h_m / (2.0 * G)

def g_from_mass_r(mass_kg, r_m):
    return G * mass_kg / r_m**2

def mass_from_g_r(g_m_s2, r_m):
    return g_m_s2 * r_m**2 / G

def circular_velocity_from_mass_r(mass_kg, r_m):
    return np.sqrt(G * mass_kg / r_m)

def mass_from_v_r(v_m_s, r_m):
    return v_m_s**2 * r_m / G

def shapiro_K_from_mass(mass_kg):
    # one-way coefficient in seconds
    return 2.0 * G * mass_kg / C**3

def mass_from_shapiro_K(K_s):
    return K_s * C**3 / (2.0 * G)

def grav_redshift_weak_from_mass_r(mass_kg, r_m):
    return G * mass_kg / (C**2 * r_m)

def mass_from_grav_redshift_weak(z_grav, r_m):
    return z_grav * C**2 * r_m / G

def lensing_alpha_weak_from_mass_b(mass_kg, b_m):
    # radians
    return 4.0 * G * mass_kg / (C**2 * b_m)

def mass_from_lensing_alpha_b(alpha_rad, b_m):
    return alpha_rad * C**2 * b_m / (4.0 * G)

cases = [
    {"object": "1 solar mass", "mass_solar": 1.0},
    {"object": "10 solar masses", "mass_solar": 10.0},
    {"object": "Sagittarius A* scale", "mass_solar": 4.3e6},
    {"object": "M87* scale", "mass_solar": 6.5e9},
]

# We test observables at radii safely outside the horizon.
radius_multipliers = [3, 10, 100, 1000, 1_000_000]

rows = []
for case in cases:
    mass_kg = case["mass_solar"] * M_SUN
    Rs = rs_from_mass(mass_kg)

    # Horizon method.
    M_horizon = mass_from_horizon_radius(Rs)
    K_shapiro = shapiro_K_from_mass(mass_kg)
    M_shapiro = mass_from_shapiro_K(K_shapiro)

    for mult in radius_multipliers:
        r = mult * Rs

        A = A_at_r(mass_kg, r)
        grad_A_abs = A / r
        g = g_from_mass_r(mass_kg, r)
        M_g = mass_from_g_r(g, r)

        v = circular_velocity_from_mass_r(mass_kg, r)
        M_v = mass_from_v_r(v, r)

        z_grav = grav_redshift_weak_from_mass_r(mass_kg, r)
        M_z = mass_from_grav_redshift_weak(z_grav, r)

        alpha = lensing_alpha_weak_from_mass_b(mass_kg, r)
        M_alpha = mass_from_lensing_alpha_b(alpha, r)

        rows.append({
            "object": case["object"],
            "mass_input_solar": case["mass_solar"],
            "mass_input_kg": mass_kg,
            "Rs_m": Rs,
            "Rs_km": Rs / 1000.0,
            "test_radius_multiple_Rs": mult,
            "test_radius_m": r,

            "A_Rs_over_r": A,
            "grad_A_abs": grad_A_abs,

            "mass_from_horizon_solar": M_horizon / M_SUN,
            "ratio_horizon": M_horizon / mass_kg,

            "g_m_s2": g,
            "mass_from_acceleration_solar": M_g / M_SUN,
            "ratio_acceleration": M_g / mass_kg,

            "circular_velocity_m_s": v,
            "circular_velocity_fraction_c": v / C,
            "mass_from_orbital_velocity_solar": M_v / M_SUN,
            "ratio_orbital_velocity": M_v / mass_kg,

            "shapiro_K_one_way_s": K_shapiro,
            "mass_from_shapiro_K_solar": M_shapiro / M_SUN,
            "ratio_shapiro": M_shapiro / mass_kg,

            "grav_redshift_weak": z_grav,
            "mass_from_grav_redshift_solar": M_z / M_SUN,
            "ratio_grav_redshift": M_z / mass_kg,

            "lensing_alpha_rad": alpha,
            "lensing_alpha_arcsec": alpha * 206264.806247,
            "mass_from_lensing_solar": M_alpha / M_SUN,
            "ratio_lensing": M_alpha / mass_kg,
        })

results = pd.DataFrame(rows)

ratio_cols = [
    "ratio_horizon",
    "ratio_acceleration",
    "ratio_orbital_velocity",
    "ratio_shapiro",
    "ratio_grav_redshift",
    "ratio_lensing",
]

summary_rows = []
for obj, group in results.groupby("object"):
    entry = {
        "object": obj,
        "mass_input_solar": group["mass_input_solar"].iloc[0],
        "Rs_km": group["Rs_km"].iloc[0],
        "rows_tested": len(group),
    }
    for col in ratio_cols:
        vals = group[col].to_numpy(dtype=float)
        entry[f"{col}_mean"] = float(np.nanmean(vals))
        entry[f"{col}_max_abs_error_from_1"] = float(np.nanmax(np.abs(vals - 1.0)))
    all_ratios = group[ratio_cols].to_numpy(dtype=float).ravel()
    entry["all_methods_mean_ratio"] = float(np.nanmean(all_ratios))
    entry["all_methods_max_abs_error_from_1"] = float(np.nanmax(np.abs(all_ratios - 1.0)))
    summary_rows.append(entry)

summary = pd.DataFrame(summary_rows)

results_path = OUT / "stam_multimethod_mass_inference_results.csv"
summary_path = OUT / "stam_multimethod_mass_inference_summary.csv"

results.to_csv(results_path, index=False)
summary.to_csv(summary_path, index=False)

print("\n=== STAM MULTI-METHOD MASS INFERENCE TEST ===")
print("Core field: A(r) = Rs/r = 2GM/(c^2 r)")
print("Mass inferred from: horizon, acceleration, orbital velocity, Shapiro coefficient, gravitational redshift, lensing.")
print("\n=== SUMMARY ===")
print(summary.to_string(index=False))
print("\nFiles written:")
print(results_path)
print(summary_path)
