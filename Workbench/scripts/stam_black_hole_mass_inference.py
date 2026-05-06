
from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# STAM black-hole mass inference test
# ============================================================
# STAM threshold:
#   A(r) = Rs/r
#   A = 1 -> r = Rs
#
# Schwarzschild radius:
#   Rs = 2GM/c^2
#
# Therefore mass from horizon radius:
#   M = c^2 * r_h / (2G)
#
# STAM gravity bridge:
#   g = (c^2/2) grad(A)
#
# In circular orbit limit:
#   v^2/r = GM/r^2
#   M = v^2 r / G
#
# Shapiro one-way coefficient:
#   K = 2GM/c^3
#   M = K c^3 / (2G)
#
# This script checks whether the STAM threshold relation and other
# gravitational observables infer consistent black-hole masses.
# ============================================================

BASE = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path(".")
OUT = BASE / "output"
OUT.mkdir(exist_ok=True)

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30

def rs_from_mass(mass_kg):
    return 2.0 * G * mass_kg / C**2

def mass_from_threshold_radius(radius_m):
    return C**2 * radius_m / (2.0 * G)

def mass_from_orbital_velocity(v_m_s, r_m):
    return v_m_s**2 * r_m / G

def circular_velocity_from_mass(mass_kg, r_m):
    return np.sqrt(G * mass_kg / r_m)

def shapiro_coefficient_from_mass(mass_kg):
    # One-way coefficient K = 2GM/c^3 = Rs/c
    return 2.0 * G * mass_kg / C**3

def mass_from_shapiro_coefficient(K_seconds):
    return K_seconds * C**3 / (2.0 * G)

# Illustrative cases.
# These are not a fitted observational catalog. They are sanity checks for the formula.
cases = [
    {"object": "1 solar mass threshold", "mass_solar_input": 1.0},
    {"object": "10 solar mass threshold", "mass_solar_input": 10.0},
    {"object": "Sagittarius A* scale example", "mass_solar_input": 4.3e6},
    {"object": "M87* scale example", "mass_solar_input": 6.5e9},
]

rows = []
for case in cases:
    mass_input_kg = case["mass_solar_input"] * M_SUN
    rs_m = rs_from_mass(mass_input_kg)
    rs_km = rs_m / 1000.0

    # STAM threshold mass from horizon radius.
    mass_from_rs_kg = mass_from_threshold_radius(rs_m)

    # Test circular orbit at 10 Rs and 100 Rs.
    for orbit_mult in [3, 10, 100, 1000]:
        orbit_r_m = orbit_mult * rs_m
        v = circular_velocity_from_mass(mass_input_kg, orbit_r_m)
        mass_from_orbit_kg = mass_from_orbital_velocity(v, orbit_r_m)

        rows.append({
            "object": case["object"],
            "mass_solar_input": case["mass_solar_input"],
            "mass_input_kg": mass_input_kg,
            "stam_threshold_A": 1.0,
            "horizon_radius_Rs_m": rs_m,
            "horizon_radius_Rs_km": rs_km,
            "mass_from_STAM_threshold_kg": mass_from_rs_kg,
            "mass_from_STAM_threshold_solar": mass_from_rs_kg / M_SUN,
            "threshold_mass_ratio_to_input": mass_from_rs_kg / mass_input_kg,
            "orbit_radius_multiple_Rs": orbit_mult,
            "orbit_radius_m": orbit_r_m,
            "circular_velocity_m_s": v,
            "circular_velocity_fraction_c": v / C,
            "mass_from_orbital_velocity_kg": mass_from_orbit_kg,
            "mass_from_orbital_velocity_solar": mass_from_orbit_kg / M_SUN,
            "orbital_mass_ratio_to_input": mass_from_orbit_kg / mass_input_kg,
            "shapiro_K_one_way_seconds": shapiro_coefficient_from_mass(mass_input_kg),
            "mass_from_shapiro_K_kg": mass_from_shapiro_coefficient(shapiro_coefficient_from_mass(mass_input_kg)),
            "mass_from_shapiro_K_solar": mass_from_shapiro_coefficient(shapiro_coefficient_from_mass(mass_input_kg)) / M_SUN,
            "shapiro_mass_ratio_to_input": mass_from_shapiro_coefficient(shapiro_coefficient_from_mass(mass_input_kg)) / mass_input_kg,
        })

results = pd.DataFrame(rows)

summary = results.groupby("object").agg(
    mass_solar_input=("mass_solar_input", "first"),
    horizon_radius_Rs_km=("horizon_radius_Rs_km", "first"),
    mass_from_STAM_threshold_solar=("mass_from_STAM_threshold_solar", "first"),
    threshold_mass_ratio_to_input=("threshold_mass_ratio_to_input", "first"),
    shapiro_K_one_way_seconds=("shapiro_K_one_way_seconds", "first"),
    mass_from_shapiro_K_solar=("mass_from_shapiro_K_solar", "first"),
    shapiro_mass_ratio_to_input=("shapiro_mass_ratio_to_input", "first"),
    orbital_mass_ratio_mean=("orbital_mass_ratio_to_input", "mean"),
    orbital_mass_ratio_max_error=("orbital_mass_ratio_to_input", lambda x: float(np.max(np.abs(x - 1.0)))),
).reset_index()

results_path = OUT / "stam_black_hole_mass_inference_results.csv"
summary_path = OUT / "stam_black_hole_mass_inference_summary.csv"

results.to_csv(results_path, index=False)
summary.to_csv(summary_path, index=False)

print("\n=== STAM BLACK-HOLE MASS INFERENCE TEST ===")
print("Core formulas:")
print("  A(r) = Rs/r")
print("  A = 1 -> r = Rs")
print("  M = c^2 * Rs / (2G)")
print("  M_orbit = v^2*r/G")
print("  K_shapiro = 2GM/c^3 -> M = K*c^3/(2G)")
print("\n=== SUMMARY ===")
print(summary.to_string(index=False))
print("\nFiles written:")
print(results_path)
print(summary_path)
