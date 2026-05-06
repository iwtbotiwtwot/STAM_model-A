from pathlib import Path
import math
import numpy as np
import pandas as pd

# ============================================================
# STAM Shapiro delay test
# ============================================================
# STAM premise:
#   A(r) = Rs/r = 2GM/(c^2 r)
#
# Path delay:
#   Δt_STAM = (1/c) ∫ A(r) ds
#            = (Rs/c) ∫ ds/r
#            = (2GM/c^3) ∫ ds/r
#
# For a near-straight path with impact parameter b:
#   r(s) = sqrt(s^2 + b^2)
#   ∫ ds/r = asinh(s2/b) + asinh(s1/b)
#
# This matches the weak-field GR Shapiro one-way delay:
#   Δt_GR = (2GM/c^3) ln((r1+r2+R)/(r1+r2-R))
#
# Outputs:
#   output/stam_shapiro_delay_test_results.csv
# ============================================================

BASE = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path(".")
OUT = BASE / "output"
OUT.mkdir(exist_ok=True)

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30
R_SUN = 696_340_000.0
AU = 149_597_870_700.0

def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2

def stam_shapiro_oneway(mass_kg: float, r1_m: float, r2_m: float, b_m: float) -> tuple[float, float]:
    """One-way STAM Shapiro delay from ∫A ds / c."""
    x1 = math.sqrt(max(r1_m*r1_m - b_m*b_m, 0.0))
    x2 = math.sqrt(max(r2_m*r2_m - b_m*b_m, 0.0))
    integral_ds_over_r = math.asinh(x1 / b_m) + math.asinh(x2 / b_m)
    delay_s = schwarzschild_radius(mass_kg) / C * integral_ds_over_r
    return delay_s, integral_ds_over_r

def gr_shapiro_oneway(mass_kg: float, r1_m: float, r2_m: float, b_m: float) -> float:
    """Standard weak-field one-way Shapiro delay using straight-line geometry."""
    x1 = math.sqrt(max(r1_m*r1_m - b_m*b_m, 0.0))
    x2 = math.sqrt(max(r2_m*r2_m - b_m*b_m, 0.0))
    R = x1 + x2
    return 2.0 * G * mass_kg / C**3 * math.log((r1_m + r2_m + R) / (r1_m + r2_m - R))

def gr_large_distance_approx_oneway(mass_kg: float, r1_m: float, r2_m: float, b_m: float) -> float:
    """Common large-distance approximation."""
    return 2.0 * G * mass_kg / C**3 * math.log(4.0 * r1_m * r2_m / (b_m*b_m))

def main():
    print("\n=== STAM SHAPIRO DELAY TEST ===")
    print("Testing: Δt_STAM = (1/c)∫A(r)ds with A(r)=Rs/r")
    print("Expected weak-field match: Δt_GR = (2GM/c^3)∫ds/r")

    cases = [
        ("Sun: Earth-Mercury grazing", AU, 0.387 * AU, R_SUN),
        ("Sun: Earth-Venus grazing", AU, 0.723 * AU, R_SUN),
        ("Sun: Earth-Mars grazing", AU, 1.524 * AU, R_SUN),
        ("Sun: Earth-Saturn grazing", AU, 9.58 * AU, R_SUN),
        ("Sun: Earth-Cassini-ish 8.43 AU grazing", AU, 8.43 * AU, R_SUN),
        ("Sun: Earth-Mercury b=5Rsun", AU, 0.387 * AU, 5 * R_SUN),
        ("Sun: Earth-Mars b=10Rsun", AU, 1.524 * AU, 10 * R_SUN),
    ]

    rows = []
    for name, r1, r2, b in cases:
        t_stam, integral = stam_shapiro_oneway(M_SUN, r1, r2, b)
        t_gr = gr_shapiro_oneway(M_SUN, r1, r2, b)
        t_approx = gr_large_distance_approx_oneway(M_SUN, r1, r2, b)

        rows.append({
            "case": name,
            "r1_AU": r1 / AU,
            "r2_AU": r2 / AU,
            "impact_b_Rsun": b / R_SUN,
            "integral_ds_over_r": integral,
            "STAM_oneway_s": t_stam,
            "GR_oneway_s": t_gr,
            "large_distance_approx_oneway_s": t_approx,
            "STAM_minus_GR_s": t_stam - t_gr,
            "relative_error_STAM_vs_GR": (t_stam - t_gr) / t_gr,
            "STAM_roundtrip_us": 2.0 * t_stam * 1e6,
            "GR_roundtrip_us": 2.0 * t_gr * 1e6,
        })

    df = pd.DataFrame(rows)
    out_path = OUT / "stam_shapiro_delay_test_results.csv"
    df.to_csv(out_path, index=False)

    print("\nResults:")
    print(df.to_string(index=False))
    print("\nWritten:")
    print(out_path)
    print("\nInterpretation:")
    print("If A(r)=Rs/r, then STAM's path delay ∫A ds/c matches the weak-field GR Shapiro delay.")
    print("This is an internal-consistency bridge, not a new deviation from GR yet.")

if __name__ == "__main__":
    main()
