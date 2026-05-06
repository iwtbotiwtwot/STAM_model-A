#!/usr/bin/env python3
"""
SF1_A_metric_schwarzschild_recovery_test.py

Accumulation Theory / STAM Model-A
SF1 — A-metric Schwarzschild recovery test

Purpose:
    Establish a disciplined strong-field baseline by rewriting the Schwarzschild
    metric in A-language.

Core definition:
    A(r) = Rs/r = 2GM/(c^2 r)

Candidate strong-field A-metric:
    ds^2 = -(1 - A)c^2 dt^2 + dr^2/(1 - A) + r^2 dΩ^2

This is Schwarzschild written in A notation, not a new strong-field metric yet.

What SF1 checks:
    1. Weak-field clock:
        sqrt(1-A) ≈ 1 - A/2

    2. Horizon:
        A = 1 -> F(A)=1-A=0

    3. Escape-speed identity:
        A = v_escape^2/c^2

    4. Strong-field landmarks:
        photon sphere r = 1.5 Rs -> A = 2/3
        ISCO          r = 3.0 Rs -> A = 1/3

    5. Radial null coordinate speed:
        dr/dt = c(1-A)

    6. Radial coordinate light-time excess:
        Δt_exact = ∫ dr/[c(1-A)] - ∫dr/c
                 = (Rs/c) ln((r2-Rs)/(r1-Rs))
        weak approximation:
        Δt_weak ≈ ∫ A dr/c = (Rs/c) ln(r2/r1)

This script is an anchor test:
    STAM must recover known Schwarzschild structure before proposing modifications.
"""

from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Normalized units:
# c = 1, Rs = 1.
# Then A = 1/r.
C = 1.0
RS = 1.0


def A_of_r(r):
    return RS / np.asarray(r, dtype=float)


def F_of_A(A):
    return 1.0 - np.asarray(A, dtype=float)


def G_of_A(A):
    return 1.0 / (1.0 - np.asarray(A, dtype=float))


def clock_factor(A):
    return np.sqrt(np.maximum(1.0 - np.asarray(A, dtype=float), 0.0))


def clock_weak(A):
    return 1.0 - 0.5*np.asarray(A, dtype=float)


def escape_v_over_c(A):
    return np.sqrt(np.asarray(A, dtype=float))


def radial_null_speed_over_c(A):
    # For ds^2=0 radial: 0 = -(1-A)c^2 dt^2 + dr^2/(1-A)
    # dr/dt = c(1-A)
    return 1.0 - np.asarray(A, dtype=float)


def exact_radial_excess_time(r1, r2):
    # Coordinate time excess over flat radial light travel in normalized units.
    # t_exact = ∫ dr/(1 - Rs/r) = ∫ r/(r-Rs) dr = (r2-r1)+Rs ln((r2-Rs)/(r1-Rs))
    # t_flat = r2-r1
    return RS * math.log((r2 - RS)/(r1 - RS)) / C


def weak_radial_excess_time(r1, r2):
    # ∫ A dr / c = Rs ln(r2/r1)/c
    return RS * math.log(r2/r1) / C


def run(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # A grid excluding horizon exactly for finite G.
    A_grid = np.linspace(0.0, 0.999, 1000)
    r_grid = RS / np.maximum(A_grid, 1e-12)
    r_grid[0] = np.inf

    df = pd.DataFrame({
        "A": A_grid,
        "F_time_factor_metric": F_of_A(A_grid),
        "G_radial_factor_metric": G_of_A(A_grid),
        "clock_exact_sqrt_1_minus_A": clock_factor(A_grid),
        "clock_weak_1_minus_A_over_2": clock_weak(A_grid),
        "clock_error_exact_minus_weak": clock_factor(A_grid) - clock_weak(A_grid),
        "v_escape_over_c": escape_v_over_c(A_grid),
        "radial_null_dr_dt_over_c": radial_null_speed_over_c(A_grid),
    })
    df.to_csv(outdir / "SF1_A_metric_grid.csv", index=False)

    # Landmark table
    landmarks = [
        {
            "name": "weak field example",
            "r_over_Rs": 1000.0,
            "A": 1/1000.0,
            "meaning": "A small; weak-field expansion should be excellent."
        },
        {
            "name": "ISCO Schwarzschild",
            "r_over_Rs": 3.0,
            "A": 1/3.0,
            "meaning": "Innermost stable circular orbit for massive particles in Schwarzschild."
        },
        {
            "name": "photon sphere Schwarzschild",
            "r_over_Rs": 1.5,
            "A": 2/3.0,
            "meaning": "Circular photon orbit in Schwarzschild."
        },
        {
            "name": "horizon threshold",
            "r_over_Rs": 1.0,
            "A": 1.0,
            "meaning": "A=1; F(A)=0; horizon."
        },
    ]
    lrows = []
    for lm in landmarks:
        A = lm["A"]
        lrows.append({
            **lm,
            "F_1_minus_A": F_of_A(A),
            "clock_sqrt_1_minus_A": clock_factor(A),
            "clock_weak_1_minus_A_over_2": clock_weak(A),
            "v_escape_over_c": escape_v_over_c(A),
            "radial_null_dr_dt_over_c": radial_null_speed_over_c(A),
        })
    landmarks_df = pd.DataFrame(lrows)
    landmarks_df.to_csv(outdir / "SF1_landmarks.csv", index=False)

    # Radial light-time comparison for several intervals
    intervals = [
        (10.0, 100.0),
        (3.0, 100.0),
        (1.5, 100.0),
        (1.1, 100.0),
        (1.01, 100.0),
    ]
    trows = []
    for r1, r2 in intervals:
        tex = exact_radial_excess_time(r1, r2)
        tw = weak_radial_excess_time(r1, r2)
        trows.append({
            "r1_over_Rs": r1,
            "r2_over_Rs": r2,
            "A_at_r1": 1/r1,
            "A_at_r2": 1/r2,
            "exact_excess_time_Rs_over_c_units": tex,
            "weak_integral_A_dr_over_c": tw,
            "exact_minus_weak": tex - tw,
            "weak_over_exact": tw/tex if tex != 0 else np.nan,
        })
    time_df = pd.DataFrame(trows)
    time_df.to_csv(outdir / "SF1_radial_light_time_excess.csv", index=False)

    # Plots
    plt.figure(figsize=(8,5))
    plt.plot(df["A"], df["clock_exact_sqrt_1_minus_A"], label="sqrt(1-A)")
    plt.plot(df["A"], df["clock_weak_1_minus_A_over_2"], "--", label="1 - A/2")
    plt.axvline(1/3, linestyle=":", label="A=1/3 ISCO")
    plt.axvline(2/3, linestyle=":", label="A=2/3 photon sphere")
    plt.axvline(1.0, linestyle=":", label="A=1 horizon")
    plt.title("SF1 clock factor: exact strong-field vs weak-field")
    plt.xlabel("A")
    plt.ylabel("clock factor relative to far observer")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF1_clock_factor.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(df["A"], df["F_time_factor_metric"], label="F(A)=1-A")
    # Clip G for visual.
    plt.plot(df["A"], np.minimum(df["G_radial_factor_metric"], 20), label="G(A)=1/(1-A), clipped at 20")
    plt.axvline(1/3, linestyle=":", label="A=1/3 ISCO")
    plt.axvline(2/3, linestyle=":", label="A=2/3 photon sphere")
    plt.axvline(1.0, linestyle=":", label="A=1 horizon")
    plt.title("SF1 A-metric factors")
    plt.xlabel("A")
    plt.ylabel("metric factor")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF1_metric_factors.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(df["A"], df["v_escape_over_c"], label="v_escape/c = sqrt(A)")
    plt.plot(df["A"], df["radial_null_dr_dt_over_c"], label="radial coordinate light speed/c = 1-A")
    plt.axvline(1, linestyle=":", label="A=1")
    plt.title("SF1 escape and radial-null coordinate behavior")
    plt.xlabel("A")
    plt.ylabel("normalized value")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF1_escape_and_radial_null.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(time_df["A_at_r1"], time_df["exact_excess_time_Rs_over_c_units"], marker="o", label="exact Schwarzschild radial excess")
    plt.plot(time_df["A_at_r1"], time_df["weak_integral_A_dr_over_c"], marker="o", label="weak ∫A dr/c")
    plt.title("SF1 radial light-time excess: exact vs weak A-integral")
    plt.xlabel("A at inner radius r1")
    plt.ylabel("excess time in Rs/c units")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "SF1_radial_light_time_excess.png", dpi=180)
    plt.close()

    # Quantitative weak-field clock errors.
    weak_checks = []
    for A in [1e-9, 1e-6, 1e-3, 1e-2, 0.1, 1/3, 2/3]:
        exact = float(clock_factor(A))
        weak = float(clock_weak(A))
        weak_checks.append({
            "A": A,
            "clock_exact": exact,
            "clock_weak": weak,
            "absolute_error": exact - weak,
            "relative_error": (exact - weak)/exact if exact != 0 else np.nan,
        })
    weak_df = pd.DataFrame(weak_checks)
    weak_df.to_csv(outdir / "SF1_weak_clock_error_checks.csv", index=False)

    verdict = {
        "script": "SF1_A_metric_schwarzschild_recovery_test.py",
        "metric": "ds^2 = -(1-A)c^2dt^2 + dr^2/(1-A) + r^2dOmega^2",
        "A_definition": "A(r)=Rs/r=2GM/(c^2 r)",
        "main_result": [
            "The candidate A-metric is exactly Schwarzschild expressed in A language.",
            "Weak-field clock factor sqrt(1-A) reduces to 1-A/2.",
            "A=1 gives F(A)=0 and radial null dr/dt=0 in Schwarzschild coordinates.",
            "A=2/3 corresponds to the Schwarzschild photon sphere.",
            "A=1/3 corresponds to the Schwarzschild ISCO.",
            "Weak path delay ∫A dr/c approximates radial Schwarzschild excess only in weak field; exact expression diverges near A=1."
        ],
        "interpretation": "SF1 provides a strong-field baseline. Any future STAM metric modification must recover these checks or explain why it departs.",
        "caution": "This is not yet new physics; it is the disciplined anchor: Schwarzschild rewritten as an A-metric.",
    }
    with open(outdir / "SF1_verdict.json", "w") as f:
        json.dump(verdict, f, indent=2)

    print("SF1 A-metric Schwarzschild recovery test complete.")
    print(f"Output directory: {outdir}")
    print()
    print("Landmarks:")
    print(landmarks_df.to_string(index=False))
    print()
    print("Radial light-time excess:")
    print(time_df.to_string(index=False))
    print()
    print("Weak clock checks:")
    print(weak_df.to_string(index=False))
    print()
    print("Main result:")
    print("  The A-metric baseline is Schwarzschild in A-language.")
    print("  A=1 horizon, A=2/3 photon sphere, A=1/3 ISCO are recovered.")
    print("  sqrt(1-A) reduces to 1-A/2 in weak field.")
    print("  ∫A dr/c is the weak path-delay limit; exact radial excess diverges near the horizon.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="results/SF1_A_metric_schwarzschild_recovery")
    args = p.parse_args()
    run(args.outdir)
