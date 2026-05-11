#!/usr/bin/env python3
"""
G22_proper_time_radial_infall.py

Proper-time of radial infall in Model-A vs Schwarzschild GR.
Author: Sean Brady / STAM Model-A
Date: 2026-05-11

Test particle released from rest at the cosmic baseline A = A_0 = 1/(12π).
Free-fall radially inward. Compute proper time τ the particle's own clock
reads as it falls to various A values.

Compare Model-A's metric k(A) = (1-A)(1-A²)² to Schwarzschild's k(A) = 1-A.
Both use h(A) = 1-A in g_tt.

Geodesic for radial timelike free-fall, released from rest at A_release:
    (dr/dτ)² = (k/h)·(A - A_release)

In natural units (c = 1, Rs = 1), A = 1/r, dr = -dA/A². Proper time:
    τ(A_end) = ∫_{A_release}^{A_end} dA / [A²·√((k/h)·(A - A_release))]

GR: k/h = 1. Model-A: k/h = (1-A²)².

This is an exploratory test. We don't predict outcomes; we compute and look.
"""

from __future__ import annotations
import sys
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

A_0 = 1.0 / (12.0 * math.pi)


# ============================================================
# Integrands
# ============================================================

def integrand_GR(A, A_release):
    if A <= A_release:
        return 0.0
    return 1.0 / (A**2 * np.sqrt(A - A_release))


def integrand_MA(A, A_release):
    if A <= A_release or A >= 1.0:
        return 0.0
    return 1.0 / (A**2 * (1.0 - A**2) * np.sqrt(A - A_release))


def proper_time_GR(A_end, A_release):
    if A_end <= A_release:
        return 0.0
    val, _ = quad(integrand_GR, A_release, A_end,
                  args=(A_release,), points=[A_release],
                  epsabs=1e-12, epsrel=1e-10, limit=400)
    return val


def proper_time_MA(A_end, A_release):
    if A_end <= A_release:
        return 0.0
    if A_end >= 1.0:
        return float("inf")
    val, _ = quad(integrand_MA, A_release, A_end,
                  args=(A_release,), points=[A_release],
                  epsabs=1e-12, epsrel=1e-10, limit=400)
    return val


# ============================================================
# Tabulation
# ============================================================

def make_table(A_release):
    landmarks = [
        ("release (A_0)",          A_release),
        ("A = 0.05",               0.05),
        ("A = 0.10",               0.10),
        ("A = 0.20",               0.20),
        ("A = 1/3 (ISCO)",         1.0/3.0),
        ("A = 0.50",               0.50),
        ("A = 2/3 (photon sph)",   2.0/3.0),
        ("A = 0.80",               0.80),
        ("A = 0.90",               0.90),
        ("A = 0.99",               0.99),
        ("A = 0.999",              0.999),
        ("A = 0.9999",             0.9999),
        ("A = 0.99999",            0.99999),
    ]
    rows = []
    for name, A_end in landmarks:
        tau_GR = proper_time_GR(A_end, A_release)
        tau_MA = proper_time_MA(A_end, A_release)
        ratio = tau_MA / tau_GR if tau_GR > 0 else float("nan")
        diff = tau_MA - tau_GR
        rows.append({
            "landmark":          name,
            "A":                 A_end,
            "r/Rs":              1.0/A_end if A_end > 0 else float("inf"),
            "tau_GR (Rs/c)":     tau_GR,
            "tau_MA (Rs/c)":     tau_MA,
            "ratio MA/GR":       ratio,
            "diff MA-GR (Rs/c)": diff,
        })
    return pd.DataFrame(rows)


# ============================================================
# Divergence coefficient (Model-A only)
# ============================================================

def divergence_coefficient(A_release):
    """Fit τ_MA(A) near A=1 to form: τ = -slope · ln(1-A) + intercept.

    Analytical leading coefficient: 1 / [2·√(1 - A_release)·1²] = 1/(2√(1-A_0))
    because near A=1: integrand ≈ 1/(2·(1-A)·√(1-A_release)).
    """
    A_test = np.array([0.99, 0.999, 0.9999, 0.99999])
    tau_test = np.array([proper_time_MA(A, A_release) for A in A_test])
    x = -np.log(1.0 - A_test)
    slope, intercept = np.polyfit(x, tau_test, 1)
    analytical = 1.0 / (2.0 * math.sqrt(1.0 - A_release))
    return {
        "fitted_slope":  slope,
        "analytical":    analytical,
        "ratio":         slope / analytical,
        "intercept":     intercept,
    }


# ============================================================
# Plots
# ============================================================

def plot_proper_time(A_release):
    A_grid = np.linspace(A_release + 1e-5, 0.999, 800)
    tau_GR_list = np.array([proper_time_GR(A, A_release) for A in A_grid])
    tau_MA_list = np.array([proper_time_MA(A, A_release) for A in A_grid])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.plot(A_grid, tau_GR_list, color="tab:gray", linewidth=2, label="GR Schwarzschild")
    ax.plot(A_grid, tau_MA_list, color="tab:blue", linewidth=2, label="Model-A")
    ax.axvline(1.0/3.0, color="tab:green", linestyle=":", alpha=0.5, label="A=1/3 (ISCO)")
    ax.axvline(2.0/3.0, color="tab:orange", linestyle=":", alpha=0.5, label="A=2/3 (photon)")
    ax.set_xlabel("A (along infall)")
    ax.set_ylabel("Proper time τ (Rs/c)")
    ax.set_title("Proper time of radial infall from rest at A_0")
    ax.legend()
    ax.grid(True, linewidth=0.3)
    ax.set_ylim(0, max(tau_MA_list[A_grid < 0.95].max() * 1.1, 200))

    ax = axes[1]
    A_grid_zoom = np.linspace(0.5, 0.99999, 500)
    tau_GR_zoom = np.array([proper_time_GR(A, A_release) for A in A_grid_zoom])
    tau_MA_zoom = np.array([proper_time_MA(A, A_release) for A in A_grid_zoom])
    ax.plot(A_grid_zoom, tau_GR_zoom, color="tab:gray", linewidth=2, label="GR (finite to horizon)")
    ax.plot(A_grid_zoom, tau_MA_zoom, color="tab:blue", linewidth=2, label="Model-A (diverges)")
    ax.axvline(2.0/3.0, color="tab:orange", linestyle=":", alpha=0.5, label="A=2/3 (photon)")
    ax.set_xlabel("A (zoom near A=1)")
    ax.set_ylabel("Proper time τ (Rs/c)")
    ax.set_title("Near saturation A=1: Model-A diverges, GR finite")
    ax.legend()
    ax.grid(True, linewidth=0.3)

    out = PLOTS / "G22_proper_time_infall.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_ratio(A_release):
    A_grid = np.linspace(A_release + 1e-4, 0.99, 800)
    ratios = []
    diffs = []
    for A in A_grid:
        tau_GR = proper_time_GR(A, A_release)
        tau_MA = proper_time_MA(A, A_release)
        if tau_GR > 0:
            ratios.append(tau_MA / tau_GR)
            diffs.append(tau_MA - tau_GR)
        else:
            ratios.append(float("nan"))
            diffs.append(float("nan"))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.plot(A_grid, ratios, color="tab:purple", linewidth=2)
    ax.axhline(1.0, color="black", linewidth=0.5)
    ax.axvline(1.0/3.0, color="tab:green", linestyle=":", alpha=0.5, label="A=1/3 (ISCO)")
    ax.axvline(2.0/3.0, color="tab:orange", linestyle=":", alpha=0.5, label="A=2/3 (photon)")
    ax.set_xlabel("A_end")
    ax.set_ylabel("τ_MA / τ_GR")
    ax.set_title("Ratio: Model-A proper time / GR proper time")
    ax.legend()
    ax.grid(True, linewidth=0.3)

    ax = axes[1]
    ax.plot(A_grid, diffs, color="tab:red", linewidth=2)
    ax.axhline(0.0, color="black", linewidth=0.5)
    ax.axvline(1.0/3.0, color="tab:green", linestyle=":", alpha=0.5)
    ax.axvline(2.0/3.0, color="tab:orange", linestyle=":", alpha=0.5)
    ax.set_xlabel("A_end")
    ax.set_ylabel("τ_MA − τ_GR (Rs/c)")
    ax.set_title("Excess Model-A proper time over GR")
    ax.grid(True, linewidth=0.3)

    out = PLOTS / "G22_ratio_and_diff.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# ============================================================
# Markdown summary
# ============================================================

def write_markdown(df, A_release, div, plot_paths):
    md = []
    md.append("# G22: Proper-time of radial infall in Model-A vs Schwarzschild\n\n")
    md.append("**Date:** 2026-05-11\n\n")
    md.append("**Setup.** Test particle released from rest at A = A_0 = 1/(12π) "
              f"≈ {A_0:.6f}. Free-fall radial infall. Compute the proper time τ "
              "the particle's own clock reads as it reaches various A values.\n\n")
    md.append("**Metric inputs (natural units c=1, Rs=1):**\n\n")
    md.append("- GR Schwarzschild: g_tt = -(1-A), g_rr = 1/(1-A); k/h = 1\n")
    md.append("- Model-A: g_tt = -(1-A), g_rr = 1/[(1-A)(1-A²)²]; k/h = (1-A²)²\n\n")
    md.append("**Geodesic (radial timelike free-fall):**\n\n")
    md.append("```\n(dr/dτ)² = (k/h)·(A − A_release)\nτ(A_end) = ∫_{A_release}^{A_end} dA / [A²·√((k/h)·(A − A_release))]\n```\n\n")
    md.append("## Numerical results (proper time in units of Rs/c)\n\n")
    md.append(df.to_markdown(index=False, floatfmt=".6f"))
    md.append("\n\n")
    md.append("## Divergence coefficient (Model-A near A=1)\n\n")
    md.append("Fit form: τ_MA(A) ≈ −slope · ln(1−A) + intercept (near A=1)\n\n")
    md.append(f"- Fitted slope:                       {div['fitted_slope']:.6f}\n")
    md.append(f"- Analytical 1/(2·√(1−A_0)):          {div['analytical']:.6f}\n")
    md.append(f"- Ratio (fit/analytical):              {div['ratio']:.6f}\n")
    md.append(f"- Fitted intercept:                    {div['intercept']:.6f}\n\n")
    md.append("## Files\n\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`\n")
    md.append("- `results/G22_proper_time_table.csv`\n")
    out = RESULTS / "G22_proper_time_radial_infall_summary.md"
    out.write_text("".join(md), encoding="utf-8")
    return out


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 78)
    print("G22: Proper-time of radial infall in Model-A vs GR Schwarzschild")
    print("=" * 78)
    print(f"A_0 = 1/(12π) = {A_0:.8f}")
    print(f"Release point: A = A_0 (rest at the cosmic baseline)")
    print(f"Units: Rs/c throughout.")
    print()

    A_release = A_0

    # Table
    df = make_table(A_release)
    df.to_csv(RESULTS / "G22_proper_time_table.csv", index=False)
    print("Proper-time table:")
    print(df.to_string(index=False, float_format='%.6f'))
    print()

    # Interval times
    A_isco = 1.0/3.0
    A_ph = 2.0/3.0

    print("Cumulative proper time from release:")
    print(f"  GR  → ISCO (A=1/3):      {proper_time_GR(A_isco, A_release):.6f}")
    print(f"  MA  → ISCO (A=1/3):      {proper_time_MA(A_isco, A_release):.6f}")
    print(f"  GR  → photon (A=2/3):    {proper_time_GR(A_ph, A_release):.6f}")
    print(f"  MA  → photon (A=2/3):    {proper_time_MA(A_ph, A_release):.6f}")
    print(f"  GR  → A=0.99:            {proper_time_GR(0.99, A_release):.6f}")
    print(f"  MA  → A=0.99:            {proper_time_MA(0.99, A_release):.6f}")
    print(f"  GR  → A=1 (horizon):     {proper_time_GR(0.9999999, A_release):.6f}")
    print(f"  MA  → A=1 (saturation):  diverges")
    print()

    print("Interval proper times (GR):")
    print(f"  release → ISCO:           {proper_time_GR(A_isco, A_release):.6f}")
    print(f"  ISCO → photon:            {proper_time_GR(A_ph, A_release) - proper_time_GR(A_isco, A_release):.6f}")
    print(f"  photon → horizon:         {proper_time_GR(0.9999999, A_release) - proper_time_GR(A_ph, A_release):.6f}")
    print()
    print("Interval proper times (Model-A):")
    print(f"  release → ISCO:           {proper_time_MA(A_isco, A_release):.6f}")
    print(f"  ISCO → photon:            {proper_time_MA(A_ph, A_release) - proper_time_MA(A_isco, A_release):.6f}")
    print(f"  photon → A=0.9:           {proper_time_MA(0.9, A_release) - proper_time_MA(A_ph, A_release):.6f}")
    print(f"  A=0.9 → A=0.99:           {proper_time_MA(0.99, A_release) - proper_time_MA(0.9, A_release):.6f}")
    print(f"  A=0.99 → A=0.999:         {proper_time_MA(0.999, A_release) - proper_time_MA(0.99, A_release):.6f}")
    print()

    # Divergence analysis
    div = divergence_coefficient(A_release)
    print("Model-A divergence near A=1:")
    print(f"  τ_MA(A) ≈ -slope · ln(1−A) + intercept")
    print(f"  Fitted slope:                  {div['fitted_slope']:.6f}")
    print(f"  Analytical 1/(2·√(1-A_0)):     {div['analytical']:.6f}")
    print(f"  Ratio fit/analytical:          {div['ratio']:.6f}")
    print(f"  Fitted intercept:              {div['intercept']:.6f}")
    print()

    # Plots
    p1 = plot_proper_time(A_release)
    p2 = plot_ratio(A_release)

    summary = write_markdown(df, A_release, div, [p1, p2])
    print(f"Summary: {summary}")
    print(f"Plot: {p1}")
    print(f"Plot: {p2}")


if __name__ == "__main__":
    main()
