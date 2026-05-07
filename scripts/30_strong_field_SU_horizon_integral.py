#!/usr/bin/env python3
"""
30_strong_field_SU_horizon_integral.py

STAM Model-A strong-field SU horizon integral test.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    Test the proposed strong-field traversal scaling:

        S_h(A) = 1 / (1 - A)
        A(r) = Rs / r
        dSU_A = dr / (1 - Rs/r)

    For radial intervals outside the horizon, compare ordinary radial distance Δr
    against the strong-field SU_A interval:

        SU_A(r_inner -> r_outer) = (r_outer - r_inner)
            + Rs * ln((r_outer - Rs) / (r_inner - Rs))

    As r_inner approaches Rs from above, Δr remains finite/small while SU_A
    diverges logarithmically.

Notes:
    - The script uses dimensionless radius x = r/Rs, so Rs = 1.
    - All SU_A results are therefore in units of Rs.
    - The direction is treated as outward from r_inner to r_outer, using positive magnitudes.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

RS = 1.0  # dimensionless units, r/Rs


def A_of_x(x: float) -> float:
    """Dimensionless accumulation field A = Rs/r = 1/x."""
    if x <= 0:
        raise ValueError("x must be positive")
    return RS / x


def S_h_of_A(A: float) -> float:
    """Strong-field scale factor S_h(A)=1/(1-A)."""
    if A >= 1:
        return math.inf
    return 1.0 / (1.0 - A)


def SU_interval_x(x_inner: float, x_outer: float) -> float:
    """
    Strong-field SU_A radial interval in Rs units.

    Valid only for x_outer > x_inner > 1.
    Formula with Rs=1:
        SU = (x_outer - x_inner) + ln((x_outer - 1)/(x_inner - 1))
    """
    if not (x_outer > x_inner > 1.0):
        raise ValueError("Require x_outer > x_inner > 1")
    return (x_outer - x_inner) + math.log((x_outer - 1.0) / (x_inner - 1.0))


def ordinary_interval_x(x_inner: float, x_outer: float) -> float:
    """Ordinary radial interval Δr in Rs units."""
    if not (x_outer > x_inner > 1.0):
        raise ValueError("Require x_outer > x_inner > 1")
    return x_outer - x_inner


def build_scale_table(A_values: Iterable[float]) -> pd.DataFrame:
    rows = []
    for A in A_values:
        rows.append({
            "A": A,
            "S_h(A)=1/(1-A)": S_h_of_A(A),
        })
    return pd.DataFrame(rows)


def build_interval_table(intervals: Iterable[tuple[float, float]]) -> pd.DataFrame:
    rows = []
    for x_outer, x_inner in intervals:
        delta_r = ordinary_interval_x(x_inner, x_outer)
        su = SU_interval_x(x_inner, x_outer)
        rows.append({
            "outer_r_over_Rs": x_outer,
            "inner_r_over_Rs": x_inner,
            "A_outer": A_of_x(x_outer),
            "A_inner": A_of_x(x_inner),
            "ordinary_delta_r_Rs": delta_r,
            "SU_A_Rs": su,
            "SU_A_over_delta_r": su / delta_r,
        })
    return pd.DataFrame(rows)


def save_plot_sh_vs_A() -> Path:
    A_vals = [1.0 - 10 ** (-p / 20) for p in range(0, 181)]  # dense near 1
    A_vals = [a for a in A_vals if 0 <= a < 1]
    S_vals = [S_h_of_A(a) for a in A_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(A_vals, S_vals)
    plt.yscale("log")
    plt.xlabel("A")
    plt.ylabel("S_h(A) = 1 / (1 - A)")
    plt.title("Strong-field scale factor diverges as A approaches 1")
    plt.grid(True, which="both", linewidth=0.4)
    out = PLOTS / "30_sh_of_A_near_horizon.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def save_plot_su_ratio(interval_df: pd.DataFrame) -> Path:
    df = interval_df.sort_values("inner_r_over_Rs", ascending=False)
    plt.figure(figsize=(8, 5))
    plt.plot(df["inner_r_over_Rs"], df["SU_A_over_delta_r"], marker="o")
    plt.yscale("log")
    plt.gca().invert_xaxis()
    plt.xlabel("Inner radius r/Rs approaching 1 from above")
    plt.ylabel("SU_A / ordinary Δr")
    plt.title("SU_A grows relative to ordinary distance near the horizon")
    plt.grid(True, which="both", linewidth=0.4)
    out = PLOTS / "30_su_over_delta_r_vs_inner_radius.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def save_plot_su_vs_delta(interval_df: pd.DataFrame) -> Path:
    df = interval_df.sort_values("inner_r_over_Rs", ascending=False)
    labels = [f"{row.outer_r_over_Rs:g}→{row.inner_r_over_Rs:g}" for row in df.itertuples()]
    x = range(len(df))

    plt.figure(figsize=(10, 5))
    plt.plot(x, df["ordinary_delta_r_Rs"], marker="o", label="ordinary Δr")
    plt.plot(x, df["SU_A_Rs"], marker="o", label="strong-field SU_A")
    plt.yscale("log")
    plt.xticks(list(x), labels, rotation=35, ha="right")
    plt.ylabel("Distance measure in Rs units")
    plt.title("Ordinary distance remains finite while SU_A grows")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.4)
    out = PLOTS / "30_ordinary_distance_vs_strong_field_su.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def save_plot_log_divergence() -> Path:
    eps_vals = [10 ** (-p / 10) for p in range(1, 91)]  # epsilon = r/Rs - 1; start below outer radius
    x_outer = 2.0
    su_vals = [SU_interval_x(1.0 + eps, x_outer) for eps in eps_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(eps_vals, su_vals, marker=".")
    plt.xscale("log")
    plt.gca().invert_xaxis()
    plt.xlabel("epsilon = r_inner/Rs - 1")
    plt.ylabel("SU_A from r_inner to 2Rs")
    plt.title("Logarithmic SU_A divergence as epsilon approaches 0")
    plt.grid(True, which="both", linewidth=0.4)
    out = PLOTS / "30_log_divergence_su_to_2Rs.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_markdown(scale_df: pd.DataFrame, interval_df: pd.DataFrame, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A Strong-Field SU Horizon Integral Test\n")
    md.append("## Purpose\n")
    md.append(
        "This test evaluates the proposed strong-field traversal layer `S_h(A)=1/(1-A)` "
        "with `A(r)=Rs/r`. The goal is to compare ordinary radial distance against the "
        "strong-field SU measure as the inner radius approaches the horizon from outside.\n"
    )
    md.append("## Core equations\n")
    md.append("```text\nA(r) = Rs/r\nS_h(A) = 1/(1-A)\ndSU_A = dr/(1-Rs/r)\nSU_A = ∫ r/(r-Rs) dr\nSU_A = r + Rs ln|r-Rs| + constant\n```\n")
    md.append("For an outward interval outside the horizon, using dimensionless units `x=r/Rs`:\n")
    md.append("```text\nSU_A(x_inner → x_outer) = (x_outer - x_inner) + ln[(x_outer - 1)/(x_inner - 1)]\n```\n")

    md.append("## Scale factor table\n")
    md.append(scale_df.to_markdown(index=False, floatfmt=".12g"))
    md.append("\n\n## Interval comparison table\n")
    md.append(interval_df.to_markdown(index=False, floatfmt=".12g"))
    md.append("\n\n## Interpretation\n")
    md.append(
        "The ordinary radial interval can become very small near the horizon, while `SU_A` remains large "
        "and increases without bound as `r_inner/Rs → 1+`. This matches the intended strong-field behavior: "
        "the A=1 horizon is already derived from `A=(v_escape/c)^2`, and the candidate `S_h(A)` layer turns "
        "that horizon threshold into a divergent outside-comparison/traversal measure.\n"
    )
    md.append(
        "This does not by itself complete the strong-field STAM claim. It does establish that the proposed "
        "local SU layer produces the expected logarithmic horizon divergence while preserving finite ordinary "
        "radial distance outside the horizon.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT)}`")
    md.append("")

    out = RESULTS / "30_strong_field_SU_horizon_integral_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    A_values = [0, 0.5, 0.9, 0.99, 0.999, 0.999999, 0.999999999]
    intervals = [
        (10.0, 2.0),
        (2.0, 1.1),
        (1.1, 1.01),
        (1.01, 1.001),
        (1.001, 1.000001),
        (1.000001, 1.000000001),
    ]

    scale_df = build_scale_table(A_values)
    interval_df = build_interval_table(intervals)

    scale_csv = RESULTS / "30_scale_factor_table.csv"
    interval_csv = RESULTS / "30_interval_comparison_table.csv"
    scale_df.to_csv(scale_csv, index=False)
    interval_df.to_csv(interval_csv, index=False)

    plot_paths = [
        save_plot_sh_vs_A(),
        save_plot_su_ratio(interval_df),
        save_plot_su_vs_delta(interval_df),
        save_plot_log_divergence(),
    ]
    summary_md = write_markdown(scale_df, interval_df, plot_paths)

    print("STAM Model-A Strong-Field SU Horizon Integral Test")
    print("=" * 62)
    print("\nScale factor table:")
    print(scale_df.to_string(index=False))
    print("\nInterval comparison table:")
    print(interval_df.to_string(index=False))
    print("\nFiles written:")
    for p in [scale_csv, interval_csv, summary_md, *plot_paths]:
        print(f"- {p}")


if __name__ == "__main__":
    main()
