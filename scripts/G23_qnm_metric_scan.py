#!/usr/bin/env python3
"""
G23_qnm_metric_scan.py

Scan k(A) = (1-A)^m · (1+A)^n forms for the QNM damping ratio at the photon sphere.
Check whether the (m=2, n=2) "balanced pair = 4 factors" structure that Sean
flagged connects to the framework's other 4-counts (Bekenstein entropy, etc).

Author: Sean Brady / STAM Model-A continuation
Date: 2026-05-11
Exploratory — no commitments, just see what falls out.
"""

from __future__ import annotations
import sys
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

A_PHOTON = 2.0 / 3.0  # photon sphere (h-determined; same in any metric with h = 1-A)


# ============================================================
# Core formulas
# ============================================================

def tau_ratio_at_photon_sphere(m: int, n: int) -> float:
    """τ_alt / τ_GR at the photon sphere for k(A) = (1-A)^m · (1+A)^n.

    Derivation:
        λ² ∝ h·k at photon sphere.
        h·k_GR = (1-A)²; h·k_alt = (1-A)^(m+1) · (1+A)^n.
        Ratio of λ_alt²/λ_GR² = (1-A)^(m-1) · (1+A)^n.
        τ ∝ 1/λ, so τ_alt/τ_GR = √(λ_GR²/λ_alt²) = (1-A)^((1-m)/2) · (1+A)^(-n/2).

    At A = 2/3: (1-A) = 1/3, (1+A) = 5/3.
        τ_ratio = (1/3)^((1-m)/2) · (5/3)^(-n/2)
                = 3^((m-1)/2) · (3/5)^(n/2)
    """
    A = A_PHOTON
    return (1.0 - A) ** ((1.0 - m) / 2.0) * (1.0 + A) ** (-n / 2.0)


def weak_field_recovers(m: int, n: int) -> bool:
    """Does k(A) = (1-A)^m · (1+A)^n recover GR weak-field?

    Radial null speed √(h·k) at small A:
        √((1-A)^(m+1)·(1+A)^n) ≈ 1 + (n - m - 1)·A/2 + O(A²)
    GR is 1 - A. Linear coefficient -1 requires n = m - 1.
    """
    return n == m - 1


# ============================================================
# Scan
# ============================================================

def scan_grid(m_range=range(1, 6), n_range=range(0, 6)):
    rows = []
    for m in m_range:
        for n in n_range:
            ratio = tau_ratio_at_photon_sphere(m, n)
            wf = weak_field_recovers(m, n)
            mod_count = (m - 1) + n  # modification factors beyond GR's h baseline
            mod_form = f"(1-A)^{m-1}·(1+A)^{n}" if (m - 1 > 0 or n > 0) else "1 (no mod)"
            rows.append({
                "m": m,
                "n": n,
                "k(A)": f"(1-A)^{m}·(1+A)^{n}",
                "modification beyond h": mod_form,
                "mod factor count": mod_count,
                "k at A=2/3": (1.0 / 3.0) ** m * (5.0 / 3.0) ** n,
                "tau_ratio": ratio,
                "WF OK": wf,
            })
    return pd.DataFrame(rows)


# ============================================================
# Plots
# ============================================================

def plot_heatmap(df, out_path: Path):
    pivot = df.pivot(index="n", columns="m", values="tau_ratio")

    fig, ax = plt.subplots(figsize=(11, 7))
    im = ax.imshow(pivot.values, cmap="RdYlBu_r", aspect="auto", origin="lower",
                   extent=(df["m"].min() - 0.5, df["m"].max() + 0.5,
                           df["n"].min() - 0.5, df["n"].max() + 0.5),
                   vmin=0.3, vmax=3.5)

    for _, row in df.iterrows():
        m, n = row["m"], row["n"]
        v = row["tau_ratio"]
        wf = row["WF OK"]
        is_ma = (m == 3 and n == 2)
        is_22 = (m == 2 and n == 2)
        color = "black"
        marker = ""
        if is_ma:
            marker = " ★"  # Model-A
        elif is_22:
            marker = " □"  # (2,2) "4" structure
        elif wf:
            marker = " ●"
        ax.text(m, n, f"{v:.3f}{marker}",
                ha="center", va="center", color=color, fontsize=9,
                weight=("bold" if (wf or is_22) else "normal"))

    diag_x = [m for m in df["m"].unique()]
    diag_y = [m - 1 for m in df["m"].unique()]
    ax.plot(diag_x, diag_y, "k--", alpha=0.4, linewidth=2,
            label="n = m−1 (weak-field-respecting)")

    ax.plot(2, 2, "kx", markersize=18, markeredgewidth=3,
            label="(m=2, n=2) — '4' structure (fails WF)")
    ax.plot(3, 2, "k*", markersize=20, markeredgewidth=2,
            label="Model-A (m=3, n=2)")

    ax.set_xlabel("m (power of (1−A) in k(A))")
    ax.set_ylabel("n (power of (1+A) in k(A))")
    ax.set_title("QNM damping ratio τ_alt/τ_GR at A=2/3 photon sphere")
    plt.colorbar(im, ax=ax, label="τ_alt / τ_GR")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_p_scaling(out_path: Path):
    """Damping ratio along the weak-field-respecting diagonal: k = h · (1-A²)^p."""
    p_range = np.array([0, 1, 2, 3, 4, 5, 6])
    ratios = np.array([tau_ratio_at_photon_sphere(p + 1, p) for p in p_range])

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(p_range, ratios, "o-", color="tab:blue", markersize=10, linewidth=2)
    for p, r in zip(p_range, ratios):
        ax.annotate(f"τ_ratio={r:.3f}",
                    xy=(p, r), xytext=(p, r + 0.15),
                    ha="center", fontsize=9)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5, label="GR (no enhancement)")
    ax.axvline(2, color="red", linestyle=":", alpha=0.5, label="p=2 → Model-A")
    ax.set_xlabel("p (number of (1−A²) factors in modification)")
    ax.set_ylabel("τ_alt / τ_GR at photon sphere")
    ax.set_title("Damping ratio along weak-field-respecting diagonal\n"
                 "k(A) = h · (1−A²)^p, m = p+1, n = p")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


# ============================================================
# Markdown
# ============================================================

def write_markdown(df, out_path: Path):
    md = []
    md.append("# G23: QNM damping ratio scan over metric forms\n\n")
    md.append("**Date:** 2026-05-11\n\n")
    md.append("Exploratory scan. Goal: see whether the (m=2, n=2) connection to "
              "the framework's structural '4' (per Sean) shows anything clean in "
              "the QNM landscape.\n\n")

    md.append("## What's scanned\n\n")
    md.append("k(A) = (1−A)^m · (1+A)^n for integer m ∈ [1,5], n ∈ [0,5]. ")
    md.append("h(A) = 1−A fixed. Photon sphere at A = 2/3.\n\n")

    md.append("## Damping ratio formula (photon sphere, A=2/3)\n\n")
    md.append("```\nτ_alt / τ_GR = 3^((m−1)/2) · (3/5)^(n/2)\n```\n\n")

    md.append("## Weak-field constraint\n\n")
    md.append("k(A) recovers GR weak-field only when **n = m − 1**. ")
    md.append("This is the 'diagonal' in the (m, n) grid.\n\n")

    wf_df = df[df["WF OK"]].copy()
    md.append("### Weak-field-respecting family\n\n")
    md.append("Along n = m−1, write k = h · (1−A²)^p with p = n, m = p+1.\n\n")
    md.append(wf_df[["m", "n", "k(A)", "tau_ratio"]].to_markdown(index=False, floatfmt=".4f"))
    md.append("\n\nScaling: **τ_ratio = (9/5)^(p/2)**. Each (1−A²) factor multiplies by √(9/5) ≈ 1.342.\n\n")

    md.append("## The (m=2, n=2) connection to '4'\n\n")
    r22 = tau_ratio_at_photon_sphere(2, 2)
    md.append(f"At (m=2, n=2): k = (1−A)²(1+A)² = (1−A²)². τ_ratio = {r22:.4f}.\n\n")
    md.append("**FAILS GR weak-field recovery** (n=2, but m−1=1). Light speed at small A "
              "goes as √(1−A) instead of (1−A) — factor-of-2 error, detectable by Cassini-class tests.\n\n")
    md.append("**However** — Model-A's actual k = h · (1−A²)². The MODIFICATION beyond h is (1−A²)², "
              "which has m_mod=2, n_mod=2 in its own exponents. **Total modification factors: 4.**\n\n")
    md.append("This is the structural '4' Sean noted. The (m=2, n=2) reading applies to the "
              "modification factor, not the full k(A).\n\n")
    md.append("**Reading**: Model-A's modification factor (1−A²)² has exactly the same 2×2 "
              "structure as Bekenstein-Hawking's α = 2-face × 2-gravity-bridge = 4 area-per-entry. "
              "Both are '2 independent pair axes meeting at the boundary.'\n\n")

    md.append("## Damping ratio along the weak-field-respecting family\n\n")
    md.append("| p (pair²-factors in mod) | (m, n) | k(A) | τ_ratio | Comment |\n")
    md.append("|---|---|---|---|---|\n")
    md.append("| 0 | (1, 0) | (1−A) | 1.0000 | GR (no modification) |\n")
    md.append("| 1 | (2, 1) | (1−A)²(1+A) | 1.3416 | minimal pair |\n")
    md.append("| 2 | (3, 2) | (1−A)³(1+A)² = (1−A)(1−A²)² | 1.8000 | **Model-A; mod has 4 factors** |\n")
    md.append("| 3 | (4, 3) | (1−A)⁴(1+A)³ | 2.4150 | |\n")
    md.append("| 4 | (5, 4) | (1−A)⁵(1+A)⁴ | 3.2404 | |\n")
    md.append("| 5 | (6, 5) | (1−A)⁶(1+A)⁵ | 4.3478 | |\n\n")

    md.append("**Model-A's p=2 is the first case along the diagonal with a fully "
              "'doubled-pair' modification (4 factors).** p=1 has 2 modification factors "
              "(one (1−A) × one (1+A)) — a single pair, not the structural 4. p=2 is the "
              "first 'two-independent-pair-axes' case.\n\n")

    md.append("## Full scan grid\n\n")
    md.append(df.to_markdown(index=False, floatfmt=".4f"))
    md.append("\n\n")

    md.append("## What this empirically surfaces\n\n")
    md.append("1. **The weak-field constraint n = m−1 is a sharp line in the grid.** "
              "Off-diagonal cases fail observational tests (e.g., (m=2, n=2) gives a factor-2 "
              "error in radial light speed at small A).\n\n")
    md.append("2. **Along the WF-respecting diagonal, the damping ratio grows monotonically with p.** "
              "Each additional (1−A²) factor in the modification multiplies the ratio by ~1.342.\n\n")
    md.append("3. **p=2 (Model-A) is the smallest WF-respecting case with 4 modification factors.** "
              "The next case (p=3) has 6 factors; p=1 has 2. p=2 is uniquely the '4-factor' case.\n\n")
    md.append("4. **If pair structure forces '4 modification factors at the boundary'** "
              "(matching Bekenstein's 2-face × 2-gravity-bridge structure), then p=2 is structurally "
              "forced, and Model-A's specific k(A) = (1−A)(1−A²)² is the natural minimal form. "
              "The 1.80 damping ratio at the photon sphere becomes a structural prediction, not "
              "a free choice.\n\n")
    md.append("5. **LIGO tension stays.** The 1.80 ratio is significantly above observed "
              "(GR-consistent) ringdown damping times. If the '4-factor' structural reading is "
              "right, the tension is a real concern for the framework. Either (a) the eikonal "
              "overshoots and exact Regge-Wheeler gives a smaller ratio, or (b) the structural "
              "4-reading is wrong, or (c) Model-A's k(A) needs revision.\n\n")

    out_path.write_text("".join(md), encoding="utf-8")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 78)
    print("G23: QNM damping ratio scan over k(A) = (1-A)^m · (1+A)^n")
    print("=" * 78)
    print()
    print(f"Photon sphere at A = 2/3 (unchanged by k modification; depends only on h)")
    print(f"Pair product (1−A)(1+A) at A=2/3 = (1/3)(5/3) = 5/9")
    print(f"Inverse pair product = 9/5 = 1.80")
    print()

    df = scan_grid()

    print("Full scan:")
    print(df.to_string(index=False, float_format='%.4f'))
    print()

    # Highlight specific points
    ma = df[(df["m"] == 3) & (df["n"] == 2)].iloc[0]
    print(f"Model-A (m=3, n=2): τ_ratio = {ma['tau_ratio']:.4f} = 9/5 = 1.80  [WF OK: {ma['WF OK']}]")
    print(f"  Modification = (1−A²)² has 4 pair factors (m_mod=2, n_mod=2)")
    print()

    bal = df[(df["m"] == 2) & (df["n"] == 2)].iloc[0]
    print(f"Balanced (m=2, n=2): τ_ratio = {bal['tau_ratio']:.4f}  [WF OK: {bal['WF OK']}]")
    print(f"  Fails WF; but 4-factor structure matches Sean's '4' intuition for the modification.")
    print()

    print("Weak-field-respecting diagonal (n = m−1):")
    for p in range(0, 6):
        m = p + 1
        n = p
        r = tau_ratio_at_photon_sphere(m, n)
        mod_count = (m - 1) + n
        print(f"  p = {p}: (m={m}, n={n}), mod factors = {mod_count}, τ_ratio = {r:.4f}")
    print()

    print("Scaling along diagonal: τ_ratio = (9/5)^(p/2)")
    print(f"  At p=2 (Model-A): (9/5)^1 = 1.80")
    print(f"  At p=4: (9/5)² = {(9/5)**2:.4f}")
    print()

    print("Structural reading:")
    print("  - p=1: modification has 2 factors (single pair)")
    print("  - p=2: modification has 4 factors (two-pair-squared) ← Model-A")
    print("  - p=3: modification has 6 factors (three-pair-squared)")
    print("  Model-A is the first case where the modification is 'doubled-pair' (4 factors).")
    print()

    # Save
    df.to_csv(RESULTS / "G23_qnm_metric_scan.csv", index=False)
    plot_heatmap(df, PLOTS / "G23_qnm_metric_scan_heatmap.png")
    plot_p_scaling(PLOTS / "G23_p_scaling.png")
    write_markdown(df, RESULTS / "G23_qnm_metric_scan_summary.md")
    print("Files written:")
    print(f"  {RESULTS / 'G23_qnm_metric_scan.csv'}")
    print(f"  {PLOTS / 'G23_qnm_metric_scan_heatmap.png'}")
    print(f"  {PLOTS / 'G23_p_scaling.png'}")
    print(f"  {RESULTS / 'G23_qnm_metric_scan_summary.md'}")


if __name__ == "__main__":
    main()
