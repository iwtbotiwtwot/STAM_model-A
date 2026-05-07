#!/usr/bin/env python3
"""
Q11_bekenstein_hawking_entropy.py

STAM Model-A: Derivation of Bekenstein-Hawking entropy S = A / (4 ell_P^2)
from the bold-STAM bubble picture and the Q10 phase-boundary temperature.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    Derive the Bekenstein-Hawking black hole entropy formula
        S = k_B A_horizon / (4 ell_P^2)        (SI units)
        S = A / 4                              (natural units)
    from STAM-native ingredients:
        - The bubble picture: all in-falling matter accumulates on the 2D
          A=1 phase boundary; there is no interior, so all degrees of freedom
          are on the surface (entropy must scale with area, not volume).
        - The Q10 phase-boundary temperature T = hbar c |grad A| / (4 pi k_B).
        - The first law of thermodynamics dE = T dS applied to a black hole
          treated as a thermal phase-boundary system with E = M c^2.

Two derivations, same answer:

    Path A — physical (bubble picture):
        - The bubble picture forces S to scale with area, not volume,
          because matter / information is on the boundary surface.
        - The proportionality constant 1/4 follows from saturation of the
          holographic bound: a 2D phase boundary with no interior must
          saturate the bound by construction.

    Path B — quantitative (first law + Q10 temperature):
        - Apply dE = T dS with E = M c^2, T = T_Q10.
        - Integrate dS / dM, using Rs = 2 G M / c^2.
        - The result is S = k_B c^3 A / (4 hbar G), exactly the
          Bekenstein-Hawking formula.

The script does Path B numerically (showing the integral works out exactly)
and walks through Path A in the markdown summary.

What this gives bold STAM:
    - A complete and self-consistent BH thermodynamics chain:
      Q8 / Q10 (temperature) -> Q11 (entropy) -> Q12 (first law verification).
    - The 1/4 prefactor in S = A / 4 comes from the same gravity-bridge
      structure that gives 1/(4 pi) in the temperature.
    - Reproduces the Bekenstein-Hawking formula exactly, while the underlying
      physical mechanism (area entropy by construction, no interior) is
      different from the standard derivation (entropy from QFT mode counting
      across the horizon).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = 2.99792458e8
HBAR = 1.054571817e-34
KB = 1.380649e-23
G = 6.67430e-11
M_SUN = 1.98847e30

# Planck length squared, ell_P^2 = hbar G / c^3
ELL_P_SQ = HBAR * G / (C ** 3)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Black hole geometric quantities ---

def Rs_of_M(M: float) -> float:
    return 2.0 * G * M / (C ** 2)


def horizon_area(M: float) -> float:
    return 4.0 * math.pi * Rs_of_M(M) ** 2


def grad_A_at_horizon(M: float) -> float:
    return 1.0 / Rs_of_M(M)


# --- Q10 / Q8 temperature ---

def T_Q10(M: float) -> float:
    """Phase-boundary temperature: k_B T = hbar c |grad A| / (4 pi)."""
    return HBAR * C * grad_A_at_horizon(M) / (4.0 * math.pi * KB)


# --- Path B: first-law integration ---

def dS_dM(M: float) -> float:
    """First law dE = T dS with E = M c^2 gives dS/dM = c^2 / T(M)."""
    return (C ** 2) / T_Q10(M)


def S_stam_from_first_law(M: float) -> float:
    """Closed-form integration of dS/dM from M' = 0 to M' = M.

    dS/dM = c^2 / T_Q10(M) = c^2 * 4 pi k_B Rs(M) / (hbar c)
          = 4 pi k_B c Rs(M) / hbar
          = 4 pi k_B c * 2 G M / (c^2 hbar)
          = 8 pi k_B G M / (hbar c)

    Integrating: S(M) = 4 pi k_B G M^2 / (hbar c).
    """
    return 4.0 * math.pi * KB * G * (M ** 2) / (HBAR * C)


def S_stam_numerical_integration(M: float, n_steps: int = 4000) -> float:
    """Numerically integrate dS/dM from 0 to M as a sanity check on the closed form."""
    M_grid = np.linspace(0.0, M, n_steps)
    # avoid singular point at M = 0
    M_grid = M_grid[1:]
    integrand = np.array([dS_dM(m) for m in M_grid])
    return float(np.trapezoid(integrand, M_grid))


# --- Bekenstein-Hawking reference ---

def S_bekenstein_hawking(M: float) -> float:
    """S_BH = k_B c^3 A / (4 hbar G) = k_B A / (4 ell_P^2)."""
    A = horizon_area(M)
    return KB * A / (4.0 * ELL_P_SQ)


# --- Tables ---

def build_entropy_table() -> pd.DataFrame:
    M_values = np.array([1e-6, 1e-3, 1.0, 1e6, 1e10]) * M_SUN
    rows = []
    for M in M_values:
        Rs = Rs_of_M(M)
        A_h = horizon_area(M)
        T = T_Q10(M)
        S_stam = S_stam_from_first_law(M)
        S_num = S_stam_numerical_integration(M)
        S_BH = S_bekenstein_hawking(M)
        S_over_A_natural = S_stam / (KB * A_h / ELL_P_SQ)  # should be 0.25
        rows.append({
            "M_over_Msun": M / M_SUN,
            "Rs_m": Rs,
            "horizon_area_m2": A_h,
            "T_Q10_K": T,
            "S_STAM_first_law_J_per_K": S_stam,
            "S_numerical_integral": S_num,
            "S_Bekenstein_Hawking": S_BH,
            "S_STAM / S_BH": S_stam / S_BH,
            "(S in natural units) / (A / 4)": S_over_A_natural / 0.25,
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_S_vs_M_comparison() -> Path:
    M_values = np.logspace(15, 35, 200)
    S_stam = np.array([S_stam_from_first_law(M) for M in M_values])
    S_BH = np.array([S_bekenstein_hawking(M) for M in M_values])

    plt.figure(figsize=(8, 5))
    plt.plot(M_values, S_stam, label="STAM Q11 (first law + Q10 T)", linewidth=2.5)
    plt.plot(M_values, S_BH, "--", label="Bekenstein-Hawking S = k_B A / (4 ell_P^2)", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Black hole mass M (kg)")
    plt.ylabel("Entropy S (J/K)")
    plt.title("STAM-derived entropy equals Bekenstein-Hawking exactly across all M")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q11_S_vs_M_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_S_vs_A() -> Path:
    """The headline plot: S is linear in A with slope k_B / (4 ell_P^2)."""
    M_values = np.logspace(15, 35, 200)
    A_values = np.array([horizon_area(M) for M in M_values])
    S_stam = np.array([S_stam_from_first_law(M) for M in M_values])

    plt.figure(figsize=(8, 5))
    plt.plot(A_values, S_stam, "o", markersize=3.5, label="STAM Q11 entropy", color="tab:blue")
    A_ref = np.logspace(np.log10(A_values.min()), np.log10(A_values.max()), 50)
    S_ref = KB * A_ref / (4.0 * ELL_P_SQ)
    plt.plot(A_ref, S_ref, "--", color="tab:red",
             label="Reference: S = k_B A / (4 ell_P^2)", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Horizon area A (m²)")
    plt.ylabel("Entropy S (J/K)")
    plt.title("S is exactly linear in A with the Bekenstein-Hawking slope — no fitted parameters")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q11_S_vs_A_linear.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_dS_dM_integration() -> Path:
    """Show the dS/dM curve and the area underneath (which IS the entropy)."""
    M_max = 5.0 * M_SUN
    M_grid = np.linspace(M_max / 5000.0, M_max, 5000)
    integrand = np.array([dS_dM(m) for m in M_grid])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(M_grid / M_SUN, 0, integrand, color="tab:blue", alpha=0.3,
                    label="Area = S(M)  (first law integration)")
    ax.plot(M_grid / M_SUN, integrand, color="tab:blue", linewidth=2,
            label="dS/dM = c^2 / T_Q10(M)  (linear in M)")
    ax.set_xlabel("Black hole mass M (solar masses)")
    ax.set_ylabel("dS/dM (J/K per kg)")
    ax.set_title("First-law integration: S(M) is the area under dS/dM")
    ax.legend()
    ax.grid(True, linewidth=0.3)

    out = PLOTS / "Q11_dS_dM_integration.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(case_df: pd.DataFrame, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A Q11: Bekenstein-Hawking Entropy from the Bubble Picture\n")
    md.append("## Purpose\n")
    md.append(
        "Derive the Bekenstein-Hawking entropy `S = k_B A / (4 ell_P^2)` for a black hole, "
        "from STAM-native ingredients: the bold-STAM bubble picture (no interior, all matter "
        "on the 2D phase boundary) and the Q10 phase-boundary temperature.\n"
    )
    md.append("## Two derivations, same answer\n")
    md.append("### Path A — physical reason (bubble picture)\n")
    md.append(
        "In bold STAM, A=1 is a 2D phase boundary between spacetime and not-spacetime. "
        "All matter that ever fell toward the black hole is on the boundary surface; there is "
        "no interior. So:\n"
        "- **No volume to fill.** The black hole has no 3D bulk for degrees of freedom to live "
        "in. Whatever entropy the system holds is on the 2D surface.\n"
        "- **Area scaling is automatic.** Entropy is proportional to area because area is the "
        "only place degrees of freedom can be.\n"
        "- **The 1/4 prefactor follows from holographic-bound saturation.** A 2D phase boundary "
        "with no interior cannot hold less entropy than the holographic bound; with no interior "
        "to share with, it must saturate the bound. The bound is `S <= k_B A / (4 ell_P^2)`, so "
        "saturation gives `S = k_B A / (4 ell_P^2)`.\n"
        "\n"
        "Author's plain-language statement: *area scaling because that's where everything is*. "
        "That sentence is the physical content of the Bekenstein-Hawking formula.\n"
    )
    md.append("### Path B — quantitative (first law + Q10 temperature)\n")
    md.append(
        "Apply the first law of thermodynamics `dE = T dS` to a black hole, with E = M c^2 and "
        "T from Q10:\n"
        "```text\n"
        "T_Q10 = hbar c |grad A| / (4 pi k_B)                     [Q8 / Q10 result]\n"
        "      = hbar c / (4 pi k_B Rs)                           [|grad A| = 1/Rs at horizon]\n"
        "\n"
        "dE = T dS  ->  dS = c^2 dM / T\n"
        "             = c^2 dM * 4 pi k_B Rs / (hbar c)\n"
        "             = 4 pi k_B c Rs / hbar  dM\n"
        "             = 8 pi k_B G M / (hbar c)  dM     [Rs = 2GM/c^2]\n"
        "\n"
        "Integrate from 0 to M:\n"
        "S(M) = 4 pi k_B G M^2 / (hbar c)\n"
        "\n"
        "Substitute A = 4 pi Rs^2 = 16 pi G^2 M^2 / c^4:\n"
        "S = k_B c^3 A / (4 hbar G)\n"
        "  = k_B A / (4 ell_P^2)                        [ell_P^2 = hbar G / c^3]\n"
        "```\n"
        "This is exactly the Bekenstein-Hawking formula. The 1/4 came out of the integration; "
        "it was not put in.\n"
    )
    md.append("## Where the 1/4 actually comes from\n")
    md.append(
        "The 1/4 in `S = A / 4` is the same structural factor as the 1/(4 pi) in the temperature "
        "`T = hbar c |grad A| / (4 pi)`. Both originate from the STAM gravity bridge "
        "`g = (c^2 / 2) grad A`, which contains the c^2/2 that propagates through to "
        "1/(2*2) = 1/4 in the entropy and 1/(2*2 pi) = 1/(4 pi) in the temperature.\n"
        "\n"
        "When the first-law integration is done, the 2 pi from the temperature combines with the "
        "2 pi from the horizon area (A = 4 pi Rs^2) to give a factor of 8 pi in dS/dM, and the "
        "integration brings in another 1/2 (since integral of M dM = M^2/2), giving 4 pi G k_B / "
        "(hbar c) * M^2. Re-expressed in terms of A, the prefactor is 1/4.\n"
    )
    md.append("## Numerical verification\n")
    md.append(case_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n")
    md.append(
        "**Reading the table:** `S_STAM / S_BH = 1.0` and "
        "`(S in natural units) / (A / 4) = 1.0` across all mass scales. Bold STAM reproduces "
        "Bekenstein-Hawking entropy exactly. The numerical integral matches the closed-form "
        "result, confirming the integration was done correctly.\n"
    )
    md.append("## What this gives bold STAM\n")
    md.append(
        "Bold STAM now has a complete black hole thermodynamics chain:\n"
        "1. **Temperature** `T = hbar c |grad A| / (4 pi k_B)` — Q8 postulate, Q10 derivation.\n"
        "2. **Spectrum** Planckian at T — Q10 derivation (Bose-Einstein on equilibrated phase boundary).\n"
        "3. **Entropy** `S = k_B A / (4 ell_P^2)` — Q11 derivation (bubble picture + first law).\n"
        "4. **Energy bookkeeping** dE = T dS — first law, used in Q11.\n"
        "\n"
        "All four are reproduced exactly relative to standard Hawking results, but for "
        "different physical reasons: phase-boundary equilibrium and bubble surface, not "
        "Schwarzschild Wick rotation and QFT mode counting.\n"
    )
    md.append(
        "**For the README:** entropy and thermodynamic-consistency status now move from open "
        "to derived. Bold STAM has a closed-form BH thermodynamics that matches the standard "
        "results numerically while having a fundamentally different mechanism (no interior, "
        "phase boundary as the only source of degrees of freedom).\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "Q11_bekenstein_hawking_entropy_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    case_df = build_entropy_table()
    case_csv = RESULTS / "Q11_entropy_table.csv"
    case_df.to_csv(case_csv, index=False)

    plot_paths = [
        plot_S_vs_M_comparison(),
        plot_S_vs_A(),
        plot_dS_dM_integration(),
    ]

    summary = write_markdown(case_df, plot_paths)

    print("STAM Model-A Q11: Bekenstein-Hawking Entropy from the Bubble Picture")
    print("=" * 72)
    print(f"\nPlanck length squared (ell_P^2): {ELL_P_SQ:.6g} m^2")
    print(f"Planck length (ell_P):           {math.sqrt(ELL_P_SQ):.6g} m")
    print()
    print("Case results (S in J/K):")
    print(case_df.to_string(index=False))
    print(f"\nS_STAM / S_BH ratio: {case_df['S_STAM / S_BH'].min():.6g} to {case_df['S_STAM / S_BH'].max():.6g}")
    print("(should be exactly 1.0 — Q11 reproduces Bekenstein-Hawking)")
    print("\nFiles written:")
    print(f"- {case_csv}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
