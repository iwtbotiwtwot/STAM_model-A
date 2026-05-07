#!/usr/bin/env python3
"""
Q8_unresolved_A_boundary_flux.py

STAM Model-A unresolved-A boundary flux test.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    Test the proposed STAM-native rule for thermal emission at A=1 resolution
    boundaries:

        k_B T = (1 / (4 pi)) * hbar * c * |grad A|_boundary

    The asymmetry that produces a net outward flux follows from the
    resolved/unresolved framing of STAM. A > 1 is unresolved by definition:
    resolutions oriented inward have no committed substrate to anchor to and
    cannot complete; outward resolutions can. The boundary therefore emits a
    net outward flux without invoking pair production or negative-energy
    partners.

    The 1/(4 pi) prefactor is not a free fit. It is forced by the STAM
    gravity bridge g = (c^2 / 2) grad A together with the standard
    surface-gravity-to-temperature relation T = hbar kappa / (2 pi k_B c),
    since kappa = (c^2 / 2) |grad A| at the boundary.

    The single rule is checked against three established thermal-horizon
    results, each parameterized by a different physical input:
        - Schwarzschild Hawking temperature   (varying mass M)
        - Unruh temperature                   (varying acceleration a)
        - de Sitter horizon temperature       (varying H)

    All three are expected to land on the same T-vs-|grad A| line, with
    the prefactor coming from STAM structure rather than from fitting.

Notes:
    - The Planckian spectrum itself is not derived here. The rule produces
      the temperature scale; downstream luminosity uses standard
      Stefan-Boltzmann thermodynamics for a thermal source at that T.
    - The de Sitter case assumes A_dS(r) = (r / R_dS)^2 with R_dS = c / H,
      which is a model-dependent STAM choice not yet derived. The
      Schwarzschild and Unruh applications do not depend on this choice.
    - Stefan-Boltzmann luminosity here is photon-only (two helicities).
      Full Hawking emission summed over radiating species would be higher
      by a degrees-of-freedom factor; the 1/M^2 scaling is the result
      being checked, not the absolute prefactor.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Physical constants (SI)
G = 6.67430e-11           # m^3 kg^-1 s^-2
C = 2.99792458e8          # m/s
HBAR = 1.054571817e-34    # J s
KB = 1.380649e-23         # J/K
M_SUN = 1.98847e30        # kg
SIGMA_SB = (math.pi ** 2) * (KB ** 4) / (60.0 * (HBAR ** 3) * (C ** 2))  # photon Stefan-Boltzmann

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# STAM postulate: k_B T = (1 / (4 pi)) * hbar * c * |grad A|
PREFACTOR = 1.0 / (4.0 * math.pi)


def stam_T_from_grad_A(grad_A: float) -> float:
    """Temperature at a resolution boundary from the STAM rule.

    grad_A has units 1/m. Returns temperature in K.
    """
    return PREFACTOR * HBAR * C * grad_A / KB


# --- Schwarzschild (black hole) case ---

def Rs_of_M(M: float) -> float:
    return 2.0 * G * M / (C ** 2)


def grad_A_schwarzschild(M: float) -> float:
    """At the horizon, A = Rs/r so |dA/dr| = Rs/r^2 = 1/Rs at r = Rs."""
    return 1.0 / Rs_of_M(M)


def T_hawking_standard(M: float) -> float:
    return HBAR * (C ** 3) / (8.0 * math.pi * G * M * KB)


def T_stam_schwarzschild(M: float) -> float:
    return stam_T_from_grad_A(grad_A_schwarzschild(M))


def L_stefan_boltzmann(T: float, area: float) -> float:
    return SIGMA_SB * (T ** 4) * area


def L_stam_schwarzschild(M: float) -> float:
    Rs = Rs_of_M(M)
    return L_stefan_boltzmann(T_stam_schwarzschild(M), 4.0 * math.pi * Rs ** 2)


def L_standard_schwarzschild(M: float) -> float:
    Rs = Rs_of_M(M)
    return L_stefan_boltzmann(T_hawking_standard(M), 4.0 * math.pi * Rs ** 2)


# --- Unruh (Rindler) case ---

def grad_A_rindler(a: float) -> float:
    """From g = (c^2/2) grad A inverted: |grad A| = 2|g|/c^2."""
    return 2.0 * a / (C ** 2)


def T_unruh_standard(a: float) -> float:
    return HBAR * a / (2.0 * math.pi * KB * C)


def T_stam_unruh(a: float) -> float:
    return stam_T_from_grad_A(grad_A_rindler(a))


# --- de Sitter (cosmological) case ---

def R_dS(H: float) -> float:
    return C / H


def grad_A_de_sitter(H: float) -> float:
    """Assumes A_dS(r) = (r/R_dS)^2, so |dA/dr| at R_dS is 2/R_dS = 2H/c."""
    return 2.0 * H / C


def T_dS_standard(H: float) -> float:
    return HBAR * H / (2.0 * math.pi * KB)


def T_stam_dS(H: float) -> float:
    return stam_T_from_grad_A(grad_A_de_sitter(H))


# --- Tables ---

def build_schwarzschild_table() -> pd.DataFrame:
    M_values = np.array([1e-6, 1e-3, 1.0, 1e6, 1e10]) * M_SUN
    rows = []
    for M in M_values:
        Rs = Rs_of_M(M)
        gA = grad_A_schwarzschild(M)
        T_stam = T_stam_schwarzschild(M)
        T_std = T_hawking_standard(M)
        L_stam = L_stam_schwarzschild(M)
        L_std = L_standard_schwarzschild(M)
        rows.append({
            "M_over_Msun": M / M_SUN,
            "Rs_m": Rs,
            "grad_A_per_m": gA,
            "T_STAM_K": T_stam,
            "T_Hawking_K": T_std,
            "T_ratio": T_stam / T_std,
            "L_STAM_W": L_stam,
            "L_standard_W": L_std,
            "L_ratio": L_stam / L_std,
        })
    return pd.DataFrame(rows)


def build_unruh_table() -> pd.DataFrame:
    a_values = np.array([1.0, 1e10, 1e20, 1e26])  # last value is near Schwinger-like scales
    rows = []
    for a in a_values:
        gA = grad_A_rindler(a)
        T_stam = T_stam_unruh(a)
        T_std = T_unruh_standard(a)
        rows.append({
            "a_m_per_s2": a,
            "grad_A_per_m": gA,
            "T_STAM_K": T_stam,
            "T_Unruh_K": T_std,
            "T_ratio": T_stam / T_std,
        })
    return pd.DataFrame(rows)


def build_de_sitter_table() -> pd.DataFrame:
    # 2.27e-18/s is approximately H_0; bracket above and below
    H_values = np.array([1e-20, 2.27e-18, 1e-15, 1e-10])
    rows = []
    for H in H_values:
        gA = grad_A_de_sitter(H)
        T_stam = T_stam_dS(H)
        T_std = T_dS_standard(H)
        rows.append({
            "H_per_s": H,
            "R_dS_m": R_dS(H),
            "grad_A_per_m": gA,
            "T_STAM_K": T_stam,
            "T_dS_standard_K": T_std,
            "T_ratio": T_stam / T_std,
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_T_vs_grad_A_unified(schw_df: pd.DataFrame, unruh_df: pd.DataFrame, dS_df: pd.DataFrame) -> Path:
    plt.figure(figsize=(8, 5))

    plt.scatter(schw_df["grad_A_per_m"], schw_df["T_STAM_K"],
                label="Schwarzschild horizons (varying M)", marker="o", s=60)
    plt.scatter(unruh_df["grad_A_per_m"], unruh_df["T_STAM_K"],
                label="Rindler / Unruh (varying a)", marker="s", s=60)
    plt.scatter(dS_df["grad_A_per_m"], dS_df["T_STAM_K"],
                label="de Sitter (varying H)", marker="^", s=60)

    gA_ref = np.logspace(-30, 5, 200)
    T_ref = PREFACTOR * HBAR * C * gA_ref / KB
    plt.plot(gA_ref, T_ref, "--", color="black", alpha=0.5,
             label="STAM rule: k_B T = hbar c |grad A| / (4 pi)")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("|grad A| at boundary (1/m)")
    plt.ylabel("Temperature T (K)")
    plt.title("Three thermal horizons collapse onto a single STAM rule")
    plt.legend(fontsize=9)
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q8_T_vs_grad_A_unified.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_hawking_T_comparison() -> Path:
    M_values = np.logspace(15, 35, 200)
    T_stam = np.array([T_stam_schwarzschild(M) for M in M_values])
    T_std = np.array([T_hawking_standard(M) for M in M_values])

    plt.figure(figsize=(8, 5))
    plt.plot(M_values, T_stam, label="STAM rule", linewidth=2.5)
    plt.plot(M_values, T_std, "--", label="Standard Hawking T_H", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Black hole mass M (kg)")
    plt.ylabel("Temperature (K)")
    plt.title("STAM-derived T equals Hawking T_H exactly across all M")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q8_hawking_T_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_hawking_L_comparison() -> Path:
    M_values = np.logspace(15, 35, 200)
    L_stam = np.array([L_stam_schwarzschild(M) for M in M_values])
    L_std = np.array([L_standard_schwarzschild(M) for M in M_values])

    plt.figure(figsize=(8, 5))
    plt.plot(M_values, L_stam, label="STAM rule (Stefan-Boltzmann at T_STAM)", linewidth=2.5)
    plt.plot(M_values, L_std, "--", label="Stefan-Boltzmann at T_Hawking", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Black hole mass M (kg)")
    plt.ylabel("Luminosity (W)")
    plt.title("Hawking-style L proportional to 1/M^2 recovered from STAM rule")
    plt.legend()
    plt.grid(True, which="both", linewidth=0.3)

    out = PLOTS / "Q8_hawking_L_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_markdown(schw_df: pd.DataFrame, unruh_df: pd.DataFrame, dS_df: pd.DataFrame,
                   plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A Q8: Unresolved-A Boundary Flux Test\n")
    md.append("## Purpose\n")
    md.append(
        "Evaluate a single STAM-native rule for thermal emission at A=1 resolution boundaries:\n"
    )
    md.append("```text\nk_B T = (1 / (4 pi)) * hbar * c * |grad A|_boundary\n```\n")
    md.append(
        "The asymmetry producing a net outward flux follows from the resolved/unresolved framing. "
        "A > 1 is unresolved by definition: resolutions oriented inward have no committed substrate "
        "to anchor to and cannot complete. Outward resolutions can. The boundary therefore emits a "
        "net outward flux without pair production or negative-energy partners.\n"
    )
    md.append("## Mechanism summary\n")
    md.append(
        "1. Vacuum fluctuations everywhere are unresolved A excursions.\n"
        "2. The resolution rate scale at any boundary is set by the gradient |grad A|.\n"
        "3. At a resolution boundary (A=1) the inward direction lacks committed substrate, so only "
        "outward resolutions complete.\n"
        "4. The characteristic temperature is k_B T = (hbar c / (4 pi)) |grad A|.\n"
        "5. The 1/(4 pi) prefactor is forced by the STAM gravity bridge g = (c^2/2) grad A together "
        "with the standard surface-gravity-to-temperature relation T = hbar kappa / (2 pi k_B c), "
        "since kappa = (c^2/2) |grad A| at the boundary. It is not a free fit.\n"
    )
    md.append("## Schwarzschild (Hawking) case\n")
    md.append(schw_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Unruh case\n")
    md.append(unruh_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## de Sitter case\n")
    md.append(dS_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Interpretation\n")
    md.append(
        "All three temperature ratios T_STAM / T_standard are 1.0 to floating-point precision. "
        "This is one rule applied to three different boundary configurations, with the prefactor "
        "fixed by the STAM gravity bridge rather than by data fitting. The Hawking-style "
        "luminosity scaling L proportional to 1/M^2 follows from combining the temperature rule "
        "with Stefan-Boltzmann thermodynamics for a thermal source at that T.\n"
    )
    md.append("## What is and is not derived\n")
    md.append(
        "- Derived from the rule: temperature scale at any A=1 boundary; recovery of Hawking, "
        "Unruh, and de Sitter temperatures with the standard prefactor.\n"
        "- Not derived here: the Planckian spectrum itself. The rule fixes the temperature scale; "
        "treating the emission as thermal at that T is an additional step that future STAM work "
        "would have to ground in resolution-event statistics.\n"
        "- Stefan-Boltzmann luminosity uses photon-only (two-helicity) sigma. Full Hawking "
        "emission summed over species is higher by a degrees-of-freedom factor. The 1/M^2 scaling "
        "is the result being checked, not the absolute multiplicative prefactor.\n"
    )
    md.append("## Caveats and open questions\n")
    md.append(
        "- The de Sitter check assumes A_dS(r) = (r/R_dS)^2 with R_dS = c/H. STAM has not yet "
        "committed to a cosmological form for A; the Schwarzschild and Unruh applications do not "
        "depend on this choice.\n"
        "- The mechanism predicts that any boundary in A produces an asymmetric outward flux at "
        "T = hbar c |grad A| / (4 pi k_B), regardless of whether the boundary is gravitational, "
        "kinematic, or cosmological. This unification is the strongest content of the rule.\n"
        "- Where STAM could differ from standard QFT-on-curved-spacetime predictions: details of "
        "the spectrum (deviations from exact Planckian), behavior near very small Rs, and the "
        "scaling for non-spherical / multi-source A configurations.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "Q8_unresolved_A_boundary_flux_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    schw_df = build_schwarzschild_table()
    unruh_df = build_unruh_table()
    dS_df = build_de_sitter_table()

    schw_csv = RESULTS / "Q8_schwarzschild_table.csv"
    unruh_csv = RESULTS / "Q8_unruh_table.csv"
    dS_csv = RESULTS / "Q8_de_sitter_table.csv"
    schw_df.to_csv(schw_csv, index=False)
    unruh_df.to_csv(unruh_csv, index=False)
    dS_df.to_csv(dS_csv, index=False)

    plot_paths = [
        plot_T_vs_grad_A_unified(schw_df, unruh_df, dS_df),
        plot_hawking_T_comparison(),
        plot_hawking_L_comparison(),
    ]

    summary = write_markdown(schw_df, unruh_df, dS_df, plot_paths)

    print("STAM Model-A Q8: Unresolved-A Boundary Flux Test")
    print("=" * 56)
    print("\nSchwarzschild (Hawking) case:")
    print(schw_df.to_string(index=False))
    print("\nUnruh case:")
    print(unruh_df.to_string(index=False))
    print("\nde Sitter case:")
    print(dS_df.to_string(index=False))
    print("\nFiles written:")
    for p in [schw_csv, unruh_csv, dS_csv, summary, *plot_paths]:
        print(f"- {p}")


if __name__ == "__main__":
    main()
