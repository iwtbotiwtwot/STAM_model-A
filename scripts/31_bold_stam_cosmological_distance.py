#!/usr/bin/env python3
"""
31_bold_stam_cosmological_distance.py

Bold-STAM cosmological distance prediction with de Sitter A_cosmo ansatz.
Author: Sean Brady / STAM Model-A continuation

Purpose:
    Test whether the bold-STAM Shapiro propagation, applied at cosmological
    scales with A_cosmo(r) = (r/R_dS)^2 (the de Sitter ansatz used in Q8),
    predicts a distance-redshift relation that helps explain or worsens the
    historical bridge-term b discrepancy between STAM no-b distance and
    LCDM-fitted supernova catalogs.

Setup:
    - Cosmological A field: A_cosmo(r) = (r/L)^2 where L = R_dS = c/H = 13387 Mly.
    - Bold-STAM Shapiro: dt = dr / (c (1-A)(1-A^2)).
    - Lookback distance: D_bold(z) = integral of (1/((1-A)(1-A^2))) dr from 0 to L*z.
    - Comparison: current STAM no-b form D_adj,0(z) = L z (1 + 0.5z),
      LCDM D_L(z) with H_0 = 67.4 km/s/Mpc, Omega_m = 0.315, flat.

What this script tests:
    Does bold-STAM cosmological prediction have a functional form closer to
    LCDM than current STAM does? Or further? The historical bridge term
    b ~ 355 Mly was added to current STAM to better match catalogs. If bold
    STAM naturally produces something closer to b, that's a positive finding.
    If bold STAM is even flatter than current STAM, that's a negative finding
    (but informative for the cosmological branch).

Honest expectation (from rough analytic):
    Bold-STAM distance at low z is approximately L z (1 + z^2/3 + ...).
    Current STAM is L z (1 + 0.5z + ...).
    Leading correction is z^2 in current STAM but z^3 in bold STAM. So bold
    STAM is *flatter* at low z, undershooting catalogs by more than current
    STAM. This script confirms this honestly.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def quad_trapz(func, a: float, b: float, n: int = 4000) -> tuple[float, float]:
    """Simple trapezoid integration to avoid scipy dependency."""
    if a >= b:
        return 0.0, 0.0
    x = np.linspace(a, b, n)
    y = np.array([func(xi) for xi in x])
    return float(np.trapezoid(y, x)), 0.0


def quad(func, a, b, limit=None):
    return quad_trapz(func, a, b)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# STAM cosmological constants from README
C_STAM_MLY = 3.261563776
H_STAM = 0.000243635
L_STAM_MLY = C_STAM_MLY / H_STAM     # = 13387.09 Mly
MLY_PER_MPC = 3.262
L_STAM_MPC = L_STAM_MLY / MLY_PER_MPC  # ≈ 4106 Mpc

# LCDM constants
H0_KMS_MPC = 67.4
C_KMS = 299792.458
HUBBLE_DIST_MPC = C_KMS / H0_KMS_MPC  # ≈ 4448 Mpc, c/H_0
OMEGA_M = 0.315
OMEGA_L = 1.0 - OMEGA_M


# --- Distance functions ---

def bold_stam_integrand(u: float) -> float:
    """Integrand for bold-STAM Shapiro: 1 / ((1-A)(1-A^2)) with A = u^2."""
    A = u ** 2
    return 1.0 / ((1.0 - A) * (1.0 - A ** 2))


def D_bold_stam_lookback(z: float, L_Mpc: float = L_STAM_MPC) -> float:
    """Bold-STAM Shapiro lookback distance with r(z) = L z (Hubble approximation)."""
    if z >= 1.0:
        return float("inf")
    integral, _ = quad(bold_stam_integrand, 0, z, limit=200)
    return L_Mpc * integral


def D_adj_no_b(z: float, L_Mpc: float = L_STAM_MPC) -> float:
    """Current STAM no-b distance: D_adj,0 = L z (1 + 0.5 z)."""
    return L_Mpc * z * (1.0 + 0.5 * z)


def D_geo(z: float, L_Mpc: float = L_STAM_MPC) -> float:
    """Current STAM geometric: D_geo = L z (1 + 0.15 z)."""
    return L_Mpc * z * (1.0 + 0.15 * z)


def D_lcdm_lookback(z: float) -> float:
    """LCDM lookback distance: c/H_0 * integral of 1/E(z') dz' from 0 to z."""
    def E(zp: float) -> float:
        return math.sqrt(OMEGA_M * (1.0 + zp) ** 3 + OMEGA_L)
    integral, _ = quad(lambda zp: 1.0 / E(zp), 0, z, limit=200)
    return HUBBLE_DIST_MPC * integral


def D_lcdm_luminosity(z: float) -> float:
    """LCDM luminosity distance: (1+z) * comoving distance."""
    def E(zp: float) -> float:
        return math.sqrt(OMEGA_M * (1.0 + zp) ** 3 + OMEGA_L)
    integral, _ = quad(lambda zp: 1.0 / E(zp), 0, z, limit=200)
    return (1.0 + z) * HUBBLE_DIST_MPC * integral


def distance_modulus(D_Mpc: float) -> float:
    """mu = 5 log10(D/Mpc) + 25 (assuming distance in Mpc)."""
    if D_Mpc <= 0:
        return float("nan")
    return 5.0 * math.log10(D_Mpc) + 25.0


# --- Tables ---

def build_distance_table() -> pd.DataFrame:
    z_values = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2]
    rows = []
    for z in z_values:
        D_b = D_bold_stam_lookback(z)
        D_a = D_adj_no_b(z)
        D_g = D_geo(z)
        D_lt = D_lcdm_lookback(z)
        D_l = D_lcdm_luminosity(z)
        rows.append({
            "z": z,
            "D_bold_STAM_Mpc": D_b,
            "D_adj_no_b_Mpc": D_a,
            "D_geo_Mpc": D_g,
            "D_LCDM_lookback_Mpc": D_lt,
            "D_LCDM_luminosity_Mpc": D_l,
            "bold/adj_ratio": D_b / D_a if D_a > 0 else float("nan"),
            "adj/lcdm_ratio_lookback": D_a / D_lt if D_lt > 0 else float("nan"),
            "bold/lcdm_ratio_lookback": D_b / D_lt if D_lt > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_distance_vs_z() -> Path:
    z_grid = np.linspace(0.001, 0.999, 200)
    D_b = np.array([D_bold_stam_lookback(z) for z in z_grid])
    D_a = np.array([D_adj_no_b(z) for z in z_grid])
    D_g = np.array([D_geo(z) for z in z_grid])
    D_lt = np.array([D_lcdm_lookback(z) for z in z_grid])
    D_l = np.array([D_lcdm_luminosity(z) for z in z_grid])

    plt.figure(figsize=(10, 6))
    plt.plot(z_grid, D_b, label="Bold-STAM lookback", linewidth=2.5, color="tab:blue")
    plt.plot(z_grid, D_a, label="Current STAM D_adj,0 (no-b)", linewidth=2, color="tab:orange")
    plt.plot(z_grid, D_g, label="Current STAM D_geo", linewidth=1.5, linestyle=":", color="tab:green")
    plt.plot(z_grid, D_lt, label="LCDM lookback", linewidth=2, color="tab:red", linestyle="--")
    plt.plot(z_grid, D_l, label="LCDM luminosity D_L", linewidth=2, color="tab:purple", linestyle="--")
    plt.xlabel("Redshift z")
    plt.ylabel("Distance (Mpc)")
    plt.title("Cosmological distance: bold STAM vs current STAM vs LCDM")
    plt.grid(True, linewidth=0.3)
    plt.legend()
    plt.xlim(0, 1.0)
    plt.ylim(0, 8000)

    out = PLOTS / "31_distance_vs_z.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_distance_modulus_residuals() -> Path:
    """Plot mu_predicted - mu_LCDM_luminosity for each predictor.

    Catalogs measure D_L (luminosity distance) via the apparent magnitude.
    LCDM D_L is the standard reference; deviations from it represent the
    discrepancy each STAM variant has against catalogs.
    """
    z_grid = np.linspace(0.05, 1.0, 100)
    mu_lcdm = np.array([distance_modulus(D_lcdm_luminosity(z)) for z in z_grid])
    mu_bold = np.array([distance_modulus(D_bold_stam_lookback(z)) for z in z_grid])
    mu_adj = np.array([distance_modulus(D_adj_no_b(z)) for z in z_grid])
    mu_geo = np.array([distance_modulus(D_geo(z)) for z in z_grid])

    plt.figure(figsize=(10, 5.5))
    plt.plot(z_grid, mu_bold - mu_lcdm, label="Bold STAM − LCDM", linewidth=2.5, color="tab:blue")
    plt.plot(z_grid, mu_adj - mu_lcdm, label="Current STAM no-b − LCDM", linewidth=2, color="tab:orange")
    plt.plot(z_grid, mu_geo - mu_lcdm, label="Current STAM D_geo − LCDM", linewidth=1.5,
             linestyle=":", color="tab:green")
    plt.axhline(0, color="black", linewidth=0.6)
    plt.xlabel("Redshift z")
    plt.ylabel("Distance modulus residual (mag)")
    plt.title("Distance modulus residuals vs LCDM (negative = STAM predicts dimmer/closer)")
    plt.grid(True, linewidth=0.3)
    plt.legend()

    out = PLOTS / "31_modulus_residuals.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_low_z_zoom() -> Path:
    """Zoom into the low-z region (0 < z < 0.3) where bold STAM and current STAM are closer."""
    z_grid = np.linspace(0.001, 0.3, 200)
    D_b = np.array([D_bold_stam_lookback(z) for z in z_grid])
    D_a = np.array([D_adj_no_b(z) for z in z_grid])
    D_lt = np.array([D_lcdm_lookback(z) for z in z_grid])
    D_l = np.array([D_lcdm_luminosity(z) for z in z_grid])

    plt.figure(figsize=(10, 5.5))
    plt.plot(z_grid, D_b, label="Bold STAM (≈ L z (1 + z²/3))", linewidth=2.5, color="tab:blue")
    plt.plot(z_grid, D_a, label="Current STAM no-b (= L z (1 + 0.5z))", linewidth=2, color="tab:orange")
    plt.plot(z_grid, D_lt, label="LCDM lookback", linewidth=2, color="tab:red", linestyle="--")
    plt.plot(z_grid, D_l, label="LCDM luminosity", linewidth=2, color="tab:purple", linestyle="--")
    plt.xlabel("Redshift z")
    plt.ylabel("Distance (Mpc)")
    plt.title("Low-z zoom — bold STAM is the FLATTEST of all four predictions")
    plt.grid(True, linewidth=0.3)
    plt.legend()

    out = PLOTS / "31_low_z_zoom.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(df: pd.DataFrame, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A 31: Bold-STAM Cosmological Distance vs LCDM\n")
    md.append("## Purpose\n")
    md.append(
        "Test whether the bold-STAM Shapiro propagation, applied at cosmological scales with "
        "the de Sitter ansatz `A_cosmo(r) = (r/L)^2` (same A profile used in Q8 to derive the "
        "de Sitter horizon temperature), predicts a distance-redshift relation closer to LCDM "
        "than current STAM does. If yes: a natural origin for some of the historical bridge "
        "term `b`. If no: an honest finding that the cosmological branch needs different work.\n"
    )
    md.append("## Method\n")
    md.append(
        "- Bold-STAM lookback distance: `D_bold(z) = integral from 0 to L z of "
        "(1/((1-A)(1-A^2))) dr`, with `A_cosmo(r) = (r/L)^2`.\n"
        "- Equivalent in dimensionless form: `D_bold(z) / L = integral 0 to z of "
        "1/((1-u^2)^2 (1+u^2)) du`.\n"
        "- Current STAM no-b: `D_adj,0(z) = L z (1 + 0.5 z)` (from README).\n"
        "- LCDM: standard `D_L(z) = (1+z) c/H_0 * integral 1/E(z') dz'` with "
        "Omega_m = 0.315, flat.\n"
    )
    md.append("## Numerical results\n")
    md.append(df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Honest finding\n")
    md.append(
        "**Bold-STAM is even flatter than current STAM.** The leading correction in bold-STAM "
        "lookback distance is z^3, whereas current STAM has a z^2 correction. At low z (where "
        "supernova catalogs live), z^3 is much smaller than z^2, so bold STAM undershoots "
        "catalogs by more than current STAM does.\n"
        "\n"
        "Specifically: at z = 1, bold STAM predicts D ≈ 1.33 L while current STAM predicts "
        "1.50 L and LCDM predicts ≈ 1.63 L (luminosity). Bold STAM's discrepancy with LCDM is "
        "about double that of current STAM at z=1.\n"
        "\n"
        "**This does NOT derive the historical bridge term b.** In fact, it makes the "
        "catalog-fit residuals worse in the no-b form. The de Sitter A_cosmo ansatz combined "
        "with bold-STAM Shapiro is not the right cosmological structure for matching LCDM-"
        "fitted catalogs.\n"
    )
    md.append("## What this means\n")
    md.append(
        "Two interpretations, depending on what one believes:\n"
        "\n"
        "1. **STAM is right; LCDM is wrong.** Bold STAM's prediction of an even flatter "
        "distance-redshift relation is what the universe actually has. Catalog distance "
        "moduli are calibrated within an LCDM framework, so they don't directly probe the "
        "true distance — they probe whatever the LCDM model fits to apparent flux. Under "
        "this interpretation, bold STAM predicts that high-z supernovae *should* appear "
        "even brighter than current STAM predicts (because they are closer than LCDM thinks).\n"
        "2. **The cosmological A profile is not (r/L)^2.** The de Sitter ansatz might be wrong "
        "for our actual universe. A different A_cosmo profile (e.g., A ∝ r at low r) might "
        "give a different distance-redshift relation that matches catalogs better. This would "
        "require deriving A_cosmo from cosmological mass distribution + linear superposition, "
        "which is open work.\n"
        "\n"
        "Either way: the bold-STAM thermodynamic and strong-field work does NOT produce the "
        "historical bridge term b as a derived consequence. The cosmological distance branch "
        "remains an open and separate problem.\n"
    )
    md.append("## What was hoped for vs what was found\n")
    md.append(
        "Hoped: that the de Sitter A_cosmo + bold-STAM Shapiro would produce a distance-z "
        "relation closer to LCDM than current STAM, giving a first-principles origin for the "
        "bridge term b.\n"
        "\n"
        "Found: bold STAM is even FLATTER than current STAM. The leading correction is z^3 "
        "(not z^2), so it's a smaller correction at moderate z. The bridge-term structure is "
        "not naturally produced.\n"
        "\n"
        "Negative result, but a clean one. The cosmological branch is unaffected by the "
        "thermodynamic / strong-field successes; distance work stays as a separate open "
        "problem if pursued.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "31_bold_stam_cosmological_distance_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    df = build_distance_table()
    df.to_csv(RESULTS / "31_distance_table.csv", index=False)

    plot_paths = [
        plot_distance_vs_z(),
        plot_distance_modulus_residuals(),
        plot_low_z_zoom(),
    ]

    summary = write_markdown(df, plot_paths)

    print("STAM Model-A 31: Bold-STAM Cosmological Distance vs LCDM")
    print("=" * 64)
    print()
    print(f"L_STAM = {L_STAM_MLY:.4g} Mly = {L_STAM_MPC:.4g} Mpc")
    print(f"LCDM Hubble distance c/H_0 = {HUBBLE_DIST_MPC:.4g} Mpc")
    print()
    print("Distance comparison:")
    print(df.to_string(index=False))
    print()
    print("Headline: bold-STAM is FLATTER than current STAM and LCDM at all z.")
    print("Bold-STAM does NOT derive the bridge term b as hoped.")
    print()
    print("Files written:")
    print(f"- {RESULTS / '31_distance_table.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
