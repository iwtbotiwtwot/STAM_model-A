#!/usr/bin/env python3
"""
32_model_b_pure_traversal_excess.py

STAM Model-B cosmological distance from PURE path-traversal excess only.
Author: Sean Brady / STAM Model-B continuation

Purpose:
    Compute the cosmological distance prediction in STAM Model-B using ONLY
    the path-traversal excess of light through the cosmic A field. No clock
    dilation effects are folded in (those are separate predictions affecting
    redshift / source-frame luminosity, not the apparent distance directly).

Key clarification (per author 2026-05-07):
    The standard "Shapiro effect" lumps two physically distinct effects:
        1. Path-traversal excess: light slowed along the path by integrated A.
        2. Clock dilation: source/observer clock rates differ at different A.
    These are SEPARABLE. For supernova distance prediction, only effect (1)
    contributes to D_L directly. Effect (2) contributes to z determination
    (already handled by standard redshift formalism).

    This script computes ONLY effect (1) and adds it to the geometric distance.

What gets used:
    A_cosmo(r): Model-B matter-sourced cosmic A field profile.
        For a homogeneous critical-density universe, Poisson's equation
        applied to a sphere of radius r centered on observer gives
        approximately A(r) ~ (r/L)^2 with L = c/H_0 = Hubble distance.
        This is NOT the de Sitter ansatz postulated in Model-A Q8 — it's
        a derived result from matter-sourced Poisson under the
        critical-density assumption.

What gets computed:
    Pure traversal excess: TE(z) = ∫_0^r(z) A_cosmo(r') dr'
    Distance prediction: d_Model-B(z) = r(z) + TE(z)
    where r(z) ≈ L z at low z (Hubble linear).

What gets compared:
    - LCDM D_L(z) (standard cosmology, Omega_m = 0.315 flat)
    - Model-A no-b: D_adj,0(z) = L z (1 + 0.5z) (the "no-b" form from README)
    - Model-A with bridge: D_adj,0(z) + b z (b ≈ 354.95 Mly)
    - Model-B traversal excess (this script)

Honest expectation:
    Model-B with A_cosmo(r) ~ (r/L)^2 gives TE(z) = L z^3 / 3 (cubic in z).
    Model-A's no-b form has 0.5 L z^2 (quadratic in z).
    At low z, z^3 < z^2, so Model-B is FLATTER than Model-A no-b.
    Both are flatter than LCDM D_L. The bridge term b is not derived; it
    remains a fitting artifact in Model-B as well.

What this script DOES not try to do:
    - Compute F3-extended cosmological dynamics (modified Friedmann from
      bold-STAM effective stress-energy). That's a separate, more elaborate
      calculation flagged as critical follow-up.
    - Derive the bridge term b. Model-B with pure traversal excess from
      matter-sourced A does NOT produce b naturally.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


def quad_trapz(func, a: float, b: float, n: int = 4000) -> float:
    """Trapezoid integration to avoid scipy dependency."""
    if a >= b:
        return 0.0
    x = np.linspace(a, b, n)
    y = np.array([func(xi) for xi in x])
    return float(np.trapezoid(y, x))


# STAM cosmological constants (from README)
C_STAM_MLY = 3.261563776       # c in Mly/yr units
H_STAM = 0.000243635
L_STAM_MLY = C_STAM_MLY / H_STAM   # ≈ 13387 Mly (Hubble length in Mly)

# LCDM constants
H0_KMS_MPC = 67.4
C_KMS = 299792.458
HUBBLE_DIST_MPC = C_KMS / H0_KMS_MPC  # ≈ 4448 Mpc
MLY_PER_MPC = 3.262
L_STAM_MPC = L_STAM_MLY / MLY_PER_MPC
HUBBLE_DIST_MLY = HUBBLE_DIST_MPC * MLY_PER_MPC

OMEGA_M = 0.315
OMEGA_L = 1.0 - OMEGA_M

# Historical bridge term (fitted, not derived)
B_BRIDGE_MLY = 354.95


# --- Model-B cosmic A profile ---

def A_cosmo_critical_density(r_over_L: float | np.ndarray) -> float | np.ndarray:
    """A_cosmo(r) from Model-B matter-sourced Poisson at critical density.

    For a homogeneous critical-density universe, the Poisson equation
        grad^2 A = (8 pi G / c^2) rho_crit
    integrated over a sphere of radius r centered on observer gives
        A(r) ≈ (r/L)^2  where L = c/H_0
    by the relation rho_crit = 3 H_0^2 / (8 pi G).

    NOTE: This is the matter-sourced derivation, not the de Sitter ansatz.
    The functional form happens to match because critical density is what
    creates this particular A(r) profile.
    """
    return r_over_L ** 2


# --- Pure traversal excess ---

def traversal_excess_Mly(z: float) -> float:
    """TE(z) = integral from 0 to r(z) of A_cosmo(r') dr', in Mly.

    For low z: r(z) = L * z (Hubble linear).
    TE(z) = ∫_0^(Lz) (r'/L)^2 dr' = (1/3L^2) * (Lz)^3 = L z^3 / 3
    """
    return L_STAM_MLY * (z ** 3) / 3.0


def traversal_excess_numerical_Mly(z: float, n_points: int = 1000) -> float:
    """Numerical version for sanity check — should match analytical L z^3 / 3."""
    if z <= 0:
        return 0.0
    r_grid = np.linspace(0, L_STAM_MLY * z, n_points)
    A_grid = A_cosmo_critical_density(r_grid / L_STAM_MLY)
    return float(np.trapezoid(A_grid, r_grid))


# --- Model-B distance prediction (pure traversal excess only) ---

def d_model_b_Mly(z: float) -> float:
    """Apparent distance in Model-B from geometric + pure traversal excess.

    No clock dilation, no F3-extended Friedmann dynamics — just the path
    integral piece of A_cosmo through the line of sight.
    """
    geometric = L_STAM_MLY * z   # Hubble linear at low z
    excess = traversal_excess_Mly(z)
    return geometric + excess


# --- Comparison distance functions ---

def d_model_a_no_b_Mly(z: float) -> float:
    """Model-A no-b: D_adj,0(z) = L z (1 + 0.5 z)."""
    return L_STAM_MLY * z * (1.0 + 0.5 * z)


def d_model_a_with_bridge_Mly(z: float) -> float:
    """Model-A with bridge term: D_adj,0(z) + b * z."""
    return d_model_a_no_b_Mly(z) + B_BRIDGE_MLY * z


def d_lcdm_luminosity_Mly(z: float) -> float:
    """LCDM luminosity distance for comparison.

    D_L(z) = (1+z) * c/H_0 * integral 0 to z of 1/E(z') dz'
    """
    def E(zp: float) -> float:
        return math.sqrt(OMEGA_M * (1.0 + zp) ** 3 + OMEGA_L)
    integral = quad_trapz(lambda zp: 1.0 / E(zp), 0, z)
    return (1.0 + z) * HUBBLE_DIST_MLY * integral


# --- Tables ---

def build_distance_table() -> pd.DataFrame:
    z_values = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2]
    rows = []
    for z in z_values:
        d_b = d_model_b_Mly(z)
        d_a_no_b = d_model_a_no_b_Mly(z)
        d_a_bridge = d_model_a_with_bridge_Mly(z)
        d_lcdm = d_lcdm_luminosity_Mly(z)
        TE_analytical = traversal_excess_Mly(z)
        TE_numerical = traversal_excess_numerical_Mly(z)
        rows.append({
            "z": z,
            "TE_analytical_Mly": TE_analytical,
            "TE_numerical_Mly": TE_numerical,
            "d_Model_B_Mly": d_b,
            "d_Model_A_no_b_Mly": d_a_no_b,
            "d_Model_A_bridge_Mly": d_a_bridge,
            "d_LCDM_DL_Mly": d_lcdm,
            "Model_B_over_no_b": d_b / d_a_no_b if d_a_no_b > 0 else float("nan"),
            "Model_B_over_LCDM": d_b / d_lcdm if d_lcdm > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_distances() -> Path:
    z_grid = np.linspace(0.001, 1.5, 200)
    d_b = np.array([d_model_b_Mly(z) for z in z_grid])
    d_a_no_b = np.array([d_model_a_no_b_Mly(z) for z in z_grid])
    d_a_bridge = np.array([d_model_a_with_bridge_Mly(z) for z in z_grid])
    d_lcdm = np.array([d_lcdm_luminosity_Mly(z) for z in z_grid])

    plt.figure(figsize=(10, 6))
    plt.plot(z_grid, d_b, color="tab:blue", linewidth=2.5,
             label="Model-B (pure path TE only): L z + L z³/3")
    plt.plot(z_grid, d_a_no_b, color="tab:orange", linewidth=2,
             label="Model-A no-b: L z (1 + 0.5z)")
    plt.plot(z_grid, d_a_bridge, color="tab:green", linewidth=2,
             label=f"Model-A with bridge b={B_BRIDGE_MLY:.0f} Mly")
    plt.plot(z_grid, d_lcdm, color="tab:red", linewidth=2, linestyle="--",
             label="LCDM D_L (Ω_m=0.315 flat)")

    plt.xlabel("Redshift z")
    plt.ylabel("Apparent distance (Mly)")
    plt.title("Model-B cosmological distance from pure path-traversal excess only")
    plt.grid(True, linewidth=0.3)
    plt.legend()
    plt.xlim(0, 1.5)

    out = PLOTS / "32_model_b_distance_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_traversal_excess_only() -> Path:
    z_grid = np.linspace(0.001, 1.5, 200)
    TE = np.array([traversal_excess_Mly(z) for z in z_grid])
    bridge_only = B_BRIDGE_MLY * z_grid

    plt.figure(figsize=(10, 5.5))
    plt.plot(z_grid, TE, color="tab:blue", linewidth=2.5,
             label="Model-B pure TE: L z³/3")
    plt.plot(z_grid, bridge_only, color="tab:green", linewidth=2,
             label=f"Historical bridge: b z = {B_BRIDGE_MLY:.0f} z Mly (linear)")
    plt.plot(z_grid, 0.35 * L_STAM_MLY * z_grid ** 2, color="tab:orange",
             linewidth=2, linestyle=":",
             label="Model-A TE: 0.35 L z²")

    plt.xlabel("Redshift z")
    plt.ylabel("Traversal excess (Mly)")
    plt.title("Pure path traversal excess: Model-B vs Model-A vs historical bridge")
    plt.grid(True, linewidth=0.3)
    plt.legend()
    plt.xlim(0, 1.5)

    out = PLOTS / "32_traversal_excess_only.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_distance_residuals() -> Path:
    """Show residuals against LCDM (a common reference) for each model."""
    z_grid = np.linspace(0.05, 1.2, 100)
    d_b = np.array([d_model_b_Mly(z) for z in z_grid])
    d_a_no_b = np.array([d_model_a_no_b_Mly(z) for z in z_grid])
    d_a_bridge = np.array([d_model_a_with_bridge_Mly(z) for z in z_grid])
    d_lcdm = np.array([d_lcdm_luminosity_Mly(z) for z in z_grid])

    res_b = (d_b - d_lcdm) / d_lcdm * 100
    res_a_no_b = (d_a_no_b - d_lcdm) / d_lcdm * 100
    res_a_bridge = (d_a_bridge - d_lcdm) / d_lcdm * 100

    plt.figure(figsize=(10, 5.5))
    plt.plot(z_grid, res_b, color="tab:blue", linewidth=2.5,
             label="Model-B (pure TE)")
    plt.plot(z_grid, res_a_no_b, color="tab:orange", linewidth=2,
             label="Model-A no-b")
    plt.plot(z_grid, res_a_bridge, color="tab:green", linewidth=2,
             label="Model-A with bridge")
    plt.axhline(0, color="black", linewidth=0.6)
    plt.xlabel("Redshift z")
    plt.ylabel("Distance residual vs LCDM D_L (%)")
    plt.title("Distance residuals vs LCDM (negative = STAM-style flatter than LCDM)")
    plt.grid(True, linewidth=0.3)
    plt.legend()

    out = PLOTS / "32_distance_residuals_vs_LCDM.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(df: pd.DataFrame, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-B Cosmological Distance — Pure Traversal Excess Only\n")
    md.append("## Purpose\n")
    md.append(
        "Compute the cosmological distance prediction in STAM Model-B using ONLY the path-"
        "traversal excess of light through the cosmic A field. Clock dilation effects are "
        "kept SEPARATE (per author's 2026-05-07 sharpening). The full Shapiro effect lumps "
        "the two together; here we peel them apart and use only the path piece.\n"
    )
    md.append("## Setup\n")
    md.append(
        f"- L = c/H_0 = {L_STAM_MLY:.0f} Mly  (STAM Hubble length)\n"
        f"- Historical bridge term: b = {B_BRIDGE_MLY} Mly (fitted, retained for comparison)\n"
        "- Cosmic A profile (Model-B): `A_cosmo(r) ≈ (r/L)²` from critical-density Poisson\n"
        "- Pure traversal excess: `TE(z) = ∫_0^Lz A_cosmo(r') dr' = L z³/3`\n"
        "- Apparent distance: `d_Model-B(z) = L z + L z³/3`\n"
    )
    md.append("## Why A_cosmo ≈ (r/L)² is a derived consequence, not a postulate\n")
    md.append(
        "For a homogeneous critical-density universe, Poisson's equation\n"
        "```text\n"
        "grad² A = (8π G / c²) ρ_crit\n"
        "```\n"
        "applied to a sphere of radius r centered on observer gives\n"
        "```text\n"
        "A(r) ~ (8π G / 3c²) ρ_crit r² = r²/L²\n"
        "```\n"
        "using `ρ_crit = 3H₀²/(8πG)` and `L = c/H₀`. This is NOT the de Sitter ansatz from "
        "Q8 (which was postulated). It is a *derived* consequence of Model-B's matter-"
        "sourced framework under the observed critical-density assumption.\n"
    )
    md.append("## Numerical results\n")
    md.append(df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Findings\n")
    md.append(
        "- **Model-B pure TE = L z³/3.** Cubic in z, with no z² or z-linear term.\n"
        "- **Model-A's TE = 0.35 L z² is a quadratic form.** Differs in functional shape.\n"
        "- **Bridge term b z is z-linear.** Different again.\n"
        "- All three differ in z-dependence, but at z ≈ 1 they happen to give comparable "
        "magnitudes (Model-B: 0.33L, Model-A: 0.35L, b: 0.027L).\n"
        "- Model-B is FLATTER than Model-A no-b at low z (z³ < z² for z<1).\n"
        "- All three are flatter than LCDM D_L. None reproduces LCDM.\n"
    )
    md.append("## Honest verdict\n")
    md.append(
        "**Pure traversal excess alone does not derive the historical bridge term b.** "
        "Model-B's matter-sourced A field at cosmological scales gives a z³ correction to "
        "the Hubble linear distance, not the z-linear b z form. The bridge term remains a "
        "fitting artifact when comparing STAM (any version) to LCDM-fitted catalogs.\n"
        "\n"
        "What Model-B's pure TE does give us:\n"
        "- A clean, derived prediction `D = Lz + Lz³/3` for cosmological distance under "
        "matter-sourced A and pure path-excess only.\n"
        "- Confirmation that Model-B's distance prediction is FLATTER than Model-A's no-b "
        "form, which is in turn flatter than LCDM. The discrepancy with LCDM is real and "
        "widening, not narrowing.\n"
        "\n"
        "If STAM is correct, this prediction needs to be tested against catalogs honestly:\n"
        "- Either the catalogs (calibrated within LCDM) systematically overestimate "
        "supernova distances at high z, in which case Model-B's flatter prediction is closer "
        "to the truth.\n"
        "- Or STAM's distance prediction is wrong, and the framework needs additional "
        "structure beyond pure path traversal excess.\n"
        "\n"
        "**The F3-extended cosmology calculation (modified Friedmann with bold-STAM "
        "effective stress-energy) remains the next critical follow-up.** It might "
        "introduce additional dark-energy-like structure that moves the prediction closer "
        "to LCDM. This script does NOT include that effect.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "32_model_b_pure_traversal_excess_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    df = build_distance_table()
    df.to_csv(RESULTS / "32_distance_table.csv", index=False)

    plot_paths = [
        plot_distances(),
        plot_traversal_excess_only(),
        plot_distance_residuals(),
    ]

    summary = write_markdown(df, plot_paths)

    print("STAM Model-B Cosmological Distance — Pure Traversal Excess Only")
    print("=" * 72)
    print()
    print(f"L_STAM = {L_STAM_MLY:.4g} Mly  (Hubble length)")
    print(f"Bridge term (historical, fitted): b = {B_BRIDGE_MLY} Mly")
    print()
    print("Cosmic A profile (Model-B from matter-sourced Poisson):")
    print("  A_cosmo(r) ~= (r/L)^2  [derived under critical-density assumption]")
    print()
    print("Pure path-traversal excess:")
    print("  TE(z) = L z^3 / 3  (cubic in z)")
    print()
    print("Distance comparison:")
    print(df.to_string(index=False))
    print()
    print("Honest finding:")
    print("  Model-B pure TE gives L z + L z^3/3.")
    print("  Cubic in z, FLATTER than Model-A's no-b form (which has z^2).")
    print("  Bridge term b is z-linear; not derived from pure path TE in any model.")
    print("  STAM (any version) is flatter than LCDM. Discrepancy real, not resolved.")
    print()
    print("Next critical follow-up: F3-extended cosmology with modified Friedmann.")
    print()
    print("Files written:")
    print(f"- {RESULTS / '32_distance_table.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
