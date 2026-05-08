#!/usr/bin/env python3
"""
34_F3_extended_cosmology.py

F3-extended cosmological calculation for STAM Model-B.
Author: Sean Brady / STAM Model-B continuation 2026-05-07

Goal:
    Compute the bold-STAM effective stress-energy at cosmological scales
    (energy density AND pressure components), determine the effective
    equation of state, and check whether it produces a derived value of
    A_0 ≈ 0.0265 — i.e., whether the F3-extended cosmology naturally
    derives the historical bridge term b ≈ 354.95 Mly from first
    principles.

Setup:
    Cosmic A profile (Model-B, matter-sourced, critical density):
        A_cosmo(r) ≈ (r/L)² for r within Hubble distance L.

    Bold-STAM metric (Framework C):
        g_tt = -(1-A) c²
        g_rr = 1/[(1-A)(1-A²)²]

    For this metric, compute the Einstein tensor components G^t_t,
    G^r_r, G^θ_θ. Interpreting these as effective stress-energy under
    Einstein gravity:
        G^t_t = -8πG ρ_eff / c²              (energy density, signed)
        G^r_r =  8πG p_r_eff / c⁴            (radial pressure)
        G^θ_θ =  8πG p_t_eff / c⁴            (tangential pressure)

    Cosmic averages over the Hubble sphere:
        <ρ_eff>, <p_r_eff>, <p_t_eff>
    Effective equation of state:
        w_eff = <p_average> / <ρc²>

    Dark-energy-like behavior corresponds to w ≈ -1.

Honest expectations:
    - F3 already showed ρ_eff has a sign change at A ≈ 0.44 (negative
      for A < 0.44, positive for A > 0.44). For cosmic profile with most
      of the volume at small A, the volume-weighted average will be
      dominated by the negative-ρ regions.
    - If the volume-weighted <ρ> is negative AND the pressure has
      appropriate sign, the effective stress-energy could behave like
      dark energy with w ≈ -1.
    - If so, this could produce a cosmological-constant-equivalent that
      modifies expansion dynamics and gives a Friedmann correction.
    - Whether this exactly reproduces A_0 = 0.0265 is the question.

Caveats (be honest):
    - This script uses the static-spherical metric formulas applied to a
      cosmic A profile. The true cosmological framework would require
      bold-STAM in FRW geometry, which has not been derived. So this is
      a "best estimate from what we have" — order-of-magnitude
      cosmological prediction, not a fully rigorous calculation.
    - The averaging assumes A_cosmo(r) = (r/L)² as a derived consequence
      of critical-density Poisson. If actual STAM cosmology has a
      different A profile, results would shift.
    - Pressure components G^r_r, G^θ_θ require numerical second
      derivatives of the metric coefficients. Numerical noise at small
      grid spacings could affect results; we use enough points to be safe.
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

# Constants
G_NEWTON = 6.67430e-11
C_LIGHT = 2.99792458e8

# Bridge term and STAM length (for comparison)
B_HISTORICAL_MLY = 354.95
L_STAM_MLY = 13387.09
A_0_EMPIRICAL = B_HISTORICAL_MLY / L_STAM_MLY  # ≈ 0.02651


# --- Bold-STAM metric coefficients ---

def h_metric(A: float) -> float:
    """g_tt = -h(A) c² in bold STAM. h(A) = 1-A."""
    return 1.0 - A


def k_metric(A: float) -> float:
    """g_rr = 1/k(A) in bold STAM. k(A) = (1-A)(1-A²)²."""
    return (1.0 - A) * (1.0 - A ** 2) ** 2


def Phi_metric(A: float) -> float:
    """Phi defined by g_tt = -e^(2Phi) c²: Phi = (1/2) ln(1-A)."""
    if A >= 1.0:
        return -math.inf
    return 0.5 * math.log(1.0 - A)


def Lambda_metric(A: float) -> float:
    """Lambda defined by g_rr = e^(2Lambda): Lambda = -(1/2) ln k(A)."""
    if A >= 1.0:
        return math.inf
    return -0.5 * math.log(k_metric(A))


# --- Cosmic A profile (Model-B, critical-density Poisson) ---

def A_cosmo(r_over_L: float) -> float:
    """A_cosmo(r) = (r/L)^2 from critical-density Poisson result.

    Bounded at A < 1, so r/L < 1 (within Hubble sphere).
    """
    return r_over_L ** 2


# --- Numerical Einstein tensor at cosmic A profile ---

def compute_einstein_tensor(u_grid: np.ndarray) -> dict:
    """Compute G^μ_ν components on grid of u = r/L values.

    Uses bold-STAM metric with A = u² cosmic profile.
    Returns dict with G_t_t, G_r_r, G_theta_theta arrays.
    All values in units where L = 1 (so 1/r² has units 1/L² in physical).
    """
    A_grid = u_grid ** 2
    Phi_grid = np.array([Phi_metric(A) for A in A_grid])
    Lambda_grid = np.array([Lambda_metric(A) for A in A_grid])

    # Coordinate is u = r/L; need d/dr = (1/L) × d/du.
    # In dimensionless units L=1, dr = du.
    du = u_grid[1] - u_grid[0]
    Phi_p = np.gradient(Phi_grid, du)         # dPhi/du = dPhi/dr × L; we'll set L=1
    Lambda_p = np.gradient(Lambda_grid, du)
    Phi_pp = np.gradient(Phi_p, du)

    e_neg_2Lambda = np.exp(-2.0 * Lambda_grid)
    r = u_grid  # in units of L

    # Avoid singularity at r=0
    r_safe = np.where(r > 1e-10, r, 1e-10)

    G_t_t = (1.0 - e_neg_2Lambda) / r_safe ** 2 + 2.0 * e_neg_2Lambda * Lambda_p / r_safe
    G_r_r = (1.0 - e_neg_2Lambda) / r_safe ** 2 - 2.0 * e_neg_2Lambda * Phi_p / r_safe
    G_theta_theta = -e_neg_2Lambda * (Phi_pp + Phi_p ** 2 - Phi_p * Lambda_p
                                      + (Phi_p - Lambda_p) / r_safe)

    return {
        "u_grid": u_grid,
        "A_grid": A_grid,
        "Phi": Phi_grid,
        "Lambda": Lambda_grid,
        "Phi_p": Phi_p,
        "Lambda_p": Lambda_p,
        "G_t_t": G_t_t,
        "G_r_r": G_r_r,
        "G_theta_theta": G_theta_theta,
    }


# --- Effective stress-energy interpretation ---

def effective_stress_energy(einstein: dict) -> dict:
    """Interpret Einstein tensor as effective T^μ_ν of dark-energy-like fluid.

    Under Einstein gravity G^μ_ν = (8πG/c⁴) T^μ_ν with T^μ_ν = diag(-ρc², p_r, p_t, p_t):
        G_t_t = -8πG ρ / c²    →  ρ_eff (relative units) ∝ -G_t_t
        G_r_r =  8πG p_r / c⁴  →  p_r_eff (relative units) ∝ G_r_r
        G_theta_theta = G_phi_phi  →  p_t_eff (relative units) ∝ G_theta_theta

    Returns relative units (proportional to true values).
    Equation of state w = p / (ρ c²) = -G_pressure / G_t_t for our convention.
    """
    rho_eff = -einstein["G_t_t"]  # proportional to ρc² (positive when G_t_t < 0)
    p_r_eff = einstein["G_r_r"]
    p_t_eff = einstein["G_theta_theta"]
    p_avg = (p_r_eff + 2.0 * p_t_eff) / 3.0  # isotropic average pressure

    # Avoid divide-by-zero where ρ_eff is small
    safe_rho = np.where(np.abs(rho_eff) > 1e-15, rho_eff, np.nan)
    w_eff = p_avg / safe_rho

    return {
        "rho_eff": rho_eff,
        "p_r_eff": p_r_eff,
        "p_t_eff": p_t_eff,
        "p_avg_eff": p_avg,
        "w_eff": w_eff,
    }


# --- Cosmic volume averaging ---

def cosmic_volume_averages(u_grid: np.ndarray, stress: dict, u_max: float = 0.95) -> dict:
    """Average effective stress-energy over Hubble sphere (volume-weighted).

    Volume element: 4π r² dr. In dimensionless u: 4π u² du.
    Average X = (3/u_max³) × ∫₀^u_max X(u) u² du.
    Avoid u very close to 1 where bold-STAM metric coefficients diverge.
    """
    mask = u_grid <= u_max
    u_sub = u_grid[mask]
    rho = stress["rho_eff"][mask]
    p_r = stress["p_r_eff"][mask]
    p_t = stress["p_t_eff"][mask]
    p_avg = stress["p_avg_eff"][mask]

    weights = u_sub ** 2  # volume weight
    norm = np.trapezoid(weights, u_sub)

    rho_avg = np.trapezoid(rho * weights, u_sub) / norm
    p_r_avg = np.trapezoid(p_r * weights, u_sub) / norm
    p_t_avg = np.trapezoid(p_t * weights, u_sub) / norm
    p_avg_avg = np.trapezoid(p_avg * weights, u_sub) / norm

    w_avg_total = p_avg_avg / rho_avg if abs(rho_avg) > 1e-15 else float("nan")

    return {
        "rho_cosmic_avg": float(rho_avg),
        "p_r_cosmic_avg": float(p_r_avg),
        "p_t_cosmic_avg": float(p_t_avg),
        "p_avg_cosmic_avg": float(p_avg_avg),
        "w_cosmic_avg": float(w_avg_total),
        "u_max_used": u_max,
    }


# --- Bridge term derivation from cosmic stress-energy ---

def bridge_from_cosmic_stress_energy(cosmic_avg: dict) -> dict:
    """Estimate the bridge-term-equivalent A_0 from cosmic average stress-energy.

    If bold-STAM cosmic average ρ_eff is negative and dark-energy-like (w ≈ -1),
    it acts as effective dark energy. The 'equivalent constant ambient A_0'
    that this would produce can be estimated.

    For a dark-energy-like component contributing fraction Ω_DE,STAM to cosmic
    energy budget:
        A_0_predicted = some-function-of(Ω_DE,STAM)

    Simplest estimate: A_0 ~ |ρ_eff_cosmic_avg| × geometric factor.

    The relative magnitude of the cosmic-averaged effective ρ (in units where
    the bracket function from F3 is order 1) gives:
        bracket_avg_volumetric = <bracket(A)> averaged over cosmic sphere
    where bracket(A) = 1 - (1-A)²(1+A)(1+A+4A²) — the F3 ρ-bracket.

    The empirical A_0 = 0.0265 should match this if F3-extended cosmology
    derives the bridge term.
    """
    # Computing A_0 from cosmic average effective stress-energy is non-trivial;
    # depends on how the effective stress-energy modifies cosmological dynamics.
    # As a rough indicator: A_0 ∝ |rho_cosmic_avg| × (some calibration factor).
    # The "size" of the dark-energy-like effect.
    abs_rho_avg = abs(cosmic_avg["rho_cosmic_avg"])
    return {
        "abs_rho_cosmic_avg": abs_rho_avg,
        "w_cosmic_avg": cosmic_avg["w_cosmic_avg"],
        "is_dark_energy_like": abs(cosmic_avg["w_cosmic_avg"] + 1.0) < 0.5,
        "comment": "abs_rho_cosmic_avg gives the relative magnitude of effective DE; "
                   "exact A_0 derivation requires modified-Friedmann solution.",
    }


# --- Tables ---

def build_pointwise_table(einstein: dict, stress: dict) -> pd.DataFrame:
    """Show G^μ_ν components and effective stress-energy at specific A values."""
    u_targets = [0.1, 0.3, 0.5, 0.7, 0.9]
    rows = []
    for u_t in u_targets:
        idx = np.argmin(np.abs(einstein["u_grid"] - u_t))
        rows.append({
            "u": einstein["u_grid"][idx],
            "A": einstein["A_grid"][idx],
            "G_t_t": einstein["G_t_t"][idx],
            "G_r_r": einstein["G_r_r"][idx],
            "G_theta_theta": einstein["G_theta_theta"][idx],
            "rho_eff": stress["rho_eff"][idx],
            "p_avg_eff": stress["p_avg_eff"][idx],
            "w_eff": stress["w_eff"][idx],
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_einstein_components(einstein: dict, stress: dict) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    u_grid = einstein["u_grid"]
    A_grid = einstein["A_grid"]
    mask = u_grid <= 0.99  # avoid the singular endpoint

    axes[0, 0].plot(u_grid[mask], stress["rho_eff"][mask], color="tab:blue", linewidth=2)
    axes[0, 0].axhline(0, color="black", linewidth=0.5)
    axes[0, 0].set_xlabel("u = r/L")
    axes[0, 0].set_ylabel("ρ_eff (relative units)")
    axes[0, 0].set_title("Effective energy density vs cosmic position")
    axes[0, 0].grid(True, linewidth=0.3)

    axes[0, 1].plot(u_grid[mask], stress["p_r_eff"][mask], color="tab:orange",
                    linewidth=2, label="p_r_eff (radial)")
    axes[0, 1].plot(u_grid[mask], stress["p_t_eff"][mask], color="tab:green",
                    linewidth=2, label="p_t_eff (tangential)")
    axes[0, 1].plot(u_grid[mask], stress["p_avg_eff"][mask], color="tab:red",
                    linewidth=1.5, linestyle="--", label="p_avg")
    axes[0, 1].axhline(0, color="black", linewidth=0.5)
    axes[0, 1].set_xlabel("u = r/L")
    axes[0, 1].set_ylabel("p_eff (relative units)")
    axes[0, 1].set_title("Effective pressure components")
    axes[0, 1].grid(True, linewidth=0.3)
    axes[0, 1].legend()

    axes[1, 0].plot(A_grid[mask], stress["rho_eff"][mask], color="tab:blue", linewidth=2,
                    label="ρ_eff")
    axes[1, 0].plot(A_grid[mask], stress["p_avg_eff"][mask], color="tab:red", linewidth=2,
                    label="p_avg_eff")
    axes[1, 0].axhline(0, color="black", linewidth=0.5)
    axes[1, 0].set_xlabel("A = u² (cosmic A field)")
    axes[1, 0].set_ylabel("Stress-energy components")
    axes[1, 0].set_title("ρ and p vs A — F3 sign-change at A ≈ 0.44 visible")
    axes[1, 0].grid(True, linewidth=0.3)
    axes[1, 0].legend()

    # w_eff (clip extreme values)
    w_clipped = np.clip(stress["w_eff"][mask], -5, 5)
    axes[1, 1].plot(u_grid[mask], w_clipped, color="tab:purple", linewidth=2)
    axes[1, 1].axhline(-1, color="tab:red", linestyle="--", linewidth=1,
                       label="w = -1 (cosmological constant)")
    axes[1, 1].axhline(0, color="black", linewidth=0.5)
    axes[1, 1].set_xlabel("u = r/L")
    axes[1, 1].set_ylabel("w_eff = p_avg / (ρ c²)")
    axes[1, 1].set_title("Effective equation of state — dark-energy-like requires w ≈ -1")
    axes[1, 1].set_ylim(-5, 5)
    axes[1, 1].grid(True, linewidth=0.3)
    axes[1, 1].legend()

    fig.suptitle("F3-extended cosmology: bold-STAM effective stress-energy", y=1.00, fontsize=13)

    out = PLOTS / "34_einstein_stress_energy.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    return out


def plot_volume_averages(cosmic_avg: dict) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axis("off")
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 5)

    ax.text(4.5, 4.5, "F3-extended cosmology — cosmic averages",
            ha="center", fontsize=14, fontweight="bold")

    ax.text(0.5, 3.7, "Cosmic-volume average over Hubble sphere "
                       f"(u_max = {cosmic_avg['u_max_used']}):",
            fontsize=11)
    ax.text(0.5, 3.1, f"<ρ_eff>     = {cosmic_avg['rho_cosmic_avg']:+.6f}",
            fontsize=11, family="monospace")
    ax.text(0.5, 2.6, f"<p_r_eff>   = {cosmic_avg['p_r_cosmic_avg']:+.6f}",
            fontsize=11, family="monospace")
    ax.text(0.5, 2.1, f"<p_t_eff>   = {cosmic_avg['p_t_cosmic_avg']:+.6f}",
            fontsize=11, family="monospace")
    ax.text(0.5, 1.6, f"<p_avg>     = {cosmic_avg['p_avg_cosmic_avg']:+.6f}",
            fontsize=11, family="monospace")
    ax.text(0.5, 1.0, f"w_cosmic    = {cosmic_avg['w_cosmic_avg']:+.6f}",
            fontsize=12, family="monospace", fontweight="bold", color="tab:purple")
    ax.text(0.5, 0.5, "Dark-energy-like requires w ≈ -1.", fontsize=10, style="italic")

    out = PLOTS / "34_cosmic_averages.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(pointwise_df: pd.DataFrame, cosmic_avg: dict, bridge_summary: dict,
                   plot_paths: list[Path]) -> Path:
    md = []
    md.append("# F3-Extended Cosmology — Bold-STAM Effective Stress-Energy at Cosmic Scale\n")
    md.append("## Purpose\n")
    md.append(
        "Compute the bold-STAM effective stress-energy (energy density AND pressure components) "
        "at cosmological scales. Determine whether the effective equation of state is dark-"
        "energy-like (w ≈ -1) and check whether this produces a derived value of `A_0 ≈ 0.0265` "
        "matching the historical bridge term b ≈ 354.95 Mly.\n"
    )
    md.append("## Setup\n")
    md.append(
        "- Cosmic A profile (Model-B): `A_cosmo(r) = (r/L)²` from critical-density Poisson.\n"
        "- Bold-STAM metric: `g_tt = -(1-A)c²`, `g_rr = 1/[(1-A)(1-A²)²]`.\n"
        "- Compute Einstein tensor components G^t_t, G^r_r, G^θ_θ for this metric.\n"
        "- Interpret as effective stress-energy under Einstein gravity:\n"
        "  ρ_eff ∝ -G^t_t, p_r_eff ∝ G^r_r, p_t_eff ∝ G^θ_θ.\n"
        "- Volume-average over Hubble sphere to get cosmic-mean stress-energy.\n"
        "- Effective equation of state: w_eff = <p_avg> / <ρc²>.\n"
    )
    md.append("## Pointwise stress-energy at sample cosmic positions\n")
    md.append(pointwise_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n## Cosmic volume averages\n")
    md.append("```text\n")
    md.append(f"<ρ_eff>      = {cosmic_avg['rho_cosmic_avg']:+.6f} (relative units)")
    md.append(f"<p_r_eff>    = {cosmic_avg['p_r_cosmic_avg']:+.6f}")
    md.append(f"<p_t_eff>    = {cosmic_avg['p_t_cosmic_avg']:+.6f}")
    md.append(f"<p_avg>      = {cosmic_avg['p_avg_cosmic_avg']:+.6f}")
    md.append(f"w_cosmic_avg = {cosmic_avg['w_cosmic_avg']:+.6f}")
    md.append(f"u_max used   = {cosmic_avg['u_max_used']} (avoiding singularity at u=1)")
    md.append("```\n")
    md.append("**Dark-energy-like reference**: w = -1 (cosmological constant).\n")
    if bridge_summary["is_dark_energy_like"]:
        md.append("**Verdict on equation of state**: dark-energy-like (within ±0.5 of w=-1). "
                  "Suggests bold-STAM cosmology naturally produces a dark-energy-like effective "
                  "fluid component.\n")
    else:
        md.append(f"**Verdict on equation of state**: NOT dark-energy-like (w_avg = "
                  f"{cosmic_avg['w_cosmic_avg']:.3f}, far from -1). Bold-STAM effective "
                  "stress-energy at cosmic scales does not reduce to a simple dark-energy "
                  "component.\n")
    md.append("\n## Bridge-term derivation status\n")
    md.append(
        "Whether this cosmic stress-energy *exactly* reproduces `A_0 = 0.0265` requires "
        "solving modified Friedmann equations with the effective T^μ_ν as input, then "
        "extracting the distance-redshift relation, then identifying the constant-ambient-A "
        "equivalent. That last step is non-trivial because the modification doesn't always "
        "reduce to a simple constant-A picture.\n"
        "\n"
        "What this script can say:\n"
        "- The bold-STAM effective stress-energy at cosmic scale has a specific equation of "
        "state w (see above).\n"
        "- The magnitude of the cosmic-average ρ_eff sets the SCALE of any dark-energy-like "
        "modification.\n"
        "- Whether the resulting Friedmann modification produces A_0 = 0.0265 specifically "
        "is a separate, more elaborate calculation.\n"
        "\n"
        "**Honest interpretation:**\n"
        f"- Cosmic-average ρ_eff magnitude: {bridge_summary['abs_rho_cosmic_avg']:.6f} (relative units)\n"
        f"- Cosmic-average w: {bridge_summary['w_cosmic_avg']:+.6f}\n"
        f"- {bridge_summary['comment']}\n"
    )
    md.append("\n## Caveats\n")
    md.append(
        "- Static-spherical-around-observer approximation. True FRW cosmology with bold-STAM "
        "would require deriving the Model-B field equations in expanding-universe context.\n"
        "- Used A_cosmo(r) = (r/L)² as derived consequence of critical-density Poisson. If "
        "actual cosmological A profile differs, results differ.\n"
        "- Volume averaging up to u_max < 1 to avoid singular endpoint. Results depend "
        "weakly on u_max choice (most volume is at moderate u).\n"
        "- Pressure components computed from numerical second derivatives — some numerical "
        "noise possible at very small/large u.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "34_F3_extended_cosmology_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    # Build grid of u = r/L values, avoiding singular endpoints
    u_grid = np.linspace(0.001, 0.99, 5000)

    einstein = compute_einstein_tensor(u_grid)
    stress = effective_stress_energy(einstein)

    cosmic_avg = cosmic_volume_averages(u_grid, stress, u_max=0.95)
    bridge_summary = bridge_from_cosmic_stress_energy(cosmic_avg)

    pointwise_df = build_pointwise_table(einstein, stress)
    pointwise_df.to_csv(RESULTS / "34_pointwise_stress_energy.csv", index=False)

    cosmic_df = pd.DataFrame([cosmic_avg])
    cosmic_df.to_csv(RESULTS / "34_cosmic_averages.csv", index=False)

    plot_paths = [
        plot_einstein_components(einstein, stress),
        plot_volume_averages(cosmic_avg),
    ]

    summary = write_markdown(pointwise_df, cosmic_avg, bridge_summary, plot_paths)

    print("F3-Extended Cosmology — Bold-STAM Effective Stress-Energy")
    print("=" * 64)
    print()
    print("Pointwise stress-energy at sample cosmic positions:")
    print(pointwise_df.to_string(index=False))
    print()
    print("Cosmic-volume averages (over u in [0, 0.95]):")
    print(f"  <rho_eff>      = {cosmic_avg['rho_cosmic_avg']:+.6f}  (relative units)")
    print(f"  <p_r_eff>      = {cosmic_avg['p_r_cosmic_avg']:+.6f}")
    print(f"  <p_t_eff>      = {cosmic_avg['p_t_cosmic_avg']:+.6f}")
    print(f"  <p_avg>        = {cosmic_avg['p_avg_cosmic_avg']:+.6f}")
    print(f"  w_cosmic_avg   = {cosmic_avg['w_cosmic_avg']:+.6f}")
    print(f"  Dark-energy-like (|w+1|<0.5)? {bridge_summary['is_dark_energy_like']}")
    print()
    print(f"Bridge term comparison:")
    print(f"  Empirical A_0 = b/L = {A_0_EMPIRICAL:.6f}")
    print(f"  Magnitude of cosmic <rho_eff>: {bridge_summary['abs_rho_cosmic_avg']:.6f}")
    print()
    print("Files written:")
    print(f"- {RESULTS / '34_pointwise_stress_energy.csv'}")
    print(f"- {RESULTS / '34_cosmic_averages.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
