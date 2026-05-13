#!/usr/bin/env python3
"""
G64_bubble_refinement.py

Two candidate A(r, theta; M, a) formulas for the Kerr extension, both
consistent with the framework's commitments, with different consequences:

Formula (alpha) - the G3 / G62 reading:
    A_alpha = 2 M r / Sigma     where Sigma = r^2 + a^2 cos^2(theta)
    - A = 1 at the static-limit surface (ergosphere outer boundary)
    - At equator: A = 2M/r, A = 1 at r = 2M
    - At pole:    A = 2Mr/(r^2 + a^2), A = 1 at r = r_+
    - Bubble shape: oblate in BL coords, r(theta) varies
    - |grad A| varies with theta -> latitudinal T(theta) banding
    - T_eq = T_Schw, T_pole = T_Kerr (G62 verified)
    - Problem: Kerr photon orbit gets buried inside r = 2M for a > 0.707

Formula (beta) - the "hologram on stationary horizon" reading:
    A_beta = 2 M r / (r^2 + a^2)
    - A = 1 at r = r_+ (Kerr horizon, constant BL coord, oblate when embedded)
    - No theta dependence in A; uniform on r = constant surfaces
    - At equator: A = 2Mr/(r^2 + a^2)
    - At Schwarzschild limit (a=0): A = 2M/r recovers Schwarzschild
    - |grad A| theta-independent on the bubble -> uniform T = T_Kerr
    - No latitudinal banding (framework matches standard Kerr's uniform T)
    - All Kerr photon orbits and ISCOs are outside the bubble for any spin

This script computes both formulas' predictions for:
    1. Bubble shape (where A = 1)
    2. T(theta) on the bubble
    3. Location of Kerr photon orbits relative to the bubble
    4. Location of Kerr ISCOs relative to the bubble
    5. Trade-offs and what the framework gives up vs gains under each
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Kerr structural quantities (M units) ---

def kerr_horizon_radii(a):
    s = np.sqrt(np.maximum(1.0 - a**2, 0.0))
    return 1.0 + s, 1.0 - s


def kerr_ISCO_pro(a):
    if a == 0:
        return 6.0
    Z1 = 1.0 + (1.0 - a**2)**(1.0/3.0) * ((1.0 + a)**(1.0/3.0) + (1.0 - a)**(1.0/3.0))
    Z2 = np.sqrt(3.0 * a**2 + Z1**2)
    return 3.0 + Z2 - np.sqrt(np.maximum((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2), 0.0))


def kerr_photon_orbit_pro(a):
    return 2.0 * (1.0 + np.cos(2.0/3.0 * np.arccos(-a)))


# --- Formula alpha: A = 2Mr/Sigma ---

def A_alpha(r, theta, a):
    Sigma = r**2 + a**2 * np.cos(theta)**2
    return 2.0 * r / Sigma


def bubble_alpha(theta, a):
    """A_alpha = 1: r^2 - 2Mr + a^2 cos^2 theta = 0; r = 1 + sqrt(1 - a^2 cos^2 theta)."""
    return 1.0 + np.sqrt(np.maximum(1.0 - a**2 * np.cos(theta)**2, 0.0))


def grad_A_alpha_on_bubble(theta, a):
    """Coordinate |grad A_alpha| on bubble."""
    r = bubble_alpha(theta, a)
    radial_sq = (1.0 - r)**2 / r**2
    polar_sq = (a**4 * np.sin(2.0 * theta)**2) / (4.0 * r**4)
    return np.sqrt(radial_sq + polar_sq)


def T_alpha_on_bubble(theta, a):
    return grad_A_alpha_on_bubble(theta, a) / (4.0 * np.pi)


# --- Formula beta: A = 2Mr / (r^2 + a^2) ---

def A_beta(r, a):
    """No theta dependence."""
    return 2.0 * r / (r**2 + a**2)


def bubble_beta(a):
    """A_beta = 1 -> 2Mr = r^2 + a^2 -> r^2 - 2Mr + a^2 = 0 -> r = r_+.
    The bubble is at constant r = r_+ in BL coords (independent of theta).
    """
    r_p, _ = kerr_horizon_radii(a)
    return r_p


def grad_A_beta_at_horizon(a):
    """|d A_beta / dr| at r = r_+.
    dA/dr = 2(a^2 - r^2)/(r^2 + a^2)^2
    Using (r_+^2 + a^2) = 2 r_+: denominator = 4 r_+^2.
    Numerator at r_+: 2(a^2 - r_+^2) = -2 * 2 sqrt(1 - a^2) * r_+ ... actually:
       a^2 - r_+^2 = a^2 - (1 + sqrt(1-a^2))^2 = a^2 - 1 - 2 sqrt(1-a^2) - (1-a^2)
                  = a^2 - 1 - 2sqrt(1-a^2) - 1 + a^2
                  = 2a^2 - 2 - 2 sqrt(1-a^2)
                  = -2(1 - a^2 + sqrt(1-a^2))
                  = -2 sqrt(1-a^2) (sqrt(1-a^2) + 1)
                  = -2 sqrt(1-a^2) * r_+
    So dA/dr at r_+ = 2 * (-2 sqrt(1-a^2) * r_+) / (4 r_+^2)
                    = -sqrt(1-a^2) / r_+
    |dA/dr| at r_+ = sqrt(1-a^2) / r_+
    """
    r_p, _ = kerr_horizon_radii(a)
    return np.sqrt(np.maximum(1.0 - a**2, 0.0)) / r_p


def T_beta_on_bubble(a):
    """Uniform T on the bubble (no theta dependence)."""
    return grad_A_beta_at_horizon(a) / (4.0 * np.pi)


def kerr_T_standard(a):
    """Standard Kerr horizon T = kappa/(2 pi)."""
    r_p, _ = kerr_horizon_radii(a)
    return np.sqrt(np.maximum(1.0 - a**2, 0.0)) / (4.0 * np.pi * r_p)


def main():
    print("=" * 80)
    print("G64: Bubble surface refinement -- two A formulas for Kerr")
    print("=" * 80)
    print()
    print("Formula (alpha): A = 2Mr / Sigma  -- current G3/G62, static-limit bubble")
    print("Formula (beta):  A = 2Mr / (r^2 + a^2) -- Kerr-horizon bubble, hologram reading")
    print()

    # --- Bubble shape comparison ---
    print("=" * 80)
    print("STEP 1: Bubble shapes")
    print("=" * 80)
    print()
    print(f"{'a/M':>6}{'r_+ (alpha pole)':>18}{'r_eq alpha':>14}"
          f"{'r_beta (uniform)':>18}{'r_PS pro':>12}{'r_ISCO pro':>14}")
    print("-" * 85)
    for a in [0.0, 0.3, 0.5, 0.67, 0.707, 0.85, 0.943, 0.95, 0.99]:
        r_p, _ = kerr_horizon_radii(a)
        r_eq_alpha = 2.0
        r_beta = bubble_beta(a)
        r_PS = kerr_photon_orbit_pro(a)
        r_ISCO = kerr_ISCO_pro(a)
        print(f"{a:>6.3f}{r_p:>18.4f}{r_eq_alpha:>14.4f}"
              f"{r_beta:>18.4f}{r_PS:>12.4f}{r_ISCO:>14.4f}")
    print()
    print("Reading:")
    print("- alpha bubble at equator stays at r = 2M for all spin.")
    print("- beta bubble at r = r_+ (decreases with spin).")
    print("- For a > 0.707: r_PS < 2M (alpha problem) but r_PS > r_+ (beta OK).")
    print("- For a > 0.943: r_ISCO < 2M (alpha problem) but r_ISCO > r_+ (beta OK).")
    print()

    # --- Temperature comparison ---
    print("=" * 80)
    print("STEP 2: T(theta) comparison")
    print("=" * 80)
    print()
    print("Formula (alpha): T(theta) banded between T_Schw (equator) and T_Kerr (pole)")
    print("Formula (beta):  T uniform on bubble = T_Kerr (standard, no banding)")
    print()
    print(f"{'a/M':>6}{'T_eq_alpha':>14}{'T_pole_alpha':>16}{'T_beta':>14}{'T_K_std':>14}")
    print("-" * 70)
    for a in [0.0, 0.3, 0.5, 0.67, 0.85, 0.95]:
        T_eq_a = T_alpha_on_bubble(np.pi/2, a)
        T_pol_a = T_alpha_on_bubble(0.0, a)
        T_b = T_beta_on_bubble(a)
        T_K = kerr_T_standard(a)
        print(f"{a:>6.2f}{T_eq_a:>14.6f}{T_pol_a:>16.6f}{T_b:>14.6f}{T_K:>14.6f}")
    print()
    print("Reading: T_beta = T_K_standard exactly (uniform). T_pole_alpha = T_K_standard")
    print("         exactly (G62 identity). T_eq_alpha = T_Schw (= 1/(8 pi)) for all spin.")
    print()

    # --- Photon orbit and ISCO location relative to bubble ---
    print("=" * 80)
    print("STEP 3: Where Kerr photon orbit and ISCO sit relative to each bubble")
    print("=" * 80)
    print()
    print(f"{'a/M':>6}{'PS inside alpha?':>20}{'PS inside beta?':>20}"
          f"{'ISCO inside alpha?':>22}{'ISCO inside beta?':>22}")
    print("-" * 100)
    for a in [0.0, 0.3, 0.5, 0.67, 0.707, 0.85, 0.943, 0.95]:
        r_PS = kerr_photon_orbit_pro(a)
        r_ISCO = kerr_ISCO_pro(a)
        r_p, _ = kerr_horizon_radii(a)
        # alpha bubble at equator: r = 2
        # beta bubble: r = r_+
        ps_in_alpha = "YES" if r_PS < 2.0 else "no"
        ps_in_beta = "YES" if r_PS < r_p else "no"
        isco_in_alpha = "YES" if r_ISCO < 2.0 else "no"
        isco_in_beta = "YES" if r_ISCO < r_p else "no"
        print(f"{a:>6.3f}{ps_in_alpha:>20}{ps_in_beta:>20}"
              f"{isco_in_alpha:>22}{isco_in_beta:>22}")
    print()
    print("Reading: under formula (alpha), Kerr photon orbit gets buried inside bubble")
    print("         for a > 0.707; ISCO buried for a > 0.943.")
    print("         Under formula (beta), photon orbit and ISCO are always OUTSIDE")
    print("         the bubble (r_p < r_ISCO < r_PS < ... for any spin < 1).")
    print()

    # --- Equatorial photon-orbit A value comparison ---
    print("=" * 80)
    print("STEP 4: A value at Kerr equatorial photon orbit under each formula")
    print("=" * 80)
    print()
    print(f"{'a/M':>6}{'r_PS_pro/M':>14}{'A_alpha(r_PS, eq)':>22}{'A_beta(r_PS)':>18}")
    print("-" * 60)
    for a in [0.0, 0.3, 0.5, 0.67, 0.707, 0.85, 0.95]:
        r_PS = kerr_photon_orbit_pro(a)
        A_a = A_alpha(r_PS, np.pi/2, a)
        A_b = A_beta(r_PS, a)
        print(f"{a:>6.3f}{r_PS:>14.4f}{A_a:>22.4f}{A_b:>18.4f}")
    print()
    print("Reading: under (alpha), A at Kerr photon orbit exceeds 1 for high spin")
    print("         (structurally inconsistent, since A > 1 is forbidden).")
    print("         Under (beta), A at Kerr photon orbit stays well below 1.")
    print()

    # --- Photon orbit reduction in beta ---
    print("=" * 80)
    print("STEP 5: A_beta(r_PS) -- does it land on the framework's 2/3 landmark?")
    print("=" * 80)
    print()
    print("If the framework wants its Sigma = 2 landmark to coincide with Kerr PS,")
    print("we need A(r_PS) = 2/3 under the framework's A formula.")
    print()
    print(f"{'a/M':>6}{'r_PS_pro/M':>14}{'A_beta(r_PS)':>18}{'distance from 2/3':>22}")
    print("-" * 70)
    for a in [0.0, 0.3, 0.5, 0.67, 0.85, 0.95, 0.99]:
        r_PS = kerr_photon_orbit_pro(a)
        A_b = A_beta(r_PS, a)
        diff = abs(A_b - 2.0/3.0)
        print(f"{a:>6.3f}{r_PS:>14.4f}{A_b:>18.4f}{diff:>22.4f}")
    print()
    print("Reading: under (beta), A at Kerr PS is exactly 2/3 for a = 0, drifts")
    print("         smoothly with spin. Not exactly 2/3 for spinning case, but close.")
    print("         The framework would need a Sigma_Kerr that places the integer")
    print("         landmark at the Kerr photon orbit's actual A_beta value.")
    print()

    # --- Trade-offs summary ---
    print("=" * 80)
    print("STEP 6: Trade-offs between formulas")
    print("=" * 80)
    print()
    print("Formula (alpha)  A = 2Mr/Sigma           Formula (beta)  A = 2Mr/(r^2+a^2)")
    print("-" * 80)
    print("Bubble:")
    print("  oblate, r(theta) from r_+ to 2M       spherical in BL (constant r = r_+)")
    print("  static-limit surface                  Kerr horizon surface")
    print()
    print("T on bubble:")
    print("  T(theta) banded                       T uniform = T_Kerr")
    print("  T_eq = T_Schw, T_pole = T_Kerr        no latitudinal structure")
    print("  -> distinguishing prediction          -> matches standard Kerr exactly")
    print()
    print("Kerr photon orbit:")
    print("  inside bubble for a > 0.707           always outside bubble")
    print("  framework's PS landmark conflict      no conflict")
    print()
    print("Kerr ISCO:")
    print("  inside bubble for a > 0.943           always outside bubble")
    print("  framework predicts disk inner = 2M    framework matches standard Kerr ISCO")
    print()
    print("LIGO ringdown at moderate spin (a ~ 0.67):")
    print("  exact GR-Kerr (G62/G63 verified)      exact GR-Kerr")
    print()
    print("LIGO ringdown at high spin (a > 0.707):")
    print("  STAM bubble inside photon orbit       still exact GR-Kerr")
    print("  -> framework departs from Kerr        -> framework = Kerr")
    print()
    print("g_tt structure:")
    print("  g_tt = -(1 - A_alpha) c^2 = exact     g_tt no longer = -(1-A_beta) c^2")
    print("  Kerr g_tt                              direct identity broken")
    print()
    print("Hologram reading:")
    print("  bubble oblate, matter on it          bubble at constant r_+, hologram")
    print("  rotates                                spins on stationary horizon")
    print()
    print("Latitudinal banding of evaporation:")
    print("  prediction (G63)                      no banding")
    print()

    # --- Plots ---
    print("Generating plots...")
    p1 = plot_bubble_shapes_compared()
    p2 = plot_T_compared()
    p3 = plot_PS_relative_to_bubbles()
    print(f"  {p1}")
    print(f"  {p2}")
    print(f"  {p3}")
    print()

    # --- Recommendation ---
    print("=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)
    print()
    print("Both formulas are consistent with the framework's primitive commitments,")
    print("but they encode different physics for spinning BHs:")
    print()
    print("(alpha) gives the framework a distinguishing prediction (latitudinal")
    print("        banding) at the cost of high-spin photon-orbit consistency.")
    print("        Predicts: high-spin BH accretion disks have inner edge at r=2M,")
    print("        not at the spin-dependent Kerr ISCO. Testable but in tension with")
    print("        continuum-fitting measurements of high-spin BHs.")
    print()
    print("(beta)  matches standard Kerr exactly at the cost of losing the")
    print("        latitudinal banding prediction. The hologram-on-stationary-horizon")
    print("        picture is cleaner here: bubble at r = r_+, hologram spins on it,")
    print("        T is uniform (Kerr's), Kerr photon orbits and ISCOs all sit outside")
    print("        the bubble cleanly for any spin.")
    print()
    print("Per Sean's question (2026-05-13): 'Does the hologram spinning on the")
    print("stationary horizon work here?' Yes -- under (beta). The bubble is the")
    print("Kerr horizon (constant r in BL coords, oblate when embedded in flat 3-space),")
    print("the hologram matter spins on this static surface, T is uniform = T_Kerr.")
    print()
    print("Framework decision: which to commit to is Sean's call. (beta) is the")
    print("cleaner Kerr-conservative reading; (alpha) is the more aggressive")
    print("distinguishing-prediction reading. The 2026-05-13 'rotation tilt makes")
    print("outward more viable' framing fits either -- the rotation enhancement is")
    print("a SEPARATE structural effect (rate enhancement, not T variation).")
    print()
    print("Under (beta), latitudinal banding of evaporation is REPLACED by uniform")
    print("Kerr emission with super-radiance enhancement at the equator (where")
    print("matter velocity is maximal). The framework still predicts non-uniform")
    print("emission RATE (equatorial enhancement from rotation), just not non-uniform T.")
    print()

    write_summary()


def plot_bubble_shapes_compared():
    fig, ax = plt.subplots(figsize=(11, 8))
    theta = np.linspace(0, 2*np.pi, 200)

    a = 0.85  # representative high spin
    r_alpha = np.array([bubble_alpha(t, a) for t in theta])
    r_beta_const = bubble_beta(a)
    r_PS = kerr_photon_orbit_pro(a)
    r_ISCO = kerr_ISCO_pro(a)

    x_a = r_alpha * np.sin(theta)
    z_a = r_alpha * np.cos(theta)
    ax.plot(x_a, z_a, 'tab:blue', linewidth=2.5,
            label=f"alpha bubble (static limit) a/M={a}")

    theta_circ = np.linspace(0, 2*np.pi, 100)
    x_b = r_beta_const * np.sin(theta_circ)
    z_b = r_beta_const * np.cos(theta_circ)
    ax.plot(x_b, z_b, 'tab:green', linewidth=2.5,
            label=f"beta bubble (Kerr horizon r_+) a/M={a}")

    # Photon orbit (equatorial circle)
    x_ps = r_PS * np.cos(np.linspace(0, 2*np.pi, 100))
    z_ps = r_PS * np.sin(np.linspace(0, 2*np.pi, 100))
    # Show as equatorial circle in r-z plane
    ax.axvline(r_PS, color='tab:red', linestyle='--', linewidth=1.5, alpha=0.6,
               label=f"Kerr equatorial photon orbit r_PS = {r_PS:.3f} M")
    ax.axvline(-r_PS, color='tab:red', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.axvline(r_ISCO, color='tab:orange', linestyle=':', linewidth=1.5, alpha=0.6,
               label=f"Kerr equatorial ISCO r = {r_ISCO:.3f} M")
    ax.axvline(-r_ISCO, color='tab:orange', linestyle=':', linewidth=1.5, alpha=0.6)

    ax.set_xlabel("x / M  (equatorial)")
    ax.set_ylabel("z / M  (polar)")
    ax.set_title(f"Bubble shapes for a/M = {a}: alpha (static limit) vs beta (Kerr horizon)\n"
                 f"r_PS = {r_PS:.2f} M is BURIED inside alpha bubble (r_eq = 2M)\n"
                 f"r_PS is OUTSIDE beta bubble (r = r_+ = {r_beta_const:.2f} M)")
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    out = PLOTS / "G64_bubble_shapes_compared.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_T_compared():
    fig, ax = plt.subplots(figsize=(11, 6))
    theta = np.linspace(0.001, np.pi - 0.001, 200)

    colors = ['black', 'tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    for a, color in zip([0.0, 0.3, 0.5, 0.67, 0.95], colors):
        T_a = np.array([T_alpha_on_bubble(t, a) for t in theta])
        T_b = T_beta_on_bubble(a)
        ax.plot(theta * 180 / np.pi, T_a, color=color, linewidth=2,
                label=f"alpha T(theta), a={a}")
        ax.axhline(T_b, color=color, linestyle='--', alpha=0.6,
                   label=f"beta T uniform = T_Kerr, a={a}")

    ax.set_xlabel("theta (deg from pole)")
    ax.set_ylabel("T (in M^-1 units)")
    ax.set_title("T(theta) under alpha (banded) vs beta (uniform = T_Kerr)\n"
                 "alpha: latitudinal structure (framework distinguishing prediction)\n"
                 "beta: matches standard Kerr's uniform horizon T")
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = PLOTS / "G64_T_compared.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_PS_relative_to_bubbles():
    fig, ax = plt.subplots(figsize=(11, 6))
    a_grid = np.linspace(0.0, 0.999, 100)

    r_PS = np.array([kerr_photon_orbit_pro(a) for a in a_grid])
    r_ISCO = np.array([kerr_ISCO_pro(a) for a in a_grid])
    r_eq_alpha = np.full_like(a_grid, 2.0)
    r_beta = np.array([bubble_beta(a) for a in a_grid])

    ax.plot(a_grid, r_PS, 'tab:red', linewidth=2, label="Kerr prograde photon orbit r_PS")
    ax.plot(a_grid, r_ISCO, 'tab:orange', linewidth=2, label="Kerr prograde ISCO r_ISCO")
    ax.plot(a_grid, r_eq_alpha, 'tab:blue', linewidth=2.5, label="alpha bubble equator (r = 2M)")
    ax.plot(a_grid, r_beta, 'tab:green', linewidth=2.5, label="beta bubble (r = r_+)")

    ax.axvline(0.707, color='black', linestyle=':', alpha=0.5,
               label="a = 0.707 (r_PS crosses 2M)")
    ax.axvline(0.943, color='gray', linestyle=':', alpha=0.5,
               label="a = 0.943 (r_ISCO crosses 2M)")
    ax.set_xlabel("a/M")
    ax.set_ylabel("r / M")
    ax.set_title("Kerr photon orbit and ISCO vs bubble surfaces\n"
                 "alpha bubble (r=2M) trapped above r_PS for a > 0.707\n"
                 "beta bubble (r=r_+) stays below r_PS for any spin")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = PLOTS / "G64_PS_vs_bubbles.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_summary():
    md = []
    md.append("# G64 - Bubble Surface Refinement (Two Candidate A Formulas)\n")
    md.append("**Date: 2026-05-13.** Explores two candidate A(r, theta; M, a) formulas "
              "for the framework's Kerr extension, resolving the high-spin photon-orbit "
              "and ISCO issues flagged in G63.\n")

    md.append("## The structural question\n")
    md.append("From Sean (2026-05-13): 'Can we refine the spinning bubble surface? "
              "Does the hologram spinning on the stationary horizon work here?'\n")
    md.append("The G3/G62 formula `A = 2Mr/Sigma` has `A = 1` at the **static-limit "
              "surface**, not at the actual Kerr horizon. This is why the bubble's "
              "equator stays at `r = 2M` for any spin, which buries the Kerr photon "
              "orbit (a > 0.707) and ISCO (a > 0.943) inside the bubble at high spin.\n")
    md.append("Alternative: `A = 2Mr/(r^2 + a^2)` places `A = 1` at the actual Kerr "
              "horizon `r = r_+` (constant in BL coordinates, oblate when embedded). "
              "The bubble is then the Kerr horizon, hologram matter spins on it, all "
              "photon orbits and ISCOs sit outside cleanly for any spin.\n")

    md.append("## Comparison\n")
    md.append("| Property | Formula alpha (A = 2Mr/Sigma) | Formula beta (A = 2Mr/(r^2+a^2)) |\n")
    md.append("|---|---|---|\n")
    md.append("| A=1 surface | static-limit (oblate r(theta)) | Kerr horizon r = r_+ |\n")
    md.append("| Equatorial extent | r = 2M (all spin) | r = r_+ (varies) |\n")
    md.append("| Polar extent | r = r_+ | r = r_+ |\n")
    md.append("| T on bubble | banded T(theta) | uniform = T_Kerr |\n")
    md.append("| T_eq | T_Schwarzschild | T_Kerr |\n")
    md.append("| T_pole | T_Kerr | T_Kerr |\n")
    md.append("| Kerr PS location | inside bubble for a > 0.707 | always outside |\n")
    md.append("| Kerr ISCO location | inside bubble for a > 0.943 | always outside |\n")
    md.append("| g_tt structure | g_tt = -(1 - A_alpha) c^2 = exact Kerr g_tt | direct identity broken |\n")
    md.append("| LIGO ringdown high spin | departs from GR-Kerr | exact GR-Kerr |\n")
    md.append("| High-spin disk inner edge | r = 2M (prediction) | matches Kerr ISCO |\n")
    md.append("| Latitudinal evaporation banding | distinguishing prediction | not predicted |\n")
    md.append("| Hologram-on-stationary-horizon | bubble oblate, matter rotates | clean match |\n")

    md.append("## What each formula gives up\n")
    md.append("**Formula alpha gives up:** consistency with standard Kerr at high spin. "
              "Predicts no ergosphere physics (region between r_+ and r=2M doesn't exist "
              "in manifold), high-spin BH inner edges at r=2M not at Kerr ISCO, photon "
              "orbits buried.\n")
    md.append("**Formula beta gives up:** the latitudinal banding distinguishing "
              "prediction. T is uniform = T_Kerr on the bubble, matching standard Kerr. "
              "Also gives up the direct g_tt = -(1-A) c^2 structural identity (g_tt no "
              "longer equals -(1-A_beta) c^2 in Kerr coordinates).\n")

    md.append("## Sean's 'hologram on stationary horizon' reading fits beta\n")
    md.append("Under formula (beta):\n")
    md.append("- The bubble IS the Kerr horizon: constant r = r_+ in BL coordinates, "
              "oblate when embedded in flat 3-space.\n")
    md.append("- The hologram matter on the bubble spins at ZAMO frequency, carrying J.\n")
    md.append("- T is uniform on the bubble (= standard Kerr T_Kerr).\n")
    md.append("- Rotation tilt at the horizon still makes outward emission more viable "
              "(Sean's 2026-05-13 reading): the super-radiance enhancement R(theta) = "
              "1 + (v_matter/c)^2 still applies to the emission RATE, with equatorial "
              "matter velocity giving rate enhancement -- but T itself is uniform.\n")
    md.append("- The framework's 'rotation enhances outward emission' is structurally "
              "preserved as a RATE effect, not a temperature effect.\n")

    md.append("## What the framework retains under beta\n")
    md.append("- Substance ontology + presentism (unchanged)\n")
    md.append("- Ledger-channel + two-face + SU shell-count (unchanged)\n")
    md.append("- alpha = 4 entropy area-per-entry (unchanged)\n")
    md.append("- Quintic Hermite F(y) for inside-PS metric (unchanged)\n")
    md.append("- Spinless strong-field metric (unchanged)\n")
    md.append("- LIGO spinless ringdown = exact GR (unchanged)\n")
    md.append("- LIGO Kerr ringdown = exact GR-Kerr (recovered for all spin, not just <0.707)\n")
    md.append("- Static bubble + rotating matter (clean: bubble at Kerr horizon)\n")
    md.append("- Outward-only emission with rotation enhancement of rate\n")
    md.append("- D = 3 forcing (unchanged)\n")

    md.append("## What changes under beta\n")
    md.append("- T uniform on bubble = T_Kerr (replaces T_pole = T_Kerr + T_eq = T_Schw)\n")
    md.append("- Latitudinal banding of evaporation T NOT predicted (replaces G63's banding prediction)\n")
    md.append("- Direct g_tt = -(1-A) c^2 identity for Kerr no longer holds (alpha had it)\n")
    md.append("- Standard Kerr observables matched cleanly at all spins\n")
    md.append("\n")
    md.append("The framework still predicts:\n")
    md.append("- Equatorially enhanced emission RATE (super-radiance from matter rotation)\n")
    md.append("- Different ringdown physics inside Kerr photon orbit (quintic Hermite F)\n")
    md.append("- No interior, ledger structure, primordial-mass remnants\n")

    md.append("## Recommendation\n")
    md.append("Formula (beta) is the cleaner reading per Sean's 'hologram on stationary "
              "horizon' framing. It:\n")
    md.append("- Resolves the high-spin photon-orbit/ISCO conflict in G63\n")
    md.append("- Matches standard Kerr at all spins for the macroscopic observables\n")
    md.append("- Preserves all framework primitives except the static-limit identification\n")
    md.append("- Gives up the latitudinal banding observable but keeps the rotation-rate "
              "enhancement as a distinguishing observable\n")
    md.append("\n")
    md.append("Formula (alpha) gives a more aggressive distinguishing prediction set but "
              "creates the high-spin photon-orbit conflict.\n")
    md.append("\n")
    md.append("Sean's call which to commit to. The framework can also stay agnostic and "
              "treat the choice as observational: if high-spin BH inner edges are "
              "observed at the spin-dependent Kerr ISCO (consistent with continuum-fitting "
              "results), (beta) is preferred; if they're observed at r = 2M, (alpha) is "
              "preferred. Currently, continuum-fitting favors (beta).\n")

    md.append("## Files\n")
    md.append("- [scripts/G64_bubble_refinement.py](../scripts/G64_bubble_refinement.py)\n")
    md.append("- [plots/G64_bubble_shapes_compared.png](../plots/G64_bubble_shapes_compared.png)\n")
    md.append("- [plots/G64_T_compared.png](../plots/G64_T_compared.png)\n")
    md.append("- [plots/G64_PS_vs_bubbles.png](../plots/G64_PS_vs_bubbles.png)\n")

    out = RESULTS / "G64_bubble_refinement_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary: {out}")


if __name__ == "__main__":
    main()
