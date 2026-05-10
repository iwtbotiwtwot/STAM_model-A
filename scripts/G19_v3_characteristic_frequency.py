#!/usr/bin/env python3
"""
G19_v3_characteristic_frequency.py

Compute the characteristic frequency of small A oscillations around the
V_3 potential minimum at A_0 = 1/(12π). This is the framework's "vacuum
tone" — the natural frequency at which the cosmic A field oscillates
around its structural-floor value.

Method:
    V_3(A) = α/A + β/(1-A), with α/β = [A_0/(1-A_0)]² fixing the
    minimum at A_0.

    For homogeneous A field with kinetic Lagrangian (1/2) M_P²(∂A)²
    (the standard scalar-field normalization with dimensionless A and
    Planck-mass kinetic prefactor), the equation of motion around
    the minimum is:

        d²(δA)/dt² = -[V''(A_0)/M_P²] × δA

    giving angular frequency ω² = V''(A_0)/M_P²_red.

    V''(A_0) computed analytically:
        V''(A_0) = 2α/A_0³ + 2β/(1-A_0)³
        With α/β = [A_0/(1-A_0)]²:
        V''(A_0) = 2β/[A_0 (1-A_0)³]

What we ask:
    What is ω numerically, and how does it compare to natural cosmic
    scales (Hubble rate H_0, Planck frequency, etc.)?

    If ω ~ H_0: the vacuum oscillation is on cosmic timescales —
        cosmologically relevant, slow.
    If ω ~ Planck: the vacuum oscillation is at quantum-gravity rates.
    If ω ~ something else: the framework picks out a specific scale.

This is an attempt to find harmonic structure in V_3 — the question
posed by the symphony reading: does V_3 produce a natural tone, and
what frequency does it sit at?
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

# Physical constants (SI)
G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
K_B = 1.380649e-23

# Cosmology at H_0 = 73 (SH0ES)
H_0_KMS_MPC = 73.04
MPC_TO_M = 3.085677581e22
H_0_SI = H_0_KMS_MPC * 1000.0 / MPC_TO_M  # in 1/s

# Critical density (SI)
RHO_CRIT_SI = 3.0 * H_0_SI**2 / (8.0 * np.pi * G)  # kg/m³
RHO_CRIT_ENERGY_SI = RHO_CRIT_SI * C**2  # J/m³

# Planck units
M_P = np.sqrt(HBAR * C / G)  # full Planck mass, kg
M_P_RED = np.sqrt(HBAR * C / (8.0 * np.pi * G))  # reduced Planck mass, kg
L_PLANCK = np.sqrt(HBAR * G / C**3)
T_PLANCK = L_PLANCK / C
OMEGA_PLANCK = 1.0 / T_PLANCK


# --- Symbolic V_3 derivatives ---

def symbolic_v3_double_prime():
    """Compute V_3''(A_0) symbolically with α/β fixed by minimum at A_0."""
    A, A_0_sym, beta = sp.symbols("A A_0 beta", positive=True)
    alpha = beta * (A_0_sym / (1 - A_0_sym))**2
    V3 = alpha/A + beta/(1 - A)
    V3_pp = sp.diff(V3, A, 2)
    V3_pp_at_A0 = sp.simplify(V3_pp.subs(A, A_0_sym))
    return V3, V3_pp_at_A0


# --- β calibrations ---

# β/ρ_crit values from prior framework work:
# Tracking ansatz (G11): β/ρ_crit = 0.6491
# Numerical KG (G8c quarantined): β/ρ_crit = 0.4265
# These differ because tracking and full dynamics give different β; we
# report both for transparency.
BETA_TILDE_TRACKING = 0.6491
BETA_TILDE_NUMERICAL_KG = 0.4265


def v_double_prime_at_A0(beta_tilde: float) -> float:
    """V''(A_0) in SI units, with β = β_tilde × ρ_crit_energy."""
    A_0 = 1.0 / (12.0 * np.pi)
    factor = 2.0 / (A_0 * (1.0 - A_0)**3)
    beta = beta_tilde * RHO_CRIT_ENERGY_SI
    return factor * beta


def omega_v3(beta_tilde: float) -> dict:
    """Compute V_3 characteristic frequency for given β calibration.

    Standard scalar-field Lagrangian with dimensionless A:
        L = (1/2) M_P²_red (∂A)² - V(A)
    gives ω² = V''(A_0) / M_P²_red.

    But M_P²_red here is in units that make the Lagrangian have units
    of energy density. The proper natural-units treatment: ω² has units
    of 1/time², V''(A_0) has units of energy density, M_P_red²_natural
    needs to be V''/ω² (i.e., energy density × time²).

    Concretely: M_P²_red (in natural units where ℏ=c=1) = 1/(8πG).
    With ℏ and c restored: M_P²_red carries units of energy × time² /
    length³ = (J·s²/m³). The ratio V''(A_0)[J/m³] / M_P²_red[J·s²/m³]
    = 1/s² as required.

    M_P²_red (with units J·s²/m³) = ℏc/(8πG) × (1/c⁴) [accounting for
    the kinetic-term normalization]... actually the cleanest derivation
    is via natural-units:

        ω² = V''(A_0) / M_P²_red    [in natural units]

    where both sides are in (energy)² and we then convert back to 1/s².
    """
    A_0 = 1.0 / (12.0 * np.pi)

    # Compute V''(A_0) in J/m³ (energy density)
    V_pp_J_per_m3 = v_double_prime_at_A0(beta_tilde)

    # Convert to natural units (eV⁴).
    # 1 J/m³ = 1 J × (1/m)³.
    # 1 J in eV: 1 J = 6.242×10¹⁸ eV
    # 1 m in eV⁻¹ (natural units): ℏc = 1.973×10⁻⁷ eV·m → 1 m = 1/(1.973×10⁻⁷) eV⁻¹ = 5.068×10⁶ eV⁻¹
    # So 1/m = 1/(5.068×10⁶ eV⁻¹) eV = 1.973×10⁻⁷ eV
    # And 1 m³ = (5.068×10⁶)³ eV⁻³
    # 1 J/m³ = (6.242×10¹⁸ eV) × (1.973×10⁻⁷ eV)³ = 6.242×10¹⁸ × 7.685×10⁻²¹ eV⁴
    #        = 4.798×10⁻² eV⁴ ... actually let me just use a single factor.
    # Numerical: 1 J/m³ = 6.241509e18 eV / (5.068e6 eV⁻¹)³ = 6.241509e18 / 1.302e20 eV⁴
    #          ≈ 4.79e-2 eV⁴
    EV_PER_J = 6.241509e18
    M_PER_INV_EV = 1.973269788e-7  # ℏc in eV·m (so 1/m = 1.973e-7 eV)
    J_PER_M3_TO_EV4 = EV_PER_J * (M_PER_INV_EV)**3
    V_pp_eV4 = V_pp_J_per_m3 * J_PER_M3_TO_EV4

    # Reduced Planck mass in eV
    M_P_RED_EV = 4.34106e-9 * EV_PER_J / (1.78266e-36 * EV_PER_J)  # this is wrong
    # Let me just use M_P_red = 2.435e27 eV (standard value)
    M_P_RED_EV = 2.435323e27

    # ω² in eV²
    omega_sq_eV2 = V_pp_eV4 / M_P_RED_EV**2
    omega_eV = np.sqrt(omega_sq_eV2)

    # Convert to 1/s
    # ℏ in eV·s = 6.582119e-16
    HBAR_EV_S = 6.582119e-16
    omega_per_s = omega_eV / HBAR_EV_S

    return {
        "beta_tilde": beta_tilde,
        "V_pp_J_per_m3": V_pp_J_per_m3,
        "V_pp_eV4": V_pp_eV4,
        "omega_eV": omega_eV,
        "omega_per_s": omega_per_s,
        "omega_over_H0": omega_per_s / H_0_SI,
        "period_s": 2.0 * np.pi / omega_per_s,
        "period_Gyr": 2.0 * np.pi / omega_per_s / (1e9 * 365.25 * 24 * 3600),
    }


def omega_v3_sanity_check_via_natural_units(beta_tilde: float) -> float:
    """Sanity check: compute ω directly in natural units to verify."""
    A_0 = 1.0 / (12.0 * np.pi)
    factor = 2.0 / (A_0 * (1.0 - A_0)**3)

    # Hubble in natural units (eV)
    HBAR_EV_S = 6.582119e-16
    H0_eV = H_0_SI * HBAR_EV_S

    # ρ_crit in eV⁴: ρ_crit = 3 H² M_P²_red (in natural units)
    M_P_RED_EV = 2.435323e27
    rho_crit_eV4 = 3.0 * H0_eV**2 * M_P_RED_EV**2

    # V''(A_0) in eV⁴
    V_pp_eV4 = factor * beta_tilde * rho_crit_eV4

    # ω in eV
    omega_eV = np.sqrt(V_pp_eV4 / M_P_RED_EV**2)

    # Convert to 1/s
    omega_per_s = omega_eV / HBAR_EV_S

    return omega_per_s


def main() -> None:
    print("G19: V_3 characteristic frequency at the structural-floor minimum")
    print("=" * 78)
    print()

    # --- Symbolic part ---
    V3_sym, V_pp_at_A0_sym = symbolic_v3_double_prime()
    print("Symbolic V_3:")
    print(f"  V_3(A) = α/A + β/(1-A),  α = β [A_0/(1-A_0)]²")
    print()
    print(f"V_3''(A_0) = 2β/[A_0(1-A_0)³]")
    print()

    A_0 = 1.0 / (12.0 * np.pi)
    factor = 2.0 / (A_0 * (1.0 - A_0)**3)
    print(f"  A_0 = 1/(12π) = {A_0:.6f}")
    print(f"  (1-A_0) = {1-A_0:.6f}")
    print(f"  (1-A_0)³ = {(1-A_0)**3:.6f}")
    print(f"  factor 2/[A_0(1-A_0)³] = {factor:.4f}")
    print(f"  V_3''(A_0) = {factor:.4f} × β")
    print()

    # --- Cosmological setup ---
    print("Cosmological setup (SI):")
    print(f"  H_0 = {H_0_KMS_MPC} km/s/Mpc = {H_0_SI:.4e} 1/s")
    print(f"  ρ_crit = {RHO_CRIT_SI:.4e} kg/m³")
    print(f"  ρ_crit × c² = {RHO_CRIT_ENERGY_SI:.4e} J/m³")
    print(f"  Reduced Planck mass M_P_red = {M_P_RED:.4e} kg")
    print(f"  Planck frequency ω_P = {OMEGA_PLANCK:.4e} 1/s")
    print()

    # --- Compute ω for both calibrations ---
    print("=" * 78)
    print("V_3 CHARACTERISTIC FREQUENCY")
    print("=" * 78)
    print()
    print("Two β calibrations available from prior framework work:")
    print()

    for label, bt in [("Tracking ansatz (G11)", BETA_TILDE_TRACKING),
                       ("Numerical KG (G8c)", BETA_TILDE_NUMERICAL_KG)]:
        result = omega_v3(bt)
        # Sanity check
        omega_check = omega_v3_sanity_check_via_natural_units(bt)

        print(f"  {label}:")
        print(f"    β/ρ_crit (β_tilde):       {bt:.4f}")
        print(f"    V''(A_0) [eV⁴]:           {result['V_pp_eV4']:.4e}")
        print(f"    ω [eV]:                   {result['omega_eV']:.4e}")
        print(f"    ω [1/s]:                  {result['omega_per_s']:.4e}")
        print(f"    ω / H_0:                  {result['omega_over_H0']:.4f}")
        print(f"    Period 2π/ω [s]:          {result['period_s']:.4e}")
        print(f"    Period [Gyr]:             {result['period_Gyr']:.4f}")
        print(f"    Sanity check ω [1/s]:     {omega_check:.4e}")
        print()

    # --- Comparisons ---
    print("=" * 78)
    print("WHAT THE NUMBER MEANS — COMPARISONS")
    print("=" * 78)
    print()

    omega_avg = 0.5 * (omega_v3(BETA_TILDE_TRACKING)["omega_per_s"] +
                       omega_v3(BETA_TILDE_NUMERICAL_KG)["omega_per_s"])

    age_universe_s = 13.8e9 * 365.25 * 24 * 3600

    print(f"Hubble rate (today):         H_0       = {H_0_SI:.4e} 1/s")
    print(f"V_3 oscillator frequency:    ω         = {omega_avg:.4e} 1/s  "
          f"(avg of calibrations)")
    print(f"Planck frequency:            ω_P       = {OMEGA_PLANCK:.4e} 1/s")
    print(f"Age of universe:             1/H_0     ≈ {1/H_0_SI:.4e} s "
          f"({1/H_0_SI / (1e9 * 365.25 * 24 * 3600):.2f} Gyr)")
    print()
    print(f"ω / H_0:                              ≈ "
          f"{omega_avg / H_0_SI:.2f}")
    print(f"ω / ω_Planck:                         ≈ "
          f"{omega_avg / OMEGA_PLANCK:.4e}")
    print(f"V_3 vacuum oscillation period:        ≈ "
          f"{2*np.pi/omega_avg / (1e9 * 365.25 * 24 * 3600):.2f} Gyr")
    print(f"  (compare to age of universe:        ≈ "
          f"{age_universe_s / (1e9 * 365.25 * 24 * 3600):.2f} Gyr)")
    print()

    # --- Verdict ---
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print()

    omega_to_H0 = omega_avg / H_0_SI
    if 1 < omega_to_H0 < 100:
        print(f"V_3's characteristic frequency sits at ω ≈ {omega_to_H0:.1f} × H_0.")
        print()
        print("This is a COSMOLOGICALLY-RELEVANT scale, not a Planck or atomic")
        print("scale. The cosmic A field's natural oscillation around its")
        print("structural-floor value happens on timescales comparable to the")
        print("age of the universe.")
        print()
        print("Specifically:")
        print(f"  - Vacuum oscillation period: ~{2*np.pi/omega_avg / (1e9 * 365.25 * 24 * 3600):.1f} Gyr")
        print(f"  - Age of universe (1/H_0):    ~{1/H_0_SI / (1e9 * 365.25 * 24 * 3600):.1f} Gyr")
        print(f"  - Ratio:                      ~{(2*np.pi/omega_avg) * H_0_SI:.2f} × Hubble time")
        print()
        print("The framework picks out a frequency that's in the same regime as")
        print("Hubble expansion — this is the cosmic A field's 'tone.' If V_3")
        print("represents the spacetime-structure potential, its tone is at the")
        print("scale of cosmic dynamics, not microscopic dynamics.")
    elif omega_to_H0 < 1:
        print(f"ω is below H_0. Vacuum oscillation is slower than Hubble expansion;")
        print("the field is effectively static on cosmic timescales.")
    elif omega_to_H0 > 100:
        print(f"ω is much larger than H_0. Vacuum oscillation is fast compared to")
        print("cosmic expansion; would average out on cosmic timescales.")
    print()

    print("=" * 78)
    print("HONEST ASSESSMENT")
    print("=" * 78)
    print()
    print("What this calculation does:")
    print("  1. Confirms V_3 has a definite harmonic-oscillator structure")
    print("     around its minimum at A_0 = 1/(12π).")
    print("  2. The characteristic frequency comes out to ω ≈ 10 × H_0,")
    print("     placing it in the cosmic-dynamics regime.")
    print("  3. This is the framework's first 'natural tone' — a specific")
    print("     frequency the structure picks out, computed from V_3's")
    print("     curvature at A_0 with the standard scalar-field kinetic")
    print("     normalization.")
    print()
    print("What this does NOT do:")
    print("  1. Prove this frequency is 'harmonic' in the symphony sense.")
    print("     One frequency alone doesn't make harmony; we'd need a")
    print("     spectrum of modes with structured ratios.")
    print("  2. Tell us whether the specific number ~10×H_0 is structural")
    print("     or just a consequence of the empirical Ω_DE calibration of")
    print("     β. If β were calibrated differently, ω would shift.")
    print("  3. Connect to the strong-field thirds (1/3, 2/3, 1) directly.")
    print("     Those orbital values come from a different mechanism")
    print("     (orbital mechanics in (1-A) potential), not from V_3's")
    print("     curvature.")
    print()
    print("What it suggests:")
    print("  V_3 produces a natural oscillation frequency at the cosmic-")
    print("  dynamics scale, not the Planck scale. If we extend this to a")
    print("  full mode analysis (spectrum of A oscillations bounded between")
    print("  A_0 and A=1, instead of just the curvature at A_0), we'd get")
    print("  a discrete spectrum. Whether those modes match the orbital")
    print("  thirds is the harmony question that this single-number")
    print("  calculation can't answer alone, but the calculation does")
    print("  establish: there IS a natural frequency in V_3, and it sits")
    print("  at a meaningful cosmic scale.")

    write_summary(omega_avg, omega_to_H0)


def write_summary(omega_avg, omega_to_H0):
    md = []
    md.append("# G19: V_3 Characteristic Frequency at the Structural-Floor Minimum\n")

    md.append("## Setup\n")
    md.append(
        "Compute the characteristic frequency of small A oscillations around "
        "V_3's minimum at A_0 = 1/(12π). The framework's vacuum-tone calculation.\n"
        "\n"
        "**V_3 form:** `V_3(A) = α/A + β/(1-A)` with `α/β = [A_0/(1-A_0)]²` "
        "fixing the minimum at A_0.\n"
        "\n"
        "**Method:** standard scalar-field normalization with dimensionless A:\n"
        "```text\n"
        "L = (1/2) M_P²_red (∂A)² − V(A)\n"
        "ω² = V''(A_0) / M_P²_red\n"
        "V''(A_0) = 2β / [A_0 (1 − A_0)³]\n"
        "```\n"
    )

    md.append("## Result\n")

    track = omega_v3(BETA_TILDE_TRACKING)
    num_kg = omega_v3(BETA_TILDE_NUMERICAL_KG)

    md.append(f"| β calibration | β/ρ_crit | ω (1/s) | ω/H_0 | Period (Gyr) |\n")
    md.append("|---|---:|---:|---:|---:|\n")
    md.append(f"| Tracking ansatz (G11) | {BETA_TILDE_TRACKING:.4f} | "
              f"{track['omega_per_s']:.4e} | {track['omega_over_H0']:.2f} | "
              f"{track['period_Gyr']:.2f} |\n")
    md.append(f"| Numerical KG (G8c) | {BETA_TILDE_NUMERICAL_KG:.4f} | "
              f"{num_kg['omega_per_s']:.4e} | {num_kg['omega_over_H0']:.2f} | "
              f"{num_kg['period_Gyr']:.2f} |\n")
    md.append("\n")

    md.append(f"**Average ω ≈ {omega_avg:.4e} 1/s ≈ {omega_to_H0:.1f} × H_0.**\n\n")

    md.append("## What this means\n")
    md.append(
        f"V_3's structural-floor minimum has a definite harmonic-oscillator "
        f"frequency. The frequency sits at **ω ≈ {omega_to_H0:.1f} × H_0** — "
        f"approximately ten Hubble rates per cycle, or a vacuum-oscillation "
        f"period of ~7-8 Gyr (compared to age of universe ~13.8 Gyr).\n"
        f"\n"
        f"The framework picks out a **cosmic-dynamics scale** for V_3's natural "
        f"tone, not a Planck or atomic scale. This is meaningful: if V_3 "
        f"represents the spacetime-structure potential, the fact that its "
        f"natural frequency falls at cosmic-dynamics scales (rather than "
        f"quantum-gravity or particle-physics scales) ties the structural floor "
        f"directly to cosmic-evolution physics rather than to microscopic physics.\n"
    )

    md.append("## What this is and isn't\n")
    md.append(
        "**Is:**\n"
        "- A clean computation of V_3's curvature at A_0 and the resulting "
        "characteristic frequency.\n"
        "- A confirmation that V_3 has a natural tone; the structure isn't flat "
        "near its minimum.\n"
        "- A first quantitative answer to 'what does the symphony sound like' "
        "at the level of fundamental frequency.\n"
        "\n"
        "**Isn't:**\n"
        "- A derivation of harmonic structure. One frequency doesn't make "
        "harmony; we'd need a spectrum of modes with structured ratios. The "
        "next step would be solving a wave equation for A bounded between "
        "A_0 and A=1, computing all eigenmodes, and seeing if their ratios "
        "match the orbital thirds (1/3, 2/3, 1) or other framework features.\n"
        "- Independent of the empirical Ω_DE calibration. β's value is set "
        "by matching observed Ω_DE ≈ 0.685; if β changes, ω scales accordingly. "
        "The specific number ~10×H_0 inherits this calibration.\n"
        "- A connection to the strong-field thirds. Those values come from "
        "orbital mechanics in (1-A) potential, not from V_3's curvature. "
        "Whether they're related is open.\n"
    )

    md.append("## What it suggests\n")
    md.append(
        "V_3 has a definite vacuum tone at cosmic-dynamics scales. The "
        "framework picks out a specific cosmologically-relevant frequency, "
        "rather than something Planck-scale or otherwise extreme. This is "
        "structurally meaningful — it connects V_3's structural-floor "
        "commitment directly to cosmic evolution physics rather than to "
        "microscopic quantum-gravity physics.\n"
        "\n"
        "Whether this constitutes 'harmonic structure' in the symphony sense "
        "depends on whether a full mode spectrum (not just the fundamental) "
        "shows structured ratios. That's the next-step calculation: solve a "
        "wave equation for A bounded between A_0 and A=1, see what eigenmode "
        "spectrum emerges. If the modes show integer ratios or match the "
        "orbital thirds, harmony is real. If not, V_3 has a single natural "
        "tone but no scale-spanning harmonic structure.\n"
    )

    out = RESULTS / "G19_v3_characteristic_frequency_summary.md"
    out.write_text("".join(md), encoding="utf-8")
    print(f"\nSummary: {out}")


if __name__ == "__main__":
    main()
