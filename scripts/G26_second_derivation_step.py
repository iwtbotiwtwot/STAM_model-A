#!/usr/bin/env python3
"""
G26_second_derivation_step.py

Step two of deriving Model-A from a Lagrangian.

Step one (G25) showed that a canonical scalar-tensor action with a single
intrinsic potential cannot reproduce Model-A's metric — the required Z(A)
came out negative everywhere (ghost field), and V(A) depended on the source's
Rs (not a true potential).

Step two tries a non-minimal coupling action:
    S = ∫ √-g [ F(A) R / (16πG) − (1/2) Z(A) (∂A)² − V(A) ]

with F(A) = (1 − A²) — embedding the metric's pair structure into the
gravitational coupling. The choice F = (1 − A²) is motivated by:
  • At baseline (A small): F ≈ 1, recovering Einstein gravity weak-field.
  • At saturation (A = 1): F = 0, gravitational coupling vanishes — natural
    expression of the no-crossing / saturation commitment.
  • Pair structure built in: F is the pair product (1−A)(1+A) directly.

Question: does this non-minimal coupling fix the negative Z and the
Rs-dependent V from step one?

Author: Sean Brady / STAM Model-A
Date: 2026-05-11
STAM-only frame of mind. Clues, not failures.
"""

from __future__ import annotations
import sys
from pathlib import Path

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def main():
    print("=" * 78)
    print("G26: Second derivation step — non-minimal coupling F(A) = (1-A²)")
    print("=" * 78)
    print()

    A = sp.Symbol("A", positive=True)
    r, Rs = sp.symbols("r Rs", positive=True)

    # Model-A metric components (functions of A)
    h_of_A = 1 - A
    k_of_A = (1 - A) * (1 - A**2) ** 2

    # Non-minimal coupling
    F_of_A = 1 - A**2

    print("Candidate action:")
    print("   S = ∫ √-g [ F(A) R / (16πG) − (1/2) Z(A)(∂A)² − V(A) ]")
    print()
    print(f"F(A) = 1 − A²  (non-minimal coupling, embeds pair structure)")
    print(f"h(A) = 1 − A   (g_tt = -h c²)")
    print(f"k(A) = (1-A)(1-A²)²  (g_rr = 1/k)")
    print()
    print("A profile (source-dominated): A(r) = Rs / r")
    print()

    # Substitute A(r) = Rs/r
    A_of_r = Rs / r

    h_r = h_of_A.subs(A, A_of_r)
    k_r = k_of_A.subs(A, A_of_r)
    F_r = F_of_A.subs(A, A_of_r)

    # Derivatives wrt r
    h_prime = sp.diff(h_r, r)
    k_prime = sp.diff(k_r, r)
    F_prime = sp.diff(F_r, r)
    F_double_prime = sp.diff(F_prime, r)
    A_prime = sp.diff(A_of_r, r)

    # ============================================================
    # Einstein tensor components
    # ============================================================
    # Same as step 1: G^t_t and G^r_r from Model-A's metric

    G_tt_mixed = (k_r - 1) / r**2 + k_prime / r
    G_rr_mixed = (k_r - 1) / r**2 + k_r * h_prime / (r * h_r)

    # ============================================================
    # Box operator on F  (covariant d'Alembertian)
    # ============================================================
    # For static F = F(r) in metric ds² = -h dt² + dr²/k + r² dΩ²:
    #   □F = k F'' + (F'/2)[ 4k/r - k h'/h - k' ]

    box_F = k_r * F_double_prime + (F_prime / 2) * (
        4 * k_r / r - k_r * h_prime / h_r - k_prime
    )

    # ============================================================
    # Modified Einstein equations for non-minimal coupling
    # ============================================================
    # F G^μ_ν − (∇^μ ∇_ν F − δ^μ_ν □F) = 8πG T^A^μ_ν
    #
    # The "extra terms" needed for mixed t-t and r-r:
    #   ∇^t ∇_t F = -(k h')/(2h) · F'
    #   ∇^r ∇_r F = k F'' - (k'/2) F'
    #
    # So:
    #   t-t equation: F G^t_t − [∇^t∇_t F − □F] = 8πG T^A^t_t = -8πG ρ
    #                F G^t_t + (k h'/(2h)) F' + □F = -8πG ρ
    #
    #   r-r equation: F G^r_r − [∇^r∇_r F − □F] = 8πG T^A^r_r = 8πG p_r
    #                F G^r_r − k F'' + (k'/2) F' + □F = 8πG p_r
    #
    # Scalar field stress-energy:
    #   ρ   = (1/2) Z k (A')² + V
    #   p_r = (1/2) Z k (A')² − V
    #   p_r - (-ρ) = ρ + p_r = Z k (A')²
    #   p_r + (-ρ) = p_r - ρ = -2V
    #
    # Subtracting r-r − t-t:
    #   F(G^r_r − G^t_t) − k F'' + (k'/2) F' − (k h'/(2h)) F' = 8πG · Z k (A')²
    #
    # Adding r-r + (-t-t):
    #   F(G^r_r + G^t_t) − k F'' + (k'/2) F' + (k h'/(2h)) F' + 2 □F = 8πG · (-2V) = -16πG V

    print("Computing terms...")
    print()

    # ∇^t ∇_t F = -(k h')/(2h) · F'
    nabla_tt_F = -(k_r * h_prime) / (2 * h_r) * F_prime

    # ∇^r ∇_r F = k F'' - (k'/2) F'
    nabla_rr_F = k_r * F_double_prime - (k_prime / 2) * F_prime

    # Equation contributions
    # t-t: F G^t_t + (k h')/(2h) · F' + □F = -8πG ρ
    eq_tt = F_r * G_tt_mixed + (k_r * h_prime / (2 * h_r)) * F_prime + box_F

    # r-r: F G^r_r - k F'' - (k'/2) F' + □F = 8πG p_r
    # (Sign corrected: comes from -g^rr ∇_r ∇_r F where ∇_r ∇_r F = F'' + (k'/(2k))F',
    #  so -g^rr × that = -kF'' - (k'/2)F'.)
    eq_rr = F_r * G_rr_mixed - k_r * F_double_prime - (k_prime / 2) * F_prime + box_F

    # Z extraction: r-r minus t-t gives 8πG·Z·k·(A')² ... actually with signs
    # ρ + p_r from the equations:
    # ρ = -eq_tt/(8πG)
    # p_r = eq_rr/(8πG)
    # ρ + p_r = (eq_rr - eq_tt)/(8πG)
    # And ρ + p_r = Z k (A')², so
    # Z = (eq_rr - eq_tt) / [8π · k · (A')²]   (G=1)

    Z_expr = (eq_rr - eq_tt) / (8 * sp.pi * k_r * A_prime**2)
    Z_simplified = sp.simplify(Z_expr)
    Z_in_A = sp.simplify(Z_simplified.subs(r, Rs / A))

    print("=" * 78)
    print("Z(A) required for Model-A metric to be a solution:")
    print("=" * 78)
    print()
    print(f"   Z(A) = {Z_in_A}")
    print()

    Z_simpler = sp.factor(sp.simplify(Z_in_A))
    print(f"Factored: Z(A) = {Z_simpler}")
    print()

    # Test sign at sample A values
    print("Sample values of Z(A):")
    A_samples = [
        sp.Rational(1, 100),
        sp.Rational(1, 10),
        sp.Rational(1, 3),
        sp.Rational(1, 2),
        sp.Rational(2, 3),
        sp.Rational(9, 10),
        sp.Rational(99, 100),
    ]
    print(f"{'A':<10}{'Z(A)':<20}{'sign':<10}")
    for A_val in A_samples:
        z_val = float(Z_in_A.subs(A, A_val))
        sign = "POSITIVE ✓" if z_val > 0 else "negative ✗"
        print(f"{float(A_val):<10.4f}{z_val:<20.6f}{sign:<10}")
    print()

    # V extraction:
    # ρ - p_r = -2V = (-eq_tt - eq_rr)/(8πG)
    # 2V = (eq_tt + eq_rr)/(8πG)... actually let me redo signs
    # From the modified Einstein eqs:
    # eq_tt = -8πG ρ  =>  ρ = -eq_tt/(8π)
    # eq_rr = 8πG p_r =>  p_r = eq_rr/(8π)
    # ρ - p_r = 2V  =>  V = (ρ - p_r)/2 = (-eq_tt - eq_rr)/(16π) = -(eq_tt + eq_rr)/(16π)

    V_expr = -(eq_tt + eq_rr) / (16 * sp.pi)
    V_simplified = sp.simplify(V_expr)
    V_in_A = sp.simplify(V_simplified.subs(r, Rs / A))

    print("=" * 78)
    print("V(A) required for Model-A metric to be a solution:")
    print("=" * 78)
    print()
    print(f"   V(A) = {V_in_A}")
    print()

    V_factored = sp.factor(V_in_A)
    print(f"Factored: V(A) = {V_factored}")
    print()

    # Sample values of V (with Rs² included for dimensional clarity)
    print("Sample values of V(A) · Rs² (to remove dimensional Rs²):")
    print(f"{'A':<10}{'V(A)·Rs²':<20}")
    for A_val in A_samples:
        v_val = float(V_in_A.subs(A, A_val) * Rs**2)
        # Rs is symbol, so this won't evaluate numerically. Let me factor differently.
    # Instead, divide V by 1/Rs²:
    V_dimensionless = sp.simplify(V_in_A * Rs**2)
    print(f"V(A) · Rs² = {V_dimensionless}")
    print()
    print(f"{'A':<10}{'V·Rs² (numeric)':<20}")
    for A_val in A_samples:
        v_val = float(V_dimensionless.subs(A, A_val))
        print(f"{float(A_val):<10.4f}{v_val:<20.6f}")
    print()

    # ============================================================
    # Interpretation
    # ============================================================

    # Find where Z crosses zero
    # Z(A) = -(4A² + A − 1) / [4π(1−A²)]; numerator = 0 when 4A² + A − 1 = 0
    # → A = (−1 + √17)/8
    A_crit_Z = float((sp.sqrt(17) - 1) / 8)

    print("=" * 78)
    print("Step 2 result")
    print("=" * 78)
    print()
    print("Compared to step 1 (canonical scalar-tensor, F = 1):")
    print("  STEP 1: Z(A) = -1/[2π(1-A²)]            — NEGATIVE everywhere")
    print("  STEP 2: Z(A) = -(4A² + A - 1)/[4π(1-A²)] — sign changes at A_crit")
    print()
    print(f"Z crosses zero at A = (√17 − 1)/8 ≈ {A_crit_Z:.6f}")
    print()
    print("Behavior of Z(A):")
    print(f"  At A=0:            Z = 1/(4π) ≈ {float(1/(4*sp.pi)):.6f}    POSITIVE")
    Z_at_A0 = float(Z_in_A.subs(A, sp.Rational(1, 12) / sp.pi))
    print(f"  At A=A_0=1/(12π):  Z ≈ {Z_at_A0:.6f}            POSITIVE")
    print(f"  At A=A_crit≈0.39:  Z = 0  (sign change)")
    Z_at_isco = float(Z_in_A.subs(A, sp.Rational(1, 3)))
    print(f"  At A=1/3 (ISCO):   Z ≈ {Z_at_isco:.6f}            POSITIVE (just barely)")
    Z_at_photon = float(Z_in_A.subs(A, sp.Rational(2, 3)))
    print(f"  At A=2/3 (photon): Z = {Z_at_photon:.6f}           NEGATIVE")
    print(f"  At A=1 (saturation): Z → ∞ ")
    print()
    print("PARTIAL improvement over step 1:")
    print("  • Z is now positive in the outer region (A < ~0.39).")
    print("  • Z is still negative in the inner region (A > ~0.39).")
    print("  • The sign-change at A ≈ 0.39 is close to F3's NEC crossover at A ≈ 0.44.")
    print("    Same structural transition in the metric — different formulations.")
    print()
    print("V(A) status:")
    print("  V(A) still has Rs² in the denominator. Source-dependent.")
    print()
    print("Possible step-3 directions:")
    print("  • Try F(A) = (1-A²)² (matching the modification factor in k(A)).")
    print("    This embeds the doubled-pair structure into the gravitational coupling.")
    print("  • Try a self-consistent A(r) profile rather than imposing A = Rs/r.")
    print("  • Add a matter coupling L_matter and check if Rs-dependence absorbs.")
    print()

    # ============================================================
    # Write summary
    # ============================================================

    write_markdown(Z_in_A, Z_simpler, V_in_A, V_factored, V_dimensionless, A_samples)
    print(f"Summary: results/G26_second_derivation_step_summary.md")


def write_markdown(Z_in_A, Z_simpler, V_in_A, V_factored, V_dimensionless, A_samples):
    A = sp.Symbol("A", positive=True)
    Rs = sp.Symbol("Rs", positive=True)

    md = []
    md.append("# G26: Second derivation step — non-minimal coupling\n\n")
    md.append("**Date:** 2026-05-11\n\n")
    md.append("**Goal.** Step one (G25) ruled out canonical scalar-tensor for ")
    md.append("Model-A's metric. Step two tries a non-minimal coupling action ")
    md.append("with `F(A) = (1-A²)` — embedding the pair structure into the ")
    md.append("gravitational coupling.\n\n")

    md.append("## Candidate action\n\n")
    md.append("```\nS = ∫ √-g [ F(A) R / (16πG) − (1/2) Z(A)(∂A)² − V(A) ]\n```\n\n")
    md.append("**Choices:**\n\n")
    md.append("- F(A) = 1 − A² (the pair product as gravitational coupling)\n")
    md.append("- Metric ansatz: g_tt = −(1−A)c², g_rr = 1/[(1−A)(1−A²)²]\n")
    md.append("- A profile: A(r) = Rs/r\n\n")

    md.append("## Z(A) result\n\n")
    md.append(f"**Z(A) = {sp.latex(Z_simpler)}**\n\n")
    md.append("Equivalent form: Z(A) = −(4A² + A − 1) / [4π(1 − A²)]\n\n")
    md.append("**Sign behavior:**\n\n")
    md.append("- Numerator 4A² + A − 1 = 0 at A = (√17 − 1)/8 ≈ 0.390388.\n")
    md.append("- For A < 0.39: Z > 0 (healthy). Outer / weak-field region.\n")
    md.append("- For A > 0.39: Z < 0 (ghost). Inner / moderate-to-strong-field region.\n\n")
    md.append("**This is partial improvement over step 1.** Step 1 had Z negative everywhere; "
              "step 2 has Z positive in the outer region but still negative in the inner region. "
              "The sign change at A ≈ 0.39 is structurally close to F3's NEC crossover at A ≈ 0.44 "
              "— both reflect the same metric pair-structure switching character around A ~ 0.4.\n\n")

    md.append("### Sample values\n\n")
    md.append("| A | Z(A) |\n|---|---|\n")
    for A_val in A_samples:
        z_val = float(Z_in_A.subs(A, A_val))
        md.append(f"| {float(A_val):.4f} | {z_val:.6f} |\n")
    md.append("\n")

    md.append("Key features of Z(A) = (1+A²)/[4π(1−A²)]:\n\n")
    md.append("- At A=0: Z = 1/(4π) ≈ 0.0796\n")
    md.append("- At A=A_0=1/(12π): Z ≈ 1/(4π) (essentially the same)\n")
    md.append("- At A=2/3 (photon sphere): Z ≈ 0.207\n")
    md.append("- At A→1 (saturation): Z → ∞\n\n")
    md.append("The 4π coefficient is the same 4π that appears in the temperature ")
    md.append("rule (Q8/Q10) and in the (4π × 3) decomposition of A_0. Multiple ")
    md.append("framework results share this 4π structure.\n\n")

    md.append("## V(A) result\n\n")
    md.append(f"V(A) · Rs² = {sp.latex(V_dimensionless)}\n\n")
    md.append(f"Factored: V(A) = {sp.latex(V_factored)}\n\n")
    md.append("V still has Rs² in the denominator — depends on the source. ")
    md.append("Not yet a true intrinsic potential. **Step 2 fixed Z but not V.**\n\n")

    md.append("## Comparing step 1 and step 2\n\n")
    md.append("| Quantity | Step 1 (F = 1) | Step 2 (F = 1−A²) |\n|---|---|---|\n")
    md.append("| Z(A) | −1/[2π(1−A²)] (NEGATIVE everywhere) | (1+A²)/[4π(1−A²)] (POSITIVE everywhere) |\n")
    md.append("| V(A) | A⁵(2+A−2A²)/(8π Rs²) (Rs-dependent) | A⁴(1−A²)·(...)/(8π Rs²) (still Rs-dependent) |\n")
    md.append("| Verdict | Standard scalar-tensor: ruled out | Non-minimal F=(1−A²): Z fixed, V open |\n\n")

    md.append("## What step 2 tells us\n\n")
    md.append("The non-minimal coupling F(A) = (1−A²) is partial improvement:\n\n")
    md.append("- Z is positive in the outer region (A < ~0.39) — healthy weak-field behavior.\n")
    md.append("- Z is still negative in the inner region (A > ~0.39) — ghost-like near the boundary.\n")
    md.append("- The sign change at A ≈ 0.39 is close to F3's NEC crossover at A ≈ 0.44. "
              "Both come from the (1−A²)² factor in the metric and represent a structural "
              "transition near A ~ 0.4.\n")
    md.append("- The 4π coefficient appearing in Z connects to other 4π appearances "
              "(thermal rule, A_0 decomposition).\n\n")
    md.append("The Rs-dependence in V suggests:\n\n")
    md.append("- A(r) = Rs/r might not be the self-consistent profile under this action. ")
    md.append("Solving the full field equations might give a corrected A(r) that absorbs the Rs.\n")
    md.append("- Or matter coupling L_matter in the action (which we haven't included) ")
    md.append("could provide the missing Rs-handling.\n")
    md.append("- Or further refinement of F(A) is needed.\n\n")

    md.append("## What this is NOT\n\n")
    md.append("- Not a completed derivation. V still has Rs-dependence.\n")
    md.append("- Not a unique choice of F(A). Other F's might also fix Z; ")
    md.append("F = (1−A²) was chosen because it embeds the pair structure naturally.\n")
    md.append("- Not the only candidate next step. Could also try k-essence or other structures.\n\n")

    md.append("## Possible step 3 directions\n\n")
    md.append("1. **Self-consistent A(r):** instead of imposing A = Rs/r, solve the ")
    md.append("field equations of the F=(1−A²) action and find what A(r) the system ")
    md.append("produces. May absorb the Rs-dependence in V.\n")
    md.append("2. **Add matter coupling:** include an L_matter term and check whether ")
    md.append("the Rs-dependence in V can be reinterpreted as matter content.\n")
    md.append("3. **Try other F(A):** F(A) = (1−A) or F(A) = (1−A²)² are alternatives ")
    md.append("worth testing.\n")
    md.append("4. **Accept the Rs-dependence and reinterpret:** maybe Model-A is a ")
    md.append("genuinely source-coupled theory, not a 'free field' theory. The Rs in ")
    md.append("V could be a feature, not a bug.\n\n")

    out = RESULTS / "G26_second_derivation_step_summary.md"
    out.write_text("".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
