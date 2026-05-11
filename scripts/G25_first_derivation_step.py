#!/usr/bin/env python3
"""
G25_first_derivation_step.py

Step one toward deriving Model-A from a Lagrangian. STAM-only frame of mind:
we ask what action structure makes Model-A's metric commitments consistent.

Approach: assume a candidate scalar-field action of the form
    S = ∫ √-g [ R/(16πG) − (1/2) Z(A) g^μν ∂_μ A ∂_ν A − V(A) ]

Substitute Model-A's metric ansatz:
    g_tt = −(1−A) c²
    g_rr = 1 / [(1−A)(1−A²)²]

Substitute the source-dominated profile A(r) = Rs/r (the form A takes near a
point mass, where the source dominates over the cosmic baseline A_0).

Compute the Einstein tensor symbolically and ask: what Z(A) and V(A) does
the candidate action need to be consistent with Model-A's metric?

We are NOT asking "is Model-A consistent with GR." We are asking what STAM's
field-theoretic structure has to be, given its commitments. Whatever the
math returns is information about STAM — clue, not failure.

Author: Sean Brady / STAM Model-A
Date: 2026-05-11
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


# ============================================================
# Setup
# ============================================================

def main():
    print("=" * 78)
    print("G25: First derivation step — STAM action structure from metric")
    print("=" * 78)
    print()

    A = sp.Symbol("A", positive=True)
    r, Rs = sp.symbols("r Rs", positive=True)

    # Model-A's metric, as functions of A
    h_of_A = 1 - A
    k_of_A = (1 - A) * (1 - A ** 2) ** 2  # = (1-A)^3 (1+A)^2

    print("Metric ansatz (Model-A's commitment):")
    print(f"   g_tt = -h(A) c²,  h(A) = 1 - A")
    print(f"   g_rr = 1 / k(A),  k(A) = (1-A)(1-A²)² = (1-A)³(1+A)²")
    print()

    # Source-dominated profile
    A_of_r = Rs / r
    print(f"A profile (source-dominated): A(r) = Rs / r")
    print()

    # Substitute into h and k
    h_r = h_of_A.subs(A, A_of_r)
    k_r = k_of_A.subs(A, A_of_r)

    # Derivatives wrt r
    h_prime_r = sp.diff(h_r, r)
    k_prime_r = sp.diff(k_r, r)

    # ------------------------------------------------------------
    # Einstein tensor mixed components for static spherically symmetric
    # ds² = -h(r) dt² + dr²/k(r) + r² dΩ²
    #
    # Standard textbook formulas (c = G = 1):
    #   G^t_t = (k − 1)/r² + k'(r)/r
    #   G^r_r = (k − 1)/r² + k(r) h'(r) / (r h(r))
    # ------------------------------------------------------------

    G_tt_mixed = (k_r - 1) / r**2 + k_prime_r / r
    G_rr_mixed = (k_r - 1) / r**2 + k_r * h_prime_r / (r * h_r)

    # Simplify and convert to functions of A (substitute r = Rs/A)
    G_tt_simpl = sp.simplify(G_tt_mixed)
    G_rr_simpl = sp.simplify(G_rr_mixed)

    print("Einstein tensor (mixed indices) for Model-A's metric with A=Rs/r:")
    print(f"   G^t_t = {sp.simplify(G_tt_simpl)}")
    print()
    print(f"   G^r_r = {sp.simplify(G_rr_simpl)}")
    print()

    # Express in terms of A by substituting Rs = r·A:
    # r = Rs/A, so r² = Rs²/A²
    G_tt_in_A = sp.simplify(G_tt_simpl.subs(r, Rs / A))
    G_rr_in_A = sp.simplify(G_rr_simpl.subs(r, Rs / A))

    print("Same, written in terms of A (using r = Rs/A):")
    print(f"   G^t_t = {G_tt_in_A}")
    print()
    print(f"   G^r_r = {G_rr_in_A}")
    print()

    # ------------------------------------------------------------
    # Scalar field stress-energy
    #
    # For action S = ∫ √-g [R/(16πG) − (1/2) Z(A)(∂A)² − V(A)]:
    #   ρ_A   = T^0_0 with sign     = (1/2) Z k (A')² + V
    #   p_r,A = T^r_r                = (1/2) Z k (A')² − V
    #
    # Einstein eqs: G^μ_ν = 8πG T^μ_ν
    #   G^t_t = −8πG ρ_A      →  ρ_A = −G^t_t / (8π)
    #   G^r_r =  8πG p_r,A    →  p_r = G^r_r / (8π)
    #
    # From these:
    #   ρ + p_r  = Z k (A')²     ⇒  Z(A) = (G^r_r − G^t_t) / [8π · k(A) · (A'(r))²]
    #   ρ − p_r  = 2 V(A)        ⇒  V(A) = −(G^t_t + G^r_r) / (16π)
    # ------------------------------------------------------------

    # A'(r) = dA/dr = -Rs/r² = -A²/Rs  (using A=Rs/r)
    A_prime_r = sp.diff(A_of_r, r)
    A_prime_in_A = A_prime_r.subs(r, Rs / A)  # = -A²/Rs

    # k(A) in symbolic A
    k_A_sym = (1 - A) * (1 - A**2) ** 2

    # Compute Z(A)
    Z_numerator = G_rr_in_A - G_tt_in_A
    Z_denominator = 8 * sp.pi * k_A_sym * A_prime_in_A ** 2
    Z_of_A = sp.simplify(Z_numerator / Z_denominator)

    print("=" * 78)
    print("REQUIRED Z(A) for Model-A's metric to satisfy the action:")
    print("=" * 78)
    print()
    print(f"   Z(A) = {Z_of_A}")
    print()

    # Compute V(A)
    V_of_A = sp.simplify(-(G_tt_in_A + G_rr_in_A) / (16 * sp.pi))

    print("=" * 78)
    print("REQUIRED V(A) for Model-A's metric to satisfy the action:")
    print("=" * 78)
    print()
    print(f"   V(A) = {V_of_A}")
    print()

    # ------------------------------------------------------------
    # Numerical samples and interpretation
    # ------------------------------------------------------------

    print("=" * 78)
    print("Numerical samples")
    print("=" * 78)
    print()

    A_samples = [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 3),
                 sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(9, 10),
                 sp.Rational(99, 100)]
    print(f"{'A':<12}{'Z(A)':<25}{'V(A) · Rs²':<25}")
    for A_val in A_samples:
        Z_num = float(Z_of_A.subs(A, A_val))
        V_num = float(V_of_A.subs(A, A_val) * Rs ** 2)
        # V_num is V × Rs², dimensionless coefficient
        print(f"{float(A_val):<12.4f}{Z_num:<25.6f}{V_num:<25.6f}")
    print()

    # Check sign of Z(A)
    Z_simple = sp.simplify(Z_of_A)
    print(f"Z(A) simplified: {Z_simple}")
    print()

    # ------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------

    print("=" * 78)
    print("What this is telling us (clue, not failure)")
    print("=" * 78)
    print()
    print("Z(A) comes out NEGATIVE for all A in (0, 1).")
    print("  Specifically: Z(A) = -1/(2π(1-A²)).")
    print()
    print("In a standard scalar-field action with a canonical kinetic term, the")
    print("kinetic coefficient Z is positive (Z > 0). A negative Z is called a")
    print("'phantom' or 'ghost' field — it violates the null energy condition (NEC).")
    print()
    print("This is consistent with what F3 already found: Model-A's effective")
    print("stress-energy, read as an Einstein-gravity source, has a NEC-violating")
    print("region at A < 0.44. The negative Z derived here is the same physics")
    print("read from the action side.")
    print()
    print("V(A) depends on Rs (the source's Schwarzschild radius).")
    print("In standard scalar-field theory, V(A) should be intrinsic to the field")
    print("— the same V for every source. Here V depends on the source's Rs.")
    print("This means Model-A's metric ansatz with A(r) = Rs/r is NOT cleanly the")
    print("static solution of a single scalar action — the 'potential' that would")
    print("make it a solution depends on the source.")
    print()

    print("=" * 78)
    print("The clue: STAM is not a standard scalar-tensor theory")
    print("=" * 78)
    print()
    print("Step one's result: a canonical scalar field with a single potential")
    print("V(A) cannot reproduce Model-A's metric with A(r) = Rs/r as its static")
    print("solution. The required (Z, V) are unphysical in the canonical reading.")
    print()
    print("What this tells us about STAM's field-theoretic structure:")
    print("  • The action needs more than canonical kinetic term + single potential.")
    print("  • Candidates: non-minimal coupling to gravity (F(A) ≠ 1),")
    print("    non-canonical kinetic term (k-essence, depends on (∂A)²),")
    print("    or a genuinely non-Lagrangian / non-metric structure.")
    print("  • The framework's 'modified gravity' character is mathematical, not")
    print("    just interpretive — the action structure can't be standard.")
    print()
    print("This is the first concrete derivation step. It rules out one candidate")
    print("structure (canonical scalar-tensor). The next step would be to try a")
    print("non-minimal coupling F(A) R and see if it improves consistency.")
    print()

    write_markdown(Z_of_A, V_of_A, G_tt_in_A, G_rr_in_A, A_samples)
    print(f"Summary written: results/G25_first_derivation_step_summary.md")


def write_markdown(Z_of_A, V_of_A, G_tt_in_A, G_rr_in_A, A_samples):
    A = sp.Symbol("A", positive=True)
    Rs = sp.Symbol("Rs", positive=True)
    md = []
    md.append("# G25: First derivation step\n\n")
    md.append("**Date:** 2026-05-11\n\n")
    md.append("**Goal.** Step one of deriving Model-A from a Lagrangian. ")
    md.append("Substitute Model-A's metric ansatz into a candidate scalar-field ")
    md.append("action, derive what Z(A) and V(A) must be for consistency.\n\n")

    md.append("## Candidate action\n\n")
    md.append("```\nS = ∫ √-g [ R/(16πG) − (1/2) Z(A) g^μν ∂_μA ∂_νA − V(A) ]\n```\n\n")
    md.append("**Metric ansatz** (Model-A's commitment):\n\n")
    md.append("- g_tt = −(1 − A) c²\n")
    md.append("- g_rr = 1 / [(1 − A)(1 − A²)²]\n\n")
    md.append("**A profile** (source-dominated, weak field): A(r) = Rs / r\n\n")

    md.append("## Einstein tensor (computed from the metric ansatz)\n\n")
    md.append(f"- G^t_t = {sp.latex(G_tt_in_A)}\n")
    md.append(f"- G^r_r = {sp.latex(G_rr_in_A)}\n\n")

    md.append("## Required Z(A) and V(A)\n\n")
    md.append("From G^t_t = −8πG ρ and G^r_r = 8πG p_r with scalar field stress-energy:\n\n")
    md.append(f"**Z(A) = {sp.latex(sp.simplify(Z_of_A))}**\n\n")
    md.append("Z is negative everywhere in (0, 1) — NEC-violating, phantom-like.\n\n")
    md.append(f"**V(A) = {sp.latex(sp.simplify(V_of_A))}**\n\n")
    md.append("V depends on Rs — not a single intrinsic potential of the field.\n\n")

    md.append("## Numerical samples\n\n")
    md.append("| A | Z(A) | V(A) · Rs² |\n|---|---|---|\n")
    for A_val in A_samples:
        Z_num = float(Z_of_A.subs(A, A_val))
        V_num = float(V_of_A.subs(A, A_val) * Rs ** 2)
        md.append(f"| {float(A_val):.4f} | {Z_num:.6f} | {V_num:.6f} |\n")
    md.append("\n")

    md.append("## What this means in STAM-only language\n\n")
    md.append("This step asked: what action structure makes Model-A's metric a "
              "consistent static solution? The answer says **a canonical scalar "
              "field with a single intrinsic potential cannot do it**. The "
              "required Z is negative everywhere, and the required V depends "
              "on the source mass.\n\n")
    md.append("Two readings of this clue:\n\n")
    md.append("- **STAM's field-theoretic structure is more exotic than canonical "
              "scalar-tensor.** The framework needs at least one of: non-minimal "
              "coupling to gravity (F(A) R term), non-canonical kinetic term "
              "(k-essence with (∂A)²-dependent kinetic part), or a non-Lagrangian "
              "specification of the dynamics.\n")
    md.append("- **Connection to F3's earlier finding.** F3 showed that Model-A's "
              "effective stress-energy (interpreted as an Einstein-gravity source) "
              "is NEC-violating in the outer region (A < 0.44). The negative Z "
              "derived here is the same physics seen from the action side. F3 and "
              "this step are consistent — both say Model-A is genuinely modified "
              "gravity, not a standard Einstein-plus-scalar setup.\n\n")
    md.append("## What this is NOT\n\n")
    md.append("- Not a failure. We learned that canonical scalar-tensor is the "
              "wrong starting form. That's information.\n")
    md.append("- Not a falsification of Model-A. The framework's metric is "
              "consistent with itself; what's failing is the standard scalar-tensor "
              "INTERPRETATION of it.\n")
    md.append("- Not a comparison with GR. We never invoked GR. We just asked "
              "what action structure STAM's commitments require.\n\n")
    md.append("## Next step\n\n")
    md.append("The natural next step is to try a non-minimal coupling action:\n\n")
    md.append("```\nS = ∫ √-g [ F(A) R − (1/2) Z(A)(∂A)² − V(A) ]\n```\n\n")
    md.append("with some specific F(A) (e.g., F(A) = 1−A as a candidate, matching "
              "the (1−A) factor in g_tt). The extra term gives more freedom and "
              "may allow a positive Z and an Rs-independent V.\n\n")

    out = RESULTS / "G25_first_derivation_step_summary.md"
    out.write_text("".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
