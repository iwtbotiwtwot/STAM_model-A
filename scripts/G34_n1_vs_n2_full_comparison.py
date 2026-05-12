#!/usr/bin/env python3
"""
G34_n1_vs_n2_full_comparison.py

Systematic comparison: what would change if Model-A committed to n=1 instead
of its current n=2 in k_n(A) = (1-A)(1-A^2)^n?

Following G33's family analysis, both n=1 and n=2 satisfy the framework's
stated structural constraints (weak-field GR, horizon proper-time divergence,
pair structure). The committed n=2 has specific predictions; n=1 would give
different predictions.

This script computes each framework signature under both n=1 and n=2, and
flags where:
  - n=1 produces a cleaner result (cleaner fraction, simpler form)
  - n=2 produces a cleaner result
  - The result is independent of n (no change)
  - A specific result requires committing to one or the other

Author: Sean Brady / STAM Model-A
Date: 2026-05-11 (evening)
"""

from __future__ import annotations
import sys
from pathlib import Path
import math

import numpy as np
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

A = sp.Symbol("A", positive=True)
A_0_VAL = 1.0 / (12.0 * math.pi)


def k_n(n_value, A_sym=A):
    return (1 - A_sym) * (1 - A_sym**2) ** n_value


def f_n(n_value, A_sym=A):
    return 1 - k_n(n_value, A_sym)


# ============================================================================
# 1. Weak-field GR recovery (both n give k(0)=1, k'(0)=-1)
# ============================================================================
def check_weak_field():
    results = {}
    for n in [1, 2]:
        k = k_n(n)
        results[n] = {
            "k(0)": sp.simplify(k.subs(A, 0)),
            "k'(0)": sp.simplify(sp.diff(k, A).subs(A, 0)),
        }
    return results


# ============================================================================
# 2. Second-order g_rr coefficient at small A
# (G17's "third integer 3" appearance was the ratio Model-A:GR = 3, n=2 specific)
# ============================================================================
def second_order_grr():
    results = {}
    for n in [0, 1, 2, 3]:
        # g_rr = 1/k(A); expand to second order in A
        grr_series = sp.series(1 / k_n(n), A, 0, 3).removeO()
        # Coefficient of A^2
        coef_A2 = grr_series.coeff(A, 2)
        results[n] = sp.simplify(coef_A2)
    return results


# ============================================================================
# 3. QNM eikonal damping ratio at photon sphere (A=2/3)
# tau_n / tau_GR = sqrt(k_GR(2/3) / k_n(2/3)) = (9/5)^(n/2)
# ============================================================================
def qnm_ratio():
    results = {}
    for n in [1, 2]:
        ratio_sq = sp.Rational(9, 5) ** n
        ratio = sp.sqrt(ratio_sq)
        results[n] = {
            "symbolic": ratio,
            "float": float(ratio),
        }
    return results


# ============================================================================
# 4. F3-style NEC crossover location
# B_n(A) = f_n(A) - A*f_n'(A); find sign change
# ============================================================================
def find_NEC_crossover_symbolic(n_value):
    """Solve B_n(A) = 0 symbolically."""
    f_expr = f_n(n_value)
    fp = sp.diff(f_expr, A)
    B = sp.expand(f_expr - A * fp)
    # Factor and find roots
    roots = sp.solve(B, A)
    # Filter to (0, 1)
    in_range = [r for r in roots if r.is_real and 0 < float(r) < 1]
    return B, in_range


# ============================================================================
# 5. Proper time of radial free-fall divergence character near A=1
# tau ~ integral 1/sqrt(k_n(A)) ~ integral 1/(1-A)^((n+1)/2) near A=1
# ============================================================================
def proper_time_character(n_value):
    p = n_value + 1
    exponent = sp.Rational(p, 2)
    if exponent < 1:
        return "finite (no divergence)"
    elif exponent == 1:
        return "logarithmic"
    elif exponent < 2:
        return f"power-law: ~(1-A)^(-{p/2 - 1}) ~ (1-A)^({1 - p/2})"
    else:
        return f"strong power-law: ~(1-A)^(-{p/2 - 1})"


# ============================================================================
# 6. Integer-3 appearances cluster (recurring framework theme)
# G17's g_rr second-order ratio Model-A:GR
# ============================================================================
def integer_3_cluster():
    """Which integer appears as the g_rr second-order ratio for each n?"""
    results = {}
    grr_coefs = second_order_grr()
    grr_GR = grr_coefs[0]
    for n in [1, 2, 3]:
        ratio = sp.simplify(grr_coefs[n] / grr_GR)
        results[n] = ratio
    return results


# ============================================================================
# 7. Pair factor count and structural decomposition reading
# ============================================================================
def pair_factor_count(n_value):
    """k_n = (1-A) * (1-A^2)^n. Count: 1 'lower' factor + n pair factors."""
    return {
        "lower_factor_count": 1,  # the (1-A)
        "pair_factor_count": n_value,  # (1-A^2)^n
        "total_factor_count": 1 + n_value,
    }


# ============================================================================
# 8. Bekenstein-Hawking entropy connection
# Standard derivation (independent of n): S = A_horizon / (4 L_P^2)
# Via thermodynamic rule k_B T = ℏc|∇A|/(4π) at horizon
# This gives S = A/(4 L_P^2) regardless of n (since T depends on |∇A|_coord
# at horizon, which is 1/Rs same for any A=Rs/r profile).
# But the (2x2)=4 structural decomposition from G18 assumed n=2.
# Under n=1: only 1 pair factor; (2x2) decomposition reads differently.
# ============================================================================
def entropy_decomposition_reading(n_value):
    if n_value == 1:
        return ("Thermodynamic derivation: S = A/(4 L_P^2). G18's (2x2)=4 "
                "decomposition becomes (2x1)=2 under n=1 -- gives S = A/(2 L_P^2), "
                "which DOES NOT match Bekenstein-Hawking. Either G18 decomposition was "
                "n=2-specific (working backward from BH), or n=1 needs different reading.")
    elif n_value == 2:
        return ("Thermodynamic derivation: S = A/(4 L_P^2). G18's (2x2)=4 "
                "decomposition matches Bekenstein-Hawking. But the (2x2) assumed n=2; "
                "this is potentially circular if used to argue for n=2.")
    else:
        return "Not currently analyzed"


# ============================================================================
# 9. PPN parameters (gamma, beta) at first order — depend only on g_tt
# Independent of n. Both pass solar-system tests identically.
# ============================================================================
def ppn_first_order(n_value):
    return {"gamma": 1, "beta": 1, "note": "depends only on g_tt; n-independent at first order"}


# ============================================================================
# Main analysis
# ============================================================================
def main():
    print("=" * 78)
    print("G34: Systematic comparison — what does n=1 clean up vs n=2 (current)?")
    print("=" * 78)
    print()

    # Section 1: weak-field GR
    print("-" * 78)
    print("1. Weak-field GR recovery (k(0), k'(0))")
    print("-" * 78)
    wf = check_weak_field()
    for n, r in wf.items():
        print(f"   n={n}: k(0) = {r['k(0)']}, k'(0) = {r[chr(34)+chr(107)+chr(39)+chr(40)+chr(48)+chr(41)+chr(34)] if False else r['k\'(0)']}")
    print(f"   -> BOTH n=1 and n=2 satisfy weak-field GR at first order")
    print(f"   -> All solar-system tests pass identically (PPN gamma = beta = 1)")
    print()

    # Section 2: Second-order g_rr coefficient
    print("-" * 78)
    print("2. Second-order g_rr coefficient ratio (G17's integer-3 finding)")
    print("-" * 78)
    so = second_order_grr()
    cluster = integer_3_cluster()
    print(f"   g_rr second-order coefficient at small A:")
    for n in [0, 1, 2, 3]:
        ratio_str = f"   ratio MA:GR = {cluster.get(n, 'N/A')}" if n > 0 else "(GR baseline)"
        print(f"      n={n}: coefficient = {so[n]}  {ratio_str}")
    print()
    print(f"   -> n=2 gives ratio 3 (the 'third integer 3' in G17)")
    print(f"   -> n=1 gives ratio 2 (the integer becomes 2, breaks the 'three 3s' cluster)")
    print()

    # Section 3: QNM eikonal damping
    print("-" * 78)
    print("3. QNM eikonal damping ratio tau_n / tau_GR at photon sphere")
    print("-" * 78)
    qnm = qnm_ratio()
    for n, r in qnm.items():
        print(f"   n={n}: ratio = {r['symbolic']} = {r['float']:.4f}")
    print()
    print(f"   LIGO current precision: tau_obs/tau_GR ~ 1.0 +/- 0.2")
    print(f"   -> n=1 gives 1.342 (MARGINAL with LIGO)")
    print(f"   -> n=2 gives 1.800 (clearly OUTSIDE current LIGO band)")
    print(f"   -> n=1 fits LIGO data better (or is less in tension)")
    print()

    # Section 4: NEC crossover
    print("-" * 78)
    print("4. F3-style NEC crossover (effective rho sign change)")
    print("-" * 78)
    for n in [1, 2]:
        B, roots = find_NEC_crossover_symbolic(n)
        print(f"   n={n}:")
        print(f"      B_{n}(A) = {sp.factor(B)}")
        print(f"      Real roots in (0,1): {[sp.nsimplify(r, rational=True) for r in roots]}")
        if roots:
            for r in roots:
                print(f"        symbolic: {r}, numerical: {float(r):.6f}")
    print()
    print(f"   -> n=1 gives EXACTLY A_crit = 1/2 (clean fraction, midpoint of saturation)")
    print(f"   -> n=2 gives A_crit ≈ 0.44 (no clean fraction)")
    print()

    # Section 5: Proper time divergence character
    print("-" * 78)
    print("5. Proper time divergence character at A=1")
    print("-" * 78)
    for n in [0, 1, 2, 3]:
        char = proper_time_character(n)
        print(f"   n={n}: {char}")
    print()
    print(f"   -> n=1 gives LOGARITHMIC divergence (matches PDF Section 5 wording)")
    print(f"   -> n=2 gives 1/sqrt(1-A) divergence (faster than logarithmic)")
    print()

    # Section 6: pair-factor structural reading
    print("-" * 78)
    print("6. Pair-factor structure and entropy decomposition")
    print("-" * 78)
    for n in [1, 2]:
        pf = pair_factor_count(n)
        ed = entropy_decomposition_reading(n)
        print(f"   n={n}:")
        print(f"      (1-A) factors: {pf['lower_factor_count']}")
        print(f"      pair factors (1-A^2): {pf['pair_factor_count']}")
        print(f"      Total: {pf['total_factor_count']}")
        print(f"      Entropy reading: {ed}")
        print()

    # Section 7: PPN
    print("-" * 78)
    print("7. PPN parameters (first-order solar system tests)")
    print("-" * 78)
    ppn = ppn_first_order(1)
    print(f"   gamma = {ppn['gamma']}, beta = {ppn['beta']}")
    print(f"   {ppn['note']}")
    print(f"   -> Both n=1 and n=2 pass solar system tests identically")
    print()

    # Section 8: cosmology and thermodynamics
    print("-" * 78)
    print("8. Cosmology and thermodynamics")
    print("-" * 78)
    print(f"   - V_3 = alpha/A + beta/(1-A): independent of k(A); same for n=1 and n=2")
    print(f"   - Bridge term b = A_0 * c/H_0: independent of k(A); same for both")
    print(f"   - Thermodynamic rule k_B T = hbar*c*|grad A|/(4pi) at boundary:")
    print(f"     gives Hawking temperature regardless of n (depends on coord |grad A| at horizon)")
    print(f"   -> Cosmological and thermodynamic results are UNCHANGED by n=1 vs n=2")
    print()

    # Section 9: summary table
    print("=" * 78)
    print("Summary table: cleaner under n=1, cleaner under n=2, or independent")
    print("=" * 78)
    print()
    print(f"   {'Feature':<35} {'n=1':<25} {'n=2 (current)':<25} {'Cleaner':<10}")
    print(f"   {'-'*35:<35} {'-'*25:<25} {'-'*25:<25} {'-'*10}")
    rows = [
        ("Weak-field GR",                "k(0)=1, k'(0)=-1",         "k(0)=1, k'(0)=-1",         "tie"),
        ("Proper-time divergence",       "logarithmic",              "(1-A)^(-1/2)",             "n=1*"),
        ("F3 NEC crossover",             "1/2 EXACTLY",              "≈0.44",                    "n=1"),
        ("QNM ratio tau/tau_GR",         "3/sqrt(5) ≈ 1.342",       "9/5 = 1.800",              "depends on obs"),
        ("g_rr 2nd-order ratio MA:GR",   "2:1 (clean)",              "3:1 (clean, in '3s' cluster)", "depends on theme"),
        ("Pair-factor count",            "1 (simpler)",              "2",                        "n=1 simpler"),
        ("Entropy decomp via (2x2)=4",   "BREAKS (gives 2, not 4)",  "matches BH 4",             "n=2"),
        ("Solar system PPN",             "gamma=beta=1",             "gamma=beta=1",             "tie"),
        ("Cosmological branch",          "unchanged",                "unchanged",                "tie"),
        ("Thermodynamic rule",           "unchanged",                "unchanged",                "tie"),
    ]
    for feat, n1, n2, win in rows:
        print(f"   {feat:<35} {n1:<25} {n2:<25} {win:<10}")
    print()
    print("   * 'logarithmic' matches PDF Section 5 wording; but both diverge, so")
    print("     this is a wording-match issue, not a derivation issue.")
    print()

    # Markdown summary
    md = []
    md.append("# G34: What does n=1 clean up vs Model-A's current n=2?")
    md.append("")
    md.append("**Date:** 2026-05-11 (evening)")
    md.append("")
    md.append("Following G33's discovery that k_n(A) = (1-A)(1-A^2)^n admits a 1-parameter "
              "family of viable choices, this script systematically compares n=1 and n=2 "
              "across the framework's predictions to identify what each commitment cleans up.")
    md.append("")
    md.append("## Comparison table")
    md.append("")
    md.append("| Feature | n=1 | n=2 (current Model-A) | Cleaner |")
    md.append("|---|---|---|---|")
    for feat, n1, n2, win in rows:
        md.append(f"| {feat} | {n1} | {n2} | {win} |")
    md.append("")
    md.append("## n=1 cleans up")
    md.append("")
    md.append("- **F3 NEC crossover at exactly A = 1/2** (symbolically: A²(1-2A) = 0 → A=1/2 exactly). "
              "n=2 gives A ≈ 0.44 with no clean fraction.")
    md.append("- **Proper-time divergence is logarithmic** at A=1 (~ ln(1/(1-A))). Matches PDF "
              "Section 5 wording. n=2 has stronger (1-A)^(-1/2) divergence.")
    md.append("- **QNM ratio tau/tau_GR = 3/sqrt(5) ≈ 1.342**, marginal with LIGO precision "
              "tau_obs/tau_GR ≈ 1.0 ± 0.2. n=2's 1.800 is clearly outside that band.")
    md.append("- **Simpler structural form**: 1 pair factor instead of 2. Fewer modification "
              "factors total.")
    md.append("")
    md.append("## n=2 cleans up (or is required for)")
    md.append("")
    md.append("- **G17's 'third integer 3' appearance** (g_rr 2nd-order ratio Model-A:GR = 3 exactly). "
              "Under n=1 this becomes 2. The 'integer 3 appears in three independent places' "
              "framework theme depends on n=2.")
    md.append("- **G18's (2×2)=4 entropy decomposition** (matching Bekenstein-Hawking S = A/(4L_P^2)). "
              "Under n=1, the structural counting gives 2, not 4. But the thermodynamic derivation "
              "(via k_B T = hbar c |grad A|/(4pi)) still gives BH entropy for both n=1 and n=2 — "
              "so the G18 (2×2) was the n=2-specific structural reading, not a derivation that "
              "FORCES n=2.")
    md.append("- **9/5 = 1.800 clean fraction** for QNM ratio. Under n=1 this becomes 3/sqrt(5) "
              "≈ 1.342 (still clean in form, less iconic).")
    md.append("")
    md.append("## Independent of n=1 vs n=2 (no change)")
    md.append("")
    md.append("- Weak-field GR (PPN gamma = beta = 1)")
    md.append("- Cosmological branch (V_3, bridge term)")
    md.append("- Thermodynamic rule and resulting Hawking/Unruh/de Sitter temperatures")
    md.append("- Strong-field landmarks (ISCO, photon sphere, horizon at A = 1/3, 2/3, 1)")
    md.append("- Quantum interpretation (resolved/unresolved A, ledger, presentism)")
    md.append("- BH ontology (two-face refinement, no interior)")
    md.append("- A_0 = 1/(12π) and the (4π × 3) decomposition")
    md.append("")
    md.append("## Reading")
    md.append("")
    md.append(
        "**n=1 cleans up several features**: cleaner NEC crossover (1/2), simpler form, "
        "logarithmic divergence matching PDF wording, QNM ratio closer to LIGO. "
        "**n=2 cleans up two structural-theme features**: integer-3 cluster (G17) and (2×2)=4 "
        "entropy decomposition (G18). But both are potentially n=2-specific readings rather "
        "than forcing arguments — they could be reinterpreted under n=1 with different counting."
    )
    md.append("")
    md.append(
        "**Honest assessment**: the framework's commitment to n=2 was made partly because of "
        "the G17 and G18 structural-theme readings. Under n=1, those readings don't disappear — "
        "they become different (integer-2 in G17, integer-2 in entropy structural decomp). "
        "The Bekenstein-Hawking 4-area-per-entry result still holds via the thermodynamic "
        "rule independent of n. So the framework can equally well commit to n=1 without "
        "breaking the central results — just with different structural-theme readings."
    )
    md.append("")
    md.append("## What this surfaces")
    md.append("")
    md.append(
        "The choice between n=1 and n=2 is a STRUCTURAL CHOICE that the framework has the "
        "freedom to make. Neither is uniquely forced by the framework's stated constraints. "
        "Choosing one over the other tilts which structural themes are foregrounded:"
    )
    md.append("")
    md.append("- **n=1 emphasizes**: clean midpoint A=1/2 in stress-energy, logarithmic horizon, "
              "simpler form, closer to LIGO data.")
    md.append("- **n=2 emphasizes**: 'three integer-3 appearances' cluster, (2×2)=4 entropy "
              "decomposition, pair-squared structure.")
    md.append("")
    md.append(
        "A real fixing point (observational measurement, or non-circular structural derivation) "
        "would settle this. Until then, the framework has a meaningful 1-parameter ambiguity "
        "in its strong-field metric commitment."
    )
    md.append("")
    md.append("## Practical recommendation for PDF v0.4")
    md.append("")
    md.append(
        "Flag the 1-parameter ambiguity explicitly. Note that Model-A currently commits to n=2 "
        "based on structural-theme considerations (integer-3 cluster, pair-squared form), but "
        "n=1 is a viable alternative that some features (NEC crossover, LIGO consistency) "
        "favor. The framework's testable predictions depend on this choice; future observational "
        "or theoretical fixing of n is an open priority."
    )

    md_path = RESULTS / "G34_n1_vs_n2_full_comparison_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Summary saved: {md_path}")


if __name__ == "__main__":
    main()
