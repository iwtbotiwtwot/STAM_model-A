#!/usr/bin/env python3
"""
G27_third_derivation_step.py

Step three of the derivation. Same machinery as G26, but F(A) = (1-A²)² —
the doubled-pair coupling, matching the (1-A²)² modification factor in k(A)
directly.

Question: does this "doubled-pair" coupling fix the Z sign-change problem
from step 2 (Z went negative for A > ~0.39)?

We're not aiming for a perfect fit. We're seeing how the math responds when
we keep increasing the "pair-strength" of the gravitational coupling.

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


def derive_Z_and_V(F_of_A, label: str):
    """Derive Z(A) and V(A) required for Model-A's metric to be a static
    spherically symmetric solution, given a non-minimal coupling F(A).

    Uses Model-A's metric ansatz (h = 1-A, k = (1-A)(1-A²)²) and the
    source-dominated profile A(r) = Rs/r.
    """
    A = sp.Symbol("A", positive=True)
    r, Rs = sp.symbols("r Rs", positive=True)

    h_of_A = 1 - A
    k_of_A = (1 - A) * (1 - A**2) ** 2

    A_of_r = Rs / r
    h_r = h_of_A.subs(A, A_of_r)
    k_r = k_of_A.subs(A, A_of_r)
    F_r = F_of_A.subs(A, A_of_r)

    h_prime = sp.diff(h_r, r)
    k_prime = sp.diff(k_r, r)
    F_prime = sp.diff(F_r, r)
    F_double_prime = sp.diff(F_prime, r)
    A_prime = sp.diff(A_of_r, r)

    # Einstein tensor mixed components
    G_tt = (k_r - 1) / r**2 + k_prime / r
    G_rr = (k_r - 1) / r**2 + k_r * h_prime / (r * h_r)

    # Box operator on F
    box_F = k_r * F_double_prime + (F_prime / 2) * (
        4 * k_r / r - k_r * h_prime / h_r - k_prime
    )

    # Modified Einstein equations:
    # t-t: F G^t_t + (k h')/(2h) F' + □F = -8π ρ
    # r-r: F G^r_r - k F'' - (k'/2) F' + □F = 8π p_r
    eq_tt = F_r * G_tt + (k_r * h_prime / (2 * h_r)) * F_prime + box_F
    eq_rr = F_r * G_rr - k_r * F_double_prime - (k_prime / 2) * F_prime + box_F

    # ρ + p_r = (eq_rr - eq_tt)/(8π) = Z k (A')²
    Z_expr = (eq_rr - eq_tt) / (8 * sp.pi * k_r * A_prime**2)
    Z = sp.simplify(Z_expr.subs(r, Rs / A))

    # ρ - p_r = 2V = (-eq_tt - eq_rr)/(8π)
    V_expr = -(eq_tt + eq_rr) / (16 * sp.pi)
    V = sp.simplify(V_expr.subs(r, Rs / A))

    return Z, V


def main():
    print("=" * 78)
    print("G27: Third derivation step — F(A) = (1-A²)² (doubled-pair coupling)")
    print("=" * 78)
    print()
    print("Candidate action:")
    print("   S = ∫ √-g [ F(A) R / (16πG) − (1/2) Z(A)(∂A)² − V(A) ]")
    print()
    print("F(A) = (1 - A²)²  (doubled-pair coupling)")
    print("h(A) = 1 - A")
    print("k(A) = (1-A)(1-A²)² = (1-A)³(1+A)²")
    print()
    print("A profile: A(r) = Rs/r")
    print()

    A = sp.Symbol("A", positive=True)
    Rs = sp.Symbol("Rs", positive=True)
    F_doubled = (1 - A**2) ** 2

    Z, V = derive_Z_and_V(F_doubled, "F = (1-A²)²")

    print("=" * 78)
    print("Z(A) result for F = (1-A²)²:")
    print("=" * 78)
    print()
    print(f"   Z(A) = {Z}")
    print()
    Z_factored = sp.factor(Z)
    print(f"   Factored: Z(A) = {Z_factored}")
    print()

    # Sample Z values
    A_samples = [
        sp.Rational(1, 100),
        sp.Rational(1, 10),
        sp.Rational(1, 3),
        sp.Rational(1, 2),
        sp.Rational(2, 3),
        sp.Rational(9, 10),
        sp.Rational(99, 100),
    ]
    print(f"{'A':<10}{'Z(A)':<25}{'sign':<10}")
    for A_val in A_samples:
        z_val = float(Z.subs(A, A_val))
        sign = "POSITIVE ✓" if z_val > 0 else "negative ✗"
        print(f"{float(A_val):<10.4f}{z_val:<25.6f}{sign:<10}")
    print()

    # Check for sign changes
    A_test = sp.symbols("A_test", positive=True)
    # Try to find roots of Z numerator
    Z_num = sp.together(Z).as_numer_denom()[0]
    Z_num_simplified = sp.simplify(Z_num)
    print(f"Z numerator: {Z_num_simplified}")
    print()
    try:
        Z_roots = sp.solve(Z_num_simplified, A)
        Z_real_roots = [r for r in Z_roots if r.is_real and 0 < r < 1]
        if Z_real_roots:
            print(f"Z = 0 at: {Z_real_roots}")
            for root in Z_real_roots:
                print(f"  Numeric: A ≈ {float(root):.6f}")
        else:
            print("No real roots of Z in (0, 1) — Z does not change sign in the manifold range.")
    except Exception as e:
        print(f"Could not solve for roots analytically: {e}")
    print()

    print("=" * 78)
    print("V(A) result for F = (1-A²)²:")
    print("=" * 78)
    print()
    V_dimensionless = sp.simplify(V * Rs**2)
    print(f"   V(A) · Rs² = {V_dimensionless}")
    print()
    V_factored = sp.factor(V)
    print(f"   Factored V(A) = {V_factored}")
    print()

    # Sample V values
    print(f"{'A':<10}{'V·Rs² (numeric)':<25}")
    for A_val in A_samples:
        v_val = float(V_dimensionless.subs(A, A_val))
        print(f"{float(A_val):<10.4f}{v_val:<25.6f}")
    print()

    # ============================================================
    # Comparison table across F choices
    # ============================================================

    print("=" * 78)
    print("Comparison across F(A) choices:")
    print("=" * 78)
    print()

    F_choices = [
        (sp.Integer(1), "F = 1 (step 1)"),
        ((1 - A**2), "F = (1-A²) (step 2)"),
        ((1 - A**2)**2, "F = (1-A²)² (step 3, this run)"),
    ]

    print(f"{'A':<8}", end="")
    for _, label in F_choices:
        print(f"{label:<32}", end="")
    print()

    for A_val in A_samples:
        print(f"{float(A_val):<8.4f}", end="")
        for F_func, _ in F_choices:
            Z_F, _ = derive_Z_and_V(F_func, "")
            try:
                z_val = float(Z_F.subs(A, A_val))
                sign_mark = "+" if z_val > 0 else "-"
                print(f"  Z={z_val:9.5f} {sign_mark}              ", end="")
            except Exception:
                print(f"  Z=ERR                          ", end="")
        print()
    print()

    print("=" * 78)
    print("Interpretation")
    print("=" * 78)
    print()
    print("Reading the pattern across F=1, F=(1-A²), F=(1-A²)²:")
    print()
    print("Watch what happens to the sign-change A_crit as we add more pair-power")
    print("to the gravitational coupling.")
    print()


if __name__ == "__main__":
    main()
