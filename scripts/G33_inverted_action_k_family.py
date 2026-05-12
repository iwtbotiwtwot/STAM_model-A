#!/usr/bin/env python3
"""
G33_inverted_action_k_family.py

Swing in the dark: treat A as fundamental, derive what k(A) can be.

Standard approach (G25-G27): g_munu fundamental, A is a scalar field on it,
derive Z(A) and V(A) from action variation. Gave ghost regions and Rs-dependent V.

Inverted approach (this script): A is fundamental. g_munu is constructed FROM
A via g_tt = -(1-A)c^2 and g_rr = 1/k(A). The question becomes: what k(A) is
consistent with the framework's structural commitments?

Constraints from the framework:
  (i) Weak-field GR at first order in A: k(A) = 1 - A + O(A^2)
      => k(0) = 1, k'(0) = -1
  (ii) Horizon at A=1 with proper-time divergence:
      k(A) ~ (1-A)^p near A=1 with p >= 2 for tau -> infinity
  (iii) Pair structure (1-A^2) appearing throughout framework results

Natural family satisfying (i)-(iii): k_n(A) = (1-A)(1-A^2)^n for n in Z_{>=1}.
  - n = 0:  k = (1-A) -- pure Schwarzschild, fails (ii) (proper time finite)
  - n = 1:  k = (1-A)(1-A^2) -- minimal STAM-like, proper time diverges
  - n = 2:  k = (1-A)(1-A^2)^2 -- Model-A's committed choice
  - n = 3:  k = (1-A)(1-A^2)^3
  - n = 4:  k = (1-A)(1-A^2)^4
  ...

What this script does:
  1. Verify constraints (i)-(iii) for each n.
  2. Compute the QNM eikonal damping ratio tau_n / tau_GR at the photon sphere.
     This scales as (9/5)^(n/2) -- so n parameterizes a continuum of predictions.
  3. Compute the F3 effective stress-energy NEC crossover for each n.
  4. Compute the proper-time divergence rate at horizon for each n.
  5. Compute G18-style ledger entropy decomposition: each n gives a different
     area-per-entry coefficient.
  6. Identify what structural argument (if any) picks out n=2 specifically.
  7. Identify what observational target would discriminate.

The expected MISS: the constraints don't uniquely force n=2. Model-A's commitment
to n=2 turns out to be a structural choice within a 1-parameter family, not a
derivation. The miss tells us where the actual fixing point would have to come from
(observation OR an additional structural principle we haven't identified).

Author: Sean Brady / STAM Model-A
Date: 2026-05-11 (evening) -- swing in the dark, eyes closed
"""

from __future__ import annotations
import sys
from pathlib import Path
from fractions import Fraction

import numpy as np
import sympy as sp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Symbolic setup
# ============================================================================
A = sp.Symbol("A", positive=True)


def k_family(n_value):
    """k_n(A) = (1 - A) * (1 - A^2)^n. Returns sympy expression."""
    return (1 - A) * (1 - A**2) ** n_value


def f_of_A(k_expr):
    """f(A) = 1 - k(A). Effective enclosed-mass function."""
    return 1 - k_expr


# ============================================================================
# Constraint checks
# ============================================================================
def check_weak_field_GR(k_expr, n_value):
    """Verify k(0) = 1 and k'(0) = -1."""
    k0 = sp.simplify(k_expr.subs(A, 0))
    kp = sp.diff(k_expr, A)
    kp0 = sp.simplify(kp.subs(A, 0))
    return k0, kp0


def check_horizon_divergence(k_expr, n_value):
    """Verify the multiplicity p of (1-A) factor at A=1, requiring p >= 2 for
    proper-time divergence."""
    eps = sp.Symbol("eps", positive=True)
    k_near = sp.series(k_expr.subs(A, 1 - eps), eps, 0, n_value + 4).removeO()
    # The leading order is eps^(n+1) for k_n = (1-A)(1-A^2)^n
    p = n_value + 1
    coef = sp.simplify(k_near.coeff(eps, p))
    return p, coef


def proper_time_divergence_rate(n_value):
    """For k_n = (1-A)(1-A^2)^n, proper time from rest at infinity to A=1:
       tau ~ integral 1/(1-A)^((n+1)/2) (modulo nonsingular factors)
    Diverges if (n+1)/2 >= 1, i.e. n >= 1.
    Returns (diverges, asymptotic_rate)."""
    rate = sp.Rational(n_value + 1, 2)
    return rate >= 1, rate


# ============================================================================
# QNM eikonal damping ratio
# ============================================================================
# Photon sphere at A = 2/3 (location independent of k, depends only on g_tt)
# Lyapunov exponent for null orbits: lambda^2 ~ k(A_ph) at fixed orbit
# tau = 1/lambda, so tau_n / tau_GR = sqrt(k_GR(A_ph) / k_n(A_ph))
# k_GR(2/3) = 1/3
# k_n(2/3) = (1/3)(5/9)^n
# tau_n / tau_GR = sqrt(1 / (5/9)^n) = (9/5)^(n/2)

def qnm_ratio(n_value):
    return sp.Rational(9, 5) ** sp.Rational(n_value, 2)


# ============================================================================
# F3-style effective stress-energy bracket
# ============================================================================
# B(A) = f(A) - A * f'(A) for the canonical F3 reading (A = Rs/r)
# Find sign-change A_crit_n

def f3_bracket(k_expr):
    f_expr = f_of_A(k_expr)
    fp = sp.diff(f_expr, A)
    return sp.simplify(f_expr - A * fp)


def find_NEC_crossover(k_expr, n_value, n_scan=10000):
    """Numerical scan for sign change of B(A) = f - A*f' on (1e-6, 0.999)."""
    B = f3_bracket(k_expr)
    B_simplified = sp.simplify(B)
    if B_simplified == 0:
        return []  # B identically zero (Schwarzschild case n=0)
    B_func = sp.lambdify(A, B, "numpy")
    A_grid = np.linspace(1e-6, 0.999, n_scan)
    B_vals = np.asarray(B_func(A_grid))
    if B_vals.ndim == 0:
        B_vals = np.full(n_scan, float(B_vals))
    sign_change = np.where(np.diff(np.sign(B_vals)) != 0)[0]
    crossings = []
    for i in sign_change:
        lo, hi = A_grid[i], A_grid[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if B_func(mid) * B_func(lo) < 0:
                hi = mid
            else:
                lo = mid
        crossings.append(0.5 * (lo + hi))
    return crossings


# ============================================================================
# Ledger-style decomposition (each n gives different area-per-entry coefficient)
# ============================================================================
# Following G18's argument: the entropy coefficient comes from how (1-A^2)
# factors decompose at the boundary. For k_n = (1-A)(1-A^2)^n:
#   - The (1-A) factor encodes "one gravity-bridge unit" (factor 2)
#   - The (1-A^2)^n factor encodes "n pair structures"
# G18 derived S = A_horizon/(4 L_P^2) via (2x2) = 4 area-per-entry.
# That assumed n=2 (two pair factors) interpreted as "two faces x pair".
# Under general n, the entropy coefficient becomes 2 x 2^n = 2^(n+1)... or similar.
# This is structural reading rather than derivation.

def ledger_area_per_entry(n_value):
    """Heuristic structural counting: 2 (gravity bridge) x 2 (per pair factor)^n
    No, that doubles per pair factor which doesn't match G18.
    G18's reading: 2 (two faces) x 2 (gravity bridge) = 4, with two pair factors
    contributing one "2" each via the two-face structure.
    Under n=1: only one pair factor, so only one "2" -> area-per-entry = 2
    Under n=2: two pair factors -> 2x2 = 4 (Bekenstein-Hawking)
    Under n=3: three pair factors -> 2x3 = 6 or 2x2x2 = 8 (ambiguous)
    The actual right answer requires a real derivation, not heuristic counting."""
    return f"2 * (factor structure for n={n_value})  -- structural reading only"


# ============================================================================
# Main analysis
# ============================================================================
def main():
    print("=" * 78)
    print("G33: Inverted action approach -- k(A) family from A-fundamental view")
    print("=" * 78)
    print()
    print("Swing in the dark: A is fundamental, derive what k(A) can be.")
    print()

    n_values = [0, 1, 2, 3, 4, 5]

    print("Family: k_n(A) = (1-A) * (1-A^2)^n")
    print()
    print(f"  {'n':>3} {'k(0)':>8} {'k_prime(0)':>11} {'p (mult at A=1)':>16} "
          f"{'tau div?':>10} {'tau/tau_GR (QNM)':>17}")
    print("  " + "-" * 78)

    results = []
    for n in n_values:
        k_expr = k_family(n)
        k0, kp0 = check_weak_field_GR(k_expr, n)
        p, coef = check_horizon_divergence(k_expr, n)
        diverges, rate = proper_time_divergence_rate(n)
        qnm = qnm_ratio(n)
        qnm_float = float(qnm)

        crossings = find_NEC_crossover(k_expr, n)
        crossing_str = (
            ", ".join(f"{c:.4f}" for c in crossings) if crossings else "none"
        )

        results.append({
            "n": n,
            "k_expr": k_expr,
            "k0": k0,
            "kp0": kp0,
            "p": p,
            "diverges": diverges,
            "qnm": qnm,
            "qnm_float": qnm_float,
            "NEC_crossings": crossings,
        })

        print(f"  {n:>3} {k0!s:>8} {kp0!s:>11} {p:>16} {('YES' if diverges else 'no'):>10} "
              f"{qnm_float:>17.4f}")

    print()
    print("All n satisfy weak-field GR (k(0)=1, k'(0)=-1). All n >= 1 satisfy")
    print("proper-time divergence at A=1. So all n >= 1 satisfy the framework's")
    print("stated constraints. n=0 is pure Schwarzschild (no STAM departure).")
    print()

    # ====================================================================
    # NEC crossover comparison
    # ====================================================================
    print("F3-style NEC crossover locations (where effective rho changes sign):")
    for r in results:
        if r["n"] == 0:
            continue
        crossings_str = ", ".join(f"{c:.4f}" for c in r["NEC_crossings"])
        print(f"  n = {r['n']}: A_crit = {crossings_str}")
    print()

    # ====================================================================
    # Observational comparison with LIGO ringdown
    # ====================================================================
    print("=" * 78)
    print("LIGO QNM ringdown comparison (eikonal, spinless BH approximation):")
    print("=" * 78)
    print("  tau_obs/tau_GR ~= 1.0 +/- 0.2  (current LIGO precision, mass-independent)")
    print()
    print(f"  {'n':>3} {'tau_n/tau_GR':>15} {'consistent with LIGO?':>25}")
    print("  " + "-" * 50)
    for r in results:
        ratio = r["qnm_float"]
        consistent = "Yes (within 1-sigma)" if abs(ratio - 1.0) < 0.2 else (
            f"NO ({ratio:.3f} far outside 1.0 +/- 0.2)" if abs(ratio - 1.0) > 0.5 else
            f"Marginal ({ratio:.3f})"
        )
        print(f"  {r['n']:>3} {ratio:>15.4f} {consistent:>25}")
    print()
    print("  Note: LIGO BHs are spinning (Kerr-like); spinless Model-A analog")
    print("  is not yet derived. So this comparison is indicative, not falsifying.")
    print()

    # ====================================================================
    # Where Model-A's n=2 sits
    # ====================================================================
    print("=" * 78)
    print("Where does Model-A's n=2 sit in this family?")
    print("=" * 78)
    print()
    print("  Constraints satisfied by n=2:")
    print("    - Weak-field GR (k(0)=1, k'(0)=-1)              shared with all n >= 0")
    print("    - Horizon at A=1 with proper-time divergence    shared with all n >= 1")
    print("    - Pair structure (1-A^2) factor                 shared with all n >= 1")
    print()
    print("  Predictions specific to n=2:")
    print("    - QNM eikonal damping tau_MA/tau_GR = 9/5 = 1.8000")
    print("    - Pair-squared form (1-A^2)^2  <-> G18's (2x2)=4 decomposition")
    print("    - 4 modification factors total in (m,n) metric scan -- pair-squared")
    print()
    print("  Other n with no specific framework anchor:")
    print(f"    - n=1: tau/tau_GR = {float(qnm_ratio(1)):.4f}, pair-1, simpler structure")
    print(f"    - n=3: tau/tau_GR = {float(qnm_ratio(3)):.4f}, pair-cubed")
    print(f"    - n=4: tau/tau_GR = {float(qnm_ratio(4)):.4f}, pair-fourth")
    print()
    print("  Heuristic structural motivation for n=2 (NOT a derivation):")
    print("    - Two-face refinement: 'two faces x pair' -> (1-A^2)^2")
    print("    - G18 ledger entropy decomposition: (2x2) area-per-entry")
    print("    - But neither uniquely forces n=2 vs n=1 or n=3")
    print()

    # ====================================================================
    # Plot the family
    # ====================================================================
    A_grid = np.linspace(0.001, 0.999, 500)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel 1: k_n(A) shapes
    ax = axes[0]
    for r in results:
        n = r["n"]
        k_func = sp.lambdify(A, r["k_expr"], "numpy")
        k_vals = k_func(A_grid)
        label = f"n={n}" + (" (Model-A)" if n == 2 else "")
        if n == 2:
            ax.plot(A_grid, k_vals, "b-", lw=3, label=label)
        elif n == 0:
            ax.plot(A_grid, k_vals, "k:", lw=2, label=f"n=0 (Schwarzschild)")
        else:
            ax.plot(A_grid, k_vals, "-", lw=1.5, alpha=0.7, label=label)
    ax.set_xlabel("A")
    ax.set_ylabel("k(A) = (1-A)(1-A^2)^n")
    ax.set_title("k(A) family shapes")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 2)

    # Panel 2: QNM ratio vs n
    ax = axes[1]
    n_arr = np.array([r["n"] for r in results])
    qnm_arr = np.array([r["qnm_float"] for r in results])
    ax.plot(n_arr, qnm_arr, "bo-", markersize=8, lw=2)
    ax.axhline(1.0, color="g", linestyle="--", lw=1.5, label="GR (tau/tau_GR = 1)")
    ax.axhline(1.0 + 0.2, color="orange", linestyle=":", lw=1.0)
    ax.axhline(1.0 - 0.2, color="orange", linestyle=":", lw=1.0, label="LIGO ~1-sigma band")
    ax.axvline(2, color="b", linestyle=":", lw=1.5, label="Model-A's n=2")
    for r in results:
        ax.text(r["n"], r["qnm_float"] + 0.1, f"{r['qnm_float']:.3f}",
                ha="center", fontsize=9)
    ax.set_xlabel("n")
    ax.set_ylabel("tau_n / tau_GR  (QNM eikonal)")
    ax.set_title("QNM damping ratio scales as (9/5)^(n/2)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(0.5, 4.0)

    fig.suptitle("G33: k(A) family from A-fundamental view -- swing in the dark")
    fig.tight_layout()
    plot_path = PLOTS / "G33_k_family.png"
    fig.savefig(plot_path, dpi=140, bbox_inches="tight")
    print(f"Plot saved: {plot_path}")
    print()

    # ====================================================================
    # The honest miss outcome
    # ====================================================================
    print("=" * 78)
    print("Outcome of the swing: HONEST MISS")
    print("=" * 78)
    print()
    print("The constraints from the framework (weak-field GR + proper-time divergence")
    print("+ pair structure) DO NOT uniquely fix n. Model-A's k(A) = (1-A)(1-A^2)^2")
    print("is one member of a 1-parameter family {k_n}_{n>=1}, all satisfying the")
    print("stated structural constraints.")
    print()
    print("The 'three principles force k(A)' claim in PDF Section 4 is overstated.")
    print("Three principles NARROW k(A) to a family; they don't uniquely specify n=2.")
    print()
    print("What we learn from the miss:")
    print("  1. The Lagrangian-gap remains. Inverting the action principle doesn't")
    print("     close it -- it just reveals that the constraints admit a family.")
    print("  2. Model-A's specific n=2 needs a fixing point we don't yet have.")
    print("     Two candidate fixing points:")
    print("       (a) An observational measurement of tau/tau_GR for a spinless BH")
    print("           (would directly pick out n).")
    print("       (b) An additional structural principle (e.g., a precise derivation")
    print("           of the G18 area-per-entry coefficient WITHOUT assuming n=2.)")
    print("  3. The pair structure (1-A^2) IS load-bearing -- it appears in all viable")
    print("     k_n forms. The number of pair factors (n) is the open question.")
    print("  4. The (9/5)^(n/2) scaling is a concrete framework prediction. Different")
    print("     n give measurably different QNM signatures.")
    print()
    print("Where this points us:")
    print("  - Push G1 toward an exact (non-eikonal) QNM computation for several n,")
    print("    not just n=2. The exact spectrum may discriminate more sharply than")
    print("    eikonal.")
    print("  - Push G18 toward a derivation of area-per-entry WITHOUT assuming the")
    print("    n=2 structure (i.e., from ledger counting at a single Planck cell).")
    print("    Would either fix n=2 structurally OR show the (2x2)=4 was n-circular.")
    print("  - Accept that the framework currently has a 1-parameter ambiguity in")
    print("    the strong-field metric, and update the PDF v0.4 accordingly.")
    print()

    # ====================================================================
    # Markdown summary
    # ====================================================================
    md = []
    md.append("# G33: Inverted action approach -- k(A) family from A-fundamental view")
    md.append("")
    md.append("**Date:** 2026-05-11 (evening)")
    md.append("")
    md.append("**Swing in the dark:** treat A as fundamental, derive what k(A) can be.")
    md.append("")
    md.append("## The setup")
    md.append("")
    md.append(
        "Standard scalar-tensor approach (G25-G27): assume g_munu fundamental, A is "
        "scalar field, derive Z(A) and V(A) from action variation. Result: ghost regions "
        "and Rs-dependent V."
    )
    md.append("")
    md.append(
        "Inverted approach (this script): A is fundamental, g_munu constructed from A "
        "via g_tt = -(1-A)c^2 and g_rr = 1/k(A). Question: what k(A) is consistent with "
        "the framework's structural commitments?"
    )
    md.append("")
    md.append("Constraints from framework:")
    md.append("- Weak-field GR at first order in A: k(0)=1, k'(0)=-1")
    md.append("- Horizon at A=1 with proper-time divergence: k(A) ~ (1-A)^p near A=1, p>=2")
    md.append("- Pair structure (1-A^2) appearing throughout framework results")
    md.append("")
    md.append("Natural family satisfying all three: **k_n(A) = (1-A)(1-A^2)^n** for n>=1.")
    md.append("")
    md.append("## Family analysis")
    md.append("")
    md.append("| n | k_n(A) | weak-field GR | proper time div | tau_n/tau_GR (QNM) | NEC crossover |")
    md.append("|---|---|---|---|---|---|")
    for r in results:
        n = r["n"]
        nec_str = ", ".join(f"{c:.4f}" for c in r["NEC_crossings"]) if r["NEC_crossings"] else "none"
        diverges_str = "YES" if r["diverges"] else "no"
        md.append(
            f"| {n} | (1-A)(1-A^2)^{n} | {r['k0']}, {r['kp0']} | {diverges_str} | "
            f"{r['qnm_float']:.4f} | {nec_str} |"
        )
    md.append("")
    md.append("**All n >= 1 satisfy the framework's stated constraints.** Model-A's n=2 is "
              "ONE MEMBER of a 1-parameter family. The 'three principles force k(A)' claim "
              "in PDF Section 4 is overstated -- three principles narrow the form to this "
              "family but don't uniquely specify n=2.")
    md.append("")
    md.append("## QNM eikonal damping ratio")
    md.append("")
    md.append("For static spherical metric with the framework's g_tt and g_rr = 1/k_n:")
    md.append("")
    md.append("```")
    md.append("tau_n / tau_GR = (9/5)^(n/2)")
    md.append("```")
    md.append("")
    md.append("This is a concrete framework prediction depending on n. Specific values:")
    md.append("")
    md.append("- n=0 (Schwarzschild): tau/tau_GR = 1 (no departure)")
    md.append("- n=1: tau/tau_GR ≈ 1.342")
    md.append("- n=2 (Model-A): tau/tau_GR = 9/5 = 1.800 -- PDF's G1 prediction")
    md.append("- n=3: tau/tau_GR ≈ 2.415")
    md.append("- n=4: tau/tau_GR = 81/25 = 3.240")
    md.append("- n=5: tau/tau_GR ≈ 4.350")
    md.append("")
    md.append("LIGO current precision is tau_obs/tau_GR ~ 1.0 +/- 0.2 (mass-independent).")
    md.append("This nominally excludes n>=1 in the eikonal approximation. However:")
    md.append("- LIGO BHs are spinning (Kerr-like); Model-A spinless analog is not yet derived")
    md.append("- Exact (non-eikonal) Regge-Wheeler computation has not been done for the family")
    md.append("- So this comparison is indicative, not falsifying")
    md.append("")
    md.append("## What picks out n=2 (current Model-A commitment)?")
    md.append("")
    md.append("Heuristic structural arguments for n=2 (NOT derivations):")
    md.append("- Two-face refinement: 'two faces x pair structure' -> (1-A^2)^2")
    md.append("- G18 ledger entropy decomposition: (2x2)=4 area-per-entry "
              "  (but G18 assumed n=2 to start, so this is potentially circular)")
    md.append("- Bekenstein-Hawking S = A/(4*L_P^2) is reproduced under n=2 via thermodynamic argument")
    md.append("- (4π × 3) decomposition of A_0 has integer 3 but doesn't directly fix n")
    md.append("")
    md.append("None of these UNIQUELY force n=2 vs n=1 or n=3. The framework's commitment to n=2 "
              "is currently structural-aesthetic, not derived.")
    md.append("")
    md.append("## What we learn from the miss")
    md.append("")
    md.append(
        "**1. The Lagrangian-gap remains.** Inverting the action principle doesn't close "
        "it -- it just reveals that the framework's constraints admit a 1-parameter family "
        "of k_n forms. The Lagrangian (if it exists) would need to fix n through some "
        "additional structural principle we haven't yet identified."
    )
    md.append("")
    md.append(
        "**2. Model-A's specific n=2 needs a fixing point.** Two candidates:"
    )
    md.append("- (a) Observational: a clean measurement of tau/tau_GR for a spinless BH ringdown")
    md.append("- (b) Theoretical: an additional structural principle that uniquely picks n=2")
    md.append("")
    md.append(
        "**3. Pair structure (1-A^2) IS load-bearing.** It appears in all viable k_n forms; "
        "n=0 (no pair factor) gives pure GR, no STAM departure. The number of pair factors n "
        "is the open question."
    )
    md.append("")
    md.append(
        "**4. The (9/5)^(n/2) scaling is a concrete framework prediction.** Different n give "
        "measurably different QNM signatures. This is the framework's most direct exposure "
        "point to falsification once the spinless-analog issue is sorted out."
    )
    md.append("")
    md.append("## Where this points us next")
    md.append("")
    md.append("- **Push G1 toward exact (non-eikonal) QNM** for the family {k_n}. The exact "
              "spectrum may discriminate more sharply than the eikonal estimate.")
    md.append("- **Push G18 toward area-per-entry derivation WITHOUT assuming n=2.** If a "
              "first-principles ledger count gives 4 area/entry, n=2 is forced; if it gives a "
              "different number, n is what that count picks out.")
    md.append("- **Update PDF v0.4** to reflect that k(A) is a 1-parameter family, not a "
              "uniquely-forced form. Model-A's n=2 is a working commitment within the family.")
    md.append("")
    md.append("## Plot")
    md.append("")
    md.append("![G33 k family](../plots/G33_k_family.png)")
    md.append("")
    md.append("Left: k_n(A) shapes for n in {0, 1, 2, 3, 4, 5}. Model-A's n=2 is bold blue. "
              "All n>=1 satisfy the framework's stated constraints; n=2 is one choice among them.")
    md.append("")
    md.append("Right: QNM eikonal damping ratio vs n. The (9/5)^(n/2) scaling parameterizes "
              "a continuum of framework predictions. Model-A's n=2 commitment gives the "
              "1.80 value cited in the PDF.")
    md.append("")
    md.append("## Reading")
    md.append("")
    md.append(
        "**The swing missed in an informative way.** We didn't derive k(A) = (1-A)(1-A^2)^2 "
        "from first principles -- we discovered that the framework's constraints admit a "
        "1-parameter family, and Model-A's commitment to n=2 is structural-aesthetic rather "
        "than forced. The miss tells us where the actual fixing point would have to come from: "
        "either an exact QNM observation, or a non-circular derivation of the entropy decomposition."
    )
    md.append("")
    md.append(
        "This is a sharpening of the strong-field-departure question (see "
        "project_strong_field_departure_question.md). The question is no longer 'what makes k(A) "
        "concrete?' but 'what picks n=2 out of the family?'."
    )

    md_path = RESULTS / "G33_inverted_action_k_family_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Summary saved: {md_path}")


if __name__ == "__main__":
    main()
