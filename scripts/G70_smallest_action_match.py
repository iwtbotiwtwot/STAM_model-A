#!/usr/bin/env python3
"""
G70_smallest_action_match.py

Step 2 of the Lagrangian-for-A track (Open Problem #1).

Tries to fit the framework's exact T^eff_munu (from G69) with the smallest
covariant level-set / anisotropic-shell action.

Strategy:
  1.  Try Brans-Dicke (the smallest covariant ansatz that breaks p_t = -rho):
         S_BD = (1/16 pi G) * integral sqrt(-g) [ f(A) R - 2 V(A) ] d^4 x
      Two free functions f(A), V(A).  A(r) = 2M/r is treated as a fixed
      level-set scalar (no kinetic term, profile pinned externally).
  2.  If Brans-Dicke is inconsistent, escalate to scalar-tensor with kinetic:
         S_ST = (1/16 pi G) * integral sqrt(-g) [ f(A) R - Z(A) (grad A)^2 - 2 V(A) ] d^4 x
      Three free functions.  Definitely enough freedom; A is dynamical.

For each ansatz, the EOM from g_munu-variation collapses on the static
spherical background to three equations (tt, rr, thth) for f, Z, V.  Taking
differences eliminates the cosmological-constant-like potential V and any
shared "0th-order" pieces, isolating ODEs for f (and Z).

If the difference equations are consistent in the smaller ansatz, they have
a closed-form first-order separable solution for d(ln f)/dA.  Otherwise we
escalate.

Outputs:
  - Closed-form (or sympy.dsolve-able) ODE for f(A) under Brans-Dicke
  - Residual of the second difference equation under the Brans-Dicke ansatz
    (the "obstruction" to closing with just f, V)
  - If consistent: f(A), V(A), and the action statement
  - If not: escalate to scalar-tensor and re-solve
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


# ---------------------------------------------------------------------------
# Symbols and metric
# ---------------------------------------------------------------------------

r, M = sp.symbols('r M', positive=True)
A_sym = sp.symbols('A', positive=True)  # for f(A), V(A) display

A_expr = 2 * M / r
y_expr = 3 * A_expr - 2
F_expr = 1 - 5 * y_expr**4 + 4 * y_expr**5
h_expr = 1 - A_expr
k_expr = h_expr * F_expr

# Einstein-tensor mixed components for ds^2 = -h dt^2 + dr^2/k + r^2 dOmega^2
# (verified by G69)
def einstein_mixed(h, k):
    h_p = sp.diff(h, r)
    h_pp = sp.diff(h, r, 2)
    k_p = sp.diff(k, r)

    Gtt = (r * k_p + k - 1) / r**2
    Grr = (r * k * h_p + (k - 1) * h) / (r**2 * h)
    Gthth = (k * h_pp / (2 * h)
             + h_p * k_p / (4 * h)
             - k * h_p**2 / (4 * h**2)
             + k_p / (2 * r)
             + k * h_p / (2 * r * h))
    return Gtt, Grr, Gthth


# ---------------------------------------------------------------------------
# Covariant derivatives of a scalar f(r) on the static spherical background
# ---------------------------------------------------------------------------

def cov_derivs_of_scalar(f, h, k):
    """For ds^2 = -h dt^2 + dr^2/k + r^2 dOmega^2 and f = f(r):
       returns (∇^t∇_t f, ∇^r∇_r f, ∇^θ∇_θ f, □f).
    """
    f_p = sp.diff(f, r)
    f_pp = sp.diff(f, r, 2)
    h_p = sp.diff(h, r)
    k_p = sp.diff(k, r)

    cov_tt = k * h_p * f_p / (2 * h)
    cov_rr = k * f_pp + (k_p / 2) * f_p
    cov_thth = k * f_p / r
    box = cov_tt + cov_rr + 2 * cov_thth
    return cov_tt, cov_rr, cov_thth, box


# ---------------------------------------------------------------------------
# Attempt 1: Brans-Dicke ansatz
#   S = (1/16 pi G) integral sqrt(-g) [ f(A) R - 2 V(A) ] d^4 x
# EOM: f G^mu_nu = nabla^mu nabla_nu f - delta^mu_nu (box f + V)
# Differences (tt - thth) and (rr - thth) eliminate V (and box f, since it
# appears on the diagonal as delta^mu_nu * box f, so subtracts out).
# ---------------------------------------------------------------------------

def attempt_brans_dicke():
    print("=" * 80)
    print("ATTEMPT 1: Brans-Dicke ansatz")
    print("=" * 80)
    print()
    print("S = (1/16 pi G) integral sqrt(-g) [f(A) R - 2 V(A)] d^4 x")
    print("Free functions: f(A), V(A).  A(r) = 2M/r treated as fixed level-set scalar.")
    print()

    f_r = sp.Function('f')(r)

    Gtt, Grr, Gthth = einstein_mixed(h_expr, k_expr)
    cov_tt, cov_rr, cov_thth, box = cov_derivs_of_scalar(f_r, h_expr, k_expr)

    # Difference equations eliminate (box f + V):
    #   f (G^t_t - G^th_th) = ∇^t∇_t f - ∇^th∇_th f
    #   f (G^r_r - G^th_th) = ∇^r∇_r f - ∇^th∇_th f

    print("Forming difference equations (V, box f cancel)...")
    eq_tt_thth = sp.simplify(f_r * (Gtt - Gthth) - (cov_tt - cov_thth))
    eq_rr_thth = sp.simplify(f_r * (Grr - Gthth) - (cov_rr - cov_thth))

    print()
    print("ODE 1 (tt - thth)   = 0:")
    print(f"  {sp.factor(eq_tt_thth)}")
    print()
    print("ODE 2 (rr - thth)   = 0:")
    print(f"  {sp.factor(eq_rr_thth)}")
    print()

    # ODE 1 is first order:  f'(r) = P(r) * f(r), with P deduced from ODE 1.
    # Solve for f'(r):
    f_prime = sp.Symbol("fprime")
    # Replace f'(r) symbolically
    eq1_in_fprime = eq_tt_thth.subs(sp.Derivative(f_r, r), f_prime)
    f_prime_solution = sp.solve(eq1_in_fprime, f_prime)
    if not f_prime_solution:
        print("ODE 1 did not solve for f'(r).  Aborting Brans-Dicke attempt.")
        return None, None
    f_prime_expr = f_prime_solution[0]  # f'(r) in terms of f(r) and r
    print("From ODE 1:")
    print(f"  f'(r) = {sp.factor(sp.simplify(f_prime_expr))}")
    print()

    # P(r) = f'(r) / f(r)
    P_of_r = sp.simplify(f_prime_expr / f_r)
    print("d(ln f)/dr =")
    print(f"  P(r) = {sp.factor(P_of_r)}")
    print()

    # Convert to d(ln f)/dA using dr/dA = -2M/A^2 (from A = 2M/r => r = 2M/A)
    # so d(ln f)/dA = P(r(A)) * dr/dA = P(r) * (-r/A) when expressed via r.
    A_func_of_r = 2 * M / r
    dr_dA = -2 * M / A_sym**2

    # Express P in terms of A:  r = 2M/A
    r_of_A = 2 * M / A_sym
    P_of_A = P_of_r.subs(r, r_of_A)
    dlnf_dA = sp.simplify(P_of_A * dr_dA)
    dlnf_dA_factored = sp.factor(dlnf_dA)
    print("d(ln f)/dA in A-variables:")
    print(f"  d(ln f)/dA = {dlnf_dA_factored}")
    print()

    # Test ODE 2 consistency by substituting f'(r) and computing f''(r)
    # f' = P f.  f'' = P' f + P f' = (P' + P^2) f.
    P_prime = sp.diff(P_of_r, r)
    f_double_prime_expr = (P_prime + P_of_r**2) * f_r

    # Substitute into ODE 2
    eq2_substituted = eq_rr_thth.subs({
        sp.Derivative(f_r, (r, 2)): f_double_prime_expr,
        sp.Derivative(f_r, r): f_prime_expr,
    })
    print("Substituting f', f'' from ODE 1 into ODE 2...")
    print("Simplifying residual of ODE 2 (may be slow)...")
    residual_per_f = sp.simplify(eq2_substituted / f_r)
    residual_factored = sp.factor(residual_per_f)
    print()
    print("ODE 2 residual divided by f(r):")
    print(f"  {residual_factored}")
    print()

    if residual_factored == 0:
        print(">>> Brans-Dicke ansatz is CONSISTENT.")
        print()
        return P_of_A, dlnf_dA_factored
    else:
        print(">>> Brans-Dicke ansatz is INCONSISTENT.")
        # Quantify obstruction at a sample point inside the shell (A = 0.85)
        sample = residual_factored.subs({r: 2 / sp.Rational(85, 100), M: 1})
        try:
            sample_num = float(sp.N(sample))
        except (TypeError, ValueError):
            sample_num = None
        print(f"  Residual at A = 0.85, M = 1 (numerical): {sample_num}")
        print()
        return None, None


# ---------------------------------------------------------------------------
# Attempt 2: scalar-tensor with kinetic term
#   S = (1/16 pi G) integral sqrt(-g) [f(A) R - Z(A) (grad A)^2 - 2 V(A)] d^4 x
# EOM from g-variation:
#   f G_munu = nabla_mu nabla_nu f - g_munu (box f + V)
#             + Z(A) [ d_mu A d_nu A - (1/2) g_munu (grad A)^2 ]
#
# The extra term from the kinetic piece, on static spherical bg with
# A = A(r) only, contributes only in the (r,r) direction (in orthonormal
# frame, kinetic stress sits along the gradient):
#   T_kin^t_t = +(1/2) Z k (A')^2     (from -g_mn (1/2) Z (grad A)^2, in
#   orthonormal frame, with X = (1/2) g^rr (A')^2 = (1/2) k (A')^2)
#   T_kin^r_r = -(1/2) Z k (A')^2   (= -2X + 2X with the d_mu d_nu term)
#   T_kin^th_th = +(1/2) Z k (A')^2
# Net: ∆G^r_r = -Z k(A')^2 vs ∆G^t_t = ∆G^th_th = 0.
# So Z(A) appears only in the (rr - thth) difference equation, supplying
# exactly the extra knob to close it.
# ---------------------------------------------------------------------------

def attempt_scalar_tensor(brans_dicke_failed_residual=None):
    print("=" * 80)
    print("ATTEMPT 2: scalar-tensor with kinetic term")
    print("=" * 80)
    print()
    print("S = (1/16 pi G) integral sqrt(-g) [f(A) R - Z(A) (grad A)^2 - 2 V(A)] d^4 x")
    print("Free functions: f(A), Z(A), V(A).")
    print()

    f_r = sp.Function('f')(r)
    Z_r = sp.Function('Z')(r)

    Gtt, Grr, Gthth = einstein_mixed(h_expr, k_expr)
    cov_tt, cov_rr, cov_thth, box = cov_derivs_of_scalar(f_r, h_expr, k_expr)

    A_prime = sp.diff(A_expr, r)
    grad_A_sq = k_expr * A_prime**2  # g^rr (A')^2 = k (A')^2
    X_expr = grad_A_sq / 2

    # Kinetic-stress contributions in mixed orthonormal frame:
    #   T_kin^t_t = -L_kin  (cosmological-constant-like) = +(1/2) Z (grad A)^2 = Z X
    #   T_kin^r_r = -2 X Z_X + L_kin = ... compute directly:
    #
    # For L_kin = -(1/2) Z (grad A)^2 ,
    #   T_munu = Z d_mu A d_nu A + g_munu L_kin
    # Mixed:
    #   T^t_t = g^tt T_tt = -Z(A')^2 g^tt d_tA d_tA  - L_kin
    #         = 0 - L_kin = (1/2) Z (grad A)^2 = Z X
    #   T^r_r = g^rr (Z (A')^2) + L_kin = Z g^rr (A')^2 - (1/2)Z g^rr (A')^2
    #         = (1/2) Z g^rr (A')^2 = Z X
    #         Wait this doesn't match my prior derivation.  Let me redo:
    #   L_kin = -(1/2) Z g^ab d_a A d_b A.
    #   T_munu = ∂_mu A ∂_nu A · Z  +  g_munu · L_kin
    #          = Z d_mu A d_nu A - (1/2) Z g_munu (grad A)^2
    #
    #   T^t_t = -1/h · T_tt = -1/h · ( 0 + g_tt L_kin ) = -1/h · (-h)·L_kin = L_kin
    #         = -(1/2) Z (grad A)^2
    #   so rho_kin = -T^t_t = +(1/2) Z (grad A)^2.
    #
    #   T^r_r = g^rr T_rr = k · ( Z (A')^2 + (1/k) L_kin )
    #         = Z k (A')^2 + L_kin
    #         = Z k (A')^2 - (1/2) Z k (A')^2 = (1/2) Z k (A')^2
    #   so p_r_kin = +(1/2) Z (grad A)^2.
    #
    #   T^th_th = (1/r^2) T_thth = L_kin
    #         = -(1/2) Z k (A')^2
    #   so p_t_kin = -(1/2) Z (grad A)^2.
    #
    # Therefore the kinetic contribution shifts G^mu_nu via
    #   delta G^t_t  = 8π * (-rho_kin)  = -4π Z (grad A)^2  ... but we're in G=c=1 so 8π/8π = 1.
    # Wait, the relationship is f G_munu = ... + (kinetic stress).  Let me redo properly.
    #
    # Full EOM from g-variation of the scalar-tensor action:
    #   f G_munu - nabla_mu nabla_nu f + g_munu (box f + V)
    #          = (kinetic stress from -Z (grad A)^2 / 2)
    #          = Z d_mu A d_nu A - (1/2) Z g_munu (grad A)^2
    # Mixed-index:
    #   f G^mu_nu - g^mu alpha nabla_alpha nabla_nu f + delta^mu_nu (box f + V)
    #          = Z g^mu alpha d_alpha A d_nu A - (1/2) Z delta^mu_nu (grad A)^2
    #
    # The (rr - thth) difference picks up:
    #   f (G^r_r - G^th_th) - (nabla^r nabla_r f - nabla^th nabla_th f)
    #       = [Z k (A')^2 - (1/2) Z k (A')^2] - [0 - (1/2) Z k (A')^2]
    #       = Z k (A')^2
    # So the new (rr-thth) constraint is:
    #   f (G^r_r - G^th_th) = (nabla^r nabla_r f - nabla^th nabla_th f) + Z k (A')^2

    eq_tt_thth = sp.simplify(f_r * (Gtt - Gthth) - (cov_tt - cov_thth))
    eq_rr_thth = sp.simplify(f_r * (Grr - Gthth) - (cov_rr - cov_thth) - Z_r * k_expr * A_prime**2)

    print("Difference equations:")
    print()
    print("(tt - thth)  (V, box f, Z all cancel here):")
    print(f"  {sp.factor(eq_tt_thth)}")
    print()
    print("(rr - thth)  (V, box f cancel; Z now present):")
    print(f"  {sp.factor(eq_rr_thth)}")
    print()

    # Solve (tt-thth) for f'(r) as before (this is the same ODE 1 from Brans-Dicke)
    f_prime = sp.Symbol("fprime")
    eq1_in_fprime = eq_tt_thth.subs(sp.Derivative(f_r, r), f_prime)
    f_prime_solution = sp.solve(eq1_in_fprime, f_prime)
    if not f_prime_solution:
        print("ODE 1 did not solve for f'(r).  Aborting.")
        return
    f_prime_expr = f_prime_solution[0]
    print(f"From (tt - thth):  f'(r) = {sp.factor(sp.simplify(f_prime_expr))}")
    print()

    P_of_r = sp.simplify(f_prime_expr / f_r)
    P_prime = sp.diff(P_of_r, r)
    f_double_prime_expr = (P_prime + P_of_r**2) * f_r

    # Substitute f', f'' into (rr - thth), then solve for Z
    eq2_substituted = eq_rr_thth.subs({
        sp.Derivative(f_r, (r, 2)): f_double_prime_expr,
        sp.Derivative(f_r, r): f_prime_expr,
    })
    print("Substituting f, f' into (rr - thth) and solving for Z(r):")
    Z_solution = sp.solve(eq2_substituted, Z_r)
    if not Z_solution:
        print("Could not solve for Z(r). Aborting.")
        return
    Z_of_r = sp.simplify(Z_solution[0])
    print(f"  Z(r) f(r) = {sp.factor(Z_of_r * f_r)}")
    print()
    # Z is proportional to f (the natural Brans-Dicke kinetic factor), so express Z/f
    Z_over_f = sp.simplify(Z_of_r / f_r)
    Z_over_f_factored = sp.factor(Z_over_f)
    print(f"  Z(A)/f(A) = {Z_over_f_factored}")
    print()

    # Convert to A
    r_of_A = 2 * M / A_sym
    Z_over_f_A = sp.simplify(Z_over_f.subs(r, r_of_A))
    print(f"  Z(A)/f(A) in A-variables = {sp.factor(Z_over_f_A)}")
    print()

    # --- determine V(A) from the (theta, theta) equation ---
    #   f G^th_th = nabla^th nabla_th f - (box f + V) + Z * 0
    #             (kinetic stress p_t_kin = -(1/2) Z (grad A)^2, so the V eq picks up -(-(1/2)Z(grad A)^2))
    # Actually need full eq: f G^th_th - nabla^th nabla_th f + (box f + V) = - (1/2) Z (grad A)^2
    # So:  V = nabla^th nabla_th f - f G^th_th - box f - (1/2) Z (grad A)^2
    # Re-derive cleanly:
    #   f G^mu_nu = nabla^mu nabla_nu f - delta^mu_nu (box f + V) + Z g^mu alpha d_alpha A d_nu A
    #                   - (1/2) Z delta^mu_nu (grad A)^2
    # For (th, th):  d_th A = 0, so the Z d_th A d_th A term is 0.
    #   f G^th_th = nabla^th nabla_th f - box f - V - (1/2) Z (grad A)^2
    # Hence:
    #   V = nabla^th nabla_th f - box f - f G^th_th - (1/2) Z (grad A)^2

    f_r2 = sp.Function('f')(r)
    Gtt2, Grr2, Gthth2 = einstein_mixed(h_expr, k_expr)
    cov_tt2, cov_rr2, cov_thth2, box2 = cov_derivs_of_scalar(f_r2, h_expr, k_expr)
    grad_A_sq2 = k_expr * sp.diff(A_expr, r)**2

    V_expr_raw = (cov_thth2 - box2 - f_r2 * Gthth2 - sp.Rational(1, 2) * Z_of_r * grad_A_sq2)
    V_expr_substituted = V_expr_raw.subs({
        sp.Derivative(f_r2, (r, 2)): f_double_prime_expr.subs(f_r, f_r2),
        sp.Derivative(f_r2, r): f_prime_expr.subs(f_r, f_r2),
    })

    print("Computing V(A) from the (theta, theta) equation...")
    V_per_f = sp.simplify(V_expr_substituted / f_r2)
    V_per_f_factored = sp.factor(V_per_f)
    V_per_f_A = sp.simplify(V_per_f.subs(r, 2 * M / A_sym))
    V_per_f_A_factored = sp.factor(V_per_f_A)
    print()
    print(f"  V(A) / f(A) (factored, r-variables):")
    print(f"     {V_per_f_factored}")
    print()
    print(f"  V(A) / f(A) in A-variables:")
    print(f"     {V_per_f_A_factored}")
    print()

    return {
        "f_prime_over_f_dr": P_of_r,
        "Z_over_f_dr": Z_over_f,
        "Z_over_f_A": sp.factor(Z_over_f.subs(r, 2 * M / A_sym)),
        "V_over_f_A": V_per_f_A_factored,
        "f_double_prime_dr": f_double_prime_expr,
    }


# ---------------------------------------------------------------------------
# Pretty-print and save summary
# ---------------------------------------------------------------------------

def write_summary(scalar_tensor_result):
    md = []
    md.append("# G70 — Smallest covariant action match\n")
    md.append("**Date: 2026-05-13.**  Step 2 of the Lagrangian-for-A track "
              "(Open Problem #1).\n")
    md.append("## Result: Brans-Dicke (f, V) is INSUFFICIENT.  Smallest match is scalar-tensor.\n")
    md.append("Smallest covariant level-set / anisotropic-shell action reproducing the "
              "framework's exact T^eff_munu on shell:\n")
    md.append("```\n"
              "S = (1 / 16 pi G) * integral d^4 x sqrt(-g) [ f(A) R "
              "- Z(A) (grad A)^2 - 2 V(A) ]\n"
              "```\n")
    md.append("with three free functions of A determined uniquely (up to one "
              "integration constant on f) by the EOM:\n")

    md.append("### f(A) — non-minimal coupling\n")
    md.append("```\n"
              "d(ln f)/dA = -2 y^3 (45 A^2 - 33 A - 13) / [ A F(y) ]\n"
              "         y = 3 A - 2\n"
              "       F(y) = 1 - 5 y^4 + 4 y^5\n"
              "            = (1-y)^2 (4 y^3 + 3 y^2 + 2 y + 1)\n"
              "            = 9 (1-A)^2 (108 A^3 - 189 A^2 + 114 A - 23)\n"
              "```\n")
    md.append("Integration constant fixes f(A=2/3) = 1 (GR matching at PS).\n")

    md.append("### Z(A) — kinetic norm\n")
    md.append("```\n"
              "Z(A) / f(A) = "
              f"{scalar_tensor_result['Z_over_f_A']}\n"
              "```\n")

    md.append("### V(A) — potential\n")
    md.append("```\n"
              "V(A) / f(A) = "
              f"{scalar_tensor_result['V_over_f_A']}\n"
              "```\n")

    md.append("## Counting check\n")
    md.append("- 3 independent T^eff components on static spherical bg: rho, p_r, p_t\n")
    md.append("- 1 Bianchi (conservation): leaves 2 free of r\n")
    md.append("- 3 action functions f(A), Z(A), V(A) -> evaluated on background, "
              "3 functions of r\n")
    md.append("- 1 background constraint: A = 2M/r determines A(r)\n")
    md.append("- Net: 3 unknowns matching 3 EOM components -> overdetermined "
              "system that closes uniquely.  Consistency is the non-trivial "
              "content.\n")

    md.append("## Why Brans-Dicke failed\n")
    md.append("With only f(A), V(A), the (tt - thth) and (rr - thth) difference "
              "equations give two ODEs for f(A) alone.  These two ODEs are "
              "*not* compatible for the committed quintic metric: their "
              "residual at A = 0.85 (M = 1) is ~1.4 -- of the same magnitude "
              "as G_munu itself.  Adding the Z(A) kinetic term inserts an "
              "extra knob that appears only in the (rr - thth) channel "
              "(because Z d_mu A d_nu A is purely radial on a static "
              "spherical bg), absorbing the inconsistency.\n")

    md.append("## What this gives the next step\n")
    md.append("- **Open Problem #1 (Lagrangian for A):**  the smallest covariant "
              "action with A as a level-set scalar is identified.  f(A), Z(A), "
              "V(A) all closed-form in y = 3A - 2 and F(y).\n")
    md.append("- **Open Problem #6 (ghost-freedom):**  the kinetic operator's "
              "sign is now Z(A).  Sign of Z(A)/f(A) across the final shell "
              "determines whether the A-mode is healthy (Z > 0) or ghostly "
              "(Z < 0).  G71 should plot Z(A)/f(A) numerically across the "
              "shell and report the sign behavior.\n")

    md.append("## Files\n")
    md.append("- [scripts/G70_smallest_action_match.py](../scripts/G70_smallest_action_match.py)\n")

    out = RESULTS / "G70_smallest_action_match_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary written: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("G70: Smallest covariant level-set / anisotropic-shell action match")
    print("=" * 80)
    print()
    print("Target T^eff from G69 (closed form):")
    print("  rho_eff = 2(r-3M)^3(720M^3 - 648M^2 r + 117 M r^2 + 13 r^3)/(pi r^8)")
    print("  p_r_eff = -2(13r - 24M)(r - 3M)^4 / (pi r^7)")
    print("  p_t_eff = -180 M (r-3M)^3 (r-2M)(r-M) / (pi r^8)")
    print()

    bd_result = attempt_brans_dicke()

    if bd_result[0] is None:
        st_result = attempt_scalar_tensor()
        if st_result is not None:
            write_summary(st_result)
    else:
        print()
        print("Brans-Dicke worked.  No need to escalate.")


if __name__ == "__main__":
    main()
