#!/usr/bin/env python3
"""
G113a_leaver_continued_fraction_schwarzschild.py

Leaver continued-fraction QNM solver for Schwarzschild axial (s = 2).

PURPOSE: calibration only.  Reproduce three published Schwarzschild Leaver
roots before anything else:

    (ell, n) =  (2, 0)   omega = 0.373672 - 0.088962 i
    (ell, n) =  (2, 1)   omega = 0.346711 - 0.273915 i
    (ell, n) =  (3, 0)   omega = 0.599443 - 0.092703 i

If the implementation reproduces these to <1e-4 absolute, we have a trusted
frequency-domain method.  G114 will then apply the same machinery (with the
recurrence adapted to V_exact) to the STAM problem.

Method
======
Leaver 1985 / Berti+Cardoso+Will 2009 review.  In M = 1 units (r_h = 2):

The Regge-Wheeler equation with QNM boundary conditions admits a series
solution

    psi(r) = (r-2)^{-2 i omega} * r^{2 i omega} * exp(i omega r)
             * sum_{n=0}^infty a_n  ((r-2)/r)^n

that satisfies horizon-ingoing and infinity-outgoing BCs simultaneously
when omega is a QNM.  The coefficients obey a three-term recurrence

    alpha_n a_{n+1} + beta_n a_n + gamma_n a_{n-1} = 0,    a_{-1} = 0

with (M = 1, axial s = 2):

    alpha_n = n^2 + (3 - 4 i omega) n + 2 - 4 i omega
    beta_n  = -2 n^2 + (-8 + 8 i omega) n - 8 omega^2 - ell(ell+1) + 4 i omega + 1
    gamma_n = n^2 + (1 - 4 i omega) n - 4 omega^2 - 4 i omega - 3

The QNM is the omega such that the continued fraction equation

    D_n(omega) = 0

with the "inverted" form for the n-th overtone:

    D_n(omega) = beta_n  -  alpha_{n-1} gamma_n / CF_down
                          -  alpha_n gamma_{n+1} / CF_up

    CF_up   = beta_{n+1} - alpha_{n+1} gamma_{n+2} /
                            (beta_{n+2} - alpha_{n+2} gamma_{n+3} / ...)
    CF_down = beta_{n-1} - alpha_{n-2} gamma_{n-1} /
                            (beta_{n-2} - alpha_{n-3} gamma_{n-2} / ...)
                            ... down to beta_0.

For n = 0, only CF_up is needed and D_0 = beta_0 - alpha_0 gamma_1 / CF_up.

Newton iteration finds omega.  Initial guesses use the Leaver-gold values.

Coefficient sign-convention check
=================================
This script first verifies that |D_n(omega_Leaver_gold)| is small for each
target mode (a self-check of the recurrence coefficients).  If that check
fails, the message is printed clearly and the script aborts BEFORE Newton
iteration -- so we know the recurrence (not the iteration) is the problem.

Output
======
    results/G113a_leaver_schwarzschild_summary.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (2, 1): complex(0.346711, -0.273915),
    (3, 0): complex(0.599443, -0.092703),
}


# ===========================================================================
# Leaver recurrence coefficients   (Schwarzschild axial s = 2,  M = 1)
# ===========================================================================

def leaver_abg(n, omega, ell, s=2):
    """Three-term recurrence coefficients alpha_n, beta_n, gamma_n.

    Convention: M = 1 (horizon at r = 2).  Axial RW potential.
    Reference: Berti+Cardoso+Will 2009 review, Section 4.2.
    """
    iw = 1j * omega
    alpha = n*n + (3.0 - 4.0*iw) * n + (2.0 - 4.0*iw)
    beta  = -2.0*n*n + (-8.0 + 8.0*iw) * n - 8.0*omega**2 - ell*(ell + 1) + 4.0*iw + 1.0
    gamma = n*n + (1.0 - 4.0*iw) * n - 4.0*omega**2 - 4.0*iw + (1.0 - s*s)
    return alpha, beta, gamma


# ===========================================================================
# Continued fraction evaluation
# ===========================================================================

def CF_up(omega, ell, n_start, N_max=300, s=2):
    """Evaluate CF_up = beta_{n_start} - alpha_{n_start} gamma_{n_start+1} / (...)
       down to N_max.  Returns the value of the truncated CF at index n_start.
    """
    a_N, b_N, g_N = leaver_abg(N_max, omega, ell, s)
    cf = b_N
    for k in range(N_max - 1, n_start - 1, -1):
        a_k, b_k, g_k = leaver_abg(k, omega, ell, s)
        a_kp1, b_kp1, g_kp1 = leaver_abg(k + 1, omega, ell, s)
        cf = b_k - a_k * g_kp1 / cf
    return cf


def CF_down(omega, ell, n_start, s=2):
    """Evaluate the 'inverted' CF starting at index 0 and going up to n_start - 1.
       Returns the value at index n_start - 1.
    """
    if n_start <= 0:
        return None
    a_0, b_0, g_0 = leaver_abg(0, omega, ell, s)
    cf = b_0
    for k in range(1, n_start):
        a_km1, b_km1, g_km1 = leaver_abg(k - 1, omega, ell, s)
        a_k, b_k, g_k = leaver_abg(k, omega, ell, s)
        cf = b_k - a_km1 * g_k / cf
    return cf


def D_n(omega, ell, n_overtone, N_max=300, s=2):
    """Inverted Leaver equation at overtone index n.  QNM => D_n(omega) = 0."""
    a_n, b_n, g_n = leaver_abg(n_overtone, omega, ell, s)
    a_np1, b_np1, g_np1 = leaver_abg(n_overtone + 1, omega, ell, s)
    cf_up_at_np1 = CF_up(omega, ell, n_overtone + 1, N_max, s)
    forward = a_n * g_np1 / cf_up_at_np1
    if n_overtone == 0:
        return b_n - forward
    a_nm1, b_nm1, g_nm1 = leaver_abg(n_overtone - 1, omega, ell, s)
    cf_down_at_nm1 = CF_down(omega, ell, n_overtone, s)
    backward = a_nm1 * g_n / cf_down_at_nm1
    return b_n - backward - forward


# ===========================================================================
# Newton iteration
# ===========================================================================

def newton_complex(f, x0, tol=1e-10, max_iter=80, h=1e-8):
    """Newton on a complex function with finite-difference derivative."""
    x = x0
    for it in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x, it, fx
        # Complex derivative via finite difference (real and imaginary stepping)
        f_dx = f(x + h)
        derivative = (f_dx - fx) / h
        if derivative == 0:
            return None, it, fx
        dx = -fx / derivative
        # Damping if step too large
        if abs(dx) > 0.5 * abs(x):
            dx = dx * 0.5 * abs(x) / abs(dx)
        x = x + dx
    return x, max_iter, f(x)


# ===========================================================================
# Coefficient self-check
# ===========================================================================

def coefficient_self_check(N_max=400):
    """Test |D_n(omega_Leaver_gold)| for each target mode.  If the
       recurrence coefficients are correctly transcribed, these should
       all be small (~1e-6 or smaller).  If they are O(1) or larger,
       the coefficients are wrong."""
    print("=" * 88)
    print("Coefficient self-check")
    print("=" * 88)
    print()
    print(f"  {'mode':>14}  {'omega_Leaver':>30}  {'|D_n(omega)|':>14}  status")
    print("-" * 78)
    all_pass = True
    for (ell, n), om_gold in LEAVER_GOLD.items():
        D_val = D_n(om_gold, ell, n, N_max=N_max)
        mag = abs(D_val)
        # |D_n| relative to typical |beta_n| scale
        a, b, g = leaver_abg(n, om_gold, ell)
        rel = mag / max(abs(b), 1.0)
        ok = "OK" if rel < 1e-3 else "FAIL"
        if ok != "OK":
            all_pass = False
        print(f"  ({ell:>1},{n:>1})       {str(om_gold):>30}  {mag:>14.4e}  rel={rel:.2e}  {ok}")
    print()
    return all_pass


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 88)
    print("G113a -- Leaver continued-fraction QNM solver for Schwarzschild axial")
    print("=" * 88)
    print()
    print("Units: M = 1 (horizon at r = 2).")
    print("Targets:")
    for (ell, n), om in LEAVER_GOLD.items():
        print(f"   (ell={ell}, n={n})   {om}")
    print()

    # First check coefficients
    coeffs_ok = coefficient_self_check(N_max=400)
    if not coeffs_ok:
        print("  >>>  COEFFICIENT SELF-CHECK FAILED  <<<")
        print("  Leaver recurrence coefficients are inconsistent with published")
        print("  Schwarzschild QNM values.  Aborting before Newton iteration.")
        print("  Falling back to G113b numerical shooting method recommended.")
        # Still write a stub summary
        md = ["# G113a -- Leaver continued fraction (Schwarzschild)\n"]
        md.append("\n## Status: COEFFICIENT SELF-CHECK FAILED\n")
        md.append("\nThe recurrence coefficients did not reproduce known Leaver roots.\n")
        md.append("See script header for the exact coefficient forms used.\n")
        md.append("Recommended: G113b numerical shooting.\n")
        (RESULTS / "G113a_leaver_schwarzschild_summary.md").write_text("".join(md), encoding="utf-8")
        return

    # Coefficients OK: run Newton on each target
    print("=" * 88)
    print("Newton iteration on continued-fraction equation")
    print("=" * 88)
    print()
    print(f"  {'mode':>14}  {'omega_solved':>30}  {'omega_Leaver':>30}  {'error':>14}  iter")
    print("-" * 110)
    results = []
    for (ell, n), om_gold in LEAVER_GOLD.items():
        f = lambda om, ell=ell, n=n: D_n(om, ell, n, N_max=400)
        # Initial guess: slightly perturbed Leaver gold so we test actual convergence
        x0 = om_gold * (1.0 + 0.01j)
        x_solved, n_iter, fx = newton_complex(f, x0, tol=1e-9, max_iter=80)
        if x_solved is None:
            print(f"  ({ell:>1},{n:>1})       Newton FAILED")
            results.append({"ell": ell, "n": n, "omega": None, "err": None})
            continue
        err = abs(x_solved - om_gold) / abs(om_gold)
        ok = "OK" if err < 1e-4 else "FAIL"
        print(f"  ({ell:>1},{n:>1})       {str(x_solved):>30}  {str(om_gold):>30}"
              f"  {err:>14.4e}  {n_iter:>4}  {ok}")
        results.append({"ell": ell, "n": n, "omega": x_solved, "err": err,
                        "n_iter": n_iter})
    print()

    # Summary
    all_calib = all(r.get("err") is not None and r["err"] < 1e-4 for r in results)
    if all_calib:
        print("CALIBRATION PASSED for all three Schwarzschild modes.")
        print("Leaver continued fraction is a trusted frequency-domain tool.")
        print("G114 may use the same machinery (with STAM-V_exact-derived recurrence")
        print("if available, otherwise via G113b numerical shooting).")
    else:
        print("CALIBRATION FAILED for at least one mode.")
        print("Recommend G113b numerical shooting as the working tool.")
    print()

    # Markdown
    md = ["# G113a -- Leaver continued-fraction QNM solver for Schwarzschild axial\n"]
    md.append("\n**Units:** M = 1 (horizon at r = 2).\n")
    md.append("\n## Coefficient self-check\n")
    md.append("Re-evaluate |D_n(omega)| at each published Leaver root with the\n"
              "implementation's coefficients.  Should be small (<1e-3 relative) if\n"
              "the recurrence is correctly transcribed.\n\n")
    md.append("(See script output for details.)\n")
    md.append("\n## Newton iteration results\n")
    md.append("| ell | n | omega_solved | omega_Leaver | rel error | iter |\n")
    md.append("|---:|---:|---|---|---:|---:|\n")
    for r in results:
        ell = r["ell"]; n = r["n"]
        if r["omega"] is None:
            md.append(f"| {ell} | {n} | FAILED | {LEAVER_GOLD[(ell,n)]} | --- | --- |\n")
        else:
            md.append(f"| {ell} | {n} | {r['omega']} | {LEAVER_GOLD[(ell,n)]} | "
                      f"{r['err']:.4e} | {r['n_iter']} |\n")
    md.append("\n## Status\n")
    if all_calib:
        md.append("PASS  --  Leaver method calibrated.  Suitable for G114 STAM work.\n")
    else:
        md.append("FAIL  --  Use G113b numerical shooting instead.\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G113a_leaver_continued_fraction_schwarzschild.py]"
              "(../scripts/G113a_leaver_continued_fraction_schwarzschild.py)\n")
    (RESULTS / "G113a_leaver_schwarzschild_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G113a_leaver_schwarzschild_summary.md'}")


if __name__ == "__main__":
    main()
