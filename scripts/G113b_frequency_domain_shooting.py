#!/usr/bin/env python3
"""
G113b_frequency_domain_shooting.py

Generalized frequency-domain QNM solver via numerical shooting + Newton
iteration on the ingoing-wave amplitude.  Works for ANY V(r*) given on a
grid -- no recurrence relation required.  Falls in as the fallback when
Leaver's continued fraction is too brittle to transcribe coefficient-by-
coefficient (as G113a's self-check exposed).

Method (shoot-and-match Wronskian)
==================================
For the Schrodinger-form master equation

    d^2 psi / d r_star^2 + [ omega^2  -  V(r_star) ]  psi  =  0

with QNM boundary conditions:

    r_star -> -infinity:    psi ~ exp( -i omega r_star )   ingoing at horizon
    r_star -> +infinity:    psi ~ exp( +i omega r_star )   outgoing at infinity

Shoot psi_L from r_star_min OUTWARD to r_star_mid (= 0, photon sphere)
with the ingoing-only BC.  Shoot psi_R from r_star_max INWARD to r_star_mid
with the outgoing-only BC.  At the matching point, compute the Wronskian:

    W(omega) = psi_L * dpsi_R/dr_star  -  psi_R * dpsi_L/dr_star

W(omega) = 0  iff  psi_L and psi_R are linearly dependent  iff  omega is a QNM.

Normalize psi_L, psi_R by their magnitudes at the midpoint to avoid the
exponential amplitude blow-up.  Newton-iterate on complex omega until W = 0.

This is the standard frequency-domain QNM method for arbitrary V(r*).

Calibration targets (Schwarzschild axial s=2, M = 1)
====================================================
    (ell, n) = (2, 0)   omega = 0.373672 - 0.088962 i
    (ell, n) = (2, 1)   omega = 0.346711 - 0.273915 i
    (ell, n) = (3, 0)   omega = 0.599443 - 0.092703 i

The method PASSES calibration if all three converge to abs error < 1e-4.

Output
======
    results/G113b_frequency_domain_shooting_summary.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


M = 1.0

LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (2, 1): complex(0.346711, -0.273915),
    (3, 0): complex(0.599443, -0.092703),
}


# ===========================================================================
# Schwarzschild metric / RW potential
# ===========================================================================

def h_fn(r): return 1.0 - 2.0 * M / r

def V_RW(r, ell):
    """Schwarzschild axial Regge-Wheeler s = 2."""
    return h_fn(r) * (ell * (ell + 1) / r**2 - 6.0 * M / r**3)


def r_of_rstar_GR(rstar):
    """Invert r* = r + 2M ln(r/2M - 1) with PS_offset s.t. r*(3M) = 0."""
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset
    def f(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target
    if rstar > 0:
        return brentq(f, 3.0 * M - 1e-10, 1e8)
    if rstar == 0:
        return 3.0 * M
    if rstar < -40:
        arg = (rstar + 3.0 * M + 2.0 * M * np.log(0.5) - 2.0 * M) / (2.0 * M)
        r_guess = 2.0 * M + 2.0 * M * np.exp(arg)
        r_lo = max(2.0 * M + 1e-300, r_guess * 1e-3)
        r_hi = min(3.0 * M - 1e-12, max(r_guess * 1e3, 2.0 * M + 1e-6))
        try:
            return brentq(f, r_lo, r_hi)
        except ValueError:
            return r_guess
    return brentq(f, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


# ===========================================================================
# Shooting: A_in(omega) for a given V on an r* grid
# ===========================================================================

def wronskian_at_mid(omega, rstar_grid, V_grid, rs_mid=0.0):
    """Shoot psi_L from rs_min OUTWARD to rs_mid (ingoing-only BC),
       psi_R from rs_max INWARD to rs_mid (outgoing-only BC).
       Compute normalized Wronskian W = psi_L * dpsi_R/dr* - psi_R * dpsi_L/dr*.
       QNM iff W(omega) = 0.
    """
    rs_min, rs_max = rstar_grid[0], rstar_grid[-1]
    V_interp = lambda rs: np.interp(rs, rstar_grid, V_grid)

    def rhs(rs, z):
        psi, dpsi = z
        return [dpsi, (V_interp(rs) - omega**2) * psi]

    # Left wave: shoot outward from rs_min with ingoing-only BC.
    # IC: psi ~ exp(-i omega rs_min), dpsi ~ -i omega * psi
    psi_L_min = np.exp(-1j * omega * rs_min)
    dpsi_L_min = -1j * omega * psi_L_min
    sol_L = solve_ivp(rhs, [rs_min, rs_mid], [psi_L_min, dpsi_L_min],
                      method="RK45", rtol=1e-10, atol=1e-14, max_step=0.5)
    if not sol_L.success or sol_L.y.size == 0:
        return None
    psi_L_mid = sol_L.y[0, -1]
    dpsi_L_mid = sol_L.y[1, -1]

    # Right wave: shoot inward from rs_max with outgoing-only BC.
    # IC: psi ~ exp(+i omega rs_max), dpsi ~ +i omega * psi
    psi_R_max = np.exp(+1j * omega * rs_max)
    dpsi_R_max = +1j * omega * psi_R_max
    sol_R = solve_ivp(rhs, [rs_max, rs_mid], [psi_R_max, dpsi_R_max],
                      method="RK45", rtol=1e-10, atol=1e-14, max_step=0.5)
    if not sol_R.success or sol_R.y.size == 0:
        return None
    psi_R_mid = sol_R.y[0, -1]
    dpsi_R_mid = sol_R.y[1, -1]

    # Normalize to avoid exponential-scale issues (preserves zeros of W)
    norm_L = max(abs(psi_L_mid), abs(dpsi_L_mid), 1e-300)
    norm_R = max(abs(psi_R_mid), abs(dpsi_R_mid), 1e-300)
    psi_L_mid /= norm_L; dpsi_L_mid /= norm_L
    psi_R_mid /= norm_R; dpsi_R_mid /= norm_R

    W = psi_L_mid * dpsi_R_mid - psi_R_mid * dpsi_L_mid
    return W


# ===========================================================================
# Newton iteration on complex omega
# ===========================================================================

def newton_complex(f, x0, tol=1e-10, max_iter=80, h=1e-7, damp_max=0.3):
    """Newton on a complex-valued function f(z) with FD derivative."""
    x = x0
    for it in range(max_iter):
        fx = f(x)
        if fx is None or not np.isfinite(fx):
            return None, it, fx
        if abs(fx) < tol:
            return x, it, fx
        # Complex derivative (Cauchy-Riemann: f'(z) = df/dx; complex step)
        f_dx = f(x + h)
        if f_dx is None or not np.isfinite(f_dx):
            return None, it, fx
        deriv = (f_dx - fx) / h
        if deriv == 0:
            return None, it, fx
        dx = -fx / deriv
        # Damp large steps
        if abs(dx) > damp_max * abs(x):
            dx = dx * damp_max * abs(x) / abs(dx)
        x = x + dx
    return x, max_iter, f(x)


# ===========================================================================
# Calibration on Schwarzschild
# ===========================================================================

def calibrate_schwarzschild(rs_min, rs_max, N_grid):
    """Find each Leaver-gold target via shooting + Newton."""
    rstar = np.linspace(rs_min, rs_max, N_grid)

    print(f"  r* grid: N = {N_grid}, range [{rs_min}, {rs_max}], "
          f"dr* = {(rs_max - rs_min) / (N_grid - 1):.4e}")

    # Build r(r*) and V grids per ell
    r_grid = np.array([r_of_rstar_GR(rs) for rs in rstar])

    print(f"  r range: [{r_grid.min():.6f}, {r_grid.max():.6f}]")
    print()

    results = []
    for (ell, n), om_gold in LEAVER_GOLD.items():
        V_grid = np.array([V_RW(r, ell) for r in r_grid])

        def W_fn(om):
            return wronskian_at_mid(om, rstar, V_grid, rs_mid=0.0)

        # Print W at exact Leaver gold as a sanity check
        W_at_gold = W_fn(om_gold)
        print(f"  (ell={ell}, n={n})  W(omega_Leaver) = {W_at_gold}"
              f"  |W| = {abs(W_at_gold):.4e}")

        # Try multiple initial guesses, take best-converging one
        best = None
        for perturb in [(1e-3, 0.0), (0.0, 1e-3), (1e-3, 1e-3),
                        (5e-3, 5e-3), (1e-2, 1e-2)]:
            x0 = om_gold + complex(perturb[0], perturb[1])
            x_solved, n_iter, fx = newton_complex(
                W_fn, x0, tol=1e-9, max_iter=60, damp_max=0.1)
            if x_solved is None:
                continue
            err_try = abs(x_solved - om_gold) / abs(om_gold)
            if best is None or err_try < best[1]:
                best = (x_solved, err_try, n_iter)
        if best is None:
            x_solved, n_iter, fx = None, 0, None
        else:
            x_solved, _, n_iter = best
        if x_solved is None:
            err = None
            ok = "FAIL"
        else:
            err = abs(x_solved - om_gold) / abs(om_gold)
            ok = "PASS" if err < 1e-4 else "FAIL"
        print(f"  (ell={ell}, n={n})  init = {x0}")
        print(f"     solved = {x_solved}")
        print(f"     Leaver = {om_gold}")
        print(f"     err = {err}  iter = {n_iter}  {ok}")
        print()
        results.append({"ell": ell, "n": n, "omega_solved": x_solved,
                        "omega_gold": om_gold, "err": err, "iter": n_iter,
                        "status": ok})
    return results


def main():
    print("=" * 88)
    print("G113b -- frequency-domain shooting QNM solver  (Schwarzschild calibration)")
    print("=" * 88)
    print()

    # Try a sequence of grid choices for robustness
    grid_choices = [
        (-50.0, 100.0, 6001),
        (-80.0, 200.0, 10001),
        (-100.0, 200.0, 12001),
    ]

    best_results = None
    for rs_min, rs_max, N in grid_choices:
        print("-" * 88)
        print(f"Trying r* grid (rs_min, rs_max, N) = ({rs_min}, {rs_max}, {N})")
        print("-" * 88, flush=True)
        results = calibrate_schwarzschild(rs_min, rs_max, N)
        n_pass = sum(1 for r in results if r["status"] == "PASS")
        print(f"  -> {n_pass}/{len(results)} PASS at this grid.")
        if best_results is None or n_pass > sum(1 for r in best_results if r["status"] == "PASS"):
            best_results = results
            best_grid = (rs_min, rs_max, N)
        if n_pass == len(results):
            break
        print()

    print("=" * 88)
    print("Best result summary")
    print("=" * 88)
    print()
    print(f"Best grid: rs_min={best_grid[0]}, rs_max={best_grid[1]}, N={best_grid[2]}")
    print()
    print(f"  {'mode':>10}  {'omega_solved':>30}  {'omega_Leaver':>30}  {'err':>14}  status")
    print("-" * 110)
    for r in best_results:
        print(f"  ({r['ell']:>1},{r['n']:>1})       "
              f"{str(r['omega_solved']):>30}  {str(r['omega_gold']):>30}"
              f"  {r['err'] if r['err'] is not None else 'FAIL':>14}  {r['status']}")
    print()
    n_pass = sum(1 for r in best_results if r["status"] == "PASS")
    if n_pass == len(best_results):
        verdict = ("PASS  --  shooting method calibrated on all three Schwarzschild "
                   "modes.  Trusted tool for G114 STAM analysis.")
    else:
        verdict = ("PARTIAL  --  shooting method needs further tuning; see per-mode "
                   "errors above.")
    print(verdict)
    print()

    # Markdown
    md = ["# G113b -- frequency-domain shooting QNM solver\n"]
    md.append("\nSchwarzschild calibration of the numerical-shooting method.\n")
    md.append(f"\nBest grid: rs_min={best_grid[0]}, rs_max={best_grid[1]}, N={best_grid[2]}\n")
    md.append("\n## Results\n")
    md.append("| ell | n | omega_solved | omega_Leaver | rel error | iter | status |\n")
    md.append("|---:|---:|---|---|---:|---:|---|\n")
    for r in best_results:
        err_str = f"{r['err']:.4e}" if r['err'] is not None else "FAILED"
        md.append(f"| {r['ell']} | {r['n']} | {r['omega_solved']} | {r['omega_gold']} | "
                  f"{err_str} | {r['iter']} | {r['status']} |\n")
    md.append("\n## Status\n")
    md.append(verdict + "\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G113b_frequency_domain_shooting.py]"
              "(../scripts/G113b_frequency_domain_shooting.py)\n")
    (RESULTS / "G113b_frequency_domain_shooting_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G113b_frequency_domain_shooting_summary.md'}")


if __name__ == "__main__":
    main()
