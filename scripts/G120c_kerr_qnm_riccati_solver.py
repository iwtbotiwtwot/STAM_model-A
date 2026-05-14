#!/usr/bin/env python3
"""
G120c_kerr_qnm_riccati_solver.py

Standalone GR Kerr QNM solver via log-derivative (Riccati) matching on
the Teukolsky radial equation.  No STAM.

Goal
====
Reproduce the G120a qnm reference table at:
    a in [0.0, 0.3, 0.5, 0.7, 0.9],  l = m = 2,  n = 0
to better than 0.5%.

Why log-derivative matching
===========================
G120b's direct amplitude shooting failed because Teukolsky R(r) for
s = -2 grows by ~10^22 across the integration; the small ingoing-mode
coefficient drowned in floating-point noise.

Log-derivative  Y = R'/R  is bounded throughout the integration domain
(no exponential amplitude growth).  The QNM condition becomes a
mismatch in Y at a matching point r_match:

    Y_left(r_match)  -  Y_right(r_match)  =  0   at the QNM omega.

This is mathematically equivalent to a Sasaki-Nakamura log-derivative
match (the numerical-robustness gain comes from log-derivative, not from
SN specifically), but uses the verified Teukolsky V_T from G120b's helper
code without re-deriving SN coefficients.

Riccati equation
================
From  d Δ R'' + 2(s+1)(r-M) R' + V_T R = 0,  with Y = R'/R:

    R'' / R  =  (R'/R)' + (R'/R)²  =  Y' + Y²

Substituting:
    Δ (Y' + Y²) + 2(s+1)(r-M) Y + V_T  =  0
=>  Y'  =  -[2(s+1)(r-M) Y + V_T] / Δ  -  Y²

For s = -2:
    Y'  =  [2(r-M) Y - V_T] / Δ  -  Y²

Boundary conditions
===================
Near horizon (r = r_+ + eps):
    R(r) ~ Delta² (r - r_+)^{-i sigma_+} (1 + a_1 (r-r_+) + ...)
    Y(r_+ + eps) ~ (2 - i sigma_+) / eps   (leading, a_1 = 0)

Near infinity (r large):
    R(r) ~ r^{3 + 2i omega M} exp(+i omega r) (1 + alpha/r + ...)
    Y(r) ~ i omega + (3 + 2 i omega M) / r   (leading)

Settings
========
eps_h   = 1e-5   (start  r_+ + eps_h  for left shoot)
r_max   = 400.0  (start  r_max for right shoot, going inward)
r_match = 25.0   (matching point; well between r_+ and r_max)

If the leading BCs hit < 0.5% calibration: done.
If close but slightly off: add first-order correction (a_1 or alpha/r²).

Outputs
=======
    results/G120c_kerr_qnm_riccati_solver_summary.md
    results/G120c_kerr_qnm_riccati_solver_table.csv
    plots/G120c_kerr_qnm_riccati_solver.png
"""

from __future__ import annotations

import sys
import json
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0
ELL, M_AZI, N_OVERTONE = 2, 2, 0
S = -2  # spin weight


# ===========================================================================
# Kerr geometry
# ===========================================================================

def r_plus(a):
    return M + np.sqrt(max(M*M - a*a, 0.0))


def r_minus(a):
    return M - np.sqrt(max(M*M - a*a, 0.0))


def Delta(r, a):
    return r*r - 2.0*M*r + a*a


def Omega_H(a):
    rp = r_plus(a)
    return a / (rp*rp + a*a)


def K_func(r, omega, a, m=M_AZI):
    return (r*r + a*a) * omega - a * m


def sigma_plus(omega, a, m=M_AZI):
    rp = r_plus(a)
    rm = r_minus(a)
    return (2.0*M*rp*omega - a*m) / (rp - rm)


# ===========================================================================
# Teukolsky V_T(r; omega, a, A_lm)  (verified in G120b)
# ===========================================================================

def V_T(r, omega, a, A_lm, m=M_AZI, s=S):
    K = K_func(r, omega, a, m)
    lam = A_lm + a*a * omega*omega - 2.0 * a * m * omega
    D = Delta(r, a)
    return K*K/D - 2j*s*(r - M)*K/D + 4j*s*omega*r - lam


# ===========================================================================
# Riccati equation for Y = R' / R
# ===========================================================================

def riccati_rhs(r, Y, omega, a, A_lm, m=M_AZI, s=S):
    """Y' = [2(r-M) Y - V_T] / Delta - Y^2     (for s = -2)

    General s:  Y' = -[2(s+1)(r-M) Y + V_T] / Delta - Y^2.
    """
    Vt = V_T(r, omega, a, A_lm, m, s)
    D = Delta(r, a)
    return -(2.0*(s+1)*(r - M)*Y + Vt) / D - Y*Y


def rhs_real(r, y_pair, omega, a, A_lm, m=M_AZI, s=S):
    """Real 2-vector form: y_pair = (Re Y, Im Y)."""
    Y = y_pair[0] + 1j * y_pair[1]
    Yp = riccati_rhs(r, Y, omega, a, A_lm, m, s)
    return [Yp.real, Yp.imag]


# ===========================================================================
# Boundary conditions
# ===========================================================================

def Y_horizon_leading(omega, a, eps_h, m=M_AZI):
    """Leading horizon log-derivative:
       Y ≈ (2 - i σ_+) / eps_h    (for s = -2)
       comes from  R ~ Delta² (r - r_+)^{-i σ_+}
       =>  d ln R / dr = 2 Delta'/Delta + (-i σ_+)/(r - r_+)
       At r = r_+ + eps:  Delta' = 2 r_+ - 2M = r_+ - r_-,  Delta ≈ eps (r_+ - r_-)
       =>  Y ≈ 2/eps - i σ_+/eps  =  (2 - i σ_+)/eps.
    """
    sp = sigma_plus(omega, a, m)
    return (2.0 - 1j * sp) / eps_h


def Y_infinity_leading(omega, r_far):
    """Leading infinity log-derivative for s = -2:
       Y(r) ~ i omega + (3 + 2 i omega M) / r.
    """
    return 1j * omega + (3.0 + 2.0j * omega * M) / r_far


# ===========================================================================
# Shoot Y from horizon outward to r_match
# ===========================================================================

def shoot_Y_left(omega, a, A_lm, eps_h, r_match, m=M_AZI, s=S):
    """Integrate Riccati from r_+ + eps_h outward to r_match.  Return Y(r_match)."""
    rp = r_plus(a)
    r0 = rp + eps_h
    Y0 = Y_horizon_leading(omega, a, eps_h, m)
    y0 = [Y0.real, Y0.imag]

    def rhs(r, y):
        return rhs_real(r, y, omega, a, A_lm, m, s)

    sol = solve_ivp(rhs, [r0, r_match], y0, method="DOP853",
                    rtol=1e-10, atol=1e-13, max_step=0.5)
    if not sol.success:
        return None
    return sol.y[0, -1] + 1j * sol.y[1, -1]


def shoot_Y_right(omega, a, A_lm, r_max, r_match, m=M_AZI, s=S):
    """Integrate Riccati from r_max inward to r_match.  Return Y(r_match)."""
    Y_max = Y_infinity_leading(omega, r_max)
    y0 = [Y_max.real, Y_max.imag]

    def rhs(r, y):
        return rhs_real(r, y, omega, a, A_lm, m, s)

    sol = solve_ivp(rhs, [r_max, r_match], y0, method="DOP853",
                    rtol=1e-10, atol=1e-13, max_step=0.5)
    if not sol.success:
        return None
    return sol.y[0, -1] + 1j * sol.y[1, -1]


def mismatch(omega, a, A_lm, eps_h, r_max, r_match, m=M_AZI, s=S):
    """Y_left(r_match) - Y_right(r_match).  QNM => mismatch = 0."""
    Y_L = shoot_Y_left(omega, a, A_lm, eps_h, r_match, m, s)
    Y_R = shoot_Y_right(omega, a, A_lm, r_max, r_match, m, s)
    if Y_L is None or Y_R is None:
        return None
    return Y_L - Y_R


# ===========================================================================
# Newton solver on complex omega
# ===========================================================================

def find_QNM(omega_guess, a, A_lm, eps_h, r_max, r_match,
              m=M_AZI, s=S, xtol=1e-9, maxfev=200):
    def f(x):
        omega = complex(x[0], x[1])
        mis = mismatch(omega, a, A_lm, eps_h, r_max, r_match, m, s)
        if mis is None or not np.isfinite(mis):
            return [1e30, 1e30]
        return [mis.real, mis.imag]

    x0 = [omega_guess.real, omega_guess.imag]
    sol = root(f, x0, method="hybr",
               options={"xtol": xtol, "maxfev": maxfev})
    if sol.success:
        return complex(sol.x[0], sol.x[1]), sol
    return None, sol


# ===========================================================================
# Main calibration
# ===========================================================================

def main():
    print("=" * 88)
    print("G120c -- Riccati / log-derivative GR Kerr QNM solver (no STAM)")
    print("=" * 88)
    print()

    ref_path = RESULTS / "G120a_gr_kerr_qnm_reference.json"
    if not ref_path.exists():
        print(f"ERROR: G120a reference not found at {ref_path}.  Run G120a first.")
        sys.exit(1)
    with ref_path.open("r", encoding="utf-8") as f:
        ref = json.load(f)
    ref_rows = [r for r in ref["data"] if r.get("status") == "OK"]
    print(f"Loaded G120a reference: {len(ref_rows)} spins")
    print(f"Modes: s={ref['modes']['s']}, l={ref['modes']['l']}, "
          f"m={ref['modes']['m']}, n={ref['modes']['n']}")
    print()

    # Calibration settings
    eps_h = 1e-5
    r_match = 25.0
    r_max = 400.0
    print(f"Settings:  eps_h = {eps_h},  r_match = {r_match},  r_max = {r_max}")
    print()

    # Stage 1: sanity at qnm omega
    print("=" * 88)
    print("Stage 1: mismatch(omega_qnm) sanity check")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'omega_qnm':>30}  {'|Y_L - Y_R|':>14}  "
          f"{'|Y_L|':>10}  {'|Y_R|':>10}")
    print("-" * 80)
    sanity_rows = []
    for ref_row in ref_rows:
        a = ref_row["a"]
        omega = complex(ref_row["omega_real"], ref_row["omega_imag"])
        A_lm = complex(ref_row["A_lm_real"], ref_row["A_lm_imag"])
        Y_L = shoot_Y_left(omega, a, A_lm, eps_h, r_match)
        Y_R = shoot_Y_right(omega, a, A_lm, r_max, r_match)
        if Y_L is None or Y_R is None:
            print(f"  {a:>5.2f}  shoot failed")
            sanity_rows.append({"a": a, "status": "fail"})
            continue
        mis = Y_L - Y_R
        print(f"  {a:>5.2f}  {str(omega):>30}  {abs(mis):>14.4e}  "
              f"{abs(Y_L):>10.4e}  {abs(Y_R):>10.4e}")
        sanity_rows.append({"a": a, "Y_L": Y_L, "Y_R": Y_R, "mismatch": mis})
    print()
    print("Reading: |Y_L - Y_R| should be small at omega_qnm.")
    print("If close to |Y_L| or |Y_R|, BC corrections are needed.")
    print()

    # Stage 2: Newton iteration
    print("=" * 88)
    print("Stage 2: Newton on mismatch(omega) = 0 (initial guess: perturbed omega_qnm)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'omega_qnm':>30}  {'omega_solver':>30}  "
          f"{'rel err':>10}  status")
    print("-" * 100)
    results = []
    for ref_row in ref_rows:
        a = ref_row["a"]
        omega_qnm = complex(ref_row["omega_real"], ref_row["omega_imag"])
        A_lm = complex(ref_row["A_lm_real"], ref_row["A_lm_imag"])

        omega_init = omega_qnm * (1.0 + 0.005 + 0.005j)
        omega_solver, sol = find_QNM(omega_init, a, A_lm,
                                       eps_h, r_max, r_match)
        if omega_solver is None:
            print(f"  {a:>5.2f}  Newton failed")
            results.append({"a": a, "status": "newton_fail",
                            "omega_qnm": omega_qnm})
            continue
        err = abs(omega_solver - omega_qnm) / abs(omega_qnm)
        status = "PASS" if err < 0.005 else "TIGHT" if err < 0.02 else "FAIL"
        print(f"  {a:>5.2f}  {str(omega_qnm):>30}  {str(omega_solver):>30}  "
              f"{err*100:>9.4f}%  {status}")
        results.append({
            "a": a, "omega_qnm": omega_qnm, "omega_solver": omega_solver,
            "rel_err": err, "status": status,
        })
    print()

    # Summary
    n_pass = sum(1 for r in results if r.get("status") == "PASS")
    n_tight = sum(1 for r in results if r.get("status") == "TIGHT")
    n_fail = sum(1 for r in results if r.get("status") == "FAIL")
    print("=" * 88)
    print("CALIBRATION SUMMARY")
    print("=" * 88)
    print()
    print(f"  PASS  (< 0.5%):   {n_pass}/{len(results)}")
    print(f"  TIGHT (0.5-2%):   {n_tight}/{len(results)}")
    print(f"  FAIL  (>= 2%):    {n_fail}/{len(results)}")
    print()
    if n_pass == len(results):
        verdict = (
            "PASS  --  All five spins calibrate to < 0.5%.  Log-derivative\n"
            "matching tames the numerical-growth problem of G120b.\n"
            "G120d (STAM Delta_STAM substitution) may now proceed."
        )
    elif n_pass + n_tight == len(results):
        verdict = (
            f"TIGHT  --  {n_pass} pass at < 0.5%, {n_tight} at 0.5-2%.\n"
            "Acceptable; add first-order BC corrections (a_1 at horizon,\n"
            "alpha/r at infinity) to tighten further."
        )
    else:
        verdict = (
            f"FAIL  --  {n_fail} spins above 2%.\n"
            "Check Riccati equation, BC implementation, r_max, r_match."
        )
    print(verdict)
    print()

    # CSV
    csv_path = RESULTS / "G120c_kerr_qnm_riccati_solver_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["a", "omega_qnm_real", "omega_qnm_imag",
                    "omega_solver_real", "omega_solver_imag",
                    "rel_error", "status"])
        for r in results:
            if r.get("status") in ("PASS", "TIGHT", "FAIL"):
                w.writerow([r["a"], r["omega_qnm"].real, r["omega_qnm"].imag,
                            r["omega_solver"].real, r["omega_solver"].imag,
                            r["rel_err"], r["status"]])
            else:
                w.writerow([r.get("a"), r.get("omega_qnm", complex()).real,
                            r.get("omega_qnm", complex()).imag,
                            "", "", "", r.get("status")])
    print(f"CSV: {csv_path}")

    # Plot
    try:
        import matplotlib.pyplot as plt
        ok = [r for r in results if r.get("status") in ("PASS", "TIGHT", "FAIL")]
        if ok:
            fig, axes = plt.subplots(1, 2, figsize=(13, 5))
            spins = np.array([r["a"] for r in ok])
            errs = np.array([r["rel_err"] * 100 for r in ok])
            colors = ["tab:green" if r["status"] == "PASS"
                      else "tab:orange" if r["status"] == "TIGHT"
                      else "tab:red" for r in ok]

            ax = axes[0]
            ax.bar(spins, errs, width=0.05, color=colors)
            ax.axhline(0.5, color="black", linestyle="--", alpha=0.5,
                       label="< 0.5% target")
            ax.axhline(2.0, color="gray", linestyle=":", alpha=0.5,
                       label="2% loose target")
            ax.set_xlabel("spin a / M")
            ax.set_ylabel("|ω_solver - ω_qnm| / |ω_qnm|  (%)")
            ax.set_title("G120c calibration error vs G120a reference")
            ax.set_yscale("log")
            ax.legend()
            ax.grid(True, alpha=0.3)

            ax = axes[1]
            re_q = np.array([r["omega_qnm"].real for r in ok])
            im_q = np.array([r["omega_qnm"].imag for r in ok])
            re_s = np.array([r["omega_solver"].real for r in ok])
            im_s = np.array([r["omega_solver"].imag for r in ok])
            ax.plot(re_q, im_q, "o-", linewidth=2, markersize=10,
                    color="tab:blue", label="ω_qnm (reference)")
            ax.plot(re_s, im_s, "s", markersize=8, color="tab:green",
                    label="ω_solver", alpha=0.7)
            for r, rq, iq in zip(ok, re_q, im_q):
                ax.annotate(f"a={r['a']:.1f}", (rq, iq),
                             xytext=(5, 5), textcoords="offset points",
                             fontsize=8)
            ax.set_xlabel("Re(ω)")
            ax.set_ylabel("Im(ω)")
            ax.set_title("QNM track: reference vs solver")
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(PLOTS / "G120c_kerr_qnm_riccati_solver.png", dpi=180)
            plt.close()
            print(f"Plot: {PLOTS / 'G120c_kerr_qnm_riccati_solver.png'}")
    except Exception as e:
        print(f"Plot failed: {e}")

    # Markdown
    md = ["# G120c -- Riccati / log-derivative GR Kerr QNM solver\n"]
    md.append("\nLog-derivative matching on the Teukolsky radial equation, "
              "with the log-derivative Y = R'/R bounded throughout integration.\n")
    md.append("\n## Settings\n")
    md.append(f"- eps_h = {eps_h}\n")
    md.append(f"- r_max = {r_max}\n")
    md.append(f"- r_match = {r_match}\n")
    md.append("- horizon BC: leading order (2 - i σ_+) / eps_h\n")
    md.append("- infinity BC: leading order i ω + (3 + 2 i ω M) / r\n")
    md.append("\n## Stage 1: mismatch at omega_qnm\n")
    md.append("| a | omega_qnm | |Y_L - Y_R| | |Y_L| | |Y_R| |\n")
    md.append("|---:|---|---:|---:|---:|\n")
    for r in sanity_rows:
        if r.get("status") == "fail":
            md.append(f"| {r['a']:.2f} | --- | shoot fail | --- | --- |\n")
            continue
        md.append(f"| {r['a']:.2f} | {complex(r['Y_L'])} | "
                  f"{abs(r['mismatch']):.4e} | {abs(r['Y_L']):.4e} | "
                  f"{abs(r['Y_R']):.4e} |\n")
    md.append("\n## Stage 2: Newton calibration\n")
    md.append("| a | omega_qnm | omega_solver | rel err | status |\n")
    md.append("|---:|---|---|---:|---|\n")
    for r in results:
        if r.get("status") in ("PASS", "TIGHT", "FAIL"):
            md.append(f"| {r['a']:.2f} | {r['omega_qnm']} | {r['omega_solver']} | "
                      f"{r['rel_err']*100:.4f}% | {r['status']} |\n")
        else:
            md.append(f"| {r.get('a','?')} | --- | --- | --- | {r.get('status')} |\n")
    md.append(f"\n## Verdict\n\n{verdict}\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G120c_kerr_qnm_riccati_solver.py]"
              "(../scripts/G120c_kerr_qnm_riccati_solver.py)\n")
    md.append("- CSV: `results/G120c_kerr_qnm_riccati_solver_table.csv`\n")
    md.append("- Plot: `plots/G120c_kerr_qnm_riccati_solver.png`\n")
    (RESULTS / "G120c_kerr_qnm_riccati_solver_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary: {RESULTS / 'G120c_kerr_qnm_riccati_solver_summary.md'}")


if __name__ == "__main__":
    main()
