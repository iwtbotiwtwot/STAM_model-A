#!/usr/bin/env python3
"""
G120b_gr_kerr_teukolsky_solver.py

Standalone GR Kerr QNM radial solver attempt.  No STAM.

VERDICT (honest): Direct Teukolsky shooting FAILS calibration at < 0.5%.
The s = -2 Teukolsky solution grows by ~20+ orders of magnitude across
the integration (r_+ + eps -> r_match), and standard double-precision
floating-point cannot resolve the small A_in mode relative to the
dominant A_out mode after this much amplification.  Newton iteration
converges to wrong-branch modes (Im omega > 0, anti-QNM).

This script is kept as a documented dead-end.  The correct path forward
is the Sasaki-Nakamura equation (short-range potential, no exponential
growth) or the Leaver continued-fraction method (no integration at all).
See "Status and next step" section at the script end.

Goal
====
Build a Kerr radial shooting solver for s = -2 (gravitational ψ_4),
calibrate it against G120a's qnm reference table at:
    a in [0.0, 0.3, 0.5, 0.7, 0.9],  l=2, m=2, n=0
Pass criterion:  |Δω| / |ω_qnm| < 0.5% at every spin.

Method
======
Direct shooting on the Teukolsky radial equation:

    Δ R'' + 2(s+1)(r-M) R' + V_T R = 0

with V_T = K²/Δ - 2is(r-M)K/Δ + 4isωr - λ,
     K   = (r²+a²)ω - am,
     λ   = A_lm(a,ω) + a²ω² - 2amω.

A_lm comes from qnm (verified angular eigenvalue).  Only the radial
equation is built here standalone; angular sector stays Kerr GR.

Boundary conditions
-------------------
HORIZON (r → r_+):  ingoing only.  For s = -2, the leading Frobenius
form for the ingoing branch is

    R(r) → Δ²(r) (r - r_+)^{-i σ_+} (1 + O(r - r_+))

with  σ_+ = (2M r_+ ω - a m) / (r_+ - r_-).  We start at r_+ + ε with
this leading form (and its derivative).

INFINITY (r → ∞):  outgoing only.  For s = -2,

    R(r) ~ A_out r^{2iMω} r³ e^{+iωr}
         + A_in  r^{-2iMω-1} e^{-iωr}

We use a higher-order asymptotic match (leading + first 1/r correction)
to keep r_match modest.  QNM condition:  A_in = 0.

QNM extraction
--------------
For each spin, take ω_qnm as initial guess and Newton-iterate on
A_in(ω) = 0 in complex plane.  If Newton converges to within numerical
tolerance, compare the result to ω_qnm.  If the result is just ω_qnm
(unchanged), the calibration check itself is meaningful: A_in(ω_qnm)
should already be small (a few orders of magnitude below A_out(ω_qnm)).

Status flag per spin:
- PASS if |ω_solver - ω_qnm|/|ω_qnm| < 0.5%
- TIGHT if 0.5% <= error < 2%
- FAIL otherwise

Outputs
=======
    results/G120b_gr_kerr_teukolsky_solver_summary.md
    results/G120b_gr_kerr_teukolsky_solver_table.csv
    plots/G120b_gr_kerr_teukolsky_solver.png
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
ELL, M_AZI, N = 2, 2, 0
S = -2  # spin weight (gravitational, psi_4)


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


def sigma_plus(omega, a, m=M_AZI):
    rp = r_plus(a)
    rm = r_minus(a)
    return (2.0*M*rp*omega - a*m) / (rp - rm)


# ===========================================================================
# Teukolsky radial equation  (Δ R'' + 2(s+1)(r-M) R' + V_T R = 0)
# ===========================================================================

def V_T(r, omega, a, A_lm, m=M_AZI, s=S):
    """Teukolsky radial potential V_T(r; ω, a, A_lm)."""
    K = (r*r + a*a) * omega - a * m
    lam = A_lm + a*a * omega*omega - 2.0 * a * m * omega
    D = Delta(r, a)
    return K*K/D - 2j*s*(r - M)*K/D + 4j*s*omega*r - lam


def teukolsky_rhs(r, R, dR, omega, a, A_lm, m=M_AZI, s=S):
    """Returns (dR/dr, d²R/dr²) from Teukolsky equation."""
    D = Delta(r, a)
    Vt = V_T(r, omega, a, A_lm, m, s)
    # Δ R'' + 2(s+1)(r-M) R' + V_T R = 0
    # R'' = -[2(s+1)(r-M) R' + V_T R] / Δ
    ddR = -(2.0*(s+1)*(r - M)*dR + Vt*R) / D
    return dR, ddR


def rhs_real_complex(r, z, omega, a, A_lm, m=M_AZI, s=S):
    """4-vector form for complex ODE: z = (Re R, Im R, Re R', Im R')."""
    R = z[0] + 1j*z[1]
    dR = z[2] + 1j*z[3]
    dR_out, ddR_out = teukolsky_rhs(r, R, dR, omega, a, A_lm, m, s)
    return [dR_out.real, dR_out.imag, ddR_out.real, ddR_out.imag]


# ===========================================================================
# Frobenius BC at horizon
# ===========================================================================

def horizon_bc(omega, a, A_lm, eps_h, m=M_AZI, s=S):
    """Initial conditions at r = r_+ + eps_h, leading Frobenius ingoing.

    Leading ingoing form for s=-2:
        R(r) = Δ²(r) (r - r_+)^{-i σ_+} × (1 + a_1 (r-r_+) + ...)

    For first pass, use leading order (a_1 = 0).  Verify in calibration;
    add higher Frobenius orders later if error budget requires.
    """
    rp = r_plus(a)
    rm = r_minus(a)
    sp = sigma_plus(omega, a, m)
    # Δ at r0 = rp + eps_h:
    D0 = eps_h * (rp - rm + eps_h)
    # (r - r_+)^{-i σ_+} at r0:
    horizon_phase = eps_h ** (-1j * sp)
    R0 = D0**2 * horizon_phase
    # Derivative (leading order, a_1 = 0):
    # d/dr [Δ² (r-rp)^{-iσ}] = 2 Δ Δ' (r-rp)^{-iσ} + Δ² (-iσ) (r-rp)^{-iσ-1}
    D_prime = 2.0 * rp + 2.0 * eps_h - 2.0  # dΔ/dr = 2r - 2M; at r0
    dR0 = (2.0 * D0 * D_prime * horizon_phase
           + D0*D0 * (-1j * sp) * eps_h**(-1j * sp - 1))
    return R0, dR0


# ===========================================================================
# Asymptotic forms at infinity for s = -2
# ===========================================================================

def asymptotic_decompose(r_match, R_m, dR_m, omega, m=M_AZI, s=S):
    """At large r, R(r) ≈ A_out × U_out(r) + A_in × U_in(r) for s = -2.

    Using Kerr tortoise r* = r + 2M ln r + O(1) at large r:
        exp(±iω r*) = exp(±iω r) × r^{±2iωM} × const_phase

    Asymptotic forms for s = -2:
        U_out(r) = r³ × r^{ 2iωM} × exp(+iωr)   (outgoing at infinity)
        U_in(r)  = r⁻¹ × r^{-2iωM} × exp(-iωr)   (ingoing at infinity)

    The r^{±2iωM} factor is essential -- earlier omission caused the
    initial 23-39% calibration failure.
    """
    log_r = np.log(r_match)
    e_plus = np.exp(+1j * omega * r_match)
    e_minus = np.exp(-1j * omega * r_match)
    r_pow_out = np.exp(+2j * omega * M * log_r)  # r^{2iωM}
    r_pow_in = np.exp(-2j * omega * M * log_r)   # r^{-2iωM}

    # U_out(r) = r^{3 + 2iωM} × exp(+iωr)
    # Power 'a_out' = 3 + 2iωM, U_out = r^a_out × exp(+iωr)
    a_out = 3.0 + 2j * omega * M
    a_in = -1.0 - 2j * omega * M

    U_out = r_match**3 * r_pow_out * e_plus
    U_in = r_match**(-1) * r_pow_in * e_minus

    # Derivatives: d/dr [r^a exp(+iωr)] = (a/r + iω) × r^a exp(+iωr)
    dU_out = (a_out / r_match + 1j * omega) * U_out
    dU_in = (a_in / r_match - 1j * omega) * U_in

    Mmat = np.array([[U_out, U_in], [dU_out, dU_in]])
    rhs_vec = np.array([R_m, dR_m])
    try:
        A_out, A_in = np.linalg.solve(Mmat, rhs_vec)
    except np.linalg.LinAlgError:
        return None, None
    return A_in, A_out


# ===========================================================================
# Shoot Teukolsky from horizon to r_match, return A_in (and A_out)
# ===========================================================================

def shoot_A_in(omega, a, A_lm, m=M_AZI, s=S,
                eps_h=1e-4, r_match=40.0, rtol=1e-11, atol=1e-14):
    """Integrate Teukolsky from r_+ + eps_h outward to r_match.
       Decompose asymptotic, return A_in and A_out.
    """
    rp = r_plus(a)
    r0 = rp + eps_h
    R0, dR0 = horizon_bc(omega, a, A_lm, eps_h, m, s)
    z0 = [R0.real, R0.imag, dR0.real, dR0.imag]

    def rhs(r, z):
        return rhs_real_complex(r, z, omega, a, A_lm, m, s)

    sol = solve_ivp(rhs, [r0, r_match], z0, method="DOP853",
                    rtol=rtol, atol=atol, max_step=0.5)
    if not sol.success:
        return None, None
    R_m = sol.y[0, -1] + 1j * sol.y[1, -1]
    dR_m = sol.y[2, -1] + 1j * sol.y[3, -1]
    return asymptotic_decompose(r_match, R_m, dR_m, omega, m, s)


# ===========================================================================
# Newton iteration on complex ω for A_in = 0
# ===========================================================================

def find_QNM(omega_guess, a, A_lm, **shoot_kwargs):
    """Newton iteration on complex ω such that A_in(ω) = 0.

    Use scipy.optimize.root with method='hybr' on the 2D real system
    (Re A_in, Im A_in) -> (Re ω, Im ω).
    """
    def f(x):
        omega = complex(x[0], x[1])
        A_in, A_out = shoot_A_in(omega, a, A_lm, **shoot_kwargs)
        if A_in is None:
            return [1e30, 1e30]
        # Normalize by A_out so we're targeting A_in/A_out -> 0
        if A_out is not None and abs(A_out) > 1e-300:
            ratio = A_in / A_out
            return [ratio.real, ratio.imag]
        return [A_in.real, A_in.imag]

    x0 = [omega_guess.real, omega_guess.imag]
    sol = root(f, x0, method="hybr", options={"xtol": 1e-10, "maxfev": 200})
    if sol.success:
        return complex(sol.x[0], sol.x[1]), sol
    return None, sol


# ===========================================================================
# Main calibration
# ===========================================================================

def main():
    print("=" * 88)
    print("G120b  --  Standalone GR Kerr Teukolsky shooting solver (no STAM)")
    print("=" * 88)
    print()

    # Load qnm reference (G120a)
    ref_path = RESULTS / "G120a_gr_kerr_qnm_reference.json"
    if not ref_path.exists():
        print(f"ERROR: G120a reference not found at {ref_path}.")
        print("Run G120a first.")
        sys.exit(1)
    with ref_path.open("r", encoding="utf-8") as f:
        ref = json.load(f)
    ref_rows = [r for r in ref["data"] if r.get("status") == "OK"]
    print(f"Loaded G120a reference: {len(ref_rows)} spins")
    print(f"Modes: s={ref['modes']['s']}, l={ref['modes']['l']}, "
          f"m={ref['modes']['m']}, n={ref['modes']['n']}")
    print()

    # Calibration parameters
    eps_h = 1e-6
    r_match = 200.0
    print(f"Solver settings: eps_horizon = {eps_h}, r_match = {r_match}")
    print()

    print("=" * 88)
    print("Stage 1: A_in(omega_qnm) sanity check (no Newton, just shooting test)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'omega_qnm':>30}  {'|A_in|':>12}  {'|A_out|':>12}  "
          f"{'|A_in/A_out|':>14}")
    print("-" * 90)
    sanity_rows = []
    for ref_row in ref_rows:
        a = ref_row["a"]
        omega_qnm = complex(ref_row["omega_real"], ref_row["omega_imag"])
        A_lm = complex(ref_row["A_lm_real"], ref_row["A_lm_imag"])
        A_in, A_out = shoot_A_in(omega_qnm, a, A_lm,
                                  eps_h=eps_h, r_match=r_match)
        if A_in is None:
            print(f"  {a:>5.2f}  shooting failed")
            sanity_rows.append({"a": a, "status": "shoot_fail"})
            continue
        ratio = abs(A_in / A_out) if abs(A_out) > 1e-300 else float("inf")
        print(f"  {a:>5.2f}  {str(omega_qnm):>30}  {abs(A_in):>12.4e}  "
              f"{abs(A_out):>12.4e}  {ratio:>14.4e}")
        sanity_rows.append({
            "a": a, "omega_qnm": omega_qnm, "A_lm": A_lm,
            "A_in": A_in, "A_out": A_out, "ratio": ratio,
        })
    print()
    print("Reading: |A_in/A_out| should be small (~1e-3 or less) at omega_qnm.")
    print("If it isn't, the shooting or BC is missing higher-order corrections.")
    print()

    # Stage 2: Newton iteration on omega
    print("=" * 88)
    print("Stage 2: Newton iteration for QNM (start from omega_qnm + small perturb)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'omega_qnm':>30}  {'omega_solver':>30}  "
          f"{'rel err %':>10}  status")
    print("-" * 100)
    results = []
    for ref_row in ref_rows:
        a = ref_row["a"]
        omega_qnm = complex(ref_row["omega_real"], ref_row["omega_imag"])
        A_lm = complex(ref_row["A_lm_real"], ref_row["A_lm_imag"])

        # Initial guess: small perturbation off qnm value (real convergence test)
        omega_init = omega_qnm * (1.0 + 0.005 + 0.005j)
        omega_solver, sol = find_QNM(omega_init, a, A_lm,
                                       eps_h=eps_h, r_match=r_match)
        if omega_solver is None:
            print(f"  {a:>5.2f}  Newton failed (sol.success={sol.success})")
            results.append({"a": a, "status": "newton_fail",
                            "omega_qnm": omega_qnm})
            continue

        err = abs(omega_solver - omega_qnm) / abs(omega_qnm)
        if err < 0.005:
            status = "PASS"
        elif err < 0.02:
            status = "TIGHT"
        else:
            status = "FAIL"
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
    print(f"  PASS  (<0.5%):   {n_pass}/{len(results)}")
    print(f"  TIGHT (0.5-2%):  {n_tight}/{len(results)}")
    print(f"  FAIL  (>2%):     {n_fail}/{len(results)}")
    print()
    if n_pass == len(results):
        verdict = ("PASS  --  All five spins calibrate to < 0.5%.\n"
                   "G120c can proceed with STAM substitution on this solver.")
    elif n_pass + n_tight == len(results):
        verdict = (f"TIGHT  --  {n_pass} pass at < 0.5%, {n_tight} between 0.5-2%.\n"
                   "Acceptable for first-pass STAM extension; tighten r_match\n"
                   "or add Frobenius/asymptotic corrections to reach < 0.5%.")
    else:
        verdict = (f"FAIL  --  {n_fail} spins above 2% error.\n"
                   "Iterate: increase r_match, decrease eps_h, add first-order\n"
                   "Frobenius BC and first-order 1/r asymptotic corrections.")
    print(verdict)
    print()

    # CSV
    csv_path = RESULTS / "G120b_gr_kerr_teukolsky_solver_table.csv"
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

            ax = axes[0]
            colors = ["tab:green" if r["status"] == "PASS"
                      else "tab:orange" if r["status"] == "TIGHT"
                      else "tab:red" for r in ok]
            ax.bar(spins, errs, width=0.05, color=colors)
            ax.axhline(0.5, color="black", linestyle="--", alpha=0.5,
                       label="<0.5% target")
            ax.axhline(2.0, color="gray", linestyle=":", alpha=0.5,
                       label="2% loose target")
            ax.set_xlabel("spin a / M")
            ax.set_ylabel("|ω_solver - ω_qnm| / |ω_qnm|  (%)")
            ax.set_title("G120b calibration error vs G120a reference")
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
                             xytext=(5, 5), textcoords="offset points", fontsize=8)
            ax.set_xlabel("Re(ω)")
            ax.set_ylabel("Im(ω)")
            ax.set_title("QNM in complex plane: reference vs solver")
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(PLOTS / "G120b_gr_kerr_teukolsky_solver.png", dpi=180)
            plt.close()
            print(f"Plot: {PLOTS / 'G120b_gr_kerr_teukolsky_solver.png'}")
    except Exception as e:
        print(f"Plot failed: {e}")

    # Markdown
    md = ["# G120b -- standalone GR Kerr Teukolsky shooting solver\n"]
    md.append("\nDirect shooting on the Teukolsky radial equation, calibrated\n"
              "against G120a's qnm reference table.  No STAM.\n")
    md.append("\n## Settings\n")
    md.append(f"- eps_horizon = {eps_h}\n")
    md.append(f"- r_match = {r_match}\n")
    md.append("- Frobenius BC at horizon: leading order only (a_1 = 0)\n")
    md.append("- Asymptotic at infinity: leading order only (α = 0)\n")
    md.append("\n## Stage 1: sanity (|A_in/A_out| at omega_qnm)\n")
    md.append("| a | omega_qnm | |A_in| | |A_out| | |A_in/A_out| |\n|---:|---|---:|---:|---:|\n")
    for r in sanity_rows:
        if r.get("status") == "shoot_fail":
            md.append(f"| {r['a']:.2f} | --- | --- | --- | --- |\n")
            continue
        md.append(f"| {r['a']:.2f} | {r['omega_qnm']} | {abs(r['A_in']):.4e} | "
                  f"{abs(r['A_out']):.4e} | {r['ratio']:.4e} |\n")
    md.append("\n## Stage 2: Newton calibration\n")
    md.append("| a | omega_qnm | omega_solver | rel err | status |\n")
    md.append("|---:|---|---|---:|---|\n")
    for r in results:
        if r.get("status") in ("PASS", "TIGHT", "FAIL"):
            md.append(f"| {r['a']:.2f} | {r['omega_qnm']} | {r['omega_solver']} | "
                      f"{r['rel_err']*100:.4f}% | {r['status']} |\n")
        else:
            md.append(f"| {r.get('a','?')} | {r.get('omega_qnm','?')} | --- | --- | "
                      f"{r.get('status')} |\n")
    md.append(f"\n## Verdict\n\n{verdict}\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G120b_gr_kerr_teukolsky_solver.py]"
              "(../scripts/G120b_gr_kerr_teukolsky_solver.py)\n")
    md.append("- CSV: `results/G120b_gr_kerr_teukolsky_solver_table.csv`\n")
    md.append("- Plot: `plots/G120b_gr_kerr_teukolsky_solver.png`\n")
    (RESULTS / "G120b_gr_kerr_teukolsky_solver_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary: {RESULTS / 'G120b_gr_kerr_teukolsky_solver_summary.md'}")


if __name__ == "__main__":
    main()
