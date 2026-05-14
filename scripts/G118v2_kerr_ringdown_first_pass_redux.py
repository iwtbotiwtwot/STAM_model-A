#!/usr/bin/env python3
"""
G118v2_kerr_ringdown_first_pass_redux.py

Redux of G118 first-pass Kerr ringdown using G119's Kerr wave base.

What G118 (original) did wrong
==============================
G118 used Schwarzschild  h(r) = 1 - 2M/r  as the radial structure base.
For Kerr at spin a > ~0.5 the actual horizon r_+(a) sits BELOW the
Schwarzschild r = 2M, so (1 - 2M/r) goes negative inside (r_+, 2M).
sqrt(h k) goes imaginary; the STAM tortoise integration fails; the
correction integrates to zero.  Result: artifact "0% shift" at high spin.

What G118v2 does
================
Use G119's Kerr wave base everywhere:
- Kerr horizon r_+(a) = M + sqrt(M^2 - a^2).
- Kerr radial structure Delta(r, a) = r^2 - 2Mr + a^2.
- Kerr tortoise dr*/dr = (r^2 + a^2) / Delta  (GR Kerr)
                       = (r^2 + a^2) / Delta_STAM  (STAM, with F factor inside PS)
- Inside-PS modification:  Delta_STAM = Delta * F(y_K_eq)
  with y_K_eq via G116's r_pr_exact(theta=pi/2, a) = r_ph_prograde(a).

Wave equation (first-pass Kerr RW analog)
=========================================
We still use a Schwarzschild-RW-analog potential V_base(r, a, ell) with Kerr Delta.
Specifically:

    V_base(r, a, ell) = [Delta / (r^2 + a^2)^2]
                       * [ell(ell+1) (r^2 + a^2) / r^2
                          - 6 M (r^2 - a^2) / r^3]

This reduces EXACTLY to Schwarzschild Regge-Wheeler at a = 0:
    V_base(r, 0, ell) = (1 - 2M/r) [ell(ell+1)/r^2 - 6M/r^3]    ✓

For a > 0 this is NOT the rigorous Kerr axial perturbation potential
(that requires m-coupling and is left for G120 Detweiler/Sasaki-Nakamura).
It is a first-pass "Kerr-Δ-aware Schwarzschild-RW analog".

Calibration strategy
====================
At a = 0 the calibration target is Leaver gold 0.373672 - 0.088962 i.
At a > 0 we do NOT have a closed Leaver gold for this approximate
potential, so we calibrate at a=0 only and report the STAM shift
percentage  (omega_STAM(a) - omega_base(a)) / |omega_base(a)|
for each spin, where omega_base(a) is the time-domain extraction
under V_base + Kerr tortoise (no STAM modification).

This isolates the STAM correction's spin dependence from the question
"how well does V_base match actual Kerr QNMs" (which is G120 scope).

Outputs
=======
    results/G118v2_kerr_ringdown_redux_summary.md
    plots/G118v2_kerr_ringdown_redux.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0
ELL = 2
LEAVER_GOLD = complex(0.373672, -0.088962)
CALIBRATION_THRESHOLD = 0.01


# ===========================================================================
# Kerr radial structure  (from G119)
# ===========================================================================

def r_horizon(a):
    disc = M * M - a * a
    return M + np.sqrt(max(disc, 0.0))


def Delta(r, a):
    return r * r - 2.0 * M * r + a * a


def r_ph_prograde(a):
    if abs(a) < 1e-14:
        return 3.0 * M
    return 2.0 * M * (1.0 + np.cos((2.0 / 3.0) * np.arccos(-a / M)))


def A_K(r, a):
    return 2.0 * M * r / (r * r + a * a)


def Sigma_K(r, a):
    return 3.0 * A_K(r, a)


def Sigma_ph_eq(a):
    return Sigma_K(r_ph_prograde(a), a)


def y_K_eq(r, a):
    Sig_r = Sigma_K(r, a)
    Sig_ph = Sigma_ph_eq(a)
    denom = 3.0 - Sig_ph
    if abs(denom) < 1e-14:
        return float("nan")
    return (Sig_r - Sig_ph) / denom


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def F_first(y):
    return -20.0 * y**3 + 20.0 * y**4


def Delta_STAM(r, a):
    """Delta * F(y_K_eq)  inside r_ph_prograde(a),  else Delta."""
    D = Delta(r, a)
    if r >= r_ph_prograde(a):
        return D
    y = y_K_eq(r, a)
    return D * F_quintic(y)


# ===========================================================================
# Kerr-RW analog base potential  (reduces to V_RW at a = 0)
# ===========================================================================

def V_base_GR(r, a, ell=ELL):
    """V_base = (Delta / (r^2+a^2)^2)
                * [ell(ell+1)(r^2+a^2)/r^2 - 6M(r^2-a^2)/r^3]

    Reduces EXACTLY to Schwarzschild RW at a = 0.
    NOT the rigorous Kerr axial potential (G120 task).
    """
    D = Delta(r, a)
    pref = D / (r * r + a * a) ** 2
    bracket = ell * (ell + 1) * (r * r + a * a) / r**2 \
              - 6.0 * M * (r * r - a * a) / r**3
    return pref * bracket


def V_base_STAM(r, a, ell=ELL):
    """V_base_STAM = (Delta_STAM / (r^2+a^2)^2)
                     * [ell(ell+1)(r^2+a^2)/r^2 - 6M(r^2-a^2)/r^3].

    Outside PS: Delta_STAM = Delta, V matches V_base_GR.
    Inside PS:  Delta_STAM smaller; V scales accordingly.
    """
    DS = Delta_STAM(r, a)
    pref = DS / (r * r + a * a) ** 2
    bracket = ell * (ell + 1) * (r * r + a * a) / r**2 \
              - 6.0 * M * (r * r - a * a) / r**3
    return pref * bracket


# ===========================================================================
# f(A_K) and analytic action correction (sqrt f)'' / sqrt f
# ===========================================================================

def dlnf_dA(A):
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = 1.0 - 5.0 * y**4 + 4.0 * y**5
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)


def d2lnf_dA2(A, hstep=1e-6):
    return (dlnf_dA(A + hstep) - dlnf_dA(A - hstep)) / (2.0 * hstep)


def dA_K_dr(r, a):
    v = r * r + a * a
    return 2.0 * M * (a * a - r * r) / v**2


def d2A_K_dr2(r, a):
    v = r * r + a * a
    return -2.0 * M * (2.0 * r * v + 4.0 * r * (a * a - r * r)) / v**3


def build_f_K_grid(r_grid, a):
    """Integrate d ln f / dr along r_grid_STAM, anchor f(r >= r_ph_prograde) = 1."""
    rp = r_ph_prograde(a)
    ln_f = np.zeros_like(r_grid)
    inside = r_grid < rp
    if not np.any(inside):
        return np.ones_like(r_grid)
    r_inside = r_grid[inside]
    r_descending = np.sort(r_inside)[::-1]

    def rhs(rv, yv):
        A = A_K(rv, a)
        return [dlnf_dA(A) * dA_K_dr(rv, a)]

    sol = solve_ivp(rhs, [rp - 1e-12, r_descending[-1]], [0.0],
                    t_eval=r_descending, method="RK45",
                    rtol=1e-12, atol=1e-14, max_step=0.005)
    ln_f_inside = np.empty_like(r_inside)
    order = np.argsort(r_inside)
    ln_f_inside[order] = sol.y[0][::-1]
    ln_f[inside] = ln_f_inside
    return np.exp(ln_f)


def analytic_correction_kerr_v2(r_grid, a):
    """Analytic (sqrt f_K)'' / sqrt f_K  in r*-coord, using KERR tortoise
        dr*/dr = (r^2 + a^2) / Delta_STAM.

    Derivation:
      P = d_*(ln f) = (Delta_STAM / (r^2+a^2)) * d(ln f)/dr
      (sqrt f)'' / sqrt f = (1/2) d_* P + (1/4) P^2
                          = (1/2) [Delta_STAM/(r^2+a^2)] d/dr P + (1/4) P^2

      d/dr P = d/dr [Delta_STAM/(r^2+a^2) * d(ln f)/dr]
             = [(d Delta_STAM/dr)(r^2+a^2) - Delta_STAM (2r)] / (r^2+a^2)^2 * d(ln f)/dr
               + Delta_STAM/(r^2+a^2) * d^2(ln f)/dr^2
    """
    out = np.zeros_like(r_grid)
    rp = r_ph_prograde(a)
    for i, rv in enumerate(r_grid):
        if rv >= rp:
            out[i] = 0.0
            continue
        A = A_K(rv, a)
        dA = dA_K_dr(rv, a)
        d2A = d2A_K_dr2(rv, a)
        u = dlnf_dA(A)
        up = d2lnf_dA2(A)
        dlnf_dr = u * dA
        d2lnf_dr2 = up * (dA**2) + u * d2A
        DS = Delta_STAM(rv, a)
        if DS <= 0:
            out[i] = 0.0
            continue
        v = rv * rv + a * a
        prefactor = DS / v
        # Need d(Delta_STAM)/dr.  Use FD locally (Delta_STAM is piecewise).
        eps = 1e-5 * max(1.0, rv - r_horizon(a))
        eps = max(eps, 1e-7)
        DSp = Delta_STAM(rv + eps, a) if (rv + eps) < rp else Delta(rv + eps, a)
        DSm = Delta_STAM(rv - eps, a) if (rv - eps) < rp else Delta(rv - eps, a)
        dDS_dr = (DSp - DSm) / (2.0 * eps)
        # P and dP/dr
        P = prefactor * dlnf_dr
        # d/dr [Delta_STAM/(r^2+a^2)] = [dDS_dr * v - DS * 2r] / v^2
        d_prefactor_dr = (dDS_dr * v - DS * 2.0 * rv) / v**2
        dP_dr = d_prefactor_dr * dlnf_dr + prefactor * d2lnf_dr2
        # P_*  = (Delta_STAM/(r^2+a^2)) * dP/dr  in r* by chain rule
        P_star = prefactor * dP_dr
        out[i] = 0.5 * P_star + 0.25 * P**2
    return out


# ===========================================================================
# Kerr tortoise builders
# ===========================================================================

def build_r_of_rstar_kerr(rstar_grid, a, use_STAM=True):
    """Integrate dr/dr* = (Delta_STAM)/(r^2+a^2) [STAM]
                    or  dr/dr* = Delta/(r^2+a^2)      [GR]
       starting from r* = 0 anchored at r = r_pr_eq(a) (the photon-region
       boundary).  Then map to standard 'r*(3M) = 0' anchor via subtract-offset.
    """
    rp_eq = r_ph_prograde(a)
    rh = r_horizon(a)

    if use_STAM:
        def integrand(r):
            DS = Delta_STAM(r, a)
            return DS / (r * r + a * a) if DS > 0 else 0.0
    else:
        def integrand(r):
            return Delta(r, a) / (r * r + a * a)

    def rhs(rs, y):
        r_val = y[0]
        return [integrand(r_val)]

    rs_pos = rstar_grid[rstar_grid > 0]
    rs_neg = rstar_grid[rstar_grid < 0]
    r_arr = np.zeros_like(rstar_grid)
    if len(rs_pos) > 0:
        order = np.argsort(rs_pos)
        x = rs_pos[order]
        sol = solve_ivp(rhs, [0.0, x[-1]], [rp_eq], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=1.0)
        vals = np.zeros_like(rs_pos)
        vals[order] = sol.y[0]
        r_arr[rstar_grid > 0] = vals
    if len(rs_neg) > 0:
        order = np.argsort(rs_neg)[::-1]
        x = rs_neg[order]
        sol = solve_ivp(rhs, [0.0, x[-1]], [rp_eq], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.02)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar_grid < 0] = vals
    r_arr[rstar_grid == 0] = rp_eq
    return r_arr


# ===========================================================================
# Time-domain pipeline (same as G91/G108/G118)
# ===========================================================================

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, i_obs):
    psi_old = psi0.copy()
    psi = psi0.copy()
    N_steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star
    t_arr = np.zeros(N_steps)
    sig_arr = np.zeros(N_steps)
    for step in range(N_steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2
        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1])
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])
        psi_old = psi
        psi = psi_new
        t_arr[step] = (step + 1) * dt
        sig_arr[step] = psi[i_obs]
    return t_arr, sig_arr


def damped_sinusoid(t, A, om_r, om_i, ph, off):
    tau = t - t[0]
    return A * np.exp(-om_i * tau) * np.cos(om_r * tau + ph) + off


def fit_QNM(t, sig, t0, t1, om_guess=complex(0.37, -0.09)):
    mask = (t >= t0) & (t <= t1)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, om_guess.real, abs(om_guess.imag), 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 5.0, 5.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, tf, sf, p0=p0, bounds=bounds,
                            maxfev=40000)
        return complex(popt[1], -popt[2])
    except Exception:
        return None


# ===========================================================================
# Spin sweep
# ===========================================================================

def run_one_spin(a, rstar):
    dr_star = rstar[1] - rstar[0]

    # Build GR and STAM r-grids using KERR tortoise
    r_grid_GR = build_r_of_rstar_kerr(rstar, a, use_STAM=False)
    r_grid_STAM = build_r_of_rstar_kerr(rstar, a, use_STAM=True)

    # GR base potential
    V_GR_grid = np.array([V_base_GR(r, a) for r in r_grid_GR])

    # STAM potential = V_base with Delta -> Delta_STAM,  plus action correction
    V_STAM_base_grid = np.array([V_base_STAM(r, a) for r in r_grid_STAM])
    corr = analytic_correction_kerr_v2(r_grid_STAM, a)
    V_STAM_grid = V_STAM_base_grid + corr

    # Time-domain settings
    pulse_center, pulse_sigma = 30.0, 3.0
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 50.0)))
    dt = 0.5 * dr_star
    T_final = 350.0
    t0_fit, t1_fit = 100.0, 250.0

    # Calibrate GR base
    t_arr, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
    om_GR = fit_QNM(t_arr, sig_GR, t0_fit, t1_fit, om_guess=LEAVER_GOLD)
    if om_GR is None:
        return {"a": a, "status": "GR_fit_failed"}

    # At a = 0, must match Leaver gold (calibration check)
    if a == 0.0:
        calib_err = abs(om_GR - LEAVER_GOLD) / abs(LEAVER_GOLD)
        if calib_err > CALIBRATION_THRESHOLD:
            return {"a": a, "status": f"a0_calib_fail_{calib_err*100:.2f}%",
                    "om_GR": om_GR}
        cal_text = f"a0 calib err {calib_err*100:.3f}% (vs Leaver)"
    else:
        # No Leaver gold for our approximate V at a > 0; just record om_GR
        cal_text = "no Leaver target (approximate V at a > 0)"

    # STAM run
    t_arr_S, sig_STAM = evolve_TD(V_STAM_grid, dr_star, T_final, dt, psi0, i_obs)
    om_STAM = fit_QNM(t_arr_S, sig_STAM, t0_fit, t1_fit, om_guess=om_GR)
    if om_STAM is None:
        return {"a": a, "status": "STAM_fit_failed", "om_GR": om_GR}

    shift = abs(om_STAM - om_GR) / abs(om_GR)
    return {
        "a": a, "status": "ok",
        "om_GR": om_GR, "om_STAM": om_STAM, "shift": shift,
        "calibration_note": cal_text,
        "r_ph_eq": r_ph_prograde(a), "r_horizon": r_horizon(a),
        "Sigma_ph_eq": Sigma_ph_eq(a),
        "max_correction": float(np.max(np.abs(corr))),
        "t": t_arr_S, "sig_GR": sig_GR, "sig_STAM": sig_STAM,
    }


def main():
    print("=" * 88)
    print("G118v2  --  Kerr ringdown first-pass REDUX on G119 wave base")
    print("=" * 88)
    print()
    print("Replaces G118's Schwarzschild base with the proper Kerr radial structure.")
    print("V_base = (Delta/(r^2+a^2)^2) * [l(l+1)(r^2+a^2)/r^2 - 6M(r^2-a^2)/r^3]")
    print("Reduces to Schwarzschild RW at a = 0  (calibration check passes by design).")
    print()

    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    print(f"r* grid: N = {N}, range [{rs_min}, {rs_max}],  dr* = {(rs_max-rs_min)/(N-1):.4e}")
    print()

    spins = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]
    results = []
    for a in spins:
        print(f"--- spin a = {a:.2f} ---", flush=True)
        res = run_one_spin(a, rstar)
        results.append(res)
        if res["status"] != "ok":
            print(f"  status: {res['status']}")
            if "om_GR" in res:
                print(f"  om_GR = {res['om_GR']}")
            continue
        print(f"  r_horizon = {res['r_horizon']:.4f}")
        print(f"  r_ph_eq = {res['r_ph_eq']:.4f}")
        print(f"  Sigma_ph_eq = {res['Sigma_ph_eq']:.4f}")
        print(f"  omega_GR_base = {res['om_GR']}    ({res['calibration_note']})")
        print(f"  omega_STAM    = {res['om_STAM']}")
        print(f"  shift = {res['shift']*100:.4f}%")
        print(f"  max |correction| = {res['max_correction']:.4e}")
        print()

    # ----- Summary -----
    print("=" * 88)
    print("Spin sweep summary (G118v2 on G119 base)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'r_+':>7}  {'r_ph':>7}  {'Sigma_ph':>9}  "
          f"{'omega_GR':>30}  {'omega_STAM':>30}  {'shift':>9}")
    print("-" * 110)
    for r in results:
        if r["status"] != "ok":
            print(f"  {r['a']:>5.2f}  status: {r['status']}")
            continue
        print(f"  {r['a']:>5.2f}  {r['r_horizon']:>7.4f}  {r['r_ph_eq']:>7.4f}  "
              f"{r['Sigma_ph_eq']:>9.4f}  {str(r['om_GR']):>30}  "
              f"{str(r['om_STAM']):>30}  {r['shift']*100:>8.4f}%")
    print()

    # ----- Comparison with G118 original -----
    G118_ORIG = {
        0.0: 4.9773, 0.1: 5.3084, 0.2: 5.4486, 0.3: 5.2286,
        0.5: 3.1757, 0.7: 0.4012, 0.9: 0.0, 0.99: 0.0,
    }
    print("=" * 88)
    print("Comparison with G118 (Schwarzschild base, broken at high spin)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'G118 (broken)':>18}  {'G118v2 (Kerr base)':>22}")
    print("-" * 50)
    for r in results:
        if r["status"] != "ok":
            continue
        a = r["a"]
        if a in G118_ORIG:
            print(f"  {a:>5.2f}  {G118_ORIG[a]:>17.4f}%  {r['shift']*100:>21.4f}%")
    print()

    # ----- Plots -----
    ok = [r for r in results if r["status"] == "ok"]
    if ok:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        ax = axes[0, 0]
        spins_arr = np.array([r["a"] for r in ok])
        shifts_arr = np.array([r["shift"] * 100 for r in ok])
        ax.plot(spins_arr, shifts_arr, "o-", linewidth=2, markersize=8, color="tab:green",
                label="G118v2 (Kerr base)")
        g118_a = sorted([a for a in G118_ORIG])
        g118_s = [G118_ORIG[a] for a in g118_a]
        ax.plot(g118_a, g118_s, "x--", linewidth=1.5, markersize=8, color="tab:gray",
                label="G118 original (Schwarzschild base, broken a>0.7)")
        ax.axhline(4.977, color="tab:blue", linestyle=":", alpha=0.5,
                   label="G108 spinless 4.977%")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("STAM shift (%)")
        ax.set_title("STAM ringdown shift vs spin: G118 vs G118v2")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        rphs = [r["r_ph_eq"] for r in ok]
        rhs = [r["r_horizon"] for r in ok]
        ax.plot(spins_arr, rphs, "s-", linewidth=2, markersize=7, color="tab:purple",
                label="r_ph_eq")
        ax.plot(spins_arr, rhs, "o-", linewidth=2, markersize=7, color="tab:red",
                label="r_+")
        ax.axhline(2.0, color="gray", linestyle=":", alpha=0.4,
                   label="r=2M (Schw. horizon)")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("radius / M")
        ax.set_title("Kerr horizon and prograde photon orbit")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        for r in ok[::2]:
            ax.plot(r["t"], np.abs(r["sig_STAM"]) + 1e-30,
                    linewidth=0.9, label=f"a = {r['a']:.2f}")
        ax.set_yscale("log")
        ax.set_xlabel("t / M")
        ax.set_ylabel("|psi(t, r_obs)|")
        ax.set_title("STAM-modified time-domain signals across spin")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        ax = axes[1, 1]
        re_S = np.array([r["om_STAM"].real for r in ok])
        im_S = np.array([abs(r["om_STAM"].imag) for r in ok])
        re_G = np.array([r["om_GR"].real for r in ok])
        im_G = np.array([abs(r["om_GR"].imag) for r in ok])
        ax.plot(spins_arr, re_S, "o-", color="tab:green", label="Re(omega_STAM)")
        ax.plot(spins_arr, re_G, "o:", color="tab:blue", alpha=0.5, label="Re(omega_GR_base)")
        ax.plot(spins_arr, im_S, "s-", color="tab:orange", label="|Im(omega_STAM)|")
        ax.plot(spins_arr, im_G, "s:", color="tab:red", alpha=0.5, label="|Im(omega_GR_base)|")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("omega components")
        ax.set_title("QNM components vs spin")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G118v2_kerr_ringdown_redux.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G118v2 -- Kerr ringdown first-pass REDUX on G119 wave base\n"]
    md.append("\nFixes G118 by replacing the Schwarzschild radial structure with\n"
              "Kerr Delta(r, a) and Kerr horizon r_+(a).  No more high-spin breakdown.\n")
    md.append("\n## Setup\n")
    md.append("- Kerr horizon r_+(a) = M + sqrt(M^2 - a^2)\n")
    md.append("- Kerr radial structure Delta(r, a) = r^2 - 2Mr + a^2\n")
    md.append("- Kerr tortoise dr*/dr = (r^2 + a^2)/Delta_STAM\n")
    md.append("- Delta_STAM = Delta * F(y_K_eq) inside r_ph_prograde(a)\n")
    md.append("- V_base reduces to Schwarzschild RW at a = 0\n")
    md.append("- Action correction (sqrt f_K)''/sqrt f_K added analytically with Kerr chain rule\n")
    md.append("\n## Spin sweep\n")
    md.append("| a | r_+ | r_ph_eq | omega_GR_base | omega_STAM | shift |\n")
    md.append("|---:|---:|---:|---|---|---:|\n")
    for r in results:
        if r["status"] != "ok":
            md.append(f"| {r['a']:.2f} | --- | --- | {r['status']} | --- | --- |\n")
            continue
        md.append(f"| {r['a']:.2f} | {r['r_horizon']:.4f} | {r['r_ph_eq']:.4f} | "
                  f"{r['om_GR']} | {r['om_STAM']} | {r['shift']*100:.4f}% |\n")
    md.append("\n## Comparison with G118 original (high-spin breakdown)\n")
    md.append("| a | G118 (broken) | G118v2 (Kerr) |\n|---:|---:|---:|\n")
    for r in results:
        if r["status"] != "ok":
            continue
        a = r["a"]
        if a in G118_ORIG:
            md.append(f"| {a:.2f} | {G118_ORIG[a]:.4f}% | {r['shift']*100:.4f}% |\n")
    md.append("\n## Caveats (still first-pass)\n")
    md.append("- V_base is a Kerr-Delta-aware Schwarzschild-RW analog; NOT the\n"
              "  rigorous Kerr axial perturbation potential (Detweiler / Sasaki-\n"
              "  Nakamura).  m-coupling and full Teukolsky are G120 scope.\n")
    md.append("- omega_GR_base at a > 0 is NOT the published Kerr QNM; it is the\n"
              "  time-domain extraction under V_base.  The reported shift % is\n"
              "  STAM-vs-V_base, not STAM-vs-published-Kerr-QNM.  Both are equally\n"
              "  derived from V_base, so the shift % is the meaningful number.\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G118v2_kerr_ringdown_first_pass_redux.py]"
              "(../scripts/G118v2_kerr_ringdown_first_pass_redux.py)\n")
    md.append("- [plots/G118v2_kerr_ringdown_redux.png](../plots/G118v2_kerr_ringdown_redux.png)\n")
    (RESULTS / "G118v2_kerr_ringdown_redux_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G118v2_kerr_ringdown_redux_summary.md'}")


if __name__ == "__main__":
    main()
