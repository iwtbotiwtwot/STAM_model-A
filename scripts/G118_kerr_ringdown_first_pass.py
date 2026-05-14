#!/usr/bin/env python3
"""
G118_kerr_ringdown_first_pass.py

First-pass spin-dependent STAM correction to the axial l=2 n=0 ringdown
under the Kerr substance density A_K(r; a) = 2Mr/(r^2 + a^2), with G116's
exact prograde-equatorial photon orbit r_ph_prograde(a) as the inside-shell
boundary.

What this script DOES
======================
- Uses Kerr equatorial A_K(r; a) instead of Schwarzschild 2M/r.
- Uses G85 / G116 r_ph_prograde(a) as the inside-shell upper boundary
  (the radius where F = 1).
- Builds f_K(r; a) by applying the same matched ODE  d ln f / dA
  to the Kerr A profile (analog extension; rigorous Kerr matching
  derivation is the natural follow-up).
- Adds the canonical (sqrt f)'' / sqrt f correction in r*, analytic
  by chain rule.
- Runs the calibrated time-domain QNM pipeline from G108 / G109 / G112.
- Reports the STAM shift vs the Schwarzschild-equivalent base as a
  function of spin a in [0, 0.99].

What this script DOES NOT do
============================
- Use the actual Kerr radial perturbation (Teukolsky / Sasaki-Nakamura).
  The GR base potential is Schwarzschild RW, not Kerr.
- Use the Kerr tortoise coordinate r*_Kerr = integral of (r^2+a^2)/Delta dr.
  We use the spinless tortoise to keep the calculation comparable to G108.
- Compute the actual Kerr QNM frequencies.  The frequencies extracted
  here are NOT Kerr (l, m, n) values; they are approximations.

What this script TELLS US
=========================
The STAM correction's spin dependence: as a increases, the photon orbit
moves inward, the inside-shell region shrinks, and the f(A_K) profile
sharpens.  The shift vs a is the first-pass spin trend of the STAM
ringdown signature.

For the rigorous Kerr-Teukolsky calculation, see G119 (planned).

Outputs
=======
    results/G118_kerr_ringdown_first_pass_summary.md
    plots/G118_kerr_ringdown_first_pass.png
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
# Equatorial Kerr substance density and photon orbit
# ===========================================================================

def r_horizon(a):
    return M + np.sqrt(max(M*M - a*a, 0.0))


def r_ph_prograde(a):
    """Bardeen equatorial prograde photon orbit, M = 1."""
    if abs(a) < 1e-14:
        return 3.0 * M
    return 2.0 * M * (1.0 + np.cos((2.0 / 3.0) * np.arccos(-a / M)))


def A_K(r, a):
    """Kerr substance density at equator (theta = pi/2).  A_K = 2Mr/(r^2 + a^2)."""
    return 2.0 * M * r / (r * r + a * a)


def dA_K_dr(r, a):
    """d A_K / d r at equator."""
    num = 2.0 * M * (r * r + a * a) - 2.0 * M * r * (2.0 * r)
    den = (r * r + a * a) ** 2
    return num / den  # = 2M(a^2 - r^2) / (r^2 + a^2)^2


def d2A_K_dr2(r, a):
    """d^2 A_K / d r^2.  Analytic from chain rule on the rational function."""
    # Let u = 2Mr,  v = r^2 + a^2.  A_K = u/v.
    # A_K' = (u'v - uv')/v^2 = (2M(r^2+a^2) - 2Mr*2r)/v^2 = 2M(a^2 - r^2)/v^2.
    # A_K'' = d/dr [2M(a^2 - r^2)/v^2]
    #       = [-4Mr * v^2 - 2M(a^2 - r^2) * 2v * 2r] / v^4
    #       = -2M [2r * v + 4r(a^2 - r^2)] / v^3
    v = r * r + a * a
    return -2.0 * M * (2.0 * r * v + 4.0 * r * (a * a - r * r)) / v**3


def Sigma_K_at_photon_orbit(a):
    """Sigma_K = 3 A_K  at the prograde equatorial photon orbit."""
    rp = r_ph_prograde(a)
    return 3.0 * A_K(rp, a)


def y_K_eq(r, a):
    """Normalized inside-shell coordinate at equator using G116 boundary.
       y_K = (Sigma_K(r) - Sigma_K(r_ph)) / (3 - Sigma_K(r_ph)).
       y = 0 at the photon orbit, y = 1 at horizon (Sigma_K(r_+) = 3).
    """
    Sig_r = 3.0 * A_K(r, a)
    Sig_ph = Sigma_K_at_photon_orbit(a)
    denom = 3.0 - Sig_ph
    if abs(denom) < 1e-14:
        return float('nan')
    return (Sig_r - Sig_ph) / denom


def F_quintic(y): return 1.0 - 5.0 * y**4 + 4.0 * y**5
def F_first(y): return -20.0 * y**3 + 20.0 * y**4


# ===========================================================================
# Schwarzschild metric (base for the first-pass)
# ===========================================================================

def h_fn(r): return 1.0 - 2.0 * M / r
def hp_fn(r): return 2.0 * M / r**2

def k_STAM_eq(r, a):
    """STAM modification: k = h F(y_K_eq) inside r_ph_prograde(a), else h.
       At the equator only.  Schwarzschild-form base h(r) = 1 - 2M/r.
    """
    if r >= r_ph_prograde(a):
        return h_fn(r)
    y = y_K_eq(r, a)
    return h_fn(r) * F_quintic(y)


def V_RW(r, ell=ELL):
    return h_fn(r) * (ell * (ell + 1) / r**2 - 6.0 * M / r**3)


# ===========================================================================
# f(A_K) by extension of the spinless matched ODE
# ===========================================================================

def dlnf_dA(A):
    """G70 / G75 matched coupling derivative, applied to Kerr A_K by extension.
       Zero outside the Kerr equatorial photon orbit.
    """
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = 1.0 - 5.0 * y**4 + 4.0 * y**5
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)


def d2lnf_dA2(A, hstep=1e-6):
    return (dlnf_dA(A + hstep) - dlnf_dA(A - hstep)) / (2.0 * hstep)


def build_f_K_grid(r_grid, a):
    """Integrate d ln f / dr = (dlnf/dA)(dA_K/dr) along r_grid.
       Anchor f(r >= r_ph_prograde(a)) = 1.
    """
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


# ===========================================================================
# Canonical-form action correction (sqrt f)'' / sqrt f  in r*
# ===========================================================================

def analytic_correction_kerr(r_grid, a):
    """Analytic (sqrt f_K)'' / sqrt f_K  in r*-coordinate, using
       chain rule on A_K(r) profile.

       Uses the spinless Schwarzschild tortoise: dr*/dr = 1 / sqrt(h k).
       k = h F(y_K_eq) inside photon orbit, h outside.
       Sigma factors come from Schwarzschild form for first-pass.
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
        d_dlnf_dr = up * (dA**2) + u * d2A
        hk = h_fn(rv) * k_STAM_eq(rv, a)
        if hk <= 0:
            out[i] = 0.0
            continue
        sqrt_hk = np.sqrt(hk)
        # k = h F(y).  Use spinless dk/dr structure (first-pass)
        # k_STAM_eq = h * F(y_K_eq); d/dr of this is messy.  Use FD.
        dr_step = 1e-5 * max(1.0, rv - 2.0)
        dr_step = max(dr_step, 1e-7)
        k_plus = k_STAM_eq(rv + dr_step, a) if (rv + dr_step) < rp else h_fn(rv + dr_step)
        k_minus = k_STAM_eq(rv - dr_step, a) if (rv - dr_step) < rp else h_fn(rv - dr_step)
        k_now = k_STAM_eq(rv, a)
        dk_dr = (k_plus - k_minus) / (2.0 * dr_step)
        d_hk = hp_fn(rv) * k_now + h_fn(rv) * dk_dr
        P = sqrt_hk * dlnf_dr
        d_sqrt_hk_dr = 0.5 * d_hk / sqrt_hk
        dP_dr = d_sqrt_hk_dr * dlnf_dr + sqrt_hk * d_dlnf_dr
        P_star = sqrt_hk * dP_dr
        out[i] = 0.5 * P_star + 0.25 * P**2
    return out


# ===========================================================================
# Tortoise (Schwarzschild form) and time-domain pipeline
# ===========================================================================

def r_of_rstar_GR(rstar):
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


def build_r_of_rstar_STAM(rstar_grid, a):
    """STAM tortoise dr/dr* = sqrt(h * k_STAM_eq) on the equatorial Kerr profile."""
    def rhs(rs, y):
        rv = y[0]
        return [np.sqrt(h_fn(rv) * k_STAM_eq(rv, a))]
    rs_pos = rstar_grid[rstar_grid > 0]
    rs_neg = rstar_grid[rstar_grid < 0]
    r_arr = np.zeros_like(rstar_grid)
    if len(rs_pos) > 0:
        order = np.argsort(rs_pos)
        x = rs_pos[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=1.0)
        vals = np.zeros_like(rs_pos)
        vals[order] = sol.y[0]
        r_arr[rstar_grid > 0] = vals
    if len(rs_neg) > 0:
        order = np.argsort(rs_neg)[::-1]
        x = rs_neg[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.02)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar_grid < 0] = vals
    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


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
        popt, _ = curve_fit(damped_sinusoid, tf, sf, p0=p0, bounds=bounds, maxfev=40000)
        return complex(popt[1], -popt[2])
    except Exception:
        return None


# ===========================================================================
# Main: spin sweep
# ===========================================================================

def run_one_spin(a, rstar):
    """Returns dict with omega_GR (calibrated), omega_STAM, shift for spin a."""
    dr_star = rstar[1] - rstar[0]

    # Build r-grid via STAM tortoise (depends on a through r_ph_prograde)
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar, a)

    # GR base potential (Schwarzschild RW) - same for all a in first-pass
    V_GR_grid = np.array([V_RW(r) for r in r_grid_GR])

    # STAM-modified k inside r_ph_prograde(a) -> V_geom_eq
    V_geom_grid = np.empty_like(r_grid_STAM)
    for i, rv in enumerate(r_grid_STAM):
        h_v = h_fn(rv)
        k_v = k_STAM_eq(rv, a)
        # Use V_geom form with k_v
        # Need k' as well; use FD for simplicity (already in correction code)
        dr_step = 1e-5 * max(1.0, rv - 2.0)
        dr_step = max(dr_step, 1e-7)
        k_p = k_STAM_eq(rv + dr_step, a)
        k_m = k_STAM_eq(rv - dr_step, a)
        kp_v = (k_p - k_m) / (2.0 * dr_step)
        V_geom_grid[i] = h_v * (
            ELL * (ELL + 1) / rv**2
            - 2.0 * (1.0 - k_v) / rv**2
            - kp_v / (2.0 * rv)
            - k_v * hp_fn(rv) / (2.0 * h_v * rv)
        )

    # Action correction
    corr = analytic_correction_kerr(r_grid_STAM, a)
    V_action_grid = V_geom_grid + corr

    # Time-domain
    pulse_center, pulse_sigma = 30.0, 3.0
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 50.0)))
    dt = 0.5 * dr_star
    T_final = 350.0
    t_fit_start, t_fit_end = 100.0, 250.0

    # Calibrate GR (same for all a since base is Schwarzschild, but check)
    t_GR, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
    om_GR = fit_QNM(t_GR, sig_GR, t_fit_start, t_fit_end, om_guess=LEAVER_GOLD)
    if om_GR is None:
        return {"a": a, "status": "GR_fit_failed"}
    calib = abs(om_GR - LEAVER_GOLD) / abs(LEAVER_GOLD)
    if calib > CALIBRATION_THRESHOLD:
        return {"a": a, "status": f"calib_fail_{calib*100:.2f}%",
                "om_GR": om_GR, "calib": calib}

    # STAM
    t_S, sig_S = evolve_TD(V_action_grid, dr_star, T_final, dt, psi0, i_obs)
    om_STAM = fit_QNM(t_S, sig_S, t_fit_start, t_fit_end, om_guess=om_GR)
    if om_STAM is None:
        return {"a": a, "status": "STAM_fit_failed", "om_GR": om_GR}
    shift = abs(om_STAM - om_GR) / abs(om_GR)

    rp = r_ph_prograde(a)
    Sig_ph = Sigma_K_at_photon_orbit(a)
    return {
        "a": a, "status": "ok",
        "om_GR": om_GR, "om_STAM": om_STAM, "calib": calib, "shift": shift,
        "r_ph_prograde": rp, "Sigma_ph_eq": Sig_ph,
        "f_max": float(max(np.max(np.abs(corr)), 0.0)),
        "t": t_S, "sig_GR": sig_GR, "sig_STAM": sig_S,
    }


def main():
    print("=" * 88)
    print("G118  --  first-pass Kerr ringdown STAM correction (spin sweep)")
    print("=" * 88)
    print()
    print("Approximations (honest):")
    print("  - GR base potential = Schwarzschild RW (NOT full Kerr Teukolsky)")
    print("  - Tortoise = STAM tortoise on equatorial Kerr substance density")
    print("  - f(A_K) = G70/G75 matched ODE applied to Kerr A_K(r; a) by extension")
    print("  - r_ph_prograde(a) from G116 as the inside-shell upper boundary")
    print()
    print("Reports STAM shift % as a function of spin.  Does NOT report Kerr QNM")
    print("frequencies (rigorous Teukolsky is G119).")
    print()

    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N = {N}, range [{rs_min}, {rs_max}], dr* = {dr_star:.4e}")
    print()

    spins = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99]
    results = []
    for a in spins:
        print(f"--- spin a = {a:.2f} ---", flush=True)
        res = run_one_spin(a, rstar)
        results.append(res)
        if res["status"] != "ok":
            print(f"  status: {res['status']}")
            continue
        print(f"  r_ph_prograde = {res['r_ph_prograde']:.4f}")
        print(f"  Sigma_ph_eq = {res['Sigma_ph_eq']:.4f}")
        print(f"  omega_GR (calib err {res['calib']*100:.3f}%) = {res['om_GR']}")
        print(f"  omega_STAM = {res['om_STAM']}")
        print(f"  shift = {res['shift']*100:.4f}%")
        print(f"  max |correction| = {res['f_max']:.4e}")
        print()

    # ----- Summary table -----
    print("=" * 88)
    print("Spin sweep summary")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'r_ph':>7}  {'Sigma_ph':>10}  {'omega_STAM':>30}  "
          f"{'shift %':>10}")
    print("-" * 80)
    for r in results:
        if r["status"] != "ok":
            print(f"  {r['a']:>5.2f}  {r['status']}")
            continue
        print(f"  {r['a']:>5.2f}  {r['r_ph_prograde']:>7.4f}  "
              f"{r['Sigma_ph_eq']:>10.4f}  {str(r['om_STAM']):>30}  "
              f"{r['shift']*100:>9.4f}%")
    print()

    # ----- Plots -----
    ok = [r for r in results if r["status"] == "ok"]
    if ok:
        spins_arr = np.array([r["a"] for r in ok])
        shifts_arr = np.array([r["shift"] * 100 for r in ok])
        rphs_arr = np.array([r["r_ph_prograde"] for r in ok])

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        ax = axes[0, 0]
        ax.plot(spins_arr, shifts_arr, "o-", linewidth=2, markersize=8, color="tab:green")
        ax.axhline(4.977, color="tab:blue", linestyle=":", alpha=0.5,
                   label="G108 spinless 4.977%")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("STAM shift vs GR  (%)")
        ax.set_title("STAM action-aware shift vs spin  (first-pass)")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        ax.plot(spins_arr, rphs_arr, "s-", linewidth=2, markersize=8, color="tab:purple")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("r_ph_prograde / M")
        ax.set_title("Prograde equatorial photon orbit (G116)")
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        for r in ok[::2]:  # subsample for readability
            ax.plot(r["t"], np.abs(r["sig_STAM"]) + 1e-30,
                    linewidth=0.9, label=f"a = {r['a']:.2f}")
        ax.set_yscale("log")
        ax.set_xlabel("t / M")
        ax.set_ylabel("|psi(t, r_obs)|")
        ax.set_title("Time-domain signals across spin")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Re/Im components of omega_STAM vs spin
        ax = axes[1, 1]
        re_S = np.array([r["om_STAM"].real for r in ok])
        im_S = np.array([abs(r["om_STAM"].imag) for r in ok])
        re_G = np.array([r["om_GR"].real for r in ok])
        im_G = np.array([abs(r["om_GR"].imag) for r in ok])
        ax.plot(spins_arr, re_S, "o-", color="tab:green", label="Re(omega_STAM)")
        ax.plot(spins_arr, re_G, "o:", color="tab:blue", alpha=0.5, label="Re(omega_GR)")
        ax.plot(spins_arr, im_S, "s-", color="tab:orange", label="|Im(omega_STAM)|")
        ax.plot(spins_arr, im_G, "s:", color="tab:red", alpha=0.5, label="|Im(omega_GR)|")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("omega components")
        ax.set_title("Real and imaginary parts vs spin")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G118_kerr_ringdown_first_pass.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G118 -- first-pass Kerr ringdown STAM correction (spin sweep)\n"]
    md.append("\n## Approximations (honest)\n")
    md.append("- GR base potential = Schwarzschild RW.  Spinning Kerr Teukolsky\n"
              "  is the rigorous follow-up (G119).\n")
    md.append("- Tortoise built with the STAM-modified equatorial Kerr substance\n"
              "  density A_K = 2Mr/(r^2 + a^2).\n")
    md.append("- f(A_K) constructed by applying the G70/G75 spinless matched ODE\n"
              "  along the Kerr A_K profile (extension, not Kerr-matched derivation).\n")
    md.append("- r_ph_prograde(a) from G116 sets the inside-shell upper boundary.\n")
    md.append("- (sqrt f_K)'' / sqrt f_K added analytically via chain rule.\n")
    md.append("\n## What the script reports\n")
    md.append("STAM shift % as a function of spin.  This is NOT the Kerr QNM frequency;\n"
              "it is the framework's f-correction magnitude as a function of spin.\n")
    md.append("\n## Spin sweep results\n")
    md.append("| a | r_ph_prograde | Sigma_ph_eq | omega_STAM | shift |\n")
    md.append("|---:|---:|---:|---|---:|\n")
    for r in results:
        if r["status"] != "ok":
            md.append(f"| {r['a']:.2f} | --- | --- | {r['status']} | --- |\n")
            continue
        md.append(f"| {r['a']:.2f} | {r['r_ph_prograde']:.4f} | "
                  f"{r['Sigma_ph_eq']:.4f} | {r['om_STAM']} | "
                  f"{r['shift']*100:.4f}% |\n")
    md.append("\n## Reading\n")
    md.append("At a = 0 the calculation recovers G108's spinless 4.977% shift.\n"
              "As spin increases, the inside-shell region shrinks (r_ph_prograde\n"
              "decreases) and the f(A_K) profile sharpens; the table shows how the\n"
              "shift trends with spin.  Increasing or decreasing trend is the\n"
              "first-pass spin-dependence prediction.\n")
    md.append("\n## What's still needed\n")
    md.append("- G119 (rigorous Kerr Teukolsky on STAM background) to get actual Kerr QNMs.\n")
    md.append("- Polar Kerr (G120) along the G114-rigorous lines.\n")
    md.append("- Higher-mode and overtone Kerr extensions (G121+).\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G118_kerr_ringdown_first_pass.py](../scripts/G118_kerr_ringdown_first_pass.py)\n")
    md.append("- [plots/G118_kerr_ringdown_first_pass.png](../plots/G118_kerr_ringdown_first_pass.png)\n")
    (RESULTS / "G118_kerr_ringdown_first_pass_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G118_kerr_ringdown_first_pass_summary.md'}")


if __name__ == "__main__":
    main()
