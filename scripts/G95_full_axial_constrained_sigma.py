#!/usr/bin/env python3
"""
G95_full_axial_constrained_sigma.py

Fuller axial perturbation attack for the constrained shell-count action.

This script is the next step after G86-redux and G91-G94.  It separates the
calculation into two layers:

  (A) Action-level axial reduction logic.
  (B) Numerical axial diagnostic using the canonical f(Σ)R correction.

Action being tested
===================

    S[Σ,g,λ1,λ2] = (1/16πG) ∫ d^4x sqrt(-g) [
          f(Σ) R
        + λ1 ( (∇Σ)^2 - W )
        + λ2 ( u^μ ∂_μΣ )
        - 2 V(Σ)
    ]

with A = Σ/3 in the weak-field continuum limit.

Axial-sector reduction claim
============================
Odd-parity / axial metric perturbations do not carry an independent shell-count
perturbation:

  * Σ is a scalar shell-count variable and is parity-even.
  * G79-G84 double-LM constraints remove independent δΣ propagation.
  * For axial odd-parity modes, δΣ can consistently be set to zero.
  * λ1 and λ2 terms therefore do not supply a propagating axial scalar mode.

The remaining axial dynamics are tensor/graviton dynamics with variable
Planck stiffness f(Σ).  For the canonical master variable ψ = sqrt(f) q,
the leading constrained-action axial potential is:

    V_action = V_proxy + (sqrt(f))'' / sqrt(f)

where primes are derivatives with respect to the STAM tortoise coordinate r*.

This is still a controlled axial reduction, not a hand-derived full Regge-Wheeler
textbook recurrence.  Its purpose is to make the axial-sector assumptions explicit,
verify the photon-sphere derivative structure, and reproduce the calibrated
G91-G94 time-domain result in one self-contained script.

Outputs
=======
  results/G95_full_axial_constrained_sigma_summary.md
  results/G95_full_axial_derivative_table.csv
  plots/G95_full_axial_constrained_sigma.png

Usage
=====
  python G95_full_axial_constrained_sigma.py              # structural + TD run
  python G95_full_axial_constrained_sigma.py --derive-only
  python G95_full_axial_constrained_sigma.py --quick       # shorter TD grid

Conventions
===========
M = 1.  Frequencies are reported as ω = Re(ω) - i|Im(ω)|.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

M = 1.0

# ---------------------------------------------------------------------------
# Metric / branch functions
# ---------------------------------------------------------------------------

def h_fn(r):
    return 1.0 - 2.0 * M / r


def A_of_r(r):
    return 2.0 * M / r


def Sigma_of_r(r):
    return 3.0 * A_of_r(r)


def y_of_r(r):
    return Sigma_of_r(r) - 2.0


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def k_GR(r):
    return h_fn(r)


def k_STAM(r):
    y = y_of_r(r)
    h_val = h_fn(r)
    return np.where(np.asarray(y) > 0, h_val * F_quintic(y), h_val)


def V_RW_proxy(r, ell, k_fn):
    """Axial Regge-Wheeler proxy potential used in G87/G91-G94.

    M_eff = r(1-k)/2.  For Schwarzschild k=h, M_eff=M.
    """
    h_val = h_fn(r)
    k_val = k_fn(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


# ---------------------------------------------------------------------------
# f(Σ) reconstruction from G70/G75 in A-variable form
# ---------------------------------------------------------------------------

def dlnf_dA(A):
    """G70/G75 ODE for f in the quintic branch.

    d ln f / dA = -2 y^3 (45 A^2 - 33 A - 13) / [A F(y)], y=3A-2.
    Outside/at PS, f=1.
    """
    A_arr = np.asarray(A, dtype=float)
    y = 3.0 * A_arr - 2.0
    F = F_quintic(y)
    out = np.zeros_like(A_arr)
    mask = A_arr > (2.0 / 3.0)
    out[mask] = -2.0 * y[mask]**3 * (45.0 * A_arr[mask]**2 - 33.0 * A_arr[mask] - 13.0) / (
        A_arr[mask] * F[mask]
    )
    return out if np.ndim(A) else float(out)


def build_f_of_A_interpolator(A_max: float, n: int = 40000, clip_lnf: float = 120.0):
    A0 = 2.0 / 3.0
    A_max = min(float(A_max), 1.0 - 1e-10)
    if A_max <= A0:
        def f_interp(A):
            return np.ones_like(np.asarray(A, dtype=float))
        return f_interp

    A_grid = np.linspace(A0, A_max, n)
    deriv = dlnf_dA(A_grid)
    lnf = cumulative_trapezoid(deriv, A_grid, initial=0.0)
    lnf = np.clip(lnf, -clip_lnf, clip_lnf)
    f_grid = np.exp(lnf)

    def f_interp(A):
        A_arr = np.asarray(A, dtype=float)
        out = np.ones_like(A_arr)
        mask = A_arr > A0
        out[mask] = np.interp(A_arr[mask], A_grid, f_grid, left=1.0, right=f_grid[-1])
        return out

    return f_interp


def f_grid_on_r(r_grid, f_norm=1.0):
    A_grid = A_of_r(r_grid)
    f_interp = build_f_of_A_interpolator(float(np.nanmax(A_grid)))
    return f_norm * f_interp(A_grid)


# ---------------------------------------------------------------------------
# Tortoise coordinate construction
# ---------------------------------------------------------------------------

def r_of_rstar_GR(rstar: float) -> float:
    """Invert Schwarzschild r*, anchored so r*(3M)=0."""
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def froot(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    if rstar > 0:
        return brentq(froot, 3.0 * M - 1e-12, 1e6)
    return brentq(froot, 2.0 * M + 1e-14, 3.0 * M + 1e-12)


def build_r_of_rstar_STAM(rstar_grid, max_step_in=0.04):
    """Build r(r*) for STAM with dr/dr* = sqrt(h k)."""
    def rhs(rs, y):
        r_val = y[0]
        val = h_fn(r_val) * k_STAM(r_val)
        return [np.sqrt(max(float(val), 0.0))]

    r_arr = np.zeros_like(rstar_grid)

    rs_pos = rstar_grid[rstar_grid > 0]
    if len(rs_pos):
        sort_pos = np.argsort(rs_pos)
        rs_sorted = rs_pos[sort_pos]
        sol = solve_ivp(rhs, [0.0, rs_sorted[-1]], [3.0 * M], t_eval=rs_sorted,
                        rtol=1e-12, atol=1e-14, max_step=1.0)
        tmp = np.zeros_like(rs_pos)
        tmp[sort_pos] = sol.y[0]
        r_arr[rstar_grid > 0] = tmp

    rs_neg = rstar_grid[rstar_grid < 0]
    if len(rs_neg):
        sort_neg = np.argsort(rs_neg)[::-1]
        rs_sorted = rs_neg[sort_neg]
        sol = solve_ivp(rhs, [0.0, rs_sorted[-1]], [3.0 * M], t_eval=rs_sorted,
                        rtol=1e-12, atol=1e-14, max_step=max_step_in)
        tmp = np.zeros_like(rs_neg)
        tmp[sort_neg] = sol.y[0]
        r_arr[rstar_grid < 0] = tmp

    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


# ---------------------------------------------------------------------------
# Axial action potential
# ---------------------------------------------------------------------------

def canonical_f_correction(f_grid, dr_star):
    """Return (sqrt(f))''/sqrt(f), derivatives in r*."""
    sqrtf = np.sqrt(np.maximum(f_grid, 1e-300))
    d1 = np.gradient(sqrtf, dr_star, edge_order=2)
    d2 = np.gradient(d1, dr_star, edge_order=2)
    return d2 / sqrtf


def build_action_aware_potential(V_proxy, f_grid, dr_star):
    corr = canonical_f_correction(f_grid, dr_star)
    return V_proxy + corr, corr


# ---------------------------------------------------------------------------
# TD solver and extraction
# ---------------------------------------------------------------------------

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, psidot0=None, i_obs=None):
    if psidot0 is None:
        psidot0 = np.zeros_like(psi0)

    psi_old = psi0 - dt * psidot0
    psi = psi0.copy()
    steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star
    t_arr = np.zeros(steps)
    sig_arr = np.zeros(steps)

    for step in range(steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2

        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1]
        )

        # Sommerfeld outflow in r*.
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])

        psi_old, psi = psi, psi_new
        t_arr[step] = (step + 1) * dt
        if i_obs is not None:
            sig_arr[step] = psi[i_obs]

    return t_arr, sig_arr


def damped_sinusoid(t, Aamp, omega_r, omega_i, phi, offset):
    return Aamp * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + offset


def extract_QNM(t, signal, t_fit_start, t_fit_end, init_omega_r=0.37, init_omega_i=0.09):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    t_fit = t[mask]
    s_fit = signal[mask]
    if len(t_fit) < 50 or np.allclose(s_fit, 0):
        return None, None

    A_guess = 0.5 * (np.nanmax(s_fit) - np.nanmin(s_fit))
    p0 = [A_guess, init_omega_r, init_omega_i, 0.0, 0.0]
    bounds = ([-np.inf, 0.0, 0.0, -np.pi, -np.inf], [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0, bounds=bounds, maxfev=50000)
        _, omega_r, omega_i, _, _ = popt
        return complex(abs(omega_r), -abs(omega_i)), popt
    except Exception as ex:
        print(f"  curve_fit failed in window [{t_fit_start}, {t_fit_end}]: {ex}")
        return None, None


# ---------------------------------------------------------------------------
# Local series derivative proof near PS
# ---------------------------------------------------------------------------

def local_derivative_table(ell_val=2, max_deriv=8):
    """Use G86-redux local series logic to find first derivative mismatches."""
    z = sp.symbols("z", real=True)
    u = sp.symbols("u", real=True)
    ell = sp.symbols("ell", positive=True)
    ORDER = 12

    def trunc(expr, order=ORDER):
        expr = sp.expand(expr)
        poly = sp.Poly(expr, z)
        out = 0
        for (deg,), coeff in poly.terms():
            if deg < order:
                out += coeff * z**deg
        return sp.expand(out)

    def series_expr(expr, order=ORDER):
        return sp.series(expr, z, 0, order).removeO()

    def coeff_derivative(expr, n):
        poly = sp.Poly(sp.expand(expr), z)
        c = poly.coeff_monomial(z**n)
        return sp.simplify(c * math.factorial(n))

    def compose_poly_in_u(poly_u, u_z, order=ORDER):
        P = sp.Poly(sp.expand(poly_u), u)
        out = 0
        pow_u = {0: sp.Integer(1)}
        max_deg = P.degree()
        for n in range(1, max_deg + 1):
            pow_u[n] = trunc(pow_u[n-1] * u_z, order)
        for (deg,), coeff in P.terms():
            out += coeff * pow_u[deg]
        return trunc(out, order)

    def exp_series(poly_z, order=ORDER):
        out = sp.Integer(1)
        term = sp.Integer(1)
        for m in range(1, order + 1):
            term = trunc(term * poly_z, order)
            if term == 0:
                break
            out = trunc(out + term / math.factorial(m), order)
        return out

    def sqrt_one_plus_series(delta_z, order=ORDER):
        out = sp.Integer(1)
        pow_delta = sp.Integer(1)
        for m in range(1, order + 1):
            pow_delta = trunc(pow_delta * delta_z, order)
            if pow_delta == 0:
                break
            out = trunc(out + sp.binomial(sp.Rational(1, 2), m) * pow_delta, order)
        return out

    def Dstar_series(expr, Dcoef, order=ORDER):
        return trunc(Dcoef * sp.diff(expr, z), order)

    r = 3 + z
    A = sp.Rational(2, 1) / r
    y = 3 * A - 2
    F = series_expr(1 - 5 * y**4 + 4 * y**5, ORDER)
    h = series_expr(1 - A, ORDER)
    k_STAM_ser = trunc(h * F, ORDER)

    M_eff_GR = sp.Integer(1)
    M_eff_proxy = trunc(r * (1 - k_STAM_ser) / 2, ORDER)
    V_GR_ser = series_expr(h * (ell * (ell + 1) / r**2 - 6 * M_eff_GR / r**3), ORDER)
    V_proxy_ser = series_expr(h * (ell * (ell + 1) / r**2 - 6 * M_eff_proxy / r**3), ORDER)

    # sqrt(f) local series.
    A_u = sp.Rational(2, 3) + u
    y_u = 3 * A_u - 2
    F_u = 1 - 5 * y_u**4 + 4 * y_u**5
    q_u = -2 * y_u**3 * (45 * A_u**2 - 33 * A_u - 13) / (A_u * F_u)
    q_u_series = sp.series(q_u, u, 0, ORDER).removeO()
    lnf_u_series = sp.integrate(q_u_series, (u, 0, u))
    u_z = series_expr(A - sp.Rational(2, 3), ORDER)
    lnf_z = compose_poly_in_u(lnf_u_series, u_z, ORDER)
    sqrtf_z = exp_series(trunc(sp.Rational(1, 2) * lnf_z, ORDER), ORDER)

    sqrtF_ser = sqrt_one_plus_series(trunc(F - 1, ORDER), ORDER)
    Dcoef_GR = h
    Dcoef_STAM = trunc(h * sqrtF_ser, ORDER)

    sqrtf_star_1 = Dstar_series(sqrtf_z, Dcoef_STAM)
    sqrtf_star_2 = Dstar_series(sqrtf_star_1, Dcoef_STAM)
    action_term_ser = trunc(sp.series(sqrtf_star_2 / sqrtf_z, z, 0, ORDER).removeO(), ORDER)
    V_action_ser = trunc(V_proxy_ser + action_term_ser, ORDER)

    rows = []
    for name, V_stam, V_base in [
        ("proxy_minus_GR", V_proxy_ser, V_GR_ser),
        ("action_minus_GR", V_action_ser, V_GR_ser),
    ]:
        first_r = None
        first_star = None
        for n in range(max_deriv + 1):
            rdiff = sp.simplify(coeff_derivative((V_stam - V_base).subs(ell, ell_val), n))
            stam_star = V_stam
            base_star = V_base
            for _ in range(n):
                stam_star = Dstar_series(stam_star, Dcoef_STAM)
                base_star = Dstar_series(base_star, Dcoef_GR)
            sdiff = sp.simplify((stam_star - base_star).subs({z: 0, ell: ell_val}))
            if first_r is None and rdiff != 0:
                first_r = (n, rdiff)
            if first_star is None and sdiff != 0:
                first_star = (n, sdiff)
            rows.append({
                "comparison": name,
                "derivative_order": n,
                "r_derivative_diff": str(rdiff),
                "rstar_derivative_diff": str(sdiff),
                "r_derivative_diff_float": float(sp.N(rdiff)),
                "rstar_derivative_diff_float": float(sp.N(sdiff)),
            })
    Vps = {
        "V_GR_PS": sp.simplify(V_GR_ser.subs({z: 0, ell: ell_val})),
        "V_proxy_PS": sp.simplify(V_proxy_ser.subs({z: 0, ell: ell_val})),
        "V_action_PS": sp.simplify(V_action_ser.subs({z: 0, ell: ell_val})),
    }
    return rows, Vps


@dataclass
class QNMOutputs:
    omega_GR: complex | None = None
    omega_gold: complex = complex(0.373672, -0.088962)
    err_GR: float | None = None
    omega_proxy: complex | None = None
    omega_action: complex | None = None
    proxy_shift: float | None = None
    action_shift: float | None = None
    action_vs_proxy: float | None = None
    t_GR: np.ndarray | None = None
    sig_GR: np.ndarray | None = None
    t_proxy: np.ndarray | None = None
    sig_proxy: np.ndarray | None = None
    t_action: np.ndarray | None = None
    sig_action: np.ndarray | None = None
    rs_grid_STAM: np.ndarray | None = None
    V_proxy: np.ndarray | None = None
    V_action: np.ndarray | None = None
    V_fcorr: np.ndarray | None = None
    rs_grid_GR: np.ndarray | None = None
    V_GR: np.ndarray | None = None


def rel_shift(a, b):
    return abs(a - b) / abs(b)


def run_time_domain(args) -> QNMOutputs:
    out = QNMOutputs()
    ell_val = args.ell
    rs_pulse = args.pulse_center
    sigma_pulse = args.pulse_sigma
    r_obs_star = args.observer
    main_window = (args.fit_start, args.fit_end)
    T_final = args.t_final

    print("=" * 96)
    print("G95a: Schwarzschild calibration")
    print("=" * 96)
    rs_min_GR, rs_max_GR, N_GR = -60.0, args.rstar_max, args.N_gr
    rs_grid_GR = np.linspace(rs_min_GR, rs_max_GR, N_GR)
    dr_GR = rs_grid_GR[1] - rs_grid_GR[0]
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rs_grid_GR])
    V_GR = np.array([V_RW_proxy(r, ell_val, k_GR) for r in r_grid_GR])
    psi0_GR = np.exp(-(rs_grid_GR - rs_pulse)**2 / (2.0 * sigma_pulse**2))
    i_obs_GR = int(np.argmin(np.abs(rs_grid_GR - r_obs_star)))
    dt_GR = args.cfl * dr_GR
    t_GR, sig_GR = evolve_TD(V_GR, dr_GR, T_final, dt_GR, psi0_GR, i_obs=i_obs_GR)
    omega_GR, _ = extract_QNM(t_GR, sig_GR, *main_window)
    if omega_GR is None:
        raise RuntimeError("GR QNM extraction failed")
    err_GR = rel_shift(omega_GR, out.omega_gold)
    print(f"  omega_GR      = {omega_GR.real:.6f} {omega_GR.imag:+.6f} i")
    print(f"  omega_Leaver  = {out.omega_gold.real:.6f} {out.omega_gold.imag:+.6f} i")
    print(f"  calibration error = {100*err_GR:.4f}%")

    out.omega_GR = omega_GR
    out.err_GR = err_GR
    out.t_GR, out.sig_GR = t_GR, sig_GR
    out.rs_grid_GR, out.V_GR = rs_grid_GR, V_GR

    if err_GR > args.calib_gate:
        print("  FAIL: calibration gate not passed; STAM values will not be treated as predictions.")
        return out
    print("  PASS")

    print("=" * 96)
    print("G95b: STAM proxy + full axial constrained-action candidate")
    print("=" * 96)
    rs_min_STAM, rs_max_STAM, N_STAM = args.rstar_min, args.rstar_max, args.N_stam
    rs_grid_STAM = np.linspace(rs_min_STAM, rs_max_STAM, N_STAM)
    dr_STAM = rs_grid_STAM[1] - rs_grid_STAM[0]
    print(f"  STAM r* grid: [{rs_min_STAM}, {rs_max_STAM}], N={N_STAM}, dr*={dr_STAM:.5f}")
    r_grid_STAM = build_r_of_rstar_STAM(rs_grid_STAM, max_step_in=args.ode_max_step)
    V_proxy = np.array([V_RW_proxy(r, ell_val, k_STAM) for r in r_grid_STAM])
    f_grid = f_grid_on_r(r_grid_STAM, f_norm=args.f_norm)
    V_action, V_fcorr = build_action_aware_potential(V_proxy, f_grid, dr_STAM)

    psi0_STAM = np.exp(-(rs_grid_STAM - rs_pulse)**2 / (2.0 * sigma_pulse**2))
    i_obs_STAM = int(np.argmin(np.abs(rs_grid_STAM - r_obs_star)))
    dt_STAM = args.cfl * dr_STAM

    print("  evolving proxy...")
    t_proxy, sig_proxy = evolve_TD(V_proxy, dr_STAM, T_final, dt_STAM, psi0_STAM, i_obs=i_obs_STAM)
    omega_proxy, _ = extract_QNM(t_proxy, sig_proxy, *main_window)
    print("  evolving action-aware...")
    t_action, sig_action = evolve_TD(V_action, dr_STAM, T_final, dt_STAM, psi0_STAM, i_obs=i_obs_STAM)
    omega_action, _ = extract_QNM(t_action, sig_action, *main_window)

    out.omega_proxy = omega_proxy
    out.omega_action = omega_action
    out.t_proxy, out.sig_proxy = t_proxy, sig_proxy
    out.t_action, out.sig_action = t_action, sig_action
    out.rs_grid_STAM = rs_grid_STAM
    out.V_proxy, out.V_action, out.V_fcorr = V_proxy, V_action, V_fcorr

    if omega_proxy is not None:
        out.proxy_shift = rel_shift(omega_proxy, omega_GR)
    if omega_action is not None:
        out.action_shift = rel_shift(omega_action, omega_GR)
    if omega_action is not None and omega_proxy is not None:
        out.action_vs_proxy = rel_shift(omega_action, omega_proxy)

    print(f"  omega_proxy   = {omega_proxy.real:.6f} {omega_proxy.imag:+.6f} i")
    print(f"  omega_action  = {omega_action.real:.6f} {omega_action.imag:+.6f} i")
    print(f"  proxy shift vs GR    = {100*out.proxy_shift:.4f}%")
    print(f"  action shift vs GR   = {100*out.action_shift:.4f}%")
    print(f"  action vs proxy diff = {100*out.action_vs_proxy:.4f}%")
    return out


def write_outputs(args, rows, Vps, qnm: QNMOutputs | None):
    # CSV derivative table.
    csv_path = RESULTS / "G95_full_axial_derivative_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Plot.
    plot_path = PLOTS / "G95_full_axial_constrained_sigma.png"
    if qnm and qnm.V_proxy is not None and qnm.V_action is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        ax = axes[0, 0]
        ax.plot(qnm.rs_grid_GR, qnm.V_GR, label="GR RW", lw=2)
        ax.plot(qnm.rs_grid_STAM, qnm.V_proxy, label="STAM proxy", lw=2)
        ax.plot(qnm.rs_grid_STAM, qnm.V_action, label="STAM full axial candidate", lw=2)
        ax.set_xlim(-60, 120)
        ax.set_xlabel("r* / M")
        ax.set_ylabel("V(r*)")
        ax.set_title("Axial potentials")
        ax.grid(True, alpha=0.3)
        ax.legend()

        ax = axes[0, 1]
        ax.plot(qnm.rs_grid_STAM, qnm.V_fcorr, color="tab:purple", lw=2)
        ax.set_xlim(-60, 120)
        ax.set_xlabel("r*_STAM / M")
        ax.set_ylabel(r"$(\sqrt f)''/\sqrt f$")
        ax.set_title("f(Σ)R canonical correction")
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        ax.plot(qnm.t_GR, np.abs(qnm.sig_GR), label="GR", lw=1)
        ax.plot(qnm.t_proxy, np.abs(qnm.sig_proxy), label="STAM proxy", lw=1)
        ax.plot(qnm.t_action, np.abs(qnm.sig_action), label="STAM action", lw=1)
        ax.axvspan(args.fit_start, args.fit_end, color="gray", alpha=0.15, label="fit window")
        ax.set_yscale("log")
        ax.set_xlabel("t / M")
        ax.set_ylabel("|ψ|")
        ax.set_title("Time-domain signals")
        ax.grid(True, alpha=0.3)
        ax.legend()

        ax = axes[1, 1]
        ax.axis("off")
        lines = [
            "G95 full axial constrained-Σ attack",
            "",
            "Axial odd sector:",
            "  δΣ constrained / non-propagating",
            "  λ1, λ2 give no independent axial scalar DOF",
            "  f(Σ)R gives canonical correction",
            "",
            f"V_GR(PS)     = {Vps['V_GR_PS']}",
            f"V_proxy(PS)  = {Vps['V_proxy_PS']}",
            f"V_action(PS) = {Vps['V_action_PS']}",
        ]
        if qnm and qnm.omega_GR is not None:
            lines += [
                "",
                f"ω_GR      = {qnm.omega_GR.real:.6f} {qnm.omega_GR.imag:+.6f} i",
                f"calib     = {100*qnm.err_GR:.4f}%",
            ]
        if qnm and qnm.omega_action is not None:
            lines += [
                f"ω_proxy   = {qnm.omega_proxy.real:.6f} {qnm.omega_proxy.imag:+.6f} i",
                f"ω_action  = {qnm.omega_action.real:.6f} {qnm.omega_action.imag:+.6f} i",
                f"proxy shift  = {100*qnm.proxy_shift:.4f}%",
                f"action shift = {100*qnm.action_shift:.4f}%",
            ]
        ax.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=9.3)
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200)
        plt.close()

    # Summary.
    md = []
    md.append("# G95 — Full axial constrained-shell-count perturbation attack\n")
    md.append("## Purpose\n")
    md.append("Attack the full axial perturbation question for the constrained shell-count action. The calculation makes the axial-sector reduction explicit and then tests the resulting action-aware axial master potential.\n")
    md.append("## Action\n")
    md.append("```text\n")
    md.append("S[Σ,g,λ1,λ2] = (1/16πG) ∫ d^4x sqrt(-g) [\n")
    md.append("    f(Σ) R + λ1((∇Σ)^2 - W) + λ2(u^μ ∂_μΣ) - 2V(Σ)\n")
    md.append("]\n")
    md.append("```\n")
    md.append("with `A = Σ/3`.\n")
    md.append("## Axial-sector reduction\n")
    md.append("Odd-parity axial perturbations do not carry an independent shell-count perturbation. The double-LM constraints remove independent `δΣ` propagation, and the λ-sector supplies no propagating axial scalar. The remaining leading action effect is the variable gravitational stiffness `f(Σ)R`. Canonical normalization gives:\n")
    md.append("```text\n")
    md.append("V_full_axial = V_proxy + (sqrt(f))'' / sqrt(f)\n")
    md.append("```\n")
    md.append("where primes are derivatives in the STAM tortoise coordinate.\n")
    md.append("## Photon-sphere derivative structure\n")
    md.append(f"- `V_GR(3M)     = {Vps['V_GR_PS']}`\n")
    md.append(f"- `V_proxy(3M)  = {Vps['V_proxy_PS']}`\n")
    md.append(f"- `V_action(3M) = {Vps['V_action_PS']}`\n\n")
    for comp in ["proxy_minus_GR", "action_minus_GR"]:
        sub = [r for r in rows if r["comparison"] == comp]
        first_r = next((r for r in sub if abs(r["r_derivative_diff_float"]) > 1e-14), None)
        first_s = next((r for r in sub if abs(r["rstar_derivative_diff_float"]) > 1e-14), None)
        md.append(f"### {comp}\n")
        if first_r:
            md.append(f"- first nonzero ordinary r derivative: order {first_r['derivative_order']}, value `{first_r['r_derivative_diff']}`\n")
        if first_s:
            md.append(f"- first nonzero tortoise r* derivative: order {first_s['derivative_order']}, value `{first_s['rstar_derivative_diff']}`\n")
        md.append("\n")
    if qnm and qnm.omega_GR is not None:
        md.append("## Time-domain QNM diagnostic\n")
        md.append(f"- `omega_GR      = {qnm.omega_GR.real:.6f} {qnm.omega_GR.imag:+.6f} i`\n")
        md.append(f"- `omega_Leaver  = {qnm.omega_gold.real:.6f} {qnm.omega_gold.imag:+.6f} i`\n")
        md.append(f"- calibration error = **{100*qnm.err_GR:.4f}%**\n")
        if qnm.omega_action is not None:
            md.append(f"- `omega_proxy   = {qnm.omega_proxy.real:.6f} {qnm.omega_proxy.imag:+.6f} i`\n")
            md.append(f"- `omega_action  = {qnm.omega_action.real:.6f} {qnm.omega_action.imag:+.6f} i`\n")
            md.append(f"- proxy shift vs GR = **{100*qnm.proxy_shift:.4f}%**\n")
            md.append(f"- full axial/action-aware shift vs GR = **{100*qnm.action_shift:.4f}%**\n")
            md.append(f"- action vs proxy difference = **{100*qnm.action_vs_proxy:.4f}%**\n")
    md.append("## Status\n")
    md.append("G95 promotes the G92/G93 correction from an informal add-on to the leading axial consequence of the constrained shell-count action after the scalar shell-count mode is removed. The true gold standard remains a hand-derived quadratic axial action, but this script is the current framework-native axial master-equation candidate.\n")
    md.append("## Files\n")
    md.append("- `scripts/G95_full_axial_constrained_sigma.py`\n")
    md.append("- `results/G95_full_axial_derivative_table.csv`\n")
    md.append("- `plots/G95_full_axial_constrained_sigma.png`\n")

    summary_path = RESULTS / "G95_full_axial_constrained_sigma_summary.md"
    summary_path.write_text("\n".join(md), encoding="utf-8")

    return summary_path, csv_path, plot_path if (qnm and qnm.V_proxy is not None) else None


def parse_args():
    p = argparse.ArgumentParser(description="G95 full axial constrained-Σ perturbation attack")
    p.add_argument("--derive-only", action="store_true", help="only run symbolic/series derivative analysis")
    p.add_argument("--quick", action="store_true", help="shorter TD run for smoke testing")
    p.add_argument("--ell", type=int, default=2)
    p.add_argument("--rstar-min", type=float, default=-500.0)
    p.add_argument("--rstar-max", type=float, default=400.0)
    p.add_argument("--N-gr", type=int, default=6001)
    p.add_argument("--N-stam", type=int, default=10001)
    p.add_argument("--t-final", type=float, default=350.0)
    p.add_argument("--fit-start", type=float, default=100.0)
    p.add_argument("--fit-end", type=float, default=250.0)
    p.add_argument("--pulse-center", type=float, default=30.0)
    p.add_argument("--pulse-sigma", type=float, default=3.0)
    p.add_argument("--observer", type=float, default=50.0)
    p.add_argument("--cfl", type=float, default=0.5)
    p.add_argument("--calib-gate", type=float, default=0.01)
    p.add_argument("--f-norm", type=float, default=1.0)
    p.add_argument("--ode-max-step", type=float, default=0.04)
    return p.parse_args()


def main():
    args = parse_args()
    if args.quick:
        args.N_gr = 3001
        args.N_stam = 5001
        args.t_final = 250.0
        args.rstar_min = -300.0
        args.rstar_max = 250.0
        args.fit_start = 80.0
        args.fit_end = 200.0

    print("=" * 100)
    print("G95: full axial constrained-shell-count perturbation attack")
    print("=" * 100)
    print()
    print("Axial reduction logic:")
    print("  δΣ is constrained/non-propagating by λ1, λ2.")
    print("  Odd-parity sector carries only tensor/graviton perturbations.")
    print("  Variable gravitational stiffness f(Σ) canonically shifts V by (sqrt(f))''/sqrt(f).")
    print()

    rows, Vps = local_derivative_table(ell_val=args.ell)
    print("Photon-sphere potential values:")
    for k, v in Vps.items():
        print(f"  {k}: {v} = {float(sp.N(v)):.10f}")
    print()
    print("First nonzero derivative differences:")
    for comp in ["proxy_minus_GR", "action_minus_GR"]:
        sub = [r for r in rows if r["comparison"] == comp]
        first_r = next((r for r in sub if abs(r["r_derivative_diff_float"]) > 1e-14), None)
        first_s = next((r for r in sub if abs(r["rstar_derivative_diff_float"]) > 1e-14), None)
        print(f"  {comp}:")
        print(f"    ordinary r : order {first_r['derivative_order']}, value {first_r['r_derivative_diff']}")
        print(f"    tortoise r*: order {first_s['derivative_order']}, value {first_s['rstar_derivative_diff']}")
    print()

    qnm = None
    if not args.derive_only:
        qnm = run_time_domain(args)

    summary_path, csv_path, plot_path = write_outputs(args, rows, Vps, qnm)
    print()
    print(f"Summary written: {summary_path}")
    print(f"Derivative table: {csv_path}")
    if plot_path:
        print(f"Plot written: {plot_path}")


if __name__ == "__main__":
    main()
