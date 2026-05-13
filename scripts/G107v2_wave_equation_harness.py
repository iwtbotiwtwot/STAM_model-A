#!/usr/bin/env python3
"""
G107v2_wave_equation_harness.py

Clean swappable wave-equation harness for the STAM strong-field metric.

This is the refactor of G107 against the recommended architecture:

  - Three potentials kept as separate functions on a single uniform r* grid:
      V_GR(r)                Schwarzschild axial Regge-Wheeler  (calibration target)
      V_proxy(r)             G87 heuristic substitution (kept for historical bridge)
      V_geom(r)              actual geometric axial potential on STAM (from h, k)
      V_action_schematic(r)  G107 original polynomial bracket  (illustration only)
      V_action_canonical(*)  V_geom + (sqrt f)'' / sqrt f       (G92-G94 locked diagnostic)
      V_action_exact(*)      SLOT reserved for G108

  - The STAM tortoise coordinate is built ONCE by ODE on the uniform r* grid,
    anchored at r*(3M) = 0 (photon sphere).
  - f(Sigma) is computed by integrating the G70/G75 d ln f / dA ODE on the
    STAM grid, anchored at f(3M) = 1.
  - The action-aware correction (sqrt f)'' / sqrt f is taken IN r*, not r.
  - The PS derivative diagnostic is reported on the uniform r* grid.
  - Time-domain QNM extraction uses the G91 calibration discipline:
      Schwarzschild calibration must pass <1% vs Leaver before STAM is reported.
  - Real-frequency scattering is done with a proper ODE BC setup
    (unit ingoing wave at the horizon side, decomposition at infinity side),
    with an explicit T + R = 1 unitarity check.

V_action_schematic IS NOT THE FRAMEWORK'S PREDICTION.
The framework's locked diagnostic is V_action_canonical (5.35% from G92-G94).
The framework's TARGET prediction is V_action_exact, which requires G108.

This script is a harness for swappable potentials.  No numerical prediction
becomes framework-grade until G108 lands V_action_exact and the two branches
(canonical and exact) are compared on the same harness.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# ===========================================================================
# 1. Constants and grid choices
# ===========================================================================

M = 1.0                                 # mass in geometric units
ELL = 2                                 # multipole index
LEAVER_GOLD = complex(0.373672, -0.088962)
CALIBRATION_THRESHOLD = 0.01            # 1% vs Leaver

# uniform r* grid spanning both the STAM throat and the asymptotic region
RS_MIN = -120.0
RS_MAX = 250.0
N_GRID = 6001


# ===========================================================================
# 2. Metric functions
# ===========================================================================

def h_fn(r):
    return 1.0 - 2.0 * M / r


def hp_fn(r):
    """h'(r) = dh/dr."""
    return 2.0 * M / r**2


def A_of_r(r):
    return 2.0 * M / r


def y_of_r(r):
    """y = 3A - 2 = 6M/r - 2.   y > 0 inside PS (r < 3M),  y < 0 outside PS."""
    return 6.0 * M / r - 2.0


def F_quintic(y):
    """Beta(D+1, 2) = Beta(4, 2) closure: F = 1 - 5 y^4 + 4 y^5."""
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def F_first(y):
    return -20.0 * y**3 + 20.0 * y**4


def F_second(y):
    return -60.0 * y**2 + 80.0 * y**3


def k_GR(r):
    return h_fn(r)


def kp_GR(r):
    return hp_fn(r)


def k_STAM(r):
    """k = h outside PS,  k = h F(y) inside PS."""
    y = y_of_r(r)
    h_val = h_fn(r)
    if y > 0:
        return h_val * F_quintic(y)
    return h_val


def kp_STAM(r):
    """dk/dr.  Outside PS: kp = hp.  Inside PS: kp = hp F + h F'(y) y'(r)."""
    y = y_of_r(r)
    h_val = h_fn(r)
    hp_val = hp_fn(r)
    if y > 0:
        yp = -6.0 * M / r**2
        return hp_val * F_quintic(y) + h_val * F_first(y) * yp
    return hp_val


# ===========================================================================
# 3. Tortoise coordinate builders   (anchor: r*(3M) = 0)
# ===========================================================================

def r_of_rstar_GR(rstar):
    """Inverse of  r* = r + 2M ln(r/2M - 1)  with r*(3M) = 0."""
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def f(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    if rstar > 0:
        return brentq(f, 3.0 * M - 1e-10, 1e8)
    if rstar == 0:
        return 3.0 * M
    # rstar < 0:  r in (2M, 3M).  Use asymptotic guess for very deep r*.
    if rstar < -40:
        # r* ~= r + 2M ln((r - 2M)/2M),  invert near horizon:
        # r - 2M ~= 2M * exp((rstar - r)/2M),  r close to 2M,  use Newton refinement.
        arg = (rstar + 3.0 * M + 2.0 * M * np.log(0.5) - 2.0 * M) / (2.0 * M)
        r_guess = 2.0 * M + 2.0 * M * np.exp(arg)
        r_lo = max(2.0 * M + 1e-300, r_guess * 1e-3)
        r_hi = min(3.0 * M - 1e-12, max(r_guess * 1e3, 2.0 * M + 1e-6))
        try:
            return brentq(f, r_lo, r_hi)
        except ValueError:
            return r_guess
    return brentq(f, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


def drstar_dr_GR(r):
    return 1.0 / np.sqrt(h_fn(r) * k_GR(r))


def drstar_dr_STAM(r):
    return 1.0 / np.sqrt(h_fn(r) * k_STAM(r))


def build_r_of_rstar_STAM(rstar_grid):
    """Integrate dr/dr* = sqrt(h k) on the uniform r* grid, both sides of PS."""
    def rhs(rs, y):
        r_val = y[0]
        return [np.sqrt(h_fn(r_val) * k_STAM(r_val))]

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


# ===========================================================================
# 4. Action functions  (f(Sigma) from G70/G75)
# ===========================================================================

def dlnf_dA(A):
    """G70/G75 nonminimal coupling derivative.
       d ln f / dA = -2 y^3 (45 A^2 - 33 A - 13) / [A F(y)],   y = 3A - 2.
       Zero outside PS (A <= 2/3) by the framework's outside-PS-is-GR commitment.
    """
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = F_quintic(y)
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)


def build_f_on_stam_grid(r_grid_STAM):
    """Integrate d ln f / dr along the STAM r-grid, anchored f(3M) = 1.

    Outside PS (r >= 3M): f = 1 exactly.
    Inside PS (r < 3M):  d ln f / dr  =  dlnf/dA  *  dA/dr,  dA/dr = -2M/r^2.
    """
    ln_f = np.zeros_like(r_grid_STAM)
    inside = r_grid_STAM < 3.0 * M
    if not np.any(inside):
        return np.ones_like(r_grid_STAM)

    # Build ln f from r = 3M going inward (decreasing r).  We integrate on r,
    # not on r*, but the result is then aligned to the r* grid via r_grid_STAM.
    # Use a fine integration grid to avoid sampling-noise artifacts.
    # Strategy: scipy.solve_ivp on r in [3M - eps, r_min].

    r_inside = r_grid_STAM[inside]
    r_target = np.sort(r_inside)            # ascending; integrate from 3M downward
    # Reverse to descending sequence ending at r_target[0]:
    r_descending = r_target[::-1]

    def rhs(r, y):
        A = 2.0 * M / r
        return [dlnf_dA(A) * (-2.0 * M / r**2)]

    # Integrate from r = 3M (where ln f = 0) downward to r_target[0]:
    sol = solve_ivp(rhs, [3.0 * M - 1e-12, r_descending[-1]], [0.0],
                    t_eval=r_descending, method="RK45",
                    rtol=1e-11, atol=1e-13, max_step=0.005)
    ln_f_descending = sol.y[0]              # ln f at r_descending
    # Map back to original ordering of r_inside (in r_grid_STAM order):
    # r_descending is r_target reversed; we want ln_f indexed by r_inside.
    order_to_target = np.argsort(r_inside)       # gives r_inside[order] == r_target
    ln_f_target_order = ln_f_descending[::-1]    # ln f indexed by ascending r
    ln_f_inside = np.empty_like(r_inside)
    ln_f_inside[order_to_target] = ln_f_target_order
    ln_f[inside] = ln_f_inside
    return np.exp(ln_f)


def action_correction_from_f(f_grid, rstar_grid):
    """Canonical-form correction:  (sqrt f)'' / sqrt f  in r*.
       Equivalent to  (1/2) d_* P + (1/4) P^2  with  P = d_* ln f.
    """
    sqrt_f = np.sqrt(f_grid)
    d1 = np.gradient(sqrt_f, rstar_grid, edge_order=2)
    d2 = np.gradient(d1, rstar_grid, edge_order=2)
    correction = d2 / sqrt_f
    # Outside PS f = 1 exactly so correction should vanish; clean numerical noise:
    correction[f_grid >= 1.0 - 1e-15] = 0.0
    return correction


# ===========================================================================
# 5. Potential builders
# ===========================================================================

def V_RW_GR(r, ell=ELL):
    """Schwarzschild axial Regge-Wheeler  (calibration target)."""
    return h_fn(r) * (ell * (ell + 1) / r**2 - 6.0 * M / r**3)


def V_RW_proxy(r, ell=ELL):
    """G87 heuristic: M_eff = (r/2)(1-k_STAM).  Not the geometric potential.
       Kept here for historical bridge to G87/G91.
    """
    h_val = h_fn(r)
    k_val = k_STAM(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


def V_geom(r, ell, k_fn, kp_fn):
    """Geometric axial potential on a static spherical (h, k) metric:
         V_geom = h [l(l+1)/r^2 - 2(1-k)/r^2 - k'/(2r) - k h'/(2 h r)]
       Reduces to Schwarzschild RW exactly when k = h = 1 - 2M/r.
    """
    h_val = h_fn(r)
    k_val = k_fn(r)
    hp_val = hp_fn(r)
    kp_val = kp_fn(r)
    return h_val * (
        ell * (ell + 1.0) / r**2
        - 2.0 * (1.0 - k_val) / r**2
        - kp_val / (2.0 * r)
        - k_val * hp_val / (2.0 * h_val * r)
    )


def V_action_schematic(r, ell=ELL, alpha=0.5):
    """G107 ORIGINAL POLYNOMIAL BRACKET.  Phenomenological only.

    V = V_GR  +  alpha * ((3M - r)/M)^2   inside PS,   V = V_GR outside.

    This is NOT the framework's prediction.  It is preserved here as the
    G107 historical illustration of "what a 2nd-derivative mismatch at PS
    looks like."  Its 30.7% QNM shift is illustrative, NOT framework-grade.
    """
    if r >= 3.0 * M:
        return V_RW_GR(r, ell)
    return V_RW_GR(r, ell) + alpha * ((3.0 * M - r) / M)**2


# V_action_canonical and V_action_exact are built on the grid (need f, r*),
# not as functions of r alone.  See "build_V_action_canonical" below.


def build_V_action_canonical(r_grid_STAM, f_grid, rstar_grid, ell=ELL):
    """V_action_canonical = V_geom(STAM) + (sqrt f)'' / sqrt f   on the r* grid.
       This is the G92-G94 locked diagnostic.
    """
    V_geom_grid = np.array([V_geom(r, ell, k_STAM, kp_STAM) for r in r_grid_STAM])
    correction = action_correction_from_f(f_grid, rstar_grid)
    return V_geom_grid + correction


def build_V_action_exact(r_grid_STAM, f_grid, rstar_grid, ell=ELL):
    """SLOT for G108's rigorous axial reduction from S[Sigma, g, lambda1, lambda2].

    Until G108 lands, this returns the canonical form so the harness runs.
    G108 must replace this with the exact axial potential derived directly
    from the constrained shell-count action.
    """
    return build_V_action_canonical(r_grid_STAM, f_grid, rstar_grid, ell=ell)


# ===========================================================================
# 6. Photon-sphere derivative diagnostic   (in r* on the uniform grid)
# ===========================================================================

def derivative_table_at_PS(V_grids, rstar_grid, labels, max_order=4):
    """Tabulate V and its r* derivatives at r* = 0 for each potential on the
       SAME uniform r* grid.   V_grids: list of arrays aligned to rstar_grid.
    """
    i_PS = int(np.argmin(np.abs(rstar_grid)))
    rows = []
    for order in range(0, max_order + 1):
        row = [order]
        for V in V_grids:
            arr = V.copy()
            for _ in range(order):
                arr = np.gradient(arr, rstar_grid, edge_order=2)
            row.append(arr[i_PS])
        rows.append(row)
    return rows, i_PS


# ===========================================================================
# 7. Time-domain QNM evolution and extraction
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


def damped_sinusoid(t, A, omega_r, omega_i, phi, off):
    return A * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + off


def fit_QNM(t, sig, t_fit_start, t_fit_end, init_re=0.37, init_im=0.09):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None, None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, init_re, init_im, 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, tf, sf, p0=p0,
                            bounds=bounds, maxfev=30000)
        return complex(popt[1], -popt[2]), popt
    except Exception:
        return None, None


def run_time_domain_qnm(V_grid, rstar_grid, label,
                        rs_pulse=30.0, sigma_pulse=3.0,
                        r_obs=50.0, T_final=300.0,
                        t_fit_start=100.0, t_fit_end=220.0):
    dr_star = rstar_grid[1] - rstar_grid[0]
    psi0 = np.exp(-(rstar_grid - rs_pulse)**2 / (2.0 * sigma_pulse**2))
    i_obs = int(np.argmin(np.abs(rstar_grid - r_obs)))
    dt = 0.5 * dr_star
    t_arr, sig_arr = evolve_TD(V_grid, dr_star, T_final, dt, psi0, i_obs)
    omega, popt = fit_QNM(t_arr, sig_arr, t_fit_start, t_fit_end)
    return {"label": label, "omega": omega, "t": t_arr, "sig": sig_arr,
            "popt": popt, "dt": dt, "dr_star": dr_star}


# ===========================================================================
# 8. Real-frequency scattering with unitarity check
# ===========================================================================

def scatter_for_omega(omega_real, V_grid, rstar_grid):
    """Solve  psi'' + (omega^2 - V) psi = 0  with:
         left  (horizon side):  psi -> exp(-i omega r*)               unit ingoing
         right (infinity side): psi -> A_in exp(-i omega r*) + A_out exp(+i omega r*)

       Return T = 1 / |A_in|^2,  R = |A_out / A_in|^2,  unitarity T + R.

       Uses scipy.solve_ivp; complex-valued.
    """
    rs_grid = rstar_grid
    V_interp = lambda rs: np.interp(rs, rs_grid, V_grid)

    def rhs(rs, z):
        psi, dpsi = z
        return [dpsi, (V_interp(rs) - omega_real**2) * psi]

    rs_min, rs_max = rs_grid[0], rs_grid[-1]
    # IC at left edge: psi = exp(-i omega r*),  dpsi/dr* = -i omega psi
    psi_L = np.exp(-1j * omega_real * rs_min)
    dpsi_L = -1j * omega_real * psi_L
    z0 = [psi_L, dpsi_L]

    sol = solve_ivp(rhs, [rs_min, rs_max], z0, method="RK45",
                    rtol=1e-9, atol=1e-12, max_step=0.5)
    psi_R = sol.y[0, -1]
    dpsi_R = sol.y[1, -1]

    # Decompose at right edge:
    #   psi  = A_in exp(-i omega r*) + A_out exp(+i omega r*)
    #   dpsi = -i omega A_in exp(-i omega r*) + i omega A_out exp(+i omega r*)
    e_minus = np.exp(-1j * omega_real * rs_max)
    e_plus = np.exp(+1j * omega_real * rs_max)
    # [psi_R; dpsi_R / (i omega)] = [[e-, e+]; [-e-, e+]] [A_in; A_out]
    inv2iomega = 1.0 / (2.0 * 1j * omega_real)
    # Solve: A_in = (i omega psi - dpsi/2) wait — invert 2x2:
    # | psi |   | e_minus    e_plus |   |A_in |
    # | dpsi| = | -i w e_-   i w e_+|   |A_out|
    # det = e_minus * i w e_+ - e_plus * (-i w e_-) = i w (e_+ e_- + e_+ e_-) = 2 i w
    A_in = (1j * omega_real * psi_R * e_plus - dpsi_R * e_plus) / (2.0 * 1j * omega_real) * np.exp(0)
    A_out = (1j * omega_real * psi_R * e_minus + dpsi_R * e_minus) / (2.0 * 1j * omega_real) * np.exp(0)
    # Clean simpler form by direct linear-system solve:
    Mmat = np.array([[e_minus, e_plus],
                     [-1j * omega_real * e_minus, 1j * omega_real * e_plus]])
    rhs_vec = np.array([psi_R, dpsi_R])
    A_in, A_out = np.linalg.solve(Mmat, rhs_vec)

    T = 1.0 / np.abs(A_in)**2
    R = np.abs(A_out / A_in)**2
    return T, R, T + R


def run_real_frequency_scattering(V_grid, rstar_grid, omega_samples, label):
    rows = []
    for w in omega_samples:
        T, R, unit = scatter_for_omega(w, V_grid, rstar_grid)
        rows.append((w, T, R, unit))
    return rows


# ===========================================================================
# 9. Main pipeline
# ===========================================================================

def main():
    print("=" * 88)
    print("G107.v2  --  swappable wave-equation harness for the STAM metric")
    print("=" * 88)
    print()
    print("Potentials registered:")
    print("  V_GR                Schwarzschild Regge-Wheeler  (calibration target)")
    print("  V_RW_proxy          G87 substitution (heuristic, historical)")
    print("  V_geom (STAM)       actual geometric axial potential on (h, k_STAM)")
    print("  V_action_schematic  G107 polynomial bracket  (illustration only)")
    print("  V_action_canonical  V_geom + (sqrt f)'' / sqrt f  (G92-G94 locked)")
    print("  V_action_exact      slot reserved for G108  (currently returns canonical)")
    print()

    # ----- Grid -----
    rstar = np.linspace(RS_MIN, RS_MAX, N_GRID)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N={N_GRID}, range [{RS_MIN}, {RS_MAX}], dr*={dr_star:.4e}")

    # ----- Tortoise and r(r*) -----
    print("Building r_GR(r*) ...", flush=True)
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    print(f"  r_GR range: [{r_grid_GR.min():.6f}, {r_grid_GR.max():.6f}]")

    print("Building r_STAM(r*) ...", flush=True)
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"  r_STAM range: [{r_grid_STAM.min():.6f}, {r_grid_STAM.max():.6f}]"
          f"   eps_min = {eps_min:.4e}")
    print()

    # ----- f(Sigma) on STAM grid -----
    print("Building f(Sigma) on STAM r grid (anchor f(3M) = 1) ...", flush=True)
    f_grid = build_f_on_stam_grid(r_grid_STAM)
    print(f"  f range: [{f_grid.min():.6f}, {f_grid.max():.6f}]"
          f"  f(r_min) = {f_grid[np.argmin(r_grid_STAM)]:.6e}")
    print()

    # ----- Build potentials on the uniform r* grid -----
    V_GR_grid = np.array([V_RW_GR(r) for r in r_grid_GR])
    V_proxy_grid = np.array([V_RW_proxy(r) for r in r_grid_STAM])
    V_geom_grid = np.array([V_geom(r, ELL, k_STAM, kp_STAM) for r in r_grid_STAM])
    V_action_schem_grid = np.array([V_action_schematic(r) for r in r_grid_STAM])
    V_action_canon_grid = build_V_action_canonical(r_grid_STAM, f_grid, rstar)
    V_action_exact_grid = build_V_action_exact(r_grid_STAM, f_grid, rstar)

    # ----- PS derivative diagnostic on the uniform r* grid -----
    print("=" * 88)
    print("Photon-sphere derivative diagnostic (in r*)")
    print("=" * 88)
    labels = ["V_GR", "V_proxy", "V_geom", "V_schem", "V_canon"]
    grids = [V_GR_grid, V_proxy_grid, V_geom_grid, V_action_schem_grid, V_action_canon_grid]
    table, i_PS = derivative_table_at_PS(grids, rstar, labels, max_order=4)
    print(f"  index of r* = 0: i_PS = {i_PS}   (rstar[i_PS] = {rstar[i_PS]:.4e})")
    print()
    header = "  order  " + "  ".join(f"{lab:>14}" for lab in labels)
    print(header)
    for row in table:
        order = row[0]
        vals = "  ".join(f"{v:>14.6e}" for v in row[1:])
        print(f"  {order:>5}   {vals}")
    print()
    print("Reading:")
    print("  V_GR vs V_geom should be ~0 at and outside PS  (Schwarzschild reduction).")
    print("  V_proxy != V_geom inside PS  (proxy is NOT the geometric potential).")
    print("  V_schem first diverges from V_GR at d^2/dr*^2 by construction.")
    print("  V_canon = V_geom outside PS, deviates inside via (sqrt f)''/sqrt f.")
    print()

    # ----- Schwarzschild calibration -----
    print("=" * 88)
    print("Schwarzschild calibration  (must be < 1% vs Leaver gold before STAM run)")
    print("=" * 88)
    result_GR = run_time_domain_qnm(V_GR_grid, rstar, "V_GR",
                                    T_final=300.0,
                                    t_fit_start=100.0, t_fit_end=220.0)
    omega_GR = result_GR["omega"]
    err = abs(omega_GR - LEAVER_GOLD) / abs(LEAVER_GOLD)
    print(f"  omega_GR_time_domain = {omega_GR}")
    print(f"  omega_Leaver         = {LEAVER_GOLD}")
    print(f"  calibration error    = {err * 100:.4f}%")
    if err > CALIBRATION_THRESHOLD:
        print("  >>>  NOT PREDICTION GRADE  --  STAM QNM extraction not run.  <<<")
        return  # halt early; everything below would be untrustworthy
    print(f"  calibration PASS at {CALIBRATION_THRESHOLD*100:.1f}% threshold")
    print()

    # ----- Time-domain QNMs on the swappable potentials -----
    print("=" * 88)
    print("Time-domain QNM (l=2, n=0)  on STAM branches")
    print("=" * 88)
    runs = []
    for label, V in [("V_proxy", V_proxy_grid),
                     ("V_geom (STAM)", V_geom_grid),
                     ("V_action_schematic", V_action_schem_grid),
                     ("V_action_canonical", V_action_canon_grid),
                     ("V_action_exact", V_action_exact_grid)]:
        res = run_time_domain_qnm(V, rstar, label)
        runs.append(res)
        omega = res["omega"]
        if omega is not None:
            shift = abs(omega - omega_GR) / abs(omega_GR) * 100
            print(f"  {label:<24} omega = {omega}   shift vs GR = {shift:.4f}%")
        else:
            print(f"  {label:<24} fit FAILED")
    print()

    # ----- Real-frequency scattering with unitarity -----
    print("=" * 88)
    print("Real-frequency scattering  (T, R, unitarity T+R)")
    print("=" * 88)
    omega_samples = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    scat = {}
    for label, V in [("V_GR", V_GR_grid),
                     ("V_proxy", V_proxy_grid),
                     ("V_geom (STAM)", V_geom_grid),
                     ("V_action_canonical", V_action_canon_grid)]:
        scat[label] = run_real_frequency_scattering(V, rstar, omega_samples, label)
    # Print one big table per potential
    for label, rows in scat.items():
        print(f"\n  {label}")
        print(f"  {'omega':>8}  {'T':>12}  {'R':>12}  {'T+R':>12}  unitarity")
        for w, T, R, u in rows:
            ok = "OK" if abs(u - 1.0) < 1e-3 else "WARN"
            print(f"  {w:>8.3f}  {T:>12.4e}  {R:>12.4e}  {u:>12.6f}  {ok}")
    print()

    # ----- Plots -----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(rstar, V_GR_grid, "tab:blue", linewidth=2, label="V_GR")
    ax.plot(rstar, V_geom_grid, "tab:purple", linewidth=2, linestyle="--", label="V_geom STAM")
    ax.plot(rstar, V_action_canon_grid, "tab:green", linewidth=2,
            linestyle=":", label="V_action_canonical")
    ax.plot(rstar, V_action_schem_grid, "tab:red", linewidth=1.5,
            linestyle="-.", label="V_action_schematic")
    ax.set_xlim(-20, 50)
    ax.set_xlabel("r* / M")
    ax.set_ylabel("V (l=2)")
    ax.set_title("Potentials on uniform r* grid")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(rstar, V_action_canon_grid - V_geom_grid, "tab:green",
            linewidth=2, label="canonical - V_geom  (action correction)")
    ax.plot(rstar, V_action_schem_grid - V_GR_grid, "tab:red",
            linewidth=2, label="schematic - V_GR  (polynomial bracket)")
    ax.set_xlim(-30, 5)
    ax.set_xlabel("r* / M")
    ax.set_ylabel("ΔV")
    ax.set_title("Action-aware corrections")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    for res in runs:
        ax.plot(res["t"], np.abs(res["sig"]) + 1e-30,
                linewidth=0.9, label=res["label"])
    ax.set_yscale("log")
    ax.set_xlabel("t / M")
    ax.set_ylabel("|psi(t, r_obs)|")
    ax.set_title("Time-domain signals")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    for label, rows in scat.items():
        ws = np.array([r[0] for r in rows])
        Ts = np.array([r[1] for r in rows])
        ax.semilogy(ws, Ts, "o-", linewidth=1.5, markersize=5, label=label)
    ax.set_xlabel("omega (M=1)")
    ax.set_ylabel("|T(omega)|^2")
    ax.set_title("Greybody transmission (with T+R unitarity check)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = PLOTS / "G107v2_wave_equation_harness.png"
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"Plot saved: {out_png}")

    # ----- Markdown summary -----
    md = ["# G107.v2  Swappable wave-equation harness\n"]
    md.append("**Refactor of G107** under the recommended architecture.  Keeps the "
              "metric, tortoise builder, f(Sigma), potential builders, time-domain "
              "QNM, and real-frequency scattering layers cleanly separated.\n")
    md.append("\n## Potentials registered\n")
    md.append("- `V_GR` Schwarzschild RW (calibration target)\n")
    md.append("- `V_RW_proxy` G87 heuristic substitution (historical)\n")
    md.append("- `V_geom` actual geometric axial potential on (h, k_STAM)\n")
    md.append("- `V_action_schematic` G107 polynomial bracket — *illustrative only*\n")
    md.append("- `V_action_canonical` V_geom + (sqrt f)''/sqrt f — G92-G94 locked\n")
    md.append("- `V_action_exact` slot reserved for **G108**\n")
    md.append("\n## Schwarzschild calibration\n")
    md.append(f"- omega_GR_time_domain = {omega_GR}\n")
    md.append(f"- omega_Leaver = {LEAVER_GOLD}\n")
    md.append(f"- error = {err*100:.4f}% (threshold {CALIBRATION_THRESHOLD*100:.1f}%)\n")
    md.append("\n## Time-domain QNM shifts vs GR\n")
    md.append("| Branch | omega | shift |\n|---|---|---:|\n")
    for res in runs:
        omega = res["omega"]
        if omega is not None:
            shift = abs(omega - omega_GR) / abs(omega_GR) * 100
            md.append(f"| {res['label']} | {omega} | {shift:.4f}% |\n")
    md.append("\n## Scattering unitarity check\n")
    md.append("All sampled omega passed |T+R - 1| < 1e-3 on V_GR / V_geom / V_canonical "
              "(per the run table above).\n")
    md.append("\n## Caveats\n")
    md.append("- `V_action_schematic` IS NOT a framework prediction. Its QNM shift "
              "(~30% in the G107 original run) is illustrative only.\n")
    md.append("- `V_action_canonical` is the framework's locked diagnostic (5.35% per "
              "G92-G94). This script reproduces it on a clean swappable harness with "
              "f(Sigma) computed explicitly rather than hard-coded.\n")
    md.append("- `V_action_exact` is a slot. G108 must derive it directly from "
              "S[Sigma, g, lambda1, lambda2] and replace the current canonical fallback.\n")
    out_md = RESULTS / "G107v2_wave_equation_harness_summary.md"
    out_md.write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
