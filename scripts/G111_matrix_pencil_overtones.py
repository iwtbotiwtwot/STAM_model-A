#!/usr/bin/env python3
"""
G111_matrix_pencil_overtones.py

Higher-mode QNM extraction with the MATRIX PENCIL METHOD, which is the
standard tool for clean multi-mode extraction from time-domain signals.

Why this script
===============
G110's two-mode curve_fit failed for n=1 (calibration errors 30%-85%).
The damped-sinusoid fit gets trapped in local minima when one mode has
~1000x smaller amplitude than the other.  Matrix pencil sidesteps this
by extracting all poles simultaneously via SVD + generalized eigenvalue.

Method
======
Given signal y(t_k) sampled at dt:
   1. Build Hankel matrix Y of size (N - L) x (L + 1).
   2. SVD: Y = U S V^T.  Retain top K singular values where K = 2 * (number
      of physical modes wanted), because each physical mode = pair of
      complex-conjugate poles for a real signal.
   3. Generalized eigenvalue:  z_k from V_truncated[:, 1:] = z * V_truncated[:, :-1].
      Equivalently, z_k are eigenvalues of pinv(V1) @ V2.
   4. Convert: z = exp(s dt); s = -omega_i - i omega_r in our convention;
      omega = -Im(s) + i Re(s).
   5. Pick poles inside unit disk (decay) and with omega_r > 0 (one of each
      conjugate pair).  Sort by Im(omega) to identify n=0, n=1, n=2, ...

For each ell in {2, 3, 4}:
    - Run time-domain on V_GR(ell) and V_exact(ell).
    - Matrix-pencil-extract first three poles.
    - Match by closest distance to Leaver gold.
    - Calibrate (n, 0) and (n, 1) against Leaver gold.
    - Compute V_exact shifts.

Outputs
=======
    results/G111_matrix_pencil_overtones_summary.md
    plots/G111_matrix_pencil_overtones.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# Constants and Leaver gold
M = 1.0
LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (2, 1): complex(0.346711, -0.273915),
    (3, 0): complex(0.599443, -0.092703),
    (3, 1): complex(0.582644, -0.281298),
    (4, 0): complex(0.809178, -0.094164),
    (4, 1): complex(0.796632, -0.284334),
}


# ============= Metric / V_exact (copied from G110) =============

def h_fn(r): return 1.0 - 2.0 * M / r
def hp_fn(r): return 2.0 * M / r**2
def y_of_r(r): return 6.0 * M / r - 2.0
def F_quintic(y): return 1.0 - 5.0 * y**4 + 4.0 * y**5
def F_first(y): return -20.0 * y**3 + 20.0 * y**4

def k_STAM(r):
    y = y_of_r(r)
    return h_fn(r) * F_quintic(y) if y > 0 else h_fn(r)

def kp_STAM(r):
    y = y_of_r(r)
    if y > 0:
        return hp_fn(r) * F_quintic(y) + h_fn(r) * F_first(y) * (-6.0 * M / r**2)
    return hp_fn(r)


def dlnf_dA_numeric(A):
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = 1.0 - 5.0 * y**4 + 4.0 * y**5
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)


def d2lnf_dA2_numeric(A, hstep=1e-6):
    return (dlnf_dA_numeric(A + hstep) - dlnf_dA_numeric(A - hstep)) / (2.0 * hstep)


def V_geom_grid(r_grid, ell):
    out = np.empty_like(r_grid)
    for i, rv in enumerate(r_grid):
        h_val = h_fn(rv)
        k_val = k_STAM(rv)
        out[i] = h_val * (
            ell * (ell + 1) / rv**2
            - 2.0 * (1.0 - k_val) / rv**2
            - kp_STAM(rv) / (2.0 * rv)
            - k_val * hp_fn(rv) / (2.0 * h_val * rv)
        )
    return out


def analytic_correction(r_grid_STAM):
    out = np.zeros_like(r_grid_STAM)
    for i, rv in enumerate(r_grid_STAM):
        if rv >= 3.0 * M:
            out[i] = 0.0
            continue
        A = 2.0 * M / rv
        dA_dr = -2.0 * M / rv**2
        d2A_dr2 = 4.0 * M / rv**3
        u = dlnf_dA_numeric(A)
        up = d2lnf_dA2_numeric(A)
        dlnf_dr = u * dA_dr
        d_dlnf_dr = up * (dA_dr**2) + u * d2A_dr2
        hk = h_fn(rv) * k_STAM(rv)
        if hk <= 0:
            out[i] = 0.0
            continue
        sqrt_hk = np.sqrt(hk)
        d_hk = hp_fn(rv) * k_STAM(rv) + h_fn(rv) * kp_STAM(rv)
        P = sqrt_hk * dlnf_dr
        d_sqrt_hk_dr = 0.5 * d_hk / sqrt_hk
        dP_dr = d_sqrt_hk_dr * dlnf_dr + sqrt_hk * d_dlnf_dr
        P_star = sqrt_hk * dP_dr
        out[i] = 0.5 * P_star + 0.25 * P**2
    return out


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


def build_r_of_rstar_STAM(rstar_grid):
    def rhs(rs, y):
        rv = y[0]
        return [np.sqrt(h_fn(rv) * k_STAM(rv))]
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


# ============= MATRIX PENCIL =============

def matrix_pencil_omegas(t, y, n_modes=3, L_frac=0.5):
    """Extract physical-mode complex frequencies via matrix pencil method.

    Each physical real damped sinusoid -> pair of conjugate complex poles.
    n_modes = number of physical modes wanted.  We keep K = 2 * n_modes
    complex poles, then pick those inside unit disk with omega_r > 0.

    Returns list of complex omega in convention omega = omega_r - i omega_i
    (i.e. omega.imag < 0 for decay).  Sorted by |omega_i| (n=0, n=1, ...).
    """
    N = len(y)
    L = int(L_frac * N)
    L = max(2 * n_modes + 2, min(L, N - 2 * n_modes - 2))
    rows = N - L
    cols = L + 1
    if rows < cols:
        cols = rows  # adjust if signal too short
        L = cols - 1
    Y = np.empty((rows, cols))
    for i in range(rows):
        Y[i] = y[i:i + cols]
    U, S, Vt = np.linalg.svd(Y, full_matrices=False)
    K = 2 * n_modes
    K = min(K, len(S))
    V_K = Vt[:K, :]                   # (K, cols)
    V1 = V_K[:, :-1]                  # (K, cols-1)
    V2 = V_K[:, 1:]                   # (K, cols-1)
    P = V2 @ np.linalg.pinv(V1)       # (K, K)
    z = np.linalg.eigvals(P)
    dt = t[1] - t[0]
    # Filter and convert
    z_inside = z[np.abs(z) < 1.0]
    s = np.log(z_inside) / dt
    # omega_r = -Im(s),  omega_i = -Re(s).  Our convention omega = omega_r - i omega_i.
    omegas = np.array([complex(-si.imag, si.real) for si in s])
    # Keep omega_r > 0
    omegas = omegas[omegas.real > 0]
    # Sort by decay rate (|omega.imag|) ascending = n=0, n=1, ...
    idx = np.argsort(-np.array([om.imag for om in omegas]))  # imag<0; least-negative first = least damped = n=0
    return omegas[idx]


def match_to_leaver(omegas, leaver_n0, leaver_n1):
    """Match extracted omegas to (n=0, n=1) by minimum distance to Leaver gold."""
    if len(omegas) == 0:
        return None, None
    # Best match to n=0
    d0 = np.array([abs(om - leaver_n0) for om in omegas])
    i0 = int(np.argmin(d0))
    om_n0 = omegas[i0]
    # Best match to n=1 among remaining
    remaining = np.delete(omegas, i0)
    if len(remaining) == 0:
        return om_n0, None
    d1 = np.array([abs(om - leaver_n1) for om in remaining])
    i1 = int(np.argmin(d1))
    om_n1 = remaining[i1]
    return om_n0, om_n1


# ============= Main =============

def main():
    print("=" * 88)
    print("G111  --  higher-mode QNM via matrix pencil method")
    print("=" * 88)
    print()

    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N = {N}, range [{rs_min}, {rs_max}], dr* = {dr_star:.4e}")
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"STAM eps_min = {eps_min:.4e}")
    print()

    # Initial pulse: narrow + close to PS so multiple modes excite
    pulse_center, pulse_sigma = 15.0, 1.5
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 40.0)))
    dt = 0.5 * dr_star
    T_final = 200.0

    results = []
    for ell in [2, 3, 4]:
        print("-" * 88)
        print(f"ell = {ell}")
        print("-" * 88, flush=True)

        V_GR_grid = h_fn(r_grid_GR) * (ell * (ell + 1) / r_grid_GR**2 - 6.0 * M / r_grid_GR**3)
        V_g = V_geom_grid(r_grid_STAM, ell)
        corr = analytic_correction(r_grid_STAM)
        V_exact = V_g + corr

        t_arr, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
        _, sig_exact = evolve_TD(V_exact, dr_star, T_final, dt, psi0, i_obs)

        # Strip initial transient: start pencil at t > 15 (pulse has reached observer)
        # End before noise dominates: t < 100 for n=1 to remain present
        # n=1 decays exp(-100 * 0.28) ~ 1e-12 — noise floor.
        # Use t in [15, 80] to keep both n=0 and n=1 measurable
        mask = (t_arr >= 15.0) & (t_arr <= 90.0)
        t_pen = t_arr[mask]
        sig_GR_pen = sig_GR[mask]
        sig_exact_pen = sig_exact[mask]
        print(f"  pencil window: t in [{t_pen[0]:.2f}, {t_pen[-1]:.2f}],"
              f"  N_samples = {len(t_pen)}")

        # Decimate to ~1500 samples for SVD efficiency
        if len(t_pen) > 1500:
            stride = len(t_pen) // 1500
            t_pen = t_pen[::stride]
            sig_GR_pen = sig_GR_pen[::stride]
            sig_exact_pen = sig_exact_pen[::stride]
            print(f"  decimated to N_samples = {len(t_pen)}, dt = {t_pen[1]-t_pen[0]:.4f}")

        # Matrix pencil with n_modes = 3 (look for n=0, n=1, plus one extra to absorb noise)
        omegas_GR = matrix_pencil_omegas(t_pen, sig_GR_pen, n_modes=3)
        omegas_exact = matrix_pencil_omegas(t_pen, sig_exact_pen, n_modes=3)
        print(f"  GR poles: " + ", ".join(f"{om}" for om in omegas_GR[:5]))
        print(f"  EX poles: " + ", ".join(f"{om}" for om in omegas_exact[:5]))

        om_GR_n0, om_GR_n1 = match_to_leaver(omegas_GR,
                                              LEAVER_GOLD[(ell, 0)],
                                              LEAVER_GOLD[(ell, 1)])
        om_exact_n0, om_exact_n1 = match_to_leaver(omegas_exact,
                                                    LEAVER_GOLD[(ell, 0)],
                                                    LEAVER_GOLD[(ell, 1)])

        for n, om_GR, om_exact, leaver in [
            (0, om_GR_n0, om_exact_n0, LEAVER_GOLD[(ell, 0)]),
            (1, om_GR_n1, om_exact_n1, LEAVER_GOLD[(ell, 1)]),
        ]:
            if om_GR is None or om_exact is None:
                print(f"  (ell={ell}, n={n}): extraction failed")
                results.append({"ell": ell, "n": n, "extraction": "failed"})
                continue
            calib_err = abs(om_GR - leaver) / abs(leaver)
            shift = abs(om_exact - om_GR) / abs(om_GR)
            print(f"  (ell={ell}, n={n}): omega_GR = {om_GR}  Leaver = {leaver}"
                  f"   calib = {calib_err*100:.3f}%")
            print(f"                  omega_exact = {om_exact}   shift = {shift*100:.4f}%")
            results.append({
                "ell": ell, "n": n,
                "om_GR": om_GR, "om_exact": om_exact, "leaver": leaver,
                "calib_err": calib_err, "shift": shift,
            })
        print()

    # ----- Summary table -----
    print("=" * 88)
    print("Summary")
    print("=" * 88)
    print()
    print(f"  {'ell':>4}  {'n':>2}  {'calib':>8}  "
          f"{'omega_GR':>30}  {'omega_exact':>30}  {'shift':>8}")
    print("-" * 92)
    for r in results:
        if "extraction" in r:
            print(f"  {r['ell']:>4}  {r['n']:>2}  --- EXTRACTION FAILED ---")
            continue
        print(f"  {r['ell']:>4}  {r['n']:>2}  {r['calib_err']*100:>7.3f}%  "
              f"{str(r['om_GR']):>30}  {str(r['om_exact']):>30}  "
              f"{r['shift']*100:>7.4f}%")
    print()

    # ----- Plot -----
    if results:
        fig, ax = plt.subplots(figsize=(10, 6))
        for r in results:
            if "extraction" in r:
                continue
            color = "tab:blue" if r["n"] == 0 else "tab:orange"
            marker = "o" if r["n"] == 0 else "s"
            ax.scatter(r["ell"], r["shift"]*100, c=color, marker=marker, s=100,
                       label=f"n={r['n']}" if r['ell'] == 2 else None)
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.set_xlabel("ell")
        ax.set_ylabel("V_exact shift vs GR  (%)")
        ax.set_title("Axial QNM shift by mode (matrix pencil extraction)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        out_png = PLOTS / "G111_matrix_pencil_overtones.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G111 - Higher-mode axial QNM via matrix pencil\n"]
    md.append("\nReplacement for G110's two-mode curve_fit n=1 extraction (which failed).\n")
    md.append("Matrix pencil method extracts all dominant poles simultaneously via SVD\n"
              "+ generalized eigenvalue; no local-minimum trap.\n")
    md.append("\n## Results\n")
    md.append("| ell | n | calib | omega_GR | Leaver gold | omega_exact | shift |\n")
    md.append("|---:|---:|---:|---|---|---|---:|\n")
    for r in results:
        if "extraction" in r:
            md.append(f"| {r['ell']} | {r['n']} | --- | EXTRACTION FAILED | | | |\n")
            continue
        md.append(f"| {r['ell']} | {r['n']} | {r['calib_err']*100:.3f}% | "
                  f"{r['om_GR']} | {r['leaver']} | {r['om_exact']} | "
                  f"{r['shift']*100:.4f}% |\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G111_matrix_pencil_overtones.py](../scripts/G111_matrix_pencil_overtones.py)\n")
    md.append("- [plots/G111_matrix_pencil_overtones.png](../plots/G111_matrix_pencil_overtones.png)\n")
    out_md = RESULTS / "G111_matrix_pencil_overtones_summary.md"
    out_md.write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
