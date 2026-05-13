#!/usr/bin/env python3
"""
G86_redux_action_aware_potential.py

Redux of G86 after G90-G94.

Purpose
-------
Old G86/G87 checked where STAM first departs from GR near the photon sphere.
The later chain changed the correct setup:
  - G90: correct STAM tortoise coordinate r*_STAM.
  - G91: calibrated time-domain proxy QNM after using r*_STAM.
  - G92/G93: action-aware f(Σ)R correction materially changes low-l QNM.
  - G94: action-aware shift is robust across calibrated runs.

This script does NOT extract QNM frequencies.  It explains structurally why
G92-G94 find a nonzero low-l shift by comparing three photon-sphere potentials:
  1. V_GR      : Schwarzschild axial Regge-Wheeler potential.
  2. V_proxy   : STAM tensor proxy potential from modified k(r).
  3. V_action  : action-aware diagnostic potential
                    V_proxy + (sqrt(f))''_{r*} / sqrt(f).

It reports the first nonzero derivative of (V_STAM - V_GR) at r=3M using:
  - ordinary r derivatives, and
  - each potential's natural tortoise derivative D_* = d/dr*.

Outputs
-------
  results/G86_redux_action_aware_potential_summary.md
  results/G86_redux_action_aware_potential_derivatives.csv
  plots/G86_redux_action_aware_potential.png

Conventions
-----------
M = 1.  Photon sphere r=3.  ell=2 by default.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Fast truncated-series helpers
# ---------------------------------------------------------------------------

z = sp.symbols("z", real=True)   # r = 3 + z, M=1
u = sp.symbols("u", real=True)   # u = A - 2/3
ell = sp.symbols("ell", positive=True)
ORDER = 12


def trunc(expr, order=ORDER):
    """Truncate a polynomial/power series in z to degree < order."""
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
    """d^n expr / dz^n at z=0 from a truncated polynomial."""
    poly = sp.Poly(sp.expand(expr), z)
    c = poly.coeff_monomial(z**n)
    return sp.simplify(c * math.factorial(n))


def compose_poly_in_u(poly_u, u_z, order=ORDER):
    """Substitute u(z) into a polynomial in u, truncating each power."""
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
    """exp(poly_z) truncated; poly_z starts at positive order."""
    out = sp.Integer(1)
    term = sp.Integer(1)
    # enough because lnf starts at z^4; m=0..3 is plenty for order 12
    for m in range(1, order + 1):
        term = trunc(term * poly_z, order)
        if term == 0:
            break
        out = trunc(out + term / math.factorial(m), order)
    return out


def sqrt_one_plus_series(delta_z, order=ORDER):
    """sqrt(1+delta) truncated. delta starts at z^4 here."""
    out = sp.Integer(1)
    pow_delta = sp.Integer(1)
    # binomial coefficients (1/2 choose m)
    coeff = sp.Integer(1)
    for m in range(1, order + 1):
        pow_delta = trunc(pow_delta * delta_z, order)
        if pow_delta == 0:
            break
        coeff = sp.binomial(sp.Rational(1, 2), m)
        out = trunc(out + coeff * pow_delta, order)
    return out

# ---------------------------------------------------------------------------
# Symbolic local series near r = 3M.
# ---------------------------------------------------------------------------

r = 3 + z
A = sp.Rational(2, 1) / r
y = 3 * A - 2
F = series_expr(1 - 5 * y**4 + 4 * y**5, ORDER)
h = series_expr(1 - A, ORDER)
k_GR_ser = h
k_STAM_ser = trunc(h * F, ORDER)

# V_GR and V_proxy as local series.
M_eff_GR = sp.Integer(1)
M_eff_proxy = trunc(r * (1 - k_STAM_ser) / 2, ORDER)
V_GR_ser = series_expr(h * (ell * (ell + 1) / r**2 - 6 * M_eff_GR / r**3), ORDER)
V_proxy_ser = series_expr(h * (ell * (ell + 1) / r**2 - 6 * M_eff_proxy / r**3), ORDER)

# sqrt(f) local series from d ln f/dA.
A_u = sp.Rational(2, 3) + u
y_u = 3 * A_u - 2
F_u = 1 - 5 * y_u**4 + 4 * y_u**5
q_u = -2 * y_u**3 * (45 * A_u**2 - 33 * A_u - 13) / (A_u * F_u)  # d ln f / dA
q_u_series = sp.series(q_u, u, 0, ORDER).removeO()
lnf_u_series = sp.integrate(q_u_series, (u, 0, u))

u_z = series_expr(A - sp.Rational(2, 3), ORDER)
lnf_z = compose_poly_in_u(lnf_u_series, u_z, ORDER)
sqrtf_z = exp_series(trunc(sp.Rational(1, 2) * lnf_z, ORDER), ORDER)

# Tortoise derivatives.
# GR: dr/dr* = h. STAM: dr/dr* = sqrt(h*k) = h sqrt(F).
sqrtF_ser = sqrt_one_plus_series(trunc(F - 1, ORDER), ORDER)
Dcoef_GR = h
Dcoef_STAM = trunc(h * sqrtF_ser, ORDER)


def Dstar_series(expr, Dcoef, order=ORDER):
    return trunc(Dcoef * sp.diff(expr, z), order)

sqrtf_star_1 = Dstar_series(sqrtf_z, Dcoef_STAM)
sqrtf_star_2 = Dstar_series(sqrtf_star_1, Dcoef_STAM)
action_term_ser = trunc(sqrtf_star_2 * exp_series(trunc(-sp.log(sqrtf_z).series(z,0,ORDER).removeO(), ORDER), ORDER), ORDER)
# The log/exp inverse above can be fragile; replace by direct series division.
action_term_ser = trunc(sp.series(sqrtf_star_2 / sqrtf_z, z, 0, ORDER).removeO(), ORDER)

V_action_ser = trunc(V_proxy_ser + action_term_ser, ORDER)

# ---------------------------------------------------------------------------
# Numerical functions for plotting.
# ---------------------------------------------------------------------------

def F_quintic_np(yv):
    return 1.0 - 5.0 * yv**4 + 4.0 * yv**5


def h_np(rv):
    return 1.0 - 2.0 / rv


def k_stam_np(rv):
    yv = 6.0 / rv - 2.0
    return np.where(yv > 0, h_np(rv) * F_quintic_np(yv), h_np(rv))


def V_proxy_np(rv, ell_val=2, stam=False):
    kval = k_stam_np(rv) if stam else h_np(rv)
    Meff = 0.5 * rv * (1.0 - kval)
    return h_np(rv) * (ell_val * (ell_val + 1) / rv**2 - 6.0 * Meff / rv**3)


def dlnf_dA_np(Aarr):
    Aarr = np.asarray(Aarr, dtype=float)
    yv = 3.0 * Aarr - 2.0
    Fv = F_quintic_np(yv)
    out = np.zeros_like(Aarr)
    mask = Aarr > 2.0 / 3.0
    out[mask] = -2.0 * yv[mask]**3 * (45.0 * Aarr[mask]**2 - 33.0 * Aarr[mask] - 13.0) / (Aarr[mask] * Fv[mask])
    return out


def action_term_numerical(rgrid):
    """Numerical (sqrt(f))''/sqrt(f) in STAM r* for plotting only."""
    Agrid = 2.0 / rgrid
    f = np.ones_like(rgrid)
    mask = Agrid > 2.0 / 3.0
    if np.any(mask):
        Amax = Agrid[mask].max()
        Avals = np.linspace(2.0 / 3.0, Amax, 30000)
        q = dlnf_dA_np(Avals)
        lnf = cumulative_trapezoid(q, Avals, initial=0.0)
        fA = np.exp(np.clip(lnf, -80, 80))
        f[mask] = np.interp(Agrid[mask], Avals, fA)
    sf = np.sqrt(f)
    Dcoef = np.sqrt(np.maximum(h_np(rgrid) * k_stam_np(rgrid), 0.0))
    dsf_dr = np.gradient(sf, rgrid)
    dsf_star = Dcoef * dsf_dr
    d2sf_star = Dcoef * np.gradient(dsf_star, rgrid)
    return d2sf_star / sf

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def first_nonzero(vals):
    for n, v in vals:
        if sp.simplify(v) != 0:
            return n, sp.simplify(v)
    return None, None


def main():
    ell_val = 2
    max_deriv = 8
    print("=" * 90)
    print("G86 redux: action-aware photon-sphere potential structure")
    print("=" * 90)
    print()
    print("Using ell=2, M=1, r=3+z.  Potentials are expanded as local series at z=0.")
    print()

    rows = []
    comparisons = [
        ("proxy_minus_GR", V_proxy_ser, V_GR_ser),
        ("action_minus_GR", V_action_ser, V_GR_ser),
    ]

    for name, V_stam, V_base in comparisons:
        r_vals = []
        star_vals = []
        for n in range(max_deriv + 1):
            rdiff = sp.simplify(coeff_derivative((V_stam - V_base).subs(ell, ell_val), n))
            # Each potential's own tortoise derivative.
            stam_star = V_stam
            base_star = V_base
            for _ in range(n):
                stam_star = Dstar_series(stam_star, Dcoef_STAM)
                base_star = Dstar_series(base_star, Dcoef_GR)
            sdiff = sp.simplify((stam_star - base_star).subs({z: 0, ell: ell_val}))
            r_vals.append((n, rdiff))
            star_vals.append((n, sdiff))
            rows.append({
                "comparison": name,
                "derivative_order": n,
                "r_derivative_diff": str(rdiff),
                "rstar_derivative_diff": str(sdiff),
                "r_derivative_diff_float": float(sp.N(rdiff)),
                "rstar_derivative_diff_float": float(sp.N(sdiff)),
            })
        n_r, v_r = first_nonzero(r_vals)
        n_s, v_s = first_nonzero(star_vals)
        print(f"{name}:")
        print(f"  first nonzero ordinary r derivative:   order {n_r}, value {v_r}")
        print(f"  first nonzero tortoise r* derivative: order {n_s}, value {v_s}")
        print()
        for n in range(max_deriv + 1):
            print(f"    n={n:2d}: d_r^n diff = {r_vals[n][1]} | d_*^n diff = {star_vals[n][1]}")
        print()

    print("Potential values at photon sphere (ell=2):")
    for label, expr in [("V_GR", V_GR_ser), ("V_proxy", V_proxy_ser), ("V_action", V_action_ser)]:
        val = sp.simplify(expr.subs({z: 0, ell: ell_val}))
        print(f"  {label:>9} = {val}  ({float(sp.N(val)):.10f})")
    print()

    csv_path = RESULTS / "G86_redux_action_aware_potential_derivatives.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Plot around PS.
    rgrid = np.linspace(2.05, 4.5, 2500)
    Vgr = V_proxy_np(rgrid, ell_val, stam=False)
    Vpr = V_proxy_np(rgrid, ell_val, stam=True)
    Vact = Vpr + action_term_numerical(rgrid)

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    axes[0].plot(rgrid, Vgr, label="GR RW", lw=2)
    axes[0].plot(rgrid, Vpr, label="STAM proxy", lw=2)
    axes[0].plot(rgrid, Vact, label="STAM action-aware", lw=2)
    axes[0].axvline(3.0, color="k", ls="--", alpha=0.5, label="PS r=3M")
    axes[0].set_ylabel("V(r), ell=2")
    axes[0].set_title("G86 redux: action-aware photon-sphere potential")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(rgrid, Vpr - Vgr, label="proxy - GR", lw=2)
    axes[1].plot(rgrid, Vact - Vgr, label="action - GR", lw=2)
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].axvline(3.0, color="k", ls="--", alpha=0.5)
    axes[1].set_xlabel("r/M")
    axes[1].set_ylabel("ΔV")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    plot_path = PLOTS / "G86_redux_action_aware_potential.png"
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)

    # Summary markdown.
    def first_from_rows(key, comp):
        for rr in rows:
            if rr["comparison"] == comp and rr[key] != "0":
                return rr["derivative_order"], rr[key]
        return None, None

    pr_n_r, pr_v_r = first_from_rows("r_derivative_diff", "proxy_minus_GR")
    ac_n_r, ac_v_r = first_from_rows("r_derivative_diff", "action_minus_GR")
    pr_n_s, pr_v_s = first_from_rows("rstar_derivative_diff", "proxy_minus_GR")
    ac_n_s, ac_v_s = first_from_rows("rstar_derivative_diff", "action_minus_GR")

    summary_path = RESULTS / "G86_redux_action_aware_potential_summary.md"
    md = f"""# G86 redux — action-aware photon-sphere potential structure

This reruns the original G86 structural question after G90-G94.
It compares GR, STAM metric-proxy, and STAM action-aware axial diagnostic potentials near the photon sphere.

## Setup

```text
V_action = V_proxy + (sqrt(f))'' / sqrt(f)
```

where primes are STAM tortoise-coordinate derivatives and `f(A)` is reconstructed locally from the G70/G75 relation with `f(2/3)=1`.

## First nonzero derivative differences at r = 3M, ell = 2

| comparison | ordinary r derivative | tortoise derivative |
|---|---:|---:|
| proxy − GR | order {pr_n_r}: `{pr_v_r}` | order {pr_n_s}: `{pr_v_s}` |
| action − GR | order {ac_n_r}: `{ac_v_r}` | order {ac_n_s}: `{ac_v_s}` |

## Interpretation

The metric-proxy and action-aware potentials both preserve the leading photon-sphere/eikonal match.
The finite-order derivative where the first mismatch appears explains why low-ell, non-eikonal QNM modes can shift even when the leading eikonal photon-orbit result is GR-exact.

The action-aware term tests whether the non-minimal coupling `f(Σ)R` introduces lower-order or larger photon-region structure than the metric proxy alone.
This script is structural only. It does not extract a QNM frequency; G91-G94 handle calibrated time-domain extraction.

## Files

- CSV: `{csv_path}`
- Plot: `{plot_path}`
"""
    summary_path.write_text(md, encoding="utf-8")

    print("Files written:")
    print(f"  {csv_path}")
    print(f"  {summary_path}")
    print(f"  {plot_path}")


if __name__ == "__main__":
    main()
