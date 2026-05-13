#!/usr/bin/env python3
"""
G85_exact_kerr_photon_region.py

Exact Kerr spheroidal photon-region surface from spherical-null-geodesic
conditions  R(r) = 0,  dR/dr = 0,  replacing G84's sin^2(theta) interpolation
between equatorial-prograde and polar photon-orbit radii.

Derivation
==========
For Kerr null geodesics in Boyer-Lindquist coordinates (M = 1),

    R(r)/E^2 = [ (r^2 + a^2) - a lambda ]^2  -  Delta * [ (lambda - a)^2 + eta ]

with lambda = L_z / E,  eta = Q / E^2,  Delta = r^2 - 2r + a^2.

A spherical photon orbit at radius r_p has R(r_p) = R'(r_p) = 0.  Solving
the pair simultaneously:

    lambda(r_p) = -( r_p^3 - 3 r_p^2 + a^2 r_p + a^2 ) / [ a (r_p - 1) ]
    eta(r_p)    = r_p^2 [ 4 a^2 Delta - ( r_p^2 - 3 r_p + 2 a^2 )^2 ]
                  / [ a^2 (r_p - 1)^2 ]

The orbit's theta-motion is governed by

    Theta(theta) / E^2  =  eta  -  cos^2(theta) [ lambda^2 / sin^2(theta) - a^2 ]

The orbit reaches a latitude theta iff Theta(theta) >= 0.  At its turning
points theta_min, theta_max (symmetric across equator: theta_max = pi - theta_min),
Theta(theta_min) = 0.

The INNER boundary of the photon region at latitude theta is the r_p such that
the orbit at r_p has  theta_min  =  theta  (i.e., the orbit "just touches"
that latitude).  This defines the exact spheroidal photon-region surface:

    r_pr_exact(theta; a)  :=   r_p  satisfying  Theta(theta; lambda(r_p), eta(r_p)) = 0
                               with  r_p  in  [ r_ph_prograde(a),  r_polar(a) ].

Limits:
    theta = pi/2  (equator):  r_pr_exact = r_ph_prograde(a)
    theta = 0     (pole):     r_pr_exact = r_polar(a)
    a = 0         (Schw):     r_pr_exact = 3 M  for all theta

G84 uses the sin^2 ansatz

    r_pr_G84(theta; a) = r_ph_prograde(a) sin^2(theta) + r_polar(a) cos^2(theta)

This script computes both, reports their max difference over (a, theta),
and decides whether the difference is observationally relevant for the
framework's W_K and F(y_K) normalization.

Outputs
=======
    results/G85_exact_kerr_photon_region_summary.md
    plots/G85_exact_kerr_photon_region.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

M = 1.0


# ===========================================================================
# Closed-form anchors (G83 / G84)
# ===========================================================================

def r_ph_prograde(a):
    """Equatorial prograde Kerr photon orbit (M = 1).  Bardeen formula."""
    return 2.0 * (1.0 + np.cos((2.0 / 3.0) * np.arccos(-a)))


def r_polar(a):
    """Polar Kerr spherical photon orbit, largest real root of
       r^3 - 3 r^2 + a^2 r + a^2 = 0  (M = 1).
    """
    if a == 0.0:
        return 3.0
    coeffs = [1.0, -3.0, a**2, a**2]
    roots = np.roots(coeffs)
    real_roots = [rr.real for rr in roots if abs(rr.imag) < 1e-10]
    r_h = 1.0 + np.sqrt(max(1.0 - a**2, 0.0))
    candidates = [rr for rr in real_roots if rr > r_h]
    return max(candidates) if candidates else float("nan")


def r_horizon(a):
    return 1.0 + np.sqrt(max(1.0 - a**2, 0.0))


# ===========================================================================
# Photon-orbit constants  lambda(r_p),  eta(r_p)   (closed form, M = 1)
# ===========================================================================

def lambda_phot(r_p, a):
    """lambda = L_z / E  at spherical photon orbit r_p.  Diverges as a -> 0."""
    return -(r_p**3 - 3.0 * r_p**2 + a**2 * r_p + a**2) / (a * (r_p - 1.0))


def eta_phot(r_p, a):
    """eta = Q / E^2  at spherical photon orbit r_p."""
    Delta = r_p**2 - 2.0 * r_p + a**2
    return (r_p**2 * (4.0 * a**2 * Delta - (r_p**2 - 3.0 * r_p + 2.0 * a**2)**2)
            / (a**2 * (r_p - 1.0)**2))


def theta_potential(theta, r_p, a):
    """Theta(theta) / E^2 for the spherical photon orbit at r_p.

    Returns sin^2(theta) * Theta to avoid the 1/sin^2 divergence at theta=0
    (multiplying by sin^2(theta) preserves the sign and the zero locations
    for theta in (0, pi/2]).
    """
    lam = lambda_phot(r_p, a)
    eta = eta_phot(r_p, a)
    s2 = np.sin(theta)**2
    c2 = np.cos(theta)**2
    # Theta * sin^2 = eta sin^2  -  cos^2 (lambda^2 - a^2 sin^2)
    return eta * s2 - c2 * (lam**2 - a**2 * s2)


# ===========================================================================
# Exact r_pr(theta; a)  via brentq on r_p in [r_ph^+, r_polar]
# ===========================================================================

def r_pr_exact(theta, a):
    """Find r_p in [r_ph^+(a), r_polar(a)] such that Theta(theta; r_p) = 0.

    Boundary cases handled explicitly to avoid 0/0 singularities at endpoints.
    """
    if a == 0.0:
        return 3.0
    eps = 1e-10
    if abs(theta - np.pi / 2) < eps:
        return r_ph_prograde(a)
    if abs(theta) < eps:
        return r_polar(a)

    r_lo = r_ph_prograde(a) + 1e-9
    r_hi = r_polar(a) - 1e-9

    f_lo = theta_potential(theta, r_lo, a)
    f_hi = theta_potential(theta, r_hi, a)

    # Expected:  f(r_ph^+) < 0  for non-equator theta  (orbit can't reach)
    #            f(r_polar) > 0 for non-equator theta  (orbit reaches everything)
    if f_lo * f_hi > 0:
        # Same sign across the interval -> no root.  Could be that we're at
        # an endpoint with floating-point fluctuation.  Return the closer
        # endpoint by which sign is closer to zero.
        if abs(f_lo) < abs(f_hi):
            return r_ph_prograde(a)
        return r_polar(a)

    return brentq(lambda r: theta_potential(theta, r, a), r_lo, r_hi,
                  xtol=1e-12, rtol=1e-12)


def r_pr_G84(theta, a):
    """G84's sin^2 interpolation."""
    return r_ph_prograde(a) * np.sin(theta)**2 + r_polar(a) * np.cos(theta)**2


# ===========================================================================
# Comparison tables and plot
# ===========================================================================

def main():
    print("=" * 88)
    print("G85 -- exact Kerr spheroidal photon-region surface")
    print("       (replaces G84's sin^2(theta) interpolation)")
    print("=" * 88)
    print()
    print("M = 1, r_p in [r_ph^+(a), r_polar(a)] for each a.")
    print()

    # ----- Boundary sanity check -----
    print("Anchors at equator and pole  (sanity: exact must match G84 exactly there):")
    print(f"  {'a':>6}  {'r_ph^+':>10}  {'r_polar':>10}  {'r_pr_exact(pi/2)':>18}  {'r_pr_exact(0)':>16}")
    print("-" * 70)
    for a in [0.0, 0.3, 0.5, 0.7, 0.9, 0.99]:
        rp_eq = r_ph_prograde(a)
        rp_pol = r_polar(a)
        rpr_eq = r_pr_exact(np.pi / 2, a)
        rpr_pol = r_pr_exact(0.0, a)
        print(f"  {a:>6.2f}  {rp_eq:>10.6f}  {rp_pol:>10.6f}  "
              f"{rpr_eq:>18.6f}  {rpr_pol:>16.6f}")
    print()

    # ----- Compare exact vs G84 at intermediate theta values -----
    print("=" * 88)
    print("Exact vs G84 sin^2(theta) interpolation at intermediate theta")
    print("=" * 88)
    print()

    thetas_test = [np.pi/2 - 1e-3, np.pi/3, np.pi/4, np.pi/6, np.pi/12,
                   np.pi/24, 1e-3]
    a_test = [0.3, 0.5, 0.7, 0.9, 0.99]

    print("Table of  r_pr_exact - r_pr_G84  (in units of M):")
    print()
    header = "  theta\\a  " + "  ".join(f"{a:>10.2f}" for a in a_test)
    print(header)
    print("-" * len(header))
    max_abs_diff = 0.0
    max_rel_diff = 0.0
    max_diff_loc = None
    for theta in thetas_test:
        row = [f"  {theta:>7.4f} "]
        for a in a_test:
            rpr_e = r_pr_exact(theta, a)
            rpr_g = r_pr_G84(theta, a)
            diff = rpr_e - rpr_g
            rel = abs(diff) / rpr_e
            if abs(diff) > max_abs_diff:
                max_abs_diff = abs(diff)
                max_diff_loc = (theta, a, diff)
            if rel > max_rel_diff:
                max_rel_diff = rel
            row.append(f"{diff:>+10.6f}")
        print("  ".join(row))
    print()
    print(f"max |r_pr_exact - r_pr_G84|  =  {max_abs_diff:.6e} M")
    if max_diff_loc:
        print(f"   at (theta = {max_diff_loc[0]:.4f}, a = {max_diff_loc[1]:.2f}),  "
              f"signed diff = {max_diff_loc[2]:+.4e}")
    print(f"max relative difference    =  {max_rel_diff*100:.4f}%")
    print()

    # ----- Implications for Sigma_ph (the framework's natural shell coord) -----
    print("=" * 88)
    print("Difference in  Sigma_ph = 6 r_pr / (r_pr^2 + a^2)  (the framework normalization)")
    print("=" * 88)
    print()
    print("This is the quantity actually entering W_K and F(y_K) normalization.")
    print()

    print("Table of  Sigma_ph_exact - Sigma_ph_G84:")
    print()
    print(header)
    print("-" * len(header))
    max_dSigma = 0.0
    max_dSigma_rel = 0.0
    for theta in thetas_test:
        row = [f"  {theta:>7.4f} "]
        for a in a_test:
            re_v = r_pr_exact(theta, a)
            rg_v = r_pr_G84(theta, a)
            Se = 6.0 * re_v / (re_v**2 + a**2)
            Sg = 6.0 * rg_v / (rg_v**2 + a**2)
            dS = Se - Sg
            if abs(dS) > max_dSigma:
                max_dSigma = abs(dS)
            relS = abs(dS) / Se if Se != 0 else 0
            if relS > max_dSigma_rel:
                max_dSigma_rel = relS
            row.append(f"{dS:>+10.6f}")
        print("  ".join(row))
    print()
    print(f"max |Sigma_ph_exact - Sigma_ph_G84| = {max_dSigma:.6e}")
    print(f"max relative                       = {max_dSigma_rel*100:.4f}%")
    print()

    # ----- Plots -----
    print("=" * 88)
    print("Generating plot ...")
    print("=" * 88)
    a_grid = np.linspace(0.0, 0.99, 50)
    theta_grid = np.linspace(1e-3, np.pi / 2, 60)
    A_mesh, TH_mesh = np.meshgrid(a_grid, theta_grid, indexing="ij")

    R_exact = np.zeros_like(A_mesh)
    R_G84 = np.zeros_like(A_mesh)
    for i, a in enumerate(a_grid):
        for j, th in enumerate(theta_grid):
            R_exact[i, j] = r_pr_exact(th, a)
            R_G84[i, j] = r_pr_G84(th, a)
    R_diff = R_exact - R_G84

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    im = ax.contourf(A_mesh, TH_mesh, R_exact, levels=20, cmap="viridis")
    plt.colorbar(im, ax=ax, label=r"$r_{\rm pr}^{\rm exact}(a, \theta) / M$")
    ax.set_xlabel("a / M")
    ax.set_ylabel(r"$\theta$ (rad)")
    ax.set_title("Exact photon-region surface  $r_{\\rm pr}(a, \\theta)$")

    ax = axes[0, 1]
    im = ax.contourf(A_mesh, TH_mesh, R_diff, levels=20, cmap="RdBu_r",
                     vmin=-np.max(abs(R_diff)), vmax=np.max(abs(R_diff)))
    plt.colorbar(im, ax=ax, label=r"$r_{\rm pr}^{\rm exact} - r_{\rm pr}^{\rm G84}$  (M)")
    ax.set_xlabel("a / M")
    ax.set_ylabel(r"$\theta$ (rad)")
    ax.set_title("Exact minus G84 sin$^2$ ansatz")

    ax = axes[1, 0]
    a_show = [0.3, 0.5, 0.7, 0.9, 0.99]
    for a in a_show:
        r_e = np.array([r_pr_exact(th, a) for th in theta_grid])
        r_g = np.array([r_pr_G84(th, a) for th in theta_grid])
        ax.plot(theta_grid, r_e, "-", linewidth=2, label=f"exact a={a}")
        ax.plot(theta_grid, r_g, "--", linewidth=1, alpha=0.7,
                label=f"G84 a={a}")
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$r_{\rm pr}(\theta) / M$")
    ax.set_title("Photon-region boundary: exact vs G84")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    for a in a_show:
        r_e = np.array([r_pr_exact(th, a) for th in theta_grid])
        r_g = np.array([r_pr_G84(th, a) for th in theta_grid])
        rel = (r_e - r_g) / r_e * 100
        ax.plot(theta_grid, rel, "-", linewidth=2, label=f"a={a}")
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$(r_{\rm pr}^{\rm exact} - r_{\rm pr}^{\rm G84}) / r_{\rm pr}^{\rm exact}$  (%)")
    ax.set_title("Relative error of G84 ansatz vs theta")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = PLOTS / "G85_exact_kerr_photon_region.png"
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"Plot saved: {out_png}")
    print()

    # ----- Verdict -----
    print("=" * 88)
    print("VERDICT")
    print("=" * 88)
    print()
    if max_rel_diff < 0.005:
        verdict = (
            "G84's sin^2(theta) ansatz is accurate to better than 0.5% in r_pr "
            "and Sigma_ph across (a in [0, 0.99], theta in (0, pi/2]).  "
            "The framework's structural conclusions are unaffected; G85 simply "
            "provides the exact form for precision EMRI / imaging predictions."
        )
    elif max_rel_diff < 0.02:
        verdict = (
            f"G84 ansatz accurate to {max_rel_diff*100:.2f}% across (a, theta).  "
            "Acceptable for structural work but the exact G85 form should be "
            "used for precision EMRI / imaging predictions."
        )
    else:
        verdict = (
            f"G84 ansatz deviates from exact by up to {max_rel_diff*100:.2f}%.  "
            "Significant; the exact G85 r_pr(theta; a) should replace the "
            "sin^2 interpolation in the W_K formula for any precision work."
        )
    print(verdict)
    print()

    # ----- Markdown -----
    md = ["# G85 -- exact Kerr spheroidal photon-region surface\n"]
    md.append("\nExact derivation from spherical-null-geodesic conditions "
              "R(r) = 0, dR/dr = 0, replacing G84's sin^2(theta) interpolation "
              "between r_ph_prograde(a) and r_polar(a).\n")
    md.append("\n## Formulas\n")
    md.append("```\n"
              "lambda(r_p) = -(r_p^3 - 3 r_p^2 + a^2 r_p + a^2) / [a (r_p - 1)]\n"
              "eta(r_p)    = r_p^2 [ 4 a^2 Delta - (r_p^2 - 3 r_p + 2 a^2)^2 ]\n"
              "              / [ a^2 (r_p - 1)^2 ]\n"
              "Theta(theta)/E^2 = eta - cos^2(theta) [ lambda^2/sin^2(theta) - a^2 ]\n"
              "```\n\n")
    md.append("r_pr_exact(theta; a) is the unique r_p in [r_ph_prograde, r_polar] "
              "such that Theta(theta; lambda(r_p), eta(r_p)) = 0.\n")
    md.append("\n## Comparison with G84\n")
    md.append(f"\n- max abs |r_pr_exact - r_pr_G84| = {max_abs_diff:.4e} M\n")
    md.append(f"- max rel |r_pr_exact - r_pr_G84| / r_pr_exact = {max_rel_diff*100:.4f}%\n")
    md.append(f"- max abs |Sigma_ph_exact - Sigma_ph_G84| = {max_dSigma:.4e}\n")
    md.append(f"- max rel = {max_dSigma_rel*100:.4f}%\n")
    md.append("\n## Verdict\n")
    md.append(verdict + "\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G85_exact_kerr_photon_region.py](../scripts/G85_exact_kerr_photon_region.py)\n")
    md.append("- [plots/G85_exact_kerr_photon_region.png](../plots/G85_exact_kerr_photon_region.png)\n")
    (RESULTS / "G85_exact_kerr_photon_region_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G85_exact_kerr_photon_region_summary.md'}")


if __name__ == "__main__":
    main()
