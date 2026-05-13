#!/usr/bin/env python3
"""
G89_leaver_QNM.py

Frequency-domain QNM extraction via direct shooting on the Regge-Wheeler
wave equation.  Three phases:

G89a: Schwarzschild Leaver-style calibration -- reproduce the known
      ℓ = 2, n = 0 axial RW QNM  ω = 0.373672 - 0.088962 i (M = 1).

G89b: STAM proxy result -- same machinery, replacing the potential with
      V_RW^proxy using the committed k(r) = h(r) outside PS, h(r) F(y)
      inside.

G89c: sensitivity / branch stability -- check robustness w.r.t. r_start
      (proximity to horizon), r_end (matching radius), integration
      tolerance, and initial-guess seeds.

Approach
========
Work in r-coordinate (avoids the cost of inverting r*(r) every step).
Wave equation:
   ψ'' + (h'/h) ψ' + (ω² - V)/h² ψ = 0       (in r-coord, dr*/dr = 1/h)

where prime is d/dr.  (We use the Schwarzschild-style tortoise dr*/dr = 1/h
for the proxy, with V = V_RW^proxy.)

Boundary conditions:
- At horizon r → 2M: ψ ~ (r - 2M)^(-2iMω)  (ingoing)
- At r → ∞:           ψ ~ e^(iωr)         (outgoing)
                      ⟹  dψ/dr ~ (iω/h) ψ  near r_end

QNM condition: shoot ingoing wave from r_start outward; the QNM
frequency makes the resulting ψ at r_end purely outgoing
(i.e. dψ/dr − (iω/h) ψ = 0 at r_end).

Newton/Powell hybrid in complex ω finds the root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Metric and potential
# ---------------------------------------------------------------------------

def h_fn(r, M=1.0):
    return 1.0 - 2.0 * M / r


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def k_Schw(r, M=1.0):
    return h_fn(r, M)


def k_STAM(r, M=1.0):
    A = 2.0 * M / r
    y = 3.0 * A - 2.0
    h_val = h_fn(r, M)
    if y > 0:
        return h_val * F_quintic(y)
    return h_val


def V_RW_proxy(r, ell, k_fn, M=1.0):
    h_val = h_fn(r, M)
    k_val = k_fn(r, M)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


# ---------------------------------------------------------------------------
# Shooting solver
# ---------------------------------------------------------------------------

def shoot(omega_complex, ell, k_fn, M=1.0,
          r_start=2.001, r_end=80.0,
          rtol=1e-10, atol=1e-12):
    """Integrate from r_start near horizon (ingoing BC) to r_end.

    Returns the residual  dψ/dr|_{r_end} - (i ω / h(r_end)) ψ(r_end),
    which vanishes when ω is a QNM.
    """
    delta = r_start - 2.0 * M
    # Leading-order ingoing wave at horizon: ψ ~ (r - 2M)^{-2iMω}
    psi0 = complex(1.0, 0.0)
    psip0 = -2.0j * M * omega_complex / delta

    def rhs(r, y):
        psi, psip = y[0], y[1]
        h_val = h_fn(r, M)
        hp_val = 2.0 * M / r**2  # h'(r)
        V = V_RW_proxy(r, ell, k_fn, M)
        psipp = -(hp_val / h_val) * psip + (V - omega_complex**2) / h_val**2 * psi
        return [psip, psipp]

    sol = solve_ivp(rhs, [r_start, r_end], [psi0, psip0],
                    method='RK45', rtol=rtol, atol=atol, max_step=0.5)
    if not sol.success:
        return complex(np.nan, np.nan)

    psi_end = sol.y[0, -1]
    psip_end = sol.y[1, -1]
    h_end = h_fn(r_end, M)
    return psip_end - 1j * omega_complex / h_end * psi_end


def find_QNM(omega_init, ell, k_fn, **kwargs):
    """Find ω as root of shooting residual via Powell's hybrid method."""
    def residual_real(x):
        omega = complex(x[0], x[1])
        r = shoot(omega, ell, k_fn, **kwargs)
        return [r.real, r.imag]

    sol = root(residual_real, [omega_init.real, omega_init.imag],
               method='hybr', tol=1e-9, options={'xtol': 1e-10})
    omega_found = complex(sol.x[0], sol.x[1])
    return omega_found, sol


# ---------------------------------------------------------------------------
# G89a: Schwarzschild calibration
# ---------------------------------------------------------------------------

print("=" * 80)
print("G89a: Schwarzschild Leaver calibration (ℓ = 2, n = 0)")
print("=" * 80, flush=True)
print()

omega_gold = complex(0.373672, -0.088962)
print(f"Leaver gold standard:  ω = {omega_gold}", flush=True)
print()

# Initial guess: slightly off from gold value to test convergence
omega_init = complex(0.37, -0.09)
print(f"Initial guess:         ω₀ = {omega_init}", flush=True)
print()

print("Shooting iteration ...", flush=True)
omega_Schw, sol_Schw = find_QNM(omega_init, ell=2, k_fn=k_Schw,
                                 r_start=2.001, r_end=80.0)
print(f"Converged:             ω_Schw = {omega_Schw}", flush=True)
print()

err_Schw = abs(omega_Schw - omega_gold) / abs(omega_gold)
print(f"Calibration error: |Δω/ω| = {err_Schw * 100:.4f} %")
if err_Schw < 0.01:
    print("PASS: Schwarzschild calibration within 1% — method is reliable.", flush=True)
elif err_Schw < 0.05:
    print("ACCEPTABLE: 1-5% calibration. Sub-percent accuracy needs more care.", flush=True)
else:
    print("FAIL: > 5% calibration error. Need to refine numerical settings.", flush=True)
print()


# ---------------------------------------------------------------------------
# G89b: STAM proxy result
# ---------------------------------------------------------------------------

print("=" * 80)
print("G89b: STAM proxy result (ℓ = 2, n = 0)")
print("=" * 80, flush=True)
print()
print("Same machinery, V_RW^proxy with committed k(r) = h(r) F(y) inside PS.", flush=True)
print()

# Initial guess: use converged Schwarzschild value as warm start
print(f"Initial guess (warm-start from Schwarzschild):  ω₀ = {omega_Schw}", flush=True)
print()
print("Shooting iteration ...", flush=True)
omega_STAM, sol_STAM = find_QNM(omega_Schw, ell=2, k_fn=k_STAM,
                                 r_start=2.001, r_end=80.0)
print(f"Converged:             ω_STAM = {omega_STAM}", flush=True)
print()

if not np.isnan(omega_STAM.real):
    dw_over_w = abs(omega_STAM - omega_Schw) / abs(omega_Schw)
    print(f"Framework shift:  |Δω/ω| = {dw_over_w * 100:.4f} %")
    print(f"  Re shift:   {(omega_STAM.real - omega_Schw.real)/omega_Schw.real*100:+.4f} %")
    print(f"  Im shift:   {(omega_STAM.imag - omega_Schw.imag)/omega_Schw.imag*100:+.4f} %")
    print()
print()


# ---------------------------------------------------------------------------
# G89c: Sensitivity / branch-stability check
# ---------------------------------------------------------------------------

print("=" * 80)
print("G89c: Sensitivity check")
print("=" * 80, flush=True)
print()

# Vary r_start (proximity to horizon)
print("(c.1) r_start variation (Schwarzschild ℓ=2, n=0):", flush=True)
print(f"  {'r_start':>10}  {'ω':>30}  {'err vs gold':>14}", flush=True)
for r_s in [2.0001, 2.001, 2.01, 2.05, 2.1]:
    try:
        omg, _ = find_QNM(omega_init, ell=2, k_fn=k_Schw,
                          r_start=r_s, r_end=80.0)
        err = abs(omg - omega_gold) / abs(omega_gold) * 100
        print(f"  {r_s:>10.4f}  {omg.real:+.5f}{omg.imag:+.5f}j  {err:>13.4f} %", flush=True)
    except Exception as e:
        print(f"  {r_s:>10.4f}  failed: {e}", flush=True)
print()

# Vary r_end (matching radius)
print("(c.2) r_end variation (Schwarzschild ℓ=2, n=0):", flush=True)
print(f"  {'r_end':>10}  {'ω':>30}  {'err vs gold':>14}", flush=True)
for r_e in [30, 50, 80, 120, 200]:
    try:
        omg, _ = find_QNM(omega_init, ell=2, k_fn=k_Schw,
                          r_start=2.001, r_end=r_e)
        err = abs(omg - omega_gold) / abs(omega_gold) * 100
        print(f"  {r_e:>10}  {omg.real:+.5f}{omg.imag:+.5f}j  {err:>13.4f} %", flush=True)
    except Exception as e:
        print(f"  {r_e:>10}  failed: {e}", flush=True)
print()

# Vary initial guess (Schwarzschild branch stability)
print("(c.3) Initial guess variation (Schwarzschild ℓ=2, n=0):", flush=True)
print(f"  {'ω_init':>22}  {'ω converged':>30}  {'err':>10}", flush=True)
for omega_seed in [complex(0.30, -0.05), complex(0.40, -0.10),
                    complex(0.37, -0.05), complex(0.45, -0.15),
                    complex(0.60, -0.20)]:
    try:
        omg, _ = find_QNM(omega_seed, ell=2, k_fn=k_Schw,
                          r_start=2.001, r_end=80.0)
        err = abs(omg - omega_gold) / abs(omega_gold) * 100
        print(f"  {omega_seed.real:+.3f}{omega_seed.imag:+.3f}j  "
              f"{omg.real:+.5f}{omg.imag:+.5f}j  {err:>9.4f} %", flush=True)
    except Exception as e:
        print(f"  seed {omega_seed} failed: {e}", flush=True)
print()

# Vary integration tolerance
print("(c.4) Integration tolerance variation (Schwarzschild):", flush=True)
print(f"  {'rtol':>10}  {'ω':>30}  {'err':>10}", flush=True)
for tol in [1e-6, 1e-8, 1e-10, 1e-12]:
    try:
        omg, _ = find_QNM(omega_init, ell=2, k_fn=k_Schw,
                          r_start=2.001, r_end=80.0,
                          rtol=tol, atol=tol * 1e-2)
        err = abs(omg - omega_gold) / abs(omega_gold) * 100
        print(f"  {tol:>10.0e}  {omg.real:+.5f}{omg.imag:+.5f}j  {err:>9.4f} %", flush=True)
    except Exception as e:
        print(f"  tol {tol} failed: {e}", flush=True)
print()


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

print("=" * 80)
print("VERDICT")
print("=" * 80, flush=True)
print()

if err_Schw < 0.01:
    print("CALIBRATION:  passed (< 1% vs Leaver gold).")
else:
    print(f"CALIBRATION:  not yet at sub-percent precision ({err_Schw*100:.2f} % vs gold).")
    print("              Numerical refinement needed for locked framework prediction.")
print()

if not np.isnan(omega_STAM.real):
    print(f"FRAMEWORK PREDICTION (with current calibration uncertainty):")
    print(f"  ω_Schw   = {omega_Schw}")
    print(f"  ω_STAM   = {omega_STAM}")
    dw = abs(omega_STAM - omega_Schw) / abs(omega_Schw)
    print(f"  |Δω/ω|   = {dw * 100:.4f} %")
    print()
    if dw > 0.10:
        print("  > 10% → near-term LIGO-facing prediction territory.")
    elif dw > 0.01:
        print("  1-10% → percent-level deviation; within LISA/precision-LIGO reach.")
    else:
        print("  < 1%  → framework matches Schwarzschild ringdown closely.")
print()


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

md = []
md.append("# G89 — Frequency-domain QNM via shooting (ℓ = 2, n = 0)\n")
md.append("**Date: 2026-05-13.**  Direct shooting solver for the Regge-Wheeler "
          "equation with Schwarzschild-style tortoise (proxy). Calibration "
          "against Leaver gold standard.\n")
md.append("## G89a — Schwarzschild calibration\n")
md.append(f"- Leaver gold:  ω = {omega_gold.real:.6f} − {abs(omega_gold.imag):.6f} i\n")
md.append(f"- Shooting:     ω = {omega_Schw.real:.6f} {'-' if omega_Schw.imag < 0 else '+'} "
          f"{abs(omega_Schw.imag):.6f} i\n")
md.append(f"- Calibration error: **{err_Schw * 100:.4f} %**\n")
md.append("\n")
md.append("## G89b — STAM proxy result\n")
if not np.isnan(omega_STAM.real):
    md.append(f"- ω_STAM:  {omega_STAM.real:.6f} {'-' if omega_STAM.imag < 0 else '+'} "
              f"{abs(omega_STAM.imag):.6f} i\n")
    dw_pct = abs(omega_STAM - omega_Schw) / abs(omega_Schw) * 100
    md.append(f"- Framework shift: **|Δω/ω| = {dw_pct:.4f} %**\n")
md.append("\n")
md.append("## G89c — Sensitivity check\n")
md.append("See script output for r_start, r_end, initial-guess, and tolerance "
          "sensitivity sweeps.\n")
md.append("\n")
md.append("## Caveats\n")
md.append("- Proxy potential: this uses V_RW^proxy = h[ℓ(ℓ+1)/r² − 6 M_eff/r³] "
          "with M_eff = (r/2)(1−k). The rigorous axial perturbation potential "
          "for the constrained-Σ action requires linearizing f(Σ)R + Lagrange-"
          "multiplier terms and is a separate computation (G90+).\n")
md.append("- Schwarzschild-style tortoise: integration uses dr*/dr = 1/h. The "
          "framework's actual tortoise has cubic-vanishing behavior inside PS "
          "due to k → 0³ at horizon; the proxy keeps the Schwarzschild log "
          "form for tractability.\n")
md.append("\n")
md.append("## Files\n")
md.append("- [scripts/G89_leaver_QNM.py](../scripts/G89_leaver_QNM.py)\n")

out_md = RESULTS / "G89_leaver_QNM_summary.md"
out_md.write_text("\n".join(md), encoding="utf-8")
print(f"Summary written: {out_md}", flush=True)
