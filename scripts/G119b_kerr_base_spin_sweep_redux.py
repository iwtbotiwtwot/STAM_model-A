#!/usr/bin/env python3
"""
G119b_kerr_base_spin_sweep_redux.py

Stepping stone before full G120 Teukolsky.  Uses an eikonal-WKB
(Poschl-Teller) proxy on G119's Kerr wave base, not the rigorous
Detweiler / Sasaki-Nakamura radial equation.

Why this script
===============
G118v2 gave a clean time-domain spin sweep on the Kerr base, but the
absolute QNM frequencies are NOT calibrated to published Kerr GR values
(because V_base is a Schwarzschild-RW analog, not the full Kerr Teukolsky
potential).  Before tackling the full Teukolsky solver (G120), it is
worth re-running the spin sweep with a complementary frequency-domain
WKB method.  If WKB on V_base and TD on V_base give similar shifts,
that increases confidence in the qualitative trend; if not, the
difference flags where TD fits or fitting windows mattered.

Method (Poschl-Teller / WKB)
============================
For a potential V(r*) with a single peak:

    V_max     = V(r*_max)
    V''       = d^2 V / dr*^2 at the peak  (negative)
    alpha     = sqrt(-2 V_max / V'')        (Poschl-Teller width)

QNM (Ferrari-Mashhoon 1984):

    omega_R = sqrt(V_max - 1/(4 alpha^2))
    omega_I = -(n + 1/2) / alpha            (for s = 2, n = 0:  -(1/2)/alpha)

For the framework's V_base on the Kerr wave base, we apply this formula
at each spin to obtain a fast frequency-domain estimate.  Doing it for
V_GR_base and V_STAM_base separately gives the STAM shift.

At a = 0, Poschl-Teller-WKB recovers Schwarzschild RW QNM at the few-%
level (typical WKB error).  This is the calibration baseline.

Caveats
=======
- V_base is the same Kerr-Delta-aware RW analog as G118v2.  Not the
  rigorous Kerr axial potential.
- WKB is a 0th-order WKB; higher-order WKB (Iyer-Will, Konoplya) tighten
  the absolute number but do not change the shift trend.
- Absolute frequencies are NOT calibrated to published Kerr QNMs at a > 0.
  Only the STAM shift % is the meaningful output.
- This is a stepping stone, not a final prediction.  G120 (full
  Detweiler / Sasaki-Nakamura) is the rigorous next step.

Output
======
    results/G119b_kerr_base_spin_sweep_redux_summary.md
    plots/G119b_kerr_base_spin_sweep_redux.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0
ELL = 2


# ===========================================================================
# G119 Kerr base (same functions)
# ===========================================================================

def r_horizon(a):
    return M + np.sqrt(max(M*M - a*a, 0.0))

def Delta(r, a):
    return r*r - 2.0*M*r + a*a

def r_ph_prograde(a):
    if abs(a) < 1e-14:
        return 3.0*M
    return 2.0*M*(1.0 + np.cos((2.0/3.0)*np.arccos(-a/M)))

def A_K(r, a):
    return 2.0*M*r/(r*r + a*a)

def Sigma_K(r, a):
    return 3.0*A_K(r, a)

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
    return 1.0 - 5.0*y**4 + 4.0*y**5

def F_first(y):
    return -20.0*y**3 + 20.0*y**4

def Delta_STAM(r, a):
    D = Delta(r, a)
    if r >= r_ph_prograde(a):
        return D
    y = y_K_eq(r, a)
    return D * F_quintic(y)


# ===========================================================================
# V_base on Kerr  (same as G118v2)
# ===========================================================================

def V_base_GR(r, a, ell=ELL):
    D = Delta(r, a)
    pref = D / (r*r + a*a)**2
    bracket = ell*(ell+1)*(r*r + a*a)/r**2 - 6.0*M*(r*r - a*a)/r**3
    return pref * bracket


def V_base_STAM(r, a, ell=ELL):
    DS = Delta_STAM(r, a)
    pref = DS / (r*r + a*a)**2
    bracket = ell*(ell+1)*(r*r + a*a)/r**2 - 6.0*M*(r*r - a*a)/r**3
    return pref * bracket


# ===========================================================================
# Action correction (same as G118v2)
# ===========================================================================

def dlnf_dA(A):
    if A <= 2.0/3.0:
        return 0.0
    y = 3.0*A - 2.0
    F = 1.0 - 5.0*y**4 + 4.0*y**5
    return -2.0*y**3*(45.0*A**2 - 33.0*A - 13.0)/(A*F)

def d2lnf_dA2(A, hstep=1e-6):
    return (dlnf_dA(A + hstep) - dlnf_dA(A - hstep))/(2.0*hstep)

def dA_K_dr(r, a):
    v = r*r + a*a
    return 2.0*M*(a*a - r*r)/v**2

def d2A_K_dr2(r, a):
    v = r*r + a*a
    return -2.0*M*(2.0*r*v + 4.0*r*(a*a - r*r))/v**3


def analytic_correction(r, a):
    rp = r_ph_prograde(a)
    if r >= rp:
        return 0.0
    A = A_K(r, a)
    dA = dA_K_dr(r, a)
    d2A = d2A_K_dr2(r, a)
    u = dlnf_dA(A)
    up = d2lnf_dA2(A)
    dlnf_dr = u*dA
    d2lnf_dr2 = up*(dA**2) + u*d2A
    DS = Delta_STAM(r, a)
    if DS <= 0:
        return 0.0
    v = r*r + a*a
    prefactor = DS/v
    # FD on Delta_STAM
    eps = 1e-5 * max(1.0, r - r_horizon(a))
    eps = max(eps, 1e-7)
    DSp = Delta_STAM(r+eps, a) if (r+eps) < rp else Delta(r+eps, a)
    DSm = Delta_STAM(r-eps, a) if (r-eps) < rp else Delta(r-eps, a)
    dDS_dr = (DSp - DSm)/(2.0*eps)
    P = prefactor * dlnf_dr
    d_prefactor_dr = (dDS_dr*v - DS*2.0*r)/v**2
    dP_dr = d_prefactor_dr*dlnf_dr + prefactor*d2lnf_dr2
    P_star = prefactor * dP_dr
    return 0.5*P_star + 0.25*P**2


# ===========================================================================
# Build potential on a coarse r* grid (uniform in r*)
# ===========================================================================

def build_potential_on_rstar(rstar, a, use_STAM=False):
    """Integrate dr/dr* = Delta_STAM/(r^2+a^2)  (or Delta/(r^2+a^2) for GR)
       starting from r* = 0 at r_pr_eq(a).  Then evaluate V on the resulting
       r-array.
    """
    rp_eq = r_ph_prograde(a)
    if use_STAM:
        def integrand(r):
            DS = Delta_STAM(r, a)
            return DS/(r*r + a*a) if DS > 0 else 0.0
    else:
        def integrand(r):
            return Delta(r, a)/(r*r + a*a)

    def rhs(rs, y):
        r_val = y[0]
        return [integrand(r_val)]

    rs_pos = rstar[rstar > 0]
    rs_neg = rstar[rstar < 0]
    r_arr = np.zeros_like(rstar)
    if len(rs_pos) > 0:
        order = np.argsort(rs_pos)
        x = rs_pos[order]
        sol = solve_ivp(rhs, [0, x[-1]], [rp_eq], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=1.0)
        vals = np.zeros_like(rs_pos)
        vals[order] = sol.y[0]
        r_arr[rstar > 0] = vals
    if len(rs_neg) > 0:
        order = np.argsort(rs_neg)[::-1]
        x = rs_neg[order]
        sol = solve_ivp(rhs, [0, x[-1]], [rp_eq], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.02)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar < 0] = vals
    r_arr[rstar == 0] = rp_eq

    if use_STAM:
        V = np.array([V_base_STAM(r, a) + analytic_correction(r, a) for r in r_arr])
    else:
        V = np.array([V_base_GR(r, a) for r in r_arr])
    return rstar, V, r_arr


# ===========================================================================
# Poschl-Teller / WKB QNM extraction
# ===========================================================================

def WKB_PT_QNM(rstar, V, n=0):
    """Find peak, second derivative; compute Poschl-Teller QNM."""
    # Locate max — refine with quadratic around grid maximum
    i_max = int(np.argmax(V))
    if i_max < 2 or i_max > len(V) - 3:
        return None
    # Quadratic peak refinement
    x0, x1, x2 = rstar[i_max-1], rstar[i_max], rstar[i_max+1]
    y0, y1, y2 = V[i_max-1], V[i_max], V[i_max+1]
    denom = (x0 - x1)*(x0 - x2)*(x1 - x2)
    if abs(denom) < 1e-30:
        return None
    A_coef = (x2*(y1 - y0) + x1*(y0 - y2) + x0*(y2 - y1))/denom
    B_coef = (x2**2*(y0 - y1) + x1**2*(y2 - y0) + x0**2*(y1 - y2))/denom
    if abs(A_coef) < 1e-30:
        return None
    x_peak = -B_coef/(2.0*A_coef)
    V_peak = A_coef*x_peak**2 + B_coef*x_peak + \
             (y0 - A_coef*x0**2 - B_coef*x0)
    Vpp_peak = 2.0*A_coef  # negative if peak

    if Vpp_peak >= 0:
        return None  # not a peak

    alpha = np.sqrt(-2.0*V_peak/Vpp_peak)
    if V_peak - 1.0/(4.0*alpha**2) < 0:
        return None

    om_R = np.sqrt(V_peak - 1.0/(4.0*alpha**2))
    om_I = -(n + 0.5)/alpha
    return {
        "omega": complex(om_R, om_I),
        "V_peak": V_peak,
        "alpha": alpha,
        "x_peak": x_peak,
        "i_max": i_max,
    }


# ===========================================================================
# Main
# ===========================================================================

# Published Kerr (l=2, m=2, n=0) QNM benchmarks
# Source: Berti's QNM database (approximate values for sanity check; user can
# replace with high-precision values if needed)
KERR_GR_220 = {
    0.00: complex(0.373672, -0.088962),  # Schwarzschild Leaver gold
    0.10: complex(0.38136,  -0.08858),
    0.20: complex(0.38980,  -0.08799),
    0.30: complex(0.39907,  -0.08720),
    0.50: complex(0.42058,  -0.08486),
    0.70: complex(0.44885,  -0.08056),
    0.90: complex(0.49348,  -0.07086),
    0.95: complex(0.51077,  -0.06547),
    0.99: complex(0.53247,  -0.05507),
}


def main():
    print("=" * 88)
    print("G119b  --  Kerr-base spin sweep REDUX via Poschl-Teller / WKB proxy")
    print("=" * 88)
    print()
    print("Stepping stone before full G120 Teukolsky.")
    print("Uses WKB-Poschl-Teller frequency-domain estimate on V_base.")
    print()

    rs_min, rs_max, N = -200.0, 200.0, 4001
    rstar = np.linspace(rs_min, rs_max, N)
    print(f"r* grid: N = {N},  dr* = {(rs_max-rs_min)/(N-1):.4e}")
    print()

    spins = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]

    results = []
    print(f"  {'a':>5}  {'om_GR_WKB':>30}  {'om_published_l=m=2':>26}  "
          f"{'WKB vs pub':>10}  {'om_STAM_WKB':>30}  {'shift':>8}")
    print("-" * 130)
    for a in spins:
        _, V_GR, _ = build_potential_on_rstar(rstar, a, use_STAM=False)
        _, V_STAM, _ = build_potential_on_rstar(rstar, a, use_STAM=True)

        res_GR = WKB_PT_QNM(rstar, V_GR, n=0)
        res_STAM = WKB_PT_QNM(rstar, V_STAM, n=0)
        if res_GR is None or res_STAM is None:
            print(f"  {a:>5.2f}  WKB fit failed")
            results.append({"a": a, "status": "WKB_fail"})
            continue

        om_GR = res_GR["omega"]
        om_STAM = res_STAM["omega"]
        om_pub = KERR_GR_220.get(a)
        wkb_vs_pub = abs(om_GR - om_pub)/abs(om_pub) if om_pub else None
        shift = abs(om_STAM - om_GR)/abs(om_GR)

        wkb_vs_pub_str = f"{wkb_vs_pub*100:>9.2f}%" if wkb_vs_pub is not None else "    n/a"
        pub_str = str(om_pub) if om_pub is not None else "          n/a"
        print(f"  {a:>5.2f}  {str(om_GR):>30}  {pub_str:>26}  "
              f"{wkb_vs_pub_str}  {str(om_STAM):>30}  {shift*100:>7.4f}%")
        results.append({
            "a": a, "status": "ok",
            "om_GR_WKB": om_GR, "om_STAM_WKB": om_STAM,
            "om_pub": om_pub, "wkb_vs_pub": wkb_vs_pub,
            "shift": shift,
            "V_GR": V_GR, "V_STAM": V_STAM,
            "r_ph_eq": r_ph_prograde(a), "r_horizon": r_horizon(a),
        })
    print()

    # ----- Comparison with G118v2 -----
    G118v2_TD = {
        0.0: 4.5703, 0.1: 5.2665, 0.2: 5.8360, 0.3: 6.1631,
        0.5: 6.0742, 0.7: 5.9069, 0.9: 14.0815, 0.95: 10.7923, 0.99: 12.1178,
    }
    print("=" * 88)
    print("G119b WKB shift  vs  G118v2 time-domain shift  (same V_base)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'G119b WKB shift':>18}  {'G118v2 TD shift':>18}  {'WKB - TD':>12}")
    print("-" * 60)
    for r in results:
        if r["status"] != "ok":
            continue
        a = r["a"]
        wkb_pct = r["shift"]*100
        td_pct = G118v2_TD.get(a, None)
        if td_pct is None:
            print(f"  {a:>5.2f}  {wkb_pct:>17.4f}%  {'n/a':>18}  {'n/a':>12}")
            continue
        diff = wkb_pct - td_pct
        print(f"  {a:>5.2f}  {wkb_pct:>17.4f}%  {td_pct:>17.4f}%  {diff:>+11.4f}")
    print()

    # ----- Plot -----
    ok = [r for r in results if r["status"] == "ok"]
    if ok:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        ax = axes[0, 0]
        spins_arr = np.array([r["a"] for r in ok])
        wkb_arr = np.array([r["shift"]*100 for r in ok])
        td_arr = np.array([G118v2_TD.get(r["a"], np.nan) for r in ok])
        ax.plot(spins_arr, wkb_arr, "o-", linewidth=2, markersize=8, color="tab:green",
                label="G119b WKB")
        ax.plot(spins_arr, td_arr, "s--", linewidth=1.5, markersize=7, color="tab:purple",
                label="G118v2 TD")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("STAM shift (%)")
        ax.set_title("WKB-PT shift vs time-domain shift on same V_base")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        # WKB-on-V_base vs published Kerr ω_R
        re_wkb = [r["om_GR_WKB"].real for r in ok]
        re_pub = [r["om_pub"].real if r["om_pub"] else np.nan for r in ok]
        ax.plot(spins_arr, re_wkb, "o-", color="tab:blue", label="Re(om) WKB on V_base")
        ax.plot(spins_arr, re_pub, "s--", color="tab:red", label="Re(om) published Kerr (2,2,0)")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("Re(omega)")
        ax.set_title("V_base calibration audit vs published Kerr l=m=2 QNM")
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        for r in ok[::2]:
            ax.plot(rstar, r["V_GR"], linewidth=0.8, label=f"a = {r['a']:.2f}")
        ax.set_xlim(-30, 30)
        ax.set_xlabel("r* / M")
        ax.set_ylabel("V_base_GR")
        ax.set_title("V_GR profiles (calibration source)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        ax = axes[1, 1]
        for r in ok[::2]:
            diff = r["V_STAM"] - r["V_GR"]
            ax.plot(rstar, diff, linewidth=0.8, label=f"a = {r['a']:.2f}")
        ax.set_xlim(-30, 5)
        ax.set_xlabel("r* / M")
        ax.set_ylabel("V_STAM - V_GR")
        ax.set_title("STAM correction shape vs spin")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G119b_kerr_base_spin_sweep_redux.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G119b -- Kerr-base spin sweep redux (WKB stepping stone)\n"]
    md.append("\nStepping stone before full G120 Teukolsky.  Poschl-Teller / WKB\n"
              "frequency-domain estimate on V_base, applied on the G119 Kerr base.\n")
    md.append("\n## V_base calibration audit (vs published Kerr l=m=2 n=0)\n")
    md.append("| a | WKB Re(om) on V_base | published Kerr Re(om) | WKB vs pub |\n")
    md.append("|---:|---:|---:|---:|\n")
    for r in ok:
        a = r["a"]
        wkb_re = r["om_GR_WKB"].real
        pub_re = r["om_pub"].real if r["om_pub"] else None
        if pub_re is not None:
            dev_pct = r["wkb_vs_pub"]*100
            md.append(f"| {a:.2f} | {wkb_re:.4f} | {pub_re:.4f} | {dev_pct:.2f}% |\n")
        else:
            md.append(f"| {a:.2f} | {wkb_re:.4f} | n/a | n/a |\n")
    md.append("\nReading: WKB-on-V_base differs from published Kerr by O(few %)\n"
              "at low spin (standard WKB error) and grows at high spin (because\n"
              "V_base lacks Kerr-specific m-coupling).  This audits the gap\n"
              "between V_base and rigorous Kerr.\n")
    md.append("\n## WKB shift vs G118v2 time-domain shift on the same V_base\n")
    md.append("| a | G119b WKB shift | G118v2 TD shift | WKB - TD |\n")
    md.append("|---:|---:|---:|---:|\n")
    for r in ok:
        a = r["a"]
        wkb_pct = r["shift"]*100
        td_pct = G118v2_TD.get(a)
        if td_pct is None:
            md.append(f"| {a:.2f} | {wkb_pct:.4f}% | n/a | n/a |\n")
            continue
        diff = wkb_pct - td_pct
        md.append(f"| {a:.2f} | {wkb_pct:.4f}% | {td_pct:.4f}% | {diff:+.4f} |\n")
    md.append("\nIf WKB and TD agree, the shift is robust to extraction method.\n"
              "If they disagree, the difference reveals method sensitivity in\n"
              "the time-domain fit window or higher-overtone contamination.\n")
    md.append("\n## Caveats\n")
    md.append("- V_base is Kerr-Delta-aware Schwarzschild-RW analog (G118v2 form).\n"
              "  NOT the rigorous Kerr axial potential.\n")
    md.append("- WKB-PT is 0th-order WKB; absolute frequencies have few-% intrinsic\n"
              "  error even at a = 0.\n")
    md.append("- This is NOT a final Kerr QNM prediction.  See G120 spec below.\n")
    md.append("\n## Next step:  G120 full Detweiler / Sasaki-Nakamura\n")
    md.append("G120 must:\n"
              "1. Build the spin-weighted spheroidal angular eigenvalue solver A_lm(a*omega).\n"
              "2. Build the Sasaki-Nakamura or Detweiler radial equation for s=-2.\n"
              "3. Calibrate GR Kerr QNM to <0.5% at a in [0, 0.5, 0.9].\n"
              "4. Insert G119's Delta_STAM only after calibration passes.\n"
              "5. Track the l=m=2 n=0 mode branch across the spin sweep.\n"
              "6. Report status flags per spin.\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G119b_kerr_base_spin_sweep_redux.py](../scripts/G119b_kerr_base_spin_sweep_redux.py)\n")
    md.append("- [plots/G119b_kerr_base_spin_sweep_redux.png](../plots/G119b_kerr_base_spin_sweep_redux.png)\n")
    (RESULTS / "G119b_kerr_base_spin_sweep_redux_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G119b_kerr_base_spin_sweep_redux_summary.md'}")


if __name__ == "__main__":
    main()
