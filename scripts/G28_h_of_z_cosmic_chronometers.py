#!/usr/bin/env python3
"""
G28_h_of_z_cosmic_chronometers.py

Model-independent H(z) test against cosmic chronometer data.

Background. Cosmic chronometers measure H(z) directly from differential
ages of passively-evolving galaxies, without going through distance-modulus
fitting. This makes them the cleanest available test of any cosmological
H(z) prediction — they do not carry the LCDM-bias that distance-based
measurements do.

Models compared at the same chronometer redshifts:
  1. STAM intrinsic H(z) = H_0 (1+z)^2 / (1 + z + 0.5 z^2), at H_0 = 73.
     This is the "no-bridge" form from the 2026-05-11 distance reframe:
     STAM's intrinsic Hubble rate, with the photon-A traversal bias removed.
  2. LCDM at H_0 = 73.04, Omega_m = 0.315 (SH0ES H_0, Planck content).
  3. LCDM at H_0 = 67.4, Omega_m = 0.315 (Planck H_0, Planck content).
  4. Einstein-de Sitter (matter-only), H_0 = 73.

Chi-squared:
  chi^2 = sum_i [(H_obs(z_i) - H_model(z_i)) / sigma_i]^2

Cosmic chronometer data hard-coded from published compilations:
  - Simon et al. 2005
  - Stern et al. 2010
  - Moresco et al. 2012
  - Zhang et al. 2014
  - Moresco 2015
  - Moresco et al. 2016
  - Ratsimbazafy et al. 2017
These are the standard 31 CC measurements used in the literature.

Author: Sean Brady / STAM Model-A
Date: 2026-05-11
STAM-only frame of mind. Clues, not failures.
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Cosmic chronometer data (z, H(z) in km/s/Mpc, sigma)
# ============================================================================
# Standard compilation. References indicated.
CC_DATA = np.array(
    [
        (0.07, 69.0, 19.6),     # Zhang+ 2014
        (0.09, 69.0, 12.0),     # Simon+ 2005
        (0.12, 68.6, 26.2),     # Zhang+ 2014
        (0.17, 83.0, 8.0),      # Simon+ 2005
        (0.179, 75.0, 4.0),     # Moresco+ 2012
        (0.199, 75.0, 5.0),     # Moresco+ 2012
        (0.20, 72.9, 29.6),     # Zhang+ 2014
        (0.27, 77.0, 14.0),     # Simon+ 2005
        (0.28, 88.8, 36.6),     # Zhang+ 2014
        (0.352, 83.0, 14.0),    # Moresco+ 2012
        (0.3802, 83.0, 13.5),   # Moresco+ 2016
        (0.40, 95.0, 17.0),     # Simon+ 2005
        (0.4004, 77.0, 10.2),   # Moresco+ 2016
        (0.4247, 87.1, 11.2),   # Moresco+ 2016
        (0.4497, 92.8, 12.9),   # Moresco+ 2016
        (0.47, 89.0, 49.6),     # Ratsimbazafy+ 2017
        (0.4783, 80.9, 9.0),    # Moresco+ 2016
        (0.48, 97.0, 62.0),     # Stern+ 2010
        (0.593, 104.0, 13.0),   # Moresco+ 2012
        (0.68, 92.0, 8.0),      # Moresco+ 2012
        (0.781, 105.0, 12.0),   # Moresco+ 2012
        (0.875, 125.0, 17.0),   # Moresco+ 2012
        (0.88, 90.0, 40.0),     # Stern+ 2010
        (0.90, 117.0, 23.0),    # Simon+ 2005
        (1.037, 154.0, 20.0),   # Moresco+ 2012
        (1.30, 168.0, 17.0),    # Simon+ 2005
        (1.363, 160.0, 33.6),   # Moresco 2015
        (1.43, 177.0, 18.0),    # Simon+ 2005
        (1.53, 140.0, 14.0),    # Simon+ 2005
        (1.75, 202.0, 40.0),    # Simon+ 2005
        (1.965, 186.5, 50.4),   # Moresco 2015
    ],
    dtype=float,
)


# ============================================================================
# Models
# ============================================================================
def H_stam_intrinsic(z: np.ndarray, H0: float) -> np.ndarray:
    """STAM intrinsic H(z), 'no-bridge' form. Bounded by 2*H0 as z -> infinity."""
    return H0 * (1.0 + z) ** 2 / (1.0 + z + 0.5 * z**2)


def H_lcdm(z: np.ndarray, H0: float, Om: float = 0.315) -> np.ndarray:
    """Standard flat LCDM with matter + dark energy (Omega_Lambda = 1 - Om)."""
    return H0 * np.sqrt(Om * (1.0 + z) ** 3 + (1.0 - Om))


def H_eds(z: np.ndarray, H0: float) -> np.ndarray:
    """Einstein-de Sitter, matter-only flat universe."""
    return H0 * (1.0 + z) ** 1.5


# ============================================================================
# chi-squared
# ============================================================================
def chi2(z: np.ndarray, Hobs: np.ndarray, sigma: np.ndarray, model_fn, *args) -> float:
    Hmodel = model_fn(z, *args)
    return float(np.sum(((Hobs - Hmodel) / sigma) ** 2))


def main():
    print("=" * 78)
    print("G28: H(z) vs cosmic chronometers — model-independent distance reframe test")
    print("=" * 78)
    print()

    z = CC_DATA[:, 0]
    Hobs = CC_DATA[:, 1]
    sigma = CC_DATA[:, 2]
    N = len(z)

    print(f"Cosmic chronometer compilation: {N} measurements, z in [{z.min():.3f}, {z.max():.3f}]")
    print()

    # ============================================================
    # Compute chi^2 for each model
    # ============================================================
    H0_local = 73.04   # SH0ES local H_0
    H0_planck = 67.4   # Planck H_0

    chi2_stam = chi2(z, Hobs, sigma, H_stam_intrinsic, H0_local)
    chi2_lcdm_73 = chi2(z, Hobs, sigma, H_lcdm, H0_local, 0.315)
    chi2_lcdm_67 = chi2(z, Hobs, sigma, H_lcdm, H0_planck, 0.315)
    chi2_eds = chi2(z, Hobs, sigma, H_eds, H0_local)

    dof = N - 1  # one parameter (H_0) at most; for STAM/EdS no free fit here

    print("Models (no free fitting — all H_0 values fixed a priori):")
    print(f"  1. STAM intrinsic H(z), H_0 = 73.04           chi^2 = {chi2_stam:8.2f}   chi^2/N = {chi2_stam/N:.3f}")
    print(f"  2. LCDM (Om=0.315),     H_0 = 73.04           chi^2 = {chi2_lcdm_73:8.2f}   chi^2/N = {chi2_lcdm_73/N:.3f}")
    print(f"  3. LCDM (Om=0.315),     H_0 = 67.4 (Planck)   chi^2 = {chi2_lcdm_67:8.2f}   chi^2/N = {chi2_lcdm_67/N:.3f}")
    print(f"  4. EdS (matter-only),   H_0 = 73.04           chi^2 = {chi2_eds:8.2f}   chi^2/N = {chi2_eds/N:.3f}")
    print()

    # ============================================================
    # Best-fit H_0 for each model (single-parameter fit)
    # ============================================================
    # For each shape, find H_0 that minimizes chi^2.
    # chi^2(H_0) = sum [(Hobs - H_0 * f(z))/sigma]^2
    # d/dH_0 = 0  =>  H_0 = sum[Hobs*f(z)/sigma^2] / sum[f(z)^2/sigma^2]

    def best_fit_H0(z, Hobs, sigma, shape_fn):
        """Find H_0 minimizing chi^2 for a model where H(z) = H_0 * shape(z)."""
        f = shape_fn(z, 1.0)  # shape only
        num = np.sum(Hobs * f / sigma**2)
        den = np.sum(f**2 / sigma**2)
        H0_best = num / den
        chi2_best = chi2(z, Hobs, sigma, lambda zz, h: h * shape_fn(zz, 1.0), H0_best)
        return H0_best, chi2_best

    def stam_shape(z, H0):
        return H_stam_intrinsic(z, H0)

    def lcdm_shape(z, H0):
        return H_lcdm(z, H0, 0.315)

    def eds_shape(z, H0):
        return H_eds(z, H0)

    H0_stam_bf, chi2_stam_bf = best_fit_H0(z, Hobs, sigma, stam_shape)
    H0_lcdm_bf, chi2_lcdm_bf = best_fit_H0(z, Hobs, sigma, lcdm_shape)
    H0_eds_bf, chi2_eds_bf = best_fit_H0(z, Hobs, sigma, eds_shape)

    print("Best-fit H_0 (each model fitted as one-parameter):")
    print(f"  1. STAM intrinsic shape:   H_0_best = {H0_stam_bf:6.2f}   chi^2 = {chi2_stam_bf:7.2f}   chi^2/dof = {chi2_stam_bf/(N-1):.3f}")
    print(f"  2. LCDM (Om=0.315) shape:  H_0_best = {H0_lcdm_bf:6.2f}   chi^2 = {chi2_lcdm_bf:7.2f}   chi^2/dof = {chi2_lcdm_bf/(N-1):.3f}")
    print(f"  3. EdS shape:              H_0_best = {H0_eds_bf:6.2f}   chi^2 = {chi2_eds_bf:7.2f}   chi^2/dof = {chi2_eds_bf/(N-1):.3f}")
    print()

    # ============================================================
    # Plot
    # ============================================================
    z_grid = np.linspace(0, 2.1, 400)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel 1: data + a priori (no-fit) models
    ax1.errorbar(z, Hobs, yerr=sigma, fmt="ko", markersize=4, capsize=2, label="Cosmic chronometers")
    ax1.plot(z_grid, H_stam_intrinsic(z_grid, H0_local), "b-", lw=2, label=f"STAM intrinsic, H_0=73.04")
    ax1.plot(z_grid, H_lcdm(z_grid, H0_local, 0.315), "r--", lw=2, label="LCDM, H_0=73.04, Om=0.315")
    ax1.plot(z_grid, H_lcdm(z_grid, H0_planck, 0.315), "g:", lw=2, label="LCDM, H_0=67.4 (Planck)")
    ax1.plot(z_grid, H_eds(z_grid, H0_local), "m-.", lw=1.5, label="EdS, H_0=73.04")
    ax1.set_xlabel("redshift z")
    ax1.set_ylabel("H(z)  [km/s/Mpc]")
    ax1.set_title("A priori models (no free fit)")
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 2.1)
    ax1.set_ylim(50, 260)

    # Panel 2: best-fit (single-parameter) versions
    ax2.errorbar(z, Hobs, yerr=sigma, fmt="ko", markersize=4, capsize=2, label="Cosmic chronometers")
    ax2.plot(z_grid, H_stam_intrinsic(z_grid, H0_stam_bf), "b-", lw=2, label=f"STAM intrinsic, H_0={H0_stam_bf:.1f} (best fit)")
    ax2.plot(z_grid, H_lcdm(z_grid, H0_lcdm_bf, 0.315), "r--", lw=2, label=f"LCDM, H_0={H0_lcdm_bf:.1f} (best fit)")
    ax2.plot(z_grid, H_eds(z_grid, H0_eds_bf), "m-.", lw=1.5, label=f"EdS, H_0={H0_eds_bf:.1f} (best fit)")
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("H(z)  [km/s/Mpc]")
    ax2.set_title("Best-fit H_0 (one-parameter fits)")
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 2.1)
    ax2.set_ylim(50, 260)

    fig.suptitle("G28: STAM intrinsic vs LCDM vs EdS against cosmic chronometers")
    fig.tight_layout()

    plot_path = RESULTS / "G28_h_of_z_cosmic_chronometers.png"
    fig.savefig(plot_path, dpi=140, bbox_inches="tight")
    print(f"Plot saved: {plot_path}")
    print()

    # ============================================================
    # Per-point residuals at high z
    # ============================================================
    print("Per-point comparison (high-z subset, z > 1.0):")
    print(f"  {'z':>6} {'H_obs':>8} {'sigma':>7}   {'STAM_pred':>10} {'res_STAM':>9}   {'LCDM73_pred':>11} {'res_LCDM':>10}")
    for zi, Hi, si in CC_DATA:
        if zi > 1.0:
            Hs = float(H_stam_intrinsic(np.array([zi]), H0_local)[0])
            Hl = float(H_lcdm(np.array([zi]), H0_local, 0.315)[0])
            res_s = (Hi - Hs) / si
            res_l = (Hi - Hl) / si
            print(f"  {zi:6.3f} {Hi:8.1f} {si:7.1f}   {Hs:10.1f} {res_s:9.2f}   {Hl:11.1f} {res_l:10.2f}")
    print()

    # ============================================================
    # Summary verdict
    # ============================================================
    print("Verdict summary:")
    print(f"  STAM intrinsic (H_0=73.04, no fit):  chi^2/N = {chi2_stam/N:.3f}")
    print(f"  STAM intrinsic (best-fit H_0):       chi^2/dof = {chi2_stam_bf/(N-1):.3f},  H_0 = {H0_stam_bf:.2f}")
    print(f"  LCDM-at-73 (no fit):                 chi^2/N = {chi2_lcdm_73/N:.3f}")
    print(f"  LCDM-at-67.4 (Planck, no fit):       chi^2/N = {chi2_lcdm_67/N:.3f}")
    print(f"  LCDM best-fit H_0:                   chi^2/dof = {chi2_lcdm_bf/(N-1):.3f},  H_0 = {H0_lcdm_bf:.2f}")
    print(f"  EdS best-fit H_0:                    chi^2/dof = {chi2_eds_bf/(N-1):.3f},  H_0 = {H0_eds_bf:.2f}")

    # Write markdown summary
    md = []
    md.append("# G28: H(z) vs cosmic chronometers — model-independent distance reframe test")
    md.append("")
    md.append("**Date:** 2026-05-11")
    md.append("")
    md.append("## Goal")
    md.append("")
    md.append(
        "Test STAM's intrinsic H(z) prediction against cosmic chronometer measurements. "
        "Chronometers measure H(z) from differential galaxy ages, model-independent of any "
        "distance-modulus fit — they do not carry the LCDM-bias the bridge term encodes. "
        "If STAM's intrinsic H(z) fits chronometer data well, the distance reframe (STAM "
        "intrinsic is right; what looks like dark energy in distance fits is LCDM measurement "
        "bias) gains empirical support."
    )
    md.append("")
    md.append("## STAM intrinsic H(z)")
    md.append("")
    md.append("```")
    md.append("H(z) = H_0 * (1+z)^2 / (1 + z + 0.5 z^2)")
    md.append("```")
    md.append("")
    md.append("Limits: H(0) = H_0; H(z) -> 2*H_0 as z -> infinity. Bounded.")
    md.append("")
    md.append("## Data")
    md.append("")
    md.append(f"{N} cosmic chronometer measurements, z in [{z.min():.3f}, {z.max():.3f}].")
    md.append("Standard compilation: Simon+05, Stern+10, Moresco+12, Zhang+14, Moresco 15, Moresco+16, Ratsimbazafy+17.")
    md.append("")
    md.append("## Results — a priori (no free fit)")
    md.append("")
    md.append("| Model | H_0 fixed | chi^2 | chi^2/N |")
    md.append("|---|---|---|---|")
    md.append(f"| STAM intrinsic | 73.04 | {chi2_stam:.2f} | {chi2_stam/N:.3f} |")
    md.append(f"| LCDM (Om=0.315) | 73.04 | {chi2_lcdm_73:.2f} | {chi2_lcdm_73/N:.3f} |")
    md.append(f"| LCDM (Om=0.315) | 67.4 (Planck) | {chi2_lcdm_67:.2f} | {chi2_lcdm_67/N:.3f} |")
    md.append(f"| EdS (matter-only) | 73.04 | {chi2_eds:.2f} | {chi2_eds/N:.3f} |")
    md.append("")
    md.append("## Results — single-parameter best-fit H_0")
    md.append("")
    md.append("| Model | best-fit H_0 | chi^2 | chi^2/dof |")
    md.append("|---|---|---|---|")
    md.append(f"| STAM intrinsic | {H0_stam_bf:.2f} | {chi2_stam_bf:.2f} | {chi2_stam_bf/(N-1):.3f} |")
    md.append(f"| LCDM (Om=0.315) | {H0_lcdm_bf:.2f} | {chi2_lcdm_bf:.2f} | {chi2_lcdm_bf/(N-1):.3f} |")
    md.append(f"| EdS | {H0_eds_bf:.2f} | {chi2_eds_bf:.2f} | {chi2_eds_bf/(N-1):.3f} |")
    md.append("")
    md.append("## Plot")
    md.append("")
    md.append("![H(z) vs chronometers](G28_h_of_z_cosmic_chronometers.png)")
    md.append("")
    md.append("Left panel: a priori predictions with H_0 fixed at SH0ES (73.04) or Planck (67.4).")
    md.append("Right panel: each model fitted with its own best-fit H_0 (one free parameter).")
    md.append("")
    md.append("## High-z residual table (z > 1.0)")
    md.append("")
    md.append("| z | H_obs | sigma | STAM_pred | (H_obs-STAM)/sigma | LCDM_pred (H_0=73) | (H_obs-LCDM)/sigma |")
    md.append("|---|---|---|---|---|---|---|")
    for zi, Hi, si in CC_DATA:
        if zi > 1.0:
            Hs = float(H_stam_intrinsic(np.array([zi]), H0_local)[0])
            Hl = float(H_lcdm(np.array([zi]), H0_local, 0.315)[0])
            res_s = (Hi - Hs) / si
            res_l = (Hi - Hl) / si
            md.append(f"| {zi:.3f} | {Hi:.1f} | {si:.1f} | {Hs:.1f} | {res_s:+.2f} | {Hl:.1f} | {res_l:+.2f} |")
    md.append("")

    # Reading
    md.append("## Reading")
    md.append("")
    if chi2_stam_bf / (N - 1) < 1.0 and chi2_lcdm_bf / (N - 1) < 1.0:
        md.append(
            "Both STAM intrinsic and LCDM are statistically acceptable fits (chi^2/dof < 1). "
            "The test does not strongly discriminate between them at current chronometer precision."
        )
    elif chi2_stam_bf < chi2_lcdm_bf:
        md.append(
            "STAM intrinsic gives a lower chi^2 than LCDM under one-parameter fits. "
            "Empirical support for the distance reframe at chronometer precision."
        )
    elif chi2_stam_bf > chi2_lcdm_bf * 1.2:
        md.append(
            "STAM intrinsic gives a substantially higher chi^2 than LCDM. "
            "The intrinsic H(z) form may need revision — the 2*H_0 high-z bound is the natural suspect."
        )
    else:
        md.append(
            "STAM intrinsic and LCDM give comparable chi^2 under one-parameter fits. "
            "Neither is strongly preferred at current chronometer precision."
        )
    md.append("")
    md.append("**Open questions raised by this run:**")
    md.append("")
    md.append("- The STAM intrinsic formula bounds H(z) <= 2*H_0. At z = 1.965, observed H ~ 186.5 +/- 50.4.")
    md.append("  With H_0 = 73, the STAM ceiling is 146 — within 1 sigma of the observation but on the low side.")
    md.append("  At higher z (if data existed), the formula would predict an even larger deficit.")
    md.append("  Possibility: the closed-form formula is a low-to-moderate-z approximation; the underlying STAM-intrinsic dynamics may differ at z > 1.")
    md.append("- Best-fit H_0 under STAM-intrinsic shape gives a different value than 73 or 67.4.")
    md.append("  This is itself a framework-internal result — the value tells us what H_0 STAM-intrinsic prefers if the formula is taken literally.")
    md.append("- The V_3 modified Friedmann redo is the natural follow-up: a derived H(z) from V_3 dynamics, not an ansatz.")
    md.append("")

    md_path = RESULTS / "G28_h_of_z_cosmic_chronometers_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Summary saved: {md_path}")


if __name__ == "__main__":
    main()
