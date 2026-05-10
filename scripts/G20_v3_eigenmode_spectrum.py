#!/usr/bin/env python3
"""
G20_v3_eigenmode_spectrum.py

Compute the eigenmode spectrum of A as a coordinate on (A_0, 1) with
V_3 potential and Dirichlet boundary conditions at A_0 and A=1.

Setup:
    Schrödinger-like equation on the bounded interval:
        -(1/2) ψ''(A) + V_3(A) ψ(A) = E ψ(A)
    with ψ(A_0) = ψ(1) = 0 (hard walls at structural boundaries).

    V_3(A) = α/A + β/(1-A) with α/β = [A_0/(1-A_0)]² so V_3 minimum at A_0.
    A_0 = 1/(12π).
    β = 1 (overall scale; eigenvalue ratios are invariant).
    Kinetic-term coefficient = 1/2 (standard).

What we compute:
    - Lowest ~10 eigenvalues E_n
    - Eigenfunctions ψ_n(A) for low n
    - Zero crossings of each ψ_n
    - Energy ratios E_n/E_0

What we compare to:
    - Orbital thirds: A = 1/3, 2/3, 1 (strong-field special values)
    - Pure particle-in-box: zero crossings at A_0 + k(1-A_0)/n
    - Harmonic-oscillator approximation around V_3 minimum:
      E_n - V_min = (n + 1/2) ω_HO with ω_HO = √V_3''(A_0)

Question: does V_3's eigenmode spectrum on (A_0, 1) produce zero crossings
that match the orbital thirds, or other framework-special values? If yes,
the harmony reading has structural content. If no, the orbital thirds
emerge from independent physics (orbital mechanics in (1-A) potential)
and aren't connected to V_3's mode structure.

Honest framing: this is treating A as a coordinate (Schrödinger setup),
which is non-standard but physically motivated by the framework's
"A bounded between A_0 and 1" structural commitment. The kinetic-term
normalization (1/2) is conventional; absolute eigenvalue scale depends
on the choice of effective mass, but ratios and zero-crossing locations
are invariant.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print("G20: V_3 eigenmode spectrum on (A_0, 1)")
    print("=" * 78)
    print()

    A_0 = 1.0 / (12.0 * np.pi)
    alpha_over_beta = (A_0 / (1.0 - A_0))**2
    beta = 1.0
    alpha = beta * alpha_over_beta

    print(f"A_0 = 1/(12π) = {A_0:.6f}")
    print(f"α/β = [A_0/(1-A_0)]² = {alpha_over_beta:.6e}")
    print(f"β = {beta} (sets overall energy scale)")
    print()

    # Grid on (A_0, 1) with hard walls at endpoints
    # Interior grid (Dirichlet boundary conditions ψ(A_0) = ψ(1) = 0)
    N = 1200
    A_full = np.linspace(A_0, 1.0, N + 2)
    A = A_full[1:-1]
    dA = A[1] - A[0]

    print(f"Grid: N = {N} interior points on ({A_0:.4f}, 1.0)")
    print(f"Spacing dA = {dA:.6f}")
    print()

    # V_3 evaluated on grid
    V = alpha / A + beta / (1.0 - A)
    V_min_at_A0 = beta / (1.0 - A_0)**2  # analytic V_3(A_0)
    V_pp_at_A0 = 2.0 * beta / (A_0 * (1.0 - A_0)**3)  # analytic V_3''(A_0)

    print(f"Analytic V_3(A_0) = β/(1-A_0)² = {V_min_at_A0:.6f}")
    print(f"Analytic V_3''(A_0) = 2β/[A_0(1-A_0)³] = {V_pp_at_A0:.4f}")
    print(f"HO approximation ω_HO = √V_3''(A_0) = {np.sqrt(V_pp_at_A0):.4f}")
    print()

    # Build tridiagonal Hamiltonian: -(1/2) ψ'' + V ψ = E ψ
    # Finite difference: -(1/2) × (-ψ_{i-1} + 2ψ_i - ψ_{i+1}) / dA²
    # → diagonal coefficient 1/dA² + V_i, off-diagonal -1/(2 dA²)
    diag = 1.0 / dA**2 + V
    off = -1.0 / (2.0 * dA**2) * np.ones(N - 1)

    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)

    print("Diagonalizing (this may take a few seconds)...")
    E, psi = np.linalg.eigh(H)
    print(f"Done. Computed {len(E)} eigenvalues.")
    print()

    # --- Eigenvalue spectrum ---
    print("=" * 78)
    print("EIGENVALUE SPECTRUM (lowest 10)")
    print("=" * 78)
    print()
    print(f"{'n':>3s}  {'E_n':>14s}  {'E_n - V_min':>14s}  "
          f"{'(E_n - V_min) / ω_HO':>22s}  {'E_n / E_0':>10s}")
    print("-" * 78)
    omega_HO = np.sqrt(V_pp_at_A0)
    for n in range(10):
        delta_E = E[n] - V_min_at_A0
        ho_units = delta_E / omega_HO
        ratio = E[n] / E[0]
        print(f"{n:>3d}  {E[n]:>14.6f}  {delta_E:>14.6f}  "
              f"{ho_units:>22.4f}  {ratio:>10.4f}")
    print()
    print("HO comparison: for pure harmonic oscillator, (E_n - V_min)/ω_HO = n + 1/2")
    print("Pure particle-in-box on (A_0, 1) (no potential): E_n ∝ n²")
    print()

    # --- Zero crossings of low modes ---
    print("=" * 78)
    print("ZERO CROSSINGS OF LOW MODES")
    print("=" * 78)
    print()
    print("Mode n has n zero crossings (counting endpoints).")
    print("Interior zero crossings — compare to framework-special values:")
    print(f"  Orbital ISCO:         A = 1/3 = {1/3:.4f}")
    print(f"  Orbital photon sphere: A = 2/3 = {2/3:.4f}")
    print(f"  A_0 (left wall):      A = {A_0:.4f}")
    print(f"  Right wall:           A = 1.0")
    print()

    crossings_per_mode = []
    for n in range(8):
        psi_n = psi[:, n]
        # Find sign changes (zero crossings)
        signs = np.sign(psi_n)
        sign_diffs = np.diff(signs)
        crossing_indices = np.where(sign_diffs != 0)[0]
        # Linear interpolation for crossing location
        crossings = []
        for idx in crossing_indices:
            if psi_n[idx + 1] - psi_n[idx] != 0:
                frac = -psi_n[idx] / (psi_n[idx + 1] - psi_n[idx])
                A_cross = A[idx] + frac * (A[idx + 1] - A[idx])
                crossings.append(A_cross)
        crossings_per_mode.append(crossings)
        if crossings:
            crossings_str = ", ".join(f"{c:.4f}" for c in crossings)
            print(f"  Mode {n} ({n} interior crossings): {crossings_str}")
        else:
            print(f"  Mode {n} (0 interior crossings): -")
    print()

    # --- Special-value matching: do any modes have crossings near orbital thirds? ---
    print("=" * 78)
    print("CHECK: do any modes have zero crossings near orbital thirds (1/3, 2/3)?")
    print("=" * 78)
    print()
    target_thirds = [1.0/3.0, 2.0/3.0]
    for n in range(8):
        crossings = crossings_per_mode[n]
        if not crossings:
            continue
        for target in target_thirds:
            distances = [abs(c - target) for c in crossings]
            min_dist = min(distances)
            if min_dist < 0.05:
                closest = crossings[distances.index(min_dist)]
                print(f"  Mode {n}: crossing at A = {closest:.4f} "
                      f"(target {target:.4f}, deviation {(closest-target)*100:+.2f}%)")
    print()

    # --- Plot low modes ---
    print("Plotting low modes...")
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Top: V_3 potential
    A_plot = np.linspace(A_0 + 0.001, 0.999, 500)
    V_plot = alpha / A_plot + beta / (1 - A_plot)
    axes[0].plot(A_plot, V_plot, 'k-', linewidth=2, label='V_3(A)')
    axes[0].axvline(A_0, color='gray', linestyle=':', label=f'A_0 = 1/(12π)')
    axes[0].axvline(1.0, color='gray', linestyle=':', label='A = 1')
    axes[0].axvline(1/3, color='tab:orange', linestyle='--', alpha=0.5, label='A = 1/3')
    axes[0].axvline(2/3, color='tab:orange', linestyle='--', alpha=0.5, label='A = 2/3')
    axes[0].set_xlabel('A')
    axes[0].set_ylabel('V_3(A)')
    axes[0].set_title('V_3 potential on (A_0, 1)')
    axes[0].set_ylim(0, 50)
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3)

    # Bottom: lowest few eigenfunctions
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, 6))
    for n in range(6):
        psi_n = psi[:, n]
        # Normalize for display
        psi_n_normed = psi_n / np.max(np.abs(psi_n))
        offset = -n * 1.5  # vertical offset for stacking
        axes[1].plot(A, psi_n_normed + offset, color=colors[n],
                     linewidth=1.5, label=f'n={n}, E={E[n]:.2f}')
        axes[1].axhline(offset, color=colors[n], linestyle=':', alpha=0.3)

    axes[1].axvline(1/3, color='tab:orange', linestyle='--', alpha=0.5, label='A=1/3')
    axes[1].axvline(2/3, color='tab:orange', linestyle='--', alpha=0.5, label='A=2/3')
    axes[1].axvline(A_0, color='gray', linestyle=':', alpha=0.5)
    axes[1].axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    axes[1].set_xlabel('A')
    axes[1].set_ylabel('ψ_n (offset by 1.5 per mode)')
    axes[1].set_title('Lowest 6 eigenfunctions of V_3 Schrödinger problem on (A_0, 1)')
    axes[1].legend(fontsize=8, loc='upper right', ncol=2)
    axes[1].grid(alpha=0.3)

    out_plot = PLOTS / "G20_v3_eigenmode_spectrum.png"
    plt.tight_layout()
    plt.savefig(out_plot, dpi=200)
    plt.close()
    print(f"Plot: {out_plot}")
    print()

    # --- Verdict ---
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print()

    # Check whether mode 3 (which would have 2 interior crossings) has them at thirds
    mode_3_crossings = crossings_per_mode[3] if len(crossings_per_mode) > 3 else []
    if len(mode_3_crossings) >= 2:
        c1, c2 = mode_3_crossings[0], mode_3_crossings[1]
        third_match_1 = abs(c1 - 1/3) / (1/3) * 100
        third_match_2 = abs(c2 - 2/3) / (2/3) * 100
        print(f"Mode 3 has {len(mode_3_crossings)} interior crossings:")
        print(f"  Crossing 1: {c1:.4f} (vs orbital ISCO at 1/3: {third_match_1:+.2f}%)")
        print(f"  Crossing 2: {c2:.4f} (vs orbital photon sphere at 2/3: {third_match_2:+.2f}%)")
        print()
        if max(third_match_1, third_match_2) < 5:
            print("Mode 3 zero crossings WITHIN 5% of orbital thirds.")
            print("→ Wave-equation eigenmode spectrum produces zero crossings at")
            print("  (approximately) the same locations as the orbital thirds.")
            print("  Structural alignment between V_3 wave dynamics and orbital")
            print("  mechanics in (1-A) metric is real, even though they emerge")
            print("  from independent calculations.")
        elif max(third_match_1, third_match_2) < 20:
            print("Mode 3 zero crossings within 20% of orbital thirds, but not 5%.")
            print("→ Suggestive structural alignment, not exact match. The orbital")
            print("  thirds and V_3 wave-mode crossings are in the same regime but")
            print("  not identical.")
        else:
            print("Mode 3 zero crossings far from orbital thirds.")
            print("→ V_3 wave-equation eigenmodes do NOT match the orbital thirds.")
            print("  These are independent physics; the (1-A) potential's orbital")
            print("  structure is not the same as V_3's wave-mode structure.")
    print()

    # Eigenvalue ratios
    print("Eigenvalue ratios (E_n/E_0) — check for harmonic structure:")
    for n in range(1, 8):
        print(f"  E_{n}/E_0 = {E[n]/E[0]:.4f}")
    print()
    print("For pure harmonic oscillator: E_n/E_0 = (2n+1) (1, 3, 5, 7, 9, 11, 13)")
    print("For pure particle-in-box: E_n/E_0 = (n+1)² (1, 4, 9, 16, 25, ...)")

    write_summary(A_0, V_min_at_A0, V_pp_at_A0, E, crossings_per_mode, [out_plot])


def write_summary(A_0, V_min, V_pp, E, crossings, plots):
    md = []
    md.append("# G20: V_3 Eigenmode Spectrum on (A_0, 1)\n")

    md.append("## Setup\n")
    md.append(
        "Schrödinger-like wave equation for A as a coordinate on the bounded "
        "interval (A_0, 1) with V_3 potential and Dirichlet boundary conditions.\n"
        "\n"
        "```text\n"
        "-(1/2) ψ''(A) + V_3(A) ψ(A) = E ψ(A)\n"
        "ψ(A_0) = ψ(1) = 0     (hard walls at structural boundaries)\n"
        "V_3(A) = α/A + β/(1-A) with α/β = [A_0/(1-A_0)]²\n"
        "A_0 = 1/(12π), β = 1 (sets scale)\n"
        "```\n"
        "\n"
        "Numerical solution: discretize on N=1200 grid points, build "
        "tridiagonal Hamiltonian, diagonalize.\n"
    )

    md.append("## Eigenvalue spectrum (lowest 10)\n")
    md.append("```text")
    md.append(f"{'n':>3s}  {'E_n':>14s}  {'E_n / E_0':>10s}")
    for n in range(10):
        md.append(f"{n:>3d}  {E[n]:>14.6f}  {E[n]/E[0]:>10.4f}")
    md.append("```\n")

    md.append("## Zero crossings\n")
    md.append("Interior zero crossings of low-order modes:\n\n")
    md.append("```text")
    for n in range(8):
        if crossings[n]:
            md.append(f"  Mode {n}: " + ", ".join(f"{c:.4f}" for c in crossings[n]))
        else:
            md.append(f"  Mode {n}: no interior crossings")
    md.append("```\n")

    md.append("Compare to orbital thirds: A = 1/3 ≈ 0.3333, A = 2/3 ≈ 0.6667.\n\n")

    if len(crossings) > 3 and len(crossings[3]) >= 2:
        c1, c2 = crossings[3][0], crossings[3][1]
        match1 = abs(c1 - 1/3) / (1/3) * 100
        match2 = abs(c2 - 2/3) / (2/3) * 100
        md.append(f"**Mode 3 crossings:** A = {c1:.4f} and A = {c2:.4f}\n\n")
        md.append(f"- vs orbital ISCO (A=1/3): {match1:+.2f}% deviation\n")
        md.append(f"- vs orbital photon sphere (A=2/3): {match2:+.2f}% deviation\n\n")

    md.append("## Generated plot\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G20_v3_eigenmode_spectrum_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"\nSummary: {out}")


if __name__ == "__main__":
    main()
