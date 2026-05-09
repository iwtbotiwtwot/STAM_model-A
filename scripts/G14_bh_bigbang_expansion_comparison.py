#!/usr/bin/env python3
"""
G14_bh_bigbang_expansion_comparison.py

Wild-theory exploration: compare the horizon-expansion process during
black hole formation to the scale-factor expansion of the Big Bang.

Conceptual setup (parking-lot wild theory):
    Star collapses -> matter "resets" at a singularity-like state -> space
    inverts -> horizon expands outward driven by initial mass.

    Question: Is this expansion mathematically similar to the Big Bang's
    scale-factor expansion?

In STAM specifically: matter doesn't reach a singularity — it accumulates
on the A=1 bubble surface. The bubble grows from r=0 to r=Rs as matter
joins. The "inside-out" picture: from outside we see the bubble grow;
from "inside" (if such a thing existed; STAM says it doesn't) it might
look like a Big Bang.

What this script computes:
    1. BH horizon growth Rs(t) during free-fall stellar collapse
    2. Big Bang scale factor a(t) for radiation, matter, and dark-energy eras
    3. Functional form comparison: power-law vs exponential vs combined
    4. Total expansion factors and timescales
    5. Plot both on normalized axes for visual comparison

Connections to mainstream theory:
    - Smolin's Cosmological Natural Selection: BH interiors are baby
      universes with re-randomized constants
    - Penrose's CCC (Conformal Cyclic Cosmology): max-extent universe
      maps to next Big Bang via conformal rescaling
    - Hartle-Hawking no-boundary: universe began from a "smooth" initial
      state analogous to a horizon
    - Big Bounce / loop quantum cosmology: BB is the time-reverse of a
      contracting BH

This is parking-lot wild theory. STAM doesn't commit to it. But the
mathematical comparison is well-defined and worth seeing.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# Constants (SI)
G = 6.67430e-11
C = 2.99792458e8
M_SUN = 1.98892e30
H_0_KMS = 73.04 * 1000.0 / (3.086e22)   # H_0 in 1/s


# --- BH horizon growth during stellar collapse ---

def Rs_meters(M_kg: float) -> float:
    return 2.0 * G * M_kg / C ** 2


def free_fall_timescale(R_initial_m: float, M_kg: float) -> float:
    """Free-fall time from R_initial to center: t = (π/2)√(R³/(2GM))"""
    return (np.pi / 2.0) * np.sqrt(R_initial_m ** 3 / (2.0 * G * M_kg))


def horizon_growth_during_collapse(M_kg: float, R_initial_m: float,
                                     n_steps: int = 200) -> dict:
    """Toy model: as star collapses, the horizon emerges and grows from
    r=0 to r=Rs. We model this as the matter falling in from R_initial
    to R=0, with the bubble surface "catching up" to where the matter
    currently is.

    Approximation: at time t, the matter at radius r(t) where r(t) follows
    pressureless free-fall:
        r(t)/R_initial = (cos²(η/2) + ...)

    For a simpler toy: bubble grows linearly from 0 to Rs over time tau_ff.
    """
    Rs = Rs_meters(M_kg)
    tau_ff = free_fall_timescale(R_initial_m, M_kg)

    t_grid = np.linspace(0.0, tau_ff, n_steps)
    # Linear growth (toy):
    R_horizon_linear = Rs * (t_grid / tau_ff)
    # Exponential approach to Rs (more realistic for bubble formation):
    # R(t) = Rs × (1 - exp(-t/tau_ff/0.3))
    R_horizon_exp = Rs * (1.0 - np.exp(-3.0 * t_grid / tau_ff))

    return {
        "M_kg": M_kg,
        "M_solar": M_kg / M_SUN,
        "R_initial_m": R_initial_m,
        "Rs_m": Rs,
        "tau_freefall_s": tau_ff,
        "t_grid_s": t_grid,
        "R_horizon_linear_m": R_horizon_linear,
        "R_horizon_exp_m": R_horizon_exp,
        "expansion_factor": Rs / max(R_initial_m, 1e-30),  # final/initial
    }


# --- Big Bang scale factor evolution ---

def scale_factor_radiation(t_s: float, t_today_s: float = 4.35e17) -> float:
    """a(t) ~ t^(1/2) in radiation era.
    Normalize so a = some reference at z = z_today_eq."""
    return np.sqrt(t_s / t_today_s)


def scale_factor_matter(t_s: float, t_today_s: float = 4.35e17) -> float:
    """a(t) ~ t^(2/3) in matter era."""
    return (t_s / t_today_s) ** (2.0 / 3.0)


def scale_factor_dark_energy(t_s: float, t_today_s: float = 4.35e17,
                              H0: float = H_0_KMS) -> float:
    """a(t) ~ exp(H_0 × (t - t_today)) in dark-energy-dominated era."""
    return np.exp(H0 * (t_s - t_today_s))


def scale_factor_combined(t_s: float, t_eq: float = 1.4e12,
                            t_today: float = 4.35e17,
                            H0: float = H_0_KMS) -> float:
    """Approximate combined evolution: radiation -> matter -> dark energy.
    Returns a(t) normalized to a(t_today) = 1."""
    if t_s < t_eq:
        a_eq = (t_eq / t_today) ** (2.0 / 3.0)
        return a_eq * np.sqrt(t_s / t_eq)
    elif t_s < t_today:
        return (t_s / t_today) ** (2.0 / 3.0)
    else:
        return np.exp(H0 * (t_s - t_today))


# --- comparison ---

def main() -> None:
    print("G14: BH horizon growth vs Big Bang scale factor expansion")
    print("=" * 78)
    print()
    print("Wild-theory exploration. STAM does not commit to BH=baby-universe.")
    print()

    # BH formation for representative cases
    print(">>> BH horizon growth during stellar collapse")
    print()
    bh_cases = [
        ("Stellar BH (10 M_sun, R = 5×10^9 m progenitor)", 10 * M_SUN, 5e9),
        ("Massive BH (40 M_sun, R = 1×10^10 m progenitor)", 40 * M_SUN, 1e10),
        ("PBH formation (10^15 g, horizon-scale at t = 10^-23 s)",
         1e12, 3e-15),
    ]

    bh_results = []
    print(f"  {'Case':<55s}  {'Rs (m)':>12s}  {'tau_ff (s)':>12s}  "
          f"{'expansion factor':>18s}")
    print("  " + "-" * 105)
    for label, M, R_i in bh_cases:
        r = horizon_growth_during_collapse(M, R_i)
        bh_results.append((label, r))
        print(f"  {label:<55s}  {r['Rs_m']:>12.2e}  "
              f"{r['tau_freefall_s']:>12.2e}  "
              f"{r['expansion_factor']:>18.3e}")
    print()

    # Big Bang expansion at various epochs
    print(">>> Big Bang scale factor expansion at various epochs")
    print()
    bb_eras = [
        ("Inflation (typical 60 e-folds)", 1e-32, 1e-30, "exp_inf"),
        ("Radiation era (10^-3 s to 10^12 s)", 1e-3, 1e12, "t^1/2"),
        ("Matter era (10^12 s to 10^17 s)", 1e12, 1e17, "t^2/3"),
        ("Dark energy era (10^17 s to far future)", 1e17, 1e25, "exp_de"),
    ]

    print(f"  {'Era':<42s}  {'t_start (s)':>12s}  {'t_end (s)':>12s}  "
          f"{'a_end/a_start':>15s}  {'form':>8s}")
    print("  " + "-" * 110)
    for label, t_s, t_e, form in bb_eras:
        if form == "exp_inf":
            ratio = np.exp(60.0)  # 60 e-folds, cap for numerical sanity
        elif form == "exp_de":
            ratio = np.exp(min(H_0_KMS * (t_e - t_s), 700.0))
        elif form == "t^1/2":
            ratio = (t_e / t_s) ** 0.5
        elif form == "t^2/3":
            ratio = (t_e / t_s) ** (2/3)
        print(f"  {label:<42s}  {t_s:>12.2e}  {t_e:>12.2e}  "
              f"{ratio:>15.3e}  {form:>8s}")
    print()

    # Mathematical comparison
    print(">>> Mathematical comparison of growth functional forms")
    print()
    print("  BH horizon (toy bubble formation):")
    print("    R_horizon(t) = Rs × (1 - exp(-3 t / tau_ff))")
    print("    Saturates at Rs in time ~ tau_ff")
    print("    Total fractional growth: Rs / R_initial -> small (Rs << R_progenitor)")
    print()
    print("  Big Bang scale factor:")
    print("    Radiation: a(t) ~ t^(1/2)        [decelerating, a_dotdot < 0]")
    print("    Matter:    a(t) ~ t^(2/3)        [decelerating, a_dotdot < 0]")
    print("    Dark energy: a(t) ~ exp(H_0 t)   [accelerating, a_dotdot > 0]")
    print("    Inflation: a(t) ~ exp(H_inf t)  [accelerating, a_dotdot > 0]")
    print()
    print("  Key differences:")
    print("    1. BH horizon SATURATES (asymptotic to Rs); BB scale factor")
    print("       grows monotonically (no upper bound until DE-driven runaway)")
    print("    2. BH formation time ~ tau_ff (~ minutes for stellar collapse)")
    print("       BB expansion happens over the age of the universe")
    print("    3. BH horizon growth is driven by accretion of matter onto")
    print("       a fixed-location boundary; BB expansion is the fabric of")
    print("       space stretching between fixed-content matter")
    print("    4. The 'inside-out' picture: BH formation looks like a Big")
    print("       Bang only if you redefine your time direction (white hole")
    print("       analog) — for an observer outside a forming BH, it's just")
    print("       horizon emergence")
    print()
    print("  Possible mathematical connection:")
    print("    Free-fall collapse to a BH with mass M produces a horizon at")
    print("    Rs in time tau_ff ~ sqrt(R^3/GM). For R = Rs (final state):")
    print("    tau_collapse ~ sqrt(Rs^3/GM) = sqrt(8 G^2 M^2 / c^6) = 2GM/c^3")
    print("    This is the 'natural' BH timescale (light-crossing time of Rs).")
    print()
    print("    For analogy with Big Bang:")
    print("    Time for visible universe to expand by factor 2: 1/H ~ Hubble time")
    print("    Both timescales depend on the relevant mass-energy content.")

    # plot
    plot_path = plot_comparison(bh_results)
    print()
    print(f"Plot: {plot_path}")
    summary_path = write_markdown(bh_results, [plot_path])
    print(f"Summary: {summary_path}")


def plot_comparison(bh_results) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Left: BH horizon growth (normalized)
    ax = axes[0]
    for label, r in bh_results:
        t_norm = r["t_grid_s"] / r["tau_freefall_s"]
        R_norm = r["R_horizon_exp_m"] / r["Rs_m"]
        ax.plot(t_norm, R_norm, linewidth=2,
                label=f"{label[:40]}")
    ax.set_xlabel("t / τ_freefall")
    ax.set_ylabel("R_horizon / R_s")
    ax.set_title("BH horizon emergence during stellar collapse\n"
                 "Saturates at Rs in time ~ τ_ff (toy bubble model)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)

    # Right: Big Bang a(t) on log-log
    ax = axes[1]
    t_BB = np.logspace(-30, 25, 500)
    a_BB = np.array([scale_factor_combined(t) for t in t_BB])
    ax.loglog(t_BB, a_BB, "purple", linewidth=2, label="a(t) — combined")

    # Annotate eras
    ax.axvline(1e-32, color="orange", linestyle=":", alpha=0.6,
               label="Inflation epoch")
    ax.axvline(1.4e12, color="blue", linestyle=":", alpha=0.6,
               label="Matter-radiation equality (~50,000 yr)")
    ax.axvline(4.35e17, color="green", linestyle=":", alpha=0.6,
               label="Today (~13.8 Gyr)")

    ax.set_xlabel("Time since BB (s)")
    ax.set_ylabel("Scale factor a(t)  (a_today = 1)")
    ax.set_title("Big Bang expansion: a(t) over cosmic history\n"
                 "Inflation (exp), radiation (t^½), matter (t^⅔), DE (exp)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    out = PLOTS / "G14_bh_bigbang_expansion.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_markdown(bh_results, plots) -> Path:
    md = []
    md.append("# G14: BH Horizon Growth vs Big Bang Expansion (Wild Theory)\n")
    md.append("## Setup\n")
    md.append(
        "User's wild theory: a stellar collapse to singularity 'resets' "
        "matter and inverts space, with the horizon expansion as something "
        "that might be analogous to the Big Bang. Connects to mainstream "
        "ideas (Smolin's CNS, Penrose's CCC, Big Bounce, white-hole "
        "cosmologies) without committing to any of them.\n"
        "\n"
        "STAM commits to: matter accumulates on the A=1 bubble surface, "
        "no interior. The bubble grows from r=0 to r=Rs as matter joins.\n"
    )

    md.append("## Numerical comparison\n")
    md.append("**BH horizon growth during stellar collapse:**\n")
    md.append("```text")
    md.append(f"{'Case':<55s}  {'Rs (m)':>12s}  {'tau_ff (s)':>12s}")
    md.append("-" * 90)
    for label, r in bh_results:
        md.append(f"{label:<55s}  {r['Rs_m']:>12.2e}  "
                  f"{r['tau_freefall_s']:>12.2e}")
    md.append("```\n")

    md.append("**Big Bang expansion eras:**\n")
    md.append("```text\n"
              "Era              Functional form    Expansion factor\n"
              "----             ---------------    ----------------\n"
              "Inflation        a ~ exp(H_inf t)   ~10^26 (60 e-folds)\n"
              "Radiation        a ~ t^(1/2)        ~10^4 (BBN to recomb)\n"
              "Matter           a ~ t^(2/3)        ~10^3 (recomb to today)\n"
              "Dark energy      a ~ exp(H_0 t)     1+ (today onward)\n"
              "```\n")

    md.append("## Mathematical similarities and differences\n")
    md.append(
        "**Similarities:**\n"
        "- Both involve a 'horizon-like' boundary that grows or recedes.\n"
        "- Both are characterized by a single mass-energy scale (BH mass M; "
        "universe content ρ_total).\n"
        "- Both have a characteristic timescale set by the mass-energy "
        "(τ_BH ~ 2GM/c³; τ_universe ~ 1/H).\n"
        "- For BH formation in time-reverse (white hole) and Big Bang "
        "looked at from outside, the mathematics partially aligns.\n"
        "\n"
        "**Differences:**\n"
        "- BH horizon SATURATES at Rs; Big Bang scale factor has no upper "
        "saturation (continues growing).\n"
        "- BH formation timescale: minutes to hours for stellar collapse; "
        "BB cosmology unfolds over ~13.8 billion years.\n"
        "- BH horizon expansion is driven by mass *accreting onto* a fixed "
        "location; BB expansion is the *spatial fabric stretching* between "
        "fixed-content matter.\n"
        "- The 'inside-out' picture only works if you redefine the time "
        "direction (white-hole analog).\n"
        "- BH horizon growth is sub-luminal (matter must fall in); BB "
        "expansion can superluminally separate distant points (no info "
        "transfer, but space itself stretches).\n"
    )

    md.append("## STAM's take\n")
    md.append(
        "STAM Model-A doesn't commit to BH-as-baby-universe. The framework's "
        "specific claim is that A=1 is the manifold edge — there is no "
        "'inside' to be a separate universe. Matter accumulates on the "
        "boundary, full stop.\n"
        "\n"
        "However, the framework is *compatible* with the wild-theory "
        "extension if one wanted to add it:\n"
        "- The bubble's outer surface (where matter accumulates from "
        "outside) could be one hologram\n"
        "- A complementary 'inner' surface (logically distinct, but "
        "ontologically denied in current STAM) would be the baby-universe "
        "interior\n"
        "- The two-hologram parking-lot idea (your 2-side-of-horizon "
        "intuition) sits on this edge.\n"
        "\n"
        "Mathematically, STAM doesn't predict any specific BH-to-BB "
        "mapping. The two processes have different drivers (matter "
        "accretion vs. cosmic expansion) and different functional forms "
        "(bubble saturation vs. monotonic scale-factor growth).\n"
    )

    md.append("## Mainstream connections\n")
    md.append(
        "**Smolin's Cosmological Natural Selection (1992):** every BH "
        "spawns a baby universe with slightly different physical constants. "
        "This is the closest mainstream proposal to the wild-theory "
        "intuition. Predicts that the constants we observe should be "
        "near-locally-optimized for BH production. Untested but interesting.\n"
        "\n"
        "**Penrose's Conformal Cyclic Cosmology (CCC):** at the universe's "
        "max extent (when only photons remain after all BHs evaporate), "
        "the conformal structure resets and becomes the next Big Bang. "
        "Specific predictions for CMB anisotropy patterns; modest empirical "
        "support, contested.\n"
        "\n"
        "**Big Bounce / loop quantum cosmology:** quantum gravity at the "
        "BB singularity replaces it with a non-singular bounce — the "
        "universe before BB was a contracting cosmos (which could be "
        "the inside of a BH from a parent universe).\n"
        "\n"
        "**White hole cosmologies:** the Big Bang as the time-reverse of "
        "BH collapse. Various proposals (Markov, Frolov, etc.).\n"
    )

    md.append("## Honest framing\n")
    md.append(
        "The wild-theory direction has interesting mainstream company and "
        "STAM is structurally compatible with it (the two-hologram idea "
        "sits naturally there). But STAM doesn't currently make specific "
        "BB-from-BH predictions, and the mathematical analogies break "
        "down quickly under detailed comparison.\n"
        "\n"
        "Worth keeping in the parking-lot file — alongside the multiverse "
        "and dark-energy-as-outer-hologram intuitions — but not in the "
        "active research line. The current framework's distinguishing "
        "predictions are at strong-field (G1), quantum (F6), and "
        "cosmological (G6/G11/G13) scales, not at BB-from-BH scales.\n"
    )

    md.append("## Generated plots\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G14_bh_bigbang_expansion_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


if __name__ == "__main__":
    main()
