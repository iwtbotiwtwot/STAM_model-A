#!/usr/bin/env python3
"""
F5_gw_polarization.py

FATALITY 5: Bold-STAM gravitational wave polarization structure.
Author: Sean Brady / STAM Model-A continuation

The opponent's attack:
    "GR predicts gravitational waves with two transverse-tensor polarizations
     (h_+, h_×). LIGO/Virgo/KAGRA observations have confirmed these and
     constrained any non-tensor polarizations (scalar / vector / longitudinal)
     to a few percent of the tensor amplitude. Pure scalar theories of
     gravity (Nordstrom 1913 and successors) predict only a longitudinal
     scalar 'breathing mode' and have been ruled out historically.
     If bold STAM is scalar gravity (Framework C), bold STAM is also ruled
     out by LIGO."

What this script computes:
    Take the bold-STAM Framework C construction rule
        g_tt = -(1-A) c²
        g_rr = 1/[(1-A)(1-A²)²]
    and apply it to a small time-dependent A perturbation around flat space.
    Determine the polarization structure of the resulting metric perturbation
    h_μν, compare to GR's TT modes (h_+, h_×), and assess against LIGO.

Honest expectation (verified below):
    The simplest scalar generalization of the construction rule predicts
    a single longitudinal/scalar 'breathing' mode. No transverse-tensor modes
    (h_+, h_×). This is the Nordstrom-type scalar-gravity polarization, which
    LIGO has ruled out.

    **Bold STAM in Framework C as specified is FALSIFIED by LIGO GW
    observations** unless the construction rule is amended to a scalar-tensor
    or higher-rank structure that produces tensor modes alongside the scalar
    A field.

This is the first fatality bold STAM does not survive. The path forward is
clear (commit to a richer field-equation structure), but it represents real
theoretical work — not a quick fix.

What this script does NOT do:
    - Compute LIGO event waveforms in detail.
    - Distinguish between possible scalar-tensor amendments to bold STAM.
    - Try to rescue bold STAM via tweaks to the construction rule.
    Each of those would be a separate piece of work. This script just
    diagnoses the problem cleanly.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = 2.99792458e8

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Linearized perturbation in bold-STAM Framework C ---

def metric_perturbation_bold_stam(delta_A: float | np.ndarray) -> dict:
    """Compute h_μν components for a small δA perturbation around flat space.

    Bold STAM construction rule:
        g_tt = -(1-A) c²
        g_rr = 1/[(1-A)(1-A²)²]
        g_θθ = r²,  g_φφ = r² sin²θ   (UNCHANGED — angular components are not
                                         functions of A in the simplest
                                         spherical-symmetric construction)

    Linearizing around A = 0:
        δg_tt = +c² δA               (clock runs slower)
        δg_rr = +δA + O(δA²)         (radial proper distance stretches)
        δg_θθ = 0                    (angular metric unchanged at this order)
        δg_φφ = 0

    For a GW propagating in z direction, the relevant components in
    Cartesian coordinates (for an observer at flat space far from source) are:
        h_tt = +c² δA   (time-time perturbation; LIGO calls this scalar mode)
        h_zz ~ +δA      (longitudinal; in direction of propagation)
        h_xx, h_yy = 0  (transverse in-plane components — NO TT modes!)
        h_xy = 0
    """
    return {
        "h_tt_over_c2": 1.0 * np.asarray(delta_A),
        "h_zz_longitudinal": 1.0 * np.asarray(delta_A),
        "h_xx_transverse": 0.0 * np.asarray(delta_A),
        "h_yy_transverse": 0.0 * np.asarray(delta_A),
        "h_xy_cross": 0.0 * np.asarray(delta_A),
        "h_plus": 0.0 * np.asarray(delta_A),  # = (h_xx - h_yy) / 2
        "h_cross": 0.0 * np.asarray(delta_A),  # = h_xy
        "trace": 1.0 * np.asarray(delta_A),    # h_xx + h_yy + h_zz - h_tt/c² (signed by gauge)
    }


def metric_perturbation_GR(h_plus: float, h_cross: float) -> dict:
    """For comparison: GR TT-gauge GW components.

    In TT gauge for vacuum GW propagating in z direction:
        h_tt = 0
        h_zz = 0
        h_xx = +h_plus
        h_yy = -h_plus
        h_xy = +h_cross
        Trace = h_xx + h_yy = 0  (TT gauge)
    """
    return {
        "h_tt_over_c2": 0.0,
        "h_zz_longitudinal": 0.0,
        "h_xx_transverse": +h_plus,
        "h_yy_transverse": -h_plus,
        "h_xy_cross": h_cross,
        "h_plus": h_plus,
        "h_cross": h_cross,
        "trace": 0.0,
    }


# --- LIGO constraints (representative published values) ---

# LIGO/Virgo bounds on non-tensor polarization fractions (representative;
# see e.g. arXiv:1809.04817 and follow-ups).
# These are the upper limits on the fraction of total GW power that could
# be in non-tensor modes for typical compact-binary mergers.
LIGO_SCALAR_FRACTION_UPPER_BOUND = 0.10   # ~10% upper bound for representative searches
LIGO_VECTOR_FRACTION_UPPER_BOUND = 0.10
LIGO_TENSOR_FRACTION_LOWER_BOUND = 0.50   # at least half must be tensor for these searches


def assess_bold_stam_against_ligo() -> dict:
    """Compute polarization fractions for bold-STAM-Framework-C GW prediction.

    In bold STAM with simplest generalization:
        Tensor (h_+, h_×) fraction: 0
        Scalar / longitudinal fraction: 100%
    Compare to LIGO bounds.
    """
    bold_pert = metric_perturbation_bold_stam(1.0)  # arbitrary normalization

    tensor_power = bold_pert["h_plus"] ** 2 + bold_pert["h_cross"] ** 2
    scalar_long_power = (bold_pert["h_tt_over_c2"] ** 2 +
                          bold_pert["h_zz_longitudinal"] ** 2)

    total_power = tensor_power + scalar_long_power
    if total_power == 0:
        return {"verdict": "no GW prediction"}
    tensor_frac = float(tensor_power / total_power)
    scalar_frac = float(scalar_long_power / total_power)

    return {
        "tensor_fraction": tensor_frac,
        "scalar_fraction": scalar_frac,
        "tensor_passes_ligo": tensor_frac >= LIGO_TENSOR_FRACTION_LOWER_BOUND,
        "scalar_passes_ligo": scalar_frac <= LIGO_SCALAR_FRACTION_UPPER_BOUND,
    }


# --- Visualizations ---

def plot_polarization_patterns() -> Path:
    """Visualize how a ring of test particles deforms under different polarizations."""
    theta = np.linspace(0, 2 * np.pi, 100)
    x_ring = np.cos(theta)
    y_ring = np.sin(theta)

    h_amplitude = 0.3
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))

    # GR h_plus mode
    x_plus = x_ring * (1 + h_amplitude)
    y_plus = y_ring * (1 - h_amplitude)
    axes[0].plot(x_plus, y_plus, color="tab:blue", linewidth=2)
    axes[0].plot(x_ring, y_ring, color="gray", linewidth=1, alpha=0.4, linestyle="--")
    axes[0].set_aspect("equal")
    axes[0].set_xlim(-1.5, 1.5)
    axes[0].set_ylim(-1.5, 1.5)
    axes[0].set_title("GR: h_+ mode\n(transverse stretch/compress)")
    axes[0].grid(True, alpha=0.3)

    # GR h_cross mode
    rotation = np.pi / 4
    cos_r = np.cos(rotation)
    sin_r = np.sin(rotation)
    x_rot = cos_r * x_ring - sin_r * y_ring
    y_rot = sin_r * x_ring + cos_r * y_ring
    x_cross = x_rot * (1 + h_amplitude)
    y_cross = y_rot * (1 - h_amplitude)
    x_back = cos_r * x_cross + sin_r * y_cross
    y_back = -sin_r * x_cross + cos_r * y_cross
    axes[1].plot(x_back, y_back, color="tab:blue", linewidth=2)
    axes[1].plot(x_ring, y_ring, color="gray", linewidth=1, alpha=0.4, linestyle="--")
    axes[1].set_aspect("equal")
    axes[1].set_xlim(-1.5, 1.5)
    axes[1].set_ylim(-1.5, 1.5)
    axes[1].set_title("GR: h_× mode\n(transverse, rotated 45°)")
    axes[1].grid(True, alpha=0.3)

    # Bold STAM scalar/breathing mode (uniform expansion of ring)
    x_breath = x_ring * (1 + h_amplitude)
    y_breath = y_ring * (1 + h_amplitude)
    axes[2].plot(x_breath, y_breath, color="tab:red", linewidth=2)
    axes[2].plot(x_ring, y_ring, color="gray", linewidth=1, alpha=0.4, linestyle="--")
    axes[2].set_aspect("equal")
    axes[2].set_xlim(-1.5, 1.5)
    axes[2].set_ylim(-1.5, 1.5)
    axes[2].set_title("Bold STAM (Framework C): scalar / breathing\n(uniform expansion / contraction)")
    axes[2].grid(True, alpha=0.3)

    # Bold STAM longitudinal mode (in z direction; viewed from side, ring becomes ellipse)
    # For a wave propagating in z, longitudinal mode would compress/stretch in z
    # Show in a side view
    z_axis = np.cos(theta) * 1
    r_axis = np.sin(theta) * 1
    z_compressed = z_axis * (1 + h_amplitude)
    r_compressed = r_axis * 1
    axes[3].plot(z_compressed, r_compressed, color="tab:orange", linewidth=2)
    axes[3].plot(z_axis, r_axis, color="gray", linewidth=1, alpha=0.4, linestyle="--")
    axes[3].set_aspect("equal")
    axes[3].set_xlim(-1.5, 1.5)
    axes[3].set_ylim(-1.5, 1.5)
    axes[3].set_xlabel("z (propagation direction)")
    axes[3].set_ylabel("r (transverse)")
    axes[3].set_title("Bold STAM (Framework C): longitudinal\n(stretch in propagation direction)")
    axes[3].grid(True, alpha=0.3)

    fig.suptitle("Gravitational wave polarization patterns — GR (h_+, h_×) vs Bold STAM (scalar/longitudinal)",
                 y=1.02, fontsize=13)

    out = PLOTS / "F5_polarization_patterns.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    return out


def plot_polarization_fractions() -> Path:
    """Show predicted polarization fractions and LIGO constraint bounds."""
    fig, ax = plt.subplots(figsize=(10, 5.5))

    bold_assessment = assess_bold_stam_against_ligo()

    categories = ["Tensor\n(h_+, h_×)", "Scalar /\nLongitudinal"]
    bold_values = [bold_assessment["tensor_fraction"], bold_assessment["scalar_fraction"]]
    gr_values = [1.0, 0.0]  # GR is pure tensor

    x = np.arange(len(categories))
    width = 0.3

    bars1 = ax.bar(x - width, gr_values, width, label="GR prediction", color="tab:blue", alpha=0.8)
    bars2 = ax.bar(x, bold_values, width, label="Bold STAM (Framework C, simplest gen.)",
                    color="tab:red", alpha=0.8)
    # LIGO bounds
    ax.bar(x[0] + width, LIGO_TENSOR_FRACTION_LOWER_BOUND, width, label="LIGO lower bound on tensor",
           color="tab:green", alpha=0.5)
    ax.bar(x[1] + width, LIGO_SCALAR_FRACTION_UPPER_BOUND, width, label="LIGO upper bound on scalar",
           color="tab:orange", alpha=0.5)

    for bar, val in zip(bars1, gr_values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02, f"{val:.2f}",
                ha="center", fontsize=10)
    for bar, val in zip(bars2, bold_values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02, f"{val:.2f}",
                ha="center", fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Fraction of GW power")
    ax.set_title("GW polarization fractions: GR vs Bold STAM (Framework C) vs LIGO bounds")
    ax.set_ylim(0, 1.15)
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    out = PLOTS / "F5_polarization_fractions.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(assessment: dict, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A F5: Gravitational Wave Polarization Structure\n")
    md.append("## The attack\n")
    md.append(
        "GR predicts gravitational waves with two transverse-tensor polarizations: `h_+` "
        "and `h_×`. LIGO/Virgo/KAGRA observations have confirmed these tensor modes and "
        "constrained any non-tensor polarizations (scalar / vector / longitudinal) to less "
        "than ~10% of the total GW power for typical compact-binary mergers (representative "
        "value; specific bounds vary by analysis).\n"
        "\n"
        "Pure scalar theories of gravity (Nordstrom 1913 and similar) predict only a "
        "longitudinal scalar 'breathing mode' and were ruled out historically. The opponent "
        "asks: bold STAM in Framework C is scalar gravity (with bold-STAM metric construction "
        "rule). Doesn't this place bold STAM in the same falsified class as Nordstrom?\n"
    )
    md.append("## The honest analysis\n")
    md.append(
        "Linearizing the bold-STAM construction rule around flat space (small `δA` "
        "perturbation):\n"
        "```text\n"
        "g_tt = -(1-A) c²       =>   δg_tt = +c² δA       (time-time component)\n"
        "g_rr = 1/[(1-A)(1-A²)²] =>   δg_rr ≈ +δA          (radial component)\n"
        "g_θθ = r²              =>   δg_θθ = 0            (angular — UNCHANGED)\n"
        "g_φφ = r² sin²θ        =>   δg_φφ = 0\n"
        "```\n"
        "\n"
        "For a GW propagating in the z direction, in Cartesian coordinates at the detector:\n"
        "```text\n"
        "h_tt = +c² δA          (scalar / time-time mode)\n"
        "h_zz = +δA             (longitudinal — in direction of propagation)\n"
        "h_xx = h_yy = 0        (transverse components vanish — NO h_+ MODE!)\n"
        "h_xy = 0               (NO h_× MODE!)\n"
        "```\n"
        "\n"
        "**The TT-gauge tensor modes (h_+, h_×) are exactly zero in bold STAM Framework C.** "
        "Bold STAM, in its simplest generalization to time-dependent perturbations, predicts "
        "only longitudinal scalar GW modes. This is the same problem Nordstrom's gravity had.\n"
    )
    md.append("## Comparison with LIGO\n")
    md.append(
        f"- Bold STAM tensor fraction: {assessment['tensor_fraction']:.2f}  "
        f"(LIGO requires ≥ {LIGO_TENSOR_FRACTION_LOWER_BOUND:.2f}). "
        f"**FAIL** ({assessment['tensor_passes_ligo']}).\n"
        f"- Bold STAM scalar fraction: {assessment['scalar_fraction']:.2f}  "
        f"(LIGO requires ≤ {LIGO_SCALAR_FRACTION_UPPER_BOUND:.2f}). "
        f"**FAIL** ({assessment['scalar_passes_ligo']}).\n"
        "\n"
        "**Verdict on F5: bold STAM in Framework C, with the simplest scalar generalization, "
        "is FALSIFIED by LIGO/Virgo/KAGRA observations of GW polarization.**\n"
    )
    md.append("## Path forward — what bold STAM has to do to survive\n")
    md.append(
        "F5 does not refute bold STAM as a *concept*. It refutes Framework C in its simplest "
        "form. The construction rule must be amended so that the metric has additional "
        "structure beyond what the scalar A provides directly. The natural amendments are:\n"
        "\n"
        "**Option I: Scalar-tensor amendment.** Promote the metric to an independent "
        "dynamical field, coupled to A through the construction rule but with its own "
        "tensor degrees of freedom. Bold STAM becomes a scalar-tensor theory of gravity, "
        "structurally similar to Brans-Dicke or its generalizations. The TT modes (h_+, "
        "h_×) come from the metric's own dynamics; the scalar A provides an additional "
        "mode that LIGO has bounded but not ruled out at the few-percent level.\n"
        "\n"
        "**Option II: A as a higher-rank field.** Promote A from a scalar to a tensor-"
        "valued field, where the construction rule maps a tensor A_μν to the metric. This "
        "is more elaborate and less STAM-native; the original framework treats A as a "
        "scalar.\n"
        "\n"
        "**Option III: Amended construction rule with derivative couplings.** Allow the "
        "metric to depend on derivatives of A (∂A, ∂²A) in addition to A itself. Could "
        "potentially produce TT modes from how A propagates, even with A as a scalar.\n"
        "\n"
        "Option I is the most direct fix. It does require giving up some of bold STAM's "
        "current cleanness — the metric is no longer a derived consequence of A alone. But "
        "it preserves all the local results (Q8-Q12 thermodynamics, SF3-SF5 strong-field "
        "structure, F1 evasion, F3 evasion) because those don't depend on the tensor mode "
        "structure of GWs. The amendment shows up only when GWs are computed.\n"
    )
    md.append("## Honest verdict\n")
    md.append(
        "**F5 is the first fatality bold STAM does not survive in its current form.** "
        "Framework C, as we specified it in SF6, predicts the wrong GW polarization "
        "structure. LIGO observations directly contradict this prediction.\n"
        "\n"
        "However, F5 is *not* a refutation of bold STAM as a research program. It is a "
        "specific failure of the simplest scalar generalization of the bold-STAM "
        "construction rule. The fix (scalar-tensor amendment) is well-understood "
        "theoretically — it just hasn't been written down for bold STAM specifically. "
        "Doing so is real research work, comparable in scope to what was needed to write "
        "down Q8-Q12 or SF3-SF5.\n"
        "\n"
        "The framework's status, after F5:\n"
        "- *Survived cleanly:* F1 (Penrose-Hawking), F2 (equivalence principle).\n"
        "- *Survived after Framework C commitment:* F3 (Birkhoff).\n"
        "- *Failed in current form, with clear path to amend:* F5 (GW polarization).\n"
        "\n"
        "**Bold STAM is unrefuted as a conceptual framework**, but it is not yet a "
        "complete theory of gravity. The completion requires committing to scalar-tensor "
        "structure (Option I above) and rebuilding the perturbation theory accordingly.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "F5_gw_polarization_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    assessment = assess_bold_stam_against_ligo()

    plot_paths = [
        plot_polarization_patterns(),
        plot_polarization_fractions(),
    ]

    summary = write_markdown(assessment, plot_paths)

    print("STAM Model-A F5: GW Polarization Fatality")
    print("=" * 64)
    print()
    print("Linearizing bold-STAM Framework C around flat space (delta A perturbation):")
    print("  h_tt = +c^2 dA   (scalar mode)")
    print("  h_zz = +dA       (longitudinal mode)")
    print("  h_xx = h_yy = 0  (NO transverse h_+ mode)")
    print("  h_xy = 0         (NO transverse h_x mode)")
    print()
    print("Polarization fractions:")
    print(f"  Tensor (h_+, h_x):  {assessment['tensor_fraction']:.2f}  "
          f"(LIGO requires >= {LIGO_TENSOR_FRACTION_LOWER_BOUND:.2f})")
    print(f"  Scalar / longitudinal: {assessment['scalar_fraction']:.2f}  "
          f"(LIGO requires <= {LIGO_SCALAR_FRACTION_UPPER_BOUND:.2f})")
    print()
    print("Verdict: Bold STAM in Framework C (simplest generalization) is FALSIFIED by LIGO.")
    print()
    print("Path forward: amend Framework C to scalar-tensor structure (Option I in summary).")
    print("This is real theoretical work but doesn't invalidate bold STAM as a research program.")
    print()
    print("Files written:")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
