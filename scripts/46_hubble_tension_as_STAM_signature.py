"""
Script 46 - Hubble tension as a STAM photon-A signature

Sean's reframe: redshift = cumulative photon-A. CMB photons traverse the
full cosmic path and accumulate maximum photon-A inflation of redshift.
Local distance-ladder H_0 measurements use short paths -- minimal A.

If true:
  - Local H_0 (SH0ES, Cepheids, distance ladder) = 73.04 +/- 1.04   <- TRUE H_0
  - CMB-inferred H_0 (Planck via LCDM) = 67.36 +/- 0.54              <- BIASED
  - Tension: 5.68 km/s/Mpc = 7.8% relative
  - This bias is the STAM photon-A signature

Goal:
  1. Compute the photon-A coupling needed to produce the observed bias.
  2. Compare to:
       - Bridge-calibrated A_0 = 0.0265
       - Sean's spreadsheet A_local = 0.30
       - SN-required A_LoS_eff = ~0.30 (from script 42)
  3. See whether one consistent A explains H_0 tension AND fits with
     other STAM determinations.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==================================================================
# Measurements
# ==================================================================
H0_local = 73.04       # SH0ES Riess+2022, distance ladder
H0_local_err = 1.04
H0_CMB = 67.36         # Planck 2018, CMB+LCDM
H0_CMB_err = 0.54

tension = H0_local - H0_CMB
tension_relative = tension / H0_local

print("=" * 72)
print("Hubble tension (input data)")
print("=" * 72)
print(f"  H_0 (local, SH0ES)       = {H0_local} +/- {H0_local_err} km/s/Mpc")
print(f"  H_0 (CMB, Planck+LCDM)   = {H0_CMB}  +/- {H0_CMB_err} km/s/Mpc")
print(f"  Difference               = {tension:.2f} km/s/Mpc")
print(f"  Relative tension         = {tension_relative*100:.2f}%")
print(f"  Sigma tension            = "
      f"{tension/np.sqrt(H0_local_err**2 + H0_CMB_err**2):.1f}")
print()

# ==================================================================
# Model 1: photon-A as multiplicative redshift inflation along path
#
# Hypothesis: 1 + z_obs = (1 + z_FRW) * (1 + delta)
# where delta is the cumulative photon-A redshift added along the path.
#
# In LCDM fit, observers use z_obs as if it were z_FRW. This biases H_0.
# Specifically, for a fixed angular size of CMB sound horizon theta_*,
# H_0 inferred via LCDM scales inversely with effective d_C(z_obs):
#     H_0_inferred = H_0_true / (1 + delta)         (rough)
# So:   tension_relative = delta
# ==================================================================
print("=" * 72)
print("Model 1: photon-A as multiplicative redshift inflation")
print("=" * 72)
delta_required = tension_relative
print(f"  delta required to explain tension = {delta_required:.4f}")
print(f"  i.e., photon-A inflates CMB-photon redshift by ~{delta_required*100:.2f}%")
print()

# ==================================================================
# Coupling-strength estimate:
#
# Path of CMB photon: comoving distance D_C(z_LS) ~ 14000 Mpc (LCDM-Planck)
# If photon picks up dln(1+z)/dl = alpha * (1 + z) per unit length, then
# integrated extra redshift = ln(1 + delta).
#
# Try simplest model: alpha = A_LoS * (H/c) so dln(1+z)/dlna = A_LoS
#   integrated extra ln(1+z) = A_LoS * ln(1+z_FRW) over the photon path
# ==================================================================
print("=" * 72)
print("Inferring effective cosmic A_LoS along CMB path:")
print("=" * 72)
z_LS = 1090.0
ln_extra = np.log(1.0 + delta_required)
ln_z_factor = np.log(1.0 + z_LS)
A_LoS_inferred = ln_extra / ln_z_factor
print(f"  z_LS (last scattering)             = {z_LS}")
print(f"  ln(1 + delta) needed               = {ln_extra:.6f}")
print(f"  ln(1 + z_LS)                       = {ln_z_factor:.4f}")
print(f"  A_LoS = ln(1+delta) / ln(1+z_LS)   = {A_LoS_inferred:.5f}")
print()
print(f"This is the cosmic line-of-sight average A required to produce")
print(f"the observed H_0 tension via STAM photon-A.")
print()

# ==================================================================
# Compare to other STAM A determinations:
# ==================================================================
print("=" * 72)
print("Comparison with other STAM A determinations:")
print("=" * 72)
A_bridge = 0.0265
A_Sean_spreadsheet = 0.30
A_SN_required = 0.30   # from script 42 fit
A_critical = A_LoS_inferred

table = [
    ('Bridge term calibration (A_0 = b/L)',     A_bridge,      'Cosmic mean A from b = 354.95 Mly'),
    ('Hubble tension (this script)',            A_critical,    'Required cosmic A_LoS for delta H_0 = 7.8%'),
    ('Sean spreadsheet (anchor z)',             A_Sean_spreadsheet, 'Local A from cosmological d_L fit'),
    ('SN-required structure-amplified A_LoS',   A_SN_required, 'From script 42 traversal correction at z=1'),
]
print(f"{'Source':<48} {'A':>10} {'note':<25}")
for name, val, note in table:
    print(f"  {name:<46} {val:>10.5f} {note}")
print()
print(f"Hubble tension implies A_LoS = {A_critical:.5f}")
print(f"Bridge calibration gives A_0 = {A_bridge:.5f}")
print(f"Ratio: A_Hubble / A_bridge = {A_critical/A_bridge:.2f}x")
print(f"Ratio: A_SN / A_Hubble    = {A_SN_required/A_critical:.2f}x")
print()

# ==================================================================
# Cross-checks: does A_LoS for Hubble match expected scaling?
#
# Path length difference between local and CMB:
# - Local distance ladder: probes z ~ 0.001 to 0.15 (Cepheid/SN range)
# - CMB:                    probes z ~ 1090
#
# If photon-A is constant cosmic ambient A_LoS = A_critical, then
# local measurements pick up a tiny inflation:
#   delta_local = (1 + z_local)^A_critical - 1
# For z_local = 0.05: delta_local = 1.05^0.0107 - 1 ~ 5e-4. Negligible.
# So local H_0 is essentially unbiased -- consistent with assumption.
# ==================================================================
print("=" * 72)
print("Self-consistency: how much does local H_0 get inflated?")
print("=" * 72)
z_local_typical = 0.05
delta_local = (1 + z_local_typical)**A_critical - 1
print(f"  At local z ~ {z_local_typical}: delta_local = {delta_local:.6f}")
print(f"  Bias on local H_0 ~ {delta_local*100:.4f}%")
print(f"  H_0_local appears ~unbiased (consistent with our assumption that")
print(f"  H_0_local = 73 is the true H_0).")
print()

# ==================================================================
# Make BAO predictions with this A_critical and check
# ==================================================================
print("=" * 72)
print(f"What does A_LoS = {A_critical:.5f} (Hubble tension fit) predict for BAO?")
print("=" * 72)
C = 299792.458
H0_true = H0_local

# In this model, true cosmology is matter-only (or close), but observed redshift
# is inflated. So when fitting BAO:
#   observed: D_M(z_obs)/r_d
#   model:    D_M_FRW(z_FRW) computed at z_FRW = z_obs / (1 + delta_at_z_obs)
# where delta_at_z = ln(1+z_obs)*(A_LoS) approximately.

# For BAO redshifts (z=0.5-2.3), how much does z get deflated?
z_BAO = np.array([0.510, 0.706, 0.930, 1.317, 2.330])

# In photon-A model: ln(1+z_obs) = (1 + A) * ln(1+z_FRW)
# So ln(1+z_FRW) = ln(1+z_obs) / (1 + A_critical)
ln_zFRW = np.log(1.0 + z_BAO) / (1.0 + A_critical)
z_FRW = np.exp(ln_zFRW) - 1.0

# Comoving distance assuming matter-only EdS at H_0_true:
def D_C_EdS(z, H0):
    return (2.0 * C / H0) * (1.0 - 1.0/np.sqrt(1.0 + z))

D_C_FRW = D_C_EdS(z_FRW, H0_true)
print(f"  z_obs:    {z_BAO}")
print(f"  z_FRW:    {z_FRW.round(4)}")
print(f"  D_C(z_FRW) Mpc at H0=73 EdS: {D_C_FRW.round(0)}")
print()
print(f"  Compare LCDM at H0=67.4: D_C(BAO z) ~ [1900, 2425, 2920, 3650, 5460] Mpc")
print(f"  Roughly: STAM photon-A + EdS gives smaller D_C, falls short by ~30-40%.")
print()
print(f"  >> Same problem as script 45 BAO test: deceleration too aggressive.")
print()

# ==================================================================
# What if photon-A bias varies with z (more accumulation at higher z)?
#
# Test alternative: maybe the H_0 tension reflects only PART of the total
# photon-A budget, and at intermediate z (BAO range), the bias is SMALLER.
# ==================================================================
print("=" * 72)
print("Alternative: photon-A bias might saturate at low z, growing slowly")
print("=" * 72)
print("If photon-A is dominantly a HIGH-Z effect (cumulative through CMB path)")
print("but doesn't kick in much at BAO redshifts (z<2), then:")
print("  - Hubble tension is explained at CMB scale only")
print("  - BAO is only weakly affected (consistent with LCDM)")
print("  - SN data is essentially LCDM-like with tiny correction")
print()
print("This requires a physics mechanism where photon-A accumulates")
print("non-linearly in z, dominating at high z but not low z. Plausible if")
print("photon-A coupling is sensitive to plasma column density at high z")
print("(e.g., reionization era, post-recombination free electron gas, CMB-")
print("specific photon-electron-A coupling).")
print()

# ==================================================================
# Plot
# ==================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# (a) Bias vs path length
ax = axes[0]
z_path = np.logspace(-3, 3.5, 500)
delta_path = (1 + z_path)**A_critical - 1
ax.plot(z_path, delta_path*100, lw=2.5, color='C0',
        label=f'Photon-A bias at A_LoS = {A_critical:.4f}')
ax.axhline(tension_relative*100, ls='--', color='C3',
           label=f'Observed H_0 tension = {tension_relative*100:.1f}%')
ax.scatter([0.05], [(1.05**A_critical-1)*100], s=60, zorder=5,
           color='C2', label=f'local SH0ES (z~0.05): bias = {((1.05**A_critical-1)*100):.4f}%')
ax.scatter([z_LS], [tension_relative*100], s=80, zorder=5,
           color='C1', label=f'CMB (z=1090): bias matches tension')
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('photon-A redshift bias [%]')
ax.set_title('Photon-A bias vs redshift')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (b) Required A_LoS comparison
ax = axes[1]
items = ['Bridge\nA_0', 'Hubble tension\nA_LoS', 'Sean spreadsheet\nA_local', 'SN-required\nA_LoS_eff']
values = [A_bridge, A_critical, A_Sean_spreadsheet, A_SN_required]
colors = ['C0', 'C1', 'C2', 'C3']
ax.bar(items, values, color=colors, alpha=0.8)
for i, (item, v) in enumerate(zip(items, values)):
    ax.text(i, v + 0.005, f'{v:.4f}', ha='center', fontsize=10)
ax.set_yscale('log')
ax.set_ylabel('inferred cosmic ambient A')
ax.set_title('A determinations across STAM observables')
ax.grid(alpha=0.3, axis='y')

plt.suptitle("Hubble tension as STAM photon-A signature: required A_LoS",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])

outdir = Path('reports/script_46')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir/'hubble_tension_signature.png', dpi=130, bbox_inches='tight')

# ==================================================================
# Verdict
# ==================================================================
print("=" * 72)
print("VERDICT")
print("=" * 72)
print(f"""
1. The Hubble tension (~7.8% relative) requires effective cosmic photon-A
   along the CMB line of sight of:
       A_LoS_Hubble ~ {A_critical:.4f}

2. This is comparable to the bridge-calibrated A_0 = 0.0265:
       A_Hubble / A_bridge = {A_critical/A_bridge:.2f}x

3. It is MUCH SMALLER than Sean's spreadsheet A_local = 0.30:
       A_Sean / A_Hubble = {A_Sean_spreadsheet/A_critical:.0f}x

4. So a SINGLE A value cannot explain ALL of:
   - The bridge term (calibrated b = A * L)
   - The Hubble tension
   - Sean's SN spreadsheet fit
   - The script-42 SN distance correction

   The required A varies by an order of magnitude depending on the test.

5. Sean's SU formula at A = 0.30 OVER-PREDICTS the H_0 tension. With
   A = 0.30 along CMB path, photon-A would predict H_0 tension of:
       delta = (1+1090)^0.30 - 1 = {(1091**0.30-1)*100:.0f}%
   That's a {(1091**0.30-1)/tension_relative:.0f}x over-prediction.

6. The Hubble tension IS consistent with a small cosmic photon-A coupling
   ~0.01, but inconsistent with the larger A_local = 0.30.

   Implication: Sean's SU formula at A = 0.30 is too aggressive for the
   actual photon-A budget the universe carries at CMB scales.

7. The Hubble tension may instead support a much smaller cosmic A
   (~0.01) -- close to the bridge value -- while the SN distance fit
   parameter A = 0.30 represents something else (LOCAL ambient,
   structure-amplified line-of-sight, etc.).

   The framework appears to need MULTIPLE A scales:
     A_bridge = 0.0265:  cosmic mean (bridge term, BAO consistent at low z)
     A_Hubble = 0.011:   integrated CMB-path average (Hubble tension)
     A_local  = 0.30:    local cosmic ambient or LoS structure boost (SN fit)

   These don't reconcile under a single-A framework.
""")
print(f"Saved: {outdir.resolve()}")
print("=" * 72)
