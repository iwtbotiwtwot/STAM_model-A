"""
Script 40 - Where do supernovae REALLY sit in STAM, and why do they LOOK different?

Sean's clean framing (2026-05-07):
  STAM should be able to say
    'This is where SN actually sit' (intrinsic prediction)
    'That is why they look different' (the adjustment, with a physical cause)
  Without smuggling dark energy into V(A) = beta/(1-A) postulate.

Test:
  1. Compute STAM's INTRINSIC prediction = Einstein-de Sitter (matter-only,
     water-tank-only, no V(A) potential, no Lambda). This is what STAM
     predicts from the water-tank principle alone.
  2. Compute LCDM's prediction (the data follows this closely).
  3. Compute the "adjustment needed" Delta_mu_needed(z) = mu_LCDM - mu_EdS.
  4. Compare to candidate physical adjustments at A_0 = 0.0265:
        - Constant Shapiro through ambient A: Delta_mu ~ 2.17 * A_0
        - Linear bridge: b*z with b = 355 Mly
        - Modified-Friedmann shape: 0.105*ln(1+z) (script 38)
  5. Determine: can STAM explain the adjustment with A_0 = 0.0265 photon-A
     mechanisms alone, or does it require a larger effective coupling?
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

C_KMS = 299792.458
H0_KMS_MPC = 67.4
A0 = 0.026514
OMEGA_M = 0.315
OMEGA_R = 9.2e-5

def cumtrapz0(y, x):
    dx = np.diff(x); midy = 0.5*(y[1:]+y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy*dx)])

def mu_from_h2(z_obs, h2_func, n_grid=8000, z_max=None):
    if z_max is None:
        z_max = max(np.max(z_obs)*1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0/(1.0+z_grid)
    inv_h = 1.0/np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS/H0_KMS_MPC)*cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0+z_obs)*DC_obs
    return 5.0*np.log10(dL_obs)+25.0

# Three cosmologies:

# 1. STAM intrinsic = Einstein-de Sitter (matter only, no DE, water-tank only)
#    Flat universe with Omega_m = 1 (this is what straight STAM predicts
#    as a flat matter-dominated cosmology with no built-in cosmological term).
#    For comparison we also compute "STAM intrinsic with current Omega_m=0.315
#    but no DE" which is open universe -- we use Omega_m=1 EdS as the cleanest
#    "no dark energy" reference.
def h2_EdS_flat(a):
    return a**-3   # Omega_m = 1, Omega_DE = 0

# 2. LCDM (data follows this closely)
def h2_LCDM(a):
    return OMEGA_M*a**-3 + OMEGA_R*a**-4 + (1.0-OMEGA_M-OMEGA_R)

# 3. STAM with V(A) = beta/(1-A) tracking (script 36) -- has DE built in
BETA_TILDE = OMEGA_M*(1.0-A0)**2
def A_eq(a): return 1.0 - (1.0-A0)*a**1.5
def Y_pot(a): return BETA_TILDE/(1.0-A_eq(a))
def K_kin(a): return (3.0/8.0)*(1.0-A0)**2*a**3
def h2_STAM_VofA_raw(a):
    return (OMEGA_M*a**-3 + OMEGA_R*a**-4 + Y_pot(a))/(1.0-K_kin(a))
H2_TODAY = h2_STAM_VofA_raw(1.0)
def h2_STAM_VofA(a):
    return h2_STAM_VofA_raw(a)/H2_TODAY

# ==================================================================
# Compute mu(z) for all three on a smooth grid
# ==================================================================
z = np.linspace(0.001, 2.5, 600)
mu_EdS_flat   = mu_from_h2(z, h2_EdS_flat,   z_max=3.0)
mu_LCDM       = mu_from_h2(z, h2_LCDM,       z_max=3.0)
mu_STAM_VofA  = mu_from_h2(z, h2_STAM_VofA,  z_max=3.0)

# The "adjustment needed" to bridge intrinsic-STAM to observed (LCDM):
Delta_mu_needed = mu_LCDM - mu_EdS_flat

# What we previously called the "derived adjustment" (modified-Friedmann gap):
Delta_mu_modFried = mu_LCDM - mu_STAM_VofA   # ~0.105 ln(1+z)

# ==================================================================
# Candidate physical-adjustment forms at A_0 = 0.0265
# ==================================================================
# (a) Constant Shapiro through ambient A_0:
#       Delta_mu_Shapiro = 5 log10(1 + A_0) ~ 2.17 * A_0
const_Shapiro = 5.0*np.log10(1.0 + A0)*np.ones_like(z)

# (b) Linear bridge term b*z with b = 354.95 Mly:
#       Delta_d_L = b z (Mly), Delta_mu = 5 log10(1 + b z / d_L_EdS)
b_Mly = 354.95
c_over_H0_Mly = (C_KMS/H0_KMS_MPC) * 3.262    # Mpc -> Mly factor
d_L_EdS_Mly = 10.0**((mu_EdS_flat - 25.0)/5.0)  # in Mpc
d_L_EdS_Mly = d_L_EdS_Mly * 3.262               # Mpc -> Mly
linear_bridge = 5.0*np.log10(1.0 + b_Mly*z/d_L_EdS_Mly)

# (c) Modified-Friedmann shape (just for reference):
modFried_shape = 0.105*np.log1p(z)

# ==================================================================
# How big is the gap each candidate covers?
# ==================================================================
# At z=1:
i_z1 = np.argmin(np.abs(z - 1.0))
print("=" * 72)
print("Adjustment needed to map STAM intrinsic (water-tank-only) -> observed")
print("=" * 72)
print(f"At z = 0.5: mu_LCDM - mu_EdS = {Delta_mu_needed[np.argmin(np.abs(z-0.5))]:+.3f} mag")
print(f"At z = 1.0: mu_LCDM - mu_EdS = {Delta_mu_needed[i_z1]:+.3f} mag")
print(f"At z = 2.0: mu_LCDM - mu_EdS = {Delta_mu_needed[np.argmin(np.abs(z-2.0))]:+.3f} mag")
print()
print("How much each candidate adjustment can produce at z = 1.0:")
print(f"  Constant Shapiro through A_0=0.0265 :  {const_Shapiro[i_z1]:+.3f} mag")
print(f"  Linear bridge b*z (b = 355 Mly)     :  {linear_bridge[i_z1]:+.3f} mag")
print(f"  Modified-Friedmann shape            :  {modFried_shape[i_z1]:+.3f} mag")
print(f"  Needed                              :  {Delta_mu_needed[i_z1]:+.3f} mag")
print()
print("Coverage ratios at z = 1.0 (candidate / needed):")
print(f"  Shapiro / needed     = {const_Shapiro[i_z1]/Delta_mu_needed[i_z1]:.3f}")
print(f"  Linear bridge / needed = {linear_bridge[i_z1]/Delta_mu_needed[i_z1]:.3f}")
print(f"  ModFried / needed    = {modFried_shape[i_z1]/Delta_mu_needed[i_z1]:.3f}")
print()

# What effective A would be needed if Shapiro alone were the mechanism?
A_eff_needed = (10.0**(Delta_mu_needed[i_z1]/5.0) - 1.0)
print(f"If single-multiplicative Shapiro covered the whole gap at z=1, "
      f"effective A would need to be {A_eff_needed:.3f}")
print(f"(Calibrated A_0 = {A0:.4f}, off by factor "
      f"{A_eff_needed/A0:.1f}x)")

# ==================================================================
# Plotting
# ==================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: mu(z) absolute, all three
ax = axes[0,0]
ax.plot(z, mu_EdS_flat, lw=2.2, color='C2',
        label='STAM intrinsic (matter-only, no V(A))')
ax.plot(z, mu_LCDM, lw=2.2, color='k', label='LCDM (data follows this)')
ax.plot(z, mu_STAM_VofA, lw=2.0, color='C1', ls='--',
        label='STAM with V(A) = beta/(1-A)')
ax.set_xlabel('redshift z')
ax.set_ylabel('mu (mag)')
ax.set_title('Three cosmologies: where SN sit')
ax.legend()
ax.grid(alpha=0.3)
ax.set_xscale('log')

# Panel B: the "adjustment needed" and what each mechanism can produce
ax = axes[0,1]
ax.plot(z, Delta_mu_needed, lw=3.0, color='k',
        label='Needed: mu_LCDM - mu_EdS')
ax.plot(z, const_Shapiro, lw=2.0, ls='--', color='C2',
        label=f'Shapiro through A_0={A0:.4f}: 2.17*A_0 ~ {2.17*A0:.3f}')
ax.plot(z, linear_bridge, lw=2.0, ls='-.', color='C3',
        label=f'Linear bridge b*z (b=355 Mly)')
ax.plot(z, modFried_shape, lw=2.0, ls=':', color='C1',
        label='Mod-Friedmann shape: 0.105*ln(1+z)')
ax.set_xlabel('redshift z')
ax.set_ylabel('Delta_mu (mag)')
ax.set_title('What is each mechanism able to bridge?')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel C: ratios -- coverage of the needed adjustment
ax = axes[1,0]
ratio_Shapiro = const_Shapiro / np.maximum(Delta_mu_needed, 1e-6)
ratio_bridge  = linear_bridge / np.maximum(Delta_mu_needed, 1e-6)
ratio_modFr   = modFried_shape / np.maximum(Delta_mu_needed, 1e-6)
ax.plot(z, ratio_Shapiro, lw=2.0, color='C2', label='Shapiro coverage')
ax.plot(z, ratio_bridge, lw=2.0, color='C3', label='Linear bridge coverage')
ax.plot(z, ratio_modFr, lw=2.0, color='C1', label='Mod-Friedmann coverage')
ax.axhline(1.0, ls=':', color='k', alpha=0.7, label='full coverage')
ax.set_xlabel('redshift z')
ax.set_ylabel('candidate / needed')
ax.set_title('Fraction of needed adjustment each mechanism covers')
ax.legend()
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.5)
ax.set_xlim(0.05, 2.5)

# Panel D: residual after subtracting each mechanism
ax = axes[1,1]
ax.plot(z, Delta_mu_needed - const_Shapiro, lw=2.0, color='C2',
        label='residual after Shapiro')
ax.plot(z, Delta_mu_needed - linear_bridge, lw=2.0, color='C3',
        label='residual after linear bridge')
ax.plot(z, Delta_mu_needed - modFried_shape, lw=2.0, color='C1',
        label='residual after mod-Friedmann shape')
ax.axhline(0, color='gray', alpha=0.7)
ax.set_xlabel('redshift z')
ax.set_ylabel('residual (mag)')
ax.set_title('What remains unexplained after each adjustment')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle(
    f"Can STAM say 'this is where SN sit, that is why they look different'?\n"
    f"At z=1: needed = {Delta_mu_needed[i_z1]:.2f} mag, "
    f"A_0 photon-A mechanisms cover only ~{100*linear_bridge[i_z1]/Delta_mu_needed[i_z1]:.0f}%",
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.95])

outdir = Path('reports/script_40')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir/'intrinsic_vs_observed.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir/'intrinsic_vs_observed.pdf', bbox_inches='tight')

print(f"\nSaved: {outdir.resolve()}")
print("=" * 72)
