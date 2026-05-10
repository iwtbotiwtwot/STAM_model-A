"""Quick diagnostic to verify Pantheon+ STAT+SYS covariance ordering matches pantheon.csv."""
import numpy as np
import pandas as pd

with open('data/Pantheon+SH0ES_STAT+SYS.cov', 'r') as f:
    N = int(f.readline().strip())
    flat = np.fromfile(f, sep='\n', dtype=np.float64, count=N * N)
C = flat.reshape(N, N)

asym = np.abs(C - C.T)
print(f"Max asymmetry: {asym.max():.3e}")
print(f"Mean asymmetry: {asym.mean():.3e}")
print(f"Symmetric to 1e-6: {np.allclose(C, C.T, atol=1e-6)}")
print(f"Symmetric to 1e-12: {np.allclose(C, C.T, atol=1e-12)}")
print()

p = pd.read_csv('data/pantheon.csv')
mask = (p['IS_CALIBRATOR'] == 0) & (p['zCMB'] > 0.05)
idx = np.where(mask)[0][:8]

print("Compare per-SN diag at z>0.05 (PV negligible):")
print(f"{'i':<4} {'CID':<18} {'z':<8} {'sigma_csv':<12} {'sqrt(Cii)':<12} {'ratio':<8}")
for i in idx:
    cid = p['CID'].iloc[i]
    z = p['zCMB'].iloc[i]
    sigma_csv = p['MU_SH0ES_ERR_DIAG'].iloc[i]
    sigma_cov = np.sqrt(C[i, i])
    print(f"{i:<4} {cid:<18} {z:<8.4f} {sigma_csv:<12.4f} {sigma_cov:<12.4f} "
          f"{sigma_cov / sigma_csv:.3f}")

print()
print("Now check whether after subtracting PV variance, CSV matches Cov diag:")
# Pantheon+ peculiar velocity uncertainty: typically 250 km/s -> sigma_mu_pv = (5/ln10) * v_pec / cz
c_kms = 299792.458
v_pec = 250.0
print(f"{'i':<4} {'CID':<18} {'z':<8} {'sigma_csv':<10} {'sigma_pv':<10} "
      f"{'diff^2':<12} {'Cii':<12}")
for i in idx:
    cid = p['CID'].iloc[i]
    z = p['zCMB'].iloc[i]
    sigma_csv = p['MU_SH0ES_ERR_DIAG'].iloc[i]
    sigma_pv = (5.0 / np.log(10.0)) * v_pec / (c_kms * z)
    sigma_csv2 = sigma_csv**2
    sigma_pv2 = sigma_pv**2
    diff2 = sigma_csv2 - sigma_pv2
    Cii = C[i, i]
    print(f"{i:<4} {cid:<18} {z:<8.4f} {sigma_csv:<10.4f} {sigma_pv:<10.4f} "
          f"{diff2:<12.4f} {Cii:<12.4f}")
