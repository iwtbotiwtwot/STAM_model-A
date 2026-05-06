# STAM Model-A b derivation from TE anchor redshift

Hypothesis:

```text
b is the linear shadow/tangent of quadratic TE.
```

Core formulas:

```text
TE(z) = 0.35Lz²
dTE/dz = 0.70Lz
b_pred = 0.70L z_anchor
```

If `z_anchor` is obtained from catalog/calibration structure instead of fitting `b`, then `b` becomes derived.

Run:

```powershell
python scripts\27_b_from_TE_anchor_redshift.py
```
