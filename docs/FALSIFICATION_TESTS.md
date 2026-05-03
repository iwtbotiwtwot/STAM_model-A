# Falsification Tests

The purpose of this document is to define tests that could count against STAM_model-A. Supportive fits are not enough.

## Test class 1 — Internal identity checks

These checks should pass before any empirical work.

### FT-01: Local gravity identity

Input:

```text
A(r)=2GM/(c^2 r)
```

Prediction:

```text
(c^2/2)dA/dr = -GM/r^2
```

Fail condition:

The implemented functions do not reproduce the identity to numerical precision.

Status: implemented in `scripts/00_validate_model_a.py` and unit tests.

### FT-02: Horizon identity

Input:

```text
r = Rs = 2GM/c^2
```

Prediction:

```text
A(r)=1
```

Fail condition:

The implemented local accumulation function does not return unity at `r=Rs`.

Status: implemented.

### FT-03: Cosmological path-integral identity

Input:

```text
D_geo=Lz(1+0.15z)
A_path,local=(b/L+0.70z)/(1+0.30z)
```

Prediction:

```text
int_0^z [1+A_path,local(u)] dD_geo/du du
= Lz(1+0.5z)+bz
```

Fail condition:

Numerical integral and analytic `D_adj` disagree beyond tolerance.

Status: implemented.

## Test class 2 — Locked supernova prediction

### FT-04: Fit one catalog, predict another

Procedure:

1. Choose one training catalog.
2. Record dataset version and preprocessing.
3. Fit `b` only on the training catalog.
4. Freeze `b` in `PRIORITY_RECORD.md` or a locked parameter file.
5. Predict one or more holdout catalogs without retuning.
6. Report residuals versus redshift.

Fail condition:

A predeclared holdout residual threshold is exceeded, or residuals show a systematic redshift trend not accounted for by the model.

Recommended first thresholds to define before running:

```text
mean residual threshold:
redshift-slope threshold:
chi-square / reduced chi-square threshold:
outlier handling rule:
```

Status: template script included; thresholds not yet locked.

### FT-05: Redshift-bin `b` stability

For each supernova or distance indicator with observed adjusted distance `D_obs(z)`, infer:

```text
b_inferred(z) = [D_obs(z) - Lz(1+0.5z)] / z
```

Then bin by redshift.

Fail condition:

`b_inferred` is not approximately catalog-constant, or the required drift is larger than the model permits.

Important:

This test is stronger than a global fit because it asks whether a single `b` actually describes the catalog.

Status: template script included.

### FT-06: Catalog-family compatibility

Known candidate values:

```text
Pantheon/Union-style b ~= 355
DES-style b            ~= 1335
```

Fail condition:

The model requires mutually incompatible `b` values and no independent calibration/path interpretation explains the discrepancy.

Status: open.

## Test class 3 — BAO / angular diameter / CMB consistency

### FT-07: Distance-duality check

If STAM changes luminosity-distance interpretation, test whether it remains compatible with:

```text
D_L = (1+z)^2 D_A
```

or specify why the relation is modified.

Fail condition:

The model fits luminosity distances but cannot produce compatible angular-diameter distances.

Status: open.

### FT-08: BAO scale prediction

Procedure:

1. Use locked STAM parameters.
2. Predict radial and transverse distance measures used by BAO.
3. Compare to BAO measurements without retuning `b`.

Fail condition:

The same geometry/path model cannot match supernova and BAO distances.

Status: open.

### FT-09: CMB acoustic-scale compatibility

Procedure:

1. Define STAM distance to last scattering.
2. Compare to the angular acoustic scale.
3. Avoid post-hoc fitting unless explicitly labeled as such.

Fail condition:

The model cannot reproduce the observed angular scale or requires parameters inconsistent with lower-redshift tests.

Status: open.

## Test class 4 — Local relativistic tests

### FT-10: Solar System Shapiro coefficient

STAM's delay form:

```text
Delta t = (2GM/c^3) int ds/r
```

must be checked against the standard weak-field coefficient used in Solar System tests.

Fail condition:

Coefficient-level disagreement beyond observational bounds, unless the model gives a principled correction that also passes data.

Status: open beyond algebraic structure.

### FT-11: Gravitational redshift / clock tests

Need an explicit clock-rate prescription from `A`, for example a weak-field relation such as:

```text
z_grav ~ A/2
```

or another derived form.

Fail condition:

The chosen prescription fails laboratory, GPS, or Solar System clock measurements.

Status: open.

### FT-12: Lensing closure

Acceleration, Shapiro delay, and lensing must be predicted by one field normalization.

Fail condition:

The model requires inconsistent factors to match dynamics and lensing.

Status: open.

## Test class 5 — Horizon regime

### FT-13: Invariant horizon criterion

The statement `A=1` must be translated into a coordinate-invariant or operational criterion.

Fail condition:

`A=1` is only a coordinate artifact or cannot reproduce known horizon behavior.

Status: open.

## Reporting rule

Every falsification script should write a machine-readable result file:

```text
results/<test_name>_<date>.json
```

with:

```text
commit_hash
input_files
parameters
fit_or_holdout
metric_values
pass_fail_status
notes
```
