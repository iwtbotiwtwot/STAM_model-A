# Parameter Locking

Parameter locking prevents accidental overfitting.

## Candidate parameters

Current historical candidate `b` values:

```text
pantheon_union_style = 354.95
original_retained    = 461.3626922
des_style            = 1335.412792
```

Model scale:

```text
C = 3.261563776
H = 0.000243635
L = C/H
```

## Rule for empirical tests

A test must be labeled one of:

```text
identity_check
fit
holdout_prediction
post_hoc_diagnostic
```

Only `holdout_prediction` tests count as prediction.

## Lock file recommendation

Before running a holdout test, create:

```text
results/locks/<YYYY-MM-DD>_<test_name>.json
```

Example:

```json
{
  "model": "STAM_model-A",
  "test": "Pantheon-trained Union3 holdout",
  "fit_catalog": "Pantheon",
  "holdout_catalog": "Union3",
  "locked_parameters": {
    "b": 354.95,
    "L": 13387.090426252385
  },
  "pass_fail_thresholds": {
    "mean_abs_mu_residual_max": null,
    "redshift_slope_abs_max": null,
    "reduced_chi2_max": null
  },
  "notes": "Fill thresholds before running."
}
```

## Do not do this

Do not fit `b` on a catalog and call the same catalog an independent validation.

Do not change preprocessing after viewing residuals unless the change is recorded as a new test version.

Do not compare DES/Pantheon/Union3 by eye only; compute the implied `b(z)` and residual trends.
