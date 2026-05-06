# STAM Model-A supernova light-curve time-dilation proxy diagnostic

## Purpose

This test is a first step toward the STAM time-interpretation question.

STAM treats time measurements as physical process-rate/traversal observables affected by `A`, not as an independent substance.

A direct supernova time-dilation test would require:

```text
observed-frame light-curve width
rest-frame/template light-curve width
redshift
```

The uploaded Pantheon and DES catalog files contain SALT2 `x1`, which is a fitted light-curve shape parameter, not raw observed-frame duration.

## Status

This script is therefore a proxy/catalog diagnostic only.

It asks:

```text
Does x1 show strong redshift structure?
Does x1 correlate with no-b Model-A distance residuals?
What raw light-curve duration schema is needed for the real test?
```

## Run

```powershell
python scripts\19_sn_lightcurve_time_proxy.py `
  --pantheon data\pantheon.csv `
  --des data\des.csv
```

Outputs:

```text
results\sn_lightcurve_time_proxy\summary.json
results\sn_lightcurve_time_proxy\sn_lightcurve_proxy_summary.csv
results\sn_lightcurve_time_proxy\future_raw_lightcurve_duration_schema.csv
```
