# STAM Model-A galactic cumulative A scratch test

## Purpose

This scratch test explores the idea:

```text
Many individually tiny A contributions can sum into a galaxy-scale accumulation field.
```

Linear Model-A version:

```text
A_total(x)=Σ 2Gm_i/(c²|x-x_i|)
```

## Important distinction

Linear cumulative A gives:

```text
A_total = sum of component A values
```

The stronger idea:

```text
A_total > linear sum
```

requires an additional collective/coherence/envelope term:

```text
A_total = A_linear + A_collective
```

The script includes this only as an exploratory toy.

## Run

```powershell
python scripts\15_galactic_cumulative_A_scratch.py
```

Outputs:

```text
results\galactic_cumulative_A\summary.json
results\galactic_cumulative_A\high_precision_A_scale_check.csv
results\galactic_cumulative_A\linear_cumulative_visible_disk_toy.csv
results\galactic_cumulative_A\collective_envelope_toy.csv
```
