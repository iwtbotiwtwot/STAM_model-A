#!/usr/bin/env python3
"""Generate STAM_model-A distance and path-accumulation tables.

This baseline script intentionally avoids plotting dependencies. Plotting can be
added later, but the first repo scripts should run in minimal environments.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from stam_model_a.constants import B_CANDIDATES, L_STAM
from stam_model_a.cosmology import A_path_average, A_path_local, D_adj, D_excess, D_geo


def main() -> None:
    out_dir = Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    z = np.linspace(0.0, 2.5, 251)
    table = pd.DataFrame({"z": z, "D_geo": D_geo(z, L=L_STAM)})

    for name, b in B_CANDIDATES.items():
        table[f"D_adj_{name}"] = D_adj(z, b=b, L=L_STAM)
        table[f"D_excess_{name}"] = D_excess(z, b=b, L=L_STAM)
        table[f"A_path_average_{name}"] = A_path_average(z, b=b, L=L_STAM)
        table[f"A_path_local_{name}"] = A_path_local(z, b=b, L=L_STAM)

    out_path = out_dir / "stam_distance_and_path_accumulation_table.csv"
    table.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
