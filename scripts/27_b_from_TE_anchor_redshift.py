#!/usr/bin/env python3
"""STAM Model-A b derivation from TE anchor redshift.

Hypothesis:
  b is the linear shadow/tangent of quadratic Traversal Excess (TE).

TE(z) = 0.35Lz^2
dTE/dz = 0.70Lz

Candidate:
  b_pred = 0.70L z_anchor

If z_anchor comes from catalog/calibration structure rather than fitting b,
then b can be derived instead of fitted.
"""

from __future__ import annotations

import json, math
from pathlib import Path
import numpy as np
import pandas as pd

C_MLY_PER_MPC=3.261563776
H_STAM=0.000243635
L_MLY=C_MLY_PER_MPC/H_STAM

HISTORICAL_B={"Pantheon_Union_style":354.95,"Original_retained":461.3626922,"DES_style":1335.412792}

def b_pred(z): return 0.70*L_MLY*z

def weighted_mean(v,w=None):
    v=np.asarray(v,dtype=float)
    if w is None: return float(np.mean(v))
    w=np.asarray(w,dtype=float); m=np.isfinite(v)&np.isfinite(w)&(w>0)
    return float(np.sum(v[m]*w[m])/np.sum(w[m]))

def load_catalog(name,path):
    df=pd.read_csv(path)
    if name=="Union3":
        z=pd.to_numeric(df["z"],errors="coerce")
        err=np.full(len(df),np.nan)
    elif name=="Pantheon":
        z=pd.to_numeric(df["zHD"],errors="coerce")
        err=pd.to_numeric(df["MU_SH0ES_ERR_DIAG"],errors="coerce")
    else:
        z=pd.to_numeric(df["zHD"],errors="coerce")
        err=pd.to_numeric(df["MUERR"],errors="coerce")
    out=pd.DataFrame({"catalog":name,"z":z,"err":err})
    out=out[np.isfinite(out.z)&(out.z>0)].copy()
    e=out.err.to_numpy(float)
    m=np.isfinite(e)&(e>0)
    out["weight"]=1.0
    if m.any():
        med=np.nanmedian(e[m])
        ef=np.where(m,e,med)
        out["weight"]=1/(ef*ef)
    return out

def main():
    out=Path("results/b_from_TE_anchor")
    out.mkdir(parents=True,exist_ok=True)
    paths={"Union3":Path("data/union3_bins.csv"),"Pantheon":Path("data/pantheon.csv"),"DES":Path("data/des.csv")}
    frames=[load_catalog(k,p) for k,p in paths.items() if p.exists()]
    if not frames:
        summary={"test":"b from TE anchor","status":"no catalog files found"}
        (out/"summary.json").write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2)); return
    cat=pd.concat(frames,ignore_index=True)
    rows=[]
    groups={k:g for k,g in cat.groupby("catalog")}
    groups["Pantheon_Union3"]=cat[cat.catalog.isin(["Pantheon","Union3"])]
    groups["All_three"]=cat
    for group,g in groups.items():
        anchors={"mean_z":float(g.z.mean()),"median_z":float(g.z.median()),"weighted_mean_z":weighted_mean(g.z,g.weight)}
        for cut in [0.03,0.05,0.075,0.1,0.15,0.2,0.3]:
            h=g[g.z<=cut]
            if len(h):
                anchors[f"mean_z_le_{cut}"]=float(h.z.mean())
                anchors[f"weighted_mean_z_le_{cut}"]=weighted_mean(h.z,h.weight)
        for an,z in anchors.items():
            for target,b in HISTORICAL_B.items():
                pred=b_pred(z)
                rows.append({"group":group,"anchor_name":an,"z_anchor":z,"b_pred_Mly":pred,"target":target,"target_b_Mly":b,"abs_error_Mly":abs(pred-b),"pct_error":100*(pred-b)/b})
    res=pd.DataFrame(rows).sort_values("abs_error_Mly")
    res.to_csv(out/"b_prediction_vs_historical_targets.csv",index=False)
    reverse=[{"b_label":k,"b_Mly":v,"z_anchor_if_b_eq_0p70Lz":v/(0.70*L_MLY)} for k,v in HISTORICAL_B.items()]
    pd.DataFrame(reverse).to_csv(out/"historical_b_reverse_anchor_redshifts.csv",index=False)
    summary={"test":"b from TE anchor redshift","formula":"b_pred=0.70L z_anchor","best_matches":res.head(20).to_dict(orient="records"),"reverse_anchors":reverse}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
