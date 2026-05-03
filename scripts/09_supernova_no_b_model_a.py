#!/usr/bin/env python3
"""Run supernova catalogs against no-b STAM Model-A distance prediction.

Primary no-b physical Model-A prediction:
    D_adj,0(z) = L z (1 + 0.5 z)

Comparison geometric spine:
    D_geo(z) = L z (1 + 0.15 z)

Historical diagnostic only:
    D_adj,b(z) = L z (1 + 0.5 z) + b z
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

C_MLY_PER_MPC = 3.261563776
H_STAM = 0.000243635
L_MLY = C_MLY_PER_MPC / H_STAM
B_HIST = 354.95

def d_geo_mly(z):
    z=np.asarray(z,dtype=float); return L_MLY*z*(1.0+0.15*z)
def d_adj0_mly(z):
    z=np.asarray(z,dtype=float); return L_MLY*z*(1.0+0.5*z)
def d_adj_b_mly(z,b=B_HIST):
    z=np.asarray(z,dtype=float); return L_MLY*z*(1.0+0.5*z)+b*z
def mu_from_mly(d_mly):
    return 5.0*np.log10(np.asarray(d_mly,dtype=float)/C_MLY_PER_MPC)+25.0
def mly_from_mu(mu):
    return C_MLY_PER_MPC*10.0**((np.asarray(mu,dtype=float)-25.0)/5.0)

def load_catalog(name: str, path: Path) -> pd.DataFrame:
    df=pd.read_csv(path)
    if name.lower()=="union3":
        zcol,mucol,errcol,idcol="z","mb",None,"bin"
    elif name.lower()=="pantheon":
        zcol,mucol,errcol,idcol="zHD","MU_SH0ES","MU_SH0ES_ERR_DIAG","CID"
    elif name.lower()=="des":
        zcol,mucol,errcol,idcol="zHD","MU","MUERR","CID"
    else:
        raise ValueError(f"unknown catalog preset: {name}")
    out=pd.DataFrame({
        "catalog":name,
        "object_id":df[idcol].astype(str) if idcol in df.columns else np.arange(len(df)).astype(str),
        "source_row":np.arange(len(df)),
        "z":pd.to_numeric(df[zcol],errors="coerce"),
        "mu_obs":pd.to_numeric(df[mucol],errors="coerce"),
    })
    out["mu_err"]=pd.to_numeric(df[errcol],errors="coerce") if errcol and errcol in df.columns else np.nan
    return out[np.isfinite(out["z"]) & np.isfinite(out["mu_obs"]) & (out["z"]>0)].copy()

def add_model_columns(df: pd.DataFrame) -> pd.DataFrame:
    out=df.copy()
    out["D_obs_Mly"]=mly_from_mu(out["mu_obs"])
    out["D_geo_Mly"]=d_geo_mly(out["z"])
    out["D_adj0_Mly"]=d_adj0_mly(out["z"])
    out["D_adj_b35495_Mly"]=d_adj_b_mly(out["z"])
    out["mu_geo"]=mu_from_mly(out["D_geo_Mly"])
    out["mu_adj0"]=mu_from_mly(out["D_adj0_Mly"])
    out["mu_adj_b35495"]=mu_from_mly(out["D_adj_b35495_Mly"])
    out["resid_mu_geo"]=out["mu_obs"]-out["mu_geo"]
    out["resid_mu_adj0"]=out["mu_obs"]-out["mu_adj0"]
    out["resid_mu_adj_b35495"]=out["mu_obs"]-out["mu_adj_b35495"]
    out["distance_ratio_obs_to_adj0"]=out["D_obs_Mly"]/out["D_adj0_Mly"]
    out["b_bridge_from_no_b"]=(out["D_obs_Mly"]-out["D_adj0_Mly"])/out["z"]
    return out

def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for cat,g in df.groupby("catalog"):
        for model,col in [("D_geo/SU","resid_mu_geo"),("Model-A no-b D_adj0","resid_mu_adj0"),("historical b=354.95 diagnostic","resid_mu_adj_b35495")]:
            vals=g[col].to_numpy(float)
            rows.append({"catalog":cat,"n":int(len(g)),"model":model,
                         "residual_mean_mag":float(np.mean(vals)),
                         "residual_median_mag":float(np.median(vals)),
                         "residual_std_mag":float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0,
                         "residual_rmse_mag":float(np.sqrt(np.mean(vals**2)))})
        bvals=g["b_bridge_from_no_b"].replace([np.inf,-np.inf],np.nan).dropna().to_numpy()
        rows.append({"catalog":cat,"n":int(len(g)),"model":"b_bridge_from_no_b",
                     "b_bridge_mean_Mly":float(np.mean(bvals)),
                     "b_bridge_median_Mly":float(np.median(bvals)),
                     "b_bridge_std_Mly":float(np.std(bvals,ddof=1)) if len(bvals)>1 else 0.0})
    return pd.DataFrame(rows)

def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--union3",default=Path("data/union3_bins.csv"),type=Path)
    parser.add_argument("--pantheon",default=Path("data/pantheon.csv"),type=Path)
    parser.add_argument("--des",default=Path("data/des.csv"),type=Path)
    parser.add_argument("--z-min",default=0.05,type=float)
    parser.add_argument("--z-max",default=1.14418,type=float)
    parser.add_argument("--output-dir",default=Path("results/no_b_supernova"),type=Path)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    all_df=pd.concat([load_catalog("Union3",args.union3),load_catalog("Pantheon",args.pantheon),load_catalog("DES",args.des)],ignore_index=True)
    all_df=add_model_columns(all_df)
    common=all_df[(all_df["z"]>=args.z_min)&(all_df["z"]<=args.z_max)].copy()
    all_df.to_csv(args.output_dir/"supernova_no_b_pointwise_all.csv",index=False)
    common.to_csv(args.output_dir/"supernova_no_b_pointwise_common_range.csv",index=False)
    summarize(common).to_csv(args.output_dir/"supernova_no_b_catalog_summary_common_range.csv",index=False)
    summary={"test":"Supernova catalogs vs no-b STAM Model-A distance prediction",
             "primary_prediction":"D_adj0(z)=Lz(1+0.5z)",
             "comparison":"D_geo(z)=Lz(1+0.15z)",
             "historical_diagnostic":"D_adj,b(z)=Lz(1+0.5z)+354.95z",
             "z_min":args.z_min,"z_max":args.z_max,
             "n_by_catalog_common_range":{k:int(v) for k,v in common.groupby("catalog").size().items()}}
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
    print("PASS: ran supernova catalogs against no-b Model-A distance prediction.")
    print(json.dumps(summary,indent=2,sort_keys=True))
    print(f"Wrote outputs to {args.output_dir}")

if __name__=="__main__":
    main()
