#!/usr/bin/env python3
"""
Script 38 — derived s-scale test

Context:
    Script 37 showed that b can be treated as a TE ledger mismatch:

        b_derived_from_TE = b_path_required - b_TE_remaining

    and that a constrained ratio worked well:

        λ_path = -s
        λ_traverse = k*s
        k ≈ 0.70

Next question:
    Can the remaining scale s be derived rather than freely fitted?

Candidate derivation:
    Original SU threshold anchor:
        z0 = 0.300
        a = 1/z0 = 3.33333

    TE/strong-field intuition:
        path-to-traversal translation carries a correction like:
            1/(1 - A)

    Use the local accumulation proxy at the anchor:
        A_local(z) = 0.70z/(1 + 0.30z)

    Derived candidate:
        s_pred(z0) = (1/z0) / (1 - A_local(z0))

    For z0=0.300:
        A_local ≈ 0.19266
        s_pred ≈ 4.129
    This is close to the Script 37 fitted k=0.70 Pantheon/P+U scale (~4.12).

This script tests:
    1. fitted k=0.70 s per catalog
    2. derived s_pred for z0=0.300, 1/3, 0.350
    3. whether derived s_pred explains b without fitting s
"""

from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import minimize_scalar
    SCIPY_OK=True
except Exception:
    SCIPY_OK=False

H=0.000243635
L_MPC=1/H

def D_geo(z): return L_MPC*z*(1+0.15*z)
def D_excess(z): return 0.35*L_MPC*z*z
def A_avg(z): return 0.35*z/(1+0.15*z)
def A_local(z): return 0.70*z/(1+0.30*z)
def mu(D):
    D=np.maximum(np.asarray(D,float),1e-12)
    return 5*np.log10(D)+25

def find_col(cols,names):
    for t in names:
        for c in cols:
            if c.lower()==t.lower(): return c
    for t in names:
        for c in cols:
            if t.lower() in c.lower(): return c
    return None

def infer(path):
    n=Path(path).stem.lower()
    if "union" in n: return "Union3"
    if "des" in n: return "DES"
    if "pantheon" in n: return "Pantheon"
    return Path(path).stem

def load(path):
    path=Path(path); df=pd.read_csv(path); cols=list(df.columns)
    zcol=find_col(cols,["z","zHD","zCMB","zcmb","redshift"])
    mucol=find_col(cols,["MU","MU_SH0ES","mu","mb","distance_modulus","distmod"])
    if zcol is None or mucol is None: return None
    out=pd.DataFrame({"catalog":infer(path),"z":pd.to_numeric(df[zcol],errors="coerce"),"mu_obs":pd.to_numeric(df[mucol],errors="coerce"),"source_file":path.name})
    out=out.replace([np.inf,-np.inf],np.nan).dropna(subset=["z","mu_obs"])
    return out[out["z"]>0].copy()

def D_path(z, lam_path):
    return D_geo(z)+lam_path*D_excess(z)

def D_TE(z, s, k=0.70):
    return D_geo(z) - s*D_excess(z) + (k*s)*D_geo(z)*A_local(z)

def required_b(y, model_mu): return float(np.mean(y-model_mu))
def rmse(y,pred): return float(np.sqrt(np.mean((y-pred)**2)))

def eval_path(df):
    z=df["z"].values; y=df["mu_obs"].values
    def eval_lam(lam):
        D=D_path(z,lam)
        if np.any(D<=0): return 1e50,0,None
        m=mu(D); b=required_b(y,m); pred=m+b
        return float(np.sum((y-pred)**2)),b,pred
    if SCIPY_OK:
        opt=minimize_scalar(lambda x: eval_lam(x)[0], bounds=(-12,12), method="bounded")
        lam=float(opt.x)
    else:
        grid=np.linspace(-12,12,2001); lam=float(grid[np.argmin([eval_lam(v)[0] for v in grid])])
    sse,b,pred=eval_lam(lam)
    return {"lambda_path":lam,"b_mu":b,"rmse":rmse(y,pred),"sse":sse}

def eval_TE_fixed_s(df,s,k=0.70):
    z=df["z"].values; y=df["mu_obs"].values
    D=D_TE(z,s,k)
    if np.any(D<=0):
        return {"s":s,"k":k,"b_mu":np.nan,"rmse":np.inf,"sse":np.inf}
    m=mu(D); b=required_b(y,m); pred=m+b
    return {"s":s,"k":k,"b_mu":b,"rmse":rmse(y,pred),"sse":float(np.sum((y-pred)**2))}

def fit_TE_s(df,k=0.70):
    z=df["z"].values; y=df["mu_obs"].values
    def eval_s(s):
        D=D_TE(z,s,k)
        if np.any(D<=0): return 1e50,0,None
        m=mu(D); b=required_b(y,m); pred=m+b
        return float(np.sum((y-pred)**2)),b,pred
    if SCIPY_OK:
        opt=minimize_scalar(lambda s: eval_s(s)[0], bounds=(-12,12), method="bounded")
        s=float(opt.x)
    else:
        grid=np.linspace(-12,12,2001); s=float(grid[np.argmin([eval_s(v)[0] for v in grid])])
    sse,b,pred=eval_s(s)
    return {"s":s,"k":k,"b_mu":b,"rmse":rmse(y,pred),"sse":sse}

def decompose(df, path_fit, te_fit):
    z=df["z"].values; y=df["mu_obs"].values
    m_path=mu(D_path(z,path_fit["lambda_path"]))
    m_te=mu(D_TE(z,te_fit["s"],te_fit["k"]))
    b_path=required_b(y,m_path)
    b_te=required_b(y,m_te)
    return {
        "b_path_required":b_path,
        "b_TE_remaining":b_te,
        "b_derived_from_TE":b_path-b_te,
        "fraction_abs_b_explained":(abs(b_path)-abs(b_te))/abs(b_path) if abs(b_path)>1e-12 else np.nan,
        "rmse_path_with_b":rmse(y,m_path+b_path),
        "rmse_TE_with_b":rmse(y,m_te+b_te),
        "rmse_improvement":rmse(y,m_path+b_path)-rmse(y,m_te+b_te),
    }

def s_pred(z0, mode="A_local"):
    a=1/z0
    A = A_local(z0) if mode=="A_local" else A_avg(z0)
    return a/(1-A)

def run(inputs,outdir,zmin=0.05):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    dfs=[load(p) for p in inputs]; dfs=[d for d in dfs if d is not None]
    all_df=pd.concat(dfs,ignore_index=True)
    all_df=all_df[all_df["z"]>=zmin].copy()
    all_df.to_csv(outdir/"script38_loaded_catalog_points.csv",index=False)

    datasets={c:all_df[all_df["catalog"]==c].copy() for c in sorted(all_df["catalog"].unique())}
    if "Pantheon" in datasets and "Union3" in datasets:
        datasets["Pantheon_plus_Union3"]=pd.concat([datasets["Pantheon"],datasets["Union3"]],ignore_index=True)
    if "DES" in datasets and "Pantheon" in datasets:
        zlo=max(datasets["DES"]["z"].min(),datasets["Pantheon"]["z"].min(),zmin)
        zhi=min(datasets["DES"]["z"].max(),datasets["Pantheon"]["z"].max())
        datasets["DES_plus_Pantheon_common"]=pd.concat([
            datasets["DES"][(datasets["DES"]["z"]>=zlo)&(datasets["DES"]["z"]<=zhi)],
            datasets["Pantheon"][(datasets["Pantheon"]["z"]>=zlo)&(datasets["Pantheon"]["z"]<=zhi)]
        ],ignore_index=True)
    datasets["All_three"]=all_df.copy()

    anchors=[("z0_0p300",0.300),("z0_1third",1/3),("z0_0p350",0.350)]
    candidate_rows=[]
    for label,z0 in anchors:
        for mode in ["A_local","A_avg"]:
            candidate_rows.append({"candidate":f"{label}_{mode}","z0":z0,"A_mode":mode,"A_at_z0":float(A_local(z0) if mode=="A_local" else A_avg(z0)),"s_pred":float(s_pred(z0,mode)),"k":0.70})
    candidates=pd.DataFrame(candidate_rows)
    candidates.to_csv(outdir/"script38_s_pred_candidates.csv",index=False)

    rows=[]; fit_rows=[]
    for name,df in datasets.items():
        if len(df)<8: continue
        pf=eval_path(df)
        fit=fit_TE_s(df,k=0.70)
        dec=decompose(df,pf,fit)
        rows.append({"dataset":name,"variant":"fit_s_k0p70","s":fit["s"],"k":fit["k"],"s_error_vs_fit":0.0,**dec})
        fit_rows.append({"dataset":name,"model":"path_fit",**pf})
        fit_rows.append({"dataset":name,"model":"TE_fit_s_k0p70",**fit})
        for _,c in candidates.iterrows():
            tf=eval_TE_fixed_s(df,c["s_pred"],k=0.70)
            dec=decompose(df,pf,tf)
            rows.append({"dataset":name,"variant":c["candidate"],"s":c["s_pred"],"k":0.70,"z0":c["z0"],"A_mode":c["A_mode"],"A_at_z0":c["A_at_z0"],"s_error_vs_fit":c["s_pred"]-fit["s"],**dec})
    results=pd.DataFrame(rows); results.to_csv(outdir/"script38_derived_s_results.csv",index=False)
    fits=pd.DataFrame(fit_rows); fits.to_csv(outdir/"script38_fit_reference.csv",index=False)

    # summarize family
    fam_names=["Pantheon","Union3","Pantheon_plus_Union3"]
    fam=results[results["dataset"].isin(fam_names)]
    summary_rows=[]
    for var,g in fam.groupby("variant"):
        summary_rows.append({
            "variant":var,
            "mean_fraction_abs_b_explained":float(g["fraction_abs_b_explained"].mean()),
            "mean_rmse_improvement":float(g["rmse_improvement"].mean()),
            "mean_s_error_vs_fit":float(g["s_error_vs_fit"].mean()),
            "mean_abs_s_error_vs_fit":float(np.abs(g["s_error_vs_fit"]).mean()),
        })
    summary=pd.DataFrame(summary_rows).sort_values("mean_abs_s_error_vs_fit")
    summary.to_csv(outdir/"script38_family_candidate_summary.csv",index=False)

    # plots
    plot=results[results["variant"].isin(["fit_s_k0p70","z0_0p300_A_local","z0_1third_A_local","z0_0p350_A_local"])].copy()
    if len(plot):
        labels=plot["dataset"]+" | "+plot["variant"]
        plt.figure(figsize=(max(12,len(plot)*0.42),5))
        x=np.arange(len(plot))
        plt.bar(x,plot["fraction_abs_b_explained"])
        plt.axhline(0.75,ls="--",label="75%")
        plt.xticks(x,labels,rotation=75,ha="right",fontsize=8)
        plt.ylabel("fraction of path-only b explained")
        plt.title("Script 38: Derived s candidates vs fitted s")
        plt.legend(); plt.tight_layout()
        plt.savefig(outdir/"script38_b_explained_candidates.png",dpi=180); plt.close()

    # fitted s by dataset vs candidate lines
    fit_s=results[results["variant"]=="fit_s_k0p70"].copy()
    plt.figure(figsize=(10,5))
    plt.bar(fit_s["dataset"],fit_s["s"])
    for _,c in candidates[candidates["A_mode"]=="A_local"].iterrows():
        plt.axhline(c["s_pred"],ls="--",label=c["candidate"])
    plt.xticks(rotation=45,ha="right")
    plt.ylabel("fitted s, k=0.70")
    plt.title("Script 38: fitted s compared to derived A_local candidates")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(outdir/"script38_fitted_s_vs_candidates.png",dpi=180); plt.close()

    verdict={
        "script":"38_derived_s_scale_test.py",
        "purpose":"Test whether the remaining TE scale s can be derived from the z0 threshold anchor and A-local strong-field correction.",
        "candidate_formula":"s_pred(z0) = (1/z0)/(1 - A_local(z0))",
        "z0_0p300_A_local_s_pred":float(candidates[(candidates["candidate"]=="z0_0p300_A_local")]["s_pred"].iloc[0]),
        "family_summary_best_by_s_error":summary.head(10).to_dict(orient="records"),
        "interpretation":[
            "If z0=0.300 with A_local gives s close to fitted Pantheon/P+U s, then the original flag plus local A correction may derive the TE scale.",
            "If this derived s keeps b-explanation high, b is closer to being naturally derived.",
            "DES remains a catalog-family stress case."
        ]
    }
    (outdir/"script38_verdict.json").write_text(json.dumps(verdict,indent=2))
    print("Script 38 derived s-scale test complete.")
    print(f"Output directory: {outdir}")
    print("\nCandidate s values:")
    print(candidates.to_string(index=False))
    print("\nResults:")
    print(results[["dataset","variant","s","s_error_vs_fit","b_path_required","b_TE_remaining","b_derived_from_TE","fraction_abs_b_explained","rmse_improvement"]].to_string(index=False))
    print("\nPantheon/Union family candidate summary:")
    print(summary.to_string(index=False))
    print("\nVerdict:")
    print(json.dumps(verdict,indent=2))
if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--outdir",default="results/script38_derived_s_scale"); p.add_argument("--zmin",type=float,default=0.05)
    a=p.parse_args(); run(a.inputs,a.outdir,zmin=a.zmin)
