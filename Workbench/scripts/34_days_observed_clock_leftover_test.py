#!/usr/bin/env python3
"""
Script 34 — days-observed / clock-leftover test

Purpose:
    Test the "leftover" clock idea.

Question:
    Observed supernova days/stretch are expected to include standard (1+z)
    time dilation. After that correction, is there leftover clock structure
    that tracks STAM A-ledger terms?

Important data limitation:
    The currently available DES/Pantheon files do not contain raw observed
    light-curve durations in days and rest-frame durations as separate columns.
    They do contain SALT2 x1, a light-curve width/shape parameter that is already
    part of SN standardization. Therefore this script treats x1 as a
    clock-leftover proxy, not as raw observed days.

Interpretation:
    If x1/residual shape correlates strongly with A_avg or A_local after normal
    catalog processing, STAM may have a measurable clock-ledger leftover.
    If distance residuals do not improve when adding x1, then the clock term is
    likely already absorbed by standard SN light-curve correction or is not
    exposed in these catalog files.

This is a diagnostic, not a final physical law.
"""

from pathlib import Path
import argparse, json, math, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.optimize import minimize, differential_evolution
    SCIPY_OK=True
except Exception:
    SCIPY_OK=False

H=0.000243635
L_MPC=1/H

def D_geo_mpc(z): return L_MPC*z*(1+0.15*z)
def D_excess_mpc(z): return 0.35*L_MPC*z*z
def A_avg(z): return 0.35*z/(1+0.15*z)
def A_local(z): return 0.70*z/(1+0.30*z)
def mu_from_mpc(D):
    D=np.maximum(np.asarray(D,dtype=float),1e-12)
    return 5*np.log10(D)+25

def find_col(cols,names):
    for t in names:
        for c in cols:
            if c.lower()==t.lower(): return c
    for t in names:
        for c in cols:
            if t.lower() in c.lower(): return c
    return None

def infer_name(p):
    n=Path(p).stem.lower()
    if "des" in n: return "DES"
    if "pantheon" in n: return "Pantheon"
    if "union" in n: return "Union3"
    return Path(p).stem

def load(path):
    path=Path(path); df=pd.read_csv(path); cols=list(df.columns)
    zcol=find_col(cols,["zHD","z","zCMB","redshift","zcmb"])
    mucol=find_col(cols,["MU","MU_SH0ES","mu","distance_modulus","distmod"])
    x1col=find_col(cols,["x1","stretch","s","width","duration","days","tobs","t_obs"])
    if zcol is None or mucol is None:
        return None
    out=pd.DataFrame({
        "catalog":infer_name(path),
        "z":pd.to_numeric(df[zcol],errors="coerce"),
        "mu":pd.to_numeric(df[mucol],errors="coerce"),
        "source_file":path.name,
    })
    if x1col:
        out["clock_proxy"]=pd.to_numeric(df[x1col],errors="coerce")
        out["clock_proxy_col"]=x1col
    else:
        out["clock_proxy"]=np.nan
        out["clock_proxy_col"]=""
    out=out.replace([np.inf,-np.inf],np.nan).dropna(subset=["z","mu"])
    out=out[out["z"]>0].copy()
    out["A_avg"]=A_avg(out["z"].values)
    out["A_local"]=A_local(out["z"].values)
    return out

def linfit(y, X):
    # X includes columns; intercept added outside by caller if desired
    mask=np.all(np.isfinite(X),axis=1)&np.isfinite(y)
    y=y[mask]; X=X[mask]
    if len(y)<5: return None
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    pred=X@beta
    resid=y-pred
    sse=float(np.sum(resid**2)); n=len(y); k=X.shape[1]
    rmse=float(np.sqrt(np.mean(resid**2)))
    aic=n*np.log(max(sse/n,1e-30))+2*k
    return {"n":n,"k":k,"beta":beta.tolist(),"sse":sse,"rmse":rmse,"aic":float(aic)}

def corr(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<5: return np.nan
    return float(np.corrcoef(a[m],b[m])[0,1])

def distance_model(z, params, use_x1=False, x1=None):
    D=D_geo_mpc(z)+params.get("lambda_path",0)*D_excess_mpc(z)+params.get("lambda_traverse",0)*D_geo_mpc(z)*A_local(z)
    mu=mu_from_mpc(D)+params.get("b_mu",0)
    if use_x1 and x1 is not None:
        mu=mu+params.get("gamma_clock_x1",0)*x1
    return mu

def fit_distance(df, use_x1=False):
    z=df["z"].values.astype(float); y=df["mu"].values.astype(float)
    x1=df["clock_proxy"].values.astype(float) if "clock_proxy" in df else np.full_like(y,np.nan)
    if use_x1:
        m=np.isfinite(z)&np.isfinite(y)&np.isfinite(x1)
    else:
        m=np.isfinite(z)&np.isfinite(y)
    z=z[m]; y=y[m]; x1=x1[m]
    if len(y)<10: return None
    names=["lambda_path","lambda_traverse"]+(["gamma_clock_x1"] if use_x1 else [])
    bounds=[(-12,12),(-12,12)]+([(-1,1)] if use_x1 else [])
    def eval_x(v):
        params={n:float(val) for n,val in zip(names,v)}
        D=D_geo_mpc(z)+params.get("lambda_path",0)*D_excess_mpc(z)+params.get("lambda_traverse",0)*D_geo_mpc(z)*A_local(z)
        if np.any(D<=0) or np.any(~np.isfinite(D)): return 1e50,None,None
        mu=mu_from_mpc(D)
        if use_x1: mu=mu+params.get("gamma_clock_x1",0)*x1
        b=float(np.mean(y-mu)); params["b_mu"]=b
        pred=mu+b
        sse=float(np.sum((y-pred)**2))
        return sse,params,pred
    if SCIPY_OK:
        de=differential_evolution(lambda v: eval_x(v)[0], bounds, seed=34, maxiter=80, tol=1e-7, polish=False)
        loc=minimize(lambda v: eval_x(v)[0], de.x, method="Nelder-Mead", options={"maxiter":1000})
        x=loc.x if loc.fun<=de.fun else de.x
    else:
        # simple fallback
        grids=[np.linspace(a,b,61) for a,b in bounds]
        best=(1e99,None)
        import itertools
        for v in itertools.product(*grids):
            s,_,_=eval_x(v)
            if s<best[0]: best=(s,v)
        x=np.array(best[1])
    sse,params,pred=eval_x(x)
    n=len(y); k=len(names)+1
    rmse=float(np.sqrt(np.mean((y-pred)**2)))
    aic=n*np.log(max(sse/n,1e-30))+2*k
    return {"n":n,"use_x1":use_x1,"rmse":rmse,"aic":float(aic),"sse":sse,**params}

def run(inputs,outdir,zmin=0.05):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    dfs=[]
    for p in inputs:
        d=load(p)
        if d is not None: dfs.append(d)
    all_df=pd.concat(dfs,ignore_index=True)
    all_df=all_df[all_df["z"]>=zmin].copy()
    all_df.to_csv(outdir/"script34_loaded_points.csv",index=False)

    datasets={c:all_df[all_df["catalog"]==c].copy() for c in sorted(all_df["catalog"].unique())}
    if "Pantheon" in datasets and "DES" in datasets:
        datasets["DES_plus_Pantheon"]=pd.concat([datasets["DES"],datasets["Pantheon"]],ignore_index=True)

    clock_rows=[]; dist_rows=[]; residual_rows=[]
    for name,df in datasets.items():
        dfx=df.dropna(subset=["clock_proxy"]).copy()
        if len(dfx)>=10:
            y=dfx["clock_proxy"].values
            models={
                "const":np.column_stack([np.ones(len(dfx))]),
                "z":np.column_stack([np.ones(len(dfx)),dfx["z"].values]),
                "A_avg":np.column_stack([np.ones(len(dfx)),dfx["A_avg"].values]),
                "A_local":np.column_stack([np.ones(len(dfx)),dfx["A_local"].values]),
                "z_plus_Alocal":np.column_stack([np.ones(len(dfx)),dfx["z"].values,dfx["A_local"].values]),
            }
            for split in [0.30,1/3,0.35]:
                hinge=np.maximum(0,dfx["z"].values-split)
                models[f"hinge_{split:.6f}"]=np.column_stack([np.ones(len(dfx)),dfx["z"].values,hinge])
            for mname,X in models.items():
                fit=linfit(y,X)
                if fit:
                    clock_rows.append({"dataset":name,"clock_proxy_col":",".join(sorted(dfx["clock_proxy_col"].unique())),"clock_model":mname,**{k:v for k,v in fit.items() if k!="beta"},"beta":json.dumps(fit["beta"])})
            clock_rows.append({"dataset":name,"clock_proxy_col":",".join(sorted(dfx["clock_proxy_col"].unique())),"clock_model":"correlations","n":len(dfx),"k":0,"sse":np.nan,"rmse":np.nan,"aic":np.nan,"beta":json.dumps({
                "corr_x1_z":corr(dfx["clock_proxy"].values,dfx["z"].values),
                "corr_x1_Aavg":corr(dfx["clock_proxy"].values,dfx["A_avg"].values),
                "corr_x1_Alocal":corr(dfx["clock_proxy"].values,dfx["A_local"].values),
            })})

            # Distance fits M3 vs M3+x1
            f0=fit_distance(dfx,use_x1=False); f1=fit_distance(dfx,use_x1=True)
            if f0 and f1:
                dist_rows.append({"dataset":name,"model":"M3_path_traverse","delta_vs_no_x1_rmse":0,**f0})
                rec={"dataset":name,"model":"M3_plus_clock_proxy_x1","delta_vs_no_x1_rmse":f1["rmse"]-f0["rmse"],"delta_vs_no_x1_aic":f1["aic"]-f0["aic"],**f1}
                dist_rows.append(rec)
                # residual correlation from no-x1 fit
                z=dfx["z"].values; ymu=dfx["mu"].values; x1=dfx["clock_proxy"].values
                pred=distance_model(z,f0,use_x1=False)
                resid=ymu-pred
                residual_rows.append({"dataset":name,"n":len(dfx),"corr_resid_x1":corr(resid,x1),"corr_resid_Alocal":corr(resid,dfx["A_local"].values),"no_x1_rmse":f0["rmse"],"with_x1_rmse":f1["rmse"],"rmse_delta":f1["rmse"]-f0["rmse"],"aic_delta":f1["aic"]-f0["aic"]})

            # Plots
            plt.figure(figsize=(9,5))
            plt.scatter(dfx["z"],dfx["clock_proxy"],s=12,alpha=0.55)
            for xv,lab in [(0.30,"0.30"),(1/3,"1/3"),(0.35,"0.35")]:
                plt.axvline(xv,ls="--",lw=1,label=f"z={lab}")
            plt.title(f"Script 34 clock-leftover proxy vs z: {name}")
            plt.xlabel("z"); plt.ylabel("clock proxy / SALT2 x1")
            plt.legend(); plt.tight_layout()
            plt.savefig(outdir/f"script34_{name}_clock_proxy_vs_z.png",dpi=170); plt.close()

    clock=pd.DataFrame(clock_rows); clock.to_csv(outdir/"script34_clock_proxy_model_fits.csv",index=False)
    dist=pd.DataFrame(dist_rows); dist.to_csv(outdir/"script34_distance_fit_with_clock_proxy.csv",index=False)
    resids=pd.DataFrame(residual_rows); resids.to_csv(outdir/"script34_residual_clock_proxy_correlations.csv",index=False)

    # Best clock model summary
    best_clock=clock[clock["clock_model"]!="correlations"].sort_values("aic").groupby("dataset",as_index=False).first() if len(clock) else pd.DataFrame()
    best_clock.to_csv(outdir/"script34_best_clock_proxy_models.csv",index=False)

    # Plot distance improvement
    if len(dist):
        dx=dist[dist["model"]=="M3_plus_clock_proxy_x1"].copy()
        if len(dx):
            plt.figure(figsize=(8,5))
            plt.bar(dx["dataset"],dx["delta_vs_no_x1_rmse"])
            plt.axhline(0,lw=1)
            plt.ylabel("RMSE change adding x1 clock proxy")
            plt.title("Script 34: Does clock proxy improve distance ledger?")
            plt.xticks(rotation=30,ha="right")
            plt.tight_layout()
            plt.savefig(outdir/"script34_x1_distance_rmse_delta.png",dpi=170); plt.close()

    verdict={
        "script":"34_days_observed_clock_leftover_test.py",
        "data_limitation":"The available files do not include raw observed days and rest-frame days. DES/Pantheon include SALT2 x1, used here as a clock-leftover/light-curve-width proxy after standard time-dilation handling.",
        "core_question":"After standard (1+z) light-curve handling, is there leftover clock structure that tracks A_avg/A_local or helps the distance ledger?",
        "main_outputs":["script34_clock_proxy_model_fits.csv","script34_distance_fit_with_clock_proxy.csv","script34_residual_clock_proxy_correlations.csv"],
        "interpretation_rules":[
            "If x1 strongly tracks A and improves residuals, clock ledger may be exposed in SN shape leftovers.",
            "If x1 does not improve distance residuals, the clock term is likely already standardized out or not exposed in these files.",
            "Raw observed days would be a better test than x1."
        ]
    }
    (outdir/"script34_verdict.json").write_text(json.dumps(verdict,indent=2))
    print("Script 34 days-observed / clock-leftover test complete.")
    print(f"Output directory: {outdir}")
    print("\nBest clock proxy models:")
    print(best_clock[["dataset","clock_model","rmse","aic","beta"]].to_string(index=False) if len(best_clock) else "No clock proxy data found.")
    print("\nDistance fits with clock proxy:")
    print(dist[["dataset","model","rmse","aic","b_mu","lambda_path","lambda_traverse","gamma_clock_x1","delta_vs_no_x1_rmse"]].to_string(index=False) if len(dist) else "No distance fits.")
    print("\nResidual/clock correlations:")
    print(resids.to_string(index=False) if len(resids) else "No residual correlations.")
    print("\nMain caution: this used SALT2 x1 as a leftover proxy, not raw observed days.")
if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--inputs",nargs="+",required=True); p.add_argument("--outdir",default="results/script34_days_observed_clock_leftover"); p.add_argument("--zmin",type=float,default=0.05)
    a=p.parse_args(); run(a.inputs,a.outdir,zmin=a.zmin)
