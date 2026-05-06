#!/usr/bin/env python3
from pathlib import Path
import argparse, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def visibility_from_A_information(E, P, Env, erase_path=False, erase_environment=False):
    # Existence receipt E does NOT resolve path-A by itself.
    P_eff = 0.0 if erase_path else P
    Env_eff = 0.0 if erase_environment else Env
    R_A = float(np.clip(max(P_eff, Env_eff), 0.0, 1.0))
    V = math.sqrt(max(0.0, 1.0 - R_A*R_A))
    return V, R_A, float(np.clip(E, 0.0, 1.0))

def double_slit_probability(x, V=1.0, wavelength=1.0, slit_sep=10.0, slit_width=1.8, L=120.0):
    slit1, slit2 = -slit_sep/2, slit_sep/2
    r1 = np.sqrt(L**2 + (x-slit1)**2)
    r2 = np.sqrt(L**2 + (x-slit2)**2)
    k = 2*math.pi/wavelength
    theta = np.arctan2(x, L)
    beta = math.pi*slit_width*np.sin(theta)/wavelength
    envelope = np.sinc(beta/math.pi)**2
    amp = np.sqrt(np.maximum(envelope,0))/math.sqrt(2)
    psi1 = amp*np.exp(1j*k*r1)
    psi2 = amp*np.exp(1j*k*r2)
    P_incoh = np.abs(psi1)**2 + np.abs(psi2)**2
    cross = 2*np.real(psi1*np.conj(psi2))
    P = np.maximum(P_incoh + V*cross, 0)
    if np.max(P) > 0: P = P/np.max(P)
    env = envelope/np.max(envelope)
    return P, env

def contrast(x, y, window=25.0):
    sub = y[np.abs(x) <= window]
    return float((np.max(sub)-np.min(sub))/(np.max(sub)+np.min(sub)))

def behavior(V):
    if V > 0.85: return "quantum-like / strong interference"
    if V > 0.35: return "partial / weakened interference"
    return "classical-like / path-resolved"

def run(outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-55,55,1200)
    cases = [
        ("A_uncataloged_particle", "No existence receipt, no path receipt", 0,0,0,False,False),
        ("B_existence_known_path_unknown", "Existence known, path-A unknown", 1,0,0,False,False),
        ("C_path_measured", "Path-A measured", 1,1,0,False,False),
        ("D_environment_keeps_receipt", "Environment keeps receipt", 1,0,1,False,False),
        ("E_weak_environment_receipt", "Weak environment receipt", 1,0,0.45,False,False),
        ("F_path_receipt_erased", "Path receipt erased", 1,1,0,True,False),
        ("G_environment_receipt_erased", "Environment receipt erased", 1,0,1,False,True),
        ("H_known_existence_and_endpoint_only", "Known emission + endpoint, path unknown", 1,0,0,False,False),
    ]
    rows=[]
    plt.figure(figsize=(11,7))
    key={"B_existence_known_path_unknown","C_path_measured","D_environment_keeps_receipt","E_weak_environment_receipt","F_path_receipt_erased"}
    for name,label,E,P,Env,eraseP,eraseEnv in cases:
        V,R_A,E2 = visibility_from_A_information(E,P,Env,eraseP,eraseEnv)
        y,env = double_slit_probability(x,V)
        rows.append({
            "case":name, "label":label, "existence_receipt_E":E2, "path_receipt_P":P,
            "environment_receipt_Env":Env, "erase_path":eraseP, "erase_environment":eraseEnv,
            "resolved_A_information_R_A":R_A, "visibility_V":V,
            "central_contrast":contrast(x,y), "predicted_behavior":behavior(V)
        })
        if name in key:
            plt.plot(x,y,label=f"{label} | V={V:.2f}")
        plt_i=plt.figure(figsize=(10,5))
        plt.plot(x,y,label=f"pattern | V={V:.3f}")
        plt.plot(x,env,"--",label="envelope")
        plt.title(f"Q6: {label}")
        plt.xlabel("screen position x"); plt.ylabel("normalized probability")
        plt.legend(); plt.tight_layout()
        plt.savefig(outdir/f"Q6_{name}.png", dpi=160)
        plt.close(plt_i)
    plt.title("Q6 A-information receipt test: existence known is not path-A known")
    plt.xlabel("screen position x"); plt.ylabel("normalized probability")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(outdir/"Q6_key_cases_overlay.png", dpi=160); plt.close()

    df=pd.DataFrame(rows); df.to_csv(outdir/"Q6_case_summary.csv", index=False)

    vals=np.linspace(0,1,101)
    grid=[]
    for P in vals:
        for Env in vals:
            V,R_A,_=visibility_from_A_information(1,P,Env)
            grid.append({"path_receipt_P":P,"environment_receipt_Env":Env,"R_A":R_A,"visibility_V":V})
    g=pd.DataFrame(grid); g.to_csv(outdir/"Q6_path_environment_receipt_grid.csv", index=False)
    pivot=g.pivot(index="environment_receipt_Env", columns="path_receipt_P", values="visibility_V")
    plt.figure(figsize=(7,6))
    plt.imshow(pivot.values, origin="lower", aspect="auto", extent=[0,1,0,1])
    plt.colorbar(label="interference visibility V")
    plt.xlabel("path-A receipt strength"); plt.ylabel("environment receipt strength")
    plt.title("Any durable receipt resolves A-information")
    plt.tight_layout(); plt.savefig(outdir/"Q6_visibility_grid_path_vs_environment.png", dpi=160); plt.close()

    e_rows=[]
    for E in vals:
        V,R_A,_=visibility_from_A_information(E,0,0)
        e_rows.append({"existence_receipt_E":E,"R_A":R_A,"visibility_V":V})
    e=pd.DataFrame(e_rows); e.to_csv(outdir/"Q6_existence_receipt_sweep.csv", index=False)
    plt.figure(figsize=(8,5))
    plt.plot(e["existence_receipt_E"], e["visibility_V"], label="visibility")
    plt.plot(e["existence_receipt_E"], e["R_A"], label="resolved A-information")
    plt.title("Existence receipt alone does not resolve path-A")
    plt.xlabel("existence receipt strength"); plt.ylabel("value")
    plt.legend(); plt.tight_layout(); plt.savefig(outdir/"Q6_existence_receipt_sweep.png", dpi=160); plt.close()

    summary={
        "claim_tested":"Quantum-like behavior is tied to unresolved A-information/path receipt, not merely unknown existence.",
        "main_result":[
            "Known existence alone leaves R_A=0 and interference survives.",
            "Durable path receipt resolves A-information and suppresses interference.",
            "Environment/influence can resolve A without human measurement.",
            "Erased/non-durable receipt restores unresolved behavior in the toy.",
            "The key bucket is A-information known vs unknown, not particle existence known vs unknown."
        ],
        "danger_cases":[
            "Durable path-specific receipt exists but full interference remains.",
            "No path or environment receipt exists anywhere but interference disappears permanently."
        ],
        "barstool_read":"Knowing a particle exists is not the same as the universe keeping the path receipt."
    }
    (outdir/"Q6_summary.json").write_text(json.dumps(summary,indent=2))
    print("Q6 A-information receipt test complete.")
    print(f"Output directory: {outdir}")
    print(df[["case","existence_receipt_E","path_receipt_P","environment_receipt_Env","erase_path","erase_environment","resolved_A_information_R_A","visibility_V","central_contrast","predicted_behavior"]].to_string(index=False))
    print("\nMain result: existence known is not path-A known; durable path/environment receipts resolve A-information.")
if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--outdir", default="results/Q6_A_information_receipt_test")
    run(p.parse_args().outdir)
