#!/usr/bin/env python3
"""
G113/G114 frequency-domain QNM solver.

G113a: Schwarzschild calibration.
G113b: generalized frequency-domain Riccati log-derivative solver.
G114: STAM V_exact roots.

This is NOT a Leaver recurrence. It is a generalized frequency-domain
Riccati matching solver for arbitrary numerical potentials:

    psi'' + [omega^2 - V(x)] psi = 0, x=r*
    y = psi'/psi
    y' = V - omega^2 - y^2

Boundary:
    left/horizon:  y=-i omega
    right/infty:   y=+i omega

Root:
    y_left(x_match) - y_right(x_match) = 0
"""

from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import root
from scipy.special import lambertw

M=1.0
LEAVER={(2,0):0.373672-0.088962j,(2,1):0.346711-0.273915j,(3,0):0.599443-0.092703j,(4,0):0.809178-0.094164j}

def h(r): return 1-2/r
def hp(r): return 2/r**2
def y_of_r(r): return 6/r-2
def F(y): return 1-5*y**4+4*y**5
def Fp(y): return -20*y**3+20*y**4
def k_stam(r):
    y=y_of_r(r)
    return h(r)*F(y) if y>0 else h(r)
def kp_stam(r):
    y=y_of_r(r)
    return hp(r)*F(y)+h(r)*Fp(y)*(-6/r**2) if y>0 else hp(r)

def r_of_x_GR(x):
    off=3+2*np.log(0.5)
    R=x+off
    u=np.real(lambertw(np.exp(R/2-1)))
    return 2*(1+u)

def build_r_stam(xgrid):
    xgrid=np.asarray(xgrid)
    r=np.zeros_like(xgrid)
    def rhs(x,y): return [np.sqrt(max(h(y[0])*k_stam(y[0]),0))]
    pos=xgrid>=0
    if np.any(pos):
        xs=np.sort(xgrid[pos])
        sol=solve_ivp(rhs,(0,float(xs[-1])),[3.0],t_eval=xs,rtol=1e-8,atol=1e-10,max_step=1.0)
        r[pos]=np.interp(xgrid[pos],xs,sol.y[0])
    neg=xgrid<0
    if np.any(neg):
        xs=np.sort(xgrid[neg])[::-1]
        sol=solve_ivp(rhs,(0,float(xs[-1])),[3.0],t_eval=xs,rtol=1e-8,atol=1e-10,max_step=0.1)
        # interp wants ascending
        r[neg]=np.interp(xgrid[neg],xs[::-1],sol.y[0][::-1])
    return r

def V_GR_r(r,ell):
    return h(r)*(ell*(ell+1)/r**2-6/r**3)

def V_geom_r(r,ell):
    hv=h(r); kv=k_stam(r)
    return hv*(ell*(ell+1)/r**2 - 2*(1-kv)/r**2 - kp_stam(r)/(2*r) - kv*hp(r)/(2*hv*r))

def dlnf_dA(A):
    if A<=2/3: return 0.0
    y=3*A-2; FF=F(y)
    return -2*y**3*(45*A*A-33*A-13)/(A*FF)

def d2lnf_dA2(A):
    hh=1e-5
    if A-hh<=2/3:
        return (dlnf_dA(A+hh)-dlnf_dA(A))/hh
    if A+hh>=0.999999:
        hh=max(1e-8,(0.999999-A)/2)
    return (dlnf_dA(A+hh)-dlnf_dA(A-hh))/(2*hh)

def corr_exact_r(r):
    if r>=3: return 0.0
    A=2/r
    dA=-2/r**2; d2A=4/r**3
    u=dlnf_dA(A); up=d2lnf_dA2(A)
    dlnf_dr=u*dA
    d_dlnf_dr=up*dA*dA+u*d2A
    hk=h(r)*k_stam(r)
    if hk<=0: return 0.0
    sq=np.sqrt(hk)
    dhk=hp(r)*k_stam(r)+h(r)*kp_stam(r)
    P=sq*dlnf_dr
    dPdr=0.5*dhk/sq*dlnf_dr+sq*d_dlnf_dr
    Pstar=sq*dPdr
    return 0.5*Pstar+0.25*P*P

def build_potential(ell,kind,xmin,xmax,N):
    x=np.linspace(xmin,xmax,N)
    if kind=="GR":
        r=r_of_x_GR(x); V=V_GR_r(r,ell)
    else:
        r=build_r_stam(x)
        V=np.array([V_geom_r(rv,ell)+(corr_exact_r(rv) if kind=="STAM_exact" else 0.0) for rv in r])
    return x,V

def mismatch(omega,x,V,xm=0.0):
    Vf=interp1d(x,V,kind="linear",fill_value=(V[0],V[-1]),bounds_error=False)
    def rhs(xv,yv):
        yy=yv[0]+1j*yv[1]
        dy=float(Vf(xv))-omega*omega-yy*yy
        return [dy.real,dy.imag]
    y0=-1j*omega
    solL=solve_ivp(rhs,(float(x[0]),xm),[y0.real,y0.imag],rtol=1e-7,atol=1e-9,max_step=0.5)
    if not solL.success: return 1e9+1e9j
    yL=solL.y[0,-1]+1j*solL.y[1,-1]
    y0=1j*omega
    solR=solve_ivp(rhs,(float(x[-1]),xm),[y0.real,y0.imag],rtol=1e-7,atol=1e-9,max_step=0.5)
    if not solR.success: return 1e9+1e9j
    yR=solR.y[0,-1]+1j*solR.y[1,-1]
    return yL-yR

def solve_root(ell,n,kind,guess,xmin,xmax,N):
    x,V=build_potential(ell,kind,xmin,xmax,N)
    def fun(z):
        m=mismatch(z[0]+1j*z[1],x,V)
        return [m.real,m.imag]
    sol=root(fun,[guess.real,guess.imag],method="hybr",tol=1e-8)
    w=sol.x[0]+1j*sol.x[1]
    res=abs(mismatch(w,x,V))
    return w,sol.success,res

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",default="/mnt/data/G113_G114_run")
    ap.add_argument("--xmin",type=float,default=-80)
    ap.add_argument("--xmax",type=float,default=180)
    ap.add_argument("--N",type=int,default=2500)
    args=ap.parse_args()
    out=Path(args.outdir); resdir=out/"results"; resdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    print("="*80)
    print("G113/G114 frequency-domain Riccati QNM solver")
    print("="*80)
    print("Method: Riccati matching, not Leaver recurrence.\n")
    for ell,n in [(2,0),(2,1),(3,0),(4,0)]:
        gold=LEAVER[(ell,n)]
        w,ok,res=solve_root(ell,n,"GR",gold,args.xmin,args.xmax,args.N)
        err=abs(w-gold)/abs(gold)
        rows.append(["GR",ell,n,w.real,w.imag,gold.real,gold.imag,err,ok,res])
        print(f"G113a GR l={ell} n={n}: omega={w.real:.9f}{w.imag:+.9f}i gold={gold.real:.9f}{gold.imag:+.9f}i err={err*100:.3f}% ok={ok} resid={res:.2e}")
    print()
    for ell,n in [(2,0),(3,0),(4,0)]:
        gold=LEAVER[(ell,n)]
        # use time-domain expected shift as rough seed for l=2, smaller for others
        seed=gold*(1+0.03)
        w,ok,res=solve_root(ell,n,"STAM_exact",seed,args.xmin,args.xmax,args.N)
        shift=abs(w-gold)/abs(gold)
        rows.append(["STAM_exact",ell,n,w.real,w.imag,gold.real,gold.imag,shift,ok,res])
        print(f"G114 STAM l={ell} n={n}: omega={w.real:.9f}{w.imag:+.9f}i shift_vs_GRgold={shift*100:.3f}% ok={ok} resid={res:.2e}")
    csvp=resdir/"G113_G114_frequency_domain_roots.csv"
    with csvp.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["kind","ell","n","omega_re","omega_im","ref_re","ref_im","rel_err_or_shift","success","residual_abs"])
        w.writerows(rows)
    md=resdir/"G113_G114_frequency_domain_summary.md"
    lines=["# G113/G114 frequency-domain QNM solver\n\n","Method: Riccati log-derivative matching on numerical potential. This is not a true Leaver recurrence.\n\n","| kind | ell | n | omega | reference | err/shift |\n","|---|---:|---:|---:|---:|---:|\n"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]:.9f}{r[4]:+.9f}i | {r[5]:.9f}{r[6]:+.9f}i | {100*r[7]:.3f}% |\n")
    lines.append("\nCalibration decides whether STAM roots are prediction-grade. Poor n=1 calibration means n=1 remains open.\n")
    md.write_text("".join(lines),encoding="utf-8")
    print(f"\nCSV written: {csvp}")
    print(f"Summary written: {md}")

if __name__=="__main__":
    main()
