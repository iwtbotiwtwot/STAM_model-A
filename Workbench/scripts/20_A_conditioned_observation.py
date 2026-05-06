#!/usr/bin/env python3
import json
from pathlib import Path
C0=299792458.0
D=10_000_000_000.0*1000
def S(A): return 1+A
def frame(A):
    s=S(A); Da=s*D; T=Da/C0
    return {"A":A,"S":s,"local_c_over_c0":1.0,"naive_external_c_over_c0":1/s,"proper_accumulated_c_over_c0":1.0,"round_trip_cross_s":2*T,"round_trip_local_s":2*T/s}
def main():
    out=Path("outputs/A_conditioned_observation"); out.mkdir(parents=True,exist_ok=True)
    rows=[frame(A) for A in [-0.5,0,0.001,0.1,0.99]]
    summary={"test":"A-conditioned observation","A_alpha":0.001,"A_epsilon":0.1,"period_ratio_epsilon_over_alpha":S(0.1)/S(0.001),"tick_rate_ratio":S(0.001)/S(0.1),"local_c_check":1.0,"rows":rows}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)); print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
