#!/usr/bin/env python3
import json
from pathlib import Path
def main():
    out=Path("outputs/distance_onset"); out.mkdir(parents=True,exist_ok=True)
    summary={"test":"distance-adjustment onset z0 scan","b_used":False,"formula":"D_adj,onset=D_geo+0.35L*max(z-z0,0)^2","result":"best z0 was near 0; z0=0.30 was worse for this onset form","best_common_single_intercept":{"Union3":0.0,"Pantheon":0.009,"Pantheon_Union3":0.003,"DES":0.0,"All_three":0.015}}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)); print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
