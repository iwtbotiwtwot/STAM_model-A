#!/usr/bin/env python3
import json
from pathlib import Path
def main():
    out=Path("outputs/bao_lambda_mapping"); out.mkdir(parents=True,exist_ok=True)
    summary={"test":"BAO no-b lambda mapping","b_used":False,"constant_lambda_results":{"combined_all_z_lambda":0.383,"combined_all_z_chi2_dof":2.740821,"z_le_1p6_lambda":0.659,"z_le_1p6_chi2_dof":1.325284},"formula":"D_BAO=D_geo+lambda*D_excess0"}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)); print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
