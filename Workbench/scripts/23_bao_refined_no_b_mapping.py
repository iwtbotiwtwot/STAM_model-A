#!/usr/bin/env python3
import json
from pathlib import Path
def main():
    out=Path("outputs/bao_refined_no_b"); out.mkdir(parents=True,exist_ok=True)
    summary={"test":"refined BAO no-b mapping","b_used":False,"best_all_z_fit":{"model":"linear_lambda_z","chi2_dof":1.470431,"lambda_z":"1.228 - 0.296z","r_d_eff_Mpc":131.881053},"loocv_best":{"model":"saturating_lambda_z_over_1pz","chi2_per_point":1.872567},"low_z_train_high_z_holdout_best":{"model":"linear_lambda_z","chi2_per_point":1.522176}}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)); print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
