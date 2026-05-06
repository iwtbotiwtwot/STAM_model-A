#!/usr/bin/env python3
import json
from pathlib import Path
def main():
    out=Path("outputs/lambda_cross_observable"); out.mkdir(parents=True,exist_ok=True)
    summary={"test":"cross-observable lambda sniff SN vs BAO","b_used":False,"bao_constant_lambda":0.3832429338743046,"bao_linear_lambda":"1.2276269984933077 - 0.29632959757613064z","best_bao_pointwise_shape":"linear_z beta0=1.241546 beta1=-0.308341","sn_result":"shape-calibrated SN lambda structured but catalog-scale sensitive; no single shared curve established yet"}
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)); print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
