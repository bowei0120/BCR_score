#!/usr/bin/env python3
"""Estimate Beta-Binomial parameters and in-sample F scales for a BCR cohort.

Pipeline (training side of predict.py):

  1. load the per-chain cover/UC/CD CDR3aa sets (fisher_sets.json, the exact
     schema written by feature_selection.py / consumed by predict.py),
  2. compute per-sample coverage features (coverage.tsv),
  3. fit the Beta-Binomial disease model per chain x key in
     {HT_UC, UC, HT_CD, CD} (params.json, same schema as params_full.json),
  4. score every sample with the log Bayes factor F (scales.tsv).

Usage:
    python3 estimate_params.py --table COHORT.tsv --fisher-json PATH --out-dir OUT
"""
import argparse
import json
from pathlib import Path

import pandas as pd

import bcr_core
from predict import load_fisher_json

BASE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(
        description="Fit Beta-Binomial params (params_full.json schema) and "
                    "in-sample F scales from the Fisher coverage features")
    parser.add_argument("--table", type=str, required=True,
                        help="Cohort TSV (columns: sample, patient, label, clono_path)")
    parser.add_argument("--fisher-json", type=str,
                        default=str(BASE / "fisher_sets.json"),
                        help="Per-chain cover/UC/CD CDR3 reference sets JSON "
                             "(default: fisher_sets.json next to this script)")
    parser.add_argument("--out-dir", type=str, required=True,
                        help="Output dir (coverage.tsv, params.json, scales.tsv)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort = bcr_core.load_cohort(args.table)
    refs = load_fisher_json(Path(args.fisher_json))

    cov = bcr_core.coverage_table(cohort, refs)
    cov.to_csv(out_dir / "coverage.tsv", sep="\t", index=None)
    print("wrote {} ({} samples)".format(out_dir / "coverage.tsv", len(cov)))

    params = bcr_core.fit_params(cov)
    with open(out_dir / "params.json", "w") as f:
        json.dump(params, f, indent=2)
    rows = [{"chain": chain, "key": key, "a": p["a"], "b": p["b"], "N": p["N"]}
            for chain in params for key, p in params[chain].items()]
    print(pd.DataFrame(rows).to_string(index=False))

    scales = bcr_core.compute_scales(cov, params)
    scale_cols = ["sample", "patient", "label"]
    for chain in bcr_core.CHAINS:
        for group in bcr_core.GROUPS:
            scale_cols.append("scale_{}_{}".format(chain, group))
    scale_cols += ["scale_UC", "scale_CD"]
    scales[scale_cols].to_csv(out_dir / "scales.tsv", sep="\t", index=None)
    print("wrote {} ({} samples)".format(out_dir / "scales.tsv", len(scales)))


if __name__ == "__main__":
    main()
