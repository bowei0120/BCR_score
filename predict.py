#!/usr/bin/env python3
"""Inference side of the BCR score: predict.py is a thin wrapper that reuses
bcr_core.coverage_table + bcr_core.compute_scales so training (estimate_params.py)
and inference share one code path; required inputs are the precomputed
fisher_sets.json / params_full.json."""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import bcr_core

CHAINS = list(bcr_core.CHAINS)
GROUPS = list(bcr_core.GROUPS)
OUT_COLS = [
    "sample", "scale_UC", "scale_CD",
    "scale_IGL_UC", "scale_IGL_CD",
    "scale_IGK_UC", "scale_IGK_CD",
    "n_IGL", "k_IGL_UC", "k_IGL_CD",
    "n_IGK", "k_IGK_UC", "k_IGK_CD",
]


def load_params(path):
    """Load per-chain Beta-Binomial model parameters (key in HT_UC/UC/HT_CD/CD)."""
    with open(path) as f:
        return json.load(f)


def load_fisher_json(path):
    """Load per-chain universe/UC/CD CDR3aa reference sets (sorted lists -> sets)."""
    with open(path) as f:
        payload = json.load(f)
    return {chain: {key: set(vals) for key, vals in payload[chain].items()}
            for chain in CHAINS}


def main():
    parser = argparse.ArgumentParser(
        description="Compute disease BCR scores (UC/CD) and covered sequences for "
                    "new samples from a sample table (columns: sample, clono_path)")
    parser.add_argument("--table", type=str, required=True,
                        help="TSV with columns 'sample' and 'clono_path'")
    parser.add_argument("--fisher-json", type=str, required=True,
                        help="Precomputed CDR3 reference sets JSON")
    parser.add_argument("--params", type=str, required=True,
                        help="Pre-trained Beta-Binomial params JSON")
    parser.add_argument("--out", type=str, default=None,
                        help="Write TSV to this file (default: stdout)")
    args = parser.parse_args()

    table = pd.read_table(args.table, comment="#")
    for col in ("sample", "clono_path"):
        if col not in table.columns:
            print("ERROR: --table must have a '{}' column".format(col), file=sys.stderr)
            sys.exit(2)

    cohort = table[["sample", "clono_path"]].copy()
    cohort["patient"] = ""
    cohort["label"] = ""

    refs = load_fisher_json(Path(args.fisher_json))
    params = load_params(args.params)
    cov = bcr_core.coverage_table(cohort, refs)
    scored = bcr_core.compute_scales(cov, params)

    rename = {"{}_universe".format(c): "n_{}".format(c) for c in CHAINS}
    for c in CHAINS:
        for g in GROUPS:
            rename["{}_{}".format(c, g)] = "k_{}_{}".format(c, g)
    out_df = scored.rename(columns=rename)[OUT_COLS]
    if args.out:
        out_df.to_csv(args.out, sep="\t", index=None)
    else:
        sys.stdout.write(out_df.to_csv(sep="\t", index=None))


if __name__ == "__main__":
    main()
