#!/usr/bin/env python3
"""Generate N random-label (permuted) files per group for permutation testing.

Input: a cohort TSV with columns  sample, patient, label, clono_path.

Domain = rows whose clono_path is non-empty AND exists on disk.

Groups (fixed):
  - CD = {CD, HT} -> random_CD/
  - UC = {UC, HT} -> random_UC/

Output files keep the SAME columns (and order) as the input table; only the
label column is shuffled, and ONLY WITHIN the domain so the case:control
ratio of domain samples stays exactly unchanged. Non-domain rows keep
original labels. A fixed RNG seed makes the whole set reproducible.

Usage:
    python3 generate_random_labels.py --table cohort.tsv \
        [--out-base-dir .] [--n-files 200] [--seed 42]
"""
import argparse
import random
import sys
from pathlib import Path

import pandas as pd

GROUPS = [("CD", {"CD", "HT"}), ("UC", {"UC", "HT"})]


def is_domain(clono_path) -> bool:
    if not isinstance(clono_path, str) or not clono_path.strip():
        return False
    return Path(clono_path).is_file()


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate random-label files.")
    ap.add_argument("--table", required=True, help="Cohort TSV")
    ap.add_argument("--out-base-dir", default=".",
                    help="Where random_<PREFIX>/ dirs are written (default: .)")
    ap.add_argument("--n-files", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    df = pd.read_csv(args.table, sep="\t", dtype=str).fillna("")
    for col in ("sample", "patient", "label", "clono_path"):
        if col not in df.columns:
            sys.exit(f"ERROR: missing column '{col}' in {args.table}")

    df["in_domain"] = df["clono_path"].map(is_domain)
    out_base = Path(args.out_base_dir)
    out_base.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)

    for prefix, allowed in GROUPS:
        sub = df[df["label"].isin(allowed)].copy()
        if sub.empty:
            print(f"WARNING: no rows for {prefix} (labels={sorted(allowed)}); skip",
                  file=sys.stderr)
            continue

        dom_mask = sub["in_domain"].to_numpy()
        dom_labels = sub.loc[dom_mask, "label"].to_numpy()
        n_dom = len(dom_labels)
        dom_counts = pd.Series(dom_labels).value_counts().to_dict()

        outdir = out_base / f"random_{prefix}"
        outdir.mkdir(exist_ok=True)

        # 除 label 外，其余列均原样保留（含列顺序）；写出时去掉辅助列
        base = sub.drop(columns=["in_domain"])
        columns = list(base.columns)
        label_idx = columns.index("label")

        for i in range(1, args.n_files + 1):
            out_df = base.copy()
            if n_dom:
                shuffled = dom_labels.copy()
                rng.shuffle(shuffled)
                # 只改 domain 行的 label
                labels = out_df["label"].to_numpy()
                labels[dom_mask] = shuffled
                out_df["label"] = labels
            out_path = outdir / f"{prefix}_{i:03d}.txt"
            out_df.to_csv(out_path, sep="\t", index=False)

        print(f"{outdir.name}: {args.n_files} files, {len(sub)} samples each, "
              f"domain={n_dom} {dom_counts} - labels shuffled within domain only")
    return 0


if __name__ == "__main__":
    sys.exit(main())