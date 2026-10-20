#!/usr/bin/env python3
"""Fisher exact test CDR3 feature selection on a BCR cohort.

For each chain (default IGL/IGK) builds one Fisher exact table over the
cohort (carrier counts per sample, disease-oriented OR vs HT), writes the
full table, then thresholds it into the per-chain cover/UC/CD CDR3aa sets
and writes them as a single JSON in the exact schema consumed by predict.py
load_fisher_json:

    {"IGL": {"cover": [..], "UC": [..], "CD": [..]}, "IGK": {...}}

Usage:
    python3 feature_selection.py --table COHORT.tsv --out-dir OUT \
        [--chains IGL IGK] [--p 0.01] [--or-thresh 1.5] [--min-count 5]
"""
import argparse
import json
from pathlib import Path

import bcr_core


def main():
    parser = argparse.ArgumentParser(
        description="Fisher exact test feature selection: per-chain CDR3 "
                    "carrier tables -> cover/UC/CD sets JSON for predict.py")
    parser.add_argument("--table", type=str, required=True,
                        help="Cohort TSV (columns: sample, patient, label, clono_path)")
    parser.add_argument("--out-dir", type=str, required=True,
                        help="Output dir ({chain}_fisher.tsv + fisher_sets.json)")
    parser.add_argument("--chains", type=str, nargs="+", default=list(bcr_core.CHAINS),
                        help="Chains to analyze (default: IGL IGK)")
    parser.add_argument("--p", type=float, default=0.01,
                        help="Max two-sided p-value for the UC/CD sets (default: 0.01)")
    parser.add_argument("--or-thresh", type=float, default=1.5,
                        help="Min disease-oriented odds ratio (default: 1.5)")
    parser.add_argument("--min-count", type=int, default=5,
                        help="Min carrier count in the disease group (default: 5)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort = bcr_core.load_cohort(args.table)
    tables = bcr_core.fisher_tables(cohort, chains=args.chains)

    payload = {}
    for chain, df in tables.items():
        tsv = out_dir / "{}_fisher.tsv".format(chain)
        df.to_csv(tsv, sep="\t", index=None)
        print("wrote {} ({} CDR3 rows)".format(tsv, len(df)))
        sets = bcr_core.build_sets(df, p_thresh=args.p, or_thresh=args.or_thresh,
                                   min_count=args.min_count)
        payload[chain] = {key: sorted(values) for key, values in sets.items()}
        print("{}: universe={} |UC|={} |CD|={}".format(
            chain, len(sets["cover"]), len(sets["UC"]), len(sets["CD"])))

    json_path = out_dir / "fisher_sets.json"
    with open(json_path, "w") as f:
        json.dump(payload, f)
    print("wrote {} ({:.1f} MB)".format(
        json_path, json_path.stat().st_size / 1e6))


if __name__ == "__main__":
    main()
