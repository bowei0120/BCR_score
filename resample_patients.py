"""Patient-level resampling of a cohort table to one specimen per patient.

Cohort tables list biopsy specimens; some patients contributed multiple
specimens, so downstream patient-level analyses would violate the
independence assumption if they treated every specimen as an independent
observation.  This script addresses that within-patient non-independence by
resampling the cohort N_ITER times: in each iteration every patient
contributes exactly ONE of their specimens, chosen uniformly at random
(single-specimen patients trivially always keep theirs).  The resulting
iter_XXX.tsv files have the exact same schema as the input table and can be
fed back to other scripts as --table inputs.

Iteration i is drawn with numpy default_rng(SEED + i), so runs are
deterministic and iterations are independent of each other.

Usage:
  python3 resample_patients.py --table COHORT.tsv --out-dir OUT
  python3 resample_patients.py --table COHORT.tsv --out-dir OUT --n-iter 1000 --seed 0
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT_COLS = ["sample", "patient", "label", "clono_path"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resample a cohort table to one random specimen per patient")
    parser.add_argument("--table", type=str, required=True,
                        help="TSV with columns sample, patient, label, clono_path")
    parser.add_argument("--out-dir", type=str, required=True,
                        help="output directory for iter_XXX.tsv files")
    parser.add_argument("--n-iter", type=int, default=1000,
                        help="number of resampling iterations (default: 1000)")
    parser.add_argument("--seed", type=int, default=0,
                        help="base seed; iteration i uses seed+i (default: 0)")
    return parser.parse_args()


def load_cohort(path):
    """Load the cohort TSV and validate its four required columns."""
    table = pd.read_table(path, dtype=str)
    for col in OUT_COLS:
        if col not in table.columns:
            print("ERROR: --table must have a '{}' column".format(col),
                  file=sys.stderr)
            sys.exit(2)
    return table[OUT_COLS]


def main():
    args = parse_args()
    table = load_cohort(args.table)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    by_sample = {row["sample"]: row
                 for _, row in table.iterrows()}
    groups = {pt: sorted(g["sample"].tolist())
              for pt, g in table.groupby("patient")}
    patients = sorted(groups)
    n_patients = len(patients)
    multi = sum(1 for pt in patients if len(groups[pt]) > 1)
    print("[resample_patients] patients={} (multi-sample={}) samples={} "
          "iters={} seed={}".format(n_patients, multi, len(table),
                                    args.n_iter, args.seed))

    for i in range(args.n_iter):
        rng = np.random.default_rng(args.seed + i)
        rows = []
        for pt in patients:
            samples = groups[pt]
            chosen = samples[0] if len(samples) == 1 \
                else samples[int(rng.integers(len(samples)))]
            rows.append(by_sample[chosen])
        df = pd.DataFrame(rows, columns=OUT_COLS)
        assert len(df) == n_patients, (i, len(df), n_patients)
        df = df.sort_values(["patient", "sample"], kind="mergesort")
        df.to_csv(out_dir / "iter_{:03d}.tsv".format(i),
                  sep="\t", index=None)
        if (i + 1) % 100 == 0:
            print("[resample_patients] {}/{} iterations written".format(
                i + 1, args.n_iter))

    print("[resample_patients] done: {} patients, {} samples per iteration, "
          "{} files -> {}".format(n_patients, n_patients, args.n_iter,
                                  out_dir))


if __name__ == "__main__":
    main()
