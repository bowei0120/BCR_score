#!/usr/bin/env python3
"""Leakage-free leave-one-sample-out cross-validation of the BCR score.

For every held-out sample both the Fisher feature selection and the
Beta-Binomial parameter estimation are re-run on the remaining samples ONLY.
Each sample is therefore scored by a model that never saw it.

Pipeline per held-out sample u:
  1. train = cohort rows whose sample != u
  2. fisher_tables(train)        [feature selection excludes u]
  3. build_sets per chain
  4. coverage_table(train)       [train coverage only]
  5. fit_params(train coverage)  [params exclude u]
  6. coverage_table(u rows) against the train sets
  7. compute_scales(u rows, train params)
  8. collect the scored rows with held_out = u

Usage:
    python3 cross_validate.py --table COHORT.tsv --out-dir OUT \
        [--chains IGL IGK] [--p 0.01] [--or-thresh 1.5] [--min-count 5]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import bcr_core

try:
    _trapz = np.trapezoid
except AttributeError:  # numpy < 2.0
    _trapz = np.trapz


def auc_roc(scores, positives):
    """Trapezoidal ROC AUC (no sklearn): sort desc, cumulative TPR/FPR.

    Returns (auc, fpr, tpr) with a (0, 0) origin prepended, or
    (nan, None, None) when either class is absent from the panel.
    """
    scores = np.asarray(scores, dtype=float)
    pos = np.asarray(positives, dtype=bool)
    n_pos = int(pos.sum())
    n_neg = int((~pos).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), None, None
    order = np.argsort(-scores, kind="mergesort")
    ranked = scores[order]
    p = pos[order]
    ends = np.append(np.where(np.diff(ranked) != 0)[0], ranked.size - 1)
    tps = np.cumsum(p)[ends]
    fps = np.cumsum(~p)[ends]
    tpr = np.concatenate(([0.0], tps / n_pos))
    fpr = np.concatenate(([0.0], fps / n_neg))
    return float(_trapz(tpr, fpr)), fpr, tpr


def warn_missing_groups(train):
    """Warn (stderr) about labels absent from a training set."""
    present = set(train["label"])
    for group in bcr_core.GROUPS:
        missing = [lbl for lbl in ("HT", group) if lbl not in present]
        if missing:
            print("WARNING: training set has no {} rows -> {} model "
                  "degenerate".format("+".join(missing), group),
                  file=sys.stderr)


def fit_and_score(train, test, args):
    """Fit the model on train only, then score the independent test rows.

    Feature selection (fisher_tables/build_sets) and parameter estimation
    (fit_params) both see only `train`, so the returned test rows are
    out-of-sample. Scoring options (chains, thresholds) come from args.
    Returns the scored test rows (coverage + scale columns).
    """
    tables = bcr_core.fisher_tables(train, chains=args.chains)
    sets_by_chain = {}
    for chain in args.chains:
        sets_by_chain[chain] = bcr_core.build_sets(
            tables[chain], p_thresh=args.p, or_thresh=args.or_thresh,
            min_count=args.min_count)
    warn_missing_groups(train)
    train_cov = bcr_core.coverage_table(train, sets_by_chain, chains=args.chains)
    params = bcr_core.fit_params(train_cov, chains=args.chains)
    test_cov = bcr_core.coverage_table(test, sets_by_chain, chains=args.chains)
    return bcr_core.compute_scales(test_cov, params, chains=args.chains)


def build_jobs(cohort):
    """List of (held_out_sample, train_df, test_df): one per held-out sample."""
    jobs = []
    for uid in pd.unique(cohort["sample"]):
        jobs.append((uid, cohort[cohort["sample"] != uid],
                     cohort[cohort["sample"] == uid]))
    return jobs


def collect_auc(result, groups):
    """One out-of-sample AUC row per group over all held-out samples."""
    auc_rows, roc_rows = [], []
    for group in groups:
        col = "scale_{}".format(group)
        sub = result[result["label"].isin(["HT", group])].dropna(subset=[col])
        n_pos = int((sub["label"] == group).sum())
        n_neg = int((sub["label"] == "HT").sum())
        auc, fpr, tpr = auc_roc(sub[col].to_numpy(),
                                (sub["label"] == group).to_numpy())
        if fpr is None:
            print("WARNING: {} vs HT: one class absent -> AUC NaN".format(group),
                  file=sys.stderr)
        auc_rows.append({"group": group, "n_pos": n_pos, "n_neg": n_neg,
                         "auc": auc})
        if fpr is not None:
            for fp, tp in zip(fpr, tpr):
                roc_rows.append({"group": group, "fpr": fp, "tpr": tp})
    return auc_rows, roc_rows


def main():
    parser = argparse.ArgumentParser(
        description="Leakage-free leave-one-sample-out CV of the BCR "
                    "Beta-Binomial score (Fisher selection and BB params "
                    "refit excluding every held-out sample)")
    parser.add_argument("--table", type=str, required=True,
                        help="Cohort TSV (columns: sample, patient, label, clono_path)")
    parser.add_argument("--out-dir", type=str, required=True,
                        help="Output dir (loo_scores.tsv, auc.tsv, roc.tsv)")
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
    if "sample" not in cohort.columns:
        print("ERROR: cohort table must have a 'sample' column", file=sys.stderr)
        sys.exit(2)

    jobs = build_jobs(cohort)
    print("cross_validate: {} held-out sample(s), chains={}"
          .format(len(jobs), " ".join(args.chains)), flush=True)

    scored_frames = []
    for held, train, test in jobs:
        scored = fit_and_score(train, test, args)
        scored["held_out"] = held
        scored_frames.append(scored)

    if not scored_frames:
        print("ERROR: no held-out samples -> nothing to score", file=sys.stderr)
        sys.exit(2)
    result = pd.concat(scored_frames, ignore_index=True)

    scale_cols = []
    for chain in args.chains:
        for group in bcr_core.GROUPS:
            scale_cols.append("scale_{}_{}".format(chain, group))
    scale_cols += ["scale_UC", "scale_CD"]
    out_cols = ["sample", "patient", "label", "held_out"] + scale_cols
    result[out_cols].to_csv(out_dir / "loo_scores.tsv", sep="\t", index=None)
    print("wrote {} ({} samples)".format(out_dir / "loo_scores.tsv", len(result)))

    auc_rows, roc_rows = collect_auc(result, bcr_core.GROUPS)
    pd.DataFrame(auc_rows).to_csv(out_dir / "auc.tsv", sep="\t", index=None)
    pd.DataFrame(roc_rows).to_csv(out_dir / "roc.tsv", sep="\t", index=None)

    print("AUC (out-of-sample):")
    for row in auc_rows:
        print("  {} vs HT: n_pos={} n_neg={} auc={:.3f}".format(
            row["group"], row["n_pos"], row["n_neg"], row["auc"]))
    print("wrote {} and {}".format(out_dir / "auc.tsv", out_dir / "roc.tsv"))


if __name__ == "__main__":
    main()
