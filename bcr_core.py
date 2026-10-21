#!/usr/bin/env python3
"""Shared library for BCR-score training: cohort IO, Fisher exact feature
selection, coverage features, and Beta-Binomial model fitting/scoring.

Semantics stay consistent with the inference side (predict.py) and the
original src_fold utilities:

  - chain membership via mixcr_type/trust_type (';'-separated gene types)
  - coverage features count ALL unique CDR3aa of a sample (chain-unfiltered,
    matching predict.py sample_features)
  - Fisher tables count carriers per sample (sample-level) and use the
    disease-oriented contingency [[k_g, n_g-k_g], [k_HT, n_HT-k_HT]] so
    OR > 1 means enriched in disease
  - score_F is the log Bayes factor used by predict.py (identical formula)
"""
import math
import sys

import numpy as np
import pandas as pd
from scipy.special import betaln, digamma
from scipy.stats import fisher_exact

CHAINS = ("IGL", "IGK")
GROUPS = ("UC", "CD")
LABELS = ("HT", "UC", "CD")
KEYS = ("HT_UC", "UC", "HT_CD", "CD")
LABEL_OF = {"HT_UC": "HT", "UC": "UC", "HT_CD": "HT", "CD": "CD"}
COMPARE_OF = {"HT_UC": "UC", "UC": "UC", "HT_CD": "CD", "CD": "CD"}
B_GRID = np.arange(1, 3001, 100, dtype=float)


def load_cohort(table_path):
    """Read the cohort TSV (columns sample, patient, label, clono_path).

    Prints ERROR to stderr and exits with code 2 if a required column is
    missing.
    """
    cohort = pd.read_table(table_path, comment="#")
    for col in ("sample", "patient", "label", "clono_path"):
        if col not in cohort.columns:
            print("ERROR: cohort table must have a '{}' column".format(col),
                  file=sys.stderr)
            sys.exit(2)
    return cohort


def chain_aa_by_sample(clono_path, chains=CHAINS):
    """Read one clonotype table -> {chain: set(CDR3aa)}.

    A row belongs to a chain if that chain string appears in mixcr_type OR
    trust_type (values like 'IGH' or 'IGK;IGL'; membership after split on
    ';'). Unique non-NaN CDR3aa only.
    """
    out = {chain: set() for chain in chains}
    df = pd.read_table(clono_path)
    if df.empty:
        return out
    df = df.dropna(subset=["CDR3aa"])
    if df.empty:
        return out
    types = (df["mixcr_type"].fillna("").astype(str) + ";"
             + df["trust_type"].fillna("").astype(str)).str.split(";")
    for chain in chains:
        mask = types.map(lambda tokens, c=chain: c in tokens)
        out[chain] = set(df.loc[mask, "CDR3aa"])
    return out


def all_aa_set(clono_path):
    """set of ALL unique non-NaN CDR3aa in the table, regardless of chain.

    Matches predict.py load_aa_set (empty file -> empty set).
    """
    df = pd.read_table(clono_path)
    if df.empty:
        return set()
    return set(df["CDR3aa"].dropna())


def _attach_fisher(df):
    """Add disease-oriented Fisher OR/P columns to one chain's count table.

    Contingency for group g in {UC, CD}: [[k_g, n_g-k_g], [k_HT, n_HT-k_HT]]
    (two-sided scipy.stats.fisher_exact), so OR > 1 = enriched in disease.
    Identical count quadruples are tested only once and merged back.
    """
    for group in GROUPS:
        key_cols = [group, "HT", group + "_", "HT_"]
        dup = df[key_cols].drop_duplicates()
        ors, ps = [], []
        for _, row in dup.iterrows():
            or_val, p_val = fisher_exact([[row[group], row[group + "_"]],
                                          [row["HT"], row["HT_"]]])
            ors.append(or_val)
            ps.append(p_val)
        dup["OR_{}".format(group)] = ors
        dup["P_{}".format(group)] = ps
        df = pd.merge(df, dup, on=key_cols, how="left")
    cols = (["CDR3aa"]
            + [c for lbl in LABELS for c in (lbl, lbl + "_")]
            + [c for g in GROUPS for c in ("OR_{}".format(g), "P_{}".format(g))])
    return df[cols].sort_values("UC", ascending=False).reset_index(drop=True)


def fisher_tables(cohort, chains=CHAINS):
    """Per chain, one Fisher exact table over the cohort (sample-level).

    k_X = number of DISTINCT SAMPLES in label X carrying that CDR3aa; n_X =
    number of samples in label X with a clonotype file. patient-level
    within-patient non-independence is handled upstream by resample_patients.py
    (one specimen per patient), so analysis scripts all run sample-level.

    Carrier counting keeps, per CDR3aa, the index of the last carrier unit
    seen while iterating the cohort sorted by (label, sample), so each unit
    is counted at most once per CDR3aa without dense matrices. Samples whose
    clonotype file cannot be read are warned to stderr and treated as carrying
    no CDR3aa (their unit still counts in n_X).

    Returns {chain: DataFrame sorted by UC desc} with columns
    CDR3aa, HT, UC, CD, HT_, UC_, CD_, OR_UC, P_UC, OR_CD, P_CD.
    """
    counts = {chain: {lbl: {} for lbl in LABELS} for chain in chains}
    last_seen = {chain: {} for chain in chains}
    n_units = {lbl: 0 for lbl in LABELS}
    prev_label, cur_unit, cur_idx = None, None, -1
    ordered = cohort.sort_values(["label", "sample"])
    for _, row in ordered.iterrows():
        label = row["label"]
        if label not in LABELS:
            continue
        unit_id = row["sample"]
        if label != prev_label or unit_id != cur_unit:
            prev_label, cur_unit = label, unit_id
            cur_idx += 1
            n_units[label] += 1
        try:
            aa_by_chain = chain_aa_by_sample(row["clono_path"], chains=chains)
        except Exception as exc:
            print("WARNING: {}: failed to read clonotypes ({}): {}".format(
                row["sample"], row["clono_path"], exc), file=sys.stderr)
            aa_by_chain = {chain: set() for chain in chains}
        for chain in chains:
            reg = counts[chain][label]
            seen = last_seen[chain]
            for cdr3 in aa_by_chain[chain]:
                if seen.get(cdr3) != cur_idx:
                    seen[cdr3] = cur_idx
                    reg[cdr3] = reg.get(cdr3, 0) + 1
    tables = {}
    for chain in chains:
        all_aa = set()
        for lbl in LABELS:
            all_aa.update(counts[chain][lbl])
        df = pd.DataFrame({"CDR3aa": list(all_aa)})
        for lbl in LABELS:
            df[lbl] = df["CDR3aa"].map(counts[chain][lbl]).fillna(0).astype(int)
            df[lbl + "_"] = n_units[lbl] - df[lbl]
        tables[chain] = _attach_fisher(df)
    return tables


def build_sets(fisher_df, p_thresh=0.01, or_thresh=1.5, min_count=5):
    """Threshold a Fisher table into the universe/UC/CD CDR3aa sets.

    universe = set of ALL CDR3aa observed for this chain in the training
    cohort (the chain universe; matches the fisher_sets.json semantics
    where universe ~2M rows).
    UC/CD = CDR3aa with (P <= p_thresh) & (disease-oriented OR > or_thresh)
    & (carrier count >= min_count).
    """
    sets = {"universe": set(fisher_df["CDR3aa"])}
    for group in GROUPS:
        mask = ((fisher_df["P_{}".format(group)] <= p_thresh)
                & (fisher_df["OR_{}".format(group)] > or_thresh)
                & (fisher_df[group] >= min_count))
        sets[group] = set(fisher_df.loc[mask, "CDR3aa"])
    return sets


def coverage_table(cohort, sets_by_chain, chains=CHAINS):
    """Per-sample coverage features against the Fisher reference sets.

    For each sample and chain: n = |all_aa_set(sample) ∩ universe|,
    k_UC = |all_aa_set ∩ UC-set|, k_CD = |all_aa_set ∩ CD-set|. IMPORTANT:
    the sample side uses ALL CDR3aa (unfiltered by chain), matching
    predict.py sample_features. Read failures and empty files: warn to
    stderr and give the sample NaN coverage (do not crash the run).

    -> DataFrame[sample, patient, label, {chain}_universe, {chain}_UC, {chain}_CD]
    """
    cols = ["sample", "patient", "label"]
    for chain in chains:
        cols += ["{}_universe".format(chain), "{}_UC".format(chain),
                 "{}_CD".format(chain)]
    rows = []
    for _, row in cohort.iterrows():
        out = {"sample": row["sample"], "patient": row["patient"],
               "label": row["label"]}
        aa = None
        try:
            aa = all_aa_set(row["clono_path"])
            if not aa:
                print("WARNING: {}: clonotype file empty -> NaN coverage".format(
                    row["sample"]), file=sys.stderr)
                aa = None
        except Exception as exc:
            print("WARNING: {}: failed to read clonotypes ({}): {}".format(
                row["sample"], row["clono_path"], exc), file=sys.stderr)
        for chain in chains:
            if aa is None:
                out["{}_universe".format(chain)] = float("nan")
                out["{}_UC".format(chain)] = float("nan")
                out["{}_CD".format(chain)] = float("nan")
            else:
                chain_sets = sets_by_chain[chain]
                out["{}_universe".format(chain)] = len(aa & chain_sets["universe"])
                out["{}_UC".format(chain)] = len(aa & chain_sets["UC"])
                out["{}_CD".format(chain)] = len(aa & chain_sets["CD"])
        rows.append(out)
    return pd.DataFrame(rows, columns=cols)


def fit_bb(n, k):
    """Beta-Binomial MLE (integer a, float b).

    VERBATIM algorithm from src_fold/util/ig_cover_model_full.py fit_bb:
    grid a in 1..29 against b_grid = arange(1, 3001, 100), then up to 2000
    gradient-ascent steps on b (digamma) for a in {a*-1, a*, a*+1}; keep the
    best log-likelihood.
    """
    n = np.asarray(n, dtype=float)
    k = np.asarray(k, dtype=float)
    N = n.size
    best = (-np.inf, 1, 1.0)
    for a in range(1, 30):
        ll = -N * betaln(a, B_GRID) + betaln(k[:, None] + a, (n - k)[:, None] + B_GRID[None, :]).sum(axis=0)
        j = int(np.argmax(ll))
        if ll[j] > best[0]:
            best = (float(ll[j]), a, float(B_GRID[j]))
    best_l, best_a, best_b = best
    for a in (max(1, best_a - 1), best_a, best_a + 1):
        b = best_b
        l_run, b_run = -np.inf, b
        for _ in range(2000):
            gb = -N * (digamma(b) - digamma(a + b)) + (digamma(n - k + b) - digamma(n + a + b)).sum()
            b_new = b + gb
            if not (1e-9 < b_new < 1e9):
                break
            b = b_new
            ll = -N * betaln(a, b) + betaln(k + a, n - k + b).sum()
            if ll > l_run:
                l_run, b_run = float(ll), b
        if l_run > best_l:
            best_l, best_a, best_b = l_run, a, b_run
    return best_a, best_b


def score_F(n, k, N0, a0, b0, N1, a1, b1):
    """Log Bayes factor: disease model (index 1) vs healthy model (index 0).

    Identical formula to predict.py score_F (accepts scalars or arrays for
    n/k via vectorized betaln).
    """
    return (math.log(N1 + 1) - math.log(N0 + 1) + betaln(a0, b0) - betaln(a1, b1)
            + betaln(k + a1, n - k + b1) - betaln(k + a0, n - k + b0))


def fit_params(cov_df, chains=CHAINS):
    """Fit per chain x key in {HT_UC, UC, HT_CD, CD} Beta-Binomial params.

    For key K: group = samples with label LABEL_OF[K]; k column is
    {chain}_UC for K in {HT_UC, UC} and {chain}_CD for K in {HT_CD, CD};
    n column is {chain}_universe. Rows with NaN coverage for that chain are
    dropped before fitting.

    -> params dict EXACTLY the schema of params_full.json:
    {chain: {key: {"a": int, "b": float, "N": int(len(group))}}}
    """
    params = {}
    for chain in chains:
        params[chain] = {}
        for key in KEYS:
            grp = cov_df.loc[cov_df["label"] == LABEL_OF[key]]
            n_col = "{}_universe".format(chain)
            k_col = "{}_{}".format(chain, COMPARE_OF[key])
            grp = grp.dropna(subset=[n_col, k_col])
            a, b = fit_bb(grp[n_col].to_numpy(), grp[k_col].to_numpy())
            params[chain][key] = {"a": a, "b": b, "N": int(len(grp))}
    return params


def compute_scores(cov_df, params, chains=CHAINS):
    """Per-chain and total BCR scores for every sample in cov_df.

    CRITICAL: each chain uses ITS OWN n/k columns with ITS OWN params.
    For group in (UC, CD): score_{chain}_{group} = score_F(n_{chain},
    k_{chain}_{group}, HT_<group> params, <group> params);
    score_{group} = sum over chains.

    -> copy of cov_df with score_{chain}_{group} and score_{group} columns.
    """
    out = cov_df.copy()
    for group in GROUPS:
        total = None
        for chain in chains:
            ht = params[chain]["HT_{}".format(group)]
            dz = params[chain][group]
            sc = score_F(out["{}_universe".format(chain)],
                         out["{}_{}".format(chain, group)],
                         ht["N"], ht["a"], ht["b"],
                         dz["N"], dz["a"], dz["b"])
            out["score_{}_{}".format(chain, group)] = sc
            total = sc if total is None else total + sc
        out["score_{}".format(group)] = total
    return out
