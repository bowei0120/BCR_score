import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd
from scipy.special import betaln

BASE = Path(__file__).resolve().parent.parent
CHAINS = ["IGL", "IGK"]
GROUPS = ["UC", "CD"]
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
    """Load per-chain cover/UC/CD CDR3aa reference sets (sorted lists -> sets)."""
    with open(path) as f:
        payload = json.load(f)
    return {chain: {key: set(vals) for key, vals in payload[chain].items()}
            for chain in CHAINS}


def load_aa_set(clono_path):
    """Unique non-NaN CDR3aa of a clonotype TSV, or None if the file is empty."""
    df = pd.read_table(clono_path)
    if df.empty:
        return None
    return set(df["CDR3aa"].dropna())


def score_F(n, k, N0, a0, b0, N1, a1, b1):
    """Log Bayes factor: disease model (index 1) vs healthy model (index 0)."""
    return (math.log(N1 + 1) - math.log(N0 + 1) + betaln(a0, b0) - betaln(a1, b1)
            + betaln(k + a1, n - k + b1) - betaln(k + a0, n - k + b0))


def chain_scale(chain, n, k, params, group):
    """F scale for one chain: healthy (HT_<group>) vs disease (<group>) model."""
    ht = params[chain]["HT_{}".format(group)]
    dz = params[chain][group]
    return score_F(n, k, ht["N"], ht["a"], ht["b"], dz["N"], dz["a"], dz["b"])


def sample_features(aa_set, refs):
    """Count sample CDR3aa hits against the cover/UC/CD reference sets per chain."""
    feats = {}
    for chain in CHAINS:
        hits = {g: sorted(aa_set & refs[chain][g]) for g in ("cover", "UC", "CD")}
        feats[chain] = {
            "n": len(hits["cover"]),
            "k_UC": len(hits["UC"]),
            "k_CD": len(hits["CD"]),
        }
    return feats


def sample_scales(feats, params):
    """Per-chain and total F scales (UC/CD) for a sample."""
    scales = {}
    for group in GROUPS:
        total = 0.0
        for chain in CHAINS:
            sc = chain_scale(chain, feats[chain]["n"], feats[chain]["k_" + group],
                             params, group)
            scales["{}_{}".format(chain, group)] = sc
            total += sc
        scales[group] = total
    return scales


def build_row(sample, feats, scales):
    """Output row dict; feats=None -> NA row (NaN counts/scales, empty hits)."""
    row = {"sample": sample}
    for chain in CHAINS:
        for group in GROUPS:
            row["scale_{}_{}".format(chain, group)] = float("nan")
        row["n_{}".format(chain)] = float("nan")
        row["k_{}_UC".format(chain)] = float("nan")
        row["k_{}_CD".format(chain)] = float("nan")
    row["scale_UC"], row["scale_CD"] = float("nan"), float("nan")
    if feats is not None:
        for chain in CHAINS:
            row["n_{}".format(chain)] = feats[chain]["n"]
            row["k_{}_UC".format(chain)] = feats[chain]["k_UC"]
            row["k_{}_CD".format(chain)] = feats[chain]["k_CD"]
        for key, val in scales.items():
            row["scale_{}".format(key)] = val
    return row


def main():
    parser = argparse.ArgumentParser(
        description="Compute disease F scales (UC/CD) and covered sequences for "
                    "new samples from a sample table (columns: sample, clono_path)")
    parser.add_argument("--table", type=str, required=True,
                        help="TSV with columns 'sample' and 'clono_path'")
    parser.add_argument("--fisher-json", type=str,
                        default=str(BASE / "predict" / "fisher_sets.json"),
                        help="Precomputed CDR3 reference sets JSON")
    parser.add_argument("--params", type=str,
                        default=str(BASE / "predict" / "params_full.json"),
                        help="Pre-trained Beta-Binomial params JSON")
    parser.add_argument("--out", type=str, default=None,
                        help="Write TSV to this file (default: stdout)")
    args = parser.parse_args()

    params = load_params(args.params)
    refs = load_fisher_json(Path(args.fisher_json))

    table = pd.read_table(args.table, comment="#")
    for col in ("sample", "clono_path"):
        if col not in table.columns:
            print("ERROR: --table must have a '{}' column".format(col), file=sys.stderr)
            sys.exit(2)

    rows = []
    for _, row in table.iterrows():
        sample = row["sample"]
        try:
            aa_set = load_aa_set(Path(row["clono_path"]))
            if aa_set is None:
                print("WARNING: {}: clonotype file empty -> NA row".format(sample),
                      file=sys.stderr)
                rows.append(build_row(sample, None, None))
                continue
            feats = sample_features(aa_set, refs)
            scales = sample_scales(feats, params)
            rows.append(build_row(sample, feats, scales))
        except Exception as exc:
            print("WARNING: {}: failed to score ({}): {}".format(
                sample, row["clono_path"], exc), file=sys.stderr)
            rows.append(build_row(sample, None, None))

    out_df = pd.DataFrame(rows, columns=OUT_COLS)
    if args.out:
        out_df.to_csv(args.out, sep="\t", index=None)
    else:
        sys.stdout.write(out_df.to_csv(sep="\t", index=None))


if __name__ == "__main__":
    main()
