#!/usr/bin/env python3

import argparse
import os
import pathlib

import pandas as pd


def join_all(series):
    """Join a column of values into a single ';'-separated string."""
    return ";".join(series)


def parse_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine the TCR CDR3 result of mixcr and trust4"
    )
    parser.add_argument(
        "-m",
        "--mixcr-path",
        required=True,
        help="path to the mixcr clonotypes file",
    )
    parser.add_argument(
        "-t",
        "--trust4-path",
        required=True,
        help="path to the trust4 report file",
    )
    parser.add_argument(
        "-i",
        "--sample-id",
        required=True,
        help="sample id",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="output directory",
    )
    return parser.parse_args()


def combine_results(mixcr_path, trust4_path):
    """Merge the TCR CDR3 results of mixcr and trust4.

    Returns:
        data: merged clonotype table
        gene_type_dict: count of clonotypes per gene type
    """
    all_gene_type = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]

    # mixcr table
    mixcr = pd.read_csv(mixcr_path, sep="\t", usecols=[
        "cloneCount", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "nSeqCDR3", "aaSeqCDR3"
    ])
    mixcr.columns = ["mixcr_count", "mixcr_V", "mixcr_D", "mixcr_J", "CDR3nt", "CDR3aa"]
    mixcr["mixcr_tpye"] = mixcr.apply(lambda x: ";".join(
        [gene_type for gene_type in all_gene_type
         if gene_type in x["mixcr_V"] or gene_type in x["mixcr_J"]]), axis=1)
    mixcr["mixcr_count"] = mixcr["mixcr_count"].astype(int).astype(str)
    mixcr = mixcr.fillna("")
    mixcr = mixcr.groupby(by=["CDR3nt", "CDR3aa"], as_index=False).agg(join_all)

    # trust table
    trust = pd.read_csv(trust4_path, sep="\t", usecols=[
        "#count", "CDR3nt", "CDR3aa", "V", "D", "J"
    ]).rename(columns={"#count": "trust_count", "V": "trust_V", "D": "trust_D", "J": "trust_J"})
    trust["trust_tpye"] = trust.apply(lambda x: ";".join(
        [gene_type for gene_type in all_gene_type
         if gene_type in x["trust_V"] or gene_type in x["trust_J"]]), axis=1)
    trust["trust_count"] = trust["trust_count"].astype(int).astype(str)
    trust = trust.fillna("")
    trust = trust.groupby(by=["CDR3nt", "CDR3aa"], as_index=False).agg(join_all)

    # merge
    data = pd.merge(mixcr, trust, how="outer", on=["CDR3nt", "CDR3aa"])
    valid_cdr3 = lambda s: all(a in "ARNDCEQGHILKMFPSTWYV" for a in s)
    data = data.loc[data["CDR3aa"].apply(valid_cdr3), :]
    data = data.fillna("")

    gene_type_dict = {}
    for gene_type in all_gene_type:
        data[gene_type] = data.apply(lambda x: 1 if (
            gene_type in x["mixcr_tpye"] or gene_type in x["trust_tpye"]) else 0, axis=1)
        gene_type_dict[gene_type] = sum(data[gene_type])
    gene_type_dict["ALLnt"] = len(data.index)
    gene_type_dict["ALLaa"] = data["CDR3aa"].nunique()

    data = data[["CDR3nt", "CDR3aa", "mixcr_V", "mixcr_D", "mixcr_J",
                 "mixcr_count", "mixcr_tpye", "trust_V", "trust_D", "trust_J",
                 "trust_count", "trust_tpye"]]
    return data, gene_type_dict


def metrics_manifest(gene_type_dict, metrics_manifest_file):
    """Write the per-gene-type clonotype counts to a metrics manifest YAML."""
    with open(metrics_manifest_file, "w") as f:
        f.write("---\nmetrics:\n")
        for gene_type, gene_num in gene_type_dict.items():
            f.write(f"\t{gene_type}:\t{gene_num}\n")


def main():
    args = parse_input()
    data, gene_type_dict = combine_results(args.mixcr_path, args.trust4_path)
    pathlib.Path(args.output_dir).mkdir(exist_ok=True, parents=True)
    data.to_csv(os.path.join(args.output_dir, f"{args.sample_id}.TCR_clonotypes.txt"), sep="\t", index=None)
    metrics_manifest(gene_type_dict, os.path.join(args.output_dir, f"{args.sample_id}_metrics_manifest.yaml"))


if __name__ == "__main__":
    main()
