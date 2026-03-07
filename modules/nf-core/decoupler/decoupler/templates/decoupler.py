#!/usr/bin/env python3
import os

# This must execute before decoupler is imported, so the env vars are visible when Numba decides whether to cache
os.environ["NUMBA_DISABLE_CACHE"] = "1"
os.environ["NUMBA_CACHE_DIR"] = "./tmp"
os.environ["MPLCONFIGDIR"] = "./tmp"

import argparse
import shlex
import sys

import decoupler as dc  # noqa: E402
import matplotlib.pyplot as plt
import pandas as pd  # noqa: E402

mat = pd.read_csv("${mat}", sep="\t", index_col=0)
net = pd.read_csv("${net}", sep="\t")
annot: str = "${annot}"


def parse_ext_args(args_string: str):
    """
    Parse external arguments passed from Nextflow's ext.args.
    Expected arguments:
      --min_n <int>
      --transpose <str> (TRUE or FALSE)
      --contrast <str> (optional, e.g., treatment_vs_control)
      --column <str> (Column name to use for transposition; default: log2FoldChange)
      --ensembl_ids <str> (TRUE to convert ENSEMBL IDs to gene symbols, FALSE to skip)
      --methods <str> (Comma-separated list of methods to use (e.g., 'mlm,ulm'))
    """
    if args_string == "null":
        args_string = ""
    args_list = shlex.split(args_string)
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_n", type=int, default=1, help="Minimum n value")
    parser.add_argument("--transpose", type=str, default="FALSE", help="Transpose DESeq2 data if TRUE")
    parser.add_argument("--column", type=str, default="log2FoldChange", help="Column name to use for transposition")
    parser.add_argument("--ensembl_ids", type=str, default="FALSE", help="Convert ENSEMBL IDs to gene symbols if TRUE")
    parser.add_argument(
        "--methods", type=str, default="ulm", help="Comma-separated list of methods to use (e.g., 'mlm,ulm')"
    )
    parser.add_argument("--features_id_col", type=str, default="gene_id", help="Column name for feature IDs")
    parser.add_argument("--features_symbol_col", type=str, default="gene_name", help="Column name for feature symbols")
    return parser.parse_args(args_list)


# Parse external arguments
raw_args = "${task.ext.args}"
parsed_args = parse_ext_args(raw_args)
methods = [m.strip() for m in parsed_args.methods.split(",") if m.strip()]

if parsed_args.ensembl_ids.upper() == "TRUE":
    try:
        if not os.path.exists("${annot}"):
            raise FileNotFoundError(f"Annotation file not found: ${annot}")

        annot_df = pd.read_csv("${annot}", sep="\t")

        required_cols = {parsed_args.features_id_col, parsed_args.features_symbol_col}

        missing = required_cols - set(annot_df.columns)
        if missing:
            raise ValueError(
                f"Missing required columns in annotation file: {missing}. Available columns: {list(annot_df.columns)}"
            )

        gene_mapping = dict(zip(annot_df[parsed_args.features_id_col], annot_df[parsed_args.features_symbol_col]))
        new_index = [gene_mapping.get(ens, None) for ens in mat.index]
        mat.index = new_index
        mat = mat[mat.index.notnull()]
        mat = mat[~mat.index.duplicated(keep="first")]
    except Exception as e:
        print(f"ERROR: Failed to process annotation file: {e}")
        sys.exit(1)

if parsed_args.transpose.upper() == "TRUE":
    mat = mat[[parsed_args.column]].T.rename(index={parsed_args.column: "${meta.id}"})

parsedargs = {"args": {}}
parsedargs["min_n"] = parsed_args.min_n


results = dc.decouple(mat=mat, net=net, methods=methods, **parsedargs)

for result in results:
    # Save table
    results[result].to_csv("${task.ext.prefix}" + "_" + result + "_decoupler.tsv", sep="\t")
    contrast_name = results[result].index[0]
    plt.figure(figsize=(8, 6))
    dc.plot_barplot(results[result], contrast_name, top=25, vertical=False)
    plt.savefig("${task.ext.prefix}" + "_" + result + "_decoupler_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

## VERSIONS FILE
with open("versions.yml", "a") as version_file:
    version_file.write('"${task.process}":' + "\\n")
    version_file.write("decoupler-py: " + dc.__version__ + "\\n")
