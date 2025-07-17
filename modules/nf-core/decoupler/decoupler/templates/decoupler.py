#!/usr/bin/env python3
import os
import argparse
import shlex
import sys
import gzip

os.environ["NUMBA_CACHE_DIR"] = "./tmp"
os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_DISABLE_CACHE"] = "1"

import numba
numba.config.DISABLE_CACHE = True

import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt


mat = pd.read_csv("${mat}", sep="\t", index_col=0)
net = pd.read_csv("${net}", sep="\t")

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
    parser.add_argument("--methods", type=str, default = "ulm", help="Comma-separated list of methods to use (e.g., 'mlm,ulm')")
    return parser.parse_args(args_list)

def parse_gtf(gtf_file: str):
    """
    Parse an optional GTF file to create a mapping of ENSEMBL gene IDs to gene symbols (required to use Progeny data).
    """
    mapping = {}
    opener = gzip.open if gtf_file.endswith('.gz') else open
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attributes_field = fields[8]
            attributes = {}
            for attr in attributes_field.split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split(" ", 1)
                if len(parts) != 2:
                    continue
                key, value = parts
                attributes[key] = value.replace('"', '').strip()
            gene_id = attributes.get("gene_id")
            gene_symbol = attributes.get("gene_name") or attributes.get("gene_symbol") or attributes.get("external_gene_name")
            if gene_id and gene_symbol:
                mapping[gene_id] = gene_symbol
    return mapping

# Parse external arguments
raw_args = "${task.ext.args}"
parsed_args = parse_ext_args(raw_args)
methods = [m.strip() for m in parsed_args.methods.split(",") if m.strip()]

if parsed_args.ensembl_ids.upper() == "TRUE":
    try:
        gene_mapping = parse_gtf("${gtf}")
        new_index = [gene_mapping.get(ens, None) for ens in mat.index]
        mat.index = new_index
        mat = mat[mat.index.notnull()]
        mat = mat[~mat.index.duplicated(keep='first')]
    except Exception as e:
        print("ERROR: Failed to parse GTF:", e)
        sys.exit(1)

if parsed_args.transpose.upper() == "TRUE":
    mat = mat[[parsed_args.column]].T.rename(index={parsed_args.column: "${meta.id}"})

parsedargs = {'args': {}}
parsedargs['min_n'] = parsed_args.min_n


results = dc.decouple(
    mat=mat,
    net=net,
    methods=methods,
    **parsedargs
)

for result in results:
    # Save table
    results[result].to_csv("${task.ext.prefix}" + "_" + result + "_decoupler.tsv", sep="\t")
    contrast_name = results[result].index[0]
    plt.figure(figsize=(8, 6))
    dc.plot_barplot(results[result], contrast_name , top=25, vertical=False)
    plt.savefig("${task.ext.prefix}" + "_" + result + "_decoupler_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

## VERSIONS FILE
with open('versions.yml', 'a') as version_file:
    version_file.write('"${task.process}":' + "\\n")
    version_file.write("decoupler-py: " + dc.__version__ + "\\n")
