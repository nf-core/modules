#!/usr/bin/env python3
import os
import argparse
import shlex
import sys

os.environ["NUMBA_CACHE_DIR"] = "./tmp"
os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_DISABLE_CACHE"] = "1"

import numba
numba.config.DISABLE_CACHE = True

import pandas as pd
import scanpy as sc
import decoupler as dc

methods = ['aucell', 'gsea', 'gsva', 'mdt', 'mlm', 'ora', 'udt',
    'ulm', 'viper', 'wmean', 'wsum']

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
      --species <str> (Species to use for BioMart annotations; default: mmusculus)
      --ensembl_ids <str> (TRUE to convert ENSEMBL IDs to gene symbols, FALSE to skip)
    """
    if args_string == "null":
        args_string = ""
    args_list = shlex.split(args_string)
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_n", type=int, default=1, help="Minimum n value")
    parser.add_argument("--transpose", type=str, default="FALSE", help="Transpose DESeq2 data if TRUE")
    parser.add_argument("--column", type=str, default="log2FoldChange", help="Column name to use for transposition")
    parser.add_argument("--species", type=str, default="mmusculus", help="Species for BioMart annotations")
    parser.add_argument("--ensembl_ids", type=str, default="FALSE", help="Convert ENSEMBL IDs to gene symbols if TRUE")
    return parser.parse_args(args_list)

# Parse external arguments
raw_args = "${task.ext.args}"
parsed_args = parse_ext_args(raw_args)

if parsed_args.ensembl_ids.upper() == "TRUE":
    try:
        annot = sc.queries.biomart_annotations(
            "mmusculus",
            ["ensembl_gene_id", "external_gene_name"],
            host="http://www.ensembl.org",
            use_cache=False
        ).set_index("ensembl_gene_id")
    except Exception as e:
        print("ERROR: Failed to retrieve annotations from BioMart:", e)
        sys.exit(1)

    new_index = [annot.loc[ens, "external_gene_name"] if ens in annot.index else None for ens in mat.index]
    mat.index = new_index

    mat = mat[mat.index.notnull()]
    mat = mat[~mat.index.duplicated(keep='first')]

if parsed_args.transpose.upper() == "TRUE":
    mat = mat[[parsed_args.column]].T.rename(index={parsed_args.column: "${meta.contrast}"})

parsedargs = {'args': {}}
parsedargs['min_n'] = parsed_args.min_n


results = dc.decouple(
    mat=mat,
    net=net,
    **parsedargs
)

for result in results:
    results[result].to_csv(result + "__decoupler.tsv", sep="\t")

## VERSIONS FILE
with open('versions.yml', 'a') as version_file:
    version_file.write('"${task.process}":' + "\\n")
    version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")