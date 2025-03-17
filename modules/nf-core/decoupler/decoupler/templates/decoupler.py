#!/usr/bin/env python3
import os
import pandas as pd
import argparse
import shlex

os.environ["NUMBA_CACHE_DIR"] = "./tmp"

import decoupler as dc

methods = ['aucell', 'gsea', 'gsva', 'mdt', 'mlm', 'ora', 'udt',
    'ulm', 'viper', 'wmean', 'wsum']

mat = pd.read_csv("${mat}", sep="\t", index_col=0)
net = pd.read_csv("${net}", sep="\t")

# Parse external arguments passed as a string, e.g. "--min_n 1"
raw_args = "${args}"
args_list = shlex.split(raw_args)
parser = argparse.ArgumentParser()
# Define expected external arguments. Here we only define min_n; add more as needed.
parser.add_argument("--min_n", type=int, default=1, help="Minimum n value")
parsed_args, unknown = parser.parse_known_args(args_list)

# Build the arguments dictionary for decoupler
# If you need method-specific arguments, adjust accordingly.
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