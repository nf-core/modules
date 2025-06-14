import argparse
import json

import scanpy as sc
import scanpy.external as sce

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--params", default="{}")

args = parser.parse_args()

# in case there are no optional arguments passed in
try:
    params = json.loads(args.params) if args.params else {}
except json.JSONDecodeError:
    params = {}

adata = sc.read_h5ad(args.input)
sce.pp.hashsolo(adata, **params)
adata.write(args.output)
