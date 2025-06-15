#!/usr/bin/env python3

import json

import scanpy as sc
import scanpy.external as sce

# in case there are no optional arguments passed in
try:
    params = json.loads("${params}") if "${params}" else {}
except json.JSONDecodeError:
    params = {}

adata = sc.read_h5ad("${input_h5ad}")
columns = "${cell_hashing_columns}"[1:-1].split(", ")
columns_str = [str(x) for x in columns]
sce.pp.hashsolo(adata, columns_str,"${params}")
adata.write("test/${prefix}_hashsolo.h5ad")
