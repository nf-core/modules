#!/usr/bin/env python3

import os
import platform
import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import scanpy.external as sce


adata = sc.read_h5ad("${input_h5ad}")
columns = "${cell_hashing_columns.join(' ')}".split()
columns_str = [str(x) for x in columns]
sce.pp.hashsolo(adata, columns_str, priors=[float(prior) for prior in "${priors.join(',')}".split(',')])

adata.write("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
