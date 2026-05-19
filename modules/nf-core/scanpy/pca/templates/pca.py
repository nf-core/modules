#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

input_file = "${h5ad}"
output_file = "${output_file}"
adata = ad.read_zarr(input_file) if input_file.endswith(".zarr") else sc.read_h5ad(input_file)
prefix = "${prefix}"
key_added = "${key_added}"

# Run PCA
sc.pp.pca(adata, random_state=0, key_added=key_added)

# Round to 8 decimal places
# This ensures hashes are stable
adata.obsm[key_added] = np.round(adata.obsm[key_added], 8)

if output_file.endswith(".zarr"):
    adata.write_zarr(output_file)
else:
    adata.write_h5ad(output_file)
df = pd.DataFrame(adata.obsm[key_added], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
