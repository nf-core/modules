#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

input_file = "${anndata}"
output_file = "${output_file}"

input_suffix = Path(input_file).suffix
if input_suffix == ".h5ad":
    adata = ad.read_h5ad(input_file)
elif input_suffix == ".zarr":
    adata = ad.read_zarr(input_file)
else:
    raise ValueError(f"Unsupported AnnData input format: {input_suffix}")

prefix = "${prefix}"
key_added = "${key_added}"

# Run PCA
sc.pp.pca(adata, random_state=0, key_added=key_added)

# Round to 8 decimal places
# This ensures hashes are stable
adata.obsm[key_added] = np.round(adata.obsm[key_added], 8)

output_suffix = Path(output_file).suffix
if output_suffix == ".h5ad":
    adata.write_h5ad(output_file)
elif output_suffix == ".zarr":
    adata.write_zarr(output_file)
else:
    raise ValueError(f"Unsupported AnnData output format: {output_suffix}")

df = pd.DataFrame(adata.obsm[key_added], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
