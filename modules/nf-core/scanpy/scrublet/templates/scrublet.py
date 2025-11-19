#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col ?: ''}"

kwargs = {}
if batch_col and adata.obs[batch_col].nunique() > 1:
    kwargs["batch_key"] = batch_col

sc.pp.scrublet(adata, **kwargs)

adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
df = adata.obs[["predicted_doublet"]]
df.columns = ["${prefix}"]
df.to_pickle("${prefix}.pkl")

adata = adata[~adata.obs["predicted_doublet"]].copy()

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
