#!/usr/bin/env python3

import os

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import platform
from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col}"

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

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
