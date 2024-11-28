#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_CACHE_DIR"] = "./tmp"

import anndata as ad
import doubletdetection


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
        spaces = "    " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = ad.read_h5ad("${h5ad}")

clf = doubletdetection.BoostClassifier()
doublets = clf.fit(adata.X).predict()
scores = clf.doublet_score()

adata.obs["doublet"] = [label == 1 for label in doublets]
adata.obs["doublet_score"] = scores

adata.write_h5ad("${prefix}.h5ad")

df = adata.obs[["doublet"]]
df.columns = ["${prefix}"]
df.to_pickle("${prefix}.pkl")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "doubletdetection": doubletdetection.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
