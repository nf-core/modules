#!/usr/bin/env python3

import platform

import anndata as ad
import pandas as pd


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


df = pd.read_csv("${barcodes}", header=None)
adata = ad.read_h5ad("${h5ad}")

adata = adata[df[0].values]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "anndata": ad.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
