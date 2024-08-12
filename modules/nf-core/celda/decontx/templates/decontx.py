#!/usr/bin/env python3

import anndata as ad
import anndata2ri
import rpy2
import rpy2.robjects as ro
import platform
import os
celda = ro.packages.importr('celda')

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

adata = ad.read_h5ad("${filtered}")
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
sce = anndata2ri.py2rpy(adata)

kwargs = {}

if '${batch_col}' in adata.obs and len(adata.obs['${batch_col}'].unique()) > 1:
    kwargs['batch'] = adata.obs['${batch_col}'].tolist()

unfiltered_path = "${unfiltered}"
if os.path.exists(unfiltered_path):
    adata_unfiltered = ad.read_h5ad(unfiltered_path)
    if "counts" not in adata_unfiltered.layers:
        adata_unfiltered.layers["counts"] = adata_unfiltered.X.copy()
    kwargs["background"] = anndata2ri.py2rpy(adata_unfiltered)

corrected = celda.decontX(sce, **kwargs)
counts = celda.decontXcounts(corrected)

adata.X = anndata2ri.rpy2py(counts).T
adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": rpy2.__version__,
        "celda": celda.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
