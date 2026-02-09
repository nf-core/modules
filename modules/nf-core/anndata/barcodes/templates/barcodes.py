#!/usr/bin/env python3

import platform

import anndata as ad
import pandas as pd
import yaml

df = pd.read_csv("${barcodes}", header=None)
adata = ad.read_h5ad("${h5ad}")

adata = adata[df[0].values]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "anndata": ad.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
